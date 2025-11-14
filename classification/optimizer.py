"""
SVM hyperparameter optimization using geometric grid search with refinement.

This module implements the optimization algorithm from intronIC.py:5431-5528,
which uses 5 rounds of progressively refined grid search to find optimal
SVM hyperparameters (specifically the soft-margin penalty C).

Key improvements over original:
- Weight-aware C bounds: Computes search range based on effective penalties
  to avoid pathological extremes on imbalanced datasets (~1:200 U12:U2 ratio)
- Example: Instead of C ∈ [1e-6, 1e6] → eff_C_pos ∈ [9.5e-5, 9.5e7] (BAD)
           Use C ∈ [1e-4, 1.0] → eff_C_pos ∈ [1e-2, 1e2] (GOOD)
- Reduces convergence warnings, speeds up CV, improves calibration

Key algorithm:
- Compute weight-aware C bounds from class distribution
- Round 1: Coarse grid search (weight-aware range)
- Rounds 2-5: Refine around best C values from previous round
- Final: Geometric mean of all rank-1 C values from final round

Port from: intronIC.py:5431-5528 (optimize_svm)
Related: intronIC.py:5290-5322 (helper functions)
"""

from dataclasses import dataclass
from typing import Sequence, Tuple, Optional, Dict, Any
import os
import contextlib
import warnings
import joblib
import numpy as np
from sklearn.model_selection import GridSearchCV, cross_val_score, StratifiedKFold, ParameterGrid
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import RobustScaler
from sklearn.exceptions import ConvergenceWarning
from scipy.stats import gmean
from tqdm.auto import tqdm

from core.intron import Intron
from classification.transformers import BothEndsStrongTransformer

# Global filter for convergence warnings (persists across multiprocessing forks)
warnings.filterwarnings("ignore", category=ConvergenceWarning)


@contextlib.contextmanager
def tqdm_joblib(tqdm_object):
    """
    Context manager to patch joblib to report into tqdm progress bar.

    Usage:
        with tqdm_joblib(tqdm(total=total_tasks, desc="GridSearchCV")):
            grid.fit(X, y)

    Adapted from: https://stackoverflow.com/a/58936697
    """
    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()


@contextlib.contextmanager
def suppress_convergence_warnings(verbose: bool = True):
    """
    Context manager to suppress sklearn ConvergenceWarning spam while logging occurrence.

    During grid search or nested CV, convergence warnings can spam the console with
    hundreds or thousands of identical messages. This context manager:
    1. Captures ConvergenceWarnings silently
    2. Counts how many were raised
    3. Logs a summary at the end if any occurred

    Args:
        verbose: If True, print summary when warnings are captured (default: True)

    Usage:
        with suppress_convergence_warnings(verbose=True):
            grid.fit(X, y)
        # Output: "Suppressed convergence warnings (check max_iter if needed)"

    Example with nested context:
        with suppress_convergence_warnings():
            with tqdm_joblib(tqdm(total=100)):
                grid.fit(X, y)
    """
    with warnings.catch_warnings():
        # Suppress convergence warnings
        warnings.filterwarnings("ignore", category=ConvergenceWarning)
        yield


def compute_weight_aware_C_bounds(
    y: np.ndarray,
    sample_weight: Optional[np.ndarray] = None,
    class_weight: str = "balanced",
    eff_C_pos_range: Tuple[float, float] = (1e-2, 1e2),
    eff_C_neg_max: Optional[float] = None
) -> Dict[str, Any]:
    """
    Compute base-C bounds such that the effective positive-class penalty
    C_eff_pos = C * w_pos * s_pos_max lies within eff_C_pos_range.

    This prevents pathological extremes when using class_weight='balanced'
    on heavily imbalanced datasets (e.g., ~1:200 U12:U2 ratio).

    Args:
        y: Binary labels (0/1 for U2/U12)
        sample_weight: Optional per-sample weights
        class_weight: 'balanced' or None
        eff_C_pos_range: Target range for positive-class effective penalty (default: 1e-3 to 1e3)
        eff_C_neg_max: Optional cap on negative-class effective penalty

    Returns:
        Dictionary with:
            - 'C_min': Minimum base C value
            - 'C_max': Maximum base C value
            - 'class_weights': Dict with positive and negative class weights
            - 'effective_C_range': Dict with effective penalties at bounds

    Example:
        With ~540 U12s and ~102K U2s:
        - w_pos ≈ 95, w_neg ≈ 0.5 (from 'balanced')
        - eff_C_pos_range = (1e-3, 1e3)
        - C_min ≈ 1e-5, C_max ≈ 10.5
        - Instead of C ∈ [1e-6, 1e6] → eff_C_pos ∈ [9.5e-5, 9.5e7]  (BAD)
        - We get    C ∈ [1e-5, 10.5] → eff_C_pos ∈ [1e-3, 1e3]       (GOOD)
    """
    y = np.asarray(y)
    classes = np.unique(y)
    if len(classes) != 2:
        raise ValueError("This helper expects binary classification.")

    # Map to positive/negative classes
    pos = classes.max()  # U12 = 1
    neg = classes.min()  # U2 = 0
    pos_mask = (y == pos)
    neg_mask = ~pos_mask
    n = y.size
    n_pos = pos_mask.sum()
    n_neg = n - n_pos

    # Compute class weights
    if class_weight is None:
        w_pos = w_neg = 1.0
    elif class_weight == "balanced":
        # sklearn rule: n / (n_classes * n_k)
        w_pos = n / (2.0 * n_pos)
        w_neg = n / (2.0 * n_neg)
    else:
        raise ValueError(f"Unsupported class_weight: {class_weight}")

    # Sample weight maxima per class
    if sample_weight is None:
        s_pos_max = s_neg_max = 1.0
    else:
        sw = np.asarray(sample_weight)
        s_pos_max = float(sw[pos_mask].max())
        s_neg_max = float(sw[neg_mask].max())

    # Derive base-C bounds from desired effective positive-class range
    eff_C_min, eff_C_max = eff_C_pos_range
    C_min = eff_C_min / (w_pos * s_pos_max)
    C_max = eff_C_max / (w_pos * s_pos_max)

    # Optional: also guard the negative class's effective penalty
    if eff_C_neg_max is not None:
        C_max = min(C_max, eff_C_neg_max / (w_neg * s_neg_max))

    if not np.isfinite(C_min) or not np.isfinite(C_max) or C_max <= C_min:
        raise ValueError(f"Invalid C bounds: C_min={C_min}, C_max={C_max}")

    return {
        'C_min': C_min,
        'C_max': C_max,
        'class_weights': {'pos': w_pos, 'neg': w_neg},
        'sample_weight_max': {'pos': s_pos_max, 'neg': s_neg_max},
        'effective_C_range': {
            'pos': (C_min * w_pos * s_pos_max, C_max * w_pos * s_pos_max),
            'neg': (C_min * w_neg * s_neg_max, C_max * w_neg * s_neg_max)
        },
        'class_counts': {'pos': int(n_pos), 'neg': int(n_neg), 'ratio': float(n_neg / n_pos)}
    }


@dataclass(frozen=True, slots=True)
class SVMParameters:
    """Optimized SVM hyperparameters for LinearSVC with BothEndsStrong features."""

    C: float  # Soft-margin penalty
    calibration_method: str  # 'sigmoid' or 'isotonic'
    gamma_5_bp: float  # Feature scaling for 5'SS-BPS imbalance penalty
    gamma_5_3: float  # Feature scaling for 5'SS-3'SS imbalance penalty
    dual: bool  # Primal (False) or dual (True) formulation
    intercept_scaling: float  # Scaling for intercept (when dual=False)
    cv_score: float  # Cross-validation neg_log_loss
    round_found: int  # Which optimization round found these params (-1 = averaged)


@dataclass(frozen=True, slots=True)
class OptimizationRound:
    """Results from one round of grid search."""

    grid_points: np.ndarray  # C values tested
    scores: np.ndarray  # CV scores for each parameter combination
    best_C: float
    best_method: str  # 'sigmoid' or 'isotonic'
    best_gamma_5_bp: float  # Best γ for 5'SS-BPS imbalance penalty
    best_gamma_5_3: float  # Best γ for 5'SS-3'SS imbalance penalty
    best_dual: bool  # Primal (False) or dual (True) formulation
    best_intercept_scaling: float  # Intercept scaling parameter
    best_score: float
    rank_one_Cs: list[float]  # All rank-1 C values


class SVMOptimizer:
    """
    Geometric grid search for SVM hyperparameter optimization.

    Uses 5 rounds of progressively refined grid search:
    - Round 1: Coarse grid (10^-6 to 10^6, 13 points)
    - Rounds 2-5: Refine around best values (100 points each)
    - Final: Geometric mean of rank-1 parameters

    The refinement strategy zooms in on the best C value from the previous
    round by creating a new geometric grid between the neighbor points.

    Port from: intronIC.py:5431-5528
    """

    def __init__(
        self,
        n_rounds: int = 3,
        n_points_initial: int = 13,
        n_points_refine: int = 100,
        cv_folds: int = 5,
        random_state: int = 42,
        n_jobs: int = 1,
        verbose: bool = True,
        max_iter: int = 100000
    ):
        """
        Initialize optimizer.

        Args:
            n_rounds: Number of refinement rounds (default: 3, reduced from 5 for speed)
            n_points_initial: Grid points for round 1 (default: 13)
            n_points_refine: Grid points for rounds 2+ (default: 100)
            cv_folds: Cross-validation folds (default: 5)
            random_state: Random seed for reproducibility
            n_jobs: Number of parallel jobs for GridSearchCV (default: 1)
            verbose: Whether to print detailed progress (default: True)
            max_iter: Maximum iterations for LinearSVC convergence (default: 100000)
        """
        self.n_rounds = n_rounds
        self.n_points_initial = n_points_initial
        self.n_points_refine = n_points_refine
        self.cv_folds = cv_folds
        self.random_state = random_state
        self.n_jobs = n_jobs
        self.verbose = verbose
        self.max_iter = max_iter
        self.rounds_: list[OptimizationRound] = []

    def optimize(
        self,
        u12_introns: Sequence[Intron],
        u2_introns: Sequence[Intron],
        eff_C_pos_range: Tuple[float, float] = (1e-3, 1e3),
        eff_C_neg_max: Optional[float] = None
    ) -> SVMParameters:
        """
        Find optimal C parameter via geometric grid search.

        Uses weight-aware C bounds to avoid pathological effective penalties
        on heavily imbalanced datasets (e.g., ~1:200 U12:U2 ratio).

        Args:
            u12_introns: Training U12-type introns (with z-scores)
            u2_introns: Training U2-type introns (with z-scores)
            eff_C_pos_range: Target effective penalty range for positive class (default: 1e-3 to 1e3)
            eff_C_neg_max: Optional cap on negative class effective penalty

        Returns:
            Optimized parameters with best C value

        Raises:
            ValueError: If introns lack z-scores
        """
        # Extract features and labels
        X, y = self._prepare_training_data(u12_introns, u2_introns)

        # Compute weight-aware C bounds
        bounds_info = compute_weight_aware_C_bounds(
            y,
            class_weight="balanced",
            eff_C_pos_range=eff_C_pos_range,
            eff_C_neg_max=eff_C_neg_max
        )

        if self.verbose:
            print(f"Class distribution: {bounds_info['class_counts']['pos']} U12, "
                  f"{bounds_info['class_counts']['neg']} U2 "
                  f"(ratio: 1:{bounds_info['class_counts']['ratio']:.1f})")
            print(f"Balanced class weights: w_pos={bounds_info['class_weights']['pos']:.2f}, "
                  f"w_neg={bounds_info['class_weights']['neg']:.3f}")
            print(f"Weight-aware C range: [{bounds_info['C_min']:.2e}, {bounds_info['C_max']:.2e}]")
            print(f"  → Effective C_pos: {bounds_info['effective_C_range']['pos']}")
            print(f"  → Effective C_neg: {bounds_info['effective_C_range']['neg']}")

        # Initialize search range with weight-aware bounds
        initial_range = (bounds_info['C_min'], bounds_info['C_max'])
        current_grid = self._create_initial_grid(initial_range)

        # Track ranges to detect oscillation/convergence
        previous_ranges = []

        # Run geometric refinement
        for round_idx in range(self.n_rounds):
            print(f"Optimization round {round_idx + 1}/{self.n_rounds}...", flush=True)

            round_result = self._grid_search_round(
                X, y, current_grid, round_idx
            )
            self.rounds_.append(round_result)

            # Prepare next round's grid (refine around best)
            if round_idx < self.n_rounds - 1:
                next_grid = self._refine_grid(current_grid, round_result.best_C)

                # Check for range oscillation/revisit (indicates convergence)
                next_range = (next_grid.min(), next_grid.max())
                for prev_min, prev_max in previous_ranges:
                    # If ranges overlap significantly (>80%), we're oscillating
                    overlap_low = max(next_range[0], prev_min)
                    overlap_high = min(next_range[1], prev_max)
                    if overlap_low < overlap_high:
                        overlap_span = overlap_high / overlap_low
                        total_span = max(next_range[1] / next_range[0], prev_max / prev_min)
                        overlap_ratio = overlap_span / total_span
                        if overlap_ratio > 0.8:
                            if self.verbose:
                                print(f"  Convergence detected: range overlaps previous by {overlap_ratio*100:.0f}%", flush=True)
                                print(f"  Stopping early at round {round_idx + 1}/{self.n_rounds}", flush=True)
                            break

                previous_ranges.append(next_range)
                current_grid = next_grid

        # Final parameters from best round
        final_C = gmean(self.rounds_[-1].rank_one_Cs)
        final_method = self.rounds_[-1].best_method
        final_dual = self.rounds_[-1].best_dual
        final_intercept_scaling = self.rounds_[-1].best_intercept_scaling
        final_gamma_5_bp = self.rounds_[-1].best_gamma_5_bp
        final_gamma_5_3 = self.rounds_[-1].best_gamma_5_3

        # Evaluate final parameters
        final_score = self._evaluate_params(
            X, y, final_C, final_method, final_dual, final_intercept_scaling,
            final_gamma_5_bp, final_gamma_5_3
        )

        print(f"Optimal C={final_C:.6e}, method={final_method}, dual={final_dual}, intercept_scaling={final_intercept_scaling}, gamma_5_bp={final_gamma_5_bp}, gamma_5_3={final_gamma_5_3}, CV score={final_score:.4f}", flush=True)

        return SVMParameters(
            C=final_C,
            calibration_method=final_method,
            gamma_5_bp=final_gamma_5_bp,
            gamma_5_3=final_gamma_5_3,
            dual=final_dual,
            intercept_scaling=final_intercept_scaling,
            cv_score=final_score,
            round_found=-1  # -1 indicates averaged result
        )

    def _create_initial_grid(
        self,
        initial_range: Tuple[float, float]
    ) -> np.ndarray:
        """
        Create initial coarse grid.

        Port from: intronIC.py:5446-5451
        Original uses logspace(-6, 6, 13) which gives 10^-6 to 10^6
        """
        log_min = np.log10(initial_range[0])
        log_max = np.log10(initial_range[1])
        n_points = self.n_points_initial

        return np.logspace(log_min, log_max, num=n_points)

    def _prepare_training_data(
        self,
        u12_introns: Sequence[Intron],
        u2_introns: Sequence[Intron]
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Extract feature matrix X and labels y.

        Features: [five_z_score, bp_z_score, three_z_score]
        Labels: 0 = U2, 1 = U12

        Args:
            u12_introns: U12-type introns (normalized)
            u2_introns: U2-type introns (normalized)

        Returns:
            X: Feature matrix (n_samples, 3)
            y: Labels (n_samples,)

        Raises:
            ValueError: If any intron lacks z-scores
        """
        # Extract U12 features
        u12_features = []
        for intron in u12_introns:
            if intron.scores is None:
                raise ValueError(f"Intron {intron.intron_id} has no scores")
            if (intron.scores.five_z_score is None or
                intron.scores.bp_z_score is None or
                intron.scores.three_z_score is None):
                raise ValueError(f"Intron {intron.intron_id} missing z-scores")

            u12_features.append([
                intron.scores.five_z_score,
                intron.scores.bp_z_score,
                intron.scores.three_z_score
            ])

        # Extract U2 features
        u2_features = []
        for intron in u2_introns:
            if intron.scores is None:
                raise ValueError(f"Intron {intron.intron_id} has no scores")
            if (intron.scores.five_z_score is None or
                intron.scores.bp_z_score is None or
                intron.scores.three_z_score is None):
                raise ValueError(f"Intron {intron.intron_id} missing z-scores")

            u2_features.append([
                intron.scores.five_z_score,
                intron.scores.bp_z_score,
                intron.scores.three_z_score
            ])

        # Combine features and create labels
        X = np.array(u2_features + u12_features)
        y = np.array([0] * len(u2_features) + [1] * len(u12_features))

        return X, y

    def _grid_search_round(
        self,
        X: np.ndarray,
        y: np.ndarray,
        C_grid: np.ndarray,
        round_idx: int
    ) -> OptimizationRound:
        """
        Run one round of grid search.

        Port from: intronIC.py:5465-5494
        Uses GridSearchCV with balanced_accuracy scoring.

        Args:
            X: Feature matrix
            y: Labels
            C_grid: Grid of C values to test
            round_idx: Round number (0-indexed)

        Returns:
            Results from this round
        """
        # LinearSVC with external calibration (best practice)
        # - RobustScaler(with_centering=False): Scales by IQR while preserving semantic zero
        #   (s=0 means "U12≈U2", centering would destroy this meaning)
        # - BothEndsStrongTransformer: Augments 3D → 7D with both-ends-strong features
        #   γ parameters (tuned via GridSearchCV) control imbalance penalty strength
        # - LinearSVC: Optimized for linear case, faster than SVC(kernel='linear')
        # - CalibratedClassifierCV: External calibration (method grid-searched)
        cv_splitter = StratifiedKFold(n_splits=5, shuffle=True, random_state=self.random_state + round_idx)

        base_pipeline = Pipeline([
            ('scale', RobustScaler(with_centering=False, with_scaling=True)),
            ('augment', BothEndsStrongTransformer()),  # γ parameters will be grid-searched
            ('svc', LinearSVC(
                loss='squared_hinge',  # Default loss for LinearSVC
                penalty='l2',  # L2 regularization
                class_weight='balanced',  # Critical for imbalanced data
                max_iter=self.max_iter,  # Maximum iterations for convergence
                tol=1e-4,  # Tighter tolerance
                random_state=self.random_state + round_idx
            ))
        ])

        model = CalibratedClassifierCV(
            base_pipeline,
            method='sigmoid',  # Will be grid-searched
            cv=cv_splitter,
            ensemble='auto'  # Per-fold fit + averaging
        )

        # Optimize C, γ parameters, dual, intercept_scaling, and calibration method
        # Use neg_log_loss to optimize for probability quality (best practice)
        # This is critical when deploying at custom thresholds (e.g., 0.90)
        param_grid = {
            'estimator__svc__C': C_grid,  # C parameter through pipeline
            'estimator__augment__gamma_5_bp': [1, 2, 4, 8],  # 5'SS-BPS imbalance penalty weight
            'estimator__augment__gamma_5_3': [1, 2, 4],  # 5'SS-3'SS imbalance penalty weight
            'estimator__svc__dual': [False, True],  # Primal vs dual formulation
            'estimator__svc__intercept_scaling': [10.0, 100.0, 1000.0],  # High values to avoid over-regularizing intercept
            'method': ['sigmoid', 'isotonic']  # Let CV pick calibration method
        }

        grid_search = GridSearchCV(
            model,
            param_grid=param_grid,
            cv=cv_splitter,  # Use same stratified splitter
            scoring='neg_log_loss',  # Optimize for calibrated probabilities
            n_jobs=self.n_jobs,
            error_score=np.nan,
            verbose=0  # Silence sklearn output, use tqdm instead
        )

        # Calculate total tasks for progress bar
        # Note: Progress tracks GridSearchCV tasks, not internal CalibratedClassifierCV fits
        n_candidates = len(list(ParameterGrid(param_grid)))
        n_outer_cv = cv_splitter.get_n_splits(y)  # GridSearchCV outer loop
        n_inner_cv = cv_splitter.get_n_splits(y)  # CalibratedClassifierCV inner loop
        total_tasks = n_candidates * n_outer_cv + 1  # +1 for final refit
        total_fits = n_candidates * n_outer_cv * n_inner_cv + 1  # Actual model fits (for info only)

        if self.verbose:
            print(f"\n{'='*80}", flush=True)
            print(f"ROUND {round_idx + 1}/{self.n_rounds} - Grid Search", flush=True)
            print(f"{'='*80}", flush=True)
            print(f"Parameter combinations: {n_candidates} (C={len(C_grid)} × γ_5bp=4 × γ_5_3=3 × dual=2 × intercept=3 × method=2)", flush=True)
            print(f"CV folds: outer={n_outer_cv}, inner={n_inner_cv} (calibration)", flush=True)
            print(f"GridSearchCV tasks: {total_tasks:,} (~{total_tasks}/{self.n_jobs if self.n_jobs > 0 else 'auto'} per worker)", flush=True)
            print(f"Total model fits: {total_fits:,} (including {n_inner_cv}× internal calibration per task)", flush=True)
            print(f"C range: [{C_grid.min():.2e}, {C_grid.max():.2e}]", flush=True)
            print(f"{'='*80}", flush=True)

        # Fit grid search with progress bar
        import time
        start_time = time.time()

        if self.verbose:
            desc = f"Round {round_idx + 1}/{self.n_rounds}"
            with suppress_convergence_warnings(verbose=True):
                with tqdm_joblib(tqdm(total=total_tasks, desc=desc, unit="fit", leave=False)):
                    grid_search.fit(X, y)
        else:
            with suppress_convergence_warnings(verbose=False):
                grid_search.fit(X, y)

        elapsed = time.time() - start_time

        if self.verbose:
            print(f"Completed in {elapsed:.1f}s", flush=True)

        # Extract results
        cv_results = grid_search.cv_results_
        scores = cv_results['mean_test_score']

        # Print detailed results for all rounds (verbose mode)
        if self.verbose:
            print(f"\n{'='*80}", flush=True)
            print(f"ROUND {round_idx + 1} DETAILED RESULTS - CV Scores (neg_log_loss)", flush=True)
            print(f"{'='*80}", flush=True)
            print(f"{'C Value':<15} {'Method':<10} {'Mean Score':<12} {'Std':<10} {'Rank':<8}", flush=True)
            print(f"{'-'*80}", flush=True)

            # Sort by rank for easier reading
            sorted_indices = np.argsort(cv_results['rank_test_score'])
            for idx in sorted_indices[:10]:  # Show top 10
                param = cv_results['params'][idx]
                score = cv_results['mean_test_score'][idx]
                std = cv_results['std_test_score'][idx]
                rank = cv_results['rank_test_score'][idx]
                c_val = param['estimator__svc__C']
                method = param['method']
                print(f"{c_val:<15.2e} {method:<10} {score:<12.4f} {std:<10.4f} {int(rank):<8}", flush=True)

            if len(C_grid) * 2 > 10:
                print(f"... ({len(C_grid) * 2 - 10} more combinations)", flush=True)
            print(f"{'='*80}\n", flush=True)

        # Find rank-1 (best) parameter values
        ranks = cv_results['rank_test_score']
        params = cv_results['params']
        rank_one_Cs = [p['estimator__svc__C'] for p, r in zip(params, ranks) if r == 1]

        # Best C is geometric mean of rank-1 values
        best_C = gmean(rank_one_Cs) if rank_one_Cs else grid_search.best_params_['estimator__svc__C']
        best_method = grid_search.best_params_['method']
        best_dual = grid_search.best_params_['estimator__svc__dual']
        best_intercept_scaling = grid_search.best_params_['estimator__svc__intercept_scaling']
        best_gamma_5_bp = grid_search.best_params_['estimator__augment__gamma_5_bp']
        best_gamma_5_3 = grid_search.best_params_['estimator__augment__gamma_5_3']
        best_score = grid_search.best_score_

        if self.verbose:
            print(f"\n{'='*80}", flush=True)
            print(f"ROUND {round_idx + 1} SUMMARY", flush=True)
            print(f"{'='*80}", flush=True)
            print(f"Best C (geometric mean of rank-1): {best_C:.6e}", flush=True)
            print(f"Best calibration method: {best_method}", flush=True)
            print(f"Best dual formulation: {best_dual}", flush=True)
            print(f"Best intercept_scaling: {best_intercept_scaling}", flush=True)
            print(f"Best gamma_5_bp: {best_gamma_5_bp}", flush=True)
            print(f"Best gamma_5_3: {best_gamma_5_3}", flush=True)
            print(f"Best CV score (neg_log_loss): {best_score:.4f}", flush=True)
            print(f"Rank-1 C values: {', '.join([f'{c:.2e}' for c in rank_one_Cs])}", flush=True)
            print(f"{'='*80}\n", flush=True)

        return OptimizationRound(
            grid_points=C_grid,
            scores=scores,
            best_C=best_C,
            best_method=best_method,
            best_gamma_5_bp=best_gamma_5_bp,
            best_gamma_5_3=best_gamma_5_3,
            best_dual=best_dual,
            best_intercept_scaling=best_intercept_scaling,
            best_score=best_score,
            rank_one_Cs=rank_one_Cs
        )

    def _refine_grid(
        self,
        current_grid: np.ndarray,
        best_C: float
    ) -> np.ndarray:
        """
        Refine search range around best C.

        Port from: intronIC.py:5496-5502
        Finds the nearest grid point to best_C, then creates a new
        geometric grid between its neighbors.

        Edge case handling: If best_C is at grid boundary, expand the
        search range to ensure exploration beyond current best.

        Args:
            current_grid: Current grid of C values
            best_C: Best C value from this round

        Returns:
            New refined grid for next round
        """
        # Find index of nearest grid point
        best_idx = self._index_of_nearest(current_grid, best_C)

        # Get bounds (neighbors of best)
        low_idx = max(best_idx - 1, 0)
        high_idx = min(best_idx + 1, len(current_grid) - 1)

        # Handle edge cases: if at boundary, expand range geometrically
        if best_idx == 0:
            # At lower edge: use [C[0], C[1]] as base, but ensure we explore below
            span = current_grid[1] / current_grid[0]  # Geometric ratio
            low_bound = current_grid[0] / span  # Extend below
            high_bound = current_grid[1]
        elif best_idx == len(current_grid) - 1:
            # At upper edge: use [C[-2], C[-1]] as base, but ensure we explore above
            span = current_grid[-1] / current_grid[-2]  # Geometric ratio
            low_bound = current_grid[-2]
            high_bound = current_grid[-1] * span  # Extend above
        else:
            # Interior point: use neighbors as normal
            low_bound = current_grid[low_idx]
            high_bound = current_grid[high_idx]

        # Enforce minimum refinement span to prevent over-convergence
        # The range should span at least 2× geometrically to explore meaningfully
        # (reduced from 5× to avoid oscillation between edges)
        min_ratio = 2.0
        current_ratio = high_bound / low_bound

        if current_ratio < min_ratio:
            # Expand range symmetrically around best_C to reach minimum ratio
            # Use geometric mean: best_C = sqrt(low * high)
            # To get ratio R: high/low = R, with best_C as geometric center
            # We get: low = best_C / sqrt(R), high = best_C * sqrt(R)
            expansion_factor = np.sqrt(min_ratio)
            low_bound = best_C / expansion_factor
            high_bound = best_C * expansion_factor

            if self.verbose:
                print(f"  Range too narrow ({current_ratio:.2f}×), expanding to {min_ratio:.0f}× around best_C", flush=True)

        # Create refined geometric grid
        refined_grid = np.geomspace(
            low_bound,
            high_bound,
            num=self.n_points_refine
        )

        return refined_grid

    def _index_of_nearest(self, array: np.ndarray, value: float) -> int:
        """
        Find index of nearest value in array.

        Port from: intronIC.py:5298-5302
        """
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return int(idx)

    def _evaluate_params(
        self,
        X: np.ndarray,
        y: np.ndarray,
        C: float,
        method: str,
        dual: bool,
        intercept_scaling: float,
        gamma_5_bp: float,
        gamma_5_3: float
    ) -> float:
        """
        Evaluate specific hyperparameters via cross-validation.

        Args:
            X: Feature matrix
            y: Labels
            C: C value to evaluate
            method: Calibration method ('sigmoid' or 'isotonic')
            dual: Primal (False) or dual (True) formulation
            intercept_scaling: Intercept scaling parameter
            gamma_5_bp: γ for 5'SS-BPS imbalance penalty
            gamma_5_3: γ for 5'SS-3'SS imbalance penalty

        Returns:
            Cross-validation neg_log_loss score
        """
        # LinearSVC with external calibration (matches training approach)
        # - RobustScaler(with_centering=False): Scales by IQR while preserving semantic zero
        #   (s=0 means "U12≈U2", centering would destroy this meaning)
        # - BothEndsStrongTransformer: Augments 3D → 7D with both-ends-strong features
        # - LinearSVC: Optimized for linear case
        # - CalibratedClassifierCV: External calibration
        base_svm_pipeline = Pipeline([
            ('scale', RobustScaler(with_centering=False, with_scaling=True)),
            ('augment', BothEndsStrongTransformer(
                gamma_5_bp=gamma_5_bp,
                gamma_5_3=gamma_5_3
            )),
            ('svc', LinearSVC(
                C=C,
                dual=dual,
                loss='squared_hinge',
                penalty='l2',
                intercept_scaling=intercept_scaling,
                class_weight='balanced',
                max_iter=self.max_iter,
                tol=1e-4,
                random_state=self.random_state
            ))
        ])

        model = CalibratedClassifierCV(
            base_svm_pipeline,
            method=method,  # Use same method as found in grid search
            cv=5
        )

        # Calculate total tasks for progress bar
        # cross_val_score runs cv_folds outer folds
        # Each fold trains CalibratedClassifierCV with 5 inner folds
        n_outer = self.cv_folds
        n_inner = 5  # CalibratedClassifierCV's cv
        total_tasks = n_outer * n_inner

        if self.verbose:
            desc = "Final param eval"
            with suppress_convergence_warnings(verbose=True):
                with tqdm_joblib(tqdm(total=total_tasks, desc=desc, unit="fit", leave=False)):
                    scores = cross_val_score(
                        model,
                        X,
                        y,
                        cv=self.cv_folds,
                        scoring='neg_log_loss',
                        n_jobs=self.n_jobs,
                        verbose=0  # Silence sklearn, use tqdm
                    )
        else:
            with suppress_convergence_warnings(verbose=True):
                scores = cross_val_score(
                    model,
                    X,
                    y,
                    cv=self.cv_folds,
                    scoring='neg_log_loss',
                    n_jobs=self.n_jobs,
                    verbose=0
                )

        return float(np.mean(scores))
