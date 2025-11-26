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
import sys
import contextlib
import warnings
import joblib
import numpy as np
from sklearn.model_selection import GridSearchCV, cross_val_score, StratifiedKFold, ParameterGrid
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
from sklearn.pipeline import Pipeline
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import make_scorer, balanced_accuracy_score, log_loss, fbeta_score
from scipy.stats import gmean
from tqdm.auto import tqdm

from sklearn.preprocessing import RobustScaler

from intronIC.core.intron import Intron
from intronIC.classification.transformers import BothEndsStrongTransformer

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
    """Optimized SVM hyperparameters for LinearSVC with BothEndsStrong features (NEW ARCHITECTURE - NO CLIPPING)."""

    C: float  # Soft-margin penalty (L2 regularization strength)
    calibration_method: str  # 'sigmoid' or 'isotonic' (isotonic preferred)
    saturate_enabled: bool  # SaturatingTransform enabled (optional log compression)
    include_max: bool  # BothEndsStrong max features (always False in new arch)
    include_pairwise_mins: bool  # BothEndsStrong pairwise mins (always False in new arch)

    # New grid-searched parameters (2025-01-19: robustness improvements)
    penalty: str  # 'l1' or 'l2'
    class_weight_multiplier: float  # Multiplier for balanced weights (0.8, 1.0, 1.2)
    loss: str  # 'hinge' or 'squared_hinge' (for L2); 'squared_hinge' only (for L1)
    gamma_imbalance: float = 1.0  # Gamma scaling for imbalance features (default: 1.0 = no scaling)

    # Fixed parameters (not grid-searched)
    dual: bool = False  # Always False (primal formulation)
    intercept_scaling: float = 1.0  # Fixed to 1.0 (sklearn default; consider 10.0 if L1 intercept issues)

    cv_score: float = 0.0  # Cross-validation balanced_accuracy score
    round_found: int = 0  # Which optimization round found these params (-1 = averaged)


@dataclass(frozen=True, slots=True)
class OptimizationRound:
    """Results from one round of grid search (NEW ARCHITECTURE - NO CLIPPING)."""

    grid_points: np.ndarray  # C values tested
    scores: np.ndarray  # CV scores for each parameter combination
    best_C: float
    best_method: str  # 'isotonic' (or 'sigmoid')
    best_saturate_enabled: bool  # SaturatingTransform enabled
    best_include_max: bool  # BothEndsStrong max features (always False)
    best_include_pairwise_mins: bool  # BothEndsStrong pairwise mins (always False)
    best_score: float  # balanced_accuracy score
    rank_one_Cs: list[float]  # All rank-1 C values

    # New grid-searched parameters (2025-01-19: robustness improvements)
    best_penalty: str = 'l2'  # 'l1' or 'l2'
    best_class_weight_multiplier: float = 1.0  # Multiplier for balanced weights
    best_loss: str = 'squared_hinge'  # 'hinge' or 'squared_hinge'
    best_gamma_imbalance: float = 1.0  # Gamma scaling for imbalance features

    # Fixed parameters (for backward compatibility, not grid-searched)
    best_dual: bool = False  # Always False (primal formulation)
    best_intercept_scaling: float = 1.0  # Fixed to 1.0 (sklearn default)


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
        cv_folds: int = 7,
        random_state: int = 42,
        n_jobs: int = 1,
        verbose: bool = True,
        max_iter: int = 100000,
        scoring_metric: str = 'balanced_accuracy',
        penalty_options: Optional[list] = None,
        loss_options: Optional[list] = None,
        class_weight_multipliers: Optional[list] = None,
        use_multiplier_tiebreaker: bool = True,
        features_list: Optional[list] = None,
        gamma_imbalance_options: Optional[list] = None,
        param_grid_override: Optional[Dict[str, list]] = None,
        progress_tracker: Optional[Any] = None
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
            scoring_metric: Metric for evaluating parameters (default: 'balanced_accuracy').
                          Options: 'balanced_accuracy', 'f_beta_0.5', 'f_beta_0.75'
            penalty_options: Penalty types to search (default: ['l1', 'l2'])
            loss_options: Loss functions to search (default: ['hinge', 'squared_hinge'])
            class_weight_multipliers: Class weight multipliers to search (default: [0.8, 1.0, 1.2])
            use_multiplier_tiebreaker: Prefer 1.0 when multipliers tied (default: True)
            features_list: List of composite feature names to include (default: None = use default 7D)
                          Available features: 'min_5_bp', 'min_5_3', 'min_all',
                                            'neg_absdiff_5_bp', 'neg_absdiff_5_3', 'neg_absdiff_bp_3',
                                            'max_5_bp', 'max_5_3'
            gamma_imbalance_options: List of gamma scaling factors to grid search (default: [1.0])
                                    Scales all neg_absdiff_* features to nudge L2 toward L1 behavior
            param_grid_override: Optional custom parameter grid for testing
                               (if None, uses default full grid). Keys: C is auto-inserted,
                               but can specify: 'estimator__augment__include_max',
                               'estimator__svc__dual', 'estimator__svc__intercept_scaling', 'method'
            progress_tracker: Optional ProgressTracker for global step counting
        """
        self.n_rounds = n_rounds
        self.n_points_initial = n_points_initial
        self.n_points_refine = n_points_refine
        self.cv_folds = cv_folds
        self.random_state = random_state
        self.n_jobs = n_jobs
        self.verbose = verbose
        self.max_iter = max_iter
        self.scoring_metric = scoring_metric
        self.penalty_options = penalty_options if penalty_options is not None else ['l1', 'l2']
        self.loss_options = loss_options if loss_options is not None else ['hinge', 'squared_hinge']
        self.class_weight_multipliers = class_weight_multipliers if class_weight_multipliers is not None else [0.8, 1.0, 1.2]
        self.use_multiplier_tiebreaker = use_multiplier_tiebreaker
        self.features_list = features_list
        self.gamma_imbalance_options = gamma_imbalance_options if gamma_imbalance_options is not None else [1.0]
        self.param_grid_override = param_grid_override
        self.progress_tracker = progress_tracker
        self.rounds_: list[OptimizationRound] = []

    def _create_scorer(self):
        """
        Create sklearn scorer based on configured scoring_metric.

        Returns:
            Scorer object for GridSearchCV

        Supported metrics:
            - 'balanced_accuracy': (TPR + TNR) / 2, treats both classes equally
            - 'f_beta_0.5': F_0.5 score, heavily weights precision (minimizes FPs)
            - 'f_beta_0.75': F_0.75 score, slightly weights precision
        """
        if self.scoring_metric == 'balanced_accuracy':
            return make_scorer(balanced_accuracy_score)
        elif self.scoring_metric.startswith('f_beta_'):
            # Extract beta value from string (e.g., 'f_beta_0.5' -> 0.5)
            beta_str = self.scoring_metric.split('_')[-1]
            try:
                beta = float(beta_str)
            except ValueError:
                raise ValueError(
                    f"Invalid scoring_metric: '{self.scoring_metric}'. "
                    f"Could not parse beta value from '{beta_str}'"
                )
            return make_scorer(fbeta_score, beta=beta)
        else:
            raise ValueError(
                f"Unsupported scoring_metric: '{self.scoring_metric}'. "
                f"Options: 'balanced_accuracy', 'f_beta_0.5', 'f_beta_0.75'"
            )

    def optimize(
        self,
        u12_introns: Sequence[Intron],
        u2_introns: Sequence[Intron],
        eff_C_pos_range: Tuple[float, float] = (1e-2, 1e4),
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
        early_stop = False
        for round_idx in range(self.n_rounds):
            # Note: Round header printed by _grid_search_round() for better organization
            round_result = self._grid_search_round(
                X, y, current_grid, round_idx, eff_C_pos_range
            )
            self.rounds_.append(round_result)

            # Update global progress
            if self.progress_tracker:
                self.progress_tracker.increment(f"Completed optimization round {round_idx + 1}/{self.n_rounds}")

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
                            early_stop = True
                            break

                # If early stop flagged, exit outer loop
                if early_stop:
                    break

                previous_ranges.append(next_range)
                current_grid = next_grid

        # STAGE 1 COMPLETE: Optimal parameters found
        final_C = gmean(self.rounds_[-1].rank_one_Cs)
        final_saturate_enabled = self.rounds_[-1].best_saturate_enabled
        final_include_max = self.rounds_[-1].best_include_max
        final_include_pairwise_mins = self.rounds_[-1].best_include_pairwise_mins
        final_penalty = self.rounds_[-1].best_penalty
        final_loss = self.rounds_[-1].best_loss
        final_class_weight_multiplier = self.rounds_[-1].best_class_weight_multiplier
        final_gamma_imbalance = self.rounds_[-1].best_gamma_imbalance

        if self.verbose:
            print(f"\n{'='*80}", flush=True)
            print(f"Stage 1 Optimization Complete", flush=True)
            print(f"{'='*80}", flush=True)
            print(f"  C = {final_C:.6e}", flush=True)
            print(f"  penalty = {final_penalty}", flush=True)
            print(f"  loss = {final_loss}", flush=True)
            print(f"  class_weight_multiplier = {final_class_weight_multiplier}", flush=True)
            print(f"  gamma_imbalance = {final_gamma_imbalance}", flush=True)
            print(f"{'='*80}\n", flush=True)

        # Calibration method selection using log-loss
        if self.verbose:
            print(f"\n{'='*80}", flush=True)
            print(f"Calibration Method Selection (log-loss comparison)", flush=True)
            print(f"{'='*80}", flush=True)

        # Evaluate sigmoid calibration
        score_sigmoid = self._evaluate_calibration_method(
            X, y, final_C, 'sigmoid',
            final_saturate_enabled, final_include_max, final_include_pairwise_mins,
            final_penalty, final_loss, final_class_weight_multiplier
        )

        # Evaluate isotonic calibration
        score_isotonic = self._evaluate_calibration_method(
            X, y, final_C, 'isotonic',
            final_saturate_enabled, final_include_max, final_include_pairwise_mins,
            final_penalty, final_loss, final_class_weight_multiplier
        )

        # Select winner (lower log-loss = better)
        if score_sigmoid < score_isotonic:
            final_method = 'sigmoid'
            final_score = score_sigmoid
            winner_margin = score_isotonic - score_sigmoid
        else:
            final_method = 'isotonic'
            final_score = score_isotonic
            winner_margin = score_sigmoid - score_isotonic

        if self.verbose:
            print(f"\nCalibration Method Comparison:", flush=True)
            print(f"  Sigmoid:  log-loss = {score_sigmoid:.6f}", flush=True)
            print(f"  Isotonic: log-loss = {score_isotonic:.6f}", flush=True)
            print(f"\n✓ Winner: {final_method} (margin: {winner_margin:.6f})", flush=True)
            print(f"{'='*80}\n", flush=True)

        print(f"Optimal parameters: C={final_C:.6e}, penalty={final_penalty}, loss={final_loss}, "
              f"class_weight_mult={final_class_weight_multiplier}, gamma={final_gamma_imbalance}, "
              f"calibration={final_method}, log-loss={final_score:.6f}", flush=True)

        return SVMParameters(
            C=final_C,
            calibration_method=final_method,
            saturate_enabled=final_saturate_enabled,
            include_max=final_include_max,
            include_pairwise_mins=final_include_pairwise_mins,
            penalty=final_penalty,
            class_weight_multiplier=final_class_weight_multiplier,
            loss=final_loss,
            gamma_imbalance=final_gamma_imbalance,
            dual=False,  # Fixed: primal formulation (n_features << n_samples)
            intercept_scaling=1.0,  # Fixed: sklearn default (works for L1 and L2)
            cv_score=final_score,  # log-loss from Stage 2
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
        round_idx: int,
        eff_C_pos_range: Tuple[float, float] = (1e-2, 1e4)
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
        # TWO-STAGE OPTIMIZATION (Expert approach):
        #
        # STAGE 1 (this function): Optimize C on UNCALIBRATED model
        # - Pipeline: z-scores → LinearSVC (NO scaler, NO calibration)
        # - Input: Pre-scaled z-scores from ScoreNormalizer
        # - Metric: balanced_accuracy (discrimination quality, handles imbalance)
        # - Goal: Find C* that maximizes classification performance
        #
        # STAGE 2 (in optimize()): Choose calibration method at fixed C*
        # - Wrap optimal model in CalibratedClassifierCV
        # - Compare method='sigmoid' vs method='isotonic'
        # - Metric: log-loss (probability calibration quality)
        # - Goal: Find calibration method that gives best-calibrated probabilities
        #
        # This separates discrimination (C) from calibration (method).
        #
        # CORRECTED: Removed RobustScaler (scaling done outside pipeline)
        # See: Expert workflow doc, SCALER_ARCHITECTURE_REVIEW.md
        cv_splitter = StratifiedKFold(n_splits=self.cv_folds, shuffle=True, random_state=self.random_state + round_idx)

        # UNCALIBRATED pipeline for Stage 1
        # Calibration will be added in Stage 2 as post-processing
        # Pipeline: z-scores → augmented features → svc

        # Debug: Log features being used
        if self.verbose and round_idx == 0:
            if self.features_list is not None:
                n_features = 3 + len(self.features_list)  # 3 base + composites
                print(f"  Using {n_features}D feature set: [s5, sBP, s3] + {self.features_list}")
            else:
                # Using transformer default (which is ['absdiff_bp_3'])
                print(f"  Using transformer default: 4D feature set [s5, sBP, s3, absdiff_bp_3]")

        model = Pipeline([
            ('transform', BothEndsStrongTransformer(
                features=self.features_list
            )),
            ('svc', LinearSVC(
                # Parameters set via grid search
                dual=False,  # Primal formulation (recommended when n_features < n_samples)
                max_iter=self.max_iter,
                tol=1e-4,
                random_state=self.random_state + round_idx
            ))
        ])

        # Compute balanced class weights with multipliers
        # Uses configured class_weight_multipliers from config
        n_samples = len(y)
        n_pos = int(np.sum(y == 1))  # U12
        n_neg = int(np.sum(y == 0))  # U2
        w_pos_base = n_samples / (2.0 * n_pos)
        w_neg_base = n_samples / (2.0 * n_neg)

        # STAGE 1 PARAMETER GRID:
        # Optimize C, penalty, loss, and class_weight_multiplier on uncalibrated model
        # No calibration parameters (that's Stage 2)
        #
        # CLASS_WEIGHT_MULTIPLIER: Controls precision-recall tradeoff
        # - Only scales w_pos (U12), keeps w_neg (U2) fixed
        # - This changes the RELATIVE penalty: w_pos/w_neg ratio varies with multiplier
        # - multiplier > 1.0: penalize U12 errors more → higher recall (find more U12s)
        # - multiplier < 1.0: penalize U12 errors less → higher precision (fewer false positives)
        #
        # IMPORTANT: Use same C_grid for all multipliers (no per-multiplier C scaling)
        # - Let effective penalties naturally vary with multiplier
        # - This allows genuine exploration of precision-recall tradeoff
        # - C bounds computed for multiplier=1.0 (balanced weights)

        if self.param_grid_override is not None:
            # Use custom parameter grid for fast testing
            param_grid = {'svc__C': C_grid}  # Direct access (no 'estimator__' prefix)
            param_grid.update(self.param_grid_override)
        else:
            # Full parameter grid for production
            # Uses configured penalty_options, loss_options, class_weight_multipliers
            #
            # SKLEARN CONSTRAINTS (dual=False):
            # - penalty='l1': only supports loss='squared_hinge'
            # - penalty='l2': only supports loss='squared_hinge' (NOT hinge when dual=False)
            #
            # Since we use dual=False, we can only use loss='squared_hinge'
            # Filter loss_options to only include 'squared_hinge'
            valid_losses = [loss for loss in self.loss_options if loss == 'squared_hinge']
            if not valid_losses:
                # Fallback to squared_hinge if user removed it from config
                valid_losses = ['squared_hinge']

            # Create list of parameter grids (one per (multiplier, gamma) combination)
            # GridSearchCV accepts a list of dicts, and will search the union of them
            # Grid over both class_weight_multipliers and gamma_imbalance_options
            param_grid_list = []
            for multiplier in self.class_weight_multipliers:
                # Create class_weight dict for this multiplier
                # CRITICAL: Only scale w_pos (U12), keep w_neg (U2) fixed
                # This creates a genuine precision-recall tradeoff
                class_weight = {
                    0: w_neg_base,               # U2 - FIXED
                    1: w_pos_base * multiplier   # U12 - SCALED by multiplier
                }

                # Create a parameter grid for this multiplier
                # Use SAME C_grid for all multipliers (no per-multiplier scaling)
                # Include gamma_imbalance in grid (scales neg_absdiff_* features)
                param_grid_list.append({
                    'svc__C': list(C_grid),  # Same C values for all multipliers
                    'svc__class_weight': [class_weight],  # Single class_weight dict (in a list)
                    'svc__penalty': self.penalty_options,
                    'svc__loss': valid_losses,
                    'transform__gamma_imbalance': self.gamma_imbalance_options  # Grid over gamma
                })

            param_grid = param_grid_list
            # Grid size: len(C_grid) × len(penalty_options) × 1 loss × len(class_weight_multipliers) × len(gamma_imbalance_options)

        # Create scorer based on configured metric
        # Options: balanced_accuracy, f_beta_0.5, f_beta_0.75
        scorer = self._create_scorer()

        grid_search = GridSearchCV(
            model,
            param_grid=param_grid,
            cv=cv_splitter,  # Use same stratified splitter
            scoring=scorer,  # Configured metric (default: balanced_accuracy)
            n_jobs=self.n_jobs,
            error_score=np.nan,
            verbose=0  # Silence sklearn output, use tqdm instead
        )

        # Calculate total tasks for progress bar
        # Stage 1: Uncalibrated model, direct CV evaluation
        # param_grid can be either a dict or list of dicts
        n_candidates = len(list(ParameterGrid(param_grid)))
        n_cv_folds = cv_splitter.get_n_splits(y)  # Cross-validation folds
        total_tasks = n_candidates * n_cv_folds + 1  # +1 for final refit
        total_fits = total_tasks  # Each task = one model fit (no calibration wrapper)

        if self.verbose:
            print(f"\n{'='*80}", flush=True)
            global_step_str = f" {self.progress_tracker.format_step()}" if self.progress_tracker else ""
            print(f"ROUND {round_idx + 1}/{self.n_rounds} - Grid Search (C Optimization){global_step_str}", flush=True)
            print(f"{'='*80}", flush=True)
            print(f"Parameter combinations: {n_candidates} (C={len(C_grid)})", flush=True)
            print(f"Model: UNCALIBRATED LinearSVC (calibration selected later)", flush=True)
            print(f"CV folds: {n_cv_folds}", flush=True)
            print(f"GridSearchCV tasks: {total_tasks:,} (~{total_tasks}/{self.n_jobs if self.n_jobs > 0 else 'auto'} per worker)", flush=True)
            print(f"Total model fits: {total_fits:,}", flush=True)
            print(f"C range: [{C_grid.min():.2e}, {C_grid.max():.2e}]", flush=True)
            print(f"Metric: {self.scoring_metric}", flush=True)
            if self.n_jobs not in [0, 1, -1] and abs(self.n_jobs) > 4:
                print(f"Note: With high parallelism (n_jobs={self.n_jobs}), progress bar may update", flush=True)
                print(f"      in large increments as worker batches complete.", flush=True)
            print(f"{'='*80}", flush=True)

        # Fit grid search with progress bar
        import time
        import threading
        start_time = time.time()

        # Progress monitoring (time-based fallback for high parallelism)
        # Note: With high n_jobs (e.g., 10+), joblib batches work and tqdm only
        # updates when batches complete, which can appear as long pauses.
        # Add periodic time updates as feedback that work is progressing.
        stop_monitor = threading.Event()
        def monitor_progress():
            """Print time-based progress updates every 30s as fallback."""
            interval = 30  # seconds between updates
            while not stop_monitor.wait(interval):
                elapsed = time.time() - start_time
                print(f"  [Still running... {elapsed:.0f}s elapsed]", flush=True)

        # Start background monitor for long-running tasks
        monitor_thread = None
        if self.verbose and total_tasks > 50:  # Only for substantial grid searches
            monitor_thread = threading.Thread(target=monitor_progress, daemon=True)
            monitor_thread.start()

        try:
            if self.verbose:
                desc = f"Round {round_idx + 1}/{self.n_rounds}"
                with suppress_convergence_warnings(verbose=True):
                    with tqdm_joblib(tqdm(
                        total=total_tasks,
                        desc=desc,
                        unit="fit",
                        leave=True,  # Keep bar visible after completion
                        mininterval=0.05,  # Faster refresh (50ms) for more responsive updates
                        miniters=1,  # Update on every increment
                        ncols=100,  # Fixed width for consistent display
                        file=sys.stdout,  # Force output to stdout (not captured by logging)
                        dynamic_ncols=False  # Disable dynamic width (can interfere with updates)
                    )):
                        grid_search.fit(X, y)
            else:
                with suppress_convergence_warnings(verbose=False):
                    grid_search.fit(X, y)
        finally:
            # Stop background monitor
            if monitor_thread is not None:
                stop_monitor.set()
                monitor_thread.join(timeout=1)

        elapsed = time.time() - start_time

        if self.verbose:
            print(f"Completed in {elapsed:.1f}s", flush=True)

        # Extract results
        cv_results = grid_search.cv_results_
        scores = cv_results['mean_test_score']

        # Print detailed results for all rounds (verbose mode)
        if self.verbose:
            print(f"\n{'='*80}", flush=True)
            print(f"ROUND {round_idx + 1} DETAILED RESULTS - CV Scores ({self.scoring_metric})", flush=True)
            print(f"{'='*80}", flush=True)
            print(f"{'C Value':<15} {'Mean Score':<12} {'Std':<10} {'Rank':<8}", flush=True)
            print(f"{'-'*80}", flush=True)

            # Sort by rank for easier reading
            sorted_indices = np.argsort(cv_results['rank_test_score'])
            for idx in sorted_indices[:10]:  # Show top 10
                param = cv_results['params'][idx]
                score = cv_results['mean_test_score'][idx]
                std = cv_results['std_test_score'][idx]
                rank = cv_results['rank_test_score'][idx]
                c_val = param['svc__C']  # Direct pipeline access (no 'estimator__' prefix)
                print(f"{c_val:<15.2e} {score:<12.4f} {std:<10.4f} {int(rank):<8}", flush=True)

            if len(C_grid) > 10:
                print(f"... ({len(C_grid) - 10} more combinations)", flush=True)
            print(f"{'='*80}\n", flush=True)

        # Find rank-1 (best) parameter values
        ranks = cv_results['rank_test_score']
        params = cv_results['params']
        rank_one_Cs = [p['svc__C'] for p, r in zip(params, ranks) if r == 1]  # Direct pipeline access

        # Best C is geometric mean of rank-1 values
        best_C = gmean(rank_one_Cs) if rank_one_Cs else grid_search.best_params_['svc__C']

        # STAGE 1: Extract hyperparameters
        # Calibration method will be selected in Stage 2 (done in optimize())
        best_method = None  # Will be determined in Stage 2 (not fixed to isotonic anymore)
        # Phase 1: No augmentation/saturation (set to False for backward compatibility)
        best_saturate_enabled = False
        best_include_max = False
        best_include_pairwise_mins = False
        best_score = grid_search.best_score_

        # Extract new parameters from best_params_
        best_penalty = grid_search.best_params_.get('svc__penalty', 'l2')
        best_loss = grid_search.best_params_.get('svc__loss', 'squared_hinge')
        best_gamma_imbalance = grid_search.best_params_.get('transform__gamma_imbalance', 1.0)

        # Map class_weight dict back to multiplier
        # Extract from best_params_ and compute which multiplier was used
        best_class_weight_dict = grid_search.best_params_.get('svc__class_weight')

        # Compute multiplier from class_weight values
        # class_weight[1] = w_pos_base * multiplier
        # So: multiplier = class_weight[1] / w_pos_base
        if best_class_weight_dict and 1 in best_class_weight_dict:
            best_class_weight_multiplier = best_class_weight_dict[1] / w_pos_base
        else:
            best_class_weight_multiplier = 1.0  # Fallback

        # DIAGNOSTIC: Verify all class_weight_multipliers were evaluated
        # Extract all tested multipliers from cv_results_ and track best score per multiplier
        tested_multipliers = set()
        multiplier_scores = {}  # {multiplier: best_score}

        for i, params in enumerate(grid_search.cv_results_['params']):
            class_weight_dict = params.get('svc__class_weight')
            if class_weight_dict and 1 in class_weight_dict:
                multiplier = class_weight_dict[1] / w_pos_base
                mult_rounded = round(multiplier, 4)
                tested_multipliers.add(mult_rounded)

                # Track best score for each multiplier
                score = grid_search.cv_results_['mean_test_score'][i]
                if mult_rounded not in multiplier_scores or score > multiplier_scores[mult_rounded]:
                    multiplier_scores[mult_rounded] = score

        # TIE-BREAKING: When multiple multipliers have the same best score, prefer 1.0
        # (if use_multiplier_tiebreaker is enabled)
        best_score_value = multiplier_scores.get(round(best_class_weight_multiplier, 4), 0.0)
        score_tolerance = 1e-6  # Consider scores within this range as tied

        tied_multipliers = [
            mult for mult, score in multiplier_scores.items()
            if abs(score - best_score_value) < score_tolerance
        ]

        tie_broken = False
        original_selection = best_class_weight_multiplier

        if self.use_multiplier_tiebreaker and len(tied_multipliers) > 1:
            # Multiple multipliers tied for best score
            if 1.0 in tied_multipliers:
                # Prefer 1.0 (standard balanced weighting) when tied
                best_class_weight_multiplier = 1.0
                tie_broken = (original_selection != 1.0)
            else:
                # If 1.0 not in tie, prefer multiplier closest to 1.0
                best_class_weight_multiplier = min(tied_multipliers, key=lambda m: abs(m - 1.0))
                tie_broken = (best_class_weight_multiplier != original_selection)

        if self.verbose:
            print(f"\nDIAGNOSTIC: Tested {len(tested_multipliers)} unique class_weight_multipliers:", flush=True)
            print(f"  Configured: {sorted(set(round(m, 4) for m in self.class_weight_multipliers))}", flush=True)
            print(f"  Actually tested: {sorted(tested_multipliers)}", flush=True)
            print(f"\n  Best score for each multiplier:", flush=True)

            for mult in sorted(multiplier_scores.keys()):
                score = multiplier_scores[mult]
                is_selected = abs(mult - round(best_class_weight_multiplier, 4)) < 0.001
                is_tied = abs(score - best_score_value) < score_tolerance

                marker = ""
                if is_selected:
                    marker = " ← SELECTED"
                elif is_tied:
                    marker = " (tied)"

                print(f"    {mult:5.2f}: {score:.6f}{marker}", flush=True)

            if len(tied_multipliers) > 1:
                print(f"\n  ℹ Note: {len(tied_multipliers)} multipliers tied at {best_score_value:.6f}", flush=True)
                if tie_broken:
                    print(f"  Tie-breaking enabled: Changed from {original_selection:.2f} to {best_class_weight_multiplier:.2f} (prefer 1.0)", flush=True)
                else:
                    if self.use_multiplier_tiebreaker:
                        print(f"  Tie-breaking enabled: Already selected {best_class_weight_multiplier:.2f} (optimal choice)", flush=True)
                    else:
                        print(f"  Tie-breaking disabled: Using GridSearchCV default {best_class_weight_multiplier:.2f}", flush=True)

        # VALIDATION: Check if selected parameters create pathological effective C
        # Effective C_pos = C × w_pos × multiplier
        eff_C_pos_selected = best_C * w_pos_base * best_class_weight_multiplier

        # Compute expected range accounting for class_weight_multiplier variation
        # The eff_C_pos_range is for multiplier=1.0 (balanced weights)
        # With other multipliers, we expect natural variation:
        #   mult < 1.0: lower eff_C_pos (prioritize precision)
        #   mult > 1.0: higher eff_C_pos (prioritize recall)
        #
        # Expand range by multiplier bounds to account for intentional variation
        mult_min = min(self.class_weight_multipliers)
        mult_max = max(self.class_weight_multipliers)

        # Expected range adjusted for multiplier variation, with 2x tolerance
        eff_C_min_expected = eff_C_pos_range[0] * mult_min / 2.0
        eff_C_max_expected = eff_C_pos_range[1] * mult_max * 2.0

        if not (eff_C_min_expected <= eff_C_pos_selected <= eff_C_max_expected):
            import warnings
            warnings.warn(
                f"Grid search selected C={best_C:.2e} with multiplier={best_class_weight_multiplier:.2f}, "
                f"giving eff_C_pos={eff_C_pos_selected:.2e}. "
                f"This is outside expected range [{eff_C_min_expected:.2e}, {eff_C_max_expected:.2e}] "
                f"(accounting for multiplier range [{mult_min:.2f}, {mult_max:.2f}]). "
                f"This may indicate optimizer instability or pathological parameter selection.",
                RuntimeWarning
            )

        if self.verbose:
            print(f"\n{'='*80}", flush=True)
            print(f"ROUND {round_idx + 1} SUMMARY", flush=True)
            print(f"{'='*80}", flush=True)
            print(f"Best C (geometric mean of rank-1): {best_C:.6e}", flush=True)
            print(f"Best penalty: {best_penalty}", flush=True)
            print(f"Best loss: {best_loss}", flush=True)
            print(f"Best class_weight_multiplier: {best_class_weight_multiplier:.2f}", flush=True)
            print(f"Effective C_pos: {eff_C_pos_selected:.2e} (target: [{eff_C_pos_range[0]:.0e}, {eff_C_pos_range[1]:.0e}])", flush=True)
            print(f"Model: UNCALIBRATED LinearSVC", flush=True)
            print(f"Best CV score ({self.scoring_metric}): {best_score:.4f}", flush=True)
            print(f"Rank-1 C values: {', '.join([f'{c:.2e}' for c in rank_one_Cs])}", flush=True)
            print(f"Note: Calibration method (sigmoid vs isotonic) will be selected after C optimization completes", flush=True)
            print(f"{'='*80}\n", flush=True)

        return OptimizationRound(
            grid_points=C_grid,
            scores=scores,
            best_C=best_C,
            best_method=best_method,  # Will be determined in Stage 2
            best_saturate_enabled=best_saturate_enabled,
            best_include_max=best_include_max,
            best_include_pairwise_mins=best_include_pairwise_mins,
            best_score=best_score,
            rank_one_Cs=rank_one_Cs,
            best_penalty=best_penalty,
            best_class_weight_multiplier=best_class_weight_multiplier,
            best_loss=best_loss,
            best_gamma_imbalance=best_gamma_imbalance
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
                print(f"\n[Next round preparation] Range too narrow ({current_ratio:.2f}×), expanding to {min_ratio:.0f}× around best_C", flush=True)

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

    def _evaluate_calibration_method(
        self,
        X: np.ndarray,
        y: np.ndarray,
        C: float,
        method: str,
        saturate_enabled: bool,
        include_max: bool,
        include_pairwise_mins: bool,
        penalty: str = 'l2',
        loss: str = 'squared_hinge',
        class_weight_multiplier: float = 1.0
    ) -> float:
        """
        Evaluate calibration method using log-loss (STAGE 2).

        Args:
            X: Feature matrix
            y: Labels
            C: C value (from Stage 1)
            method: Calibration method ('sigmoid' or 'isotonic')
            saturate_enabled: Ignored in Phase 1
            include_max: Ignored in Phase 1
            include_pairwise_mins: Ignored in Phase 1
            penalty: Penalty type ('l1' or 'l2')
            loss: Loss function ('hinge' or 'squared_hinge')
            class_weight_multiplier: Multiplier for balanced class weights

        Returns:
            Cross-validation log-loss (lower = better calibration)
        """
        # Compute balanced class weights with multiplier
        n_samples = len(y)
        n_pos = int(np.sum(y == 1))  # U12
        n_neg = int(np.sum(y == 0))  # U2
        w_pos = (n_samples / (2.0 * n_pos)) * class_weight_multiplier
        w_neg = (n_samples / (2.0 * n_neg)) * class_weight_multiplier
        class_weight = {0: w_neg, 1: w_pos}

        # CORRECTED: Removed RobustScaler (scaling done outside pipeline)
        # Pipeline: z-scores → augmented features → svc
        base_svm_pipeline = Pipeline([
            ('transform', BothEndsStrongTransformer(
                features=self.features_list
            )),
            ('svc', LinearSVC(
                C=C,
                penalty=penalty,
                dual=False,
                loss=loss,
                class_weight=class_weight,
                max_iter=self.max_iter,
                tol=1e-4,
                random_state=self.random_state
            ))
        ])

        model = CalibratedClassifierCV(
            base_svm_pipeline,
            method=method,  # The method we're evaluating
            cv=self.cv_folds,
            ensemble='auto'
        )

        # Use sklearn's built-in neg_log_loss scorer
        # (automatically handles probability predictions)
        # Returns negative log-loss, so we negate to get positive log-loss
        with suppress_convergence_warnings(verbose=False):
            scores = cross_val_score(
                model,
                X,
                y,
                cv=self.cv_folds,
                scoring='neg_log_loss',  # Built-in scorer for log-loss
                n_jobs=self.n_jobs,
                verbose=0
            )

        # Return mean log-loss (negate the negative to get positive)
        return float(-np.mean(scores))  # Return positive log-loss

    def _evaluate_params(
        self,
        X: np.ndarray,
        y: np.ndarray,
        C: float,
        method: str,
        saturate_enabled: bool,
        include_max: bool,
        include_pairwise_mins: bool
    ) -> float:
        """
        Evaluate specific hyperparameters via cross-validation (PHASE 1: CENTERED SCALING).

        Args:
            X: Feature matrix
            y: Labels
            C: C value to evaluate
            method: Calibration method ('sigmoid' or 'isotonic')
            saturate_enabled: Ignored in Phase 1 (backward compatibility)
            include_max: Ignored in Phase 1 (backward compatibility)
            include_pairwise_mins: Ignored in Phase 1 (backward compatibility)

        Returns:
            Cross-validation balanced_accuracy score
        """
        # CORRECTED: Removed RobustScaler (scaling done outside pipeline)
        # See: Expert workflow doc, SCALER_ARCHITECTURE_REVIEW.md
        #
        # Pipeline: z-scores → augmented features → svc → calibration
        # - Data arrives as z-scores (already scaled by ScoreNormalizer)
        # - BothEndsStrongTransformer: Augmented features from z-scores
        # - LinearSVC: L2-regularized linear classifier
        # - CalibratedClassifierCV: External calibration (isotonic or sigmoid)
        #
        # Scoring: balanced_accuracy (designed for imbalanced data)
        base_svm_pipeline = Pipeline([
            ('transform', BothEndsStrongTransformer(
                features=self.features_list
            )),
            ('svc', LinearSVC(
                C=C,
                penalty='l2',
                dual=False,
                loss='squared_hinge',
                class_weight='balanced',
                max_iter=self.max_iter,
                tol=1e-4,
                random_state=self.random_state
            ))
        ])

        model = CalibratedClassifierCV(
            base_svm_pipeline,
            method=method,  # Use same method as found in grid search
            cv=self.cv_folds,
            ensemble='auto'  # Per-fold fit + averaging
        )

        # balanced_accuracy scorer: (TPR + TNR) / 2
        # Designed for imbalanced data
        scorer = make_scorer(balanced_accuracy_score)

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
                        scoring=scorer,  # balanced_accuracy
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
                    scoring=scorer,  # balanced_accuracy
                    n_jobs=self.n_jobs,
                    verbose=0
                )

        return float(np.mean(scores))
