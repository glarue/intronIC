"""
SVM hyperparameter optimization using geometric grid search with refinement.

This module implements the optimization algorithm from intronIC.py:5431-5528,
which uses 5 rounds of progressively refined grid search to find optimal
SVM hyperparameters (specifically the soft-margin penalty C).

Key algorithm:
- Round 1: Coarse grid search (10^-6 to 10^6)
- Rounds 2-5: Refine around best C values from previous round
- Final: Geometric mean of all rank-1 C values from final round

Port from: intronIC.py:5431-5528 (optimize_svm)
Related: intronIC.py:5290-5322 (helper functions)
"""

from dataclasses import dataclass
from typing import Sequence, Tuple, Optional
import os
import numpy as np
from sklearn.model_selection import GridSearchCV, cross_val_score, train_test_split, StratifiedKFold
from sklearn.svm import SVC
from sklearn.calibration import CalibratedClassifierCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import RobustScaler
from scipy.stats import gmean

from core.intron import Intron


@dataclass(frozen=True, slots=True)
class SVMParameters:
    """Optimized SVM hyperparameters."""

    C: float  # Soft-margin penalty
    calibration_method: str  # 'sigmoid' or 'isotonic'
    cv_score: float  # Cross-validation neg_log_loss
    round_found: int  # Which optimization round found these params (-1 = averaged)


@dataclass(frozen=True, slots=True)
class OptimizationRound:
    """Results from one round of grid search."""

    grid_points: np.ndarray  # C values tested
    scores: np.ndarray  # CV scores for each C
    best_C: float
    best_method: str  # 'sigmoid' or 'isotonic'
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
        verbose: bool = True
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
        """
        self.n_rounds = n_rounds
        self.n_points_initial = n_points_initial
        self.n_points_refine = n_points_refine
        self.cv_folds = cv_folds
        self.random_state = random_state
        self.n_jobs = n_jobs
        self.verbose = verbose
        self.rounds_: list[OptimizationRound] = []

    def optimize(
        self,
        u12_introns: Sequence[Intron],
        u2_introns: Sequence[Intron],
        initial_range: Tuple[float, float] = (1e-6, 1e6)
    ) -> SVMParameters:
        """
        Find optimal C parameter via geometric grid search.

        Args:
            u12_introns: Training U12-type introns (with z-scores)
            u2_introns: Training U2-type introns (with z-scores)
            initial_range: Initial search range for C (min, max)

        Returns:
            Optimized parameters with best C value

        Raises:
            ValueError: If introns lack z-scores
        """
        # Extract features and labels
        X, y = self._prepare_training_data(u12_introns, u2_introns)

        # Initialize search range
        current_grid = self._create_initial_grid(initial_range)

        # Run geometric refinement
        for round_idx in range(self.n_rounds):
            print(f"Optimization round {round_idx + 1}/{self.n_rounds}...", flush=True)

            round_result = self._grid_search_round(
                X, y, current_grid, round_idx
            )
            self.rounds_.append(round_result)

            # Prepare next round's grid (refine around best)
            if round_idx < self.n_rounds - 1:
                current_grid = self._refine_grid(current_grid, round_result.best_C)

        # Final parameter: geometric mean of final round's rank-1 values
        final_C = gmean(self.rounds_[-1].rank_one_Cs)
        final_method = self.rounds_[-1].best_method

        # Evaluate final C with final method
        final_score = self._evaluate_C(X, y, final_C, final_method)

        print(f"Optimal C={final_C:.6e}, method={final_method}, CV score={final_score:.4f}", flush=True)

        return SVMParameters(
            C=final_C,
            calibration_method=final_method,
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
        # SVC with external calibration (best practice)
        # - RobustScaler: Matches legacy median/IQR scaling
        # - SVC(probability=False): Linear margin without slow internal Platt
        # - CalibratedClassifierCV: External calibration (method grid-searched)
        cv_splitter = StratifiedKFold(n_splits=5, shuffle=True, random_state=self.random_state + round_idx)

        base_pipeline = Pipeline([
            ('scale', RobustScaler(with_centering=True, with_scaling=True)),
            ('svc', SVC(
                kernel='linear',
                probability=False,  # We calibrate externally
                class_weight='balanced',
                cache_size=1000,
                random_state=self.random_state + round_idx
            ))
        ])

        model = CalibratedClassifierCV(
            base_pipeline,
            method='sigmoid',  # Will be grid-searched
            cv=cv_splitter,
            ensemble='auto'  # Per-fold fit + averaging
        )

        # Use 80% train split for GridSearchCV (matches original)
        X_train, X_test, y_train, y_test = train_test_split(
            X, y,
            train_size=0.80,
            stratify=y,
            random_state=self.random_state + round_idx
        )

        # Optimize C and calibration method
        # Use neg_log_loss to optimize for probability quality (best practice)
        # This is critical when deploying at custom thresholds (e.g., 0.90)
        param_grid = {
            'estimator__svc__C': C_grid,  # C parameter through pipeline
            'method': ['sigmoid', 'isotonic']  # Let CV pick calibration method
        }

        grid_search = GridSearchCV(
            model,
            param_grid=param_grid,
            cv=cv_splitter,  # Use same stratified splitter
            scoring='neg_log_loss',  # Optimize for calibrated probabilities
            n_jobs=self.n_jobs,
            error_score=np.nan,
            verbose=10 if self.verbose else 0  # Max verbosity to force output
        )

        if self.verbose:
            print(f"\n{'='*80}", flush=True)
            print(f"ROUND {round_idx + 1}/{self.n_rounds} - Grid Search", flush=True)
            print(f"{'='*80}", flush=True)
            print(f"Testing {len(C_grid)} C values × 2 calibration methods = {len(C_grid) * 2} total fits", flush=True)
            print(f"Each fit trains {self.cv_folds} CV folds + final model", flush=True)
            print(f"C range: [{C_grid.min():.2e}, {C_grid.max():.2e}]", flush=True)
            print(f"CV folds: {self.cv_folds}, Jobs: {self.n_jobs}", flush=True)
            print(f"{'='*80}", flush=True)
            print(f"Starting grid search... (this may take a few minutes)", flush=True)
            import sys
            sys.stdout.flush()

        # Fit on 80% subset for faster optimization
        import time
        start_time = time.time()
        grid_search.fit(X_train, y_train)
        elapsed = time.time() - start_time

        if self.verbose:
            print(f"\nGrid search completed in {elapsed:.1f} seconds", flush=True)

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

        # Find rank-1 (best) C values and method
        ranks = cv_results['rank_test_score']
        params = cv_results['params']
        rank_one_Cs = [p['estimator__svc__C'] for p, r in zip(params, ranks) if r == 1]

        # Best C is geometric mean of rank-1 values
        best_C = gmean(rank_one_Cs) if rank_one_Cs else grid_search.best_params_['estimator__svc__C']
        best_method = grid_search.best_params_['method']
        best_score = grid_search.best_score_

        if self.verbose:
            print(f"\n{'='*80}", flush=True)
            print(f"ROUND {round_idx + 1} SUMMARY", flush=True)
            print(f"{'='*80}", flush=True)
            print(f"Best C (geometric mean of rank-1): {best_C:.6e}", flush=True)
            print(f"Best calibration method: {best_method}", flush=True)
            print(f"Best CV score (neg_log_loss): {best_score:.4f}", flush=True)
            print(f"Rank-1 C values: {', '.join([f'{c:.2e}' for c in rank_one_Cs])}", flush=True)
            print(f"{'='*80}\n", flush=True)

        return OptimizationRound(
            grid_points=C_grid,
            scores=scores,
            best_C=best_C,
            best_method=best_method,
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

        low_bound = current_grid[low_idx]
        high_bound = current_grid[high_idx]

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

    def _evaluate_C(self, X: np.ndarray, y: np.ndarray, C: float, method: str) -> float:
        """
        Evaluate a specific C value and calibration method via cross-validation.

        Args:
            X: Feature matrix
            y: Labels
            C: C value to evaluate
            method: Calibration method ('sigmoid' or 'isotonic')

        Returns:
            Cross-validation neg_log_loss score
        """
        # SVC with external calibration (matches training approach)
        # - RobustScaler: Matches legacy median/IQR scaling
        # - SVC(probability=False): Linear margin without slow internal Platt
        # - CalibratedClassifierCV: External calibration
        base_svm_pipeline = Pipeline([
            ('scale', RobustScaler(with_centering=True, with_scaling=True)),
            ('svc', SVC(
                C=C,
                kernel='linear',
                probability=False,  # We calibrate externally
                class_weight='balanced',
                cache_size=1000,
                random_state=self.random_state
            ))
        ])

        model = CalibratedClassifierCV(
            base_svm_pipeline,
            method=method,  # Use same method as found in grid search
            cv=5
        )

        scores = cross_val_score(
            model,
            X,
            y,
            cv=self.cv_folds,
            scoring='neg_log_loss',  # Evaluate probability quality
            n_jobs=self.n_jobs,  # Parallelize CV folds
            verbose=2  # Show real-time progress
        )

        return float(np.mean(scores))
