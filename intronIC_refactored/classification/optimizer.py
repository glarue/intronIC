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
from sklearn.model_selection import GridSearchCV, cross_val_score, train_test_split
from sklearn.svm import SVC
from scipy.stats import gmean

from core.intron import Intron


@dataclass(frozen=True, slots=True)
class SVMParameters:
    """Optimized SVM hyperparameters."""

    C: float  # Soft-margin penalty
    cv_score: float  # Cross-validation balanced accuracy
    round_found: int  # Which optimization round found these params (-1 = averaged)


@dataclass(frozen=True, slots=True)
class OptimizationRound:
    """Results from one round of grid search."""

    grid_points: np.ndarray  # C values tested
    scores: np.ndarray  # CV scores for each C
    best_C: float
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
        n_rounds: int = 5,
        n_points_initial: int = 13,
        n_points_refine: int = 100,
        cv_folds: int = 5,
        random_state: int = 42,
        n_jobs: int = 1
    ):
        """
        Initialize optimizer.

        Args:
            n_rounds: Number of refinement rounds (default: 5)
            n_points_initial: Grid points for round 1 (default: 13)
            n_points_refine: Grid points for rounds 2-5 (default: 100)
            cv_folds: Cross-validation folds (default: 5)
            random_state: Random seed for reproducibility
            n_jobs: Number of parallel jobs for GridSearchCV (default: 1)
        """
        self.n_rounds = n_rounds
        self.n_points_initial = n_points_initial
        self.n_points_refine = n_points_refine
        self.cv_folds = cv_folds
        self.random_state = random_state
        self.n_jobs = n_jobs
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
            print(f"Optimization round {round_idx + 1}/{self.n_rounds}...")

            round_result = self._grid_search_round(
                X, y, current_grid, round_idx
            )
            self.rounds_.append(round_result)

            # Prepare next round's grid (refine around best)
            if round_idx < self.n_rounds - 1:
                current_grid = self._refine_grid(current_grid, round_result.best_C)

        # Final parameter: geometric mean of final round's rank-1 values
        final_C = gmean(self.rounds_[-1].rank_one_Cs)

        # Evaluate final C
        final_score = self._evaluate_C(X, y, final_C)

        print(f"Optimal C={final_C:.6e}, CV score={final_score:.4f}")

        return SVMParameters(
            C=final_C,
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
        # Create base model - MUST have class_weight='balanced' for imbalanced classes!
        base_svm = SVC(
            kernel='linear',
            class_weight='balanced',  # CRITICAL: 990 U2 vs 97 U12 = 10:1 imbalance
            cache_size=1000,  # MB - critical for performance with large datasets
            random_state=self.random_state + round_idx
        )

        # Use 80% train split for GridSearchCV (matches original)
        # Reduces dataset size to speed up optimization
        # Testing showed this helps with larger datasets (1000+ samples)
        X_train, X_test, y_train, y_test = train_test_split(
            X, y,
            train_size=0.80,
            stratify=y,
            random_state=self.random_state + round_idx
        )

        grid_search = GridSearchCV(
            base_svm,
            param_grid={'C': C_grid},
            cv=self.cv_folds,
            scoring='balanced_accuracy',
            n_jobs=self.n_jobs,  # Parallelize CV folds
            error_score=np.nan,
            verbose=2  # Show real-time progress
        )

        # Fit on 80% subset for faster optimization
        grid_search.fit(X_train, y_train)

        # Extract results
        cv_results = grid_search.cv_results_
        scores = cv_results['mean_test_score']

        # Find rank-1 (best) C values
        # Following original logic in rank_ones() and rank1_param_avg()
        ranks = cv_results['rank_test_score']
        params = cv_results['params']
        rank_one_Cs = [p['C'] for p, r in zip(params, ranks) if r == 1]

        # Best C is geometric mean of rank-1 values
        best_C = gmean(rank_one_Cs) if rank_one_Cs else grid_search.best_params_['C']
        best_score = grid_search.best_score_

        return OptimizationRound(
            grid_points=C_grid,
            scores=scores,
            best_C=best_C,
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

    def _evaluate_C(self, X: np.ndarray, y: np.ndarray, C: float) -> float:
        """
        Evaluate a specific C value via cross-validation.

        Args:
            X: Feature matrix
            y: Labels
            C: C value to evaluate

        Returns:
            Cross-validation balanced accuracy score
        """
        svm = SVC(
            C=C,
            kernel='linear',
            class_weight='balanced',
            cache_size=1000,  # MB - critical for performance with large datasets
            random_state=self.random_state
        )

        scores = cross_val_score(
            svm,
            X,
            y,
            cv=self.cv_folds,
            scoring='balanced_accuracy',
            n_jobs=self.n_jobs,  # Parallelize CV folds
            verbose=2  # Show real-time progress
        )

        return float(np.mean(scores))
