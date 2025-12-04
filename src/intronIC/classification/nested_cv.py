"""
Nested cross-validation for honest model evaluation.

This module implements nested CV to get unbiased estimates of model performance
before training the final production model on all reference data.

Key features:
- Outer loop: Stratified K-fold splits for evaluation
- Inner loop: Hyperparameter optimization on each training fold
- Honest metrics: Test folds never seen during training or optimization
- Aggregated results: Mean and std across all folds
"""

from dataclasses import dataclass
from typing import List, Sequence, Tuple

import numpy as np
from sklearn.metrics import average_precision_score, f1_score, precision_recall_curve
from sklearn.model_selection import StratifiedKFold

from intronIC.classification.optimizer import SVMOptimizer
from intronIC.classification.predictor import SVMPredictor
from intronIC.classification.trainer import SVMTrainer
from intronIC.core.intron import Intron


@dataclass(frozen=True, slots=True)
class FoldResult:
    """Results from one CV fold."""

    fold_idx: int
    f1_score: float
    pr_auc: float
    precision: np.ndarray  # Precision values for PR curve
    recall: np.ndarray  # Recall values for PR curve
    n_u12_train: int
    n_u2_train: int
    n_u12_test: int
    n_u2_test: int
    optimized_C: float
    calibration_method: str


@dataclass(frozen=True, slots=True)
class NestedCVResult:
    """Aggregated nested CV results."""

    fold_results: Sequence[FoldResult]
    mean_f1: float
    std_f1: float
    mean_pr_auc: float
    std_pr_auc: float
    n_folds: int
    pr_curves: List[
        Tuple[np.ndarray, np.ndarray]
    ]  # List of (precision, recall) tuples from all folds

    def __str__(self) -> str:
        """Format results for display."""
        lines = [
            "\n" + "=" * 80,
            "Nested Cross-Validation Results (Honest Evaluation)",
            "=" * 80,
        ]

        # Individual fold results
        for fold in self.fold_results:
            lines.append(
                f"Fold {fold.fold_idx + 1}/{self.n_folds}: "
                f"F1={fold.f1_score:.4f}, PR-AUC={fold.pr_auc:.4f}"
            )

        # Summary statistics
        lines.extend(
            [
                "",
                "Summary:",
                f"  Mean F1:     {self.mean_f1:.4f} ± {self.std_f1:.4f}",
                f"  Mean PR-AUC: {self.mean_pr_auc:.4f} ± {self.std_pr_auc:.4f}",
                "=" * 80,
                "",
            ]
        )

        return "\n".join(lines)


class NestedCVEvaluator:
    """
    Nested cross-validation evaluator for honest performance estimation.

    Outer loop: Stratified K-fold on reference data
    Inner loop: Hyperparameter optimization + training on each fold

    This provides unbiased estimates of model performance before
    training the final production model on all data.
    """

    def __init__(
        self,
        n_folds: int = 7,
        n_optimization_rounds: int = 3,
        n_ensemble_models: int = 1,
        classification_threshold: float = 90.0,
        subsample_u2: bool = False,
        subsample_ratio: float = 0.8,
        random_state: int = 42,
        n_jobs: int = 1,
        max_iter: int = 100000,
        verbose: bool = True,
        optimize_c: bool = True,
        fixed_c: float | None = None,
        cv_folds: int = 5,
        n_points_initial: int = 13,
        scoring_metric: str = "balanced_accuracy",
        penalty_options: list | None = None,
        loss_options: list | None = None,
        class_weight_multipliers: list | None = None,
        use_multiplier_tiebreaker: bool = True,
        param_grid_override: dict | None = None,
        eff_C_pos_range: tuple = (1e-3, 1e3),
        eff_C_neg_max: float | None = None,
        progress_tracker=None,
        features_list: list | None = None,
    ):
        """
        Initialize nested CV evaluator.

        Args:
            n_folds: Number of CV folds (default: 5)
            n_optimization_rounds: Hyperparameter search rounds (default: 3)
            n_ensemble_models: Models per ensemble (default: 1 for speed)
            classification_threshold: U12 probability threshold (default: 90.0)
            subsample_u2: Whether to subsample U2 for ensemble (default: False for speed)
            subsample_ratio: U2 subsample ratio (default: 0.8)
            random_state: Random seed
            n_jobs: Parallel jobs for optimization/prediction
            max_iter: Max iterations for LinearSVC
            verbose: Print progress
            optimize_c: Whether to optimize C parameter (default: True)
            fixed_c: Fixed C value if not optimizing (default: None)
            cv_folds: Cross-validation folds for GridSearchCV (default: 5)
            n_points_initial: Initial grid points for optimization (default: 13)
            scoring_metric: Metric for hyperparameter optimization (default: 'balanced_accuracy')
            penalty_options: Penalty types to search (default: None -> ['l2'])
            loss_options: Loss functions to search (default: None -> ['squared_hinge'])
            class_weight_multipliers: Class weight multipliers (default: None -> [1.0])
            use_multiplier_tiebreaker: Prefer 1.0 when multipliers tied (default: True)
            param_grid_override: Optional custom parameter grid (default: None)
            eff_C_pos_range: Target effective penalty range for positive class (default: 1e-3 to 1e3)
            eff_C_neg_max: Optional cap on negative class effective penalty (default: None)
            progress_tracker: Optional ProgressTracker for global step counting
            features_list: List of composite feature names for BothEndsStrongTransformer (default: None)
        """
        self.n_folds = n_folds
        self.n_optimization_rounds = n_optimization_rounds
        self.n_ensemble_models = n_ensemble_models
        self.classification_threshold = classification_threshold
        self.subsample_u2 = subsample_u2
        self.subsample_ratio = subsample_ratio
        self.random_state = random_state
        self.n_jobs = n_jobs
        self.max_iter = max_iter
        self.verbose = verbose
        self.optimize_c = optimize_c
        self.fixed_c = fixed_c
        self.cv_folds = cv_folds
        self.n_points_initial = n_points_initial
        self.scoring_metric = scoring_metric
        self.penalty_options = penalty_options
        self.loss_options = loss_options
        self.class_weight_multipliers = class_weight_multipliers
        self.use_multiplier_tiebreaker = use_multiplier_tiebreaker
        self.param_grid_override = param_grid_override
        self.eff_C_pos_range = eff_C_pos_range
        self.eff_C_neg_max = eff_C_neg_max
        self.progress_tracker = progress_tracker
        self.features_list = features_list

    def evaluate(
        self, u12_reference: Sequence[Intron], u2_reference: Sequence[Intron]
    ) -> NestedCVResult:
        """
        Run nested cross-validation evaluation.

        Args:
            u12_reference: Reference U12-type introns (normalized)
            u2_reference: Reference U2-type introns (normalized)

        Returns:
            NestedCVResult with honest performance estimates
        """
        if self.verbose:
            print("\n" + "=" * 80)
            print(f"Nested Cross-Validation Evaluation ({self.n_folds} folds)")
            print("=" * 80)
            print(f"Reference data: {len(u12_reference)} U12, {len(u2_reference)} U2")
            print("=" * 80 + "\n")

        # Prepare data for stratified splitting
        # Combine U2 and U12, create labels
        all_introns = list(u2_reference) + list(u12_reference)
        labels = np.array([0] * len(u2_reference) + [1] * len(u12_reference))

        # Outer loop: Stratified K-fold
        cv_splitter = StratifiedKFold(
            n_splits=self.n_folds, shuffle=True, random_state=self.random_state
        )

        fold_results = []

        for fold_idx, (train_indices, test_indices) in enumerate(
            cv_splitter.split(all_introns, labels)
        ):
            if self.verbose:
                print(f"\n{'=' * 80}")
                print(f"FOLD {fold_idx + 1}/{self.n_folds}")
                print(f"{'=' * 80}")

            # Split into train and test for this fold
            train_introns = [all_introns[i] for i in train_indices]
            test_introns = [all_introns[i] for i in test_indices]
            train_labels = labels[train_indices]
            test_labels = labels[test_indices]

            # Separate train into U2 and U12
            train_u2 = [
                intron
                for intron, label in zip(train_introns, train_labels)
                if label == 0
            ]
            train_u12 = [
                intron
                for intron, label in zip(train_introns, train_labels)
                if label == 1
            ]

            n_u12_train = len(train_u12)
            n_u2_train = len(train_u2)
            n_u12_test = int(np.sum(test_labels == 1))
            n_u2_test = int(np.sum(test_labels == 0))

            if self.verbose:
                print(f"Train: {n_u2_train} U2, {n_u12_train} U12")
                print(f"Test:  {n_u2_test} U2, {n_u12_test} U12")

            # Hyperparameter optimization for this fold
            # Even with fixed C, we optimize include_max/dual/intercept_scaling/calibration_method
            if self.verbose:
                print(f"\n{'─' * 80}")
                print(
                    f"FOLD {fold_idx + 1}/{self.n_folds} - Hyperparameter Optimization"
                )
                print(f"{'─' * 80}")

            # If C is fixed, constrain the grid to that single value
            param_grid = (
                self.param_grid_override.copy() if self.param_grid_override else {}
            )
            if not self.optimize_c:
                param_grid["estimator__svc__C"] = [self.fixed_c]
                if self.verbose:
                    print(
                        f"C fixed at {self.fixed_c:.6e}, optimizing include_max/dual/intercept_scaling/calibration_method"
                    )

            optimizer = SVMOptimizer(
                n_rounds=self.n_optimization_rounds,
                n_points_initial=self.n_points_initial,
                cv_folds=self.cv_folds,
                random_state=self.random_state + fold_idx,
                n_jobs=self.n_jobs,
                max_iter=self.max_iter,
                scoring_metric=self.scoring_metric,
                penalty_options=self.penalty_options,
                loss_options=self.loss_options,
                class_weight_multipliers=self.class_weight_multipliers,
                param_grid_override=param_grid if param_grid else None,
                verbose=self.verbose,
                progress_tracker=self.progress_tracker,
                features_list=self.features_list,
            )
            parameters = optimizer.optimize(
                train_u12,
                train_u2,
                eff_C_pos_range=self.eff_C_pos_range,
                eff_C_neg_max=self.eff_C_neg_max,
            )

            # Train ensemble on training fold
            if self.verbose:
                print(f"\n{'─' * 80}")
                print(f"FOLD {fold_idx + 1}/{self.n_folds} - Model Training")
                print(f"{'─' * 80}")

            trainer = SVMTrainer(
                n_models=self.n_ensemble_models,
                random_state=self.random_state + fold_idx,
                max_iter=self.max_iter,
                progress_tracker=self.progress_tracker,
                features_list=self.features_list,
            )
            ensemble = trainer.train_ensemble(
                train_u12,
                train_u2,
                parameters,
                subsample_u2=self.subsample_u2,
                subsample_ratio=self.subsample_ratio,
            )

            # Predict on test fold (honest evaluation)
            if self.verbose:
                print(f"\n{'─' * 80}")
                print(f"FOLD {fold_idx + 1}/{self.n_folds} - Test Evaluation")
                print(f"{'─' * 80}")

            predictor = SVMPredictor(
                threshold=self.classification_threshold, n_jobs=self.n_jobs
            )
            predicted_introns = list(predictor.predict(ensemble, test_introns))

            # Extract predictions and probabilities
            y_pred = []
            y_proba = []
            for intron in predicted_introns:
                if intron.scores and intron.scores.svm_score is not None:
                    y_proba.append(
                        intron.scores.svm_score / 100.0
                    )  # Convert to 0-1 range
                    y_pred.append(1 if intron.type_id == "u12" else 0)
                else:
                    # Shouldn't happen, but handle gracefully
                    y_proba.append(0.0)
                    y_pred.append(0)

            y_pred = np.array(y_pred)
            y_proba = np.array(y_proba)

            # Calculate metrics
            f1 = f1_score(test_labels, y_pred, pos_label=1)
            pr_auc = average_precision_score(test_labels, y_proba)

            # Compute precision-recall curve for plotting
            precision, recall, _ = precision_recall_curve(
                test_labels, y_proba, pos_label=1
            )

            if self.verbose:
                print(f"Fold {fold_idx + 1} Results: F1={f1:.4f}, PR-AUC={pr_auc:.4f}")

            fold_results.append(
                FoldResult(
                    fold_idx=fold_idx,
                    f1_score=float(f1),
                    pr_auc=float(pr_auc),
                    precision=precision,
                    recall=recall,
                    n_u12_train=n_u12_train,
                    n_u2_train=n_u2_train,
                    n_u12_test=n_u12_test,
                    n_u2_test=n_u2_test,
                    optimized_C=parameters.C,
                    calibration_method=parameters.calibration_method,
                )
            )

            # Update global progress after fold evaluation
            if self.progress_tracker:
                self.progress_tracker.increment(
                    f"Completed fold {fold_idx + 1}/{self.n_folds} evaluation"
                )

        # Aggregate results
        f1_scores = [fold.f1_score for fold in fold_results]
        pr_aucs = [fold.pr_auc for fold in fold_results]
        pr_curves = [(fold.precision, fold.recall) for fold in fold_results]

        result = NestedCVResult(
            fold_results=fold_results,
            mean_f1=float(np.mean(f1_scores)),
            std_f1=float(np.std(f1_scores)),
            mean_pr_auc=float(np.mean(pr_aucs)),
            std_pr_auc=float(np.std(pr_aucs)),
            n_folds=self.n_folds,
            pr_curves=pr_curves,
        )

        if self.verbose:
            print(result)

        return result
        return result
