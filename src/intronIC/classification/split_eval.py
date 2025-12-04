"""
Train/validation/test split evaluation for model assessment.

This module implements a simpler alternative to nested CV using a single
train/validation/test split for honest model evaluation.

Key features:
- Single stratified split: train/validation/test
- Train on training set
- Optimize hyperparameters (GridSearchCV uses internal CV on train)
- Evaluate on held-out test set
- Faster than nested CV (single split instead of K folds)
"""

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from sklearn.metrics import average_precision_score, f1_score, precision_recall_curve
from sklearn.model_selection import train_test_split

from intronIC.classification.optimizer import SVMOptimizer
from intronIC.classification.predictor import SVMPredictor
from intronIC.classification.trainer import SVMTrainer
from intronIC.core.intron import Intron


@dataclass(frozen=True, slots=True)
class SplitEvalResult:
    """Results from train/val/test evaluation."""

    test_f1: float
    test_pr_auc: float
    precision: np.ndarray  # Precision values for PR curve
    recall: np.ndarray  # Recall values for PR curve
    n_u12_train: int
    n_u2_train: int
    n_u12_val: int
    n_u2_val: int
    n_u12_test: int
    n_u2_test: int
    train_fraction: float
    val_fraction: float
    test_fraction: float
    optimized_C: float
    calibration_method: str

    def __str__(self) -> str:
        """Format results for display."""
        lines = [
            "\n" + "=" * 80,
            "Train/Validation/Test Split Evaluation",
            "=" * 80,
            "",
            "Data Split:",
            f"  Train: {self.n_u2_train} U2, {self.n_u12_train} U12 "
            f"({self.train_fraction * 100:.0f}%)",
            f"  Val:   {self.n_u2_val} U2, {self.n_u12_val} U12 "
            f"({self.val_fraction * 100:.0f}%)",
            f"  Test:  {self.n_u2_test} U2, {self.n_u12_test} U12 "
            f"({self.test_fraction * 100:.0f}%)",
            "",
            "Hyperparameters:",
            f"  Optimized C: {self.optimized_C:.6e}",
            f"  Calibration: {self.calibration_method}",
            "",
            "Test Set Performance (Honest Evaluation):",
            f"  F1 Score:    {self.test_f1:.4f}",
            f"  PR-AUC:      {self.test_pr_auc:.4f}",
            "=" * 80,
            "",
        ]

        return "\n".join(lines)


class SplitEvaluator:
    """
    Train/validation/test split evaluator for model assessment.

    Simpler and faster alternative to nested CV. Uses a single split
    to get honest performance estimates before training production model.
    """

    def __init__(
        self,
        test_fraction: float = 0.2,
        val_fraction: float = 0.2,
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
        Initialize split evaluator.

        Args:
            test_fraction: Fraction for test set (default: 0.2)
            val_fraction: Fraction for validation set (default: 0.2)
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
        """
        self.test_fraction = test_fraction
        self.val_fraction = val_fraction
        self.train_fraction = 1.0 - test_fraction - val_fraction

        if self.train_fraction <= 0:
            raise ValueError(
                f"test_fraction ({test_fraction}) + val_fraction ({val_fraction}) "
                "must be less than 1.0"
            )

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
    ) -> SplitEvalResult:
        """
        Run train/val/test split evaluation.

        Args:
            u12_reference: Reference U12-type introns (normalized)
            u2_reference: Reference U2-type introns (normalized)

        Returns:
            SplitEvalResult with honest test set performance
        """
        if self.verbose:
            print("\n" + "=" * 80)
            print("Train/Validation/Test Split Evaluation")
            print("=" * 80)
            print(f"Reference data: {len(u12_reference)} U12, {len(u2_reference)} U2")
            print("=" * 80 + "\n")

        # Prepare data for stratified splitting
        all_introns = list(u2_reference) + list(u12_reference)
        labels = np.array([0] * len(u2_reference) + [1] * len(u12_reference))

        # First split: separate test set
        train_val_introns, test_introns, train_val_labels, test_labels = (
            train_test_split(
                all_introns,
                labels,
                test_size=self.test_fraction,
                stratify=labels,
                random_state=self.random_state,
            )
        )

        # Second split: separate validation from training
        # Adjust val fraction to account for already-removed test set
        val_fraction_adjusted = self.val_fraction / (1.0 - self.test_fraction)

        train_introns, val_introns, train_labels, val_labels = train_test_split(
            train_val_introns,
            train_val_labels,
            test_size=val_fraction_adjusted,
            stratify=train_val_labels,
            random_state=self.random_state + 1,
        )

        # Separate each set into U2 and U12
        train_u2 = [
            intron for intron, label in zip(train_introns, train_labels) if label == 0
        ]
        train_u12 = [
            intron for intron, label in zip(train_introns, train_labels) if label == 1
        ]

        val_u2 = [
            intron for intron, label in zip(val_introns, val_labels) if label == 0
        ]
        val_u12 = [
            intron for intron, label in zip(val_introns, val_labels) if label == 1
        ]

        n_u12_train = len(train_u12)
        n_u2_train = len(train_u2)
        n_u12_val = len(val_u12)
        n_u2_val = len(val_u2)
        n_u12_test = int(np.sum(test_labels == 1))
        n_u2_test = int(np.sum(test_labels == 0))

        if self.verbose:
            print(
                f"Train: {n_u2_train} U2, {n_u12_train} U12 ({self.train_fraction * 100:.0f}%)"
            )
            print(
                f"Val:   {n_u2_val} U2, {n_u12_val} U12 ({self.val_fraction * 100:.0f}%)"
            )
            print(
                f"Test:  {n_u2_test} U2, {n_u12_test} U12 ({self.test_fraction * 100:.0f}%)"
            )

        # Optimize hyperparameters on training set
        # Note: We could use train+val here, but using only train gives
        # a more conservative estimate
        # Even with fixed C, we optimize include_max/dual/intercept_scaling/calibration_method
        if self.verbose:
            print(f"\n{'=' * 80}")
            print("Hyperparameter Optimization (training set)")
            print(f"{'=' * 80}")

        # If C is fixed, constrain the grid to that single value
        param_grid = self.param_grid_override.copy() if self.param_grid_override else {}
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
            random_state=self.random_state,
            n_jobs=self.n_jobs,
            max_iter=self.max_iter,
            scoring_metric=self.scoring_metric,
            penalty_options=self.penalty_options,
            loss_options=self.loss_options,
            class_weight_multipliers=self.class_weight_multipliers,
            use_multiplier_tiebreaker=self.use_multiplier_tiebreaker,
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

        # Train ensemble on training data
        if self.verbose:
            print(f"\n{'=' * 80}")
            print("Model Training (training set)")
            print(f"{'=' * 80}")

        trainer = SVMTrainer(
            n_models=self.n_ensemble_models,
            random_state=self.random_state,
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

        # Predict on test set (honest evaluation)
        if self.verbose:
            print(f"\n{'=' * 80}")
            print("Test Evaluation (test set)")
            print(f"{'=' * 80}")

        predictor = SVMPredictor(
            threshold=self.classification_threshold, n_jobs=self.n_jobs
        )
        predicted_introns = list(predictor.predict(ensemble, test_introns))

        # Extract predictions and probabilities
        y_pred = []
        y_proba = []
        for intron in predicted_introns:
            if intron.scores and intron.scores.svm_score is not None:
                y_proba.append(intron.scores.svm_score / 100.0)  # Convert to 0-1 range
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
        precision, recall, _ = precision_recall_curve(test_labels, y_proba, pos_label=1)

        result = SplitEvalResult(
            test_f1=float(f1),
            test_pr_auc=float(pr_auc),
            precision=precision,
            recall=recall,
            n_u12_train=n_u12_train,
            n_u2_train=n_u2_train,
            n_u12_val=n_u12_val,
            n_u2_val=n_u2_val,
            n_u12_test=n_u12_test,
            n_u2_test=n_u2_test,
            train_fraction=self.train_fraction,
            val_fraction=self.val_fraction,
            test_fraction=self.test_fraction,
            optimized_C=parameters.C,
            calibration_method=parameters.calibration_method,
        )

        if self.verbose:
            print(result)

        return result
        return result
