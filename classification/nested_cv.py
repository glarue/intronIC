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
from typing import Sequence, List, Tuple
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import f1_score, average_precision_score, precision_recall_curve

from core.intron import Intron
from classification.optimizer import SVMOptimizer
from classification.trainer import SVMTrainer
from classification.predictor import SVMPredictor


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
    pr_curves: List[Tuple[np.ndarray, np.ndarray]]  # List of (precision, recall) tuples from all folds

    def __str__(self) -> str:
        """Format results for display."""
        lines = [
            "\n" + "="*80,
            "Nested Cross-Validation Results (Honest Evaluation)",
            "="*80,
        ]

        # Individual fold results
        for fold in self.fold_results:
            lines.append(
                f"Fold {fold.fold_idx + 1}/{self.n_folds}: "
                f"F1={fold.f1_score:.4f}, PR-AUC={fold.pr_auc:.4f}"
            )

        # Summary statistics
        lines.extend([
            "",
            "Summary:",
            f"  Mean F1:     {self.mean_f1:.4f} ± {self.std_f1:.4f}",
            f"  Mean PR-AUC: {self.mean_pr_auc:.4f} ± {self.std_pr_auc:.4f}",
            "="*80,
            ""
        ])

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
        n_folds: int = 5,
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
        fixed_c: float | None = None
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

    def evaluate(
        self,
        u12_reference: Sequence[Intron],
        u2_reference: Sequence[Intron]
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
            print("\n" + "="*80)
            print(f"Nested Cross-Validation Evaluation ({self.n_folds} folds)")
            print("="*80)
            print(f"Reference data: {len(u12_reference)} U12, {len(u2_reference)} U2")
            print("="*80 + "\n")

        # Prepare data for stratified splitting
        # Combine U2 and U12, create labels
        all_introns = list(u2_reference) + list(u12_reference)
        labels = np.array([0] * len(u2_reference) + [1] * len(u12_reference))

        # Outer loop: Stratified K-fold
        cv_splitter = StratifiedKFold(
            n_splits=self.n_folds,
            shuffle=True,
            random_state=self.random_state
        )

        fold_results = []

        for fold_idx, (train_indices, test_indices) in enumerate(cv_splitter.split(all_introns, labels)):
            if self.verbose:
                print(f"\n{'='*80}")
                print(f"FOLD {fold_idx + 1}/{self.n_folds}")
                print(f"{'='*80}")

            # Split into train and test for this fold
            train_introns = [all_introns[i] for i in train_indices]
            test_introns = [all_introns[i] for i in test_indices]
            train_labels = labels[train_indices]
            test_labels = labels[test_indices]

            # Separate train into U2 and U12
            train_u2 = [intron for intron, label in zip(train_introns, train_labels) if label == 0]
            train_u12 = [intron for intron, label in zip(train_introns, train_labels) if label == 1]

            n_u12_train = len(train_u12)
            n_u2_train = len(train_u2)
            n_u12_test = int(np.sum(test_labels == 1))
            n_u2_test = int(np.sum(test_labels == 0))

            if self.verbose:
                print(f"Train: {n_u2_train} U2, {n_u12_train} U12")
                print(f"Test:  {n_u2_test} U2, {n_u12_test} U12")

            # Stage 1: Optimize hyperparameters or use fixed C
            if self.optimize_c:
                if self.verbose:
                    print(f"\nStage 1: Hyperparameter Optimization (fold {fold_idx + 1})")

                optimizer = SVMOptimizer(
                    n_rounds=self.n_optimization_rounds,
                    random_state=self.random_state + fold_idx,
                    n_jobs=self.n_jobs,
                    max_iter=self.max_iter,
                    verbose=self.verbose
                )
                parameters = optimizer.optimize(train_u12, train_u2)
            else:
                if self.verbose:
                    print(f"\nStage 1: Using Fixed C={self.fixed_c:.6e} (fold {fold_idx + 1})")

                from classification.optimizer import SVMParameters
                parameters = SVMParameters(
                    C=self.fixed_c,
                    calibration_method='sigmoid',
                    dual=False,
                    intercept_scaling=1000.0,
                    cv_score=0.0,  # Not computed when using fixed C
                    round_found=0   # Fixed, not optimized
                )

            # Stage 2: Train ensemble on training fold
            if self.verbose:
                print(f"\nStage 2: Training Ensemble (fold {fold_idx + 1})")

            trainer = SVMTrainer(
                n_models=self.n_ensemble_models,
                random_state=self.random_state + fold_idx,
                max_iter=self.max_iter
            )
            ensemble = trainer.train_ensemble(
                train_u12,
                train_u2,
                parameters,
                subsample_u2=self.subsample_u2,
                subsample_ratio=self.subsample_ratio
            )

            # Stage 3: Predict on test fold (honest evaluation)
            if self.verbose:
                print(f"\nStage 3: Evaluating on Test Fold (fold {fold_idx + 1})")

            predictor = SVMPredictor(
                threshold=self.classification_threshold,
                n_jobs=self.n_jobs
            )
            predicted_introns = list(predictor.predict(ensemble, test_introns))

            # Extract predictions and probabilities
            y_pred = []
            y_proba = []
            for intron in predicted_introns:
                if intron.scores and intron.scores.svm_score is not None:
                    y_proba.append(intron.scores.svm_score / 100.0)  # Convert to 0-1 range
                    y_pred.append(1 if intron.type_id == 'u12' else 0)
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

            if self.verbose:
                print(f"Fold {fold_idx + 1} Results: F1={f1:.4f}, PR-AUC={pr_auc:.4f}")

            fold_results.append(FoldResult(
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
                calibration_method=parameters.calibration_method
            ))

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
            pr_curves=pr_curves
        )

        if self.verbose:
            print(result)

        return result
