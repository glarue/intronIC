"""
SVM ensemble training with U2 subsampling for class imbalance.

This module implements the training algorithm from intronIC.py:5345-5430,
which trains multiple SVM models with different U2 subsamples to create
a diverse ensemble for robust classification.

Key features:
- Stratified train/test splits to maintain class proportions
- Balanced class weights to handle U12/U2 imbalance
- Multiple models with different U2 subsamples for diversity
- F1 score and Precision-Recall AUC for evaluation

Port from: intronIC.py:5345-5430 (train_svm)
"""

from dataclasses import dataclass
from typing import Sequence, Tuple, Optional
import numpy as np
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score, average_precision_score
from sklearn.base import clone

from core.intron import Intron
from classification.optimizer import SVMParameters


@dataclass(frozen=True, slots=True)
class SVMModel:
    """Trained SVM model with metadata."""

    model: SVC  # Trained sklearn SVC
    f1_score: float  # F1 on test set
    precision_recall_auc: float  # PR-AUC on test set
    train_size: int  # Number of training samples
    test_size: int  # Number of test samples
    u12_count: int  # U12 introns in training
    u2_count: int  # U2 introns in training
    parameters: SVMParameters  # Hyperparameters used


@dataclass(frozen=True, slots=True)
class SVMEnsemble:
    """Collection of trained models for ensemble prediction."""

    models: Sequence[SVMModel]
    mean_f1: float
    mean_pr_auc: float

    def __len__(self) -> int:
        return len(self.models)


class SVMTrainer:
    """
    Train ensemble of SVM classifiers with U2 subsampling.

    Handles class imbalance by:
    - Training multiple models with different U2 subsamples
    - Using balanced class weights
    - Stratified train/test splits

    Port from: intronIC.py:5345-5430
    """

    def __init__(
        self,
        n_models: int = 3,
        test_size: float = 0.2,
        random_state: int = 42,
        kernel: str = 'linear'
    ):
        """
        Initialize trainer.

        Args:
            n_models: Number of models in ensemble (default: 3)
            test_size: Fraction for test set (default: 0.2)
            random_state: Random seed
            kernel: SVM kernel type (default: 'linear')
        """
        self.n_models = n_models
        self.test_size = test_size
        self.random_state = random_state
        self.kernel = kernel

    def train_ensemble(
        self,
        u12_introns: Sequence[Intron],
        u2_introns: Sequence[Intron],
        parameters: SVMParameters,
        subsample_u2: bool = True,
        subsample_ratio: float = 0.8
    ) -> SVMEnsemble:
        """
        Train ensemble of SVM models.

        Args:
            u12_introns: Reference U12-type introns (normalized)
            u2_introns: Reference U2-type introns (normalized)
            parameters: Optimized hyperparameters
            subsample_u2: Whether to subsample U2 set for each model
            subsample_ratio: Fraction of U2s to use per model (default: 0.8)

        Returns:
            Ensemble of trained models
        """
        models = []

        for i in range(self.n_models):
            print(f"Training model {i+1}/{self.n_models}...")

            # Subsample U2 if requested (for diversity)
            if subsample_u2 and self.n_models > 1:
                u2_sample = self._subsample_u2(
                    u2_introns,
                    seed=self.random_state + i,
                    subsample_ratio=subsample_ratio
                )
            else:
                u2_sample = u2_introns

            # Train single model
            model = self._train_single_model(
                u12_introns,
                u2_sample,
                parameters,
                seed=self.random_state + i
            )
            models.append(model)

        # Calculate ensemble metrics
        mean_f1 = float(np.mean([m.f1_score for m in models]))
        mean_pr_auc = float(np.mean([m.precision_recall_auc for m in models]))

        print(f"Ensemble training complete: mean F1={mean_f1:.4f}, mean PR-AUC={mean_pr_auc:.4f}")

        return SVMEnsemble(
            models=models,
            mean_f1=mean_f1,
            mean_pr_auc=mean_pr_auc
        )

    def _train_single_model(
        self,
        u12_introns: Sequence[Intron],
        u2_introns: Sequence[Intron],
        parameters: SVMParameters,
        seed: int
    ) -> SVMModel:
        """
        Train a single SVM model.

        Port from: intronIC.py:5353-5428
        """
        # Prepare data
        X, y = self._prepare_training_data(u12_introns, u2_introns)

        # Stratified train/test split
        X_train, X_test, y_train, y_test = train_test_split(
            X, y,
            test_size=self.test_size,
            stratify=y,
            random_state=seed
        )

        # Train SVM with balanced class weights
        svm = SVC(
            C=parameters.C,
            kernel=self.kernel,
            class_weight='balanced',
            probability=True,  # Enable predict_proba
            cache_size=1000,  # MB - critical for performance with large datasets
            random_state=seed
        )
        svm.fit(X_train, y_train)

        # Evaluate on test set
        y_pred = svm.predict(X_test)
        y_proba = svm.predict_proba(X_test)[:, 1]  # U12 probabilities

        f1 = f1_score(y_test, y_pred, pos_label=1)  # 1 = U12
        pr_auc = average_precision_score(y_test, y_proba)

        return SVMModel(
            model=svm,
            f1_score=float(f1),
            precision_recall_auc=float(pr_auc),
            train_size=len(X_train),
            test_size=len(X_test),
            u12_count=int(np.sum(y == 1)),
            u2_count=int(np.sum(y == 0)),
            parameters=parameters
        )

    def _subsample_u2(
        self,
        u2_introns: Sequence[Intron],
        seed: int,
        subsample_ratio: float = 0.8
    ) -> Sequence[Intron]:
        """
        Randomly subsample U2 introns for diversity.

        Port from: intronIC.py:5356-5366 (subsetting logic)
        """
        np.random.seed(seed)

        n_samples = int(len(u2_introns) * subsample_ratio)
        indices = np.random.choice(
            len(u2_introns),
            size=n_samples,
            replace=False
        )

        return [u2_introns[i] for i in indices]

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
        # Port from: intronIC.py:5380, 5397
        X = np.array(u2_features + u12_features)
        y = np.array([0] * len(u2_features) + [1] * len(u12_features))

        return X, y
