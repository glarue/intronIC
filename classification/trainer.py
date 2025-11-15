"""
SVM ensemble training with U2 subsampling for class imbalance.

This module implements the training algorithm from intronIC.py:5345-5430,
which trains multiple SVM models with different U2 subsamples to create
a diverse ensemble for robust classification.

Key features:
- Balanced class weights to handle U12/U2 imbalance
- Multiple models with different U2 subsamples for diversity
- External calibration for probability estimation

Evaluation metrics are computed separately via nested CV or split evaluation modules.

Port from: intronIC.py:5345-5430 (train_svm)
"""

from dataclasses import dataclass
from typing import Sequence, Tuple, Optional, Any
import warnings
import contextlib
import numpy as np
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
from sklearn.base import clone
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import RobustScaler
from sklearn.exceptions import ConvergenceWarning

from core.intron import Intron
from classification.transformers import BothEndsStrongTransformer
from classification.optimizer import SVMParameters

# Global filter for convergence warnings (persists across multiprocessing forks)
warnings.filterwarnings("ignore", category=ConvergenceWarning)


@contextlib.contextmanager
def suppress_convergence_warnings(verbose: bool = True):
    """
    Context manager to suppress sklearn ConvergenceWarning spam.

    During model training with CalibratedClassifierCV (5 folds), convergence warnings
    can spam the console with many identical messages. This context manager suppresses
    them to keep output clean.

    Args:
        verbose: If True, print note about suppression (default: True)

    Usage:
        with suppress_convergence_warnings(verbose=True):
            model.fit(X, y)
    """
    with warnings.catch_warnings():
        # Suppress convergence warnings
        warnings.filterwarnings("ignore", category=ConvergenceWarning)
        yield


@dataclass(frozen=True, slots=True)
class SVMModel:
    """Trained SVM model with metadata."""

    model: Any  # Trained sklearn model (LinearSVC with external calibration)
    train_size: int  # Number of training samples
    u12_count: int  # U12 introns in training
    u2_count: int  # U2 introns in training
    parameters: SVMParameters  # Hyperparameters used


@dataclass(frozen=True, slots=True)
class SVMEnsemble:
    """Collection of trained models for ensemble prediction."""

    models: Sequence[SVMModel]

    def __len__(self) -> int:
        return len(self.models)


class SVMTrainer:
    """
    Train ensemble of SVM classifiers with U2 subsampling.

    Handles class imbalance by:
    - Training multiple models with different U2 subsamples
    - Using balanced class weights

    Evaluation metrics are computed separately via nested CV or split evaluation.

    Port from: intronIC.py:5345-5430
    """

    def __init__(
        self,
        n_models: int = 3,
        random_state: int = 42,
        kernel: str = 'linear',
        max_iter: int = 20000
    ):
        """
        Initialize trainer.

        Args:
            n_models: Number of models in ensemble (default: 3)
            random_state: Random seed
            kernel: SVM kernel type (default: 'linear')
            max_iter: Maximum iterations for LinearSVC convergence (default: 20000)
        """
        self.n_models = n_models
        self.random_state = random_state
        self.kernel = kernel
        self.max_iter = max_iter

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

        print(f"Ensemble training complete: {len(models)} models trained")

        return SVMEnsemble(models=models)

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

        # Train LinearSVC with external calibration
        # Following sklearn best practices for rare-class classification:
        # - RobustScaler(with_centering=False): Scales by IQR while preserving semantic zero
        #   (s=0 means "U12≈U2", centering would destroy this meaning)
        # - BothEndsStrongTransformer: Augments 3D → 5D (or 7D) with both-ends-strong features
        #   Adds min/max features for 5'SS-BPS and 5'SS-3'SS pairs
        #   Min features capture "both must be strong" (expert: "model will mostly weight min")
        # - LinearSVC: Faster than SVC(kernel='linear'), optimized for linear case
        #   - dual=False: Primal formulation for n_samples >> n_features (21k >> 5-7)
        #   - intercept_scaling=1000: High value to avoid over-regularizing intercept
        #   - max_iter: Configurable iteration limit (default: 100000)
        #   - tol=1e-4: Tight convergence tolerance
        # - class_weight='balanced': Handle 1.8% positive class
        # - CalibratedClassifierCV: Add external calibration (sigmoid or isotonic)
        #   Calibration method chosen by optimizer via grid search
        base_pipeline = Pipeline([
            ('scale', RobustScaler(with_centering=False, with_scaling=True)),
            ('augment', BothEndsStrongTransformer(
                include_max=parameters.include_max
            )),
            ('svc', LinearSVC(
                C=parameters.C,
                dual=parameters.dual,  # From optimizer (typically False for our data)
                penalty=parameters.penalty,  # From optimizer: 'l1' or 'l2' (L1 prunes redundant features)
                loss=parameters.loss,  # From optimizer: 'squared_hinge' (required for L1)
                intercept_scaling=parameters.intercept_scaling,  # High value when dual=False
                class_weight='balanced',  # Critical for imbalanced data
                max_iter=self.max_iter,  # Maximum iterations for convergence
                tol=1e-4,  # Tighter tolerance
                random_state=seed
            ))
        ])

        # External calibration wrapper
        # Method (sigmoid vs isotonic) chosen by hyperparameter optimization
        # ensemble='auto' uses per-fold fit + averaging for stable tails
        svm = CalibratedClassifierCV(
            base_pipeline,
            method=parameters.calibration_method,  # From optimizer: 'sigmoid' or 'isotonic'
            cv=5,  # Stratified 5-fold
            ensemble='auto'  # Per-fold calibrators averaged
        )

        # Suppress convergence warning spam but log summary
        with suppress_convergence_warnings(verbose=True):
            svm.fit(X, y)

        return SVMModel(
            model=svm,
            train_size=len(X),
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
