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
from sklearn.exceptions import ConvergenceWarning
from sklearn.preprocessing import RobustScaler

from intronIC.core.intron import Intron
from intronIC.classification.transformers import BothEndsStrongTransformer
from intronIC.classification.optimizer import SVMParameters

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
    subsample_ratio: float = 0.85  # Fraction of U2s used per model (for display)

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
        max_iter: int = 20000,
        features_list: Optional[list] = None,
        gamma_imbalance_options: Optional[list] = None,
        progress_tracker: Optional[Any] = None
    ):
        """
        Initialize trainer.

        Args:
            n_models: Number of models in ensemble (default: 3)
            random_state: Random seed
            kernel: SVM kernel type (default: 'linear')
            max_iter: Maximum iterations for LinearSVC convergence (default: 20000)
            features_list: List of composite feature names to include (default: None = use default 7D)
            gamma_imbalance_options: Gamma scaling factors (default: None, uses parameters.gamma_imbalance from SVMParameters)
            progress_tracker: Optional ProgressTracker for global step counting
        """
        self.n_models = n_models
        self.random_state = random_state
        self.kernel = kernel
        self.max_iter = max_iter
        self.features_list = features_list
        self.gamma_imbalance_options = gamma_imbalance_options
        self.progress_tracker = progress_tracker

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

        # Print training header (no stage number - context-dependent)
        print(f"\n{'='*80}", flush=True)
        if self.n_models == 1:
            print(f"Model Training", flush=True)
        else:
            print(f"Ensemble Training ({self.n_models} models)", flush=True)
        print(f"{'='*80}\n", flush=True)

        for i in range(self.n_models):
            print(f"\n{'─'*80}", flush=True)
            global_step_str = f" {self.progress_tracker.format_step()}" if self.progress_tracker else ""
            if self.n_models == 1:
                print(f"Training model...{global_step_str}", flush=True)
            else:
                print(f"MODEL {i+1}/{self.n_models}: Training ensemble model...{global_step_str}", flush=True)
            print(f"{'─'*80}", flush=True)

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

            # Update global progress
            if self.progress_tracker:
                self.progress_tracker.increment(f"Completed model {i+1}/{self.n_models}")

        print(f"\n{'='*80}", flush=True)
        if self.n_models == 1:
            print(f"Model training complete", flush=True)
        else:
            print(f"Ensemble training complete ({len(models)} models trained)", flush=True)
        print(f"{'='*80}\n", flush=True)

        return SVMEnsemble(models=models, subsample_ratio=subsample_ratio)

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

        # Compute balanced class weights with multiplier (2025-01-19: robustness improvements)
        n_samples = len(y)
        n_pos = int(np.sum(y == 1))  # U12
        n_neg = int(np.sum(y == 0))  # U2
        w_pos = (n_samples / (2.0 * n_pos)) * parameters.class_weight_multiplier
        w_neg = (n_samples / (2.0 * n_neg)) * parameters.class_weight_multiplier
        class_weight = {0: w_neg, 1: w_pos}

        # CORRECTED ARCHITECTURE (2025 - Expert guidance)
        #
        # Pipeline: z-scores → augmented features → svc
        #
        # Scaling happens OUTSIDE pipeline via ScoreNormalizer:
        # - ScoreNormalizer: RobustScaler(with_centering=True) fitted on reference LLRs
        #   Domain adaptation: Refit per-species (or reuse human for cross-species)
        # - Trainer receives pre-scaled z-scores: [five_z_score, bp_z_score, three_z_score]
        #
        # Pipeline steps:
        # - BothEndsStrongTransformer: Augmented features from z-scores
        #   → min_all, neg_absdiff_5_bp, neg_absdiff_5_3 (suppress one-end-strong FPs)
        # - LinearSVC: With penalty ∈ {l1, l2}, loss ∈ {hinge, squared_hinge}
        # - CalibratedClassifierCV: External calibration (isotonic or sigmoid)
        #
        # Key principle: Single scaling step (NOT double-scaling)
        # See: Expert workflow doc, SCALER_ARCHITECTURE_REVIEW.md

        base_pipeline = Pipeline([
            ('transform', BothEndsStrongTransformer(
                features=self.features_list,
                gamma_imbalance=parameters.gamma_imbalance
            )),
            ('svc', LinearSVC(
                C=parameters.C,
                penalty=parameters.penalty,
                dual=False,
                loss=parameters.loss,
                class_weight=class_weight,
                max_iter=self.max_iter,
                tol=1e-4,
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
        Stratified subsample of U2 introns by length and GC content (2025-01-19: robustness).

        Creates 2D bins (length × GC content) and samples proportionally from each bin
        to ensure diverse representation across both dimensions. This prevents models
        from being biased toward specific length ranges or GC compositions.

        Length bins:
        - <100 bp, 100-500 bp, 500-1000 bp, 1000-5000 bp, ≥5000 bp

        GC content bins:
        - Low: <35% GC
        - Medium-Low: 35-45% GC
        - Medium: 45-55% GC
        - Medium-High: 55-65% GC
        - High: ≥65% GC

        Args:
            u2_introns: Full U2 intron set
            seed: Random seed for reproducibility
            subsample_ratio: Fraction of U2s to use (default: 0.8)

        Returns:
            Stratified sample of U2 introns

        Port from: intronIC.py:5356-5366 (subsetting logic)
        """
        np.random.seed(seed)

        # Define length bins (in bp)
        length_bins = [
            (0, 100),
            (100, 500),
            (500, 1000),
            (1000, 5000),
            (5000, float('inf'))
        ]

        # Define GC content bins (percentage)
        gc_bins = [
            (0.0, 0.35),   # Low GC
            (0.35, 0.45),  # Medium-Low GC
            (0.45, 0.55),  # Medium GC
            (0.55, 0.65),  # Medium-High GC
            (0.65, 1.0)    # High GC
        ]

        # Helper function to calculate GC content
        def calc_gc_content(intron: Intron) -> float:
            """Calculate GC content from intron sequence."""
            if intron.sequences is None or intron.sequences.seq is None:
                return 0.5  # Default to medium GC if no sequence available
            seq = intron.sequences.seq.upper()
            if len(seq) == 0:
                return 0.5
            gc_count = seq.count('G') + seq.count('C')
            return gc_count / len(seq)

        # Bin introns by length AND GC content (2D binning)
        # binned_introns[length_idx][gc_idx] = list of introns
        binned_introns = [[[] for _ in gc_bins] for _ in length_bins]

        for intron in u2_introns:
            length = intron.length
            gc_content = calc_gc_content(intron)

            # Find length bin
            length_idx = None
            for idx, (min_len, max_len) in enumerate(length_bins):
                if min_len <= length < max_len:
                    length_idx = idx
                    break

            # Find GC bin
            gc_idx = None
            for idx, (min_gc, max_gc) in enumerate(gc_bins):
                if min_gc <= gc_content < max_gc:
                    gc_idx = idx
                    break

            # Add to 2D bin (if both indices found)
            if length_idx is not None and gc_idx is not None:
                binned_introns[length_idx][gc_idx].append(intron)

        # Sample proportionally from each 2D bin
        sampled_introns = []
        for length_idx in range(len(length_bins)):
            for gc_idx in range(len(gc_bins)):
                bin_introns = binned_introns[length_idx][gc_idx]
                if len(bin_introns) == 0:
                    continue

                # Calculate how many to sample from this bin (proportional to bin size)
                bin_fraction = len(bin_introns) / len(u2_introns)
                n_bin_samples = max(1, int(bin_fraction * len(u2_introns) * subsample_ratio))

                # Don't sample more than available
                n_bin_samples = min(n_bin_samples, len(bin_introns))

                # Sample from this bin
                if n_bin_samples < len(bin_introns):
                    indices = np.random.choice(
                        len(bin_introns),
                        size=n_bin_samples,
                        replace=False
                    )
                    sampled_introns.extend([bin_introns[i] for i in indices])
                else:
                    # Use all introns if bin is small
                    sampled_introns.extend(bin_introns)

        return sampled_introns

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
