"""
High-level intron classification pipeline.

This module provides the IntronClassifier class, which orchestrates the
complete U2/U12 classification workflow:

1. Hyperparameter optimization (SVMOptimizer)
2. Ensemble training (SVMTrainer)
3. Classification (SVMPredictor)

CRITICAL: This module does NOT re-normalize z-scores after classification.
Z-scores must be computed from reference data before classification to
prevent data leakage (Issue #1 fix).

Pipeline stages:
- Input: Reference introns + experimental introns (all with z-scores)
- Optimize: Find best SVM hyperparameters via cross-validation
- Train: Build ensemble of SVM models with U2 subsampling
- Classify: Apply ensemble to experimental introns
- Output: Classified introns with svm_score, relative_score, type_id

Port from: intronIC.py:5038-5900 (main pipeline)
"""

from typing import Sequence, Optional, Tuple
from dataclasses import dataclass

from core.intron import Intron
from classification.optimizer import SVMOptimizer, SVMParameters
from classification.trainer import SVMTrainer, SVMEnsemble
from classification.predictor import SVMPredictor


@dataclass(frozen=True, slots=True)
class ClassificationResult:
    """
    Result of classification pipeline.

    Contains the trained ensemble and classified introns.

    Attributes:
        classified_introns: Experimental introns with classification scores
        ensemble: Trained SVM ensemble
        parameters: Optimized hyperparameters
        n_u12_reference: Number of U12 reference introns
        n_u2_reference: Number of U2 reference introns
    """

    classified_introns: Sequence[Intron]
    ensemble: SVMEnsemble
    parameters: SVMParameters
    n_u12_reference: int
    n_u2_reference: int

    def get_u12_predictions(self, threshold: float = 90.0) -> Sequence[Intron]:
        """
        Get introns classified as U12-type.

        Args:
            threshold: SVM score threshold (default: 90.0)

        Returns:
            Introns with metadata.type_id == 'u12' and svm_score >= threshold
        """
        return [
            i for i in self.classified_introns
            if i.metadata and i.metadata.type_id == 'u12'
            and i.scores and i.scores.svm_score and i.scores.svm_score >= threshold
        ]

    def get_u2_predictions(self, threshold: float = 90.0) -> Sequence[Intron]:
        """
        Get introns classified as U2-type.

        Args:
            threshold: SVM score threshold (default: 90.0)

        Returns:
            Introns with metadata.type_id == 'u2' and svm_score < threshold
        """
        return [
            i for i in self.classified_introns
            if i.metadata and i.metadata.type_id == 'u2'
            and i.scores and i.scores.svm_score and i.scores.svm_score < threshold
        ]


class IntronClassifier:
    """
    High-level orchestrator for U2/U12 intron classification.

    Integrates hyperparameter optimization, ensemble training, and
    classification into a single pipeline.

    CRITICAL: Does NOT re-normalize z-scores (prevents data leakage).
    All introns must have z-scores computed from reference data before
    being passed to classify().

    Example:
        >>> # After scoring and normalization
        >>> classifier = IntronClassifier(
        ...     n_optimization_rounds=5,
        ...     n_ensemble_models=3,
        ...     classification_threshold=90.0
        ... )
        >>> result = classifier.classify(
        ...     u12_reference=u12_introns,  # Must have z-scores
        ...     u2_reference=u2_introns,    # Must have z-scores
        ...     experimental=experimental   # Must have z-scores
        ... )
        >>> u12_predictions = result.get_u12_predictions()

    Port from: intronIC.py:5038-5900
    """

    def __init__(
        self,
        n_optimization_rounds: int = 5,
        n_ensemble_models: int = 3,
        classification_threshold: float = 90.0,
        subsample_u2: bool = True,
        subsample_ratio: float = 0.8,
        random_state: int = 42,
        optimize_c: bool = True,
        fixed_c: Optional[float] = None,
        cv_processes: int = 1,
        classification_processes: int = 1
    ):
        """
        Initialize classifier.

        Args:
            n_optimization_rounds: Number of grid search refinement rounds (default: 5)
            n_ensemble_models: Number of models in ensemble (default: 3)
            classification_threshold: U12 probability threshold 0-100 (default: 90.0)
            subsample_u2: Whether to subsample U2 for ensemble diversity (default: True)
            subsample_ratio: Fraction of U2s per model (default: 0.8)
            random_state: Random seed for reproducibility (default: 42)
            optimize_c: Whether to optimize C parameter (default: True)
            fixed_c: Fixed C value if not optimizing (default: None)
            cv_processes: Number of parallel jobs for cross-validation (default: 1)
            classification_processes: Number of parallel jobs for classification (default: 1)
        """
        self.n_optimization_rounds = n_optimization_rounds
        self.n_ensemble_models = n_ensemble_models
        self.classification_threshold = classification_threshold
        self.subsample_u2 = subsample_u2
        self.subsample_ratio = subsample_ratio
        self.random_state = random_state
        self.optimize_c = optimize_c
        self.fixed_c = fixed_c
        self.cv_processes = cv_processes
        self.classification_processes = classification_processes

        # Validate parameters
        if not 0 <= classification_threshold <= 100:
            raise ValueError(
                f"classification_threshold must be 0-100, got {classification_threshold}"
            )
        if not optimize_c and fixed_c is None:
            raise ValueError("Must provide fixed_c if optimize_c is False")

    def classify(
        self,
        u12_reference: Sequence[Intron],
        u2_reference: Sequence[Intron],
        experimental: Sequence[Intron]
    ) -> ClassificationResult:
        """
        Run complete classification pipeline.

        CRITICAL: All introns must already have z-scores computed from
        reference data. This method does NOT re-normalize.

        Args:
            u12_reference: Reference U12 introns (with z-scores)
            u2_reference: Reference U2 introns (with z-scores)
            experimental: Introns to classify (with z-scores)

        Returns:
            ClassificationResult with classified introns and trained ensemble

        Raises:
            ValueError: If introns lack z-scores
        """
        # Validate inputs
        self._validate_introns_have_zscores(u12_reference, "u12_reference")
        self._validate_introns_have_zscores(u2_reference, "u2_reference")
        self._validate_introns_have_zscores(experimental, "experimental")

        print(f"Classification pipeline starting...")
        print(f"  Reference: {len(u12_reference)} U12, {len(u2_reference)} U2")
        print(f"  Experimental: {len(experimental)} introns")

        # Stage 1: Optimize hyperparameters
        if self.optimize_c:
            print("\n=== Stage 1: Hyperparameter Optimization ===")
            optimizer = SVMOptimizer(
                n_rounds=self.n_optimization_rounds,
                random_state=self.random_state,
                n_jobs=self.cv_processes
            )
            parameters = optimizer.optimize(u12_reference, u2_reference)
            print(f"Optimized C={parameters.C:.6e}, CV score={parameters.cv_score:.4f}")
        else:
            print("\n=== Stage 1: Using Fixed C Parameter ===")
            parameters = SVMParameters(
                C=self.fixed_c,
                cv_score=0.0,  # Not computed
                round_found=0   # Fixed, not optimized
            )
            print(f"Using fixed C={parameters.C:.6e}")

        # Stage 2: Train ensemble
        print("\n=== Stage 2: Ensemble Training ===")
        trainer = SVMTrainer(
            n_models=self.n_ensemble_models,
            random_state=self.random_state
        )
        ensemble = trainer.train_ensemble(
            u12_reference,
            u2_reference,
            parameters,
            subsample_u2=self.subsample_u2,
            subsample_ratio=self.subsample_ratio
        )
        print(f"Ensemble trained: {len(ensemble.models)} models")
        print(f"  Mean F1: {ensemble.mean_f1:.4f}")
        print(f"  Mean PR-AUC: {ensemble.mean_pr_auc:.4f}")

        # Stage 3: Classify experimental introns
        print("\n=== Stage 3: Classification ===")
        predictor = SVMPredictor(
            threshold=self.classification_threshold,
            n_jobs=self.classification_processes
        )
        classified = predictor.predict(ensemble, experimental)

        # Count classifications
        n_u12 = sum(
            1 for i in classified
            if i.metadata and i.metadata.type_id == 'u12'
        )
        n_u2 = len(classified) - n_u12
        print(f"Classification complete:")
        print(f"  U12: {n_u12} ({100*n_u12/len(classified):.1f}%)")
        print(f"  U2: {n_u2} ({100*n_u2/len(classified):.1f}%)")

        return ClassificationResult(
            classified_introns=classified,
            ensemble=ensemble,
            parameters=parameters,
            n_u12_reference=len(u12_reference),
            n_u2_reference=len(u2_reference)
        )

    def classify_batch(
        self,
        u12_reference: Sequence[Intron],
        u2_reference: Sequence[Intron],
        experimental: Sequence[Intron],
        batch_size: int = 10000
    ) -> ClassificationResult:
        """
        Run classification pipeline with batch processing.

        Useful for very large experimental datasets.

        Args:
            u12_reference: Reference U12 introns (with z-scores)
            u2_reference: Reference U2 introns (with z-scores)
            experimental: Introns to classify (with z-scores)
            batch_size: Classification batch size (default: 10000)

        Returns:
            ClassificationResult with classified introns
        """
        # Validate inputs
        self._validate_introns_have_zscores(u12_reference, "u12_reference")
        self._validate_introns_have_zscores(u2_reference, "u2_reference")
        self._validate_introns_have_zscores(experimental, "experimental")

        print(f"Classification pipeline starting (batch mode)...")
        print(f"  Reference: {len(u12_reference)} U12, {len(u2_reference)} U2")
        print(f"  Experimental: {len(experimental)} introns")
        print(f"  Batch size: {batch_size}")

        # Stages 1-2: Optimize and train (same as regular classify)
        if self.optimize_c:
            print("\n=== Stage 1: Hyperparameter Optimization ===")
            optimizer = SVMOptimizer(
                n_rounds=self.n_optimization_rounds,
                random_state=self.random_state,
                n_jobs=self.cv_processes
            )
            parameters = optimizer.optimize(u12_reference, u2_reference)
        else:
            print("\n=== Stage 1: Using Fixed C Parameter ===")
            parameters = SVMParameters(
                C=self.fixed_c,
                cv_score=0.0,
                round_found=0
            )

        print("\n=== Stage 2: Ensemble Training ===")
        trainer = SVMTrainer(
            n_models=self.n_ensemble_models,
            random_state=self.random_state
        )
        ensemble = trainer.train_ensemble(
            u12_reference,
            u2_reference,
            parameters,
            subsample_u2=self.subsample_u2,
            subsample_ratio=self.subsample_ratio
        )

        # Stage 3: Classify in batches
        print("\n=== Stage 3: Classification (Batch Mode) ===")
        predictor = SVMPredictor(
            threshold=self.classification_threshold,
            n_jobs=self.classification_processes
        )
        classified = predictor.predict_batch(ensemble, experimental, batch_size)

        n_u12 = sum(
            1 for i in classified
            if i.metadata and i.metadata.type_id == 'u12'
        )
        n_u2 = len(classified) - n_u12
        print(f"Classification complete:")
        print(f"  U12: {n_u12} ({100*n_u12/len(classified):.1f}%)")
        print(f"  U2: {n_u2} ({100*n_u2/len(classified):.1f}%)")

        return ClassificationResult(
            classified_introns=classified,
            ensemble=ensemble,
            parameters=parameters,
            n_u12_reference=len(u12_reference),
            n_u2_reference=len(u2_reference)
        )

    def _validate_introns_have_zscores(
        self,
        introns: Sequence[Intron],
        dataset_name: str
    ) -> None:
        """
        Validate that all introns have z-scores.

        CRITICAL: This check ensures no data leakage. Z-scores must be
        computed from reference data BEFORE classification.

        Args:
            introns: Introns to validate
            dataset_name: Name of dataset for error messages

        Raises:
            ValueError: If any intron lacks z-scores
        """
        for i, intron in enumerate(introns):
            if intron.scores is None:
                raise ValueError(
                    f"{dataset_name}[{i}] ({intron.intron_id}): "
                    f"No scores object. Run scoring pipeline first."
                )

            if (intron.scores.five_z_score is None or
                intron.scores.bp_z_score is None or
                intron.scores.three_z_score is None):
                raise ValueError(
                    f"{dataset_name}[{i}] ({intron.intron_id}): "
                    f"Missing z-scores. Must compute z-scores from reference "
                    f"data before classification."
                )
