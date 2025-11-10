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

from typing import Sequence, Optional, Tuple, Any
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
        eval_result: Optional evaluation results (NestedCVResult or SplitEvalResult)
    """

    classified_introns: Sequence[Intron]
    ensemble: SVMEnsemble
    parameters: SVMParameters
    n_u12_reference: int
    n_u2_reference: int
    eval_result: Optional[Any] = None

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
        ...     n_optimization_rounds=3,
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
        n_optimization_rounds: int = 3,
        n_ensemble_models: int = 3,
        classification_threshold: float = 90.0,
        subsample_u2: bool = True,
        subsample_ratio: float = 0.8,
        random_state: int = 42,
        optimize_c: bool = True,
        fixed_c: Optional[float] = None,
        cv_processes: int = 1,
        classification_processes: int = 1,
        max_iter: int = 20000,
        eval_mode: str = 'nested_cv',
        n_cv_folds: int = 5,
        test_fraction: float = 0.2
    ):
        """
        Initialize classifier.

        Args:
            n_optimization_rounds: Number of grid search refinement rounds (default: 3)
            n_ensemble_models: Number of models in ensemble (default: 3)
            classification_threshold: U12 probability threshold 0-100 (default: 90.0)
            subsample_u2: Whether to subsample U2 for ensemble diversity (default: True)
            subsample_ratio: Fraction of U2s per model (default: 0.8)
            random_state: Random seed for reproducibility (default: 42)
            optimize_c: Whether to optimize C parameter (default: True)
            fixed_c: Fixed C value if not optimizing (default: None)
            cv_processes: Number of parallel jobs for cross-validation (default: 1)
            classification_processes: Number of parallel jobs for classification (default: 1)
            max_iter: Maximum iterations for LinearSVC convergence (default: 20000)
            eval_mode: Evaluation mode: 'nested_cv', 'split', or 'none' (default: 'nested_cv')
            n_cv_folds: Number of CV folds for nested CV (default: 5)
            test_fraction: Test set fraction for split mode (default: 0.2)
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
        self.max_iter = max_iter
        self.n_cv_folds = n_cv_folds
        self.test_fraction = test_fraction

        # Auto-skip evaluation when using fixed C
        # Rationale: When C is pre-specified, evaluation metrics aren't useful
        # since we're not comparing different hyperparameters
        if not optimize_c and fixed_c is not None:
            if eval_mode != 'none':
                print(f"Using fixed C={fixed_c:.6e} - automatically skipping evaluation phase")
                print("(Override with --eval-mode if you want to evaluate performance)")
            self.eval_mode = 'none'
        else:
            self.eval_mode = eval_mode

        # Validate parameters
        if not 0 <= classification_threshold <= 100:
            raise ValueError(
                f"classification_threshold must be 0-100, got {classification_threshold}"
            )
        if not optimize_c and fixed_c is None:
            raise ValueError("Must provide fixed_c if optimize_c is False")
        if eval_mode not in ['nested_cv', 'split', 'none']:
            raise ValueError(
                f"eval_mode must be 'nested_cv', 'split', or 'none', got {eval_mode}"
            )

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

        # ====================================================================
        # PHASE 1: EVALUATION (Honest Performance Assessment)
        # ====================================================================
        eval_result = None

        if self.eval_mode == 'nested_cv':
            from classification.nested_cv import NestedCVEvaluator

            evaluator = NestedCVEvaluator(
                n_folds=self.n_cv_folds,
                n_optimization_rounds=self.n_optimization_rounds,
                n_ensemble_models=1,  # Use 1 for speed in CV
                classification_threshold=self.classification_threshold,
                subsample_u2=False,  # Disable for speed in CV
                random_state=self.random_state,
                n_jobs=self.cv_processes,
                max_iter=self.max_iter,
                verbose=True,
                optimize_c=self.optimize_c,
                fixed_c=self.fixed_c
            )
            eval_result = evaluator.evaluate(u12_reference, u2_reference)

        elif self.eval_mode == 'split':
            from classification.split_eval import SplitEvaluator

            evaluator = SplitEvaluator(
                test_fraction=self.test_fraction,
                val_fraction=0.2,  # Fixed validation fraction
                n_optimization_rounds=self.n_optimization_rounds,
                n_ensemble_models=1,  # Use 1 for speed in evaluation
                classification_threshold=self.classification_threshold,
                subsample_u2=False,  # Disable for speed in evaluation
                random_state=self.random_state,
                n_jobs=self.cv_processes,
                max_iter=self.max_iter,
                verbose=True,
                optimize_c=self.optimize_c,
                fixed_c=self.fixed_c
            )
            eval_result = evaluator.evaluate(u12_reference, u2_reference)

        # If eval_mode == 'none', skip evaluation entirely

        # ====================================================================
        # PHASE 2: PRODUCTION MODEL (Train on ALL reference data)
        # ====================================================================
        if self.eval_mode != 'none':
            print("\n" + "="*80)
            print("Production Model Training (all reference data)")
            print("="*80)

        # Stage 1: Optimize hyperparameters
        if self.optimize_c:
            print("\n=== Stage 1: Hyperparameter Optimization ===")
            optimizer = SVMOptimizer(
                n_rounds=self.n_optimization_rounds,
                random_state=self.random_state,
                n_jobs=self.cv_processes,
                max_iter=self.max_iter
            )
            parameters = optimizer.optimize(u12_reference, u2_reference)
            print(f"Optimized C={parameters.C:.6e}, CV score={parameters.cv_score:.4f}")
        else:
            print("\n=== Stage 1: Using Fixed C Parameter ===")
            parameters = SVMParameters(
                C=self.fixed_c,
                calibration_method='sigmoid',  # Default calibration method
                dual=False,  # Primal formulation (n_samples >> n_features: 21k >> 3)
                intercept_scaling=1000.0,  # High value to avoid over-regularizing intercept
                cv_score=0.0,  # Not computed
                round_found=0   # Fixed, not optimized
            )
            print(f"Using fixed C={parameters.C:.6e}, dual={parameters.dual}, intercept_scaling={parameters.intercept_scaling}")

        # Stage 2: Train ensemble
        print("\n=== Stage 2: Ensemble Training ===")
        trainer = SVMTrainer(
            n_models=self.n_ensemble_models,
            random_state=self.random_state,
            max_iter=self.max_iter
        )
        ensemble = trainer.train_ensemble(
            u12_reference,
            u2_reference,
            parameters,
            subsample_u2=self.subsample_u2,
            subsample_ratio=self.subsample_ratio
        )
        print(f"Ensemble trained: {len(ensemble.models)} models")

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
            n_u2_reference=len(u2_reference),
            eval_result=eval_result
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

        # ====================================================================
        # PHASE 1: EVALUATION (Honest Performance Assessment)
        # ====================================================================
        eval_result = None

        if self.eval_mode == 'nested_cv':
            from classification.nested_cv import NestedCVEvaluator

            evaluator = NestedCVEvaluator(
                n_folds=self.n_cv_folds,
                n_optimization_rounds=self.n_optimization_rounds,
                n_ensemble_models=1,  # Use 1 for speed in CV
                classification_threshold=self.classification_threshold,
                subsample_u2=False,  # Disable for speed in CV
                random_state=self.random_state,
                n_jobs=self.cv_processes,
                max_iter=self.max_iter,
                verbose=True,
                optimize_c=self.optimize_c,
                fixed_c=self.fixed_c
            )
            eval_result = evaluator.evaluate(u12_reference, u2_reference)

        elif self.eval_mode == 'split':
            from classification.split_eval import SplitEvaluator

            evaluator = SplitEvaluator(
                test_fraction=self.test_fraction,
                val_fraction=0.2,  # Fixed validation fraction
                n_optimization_rounds=self.n_optimization_rounds,
                n_ensemble_models=1,  # Use 1 for speed in evaluation
                classification_threshold=self.classification_threshold,
                subsample_u2=False,  # Disable for speed in evaluation
                random_state=self.random_state,
                n_jobs=self.cv_processes,
                max_iter=self.max_iter,
                verbose=True,
                optimize_c=self.optimize_c,
                fixed_c=self.fixed_c
            )
            eval_result = evaluator.evaluate(u12_reference, u2_reference)

        # If eval_mode == 'none', skip evaluation entirely

        # ====================================================================
        # PHASE 2: PRODUCTION MODEL (Train on ALL reference data)
        # ====================================================================
        if self.eval_mode != 'none':
            print("\n" + "="*80)
            print("Production Model Training (all reference data)")
            print("="*80)

        # Stage 1: Optimize hyperparameters
        if self.optimize_c:
            print("\n=== Stage 1: Hyperparameter Optimization ===")
            optimizer = SVMOptimizer(
                n_rounds=self.n_optimization_rounds,
                random_state=self.random_state,
                n_jobs=self.cv_processes,
                max_iter=self.max_iter
            )
            parameters = optimizer.optimize(u12_reference, u2_reference)
            print(f"Optimized C={parameters.C:.6e}, CV score={parameters.cv_score:.4f}")
        else:
            print("\n=== Stage 1: Using Fixed C Parameter ===")
            parameters = SVMParameters(
                C=self.fixed_c,
                calibration_method='sigmoid',  # Default calibration method
                dual=False,  # Primal formulation (n_samples >> n_features: 21k >> 3)
                intercept_scaling=1000.0,  # High value to avoid over-regularizing intercept
                cv_score=0.0,  # Not computed
                round_found=0   # Fixed, not optimized
            )
            print(f"Using fixed C={parameters.C:.6e}, dual={parameters.dual}, intercept_scaling={parameters.intercept_scaling}")

        # Stage 2: Train ensemble
        print("\n=== Stage 2: Ensemble Training ===")
        trainer = SVMTrainer(
            n_models=self.n_ensemble_models,
            random_state=self.random_state,
            max_iter=self.max_iter
        )
        ensemble = trainer.train_ensemble(
            u12_reference,
            u2_reference,
            parameters,
            subsample_u2=self.subsample_u2,
            subsample_ratio=self.subsample_ratio
        )
        print(f"Ensemble trained: {len(ensemble.models)} models")

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
            n_u2_reference=len(u2_reference),
            eval_result=eval_result
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
