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

from collections import Counter
from dataclasses import dataclass
from typing import Any, Optional, Sequence, Tuple

import numpy as np

from intronIC.classification.optimizer import SVMOptimizer, SVMParameters
from intronIC.classification.predictor import SVMPredictor
from intronIC.classification.trainer import SVMEnsemble, SVMTrainer
from intronIC.core.intron import Intron


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
            i
            for i in self.classified_introns
            if i.metadata
            and i.metadata.type_id == "u12"
            and i.scores
            and i.scores.svm_score
            and i.scores.svm_score >= threshold
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
            i
            for i in self.classified_introns
            if i.metadata
            and i.metadata.type_id == "u2"
            and i.scores
            and i.scores.svm_score
            and i.scores.svm_score < threshold
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
        eval_mode: str = "nested_cv",
        n_cv_folds: int = 5,
        test_fraction: float = 0.2,
        param_grid_override: Optional[dict] = None,
        n_points_initial: int = 13,
        eff_C_pos_range: tuple = (3e-4, 1e2),  # Tightened to reduce FPR
        eff_C_neg_max: Optional[float] = None,
        use_fold_averaged_params: bool = False,
        scoring_metric: str = "balanced_accuracy",
        penalty_options: Optional[list] = None,
        loss_options: Optional[list] = None,
        class_weight_multipliers: Optional[list] = None,
        use_multiplier_tiebreaker: bool = True,
        features_list: Optional[list] = None,
        gamma_imbalance_options: Optional[list] = None,
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
            param_grid_override: Optional custom parameter grid for optimizer (default: None)
            n_points_initial: Initial grid points for round 1 optimization (default: 13)
            eff_C_pos_range: Target effective penalty range for positive class
                           (default: 3e-4 to 1e+2, tightened to reduce FPR)
                           Expert: "larger C tends to raise FPR"
            eff_C_neg_max: Optional cap on negative class effective penalty (default: None)
            use_fold_averaged_params: Use fold-averaged hyperparameters from nested CV instead of
                                     re-optimizing on full dataset (default: False). When True,
                                     uses geometric mean of fold C values and majority-vote calibration.
                                     Recommended for cross-species applications. Only applies when
                                     eval_mode='nested_cv'.
            scoring_metric: Metric for hyperparameter optimization (default: 'balanced_accuracy').
                          Options: 'balanced_accuracy', 'f_beta_0.5', 'f_beta_0.75'.
                          Use 'f_beta_0.5' for precision-focused training (minimizes false positives).
            penalty_options: Penalty types to search (default: ['l2']). Options: ['l1', 'l2'].
                           WARNING: L1 can be 10-20× slower than L2.
            loss_options: Loss functions to search (default: ['squared_hinge']). Options: ['hinge', 'squared_hinge'].
                         NOTE: With dual=False, only 'squared_hinge' is valid.
            class_weight_multipliers: Class weight multipliers for precision/recall tradeoffs (default: [1.0]).
                                     Options: [0.8, 1.0, 1.2] - lower values are more conservative.
            use_multiplier_tiebreaker: When multiple multipliers are tied for best score, prefer 1.0 (default: True).
                                      Set to False to use GridSearchCV's default first-tie behavior.
            features_list: List of composite feature names to include in BothEndsStrongTransformer (default: None).
                          If None, uses default 4D feature set (base 3 z-scores + absdiff_bp_3).
                          Available features: 'min_5_bp', 'min_5_3', 'min_all',
                                            'absdiff_5_bp', 'absdiff_5_3', 'absdiff_bp_3',
                                            'max_5_bp', 'max_5_3'
                          Example: ['absdiff_bp_3'] for minimal 4D space (default, based on L1 analysis)
            gamma_imbalance_options: List of gamma scaling factors to grid search (default: None).
                                    Scales all neg_absdiff_* features to nudge L2 toward L1 behavior.
                                    If None, uses gamma=1.0 (no scaling).
                                    Higher values make imbalance more costly in the margin.
                                    Example: [1.0, 2.0, 4.0] to search over scaling factors
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
        self.param_grid_override = param_grid_override
        self.n_points_initial = n_points_initial
        self.eff_C_pos_range = eff_C_pos_range
        self.eff_C_neg_max = eff_C_neg_max
        self.use_fold_averaged_params = use_fold_averaged_params
        self.scoring_metric = scoring_metric
        self.penalty_options = (
            penalty_options if penalty_options is not None else ["l2"]
        )
        self.loss_options = (
            loss_options if loss_options is not None else ["squared_hinge"]
        )
        self.class_weight_multipliers = (
            class_weight_multipliers if class_weight_multipliers is not None else [1.0]
        )
        self.use_multiplier_tiebreaker = use_multiplier_tiebreaker
        self.features_list = features_list
        self.gamma_imbalance_options = (
            gamma_imbalance_options if gamma_imbalance_options is not None else [1.0]
        )

        # Debug: Log features being used
        if features_list is not None:
            print(
                f"IntronClassifier initialized with explicit feature list: {features_list}"
            )
        else:
            print(f"IntronClassifier initialized with default 4D feature set")

        # Auto-skip evaluation when using fixed C
        # Rationale: When C is pre-specified, evaluation metrics aren't useful
        # since we're not comparing different hyperparameters
        if not optimize_c and fixed_c is not None:
            if eval_mode != "none":
                print(
                    f"Using fixed C={fixed_c:.6e} - automatically skipping evaluation phase"
                )
                print("(Override with --eval-mode if you want to evaluate performance)")
            self.eval_mode = "none"
        else:
            self.eval_mode = eval_mode

        # Validate parameters
        if not 0 <= classification_threshold <= 100:
            raise ValueError(
                f"classification_threshold must be 0-100, got {classification_threshold}"
            )
        if not optimize_c and fixed_c is None:
            raise ValueError("Must provide fixed_c if optimize_c is False")
        if eval_mode not in ["nested_cv", "split", "none"]:
            raise ValueError(
                f"eval_mode must be 'nested_cv', 'split', or 'none', got {eval_mode}"
            )

    def _compute_fold_averaged_params(self, nested_cv_result) -> SVMParameters:
        """
        Compute fold-averaged hyperparameters from nested CV results.

        Uses geometric mean for C (better for log-scale parameters) and
        majority vote for calibration method.

        Args:
            nested_cv_result: NestedCVResult with fold-specific parameters

        Returns:
            SVMParameters with fold-averaged C and calibration_method

        Note:
            This is more conservative than re-optimizing on full dataset,
            favoring generalization over training-set fit. Recommended for
            cross-species applications.
        """
        from intronIC.classification.nested_cv import NestedCVResult

        fold_results = nested_cv_result.fold_results

        # Extract fold-specific C values and calibration methods
        c_values = np.array([fold.optimized_C for fold in fold_results])
        calibration_methods = [fold.calibration_method for fold in fold_results]

        # Geometric mean of C values (better for log-scale parameters)
        geometric_mean_C = float(np.exp(np.mean(np.log(c_values))))

        # Majority vote for calibration method
        method_counts = Counter(calibration_methods)
        majority_calibration = method_counts.most_common(1)[0][0]

        # Log the fold-specific values for transparency
        print(f"\nFold-averaged hyperparameters:")
        print(f"  Fold-specific C values: {[f'{c:.2e}' for c in c_values]}")
        print(f"  Geometric mean C: {geometric_mean_C:.6e}")
        print(f"  Fold-specific calibration: {calibration_methods}")
        print(f"  Majority-vote calibration: {majority_calibration}")
        print(
            f"  Rationale: Using conservative fold-averaged params for better cross-species generalization"
        )

        # Get fixed parameters from first fold's result (these don't vary across folds)
        # Note: FoldResult doesn't store penalty/loss/class_weight_multiplier yet
        # TODO: Enhance FoldResult to store all optimized parameters for better averaging
        # For now, use defaults that match the configured search space
        return SVMParameters(
            C=geometric_mean_C,
            calibration_method=majority_calibration,
            saturate_enabled=False,  # Corrected arch: no saturation
            include_max=False,  # Corrected arch: no max features
            include_pairwise_mins=False,  # Corrected arch: no pairwise mins
            penalty="l2",  # Default: L2 (most common in search space)
            class_weight_multiplier=1.0,  # Default: balanced (middle of search range)
            loss="squared_hinge",  # Fixed: only valid option for dual=False
            dual=False,  # Corrected arch: primal formulation
            intercept_scaling=1.0,  # Corrected arch: fixed
            cv_score=nested_cv_result.mean_f1,  # Use mean F1 from nested CV
            round_found=-1,  # -1 indicates fold-averaged (not from specific round)
        )

    def classify(
        self,
        u12_reference: Sequence[Intron],
        u2_reference: Sequence[Intron],
        experimental: Sequence[Intron],
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

        # Initialize global progress tracker
        from intronIC.classification.progress_tracker import ProgressTracker

        skip_final_opt = self.use_fold_averaged_params and self.eval_mode == "nested_cv"
        total_steps = ProgressTracker.calculate_total_steps(
            eval_mode=self.eval_mode,
            n_cv_folds=self.n_cv_folds,
            n_optimization_rounds=self.n_optimization_rounds,
            n_ensemble_models=self.n_ensemble_models,
            skip_final_optimization=skip_final_opt,
        )
        progress_tracker = ProgressTracker(total_steps=total_steps, verbose=True)
        print(f"\n[Starting training pipeline]\n")

        # ====================================================================
        # PHASE 1: EVALUATION (Honest Performance Assessment)
        # ====================================================================
        eval_result = None

        if self.eval_mode == "nested_cv":
            from intronIC.classification.nested_cv import NestedCVEvaluator

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
                fixed_c=self.fixed_c,
                cv_folds=self.n_cv_folds,
                n_points_initial=self.n_points_initial,
                scoring_metric=self.scoring_metric,
                penalty_options=self.penalty_options,
                loss_options=self.loss_options,
                class_weight_multipliers=self.class_weight_multipliers,
                use_multiplier_tiebreaker=self.use_multiplier_tiebreaker,
                param_grid_override=self.param_grid_override,
                eff_C_pos_range=self.eff_C_pos_range,
                eff_C_neg_max=self.eff_C_neg_max,
                progress_tracker=progress_tracker,
                features_list=self.features_list,
            )
            eval_result = evaluator.evaluate(u12_reference, u2_reference)

        elif self.eval_mode == "split":
            from intronIC.classification.split_eval import SplitEvaluator

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
                fixed_c=self.fixed_c,
                cv_folds=self.n_cv_folds,
                n_points_initial=self.n_points_initial,
                scoring_metric=self.scoring_metric,
                penalty_options=self.penalty_options,
                loss_options=self.loss_options,
                class_weight_multipliers=self.class_weight_multipliers,
                use_multiplier_tiebreaker=self.use_multiplier_tiebreaker,
                param_grid_override=self.param_grid_override,
                eff_C_pos_range=self.eff_C_pos_range,
                eff_C_neg_max=self.eff_C_neg_max,
                progress_tracker=progress_tracker,
                features_list=self.features_list,
            )
            eval_result = evaluator.evaluate(u12_reference, u2_reference)

        # If eval_mode == 'none', skip evaluation entirely

        # ====================================================================
        # PHASE 2: PRODUCTION MODEL (Train on ALL reference data)
        # ====================================================================
        if self.eval_mode != "none":
            print("\n" + "=" * 80)
            print("Production Model Training (all reference data)")
            print("=" * 80)

        # Stage 1: Hyperparameter Optimization
        # Choose between fold-averaged params (conservative) or re-optimization (aggressive)
        print("\n=== Stage 1: Hyperparameter Optimization ===")

        # Check if we should use fold-averaged parameters from nested CV
        if (
            self.use_fold_averaged_params
            and self.eval_mode == "nested_cv"
            and eval_result is not None
        ):
            # Use fold-averaged hyperparameters (better cross-species generalization)
            print("Using fold-averaged hyperparameters from nested CV")
            print(
                "(Skipping re-optimization on full dataset for better cross-species generalization)"
            )
            parameters = self._compute_fold_averaged_params(eval_result)
        else:
            # Standard approach: re-optimize on full dataset
            if self.use_fold_averaged_params and self.eval_mode == "nested_cv":
                print(
                    "Note: use_fold_averaged_params=True but nested CV result not available"
                )
                print("Falling back to re-optimization on full dataset")
            elif self.use_fold_averaged_params:
                print(
                    "Note: use_fold_averaged_params=True but eval_mode is not 'nested_cv'"
                )
                print("Falling back to re-optimization on full dataset")

            # If C is fixed, constrain the grid to that single value
            param_grid = (
                self.param_grid_override.copy() if self.param_grid_override else {}
            )
            if not self.optimize_c:
                # Force C to the fixed value, but still search other parameters
                param_grid["svc__C"] = [self.fixed_c]
                print(
                    f"C fixed at {self.fixed_c:.6e}, optimizing include_max/dual/intercept_scaling/calibration_method"
                )
            else:
                print(
                    "Optimizing C, include_max, dual, intercept_scaling, and calibration_method"
                )

            optimizer = SVMOptimizer(
                n_rounds=self.n_optimization_rounds,
                n_points_initial=self.n_points_initial,
                cv_folds=self.n_cv_folds,
                random_state=self.random_state,
                n_jobs=self.cv_processes,
                max_iter=self.max_iter,
                scoring_metric=self.scoring_metric,
                penalty_options=self.penalty_options,
                loss_options=self.loss_options,
                class_weight_multipliers=self.class_weight_multipliers,
                features_list=self.features_list,
                gamma_imbalance_options=self.gamma_imbalance_options,
                param_grid_override=param_grid if param_grid else None,
                progress_tracker=progress_tracker,
            )
            parameters = optimizer.optimize(
                u12_reference,
                u2_reference,
                eff_C_pos_range=self.eff_C_pos_range,
                eff_C_neg_max=self.eff_C_neg_max,
            )
            print(
                f"Best parameters: C={parameters.C:.6e}, include_max={parameters.include_max}, CV score={parameters.cv_score:.4f}"
            )

        # Stage 2: Train ensemble
        print("\n=== Stage 2: Ensemble Training ===")
        trainer = SVMTrainer(
            n_models=self.n_ensemble_models,
            random_state=self.random_state,
            max_iter=self.max_iter,
            features_list=self.features_list,
            gamma_imbalance_options=self.gamma_imbalance_options,
            progress_tracker=progress_tracker,
        )
        ensemble = trainer.train_ensemble(
            u12_reference,
            u2_reference,
            parameters,
            subsample_u2=self.subsample_u2,
            subsample_ratio=self.subsample_ratio,
        )
        print(f"Ensemble trained: {len(ensemble.models)} models")

        # Sanity check: Verify imbalance penalty weights are negative
        from intronIC.classification.model_inspector import inspect_ensemble_weights

        warnings = inspect_ensemble_weights(ensemble, verbose=True)
        if warnings:
            print("\n⚠️  WARNING: Model coefficients may not be working as expected!")
            print("Review the coefficient analysis above for details.")

        # Stage 3: Classify experimental introns
        print("\n=== Stage 3: Classification ===")

        if len(experimental) == 0:
            print("No experimental introns to classify (training-only mode)")
            classified = []
        else:
            predictor = SVMPredictor(
                threshold=self.classification_threshold,
                n_jobs=self.classification_processes,
            )
            classified = predictor.predict(ensemble, experimental)

            # Count classifications
            n_u12 = sum(
                1 for i in classified if i.metadata and i.metadata.type_id == "u12"
            )
            n_u2 = len(classified) - n_u12
            print(f"Classification complete:")
            print(f"  U12: {n_u12} ({100 * n_u12 / len(classified):.1f}%)")
            print(f"  U2: {n_u2} ({100 * n_u2 / len(classified):.1f}%)")

        return ClassificationResult(
            classified_introns=classified,
            ensemble=ensemble,
            parameters=parameters,
            n_u12_reference=len(u12_reference),
            n_u2_reference=len(u2_reference),
            eval_result=eval_result,
        )

    def classify_batch(
        self,
        u12_reference: Sequence[Intron],
        u2_reference: Sequence[Intron],
        experimental: Sequence[Intron],
        batch_size: int = 10000,
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

        # Initialize global progress tracker
        from intronIC.classification.progress_tracker import ProgressTracker

        skip_final_opt = self.use_fold_averaged_params and self.eval_mode == "nested_cv"
        total_steps = ProgressTracker.calculate_total_steps(
            eval_mode=self.eval_mode,
            n_cv_folds=self.n_cv_folds,
            n_optimization_rounds=self.n_optimization_rounds,
            n_ensemble_models=self.n_ensemble_models,
            skip_final_optimization=skip_final_opt,
        )
        progress_tracker = ProgressTracker(total_steps=total_steps, verbose=True)
        print(f"\n[Starting training pipeline (batch mode)]\n")

        # ====================================================================
        # PHASE 1: EVALUATION (Honest Performance Assessment)
        # ====================================================================
        eval_result = None

        if self.eval_mode == "nested_cv":
            from intronIC.classification.nested_cv import NestedCVEvaluator

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
                fixed_c=self.fixed_c,
                cv_folds=self.n_cv_folds,
                n_points_initial=self.n_points_initial,
                scoring_metric=self.scoring_metric,
                penalty_options=self.penalty_options,
                loss_options=self.loss_options,
                class_weight_multipliers=self.class_weight_multipliers,
                use_multiplier_tiebreaker=self.use_multiplier_tiebreaker,
                param_grid_override=self.param_grid_override,
                eff_C_pos_range=self.eff_C_pos_range,
                eff_C_neg_max=self.eff_C_neg_max,
                progress_tracker=progress_tracker,
                features_list=self.features_list,
            )
            eval_result = evaluator.evaluate(u12_reference, u2_reference)

        elif self.eval_mode == "split":
            from intronIC.classification.split_eval import SplitEvaluator

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
                fixed_c=self.fixed_c,
                cv_folds=self.n_cv_folds,
                n_points_initial=self.n_points_initial,
                scoring_metric=self.scoring_metric,
                penalty_options=self.penalty_options,
                loss_options=self.loss_options,
                class_weight_multipliers=self.class_weight_multipliers,
                use_multiplier_tiebreaker=self.use_multiplier_tiebreaker,
                param_grid_override=self.param_grid_override,
                eff_C_pos_range=self.eff_C_pos_range,
                eff_C_neg_max=self.eff_C_neg_max,
                progress_tracker=progress_tracker,
                features_list=self.features_list,
            )
            eval_result = evaluator.evaluate(u12_reference, u2_reference)

        # If eval_mode == 'none', skip evaluation entirely

        # ====================================================================
        # PHASE 2: PRODUCTION MODEL (Train on ALL reference data)
        # ====================================================================
        if self.eval_mode != "none":
            print("\n" + "=" * 80)
            print("Production Model Training (all reference data)")
            print("=" * 80)

        # Stage 1: Hyperparameter Optimization
        # Choose between fold-averaged params (conservative) or re-optimization (aggressive)
        print("\n=== Stage 1: Hyperparameter Optimization ===")

        # Check if we should use fold-averaged parameters from nested CV
        if (
            self.use_fold_averaged_params
            and self.eval_mode == "nested_cv"
            and eval_result is not None
        ):
            # Use fold-averaged hyperparameters (better cross-species generalization)
            print("Using fold-averaged hyperparameters from nested CV")
            print(
                "(Skipping re-optimization on full dataset for better cross-species generalization)"
            )
            parameters = self._compute_fold_averaged_params(eval_result)
        else:
            # Standard approach: re-optimize on full dataset
            if self.use_fold_averaged_params and self.eval_mode == "nested_cv":
                print(
                    "Note: use_fold_averaged_params=True but nested CV result not available"
                )
                print("Falling back to re-optimization on full dataset")
            elif self.use_fold_averaged_params:
                print(
                    "Note: use_fold_averaged_params=True but eval_mode is not 'nested_cv'"
                )
                print("Falling back to re-optimization on full dataset")

            # If C is fixed, constrain the grid to that single value
            param_grid = (
                self.param_grid_override.copy() if self.param_grid_override else {}
            )
            if not self.optimize_c:
                # Force C to the fixed value, but still search other parameters
                param_grid["svc__C"] = [self.fixed_c]
                print(
                    f"C fixed at {self.fixed_c:.6e}, optimizing include_max/dual/intercept_scaling/calibration_method"
                )
            else:
                print(
                    "Optimizing C, include_max, dual, intercept_scaling, and calibration_method"
                )

            optimizer = SVMOptimizer(
                n_rounds=self.n_optimization_rounds,
                n_points_initial=self.n_points_initial,
                cv_folds=self.n_cv_folds,
                random_state=self.random_state,
                n_jobs=self.cv_processes,
                max_iter=self.max_iter,
                scoring_metric=self.scoring_metric,
                penalty_options=self.penalty_options,
                loss_options=self.loss_options,
                class_weight_multipliers=self.class_weight_multipliers,
                features_list=self.features_list,
                gamma_imbalance_options=self.gamma_imbalance_options,
                param_grid_override=param_grid if param_grid else None,
                progress_tracker=progress_tracker,
            )
            parameters = optimizer.optimize(
                u12_reference,
                u2_reference,
                eff_C_pos_range=self.eff_C_pos_range,
                eff_C_neg_max=self.eff_C_neg_max,
            )
            print(
                f"Best parameters: C={parameters.C:.6e}, include_max={parameters.include_max}, CV score={parameters.cv_score:.4f}"
            )

        # Stage 2: Train ensemble
        print("\n=== Stage 2: Ensemble Training ===")
        trainer = SVMTrainer(
            n_models=self.n_ensemble_models,
            random_state=self.random_state,
            max_iter=self.max_iter,
            features_list=self.features_list,
            gamma_imbalance_options=self.gamma_imbalance_options,
            progress_tracker=progress_tracker,
        )
        ensemble = trainer.train_ensemble(
            u12_reference,
            u2_reference,
            parameters,
            subsample_u2=self.subsample_u2,
            subsample_ratio=self.subsample_ratio,
        )
        print(f"Ensemble trained: {len(ensemble.models)} models")

        # Sanity check: Verify imbalance penalty weights are negative
        from intronIC.classification.model_inspector import inspect_ensemble_weights

        warnings = inspect_ensemble_weights(ensemble, verbose=True)
        if warnings:
            print("\n⚠️  WARNING: Model coefficients may not be working as expected!")
            print("Review the coefficient analysis above for details.")

        # Stage 3: Classify in batches
        print("\n=== Stage 3: Classification (Batch Mode) ===")

        if len(experimental) == 0:
            print("No experimental introns to classify (training-only mode)")
            classified = []
        else:
            predictor = SVMPredictor(
                threshold=self.classification_threshold,
                n_jobs=self.classification_processes,
            )
            classified = predictor.predict_batch(ensemble, experimental, batch_size)

            n_u12 = sum(
                1 for i in classified if i.metadata and i.metadata.type_id == "u12"
            )
            n_u2 = len(classified) - n_u12
            print(f"Classification complete:")
            print(f"  U12: {n_u12} ({100 * n_u12 / len(classified):.1f}%)")
            print(f"  U2: {n_u2} ({100 * n_u2 / len(classified):.1f}%)")

        return ClassificationResult(
            classified_introns=classified,
            ensemble=ensemble,
            parameters=parameters,
            n_u12_reference=len(u12_reference),
            n_u2_reference=len(u2_reference),
            eval_result=eval_result,
        )

    def _validate_introns_have_zscores(
        self, introns: Sequence[Intron], dataset_name: str
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

            if (
                intron.scores.five_z_score is None
                or intron.scores.bp_z_score is None
                or intron.scores.three_z_score is None
            ):
                raise ValueError(
                    f"{dataset_name}[{i}] ({intron.intron_id}): "
                    f"Missing z-scores. Must compute z-scores from reference "
                    f"data before classification."
                )
