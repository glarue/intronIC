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

from dataclasses import dataclass, fields, MISSING, replace
from typing import Sequence, Tuple, Optional, Any
import warnings
import contextlib
import numpy as np
from sklearn.svm import LinearSVC, SVC
from sklearn.calibration import CalibratedClassifierCV
from sklearn.base import clone
from sklearn.pipeline import Pipeline
from sklearn.exceptions import ConvergenceWarning

from intronIC.core.intron import Intron
from intronIC.classification.transformers import BothEndsStrongTransformer
from intronIC.classification.svm_factory import create_svm, make_scaler_step
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
    dropped_feature: Optional[int] = None  # Index of feature dropped during training (None = no dropout)
    feature_median: Optional[float] = None  # Median value used to mask the dropped feature


def _svmmodel_getstate(self):
    """Backward-compatible pickle for SVMModel.

    Handles models pickled before dropped_feature/feature_median fields
    were added. Uninitialized slots return their default values.
    """
    state = []
    for f in fields(self):
        try:
            state.append(getattr(self, f.name))
        except AttributeError:
            state.append(f.default if f.default is not MISSING else None)
    return state


def _svmmodel_setstate(self, state):
    """Backward-compatible unpickle for SVMModel."""
    for f, val in zip(fields(self), state):
        object.__setattr__(self, f.name, val)
    # Fill missing fields from models with fewer fields
    for f in fields(self)[len(state):]:
        object.__setattr__(self, f.name,
                           f.default if f.default is not MISSING else None)


# Override the auto-generated pickle methods to handle backward compatibility
SVMModel.__getstate__ = _svmmodel_getstate
SVMModel.__setstate__ = _svmmodel_setstate


@dataclass(frozen=True, slots=True)
class SVMEnsemble:
    """Collection of trained models for ensemble prediction."""

    models: Sequence[SVMModel]
    subsample_ratio: float = 0.85  # Fraction of U2s used per model (for display)

    def __len__(self) -> int:
        return len(self.models)


#: Base motif features (background-corrected log-odds). The z-normalized variant was
#: removed in the supplant 2b refactor; scoring operates on raw motif scores only.
RAW_BASE_FEATURES = ("five_raw_score", "bp_raw_score", "three_raw_score")


def _extract_feature_vector(intron: Intron, extra_names: list,
                            base_names: tuple = RAW_BASE_FEATURES) -> list:
    """Extract feature vector from an intron for training or prediction.

    Base features: the three motif scores named in ``base_names`` (the z triple by
    default; the raw log-odds triple for the ``raw_gated`` scoring mode).
    Extra features: additional IntronScores fields by name (``extra_names``).

    Args:
        intron: Intron with scores populated
        extra_names: List of IntronScores attribute names to append
        base_names: The three base motif attribute names (default: z-scores). The
            raw_gated mode passes :data:`RAW_BASE_FEATURES`.

    Returns:
        List of floats [base1, base2, base3, ...extra]

    Raises:
        ValueError: If intron lacks scores or any base feature is missing.
    """
    if intron.scores is None:
        raise ValueError(f"Intron {intron.intron_id} has no scores")
    base = []
    for name in base_names:
        val = getattr(intron.scores, name, None)
        if val is None:
            raise ValueError(f"Intron {intron.intron_id} missing base feature {name}")
        base.append(float(val))
    for name in extra_names:
        val = getattr(intron.scores, name, None)
        base.append(float(val) if val is not None else 0.0)
    return base


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
        progress_tracker: Optional[Any] = None,
        extra_feature_names: Optional[list] = None,
        feature_dropout: int = 0,
        feature_dropout_fraction: float = 1.0,
        base_features: tuple = RAW_BASE_FEATURES,
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
            extra_feature_names: Additional IntronScores fields to include as classifier features
            feature_dropout: Number of input features to drop per model (default: 0 = disabled).
                When > 0, each model is trained with a different randomly selected feature
                masked to its population median. This creates ensemble diversity that
                produces higher variance for introns dependent on a single feature (FPs)
                than for introns strong across all features (real U12s).
            feature_dropout_fraction: Fraction of models that use dropout (default: 1.0 = all).
                When < 1.0, creates a mixed ensemble: the first n_models*(1-fraction)
                models are full-feature, the rest use dropout. E.g. n_models=50,
                feature_dropout=1, fraction=0.5 → 25 full + 25 dropout.
        """
        self.n_models = n_models
        self.random_state = random_state
        self.kernel = kernel
        self.max_iter = max_iter
        self.features_list = features_list
        self.gamma_imbalance_options = gamma_imbalance_options
        self.progress_tracker = progress_tracker
        self.extra_feature_names = extra_feature_names or []
        self.feature_dropout = feature_dropout
        self.feature_dropout_fraction = feature_dropout_fraction
        self.base_features = tuple(base_features)

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
        # Print training header (no stage number - context-dependent)
        print(f"\n{'='*80}", flush=True)
        if self.n_models == 1:
            print(f"Model Training", flush=True)
        else:
            print(f"Ensemble Training ({self.n_models} models)", flush=True)
        print(f"{'='*80}\n", flush=True)

        # Prepare per-model configurations (subsample + dropout selection)
        # Mixed ensemble: first n_full models are full-feature, rest use dropout
        n_full = self.n_models
        if self.feature_dropout > 0:
            n_full = int(self.n_models * (1 - self.feature_dropout_fraction))

        model_configs = []
        for i in range(self.n_models):
            if subsample_u2 and self.n_models > 1:
                u2_sample = self._subsample_u2(
                    u2_introns,
                    seed=self.random_state + i,
                    subsample_ratio=subsample_ratio
                )
            else:
                u2_sample = u2_introns

            if self.feature_dropout > 0 and i >= n_full:
                rng = np.random.RandomState(self.random_state + i + 1000)
                n_input_features = 3 + len(self.extra_feature_names)
                dropped = int(rng.randint(0, n_input_features))
            else:
                dropped = None

            model_configs.append((i, u2_sample, dropped))

        # Train models — parallel if n_models > 1 and joblib available
        n_parallel = min(self.n_models, 6)  # Cap at 6 parallel workers
        use_parallel = self.n_models > 1 and n_parallel > 1

        if use_parallel:
            print(f"Training {self.n_models} models in parallel "
                  f"({n_parallel} workers)\n", flush=True)
            models = self._train_ensemble_parallel(
                u12_introns, model_configs, parameters, n_parallel
            )
        else:
            models = self._train_ensemble_sequential(
                u12_introns, model_configs, parameters
            )

        print(f"\n{'='*80}", flush=True)
        if self.n_models == 1:
            print(f"Model training complete", flush=True)
        else:
            print(f"Ensemble training complete ({len(models)} models trained)", flush=True)
        print(f"{'='*80}\n", flush=True)

        return SVMEnsemble(models=models, subsample_ratio=subsample_ratio)

    def _train_ensemble_sequential(
        self,
        u12_introns: Sequence[Intron],
        model_configs: list,
        parameters: SVMParameters,
    ) -> list:
        """Train models sequentially with per-model output."""
        models = []
        for i, u2_sample, dropped in model_configs:
            print(f"\n{'─'*80}", flush=True)
            global_step_str = f" {self.progress_tracker.format_step()}" if self.progress_tracker else ""
            if self.n_models == 1:
                print(f"Training model...{global_step_str}", flush=True)
            else:
                print(f"MODEL {i+1}/{self.n_models}: Training ensemble model...{global_step_str}", flush=True)
            print(f"{'─'*80}", flush=True)

            model = self._train_single_model(
                u12_introns, u2_sample, parameters,
                seed=self.random_state + i, drop_feature=dropped,
            )
            models.append(model)

            if self.progress_tracker:
                self.progress_tracker.increment(f"Completed model {i+1}/{self.n_models}")

        return models

    def _train_ensemble_parallel(
        self,
        u12_introns: Sequence[Intron],
        model_configs: list,
        parameters: SVMParameters,
        n_parallel: int,
    ) -> list:
        """Train models in parallel with live progress updates.

        Uses joblib with threading backend. Each completed model updates
        a shared counter and prints a progress line. Per-model sklearn
        output is suppressed to avoid interleaving.

        Falls back to sequential training if joblib is unavailable.
        """
        try:
            from joblib import Parallel, delayed
        except ImportError:
            print("joblib not available, falling back to sequential training",
                  flush=True)
            return self._train_ensemble_sequential(
                u12_introns, model_configs, parameters
            )

        import io
        import sys
        import warnings
        import threading

        # Shared progress counter
        lock = threading.Lock()
        completed = [0]  # mutable for closure access
        n_total = self.n_models

        def _train_with_progress(trainer_self, u12, u2_sample, params,
                                 seed, drop_feature, model_idx):
            """Train one model, suppress output, update progress on completion."""
            old_stdout = sys.stdout
            sys.stdout = io.StringIO()
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    model = trainer_self._train_single_model(
                        u12, u2_sample, params,
                        seed=seed, drop_feature=drop_feature,
                    )
            finally:
                sys.stdout = old_stdout

            # Thread-safe progress update
            with lock:
                completed[0] += 1
                n_done = completed[0]
                drop_info = (f" (dropped feature {model.dropped_feature})"
                             if model.dropped_feature is not None else "")
                old_stdout.write(
                    f"\r  Models trained: {n_done}/{n_total}"
                    f"  [latest: model {model_idx+1}{drop_info}]"
                    + " " * 20
                )
                old_stdout.flush()

            return model

        # Use threads to avoid pickling — GIL is released during sklearn C fits
        raw_results = Parallel(n_jobs=n_parallel, prefer="threads")(
            delayed(_train_with_progress)(
                self, u12_introns, u2_sample, parameters,
                self.random_state + i, dropped, i,
            )
            for i, u2_sample, dropped in model_configs
        )

        # Final newline after progress counter
        print(flush=True)

        # Print completion summary
        for idx, model in enumerate(raw_results):
            if self.progress_tracker:
                self.progress_tracker.increment(f"Completed model {idx+1}/{self.n_models}")

        return list(raw_results)

    def _train_single_model(
        self,
        u12_introns: Sequence[Intron],
        u2_introns: Sequence[Intron],
        parameters: SVMParameters,
        seed: int,
        drop_feature: Optional[int] = None,
    ) -> SVMModel:
        """
        Train a single SVM model.

        Args:
            drop_feature: If not None, mask this feature column to its median
                before training. The model learns to classify without this feature.

        Port from: intronIC.py:5353-5428
        """
        # Prepare data
        X, y = self._prepare_training_data(u12_introns, u2_introns)

        # Apply feature dropout: mask one feature to its population median
        feature_median = None
        if drop_feature is not None and 0 <= drop_feature < X.shape[1]:
            feature_median = float(np.median(X[:, drop_feature]))
            X = X.copy()  # Don't modify shared data
            X[:, drop_feature] = feature_median

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

        # Resolve gamma: 0.0 sentinel means 'scale' was used during optimization
        gamma = parameters.gamma if parameters.gamma != 0.0 else 'scale'

        # In-pipeline GLOBAL scaler resolved from the bundle's stamped choice (default 'standard').
        # 'none' drops the 'scale' step entirely (raw features feed the kernel directly).
        steps = [
            ('transform', BothEndsStrongTransformer(
                features=self.features_list,
                gamma_imbalance=parameters.gamma_imbalance,
                extra_feature_names=self.extra_feature_names,
            )),
        ]
        scaler_step = make_scaler_step(getattr(parameters, 'scaler', 'standard'))
        if scaler_step is not None:
            steps.append(('scale', scaler_step))
        steps.append(('svc', create_svm(
            kernel=parameters.kernel,
            C=parameters.C,
            gamma=gamma,
            class_weight=class_weight,
            penalty=parameters.penalty,
            loss=parameters.loss,
            max_iter=self.max_iter,
            random_state=seed,
        )))
        base_pipeline = Pipeline(steps)

        # External calibration wrapper
        # Method (sigmoid vs isotonic) chosen by hyperparameter optimization
        # ensemble='auto' uses per-fold fit + averaging for stable tails
        # n_jobs: -1 when training sequentially (parallelize CV folds),
        #         1 when training in parallel (avoid nested parallelism)
        cal_n_jobs = 1 if self.n_models > 1 else -1
        svm = CalibratedClassifierCV(
            base_pipeline,
            method=parameters.calibration_method,  # From optimizer: 'sigmoid' or 'isotonic'
            cv=5,  # Stratified 5-fold
            ensemble='auto',  # Per-fold calibrators averaged
            n_jobs=cal_n_jobs,
        )

        # Suppress convergence warning spam but log summary
        with suppress_convergence_warnings(verbose=True):
            svm.fit(X, y)

        return SVMModel(
            model=svm,
            train_size=len(X),
            u12_count=int(np.sum(y == 1)),
            u2_count=int(np.sum(y == 0)),
            # Record the model's feature space so the predictor reconstructs it (raw vs z).
            parameters=replace(parameters, base_features=tuple(self.base_features)),
            dropped_feature=drop_feature,
            feature_median=feature_median,
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

        Features: [five_z_score, bp_z_score, three_z_score, ...extra_features]
        Labels: 0 = U2, 1 = U12

        Args:
            u12_introns: U12-type introns (normalized)
            u2_introns: U2-type introns (normalized)

        Returns:
            X: Feature matrix (n_samples, 3 + len(extra_features))
            y: Labels (n_samples,)

        Raises:
            ValueError: If any intron lacks z-scores
        """
        # Extract U12 features
        u12_features = []
        for intron in u12_introns:
            u12_features.append(_extract_feature_vector(intron, self.extra_feature_names,
                                                        base_names=self.base_features))

        # Extract U2 features
        u2_features = []
        for intron in u2_introns:
            u2_features.append(_extract_feature_vector(intron, self.extra_feature_names,
                                                       base_names=self.base_features))

        # Combine features and create labels
        # Port from: intronIC.py:5380, 5397
        X = np.array(u2_features + u12_features)
        y = np.array([0] * len(u2_features) + [1] * len(u12_features))

        return X, y
