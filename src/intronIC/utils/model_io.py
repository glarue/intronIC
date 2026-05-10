"""
Model save/load functions for intronIC.

Handles serialization of trained models with metadata.
"""

import json
import joblib
from pathlib import Path
from typing import Dict, Any, Optional, Tuple
from datetime import datetime


# ── v3 multispecies bundle support ────────────────────────────────────
# v3 bundles ship the model as a single self-describing dict that pairs
# the trained sub-models with the config + training metadata used to
# produce them, instead of the legacy v2.3 flat-dict-of-handles. The
# loader normalizes both formats into the same runtime shape so the rest
# of the pipeline stays format-agnostic.
V3_VERSION = "v3"
V3_DEFAULT_THRESHOLD = 50.0


def _build_v3_ensemble(bundle: Dict[str, Any]) -> Tuple[Any, Tuple[str, ...]]:
    """Construct a runtime SVMEnsemble from a v3 bundle's seeds dict.

    Returns (ensemble, extra_features). Imports of SVMEnsemble/SVMModel/
    SVMParameters are local because importing trainer/optimizer at module
    scope would create a cycle through scoring/cluster_validation.
    """
    from intronIC.classification.optimizer import SVMParameters
    from intronIC.classification.trainer import SVMEnsemble, SVMModel

    config = bundle.get("config", {})
    input_features = list(config.get(
        "input_features",
        ["five_z_score", "bp_z_score", "three_z_score"],
    ))
    if len(input_features) < 3:
        raise ValueError(
            f"v3 bundle config.input_features must list at least the three "
            f"base z-scores; got {input_features!r}"
        )
    extra_features = tuple(input_features[3:])

    # Synthesize SVMParameters that describe the trained models.
    # Kernel/C/gamma are recorded for diagnostics; production reads
    # extra_features and include_max from this struct.
    params = SVMParameters(
        C=float(config.get("C", 200.0)),
        calibration_method=str(config.get("calibration_method", "isotonic")),
        saturate_enabled=False,
        include_max=False,
        include_pairwise_mins=False,
        penalty="l2",
        class_weight_multiplier=1.0,
        loss="squared_hinge",
        gamma_imbalance=1.0,
        kernel=str(config.get("kernel", "rbf")),
        gamma=float(config.get("gamma", 0.0)),
        extra_features=extra_features,
    )

    training = bundle.get("training", {})
    n_train = int(training.get("n_train", 0))
    u12_count = int(training.get("u12_positives", 0))
    u2_count = max(0, n_train - u12_count)

    sub_models = []
    for seed_id in config.get("seeds", []):
        seed_block = bundle["seeds"][seed_id]
        for sub in seed_block["models"]:
            sub_models.append(SVMModel(
                model=sub,
                train_size=n_train,
                u12_count=u12_count,
                u2_count=u2_count,
                parameters=params,
                dropped_feature=None,
                feature_median=None,
            ))

    if not sub_models:
        raise ValueError("v3 bundle contains no sub-models in config.seeds")

    ensemble = SVMEnsemble(
        models=tuple(sub_models),
        subsample_ratio=float(config.get("easy_fraction", 0.85)),
    )
    return ensemble, extra_features


def _v3_to_runtime(bundle: Dict[str, Any]) -> Dict[str, Any]:
    """Translate a v3 bundle into the v2.3-style runtime dict.

    Preserves the original v3 metadata under a private key so downstream
    code can still introspect model_id / training stats if needed.

    Adaptive vs fallback normalizer:
      - "normalizer": absent (None) by design — v3 uses adaptive RobustScaler
        fitting on the experimental data, fit on a per-contig pre-pass.
      - "fallback_normalizer": optional. When adaptive fit can't be done
        (input too small, e.g. single-intron scoring), the streaming
        classify path falls through to this pre-fit scaler. Bundled as
        v3_fallback_normalizer.pkl sibling and loaded here if present.
    """
    ensemble, _ = _build_v3_ensemble(bundle)
    training = bundle.get("training", {})
    n_train = int(training.get("n_train", 0))
    u12_total = int(training.get("u12_positives", 0))
    # u12_positives is the corpus-wide count; subtract held-out positives
    # so the prior reflects what the SVM actually trained against.
    u12_train = u12_total - int(training.get("test_positives", 0))
    if n_train > 0 and 0 < u12_train <= n_train:
        training_prior = u12_train / n_train
    else:
        training_prior = 0.5

    # Look for a bundled fallback normalizer next to the model file. Kept
    # as a sibling .pkl (rather than embedded in the bundle dict) so it
    # can be updated independently without re-pickling the SVM ensemble.
    fallback_normalizer = _load_fallback_normalizer(bundle)

    return {
        "ensemble": ensemble,
        "normalizer": None,            # v3 ships pre-z-scored features; adaptive at runtime
        "fallback_normalizer": fallback_normalizer,
        "threshold": V3_DEFAULT_THRESHOLD,
        "training_prior": training_prior,
        "human_negative_stats": None,  # not used by current production path
        # Provenance — preserved so downstream can log/inspect; ignored otherwise.
        "_v3_bundle": bundle,
    }


def _load_fallback_normalizer(bundle: Dict[str, Any]) -> Any:
    """Load v3_fallback_normalizer.pkl from the same directory as the bundle.

    Returns the unpickled ScoreNormalizer instance, or None if the file
    isn't present or can't be loaded. Bundle path is recovered from
    bundle's ``_source_path`` if set; otherwise checks the package data dir.
    """
    candidates = []
    src = bundle.get("_source_path")
    if src:
        candidates.append(Path(src).parent / "v3_fallback_normalizer.pkl")
    # Fall back to the package's bundled data dir
    candidates.append(
        Path(__file__).parent.parent / "data" / "v3_fallback_normalizer.pkl"
    )
    for path in candidates:
        if path.exists():
            try:
                return joblib.load(path)
            except Exception:
                continue
    return None


def normalize_model_bundle(model_data: Any) -> Any:
    """Coerce a loaded model object into the runtime shape callers expect.

    Pass-through for v2.3-style dicts and bare SVMEnsemble objects.
    Translates v3 multispecies bundles ({"version": "v3", ...}) into
    the runtime dict shape with `ensemble`, `normalizer`, `threshold`,
    `training_prior`.
    """
    if isinstance(model_data, dict) and model_data.get("version") == V3_VERSION:
        return _v3_to_runtime(model_data)
    return model_data


def save_model(
    model: Any,
    output_dir: Path,
    base_filename: str,
    species_name: str,
    metrics: Dict[str, Any],
    config: Optional[Dict[str, Any]] = None
) -> Path:
    """Save trained model to disk with metadata.

    Args:
        model: Trained classifier (or ensemble)
        output_dir: Directory to save model
        base_filename: Base name for output files
        species_name: Species name
        metrics: Training metrics (F1, precision, recall, etc.)
        config: Optional training configuration

    Returns:
        Path to saved model file

    Saves:
        - {base_filename}.model.pkl - Pickled model
        - {base_filename}.model.metadata.json - Training metadata
    """
    # Save model
    model_path = output_dir / f"{base_filename}.model.pkl"
    joblib.dump(model, model_path)

    # Create metadata
    metadata = {
        'species_name': species_name,
        'trained_date': datetime.now().isoformat(),
        'intronIC_version': '2.0.0',
        'pipeline_architecture': 'single_scaler_v2025',  # NEW: Track architecture version
        'metrics': metrics,
    }

    # Add config if provided
    if config:
        metadata['training_config'] = config

    # Save metadata
    metadata_path = output_dir / f"{base_filename}.model.metadata.json"
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)

    return model_path


def load_model(model_path: Path) -> Any:
    """Load trained model from disk.

    Args:
        model_path: Path to .model.pkl file

    Returns:
        Loaded model

    Raises:
        FileNotFoundError: If model file not found
        Exception: If model loading fails
    """
    if not model_path.exists():
        raise FileNotFoundError(f"Model file not found: {model_path}")

    try:
        # Suppress sklearn version warnings - LinearSVC is stable across minor versions
        import warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*sklearn.*")
            model = joblib.load(model_path)
        return model
    except Exception as e:
        raise Exception(f"Failed to load model from {model_path}: {str(e)}")


def load_model_metadata(model_path: Path) -> Optional[Dict[str, Any]]:
    """Load model metadata if available.

    Args:
        model_path: Path to .model.pkl file

    Returns:
        Metadata dictionary or None if not found
    """
    metadata_path = model_path.with_suffix('.metadata.json')
    if metadata_path.exists():
        try:
            # Try to load as plain JSON first
            with open(metadata_path, 'r') as f:
                return json.load(f)
        except UnicodeDecodeError:
            # File might be compressed - try with gzip
            import gzip
            try:
                with gzip.open(metadata_path, 'rt') as f:
                    return json.load(f)
            except Exception as e:
                # If both fail, return None and log warning
                import warnings
                warnings.warn(f"Failed to load metadata from {metadata_path}: {e}")
                return None
    return None
