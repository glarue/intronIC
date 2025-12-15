"""
Model save/load functions for intronIC.

Handles serialization of trained models with metadata.
"""

import json
import joblib
from pathlib import Path
from typing import Dict, Any, Optional
from datetime import datetime


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
