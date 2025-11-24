"""
Model metadata generation for intronIC.

Generates comprehensive metadata files documenting model training,
including dataset information, hyperparameters, and evaluation metrics.

Author: intronIC refactoring project
Date: 2025-11-12
"""

import hashlib
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, Optional, Sequence
from intronIC.core.intron import Intron


def compute_file_md5(file_path: Path, chunk_size: int = 8192) -> str:
    """
    Compute MD5 hash of a file.

    Args:
        file_path: Path to file
        chunk_size: Size of chunks to read (default: 8KB)

    Returns:
        MD5 hash as hexadecimal string

    Examples:
        >>> from pathlib import Path
        >>> # compute_file_md5(Path("data/reference.iic"))
        'a1b2c3d4e5f6...'
    """
    md5 = hashlib.md5()

    with open(file_path, 'rb') as f:
        while chunk := f.read(chunk_size):
            md5.update(chunk)

    return md5.hexdigest()


def generate_training_metadata(
    model_name: str,
    u12_reference_path: Optional[Path],
    u2_reference_path: Optional[Path],
    u12_introns: Sequence[Intron],
    u2_introns: Sequence[Intron],
    optimized_C: float,
    calibration_method: str,
    cv_score: float,
    n_models: int,
    threshold: float,
    eval_result: Optional[Any] = None,
    max_iter: int = 100000,
    kernel: str = 'linear',
    seed: int = 42
) -> Dict[str, Any]:
    """
    Generate comprehensive training metadata.

    Args:
        model_name: Name of the model
        u12_reference_path: Path to U12 reference file (None if default)
        u2_reference_path: Path to U2 reference file (None if default)
        u12_introns: U12 reference introns used in training
        u2_introns: U2 reference introns used in training
        optimized_C: Optimized C parameter
        calibration_method: Calibration method used
        cv_score: Cross-validation score
        n_models: Number of ensemble models
        threshold: Classification threshold
        eval_result: Optional evaluation results (NestedCVResult or SplitEvalResult)
        max_iter: Maximum SVM iterations
        kernel: SVM kernel type
        seed: Random seed

    Returns:
        Dictionary with complete training metadata
    """
    # Determine reference file paths (use defaults if not specified)
    if u12_reference_path is None:
        u12_reference_path = Path(__file__).parent.parent / "data" / "u12_reference.introns.iic.gz"
    if u2_reference_path is None:
        u2_reference_path = Path(__file__).parent.parent / "data" / "u2_reference.introns.iic.gz"

    # Compute MD5 hashes for reproducibility
    u12_md5 = compute_file_md5(u12_reference_path) if u12_reference_path.exists() else None
    u2_md5 = compute_file_md5(u2_reference_path) if u2_reference_path.exists() else None

    metadata = {
        "model_name": model_name,
        "training_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "reference_data": {
            "u12": {
                "file": str(u12_reference_path),
                "md5": u12_md5,
                "count": len(u12_introns)
            },
            "u2": {
                "file": str(u2_reference_path),
                "md5": u2_md5,
                "count": len(u2_introns)
            },
            "total": len(u12_introns) + len(u2_introns)
        },
        "hyperparameters": {
            "C": optimized_C,
            "cv_score": cv_score,
            "kernel": kernel,
            "calibration_method": calibration_method,
            "max_iter": max_iter,
            "n_models": n_models,
            "threshold": threshold,
            "random_seed": seed
        }
    }

    # Add evaluation metrics if available
    if eval_result is not None:
        if hasattr(eval_result, 'mean_f1'):
            # Nested CV result
            metadata["evaluation"] = {
                "method": "nested_cv",
                "n_folds": eval_result.n_folds,
                "mean_f1": round(eval_result.mean_f1, 4),
                "std_f1": round(eval_result.std_f1, 4),
                "mean_pr_auc": round(eval_result.mean_pr_auc, 4),
                "std_pr_auc": round(eval_result.std_pr_auc, 4)
            }
        elif hasattr(eval_result, 'test_f1'):
            # Split evaluation result
            metadata["evaluation"] = {
                "method": "train_val_test_split",
                "test_f1": round(eval_result.test_f1, 4),
                "test_pr_auc": round(eval_result.test_pr_auc, 4),
                "split": {
                    "train": {
                        "u12": eval_result.n_u12_train,
                        "u2": eval_result.n_u2_train,
                        "fraction": eval_result.train_fraction
                    },
                    "validation": {
                        "u12": eval_result.n_u12_val,
                        "u2": eval_result.n_u2_val,
                        "fraction": eval_result.val_fraction
                    },
                    "test": {
                        "u12": eval_result.n_u12_test,
                        "u2": eval_result.n_u2_test,
                        "fraction": eval_result.test_fraction
                    }
                }
            }
        else:
            # No evaluation performed
            metadata["evaluation"] = None
    else:
        metadata["evaluation"] = None

    return metadata


def generate_pretrained_metadata(
    model_path: Path,
    threshold: float
) -> Dict[str, Any]:
    """
    Generate minimal metadata for runs using pretrained models.

    Args:
        model_path: Path to pretrained model file
        threshold: Classification threshold used

    Returns:
        Dictionary with pretrained model metadata
    """
    # Try to load accompanying metadata if it exists
    metadata_path = model_path.with_suffix('.metadata.json')

    if metadata_path.exists():
        # Load existing metadata and add usage info
        with open(metadata_path, 'r') as f:
            existing = json.load(f)

        return {
            "model_name": model_path.stem,
            "model_source": "pretrained",
            "model_path": str(model_path),
            "model_md5": compute_file_md5(model_path),
            "threshold": threshold,
            "pretrained_metadata": existing,
            "run_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
    else:
        # No existing metadata - just record what we used
        return {
            "model_name": model_path.stem,
            "model_source": "pretrained",
            "model_path": str(model_path),
            "model_md5": compute_file_md5(model_path),
            "threshold": threshold,
            "run_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "note": "No training metadata available for this pretrained model"
        }


def write_metadata(
    metadata: Dict[str, Any],
    output_path: Path
) -> None:
    """
    Write metadata to JSON file.

    Args:
        metadata: Metadata dictionary
        output_path: Path to output JSON file
    """
    with open(output_path, 'w') as f:
        json.dump(metadata, f, indent=2)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
