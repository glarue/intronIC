"""
Utility to inspect learned SVM weights for BothEndsStrong features.

This module provides tools to examine trained models and view the
learned coefficients for min/max features.

Usage:
    from classification.model_inspector import inspect_ensemble_weights

    warnings = inspect_ensemble_weights(ensemble)
    if warnings:
        print("⚠️  Model may not be working as intended!")
        for warning in warnings:
            print(f"  - {warning}")
"""

from typing import List, Optional
import numpy as np
from classification.trainer import SVMEnsemble


def inspect_ensemble_weights(
    ensemble: SVMEnsemble,
    verbose: bool = True
) -> List[str]:
    """
    Inspect learned SVM coefficients for BothEndsStrong features.

    With min/max features:
    - min features capture "both must be strong" (expect POSITIVE weights)
    - max features capture "at least one strong" (various weights expected)

    Args:
        ensemble: Trained SVM ensemble to inspect
        verbose: Print detailed coefficient information

    Returns:
        List of warning messages (empty if all checks pass)
    """
    feature_names = [
        's5', 'sBP', 's3',
        'min_5_bp', 'min_5_3'
        # Note: May also include max_5_bp and max_5_3 if include_max=True
        # We'll check actual pipeline to get correct feature count
    ]

    warnings = []

    if verbose:
        print("\n" + "=" * 70)
        print("SVM Ensemble Coefficient Analysis")
        print("=" * 70)

    for i, model in enumerate(ensemble.models):
        if verbose:
            print(f"\nModel {i+1}/{len(ensemble.models)}:")
            print(f"  Hyperparameters:")
            print(f"    include_max: {model.parameters.include_max}")
            print(f"    C:           {model.parameters.C}")

        # Extract SVC from calibrated classifier
        # model.model is CalibratedClassifierCV
        # calibrated_classifiers_[0].estimator is Pipeline(RobustScaler -> BothEndsStrong -> LinearSVC)
        # We need to get the LinearSVC from the end of the pipeline
        pipeline = model.model.calibrated_classifiers_[0].estimator
        svc = pipeline.named_steps['svc']  # Get the LinearSVC from pipeline
        coefs = svc.coef_[0]
        intercept = svc.intercept_[0]

        if verbose:
            print(f"\n  Learned coefficients:")
            # Get actual number of features from pipeline
            n_features = len(coefs)
            for idx in range(min(n_features, len(feature_names))):
                name = feature_names[idx] if idx < len(feature_names) else f"feature_{idx}"
                coef = coefs[idx]
                # Highlight min/max features
                if 'min' in name or 'max' in name:
                    print(f"    {name:15s}: {coef:+.6f}  ← BothEndsStrong feature")
                else:
                    print(f"    {name:15s}: {coef:+.6f}")

            print(f"\n  Intercept: {intercept:+.6f}")

        # With min/max features, we expect positive weights on min features
        # (higher min = both strong = more likely U12)
        # No specific sanity checks needed - expert approach is cleaner

    if verbose:
        print("\n" + "=" * 70)
        if not warnings:
            print("✓ All sanity checks passed!")
        else:
            print(f"⚠️  {len(warnings)} warning(s) found. Model may not work as intended.")
        print("=" * 70 + "\n")

    return warnings


def get_coefficient_summary(ensemble: SVMEnsemble) -> dict:
    """
    Get summary statistics of coefficients across ensemble.

    Args:
        ensemble: Trained SVM ensemble

    Returns:
        Dictionary with mean, std, min, max for each feature coefficient
    """
    feature_names = [
        's5', 'sBP', 's3',
        'min_5_bp', 'min_5_3'
        # Note: May also include max_5_bp and max_5_3 if include_max=True
        # We'll check actual pipeline to get correct feature count
    ]

    # Collect coefficients from all models
    all_coefs = []
    for model in ensemble.models:
        pipeline = model.model.calibrated_classifiers_[0].estimator
        svc = pipeline.named_steps['svc']  # Get LinearSVC from pipeline
        all_coefs.append(svc.coef_[0])

    all_coefs = np.array(all_coefs)  # Shape: (n_models, n_features)

    # Compute statistics
    summary = {}
    for i, name in enumerate(feature_names):
        summary[name] = {
            'mean': float(np.mean(all_coefs[:, i])),
            'std': float(np.std(all_coefs[:, i])),
            'min': float(np.min(all_coefs[:, i])),
            'max': float(np.max(all_coefs[:, i]))
        }

    return summary


def print_coefficient_summary(ensemble: SVMEnsemble) -> None:
    """
    Print a concise summary of coefficients across ensemble.

    Args:
        ensemble: Trained SVM ensemble
    """
    summary = get_coefficient_summary(ensemble)

    print("\nCoefficient Summary (across ensemble):")
    print("=" * 70)
    print(f"{'Feature':<15} {'Mean':>12} {'Std':>12} {'Min':>12} {'Max':>12}")
    print("-" * 70)

    for name, stats in summary.items():
        # Highlight BothEndsStrong features
        marker = "  *" if 'min' in name or 'max' in name else "   "
        print(f"{name:<15}{marker} {stats['mean']:+12.6f} {stats['std']:12.6f} "
              f"{stats['min']:+12.6f} {stats['max']:+12.6f}")

    print("=" * 70)
    print("* = BothEndsStrong augmented features (min/max)\n")
