"""
Utility to inspect learned SVM weights and verify BothEndsStrong penalties.

This module provides tools to examine trained models and verify that the
imbalance penalty features (absdiff_5_bp, absdiff_5_3) have negative weights
as intended.

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
    Inspect learned SVM coefficients and verify penalties work as intended.

    According to expert guidance:
    - absdiff_5_bp and absdiff_5_3 should have NEGATIVE coefficients
    - If not, the imbalance penalty is not working correctly
    - Solution: increase gamma or add min_* features

    Args:
        ensemble: Trained SVM ensemble to inspect
        verbose: Print detailed coefficient information

    Returns:
        List of warning messages (empty if all checks pass)
    """
    feature_names = [
        's5', 'sBP', 's3',
        'sum_5_bp', 'absdiff_5_bp',
        'sum_5_3', 'absdiff_5_3'
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
            print(f"    gamma_5_bp: {model.parameters.gamma_5_bp}")
            print(f"    gamma_5_3:  {model.parameters.gamma_5_3}")
            print(f"    C:          {model.parameters.C}")

        # Extract SVC from calibrated classifier
        # Pipeline: RobustScaler -> BothEndsStrong -> LinearSVC -> CalibratedClassifierCV
        svc = model.model.calibrated_classifiers_[0].estimator
        coefs = svc.coef_[0]
        intercept = svc.intercept_[0]

        if verbose:
            print(f"\n  Learned coefficients:")
            for name, coef in zip(feature_names, coefs):
                # Highlight the penalty features
                if 'absdiff' in name:
                    status = "✓ NEGATIVE (penalty)" if coef < 0 else "⚠️  POSITIVE (NOT penalizing!)"
                    print(f"    {name:15s}: {coef:+.6f}  {status}")
                else:
                    print(f"    {name:15s}: {coef:+.6f}")

            print(f"\n  Intercept: {intercept:+.6f}")

        # Sanity checks (per expert guidance)
        absdiff_5_bp_coef = coefs[4]  # Index 4 is absdiff_5_bp
        absdiff_5_3_coef = coefs[6]   # Index 6 is absdiff_5_3

        if absdiff_5_bp_coef > 0:
            warning = (
                f"Model {i+1}: absdiff_5_bp has POSITIVE coefficient ({absdiff_5_bp_coef:+.6f}). "
                f"This means imbalanced 5'-BP introns are MORE likely to be classified as U12! "
                f"Increase gamma_5_bp or add min_* features."
            )
            warnings.append(warning)
            if verbose:
                print(f"\n  ⚠️  WARNING: {warning}")

        if absdiff_5_3_coef > 0:
            warning = (
                f"Model {i+1}: absdiff_5_3 has POSITIVE coefficient ({absdiff_5_3_coef:+.6f}). "
                f"This means imbalanced 5'-3' introns are MORE likely to be classified as U12! "
                f"Increase gamma_5_3 or add min_* features."
            )
            warnings.append(warning)
            if verbose:
                print(f"\n  ⚠️  WARNING: {warning}")

        # Additional insight: check if penalties are "strong enough"
        # A very small negative coefficient means the penalty is weak
        if verbose and absdiff_5_bp_coef < 0 and abs(absdiff_5_bp_coef) < 0.01:
            print(f"  ℹ️  Note: absdiff_5_bp coefficient is very small ({absdiff_5_bp_coef:.6f}). "
                  f"Penalty may be weak. Consider increasing gamma_5_bp.")

        if verbose and absdiff_5_3_coef < 0 and abs(absdiff_5_3_coef) < 0.01:
            print(f"  ℹ️  Note: absdiff_5_3 coefficient is very small ({absdiff_5_3_coef:.6f}). "
                  f"Penalty may be weak. Consider increasing gamma_5_3.")

    if verbose:
        print("\n" + "=" * 70)
        if not warnings:
            print("✓ All sanity checks passed! Imbalance penalties are negative.")
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
        'sum_5_bp', 'absdiff_5_bp',
        'sum_5_3', 'absdiff_5_3'
    ]

    # Collect coefficients from all models
    all_coefs = []
    for model in ensemble.models:
        svc = model.model.calibrated_classifiers_[0].estimator
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
        # Highlight penalty features
        marker = "  *" if 'absdiff' in name else "   "
        print(f"{name:<15}{marker} {stats['mean']:+12.6f} {stats['std']:12.6f} "
              f"{stats['min']:+12.6f} {stats['max']:+12.6f}")

    print("=" * 70)
    print("* = Imbalance penalty features (should be negative)\n")
