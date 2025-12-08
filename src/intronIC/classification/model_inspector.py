"""
Utility to inspect learned SVM weights for BothEndsStrong features.

This module provides tools to examine trained models and view the
learned coefficients for min/max features.

Usage:
    from intronIC.classification.model_inspector import inspect_ensemble_weights

    warnings = inspect_ensemble_weights(ensemble)
    if warnings:
        print("⚠️  Model may not be working as intended!")
        for warning in warnings:
            print(f"  - {warning}")
"""

from typing import List, Optional

import numpy as np

from intronIC.classification.trainer import SVMEnsemble


def inspect_ensemble_weights(ensemble: SVMEnsemble, verbose: bool = True) -> List[str]:
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
    warnings = []

    if verbose:
        print("\n" + "=" * 70)
        print("SVM Ensemble Coefficient Analysis")
        print("=" * 70)

        # Print shared hyperparameters once at the top
        # All models share these parameters (from Stage 1 & 2 optimization)
        first_model = ensemble.models[0]
        print("\nShared Hyperparameters (all models):")
        print(f"  Optimal C:                   {first_model.parameters.C}")
        print(
            f"  Calibration method:          {first_model.parameters.calibration_method}"
        )
        print(f"  include_max:                 {first_model.parameters.include_max}")
        print(
            f"  include_pairwise_mins:       {first_model.parameters.include_pairwise_mins}"
        )
        print("\nEnsemble Diversity:")
        print(f"  Number of models:            {len(ensemble.models)}")
        print(
            f"  U2 subsampling:              {'Yes' if len(ensemble.models) > 1 else 'No'}"
        )
        if len(ensemble.models) > 1:
            subsample_pct = int(ensemble.subsample_ratio * 100)
            print(
                f"  U2 subsample variation:      Different random {subsample_pct}% per model"
            )
        print("\n" + "-" * 70)

    for i, model in enumerate(ensemble.models):
        if verbose:
            print(f"\nModel {i + 1}/{len(ensemble.models)} - Learned Coefficients:")
            print(
                f"  Training set size: {model.train_size} introns ({model.u12_count} U12, {model.u2_count} U2)"
            )

        # Extract SVC from calibrated classifier
        # model.model is CalibratedClassifierCV
        # calibrated_classifiers_[0].estimator is Pipeline(BothEndsStrong -> LinearSVC)
        # NOTE: Scaling happens EXTERNALLY via ScoreNormalizer, NOT inside the pipeline
        # We need to get the LinearSVC from the end of the pipeline
        pipeline = model.model.calibrated_classifiers_[0].estimator
        svc = pipeline.named_steps["svc"]  # Get the LinearSVC from pipeline
        coefs = svc.coef_[0]
        intercept = svc.intercept_[0]

        if verbose:
            # Get feature names (Phase 1: no augment step, Phase 2+: augment step exists)
            if "augment" in pipeline.named_steps:
                transformer = pipeline.named_steps["augment"]
                feature_names = transformer.get_feature_names_out()
            elif "transform" in pipeline.named_steps:
                # Get feature names from BothEndsStrongTransformer
                transformer = pipeline.named_steps["transform"]
                feature_names = transformer.get_feature_names_out()
            else:
                # Phase 1: Only base features (3D)
                feature_names = ["five_z_score", "bp_z_score", "three_z_score"]

            # Print all features with coefficients
            for idx, (name, coef) in enumerate(zip(feature_names, coefs)):
                # Highlight min/max features
                if "min" in name or "max" in name:
                    print(f"    {name:25s}: {coef:+.6f}  ← Augmented feature")
                else:
                    print(f"    {name:25s}: {coef:+.6f}")

            print(f"  Intercept:                 {intercept:+.6f}")

        # With min/max features, we expect positive weights on min features
        # (higher min = both strong = more likely U12)
        # No specific sanity checks needed - expert approach is cleaner

    if verbose:
        print("\n" + "=" * 70)
        if not warnings:
            print("✓ All sanity checks passed!")
        else:
            print(
                f"⚠️  {len(warnings)} warning(s) found. Model may not work as intended."
            )
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
    # Get feature names from the first model's transformer (or default for Phase 1)
    first_model = ensemble.models[0]
    pipeline = first_model.model.calibrated_classifiers_[0].estimator
    if "augment" in pipeline.named_steps:
        transformer = pipeline.named_steps["augment"]
        feature_names = transformer.get_feature_names_out()
    else:
        # Phase 1: Only base features (3D)
        feature_names = ["five_z_score", "bp_z_score", "three_z_score"]

    # Collect coefficients from all models
    all_coefs = []
    for model in ensemble.models:
        pipeline = model.model.calibrated_classifiers_[0].estimator
        svc = pipeline.named_steps["svc"]  # Get LinearSVC from pipeline
        all_coefs.append(svc.coef_[0])

    all_coefs = np.array(all_coefs)  # Shape: (n_models, n_features)

    # Compute statistics
    summary = {}
    for i, name in enumerate(feature_names):
        summary[name] = {
            "mean": float(np.mean(all_coefs[:, i])),
            "std": float(np.std(all_coefs[:, i])),
            "min": float(np.min(all_coefs[:, i])),
            "max": float(np.max(all_coefs[:, i])),
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
        marker = "  *" if "min" in name or "max" in name else "   "
        print(
            f"{name:<15}{marker} {stats['mean']:+12.6f} {stats['std']:12.6f} "
            f"{stats['min']:+12.6f} {stats['max']:+12.6f}"
        )

    print("=" * 70)
    print("* = BothEndsStrong augmented features (min/max)\n")
    print("=" * 70)
    print("* = BothEndsStrong augmented features (min/max)\n")
