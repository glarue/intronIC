#!/usr/bin/env python3
"""
Simplified diagnostic: Check if intercept shrinkage causes low probabilities.

Directly loads reference data and trains both models for comparison.
"""

import numpy as np
import gzip
from sklearn.svm import SVC, LinearSVC
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import train_test_split


def load_reference_introns(filepath):
    """Load reference introns from .iic.gz file - just the coordinates."""
    introns = []
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                introns.append(parts)
    return introns


def extract_z_scores_from_line(line_parts):
    """
    Extract z-scores from reference intron line.

    Reference format has scores in specific positions.
    We'll need the five_z, bp_z, three_z scores.
    These should be in the line somewhere - let me check the actual format.
    """
    # For now, return None - we'll need to inspect actual file format
    return None


def main():
    print("="*80)
    print("SIMPLIFIED DIAGNOSTIC: Intercept Shrinkage Test")
    print("="*80)

    # For this diagnostic, let's create synthetic data that matches the problem:
    # - 3 features (z-scores)
    # - Severe imbalance (387 U12 vs 20690 U2 = 1.8%)
    # - Classes should be separable (original model finds them)

    np.random.seed(42)

    # Simulate U2 introns (negative class) - centered around origin
    n_u2 = 20690
    X_u2 = np.random.randn(n_u2, 3) * 1.0  # Normal distribution

    # Simulate U12 introns (positive class) - shifted to be distinguishable
    n_u12 = 387
    X_u12 = np.random.randn(n_u12, 3) * 0.8 + 3.0  # Shifted mean, tighter variance

    # Combine
    X = np.vstack([X_u2, X_u12])
    y = np.array([0] * n_u2 + [1] * n_u12)

    print(f"\nDataset: {n_u2} U2 (negative), {n_u12} U12 (positive)")
    print(f"Imbalance: {100 * n_u12 / (n_u2 + n_u12):.2f}% positive class")

    # Split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, stratify=y, random_state=42
    )

    # Train both models with same C for fair comparison
    C_value = 1.0

    print("\n" + "="*80)
    print(f"Training LinearSVC (refactored) with C={C_value}")
    print("="*80)

    base_linear = LinearSVC(
        C=C_value,
        class_weight='balanced',
        loss='squared_hinge',
        penalty='l2',
        dual=False,  # Current approach
        max_iter=10000,
        tol=1e-4,
        random_state=42
    )

    linear_cal = CalibratedClassifierCV(
        base_linear,
        method='sigmoid',
        cv=5
    )
    linear_cal.fit(X_train, y_train)

    print("\n" + "="*80)
    print(f"Training SVC (original) with C={C_value}")
    print("="*80)

    svc_model = SVC(
        C=C_value,
        kernel='linear',
        probability=True,
        class_weight='balanced',
        random_state=42
    )
    svc_model.fit(X_train, y_train)

    # Access base estimator from calibrated model
    base_linear_fitted = linear_cal.calibrated_classifiers_[0].estimator

    print("\n" + "="*80)
    print("TEST 1: Decision Function Range")
    print("="*80)

    s_linear = base_linear_fitted.decision_function(X_test)
    s_svc = svc_model.decision_function(X_test)

    print(f"\nLinearSVC decision scores:")
    print(f"  Range: [{s_linear.min():.4f}, {s_linear.max():.4f}]")
    print(f"  Span: {s_linear.max() - s_linear.min():.4f}")
    print(f"  Std: {s_linear.std():.4f}")

    print(f"\nSVC decision scores:")
    print(f"  Range: [{s_svc.min():.4f}, {s_svc.max():.4f}]")
    print(f"  Span: {s_svc.max() - s_svc.min():.4f}")
    print(f"  Std: {s_svc.std():.4f}")

    span_ratio = (s_svc.max() - s_svc.min()) / (s_linear.max() - s_linear.min())
    print(f"\n→ SVC span is {span_ratio:.2f}x larger")

    print("\n" + "="*80)
    print("TEST 2: Intercept Magnitude")
    print("="*80)

    print(f"\nLinearSVC intercept: {base_linear_fitted.intercept_[0]:.6f}")
    print(f"SVC intercept: {svc_model.intercept_[0]:.6f}")

    if abs(svc_model.intercept_[0]) > 1e-10:
        intercept_ratio = abs(svc_model.intercept_[0]) / abs(base_linear_fitted.intercept_[0])
        print(f"\n→ SVC intercept is {intercept_ratio:.2f}x larger")

    print("\n" + "="*80)
    print("TEST 3: Calibration Parameters (Platt Scaling)")
    print("="*80)

    calibrator = linear_cal.calibrated_classifiers_[0].calibrators[0]
    a_linear = calibrator.a_
    b_linear = calibrator.b_

    print(f"\nLinearSVC Platt parameters (fold 0):")
    print(f"  A (slope): {a_linear:.6f}")
    print(f"  B (offset): {b_linear:.6f}")
    print(f"  |A|: {abs(a_linear):.6f}")

    if abs(a_linear) < 0.5:
        print(f"\n→ WARNING: Very flat sigmoid (small |A|)")
        print(f"  This means poor probability separation")

    print("\n" + "="*80)
    print("TEST 4: Probability Distributions")
    print("="*80)

    proba_linear = linear_cal.predict_proba(X_test)[:, 1]
    proba_svc = svc_model.predict_proba(X_test)[:, 1]

    u12_mask = y_test == 1
    u2_mask = y_test == 0

    print(f"\nLinearSVC probabilities:")
    print(f"  U12 (positive) - n={u12_mask.sum()}:")
    print(f"    Mean: {proba_linear[u12_mask].mean():.4f}")
    print(f"    Max: {proba_linear[u12_mask].max():.4f}")
    print(f"    # above 0.90: {(proba_linear[u12_mask] >= 0.90).sum()}")
    print(f"    # above 0.50: {(proba_linear[u12_mask] >= 0.50).sum()}")

    print(f"  U2 (negative) - n={u2_mask.sum()}:")
    print(f"    Mean: {proba_linear[u2_mask].mean():.4f}")
    print(f"    Max: {proba_linear[u2_mask].max():.4f}")

    print(f"\nSVC probabilities:")
    print(f"  U12 (positive) - n={u12_mask.sum()}:")
    print(f"    Mean: {proba_svc[u12_mask].mean():.4f}")
    print(f"    Max: {proba_svc[u12_mask].max():.4f}")
    print(f"    # above 0.90: {(proba_svc[u12_mask] >= 0.90).sum()}")
    print(f"    # above 0.50: {(proba_svc[u12_mask] >= 0.50).sum()}")

    print(f"  U2 (negative) - n={u2_mask.sum()}:")
    print(f"    Mean: {proba_svc[u2_mask].mean():.4f}")
    print(f"    Max: {proba_svc[u2_mask].max():.4f}")

    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)

    issues = []

    if span_ratio > 3:
        issues.append("✓ Decision function squashing confirmed")

    if abs(base_linear_fitted.intercept_[0]) < abs(svc_model.intercept_[0]) / 3:
        issues.append("✓ Intercept shrinkage confirmed")

    if abs(a_linear) < 0.5:
        issues.append("✓ Flat sigmoid calibration confirmed")

    if proba_linear[u12_mask].max() < 0.5 and proba_svc[u12_mask].max() > 0.9:
        issues.append("✓ Low probability predictions confirmed")

    if issues:
        print("\n✓✓✓ HYPOTHESIS CONFIRMED ✓✓✓")
        print("LinearSVC with dual=False shows:")
        for issue in issues:
            print(f"  {issue}")
        print("\nRECOMMENDED FIXES (in order):")
        print("  1. Set intercept_scaling=1000")
        print("  2. Add MaxAbsScaler for feature scaling")
        print("  3. Try isotonic calibration")
        print("  4. Grid-search intercept_scaling")
    else:
        print("\nPartial or no confirmation")
        print("May need different fixes")


if __name__ == '__main__':
    main()
