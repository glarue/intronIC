#!/usr/bin/env python3
"""
Diagnostic script to check if intercept shrinkage is causing low probabilities.

Tests:
1. Decision function range: LinearSVC vs SVC
2. Intercept magnitude: LinearSVC vs SVC
3. Calibration slope: Platt A, B parameters

Hypothesis: LinearSVC with dual=False has shrunken intercept → squashed scores → low probabilities
"""

import sys
import numpy as np
from pathlib import Path

# Add refactored code to path
sys.path.insert(0, str(Path(__file__).parent / 'intronIC_refactored'))

from sklearn.svm import SVC, LinearSVC
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report

# Import data loading utilities
from file_io.parsers import SequenceParser
from scoring.pwm import PWMLoader
from scoring.scorer import IntronScorer
from scoring.normalizer import ScoreNormalizer


def load_reference_data():
    """Load and prepare reference U12/U2 data."""
    print("Loading reference data...")

    u12_path = Path('intronIC_refactored/intronIC/data/u12_reference.introns.iic.gz')
    u2_path = Path('intronIC_refactored/intronIC/data/u2_reference.introns.iic.gz')

    parser = SequenceParser()
    u12_introns = parser.parse(u12_path)
    u2_introns = parser.parse(u2_path)

    print(f"Loaded {len(u12_introns)} U12, {len(u2_introns)} U2 reference introns")

    # Score them
    pwm_path = Path('intronIC_refactored/intronIC/data/scoring_matrices.fasta.iic')
    pwm_loader = PWMLoader(pwm_path)
    matrices = pwm_loader.load()

    scorer = IntronScorer(matrices)
    u12_introns = scorer.score_introns(u12_introns)
    u2_introns = scorer.score_introns(u2_introns)

    # Normalize
    normalizer = ScoreNormalizer()
    normalizer.fit(u12_introns, u2_introns)
    u12_introns = normalizer.transform(u12_introns)
    u2_introns = normalizer.transform(u2_introns)

    # Extract features
    X_u12 = np.array([[i.scores.five_z_score, i.scores.bp_z_score, i.scores.three_z_score]
                      for i in u12_introns])
    X_u2 = np.array([[i.scores.five_z_score, i.scores.bp_z_score, i.scores.three_z_score]
                     for i in u2_introns])

    X = np.vstack([X_u2, X_u12])
    y = np.array([0] * len(u2_introns) + [1] * len(u12_introns))

    return X, y


def train_models(X, y):
    """Train both LinearSVC and SVC models."""
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, stratify=y, random_state=42
    )

    print("\n" + "="*80)
    print("Training LinearSVC (current refactored approach)...")
    print("="*80)

    # Current refactored approach
    base_linear = LinearSVC(
        C=1.0,  # Use same C for fair comparison
        class_weight='balanced',
        loss='squared_hinge',
        penalty='l2',
        dual=False,
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
    print("Training SVC (original working approach)...")
    print("="*80)

    # Original working approach
    svc_model = SVC(
        C=1.0,  # Use same C for fair comparison
        kernel='linear',
        probability=True,
        class_weight='balanced',
        random_state=42
    )
    svc_model.fit(X_train, y_train)

    return linear_cal, svc_model, X_train, X_test, y_train, y_test


def run_diagnostics():
    """Run all diagnostic tests."""
    print("="*80)
    print("DIAGNOSTIC: Intercept Shrinkage Investigation")
    print("="*80)

    X, y = load_reference_data()
    linear_cal, svc_model, X_train, X_test, y_train, y_test = train_models(X, y)

    # Access the base estimators
    base_linear = linear_cal.calibrated_classifiers_[0].estimator

    print("\n" + "="*80)
    print("TEST 1: Decision Function Range")
    print("="*80)
    print("Hypothesis: LinearSVC scores are squashed (smaller range)")

    s_linear = base_linear.decision_function(X_test)
    s_svc = svc_model.decision_function(X_test)

    print(f"\nLinearSVC decision function:")
    print(f"  Range: [{s_linear.min():.4f}, {s_linear.max():.4f}]")
    print(f"  Span: {s_linear.max() - s_linear.min():.4f}")
    print(f"  Percentiles [1%, 99%]: [{np.percentile(s_linear, 1):.4f}, {np.percentile(s_linear, 99):.4f}]")
    print(f"  Std: {s_linear.std():.4f}")

    print(f"\nSVC decision function:")
    print(f"  Range: [{s_svc.min():.4f}, {s_svc.max():.4f}]")
    print(f"  Span: {s_svc.max() - s_svc.min():.4f}")
    print(f"  Percentiles [1%, 99%]: [{np.percentile(s_svc, 1):.4f}, {np.percentile(s_svc, 99):.4f}]")
    print(f"  Std: {s_svc.std():.4f}")

    span_ratio = (s_svc.max() - s_svc.min()) / (s_linear.max() - s_linear.min())
    print(f"\n  → SVC span is {span_ratio:.2f}x larger than LinearSVC")
    if span_ratio > 2:
        print("  ✓ CONFIRMED: LinearSVC scores are squashed")

    print("\n" + "="*80)
    print("TEST 2: Intercept Magnitude")
    print("="*80)
    print("Hypothesis: LinearSVC intercept is shrunk toward 0")

    print(f"\nLinearSVC intercept: {base_linear.intercept_[0]:.6f}")
    print(f"SVC intercept: {svc_model.intercept_[0]:.6f}")

    intercept_ratio = abs(svc_model.intercept_[0]) / abs(base_linear.intercept_[0])
    print(f"\n  → SVC intercept is {intercept_ratio:.2f}x larger")
    if intercept_ratio > 2:
        print("  ✓ CONFIRMED: LinearSVC intercept is overly regularized")

    print("\n" + "="*80)
    print("TEST 3: Calibration Slope (Platt Parameters)")
    print("="*80)
    print("Hypothesis: Flat sigmoid (small |A|) → probabilities stay near base rate")

    print("\nLinearSVC calibration (first fold):")
    calibrator = linear_cal.calibrated_classifiers_[0].calibrators_[0]
    a_linear = calibrator.a_
    b_linear = calibrator.b_
    print(f"  Platt A: {a_linear:.6f}")
    print(f"  Platt B: {b_linear:.6f}")
    print(f"  |A|: {abs(a_linear):.6f}")

    # SVC already has internal calibration, but we can check its probability behavior
    print("\nSVC internal calibration:")
    print("  (Parameters not directly accessible, but behavior is different)")

    if abs(a_linear) < 0.1:
        print("\n  ✓ CONFIRMED: Very flat sigmoid (small |A|) → poor probability separation")

    print("\n" + "="*80)
    print("TEST 4: Actual Probability Distributions")
    print("="*80)

    # Get probabilities on test set
    proba_linear = linear_cal.predict_proba(X_test)[:, 1]
    proba_svc = svc_model.predict_proba(X_test)[:, 1]

    # Separate by true class
    u12_mask = y_test == 1
    u2_mask = y_test == 0

    print("\nLinearSVC probabilities:")
    print(f"  U12 introns (positive class):")
    print(f"    Mean: {proba_linear[u12_mask].mean():.4f}")
    print(f"    Max: {proba_linear[u12_mask].max():.4f}")
    print(f"    Min: {proba_linear[u12_mask].min():.4f}")
    print(f"    # above 0.90: {(proba_linear[u12_mask] >= 0.90).sum()}/{u12_mask.sum()}")

    print(f"  U2 introns (negative class):")
    print(f"    Mean: {proba_linear[u2_mask].mean():.4f}")
    print(f"    Max: {proba_linear[u2_mask].max():.4f}")

    print("\nSVC probabilities:")
    print(f"  U12 introns (positive class):")
    print(f"    Mean: {proba_svc[u12_mask].mean():.4f}")
    print(f"    Max: {proba_svc[u12_mask].max():.4f}")
    print(f"    Min: {proba_svc[u12_mask].min():.4f}")
    print(f"    # above 0.90: {(proba_svc[u12_mask] >= 0.90).sum()}/{u12_mask.sum()}")

    print(f"  U2 introns (negative class):")
    print(f"    Mean: {proba_svc[u2_mask].mean():.4f}")
    print(f"    Max: {proba_svc[u2_mask].max():.4f}")

    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)

    max_u12_linear = proba_linear[u12_mask].max()
    max_u12_svc = proba_svc[u12_mask].max()

    if max_u12_linear < 0.5 and max_u12_svc > 0.9:
        print("\n✓✓✓ HYPOTHESIS CONFIRMED ✓✓✓")
        print("LinearSVC + external calibration produces systematically low probabilities")
        print("due to intercept shrinkage and score squashing.")
        print("\nRECOMMENDED FIXES:")
        print("1. Set intercept_scaling=1000 (or switch dual=True)")
        print("2. Add feature scaling (MaxAbsScaler)")
        print("3. Try isotonic calibration (more flexible)")
    else:
        print("\n❌ HYPOTHESIS NOT CONFIRMED")
        print("Different issue causing low probabilities")


if __name__ == '__main__':
    run_diagnostics()
