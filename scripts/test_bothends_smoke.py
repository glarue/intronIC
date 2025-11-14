#!/usr/bin/env python3
"""
Smoke test for BothEndsStrongTransformer implementation.

Quick validation that the feature augmentation is working correctly
without requiring full dataset training.

Usage:
    python scripts/test_bothends_smoke.py
"""

import sys
import numpy as np
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from classification.transformers import BothEndsStrongTransformer


def test_basic_transformation():
    """Test that 3D → 7D transformation works."""
    print("=" * 80)
    print("TEST 1: Basic Transformation")
    print("=" * 80)

    transformer = BothEndsStrongTransformer(gamma_5_bp=4, gamma_5_3=2)

    # Single sample
    X = np.array([[2.0, 1.5, -0.5]])
    X_aug = transformer.fit_transform(X)

    print(f"Input shape: {X.shape} (expected: (1, 3))")
    print(f"Output shape: {X_aug.shape} (expected: (1, 7))")

    assert X_aug.shape == (1, 7), f"Expected (1, 7), got {X_aug.shape}"
    print("✓ Shape transformation correct: 3D → 7D")

    # Multiple samples
    X_multi = np.random.randn(100, 3)
    X_multi_aug = transformer.transform(X_multi)
    assert X_multi_aug.shape == (100, 7), f"Expected (100, 7), got {X_multi_aug.shape}"
    print("✓ Batch transformation works")
    print()

    return True


def test_feature_values():
    """Test that feature calculations are correct."""
    print("=" * 80)
    print("TEST 2: Feature Value Calculations")
    print("=" * 80)

    transformer = BothEndsStrongTransformer(gamma_5_bp=4, gamma_5_3=2)

    # Test case: s5=2.0, sBP=1.5, s3=-0.5
    X = np.array([[2.0, 1.5, -0.5]])
    X_aug = transformer.transform(X)

    print(f"Input: s5={X[0,0]}, sBP={X[0,1]}, s3={X[0,2]}")
    print(f"\nExpected output features:")
    print(f"  [0] s5 (pass-through):        {X[0,0]:.3f}")
    print(f"  [1] sBP (pass-through):       {X[0,1]:.3f}")
    print(f"  [2] s3 (pass-through):        {X[0,2]:.3f}")
    print(f"  [3] sum_5_bp (s5 + sBP):      {2.0 + 1.5:.3f}")
    print(f"  [4] absdiff_5_bp (|s5-sBP|×4): {abs(2.0 - 1.5) * 4:.3f}")
    print(f"  [5] sum_5_3 (s5 + s3):        {2.0 + (-0.5):.3f}")
    print(f"  [6] absdiff_5_3 (|s5-s3|×2):  {abs(2.0 - (-0.5)) * 2:.3f}")

    print(f"\nActual output features:")
    for i, val in enumerate(X_aug[0]):
        print(f"  [{i}] {val:.3f}")

    # Verify calculations
    expected = np.array([
        2.0,                           # s5
        1.5,                           # sBP
        -0.5,                          # s3
        2.0 + 1.5,                     # sum_5_bp
        abs(2.0 - 1.5) * 4,            # absdiff_5_bp × gamma_5_bp
        2.0 + (-0.5),                  # sum_5_3
        abs(2.0 - (-0.5)) * 2          # absdiff_5_3 × gamma_5_3
    ])

    assert np.allclose(X_aug[0], expected), f"Feature mismatch!\nExpected: {expected}\nGot: {X_aug[0]}"
    print("\n✓ All feature calculations correct")
    print()

    return True


def test_gamma_scaling():
    """Test that gamma parameters correctly scale the penalty features."""
    print("=" * 80)
    print("TEST 3: Gamma Scaling")
    print("=" * 80)

    # Same input, different gamma values
    X = np.array([[2.0, 1.0, -0.5]])

    # Low gamma
    transformer_low = BothEndsStrongTransformer(gamma_5_bp=1, gamma_5_3=1)
    X_low = transformer_low.transform(X)

    # High gamma
    transformer_high = BothEndsStrongTransformer(gamma_5_bp=8, gamma_5_3=4)
    X_high = transformer_high.transform(X)

    print(f"Input: s5={X[0,0]}, sBP={X[0,1]}, s3={X[0,2]}")
    print(f"\nWith gamma_5_bp=1, gamma_5_3=1:")
    print(f"  absdiff_5_bp = {X_low[0, 4]:.3f} (expected: {abs(2.0-1.0)*1:.3f})")
    print(f"  absdiff_5_3 = {X_low[0, 6]:.3f} (expected: {abs(2.0-(-0.5))*1:.3f})")

    print(f"\nWith gamma_5_bp=8, gamma_5_3=4:")
    print(f"  absdiff_5_bp = {X_high[0, 4]:.3f} (expected: {abs(2.0-1.0)*8:.3f})")
    print(f"  absdiff_5_3 = {X_high[0, 6]:.3f} (expected: {abs(2.0-(-0.5))*4:.3f})")

    # Verify scaling
    assert np.isclose(X_low[0, 4], abs(2.0-1.0)*1), "gamma_5_bp=1 scaling failed"
    assert np.isclose(X_low[0, 6], abs(2.0-(-0.5))*1), "gamma_5_3=1 scaling failed"
    assert np.isclose(X_high[0, 4], abs(2.0-1.0)*8), "gamma_5_bp=8 scaling failed"
    assert np.isclose(X_high[0, 6], abs(2.0-(-0.5))*4), "gamma_5_3=4 scaling failed"

    print("\n✓ Gamma scaling works correctly")
    print()

    return True


def test_one_end_strong_detection():
    """Test that the features can detect one-end-strong patterns."""
    print("=" * 80)
    print("TEST 4: One-End-Strong Pattern Detection")
    print("=" * 80)

    transformer = BothEndsStrongTransformer(gamma_5_bp=4, gamma_5_3=2)

    # Balanced strong pattern (genuine U12)
    X_balanced = np.array([[3.0, 2.5, 2.0]])  # All positive, close values

    # One-end-strong pattern (false positive)
    X_imbalanced = np.array([[5.0, -1.0, -2.0]])  # Only 5'SS strong

    X_balanced_aug = transformer.transform(X_balanced)
    X_imbalanced_aug = transformer.transform(X_imbalanced)

    print("Balanced pattern (genuine U12): s5=3.0, sBP=2.5, s3=2.0")
    print(f"  sum_5_bp = {X_balanced_aug[0, 3]:.3f}")
    print(f"  absdiff_5_bp × γ = {X_balanced_aug[0, 4]:.3f}")
    print(f"  → Small imbalance penalty (good!)")

    print("\nOne-end-strong pattern (FP): s5=5.0, sBP=-1.0, s3=-2.0")
    print(f"  sum_5_bp = {X_imbalanced_aug[0, 3]:.3f}")
    print(f"  absdiff_5_bp × γ = {X_imbalanced_aug[0, 4]:.3f}")
    print(f"  → Large imbalance penalty (should reject!)")

    # Verify that imbalanced pattern gets larger penalty
    balanced_penalty = X_balanced_aug[0, 4]
    imbalanced_penalty = X_imbalanced_aug[0, 4]

    assert imbalanced_penalty > balanced_penalty * 5, \
        f"Imbalanced penalty ({imbalanced_penalty:.2f}) should be much larger than balanced ({balanced_penalty:.2f})"

    print(f"\n✓ One-end-strong patterns get {imbalanced_penalty/balanced_penalty:.1f}× larger penalty")
    print()

    return True


def test_sklearn_compatibility():
    """Test that transformer is compatible with sklearn pipelines."""
    print("=" * 80)
    print("TEST 5: Sklearn Pipeline Compatibility")
    print("=" * 80)

    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import RobustScaler
    from sklearn.svm import LinearSVC

    # Create pipeline matching our training setup
    pipeline = Pipeline([
        ('scale', RobustScaler(with_centering=False, with_scaling=True)),
        ('augment', BothEndsStrongTransformer(gamma_5_bp=4, gamma_5_3=2)),
        ('svc', LinearSVC(class_weight='balanced', random_state=42, max_iter=1000))
    ])

    # Generate simple synthetic data
    np.random.seed(42)
    X = np.random.randn(100, 3)
    y = np.random.randint(0, 2, size=100)

    try:
        pipeline.fit(X, y)
        predictions = pipeline.predict(X)
        print(f"✓ Pipeline fit successful")
        print(f"✓ Predictions shape: {predictions.shape}")
        print(f"✓ Unique predictions: {np.unique(predictions)}")
        print()
        return True
    except Exception as e:
        print(f"✗ Pipeline failed: {e}")
        print()
        return False


def test_optional_features():
    """Test optional min and total_imbalance features."""
    print("=" * 80)
    print("TEST 6: Optional Features")
    print("=" * 80)

    # Default: 7 features
    transformer_default = BothEndsStrongTransformer()
    X = np.array([[2.0, 1.5, -0.5]])
    X_default = transformer_default.transform(X)
    print(f"Default features: {X_default.shape[1]} (expected: 7)")
    assert X_default.shape[1] == 7

    # With min: 8 features
    transformer_min = BothEndsStrongTransformer(include_min=True)
    X_min = transformer_min.transform(X)
    print(f"With include_min: {X_min.shape[1]} (expected: 8)")
    assert X_min.shape[1] == 8

    # With total imbalance: 8 features
    transformer_imb = BothEndsStrongTransformer(include_total_imbalance=True)
    X_imb = transformer_imb.transform(X)
    print(f"With include_total_imbalance: {X_imb.shape[1]} (expected: 8)")
    assert X_imb.shape[1] == 8

    # With both: 9 features
    transformer_both = BothEndsStrongTransformer(include_min=True, include_total_imbalance=True)
    X_both = transformer_both.transform(X)
    print(f"With both optional features: {X_both.shape[1]} (expected: 9)")
    assert X_both.shape[1] == 9

    print("✓ All optional feature configurations work")
    print()

    return True


def main():
    """Run all smoke tests."""
    print("\n" + "=" * 80)
    print("BOTHENDSSTRONGTRANSFORMER SMOKE TESTS")
    print("=" * 80)
    print()

    tests = [
        ("Basic Transformation", test_basic_transformation),
        ("Feature Values", test_feature_values),
        ("Gamma Scaling", test_gamma_scaling),
        ("One-End-Strong Detection", test_one_end_strong_detection),
        ("Sklearn Compatibility", test_sklearn_compatibility),
        ("Optional Features", test_optional_features),
    ]

    results = []
    for name, test_func in tests:
        try:
            success = test_func()
            results.append((name, success))
        except Exception as e:
            print(f"✗ TEST FAILED: {name}")
            print(f"  Error: {e}")
            import traceback
            traceback.print_exc()
            results.append((name, False))

    # Summary
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)

    passed = sum(1 for _, success in results if success)
    total = len(results)

    for name, success in results:
        status = "✓ PASS" if success else "✗ FAIL"
        print(f"{status}: {name}")

    print(f"\nTotal: {passed}/{total} tests passed")

    if passed == total:
        print("\n🎉 All smoke tests passed! BothEndsStrongTransformer is working correctly.")
        return 0
    else:
        print(f"\n⚠️  {total - passed} test(s) failed. Please review the errors above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
