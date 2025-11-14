#!/usr/bin/env python3
"""
Smoke test for BothEndsStrongTransformer implementation.

Quick validation that the min/max feature augmentation is working correctly
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
    """Test that 3D → 5D/7D transformation works."""
    print("=" * 80)
    print("TEST 1: Basic Transformation")
    print("=" * 80)

    # Default: min features only (5D)
    transformer = BothEndsStrongTransformer(include_max=False)

    # Single sample
    X = np.array([[2.0, 1.5, -0.5]])
    X_aug = transformer.fit_transform(X)

    print(f"Input shape: {X.shape} (expected: (1, 3))")
    print(f"Output shape (min only): {X_aug.shape} (expected: (1, 5))")

    assert X_aug.shape == (1, 5), f"Expected (1, 5), got {X_aug.shape}"
    print("✓ Shape transformation correct: 3D → 5D (min features)")

    # With max features (7D)
    transformer_max = BothEndsStrongTransformer(include_max=True)
    X_aug_max = transformer_max.transform(X)
    print(f"Output shape (min+max): {X_aug_max.shape} (expected: (1, 7))")
    assert X_aug_max.shape == (1, 7), f"Expected (1, 7), got {X_aug_max.shape}"
    print("✓ Shape transformation correct: 3D → 7D (min+max features)")

    # Multiple samples
    X_multi = np.random.randn(100, 3)
    X_multi_aug = transformer.transform(X_multi)
    assert X_multi_aug.shape == (100, 5), f"Expected (100, 5), got {X_multi_aug.shape}"
    print("✓ Batch transformation works")
    print()

    return True


def test_feature_values():
    """Test that feature calculations are correct."""
    print("=" * 80)
    print("TEST 2: Feature Value Calculations (Min/Max)")
    print("=" * 80)

    transformer = BothEndsStrongTransformer(include_max=True)

    # Test case: s5=2.0, sBP=1.5, s3=-0.5
    X = np.array([[2.0, 1.5, -0.5]])
    X_aug = transformer.transform(X)

    print(f"Input: s5={X[0,0]}, sBP={X[0,1]}, s3={X[0,2]}")
    print(f"\nExpected output features:")
    print(f"  [0] s5 (pass-through):        {X[0,0]:.3f}")
    print(f"  [1] sBP (pass-through):       {X[0,1]:.3f}")
    print(f"  [2] s3 (pass-through):        {X[0,2]:.3f}")

    # min(a, b) = 0.5 * ((a + b) - |a - b|)
    # max(a, b) = 0.5 * ((a + b) + |a - b|)
    min_5_bp = 0.5 * ((2.0 + 1.5) - abs(2.0 - 1.5))
    min_5_3 = 0.5 * ((2.0 + (-0.5)) - abs(2.0 - (-0.5)))
    max_5_bp = 0.5 * ((2.0 + 1.5) + abs(2.0 - 1.5))
    max_5_3 = 0.5 * ((2.0 + (-0.5)) + abs(2.0 - (-0.5)))

    print(f"  [3] min(s5,sBP):              {min_5_bp:.3f}")
    print(f"  [4] min(s5,s3):               {min_5_3:.3f}")
    print(f"  [5] max(s5,sBP):              {max_5_bp:.3f}")
    print(f"  [6] max(s5,s3):               {max_5_3:.3f}")

    print(f"\nActual output features:")
    for i, val in enumerate(X_aug[0]):
        print(f"  [{i}] {val:.3f}")

    # Verify calculations
    expected = np.array([
        2.0,                           # s5
        1.5,                           # sBP
        -0.5,                          # s3
        min_5_bp,                      # min(s5, sBP)
        min_5_3,                       # min(s5, s3)
        max_5_bp,                      # max(s5, sBP)
        max_5_3                        # max(s5, s3)
    ])

    assert np.allclose(X_aug[0], expected), f"Feature mismatch!\nExpected: {expected}\nGot: {X_aug[0]}"
    print("\n✓ All feature calculations correct")
    print()

    return True


def test_min_max_properties():
    """Test that min/max features have correct mathematical properties."""
    print("=" * 80)
    print("TEST 3: Min/Max Properties")
    print("=" * 80)

    transformer = BothEndsStrongTransformer(include_max=True)

    # Test various inputs
    test_cases = [
        ([2.0, 1.5], "Close values"),
        ([3.0, 3.0], "Equal values"),
        ([5.0, -1.0], "Large difference"),
    ]

    for (a, b), description in test_cases:
        X = np.array([[a, b, 0.0]])  # s3 doesn't matter for this test
        X_aug = transformer.transform(X)

        min_val = X_aug[0, 3]  # min(s5, sBP)
        max_val = X_aug[0, 5]  # max(s5, sBP)

        print(f"\n{description}: s5={a}, sBP={b}")
        print(f"  min(s5,sBP) = {min_val:.3f}")
        print(f"  max(s5,sBP) = {max_val:.3f}")

        # Verify mathematical properties
        assert min_val <= min(a, b) + 1e-10, "min should be ≤ minimum input"
        assert max_val >= max(a, b) - 1e-10, "max should be ≥ maximum input"
        assert np.isclose(min_val, min(a, b)), f"min formula incorrect: {min_val} vs {min(a, b)}"
        assert np.isclose(max_val, max(a, b)), f"max formula incorrect: {max_val} vs {max(a, b)}"

    print("\n✓ Min/max formulas are mathematically correct")
    print()

    return True


def test_one_end_strong_detection():
    """Test that min features can detect one-end-strong patterns."""
    print("=" * 80)
    print("TEST 4: One-End-Strong Pattern Detection")
    print("=" * 80)

    transformer = BothEndsStrongTransformer(include_max=False)

    # Balanced strong pattern (genuine U12)
    X_balanced = np.array([[3.0, 2.5, 2.0]])  # All positive, close values

    # One-end-strong pattern (false positive)
    X_imbalanced = np.array([[5.0, -1.0, -2.0]])  # Only 5'SS strong

    X_balanced_aug = transformer.transform(X_balanced)
    X_imbalanced_aug = transformer.transform(X_imbalanced)

    print("Balanced pattern (genuine U12): s5=3.0, sBP=2.5, s3=2.0")
    print(f"  min(s5,sBP) = {X_balanced_aug[0, 3]:.3f}")
    print(f"  → High min value (both strong!)")

    print("\nOne-end-strong pattern (FP): s5=5.0, sBP=-1.0, s3=-2.0")
    print(f"  min(s5,sBP) = {X_imbalanced_aug[0, 3]:.3f}")
    print(f"  → Low/negative min value (one weak!)")

    # Verify that imbalanced pattern has much lower min
    balanced_min = X_balanced_aug[0, 3]
    imbalanced_min = X_imbalanced_aug[0, 3]

    assert balanced_min > 2.0, f"Balanced min should be high (both strong)"
    assert imbalanced_min < 0, f"Imbalanced min should be negative (one weak)"

    print(f"\n✓ Min features correctly distinguish balanced vs one-end-strong")
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
        ('augment', BothEndsStrongTransformer(include_max=False)),
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


def test_feature_names():
    """Test get_feature_names_out() method."""
    print("=" * 80)
    print("TEST 6: Feature Names")
    print("=" * 80)

    # Min only
    transformer_min = BothEndsStrongTransformer(include_max=False)
    names_min = transformer_min.get_feature_names_out()
    print(f"Features (min only): {list(names_min)}")
    assert len(names_min) == 5
    assert list(names_min) == ['s5', 'sBP', 's3', 'min_5_bp', 'min_5_3']
    print("✓ Min-only feature names correct")

    # Min + max
    transformer_max = BothEndsStrongTransformer(include_max=True)
    names_max = transformer_max.get_feature_names_out()
    print(f"Features (min+max): {list(names_max)}")
    assert len(names_max) == 7
    assert list(names_max) == ['s5', 'sBP', 's3', 'min_5_bp', 'min_5_3', 'max_5_bp', 'max_5_3']
    print("✓ Min+max feature names correct")
    print()

    return True


def main():
    """Run all smoke tests."""
    print("\n" + "=" * 80)
    print("BOTHENDSSTRONGTRANSFORMER SMOKE TESTS (Min/Max Approach)")
    print("=" * 80)
    print()

    tests = [
        ("Basic Transformation", test_basic_transformation),
        ("Feature Values", test_feature_values),
        ("Min/Max Properties", test_min_max_properties),
        ("One-End-Strong Detection", test_one_end_strong_detection),
        ("Sklearn Compatibility", test_sklearn_compatibility),
        ("Feature Names", test_feature_names),
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
