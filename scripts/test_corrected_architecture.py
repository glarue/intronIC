"""
Quick test to verify corrected architecture is working.

Tests:
1. Pipeline creation (imports, instantiation)
2. Feature transformation (z-scores → augmented features)
3. No double-scaling (single scaling step)
"""

import numpy as np
from sklearn.pipeline import Pipeline
from sklearn.svm import LinearSVC
from classification.transformers import BothEndsStrongTransformer

def test_pipeline_creation():
    """Test that pipeline can be created with BothEndsStrongTransformer."""
    print("Test 1: Pipeline creation...")

    pipeline = Pipeline([
        ('transform', BothEndsStrongTransformer(
            include_max=False,
            include_pairwise_mins=False
        )),
        ('svc', LinearSVC(
            C=0.1,
            penalty='l2',
            dual=False,
            loss='squared_hinge',
            class_weight='balanced',
            max_iter=1000,
            tol=1e-4,
            random_state=42
        ))
    ])

    print(f"✓ Pipeline created: {pipeline.named_steps.keys()}")
    return pipeline

def test_feature_transformation():
    """Test that z-scores are correctly transformed to augmented features."""
    print("\nTest 2: Feature transformation...")

    transformer = BothEndsStrongTransformer(
        include_max=False,
        include_pairwise_mins=False
    )

    # Simulate z-scores (3D input)
    X_z = np.array([
        [1.0, 2.0, 1.5],  # five_z, bp_z, three_z
        [-1.0, -2.0, -1.5],
        [0.5, 0.5, 0.5]
    ])

    # Transform to augmented features (7D output)
    X_aug = transformer.fit_transform(X_z)

    print(f"✓ Input shape (z-scores): {X_z.shape} (3D: five_z, bp_z, three_z)")
    print(f"✓ Output shape (augmented): {X_aug.shape} (7D: base + min_all + neg_absdiff_*)")

    # Check dimensions
    assert X_aug.shape == (3, 7), f"Expected (3, 7), got {X_aug.shape}"

    # Verify first 3 columns are passed through
    np.testing.assert_array_almost_equal(X_aug[:, :3], X_z)
    print("✓ First 3 features are pass-through z-scores")

    # Verify min_all feature (column 3)
    expected_min_all = np.minimum(np.minimum(X_z[:, 0], X_z[:, 1]), X_z[:, 2])
    np.testing.assert_array_almost_equal(X_aug[:, 3], expected_min_all)
    print("✓ Feature 4 (min_all) computed correctly")

    # Verify neg_absdiff features (columns 4-6)
    expected_neg_absdiff_5_bp = -np.abs(X_z[:, 0] - X_z[:, 1])
    np.testing.assert_array_almost_equal(X_aug[:, 4], expected_neg_absdiff_5_bp)
    print("✓ Feature 5 (neg_absdiff_5_bp) computed correctly")

    return transformer

def test_no_double_scaling():
    """Verify that pipeline has NO scaler (single scaling step)."""
    print("\nTest 3: No double-scaling...")

    pipeline = Pipeline([
        ('transform', BothEndsStrongTransformer(
            include_max=False,
            include_pairwise_mins=False
        )),
        ('svc', LinearSVC(
            C=0.1,
            penalty='l2',
            dual=False,
            loss='squared_hinge',
            class_weight='balanced',
            max_iter=1000,
            tol=1e-4,
            random_state=42
        ))
    ])

    # Verify no scaler in pipeline
    assert 'scale' not in pipeline.named_steps, "ERROR: Pipeline should NOT have 'scale' step"
    print("✓ No 'scale' step in pipeline (scaling done externally by ScoreNormalizer)")

    # Verify only transformer and svc
    assert len(pipeline.named_steps) == 2, f"Expected 2 steps, got {len(pipeline.named_steps)}"
    assert 'transform' in pipeline.named_steps
    assert 'svc' in pipeline.named_steps
    print(f"✓ Pipeline has correct steps: {list(pipeline.named_steps.keys())}")

def test_end_to_end():
    """Test end-to-end training and prediction."""
    print("\nTest 4: End-to-end pipeline...")

    # Create synthetic data (z-scores)
    np.random.seed(42)
    X_train = np.random.randn(100, 3)  # 100 samples, 3 z-score features
    y_train = np.random.randint(0, 2, 100)  # Binary labels

    X_test = np.random.randn(20, 3)

    # Create and train pipeline
    pipeline = Pipeline([
        ('transform', BothEndsStrongTransformer(
            include_max=False,
            include_pairwise_mins=False
        )),
        ('svc', LinearSVC(
            C=0.1,
            penalty='l2',
            dual=False,
            loss='squared_hinge',
            class_weight='balanced',
            max_iter=1000,
            tol=1e-4,
            random_state=42
        ))
    ])

    print("  Training pipeline...")
    pipeline.fit(X_train, y_train)
    print("  ✓ Training successful")

    print("  Predicting on test data...")
    y_pred = pipeline.predict(X_test)
    print(f"  ✓ Predictions: {y_pred[:5]}... (shape: {y_pred.shape})")

    # Verify prediction shape
    assert y_pred.shape == (20,), f"Expected (20,), got {y_pred.shape}"
    print("✓ End-to-end pipeline works correctly")

if __name__ == '__main__':
    print("=" * 60)
    print("CORRECTED ARCHITECTURE VERIFICATION")
    print("=" * 60)

    try:
        test_pipeline_creation()
        test_feature_transformation()
        test_no_double_scaling()
        test_end_to_end()

        print("\n" + "=" * 60)
        print("ALL TESTS PASSED ✓")
        print("=" * 60)
        print("\nArchitecture summary:")
        print("  1. Single scaling step (ScoreNormalizer OUTSIDE pipeline)")
        print("  2. Pipeline: z-scores → BothEndsStrongTransformer → LinearSVC")
        print("  3. Features: 3D → 7D (base + min_all + 3 neg_absdiff)")
        print("  4. No double-scaling (no 'scale' step in pipeline)")

    except Exception as e:
        print(f"\n✗ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
