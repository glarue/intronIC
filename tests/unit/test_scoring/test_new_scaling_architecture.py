"""
Tests for new scaling architecture (2025 redesign).

Tests the core components that fix cross-species false positives:
- SymmetricClipper: Aggressive outlier control in z-space
- SaturatingTransform: Optional log compression for extremes
- Integration: Full pipeline behavior

Design: SCALER_ARCHITECTURE_REVIEW.md, SCALING_REDESIGN_PLAN.md
"""

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

from intronIC.classification.clipping import SymmetricClipper
from intronIC.classification.saturating import SaturatingTransform
from intronIC.scoring.normalizer import ZeroAnchoredRobustScaler

# ==============================================================================
# SymmetricClipper Tests
# ==============================================================================


def test_symmetric_clipper_preserves_zero():
    """Zero values should pass through unchanged."""
    clipper = SymmetricClipper(quantile=0.95)
    X = np.array([[0.0, 0.0, 0.0]])

    clipper.fit(X)
    X_clipped = clipper.transform(X)

    assert_array_almost_equal(X_clipped, [[0.0, 0.0, 0.0]])


def test_symmetric_clipper_clips_extreme_values():
    """Extreme values should be clipped to quantile threshold."""
    # Training data: mostly moderate values + one extreme
    X_train = np.array(
        [
            [1.0, 2.0, 3.0],
            [1.5, 2.5, 3.5],
            [2.0, 3.0, 4.0],
            [0.5, 1.0, 1.5],
            # Extreme outlier
            [10.0, 15.0, 20.0],
        ]
    )

    # Use domain_adaptive=False to test fixed-cap behavior
    clipper = SymmetricClipper(quantile=0.8, domain_adaptive=False)
    clipper.fit(X_train)

    # Apply to extreme value (should be clipped)
    X_extreme = np.array([[100.0, 150.0, 200.0]])
    X_clipped = clipper.transform(X_extreme)

    # Each feature should be clipped to its learned cap
    # Caps are at 80th percentile of absolute values
    # For first feature: [1.0, 1.5, 2.0, 0.5, 10.0] → 80th % ≈ 2.0
    # Extreme value (100.0) should be clipped to cap
    assert X_clipped[0, 0] < 100.0  # Definitely clipped
    assert X_clipped[0, 0] > 0.0  # Not zeroed


def test_symmetric_clipper_clips_both_positive_and_negative():
    """Clipping should be symmetric for positive and negative values."""
    X_train = np.array([[1.0], [2.0], [3.0], [4.0], [5.0]])

    clipper = SymmetricClipper(quantile=0.8)
    clipper.fit(X_train)

    # Test symmetric extremes
    X_test = np.array([[100.0], [-100.0]])
    X_clipped = clipper.transform(X_test)

    # Both should be clipped to same absolute value (symmetric)
    assert abs(X_clipped[0, 0]) == abs(X_clipped[1, 0])
    assert X_clipped[0, 0] > 0  # Positive stays positive
    assert X_clipped[1, 0] < 0  # Negative stays negative


def test_symmetric_clipper_cross_species_use_case():
    """
    Test the actual cross-species scenario: C. elegans extreme value.

    Scenario with domain_adaptive=False (legacy behavior):
    - Train on human z-scores (moderate range)
    - Apply to C. elegans with z=11.45 from composition bias
    - Should clip to ~3-4σ, preventing FP

    NOTE: With domain_adaptive=True (default), the clipper recomputes
    caps from the input data, so this test uses domain_adaptive=False.
    """
    # Human training z-scores (typical range: -3 to +3)
    human_z = np.array(
        [
            [-2.0, 1.5, 0.8],
            [1.0, -1.2, 2.1],
            [0.5, 0.8, -0.5],
            [2.5, -2.3, 1.8],
            [1.8, 1.9, 2.8],
        ]
    )

    # Use domain_adaptive=False to test fixed-cap behavior
    clipper = SymmetricClipper(quantile=0.975, domain_adaptive=False)
    clipper.fit(human_z)

    # C. elegans with extreme 3'SS z-score
    celegans_z = np.array([[-0.15, 0.30, 11.45]])
    celegans_clipped = clipper.transform(celegans_z)

    # The extreme value should be clipped
    assert celegans_clipped[0, 2] < 11.45, "Extreme value should be clipped"
    assert celegans_clipped[0, 2] > 0, "Should remain positive"
    # Clipped value should be reasonable (roughly 3-4σ)
    assert 2.5 < celegans_clipped[0, 2] < 5.0, (
        f"Expected ~3-4σ, got {celegans_clipped[0, 2]}"
    )

    # Moderate values should pass through
    assert celegans_clipped[0, 0] == pytest.approx(-0.15, abs=0.01)
    assert celegans_clipped[0, 1] == pytest.approx(0.30, abs=0.01)


# ==============================================================================
# SaturatingTransform Tests
# ==============================================================================


def test_saturating_transform_disabled_is_identity():
    """When disabled, transform should be identity (pass through)."""
    saturate = SaturatingTransform(enabled=False)
    X = np.array([[1.5, -2.0, 10.0]])

    X_out = saturate.fit_transform(X)

    assert_array_almost_equal(X_out, X)


def test_saturating_transform_preserves_zero():
    """f(0) = 0 (zero stays zero)."""
    saturate = SaturatingTransform(enabled=True)
    X = np.array([[0.0, 0.0, 0.0]])

    X_out = saturate.fit_transform(X)

    assert_array_almost_equal(X_out, [[0.0, 0.0, 0.0]])


def test_saturating_transform_preserves_sign():
    """Positive stays positive, negative stays negative."""
    saturate = SaturatingTransform(enabled=True)
    X = np.array([[5.0, -5.0, 0.0]])

    X_out = saturate.fit_transform(X)

    assert X_out[0, 0] > 0, "Positive should stay positive"
    assert X_out[0, 1] < 0, "Negative should stay negative"
    assert X_out[0, 2] == 0, "Zero should stay zero"


def test_saturating_transform_known_values():
    """Test f(z) = sign(z) * log(1 + |z|) with known values."""
    saturate = SaturatingTransform(enabled=True)

    # Test specific values we can verify by hand
    # f(0) = 0
    # f(1) = log(2) ≈ 0.693
    # f(2) = log(3) ≈ 1.099
    # f(4) = log(5) ≈ 1.609
    # f(-2) = -log(3) ≈ -1.099
    X = np.array([[0.0, 1.0, 2.0, 4.0, -2.0]])

    X_out = saturate.fit_transform(X)

    assert X_out[0, 0] == pytest.approx(0.0, abs=0.01)
    assert X_out[0, 1] == pytest.approx(0.693, abs=0.01)
    assert X_out[0, 2] == pytest.approx(1.099, abs=0.01)
    assert X_out[0, 3] == pytest.approx(1.609, abs=0.01)
    assert X_out[0, 4] == pytest.approx(-1.099, abs=0.01)


def test_saturating_transform_compresses_extremes():
    """
    Test key property: Compresses gap between extreme and moderate values.

    Example from design doc:
    Before: 11.45 - 4 = 7.45 (huge gap)
    After:  2.52 - 1.61 = 0.91 (moderate gap)
    """
    saturate = SaturatingTransform(enabled=True)
    X = np.array([[4.0, 11.45]])

    X_out = saturate.fit_transform(X)

    # Original gap
    original_gap = 11.45 - 4.0
    assert original_gap == pytest.approx(7.45, abs=0.01)

    # Transformed gap (should be much smaller)
    transformed_gap = X_out[0, 1] - X_out[0, 0]

    # f(4) ≈ 1.609, f(11.45) ≈ 2.519
    assert X_out[0, 0] == pytest.approx(1.609, abs=0.01)
    assert X_out[0, 1] == pytest.approx(2.519, abs=0.01)

    # Gap should be compressed (from 7.45 to ~0.91)
    assert transformed_gap < 1.5, f"Gap should be compressed, got {transformed_gap}"
    assert transformed_gap == pytest.approx(0.91, abs=0.05)


def test_saturating_transform_inverse():
    """Inverse transform should recover original values."""
    saturate = SaturatingTransform(enabled=True)
    X_original = np.array([[4.0, -2.0, 11.45]])

    # Forward transform
    X_saturated = saturate.fit_transform(X_original)

    # Inverse transform
    X_recovered = saturate.inverse_transform(X_saturated)

    # Should recover original (within numerical precision)
    assert_array_almost_equal(X_recovered, X_original, decimal=10)


# ==============================================================================
# Integration Tests - Full Pipeline
# ==============================================================================


def test_full_pipeline_cross_species_scenario():
    """
    Integration test: Full pipeline on cross-species use case.

    Pipeline: ZeroAnchoredRobustScaler → SymmetricClipper → SaturatingTransform

    This simulates:
    1. Training on human raw LLRs
    2. Applying to C. elegans with composition bias
    3. Verifying extreme values are controlled
    """
    # Human training data (raw LLRs)
    human_raw = np.array(
        [
            [2.0, 1.5, 1.8],
            [1.8, 1.2, 2.1],
            [2.2, 1.8, 1.5],
            [1.5, 1.3, 2.0],
            [2.1, 1.6, 1.9],
        ]
    )

    # Build pipeline
    from sklearn.pipeline import Pipeline

    pipeline = Pipeline(
        [
            ("scale", ZeroAnchoredRobustScaler()),
            ("clip", SymmetricClipper(quantile=0.975)),
            ("saturate", SaturatingTransform(enabled=True)),
        ]
    )

    # Fit on human data
    pipeline.fit(human_raw)

    # C. elegans with extreme value
    # Assume median scale ≈ 1.8, so 14.46 / 1.8 ≈ 8.0 (z-score)
    # After clipping and saturation, should be much smaller
    celegans_raw = np.array([[1.5, 1.2, 14.46]])

    celegans_transformed = pipeline.transform(celegans_raw)

    # First two features should be moderate (similar to training range)
    assert abs(celegans_transformed[0, 0]) < 3.0
    assert abs(celegans_transformed[0, 1]) < 3.0

    # Third feature (extreme) should be controlled
    # After scaling: ~8σ
    # After clipping: ~3-4σ
    # After saturation: ~1.6
    assert celegans_transformed[0, 2] < 3.0, (
        f"Extreme value not controlled: {celegans_transformed[0, 2]}"
    )


def test_pipeline_preserves_zero_semantic():
    """
    Critical property: Zero LLR should stay zero through entire pipeline.

    s=0 means "U12≈U2" - this semantic meaning must be preserved.
    """
    # Training data with zeros
    X_train = np.array(
        [[0.0, 1.5, 2.0], [1.0, 0.0, 1.8], [2.0, 1.2, 0.0], [1.5, 1.8, 1.5]]
    )

    from sklearn.pipeline import Pipeline

    pipeline = Pipeline(
        [
            ("scale", ZeroAnchoredRobustScaler()),
            ("clip", SymmetricClipper(quantile=0.95)),
            ("saturate", SaturatingTransform(enabled=True)),
        ]
    )

    pipeline.fit(X_train)

    # Test that zero stays zero
    X_test = np.array([[0.0, 0.0, 0.0]])
    X_out = pipeline.transform(X_test)

    assert_array_almost_equal(X_out, [[0.0, 0.0, 0.0]], decimal=10)


def test_pipeline_steps_accessible():
    """
    Verify we can extract fitted components from pipeline.

    This is needed by predictor.py to compute z-scores for output files.
    """
    from sklearn.pipeline import Pipeline

    X_train = np.array([[1.0, 2.0, 3.0], [1.5, 2.5, 3.5]])

    pipeline = Pipeline(
        [
            ("scale", ZeroAnchoredRobustScaler()),
            ("clip", SymmetricClipper(quantile=0.95)),
            ("saturate", SaturatingTransform(enabled=True)),
        ]
    )

    pipeline.fit(X_train)

    # Should be able to extract each component
    scaler = pipeline.named_steps["scale"]
    clipper = pipeline.named_steps["clip"]
    saturator = pipeline.named_steps["saturate"]

    # Should be fitted
    assert scaler.scales_ is not None, "Scaler should be fitted"
    assert clipper.caps_ is not None, "Clipper should be fitted"

    # Should be able to use them independently
    X_test = np.array([[2.0, 3.0, 4.0]])
    z_scores = scaler.transform(X_test)  # What predictor.py does

    assert z_scores.shape == X_test.shape
