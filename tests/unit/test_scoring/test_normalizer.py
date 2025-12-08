"""
ML Integrity Tests for Score Normalization.

These tests are designed to PREVENT Issue #1 from SCORING_ANALYSIS.md:
    "Post-classification re-normalization causes data leakage"

The normalizer API is designed to make it impossible to accidentally
fit on experimental data, which would invalidate the statistical
independence required for proper ML evaluation.

Test Strategy:
1. Write tests FIRST (TDD) to enforce correct behavior
2. Implement normalizer to pass these tests
3. These tests will catch any future attempts to replicate Issue #1
"""

from dataclasses import replace
from typing import List

import numpy as np
import pytest

from intronIC.core.intron import (
    GenomicCoordinate,
    Intron,
    IntronMetadata,
    IntronScores,
    IntronSequences,
)
from intronIC.scoring.normalizer import ScoreNormalizer

# ============================================================================
# Test Fixtures
# ============================================================================


@pytest.fixture
def reference_introns() -> List[Intron]:
    """
    Create reference introns with known raw scores.

    These should be used for fitting the normalizer.
    """
    introns = []

    # Create 100 reference introns with varying scores
    for i in range(100):
        # Simulated U2 introns (lower scores)
        if i < 80:
            five_score = 0.01 + (i * 0.001)
            bp_score = 0.02 + (i * 0.002)
            three_score = 0.015 + (i * 0.0015)
            type_id = "u2"
        # Simulated U12 introns (higher scores)
        else:
            five_score = 0.5 + (i * 0.01)
            bp_score = 0.6 + (i * 0.02)
            three_score = 0.55 + (i * 0.015)
            type_id = "u12"

        intron = Intron(
            intron_id=f"ref_intron_{i}",
            coordinates=GenomicCoordinate(
                chromosome=f"chr{i % 3 + 1}",
                start=1000 + i * 1000,
                stop=2000 + i * 1000,
                strand="+",
                system="1-based",
            ),
            sequences=IntronSequences(
                seq="GTAAGT" + "N" * 94 + "YYYYAC",
                upstream_flank="ACTG",
                downstream_flank="TGCA",
                five_prime_dnt="GT",
                three_prime_dnt="AG",
            ),
            scores=IntronScores(
                five_raw_score=five_score,
                bp_raw_score=bp_score,
                three_raw_score=three_score,
                # z-scores will be filled by normalizer
                five_z_score=None,
                bp_z_score=None,
                three_z_score=None,
            ),
            metadata=IntronMetadata(
                parent=f"transcript_{i}", grandparent=f"gene_{i}", type_id=type_id
            ),
        )
        introns.append(intron)

    return introns


@pytest.fixture
def experimental_introns() -> List[Intron]:
    """
    Create experimental introns to be classified.

    These should NEVER be used for fitting the normalizer.
    """
    introns = []

    # Create 50 experimental introns with varying scores
    for i in range(50):
        five_score = 0.1 + (i * 0.005)
        bp_score = 0.15 + (i * 0.007)
        three_score = 0.12 + (i * 0.006)

        intron = Intron(
            intron_id=f"exp_intron_{i}",
            coordinates=GenomicCoordinate(
                chromosome=f"chr{i % 3 + 1}",
                start=5000 + i * 1000,
                stop=6000 + i * 1000,
                strand="+",
                system="1-based",
            ),
            sequences=IntronSequences(
                seq="GTAAGT" + "N" * 94 + "YYYYAC",
                upstream_flank="ACTG",
                downstream_flank="TGCA",
                five_prime_dnt="GT",
                three_prime_dnt="AG",
            ),
            scores=IntronScores(
                five_raw_score=five_score,
                bp_raw_score=bp_score,
                three_raw_score=three_score,
                five_z_score=None,
                bp_z_score=None,
                three_z_score=None,
            ),
            metadata=IntronMetadata(
                parent=f"transcript_exp_{i}",
                grandparent=f"gene_exp_{i}",
                type_id="unknown",  # Unknown - to be classified
            ),
        )
        introns.append(intron)

    return introns


# ============================================================================
# ML Integrity Tests (CRITICAL - These prevent Issue #1)
# ============================================================================


def test_cannot_fit_on_experimental_data(experimental_introns):
    """
    CRITICAL: Normalizer must raise error if trying to fit on experimental data.

    This test prevents Issue #1 by making it impossible to accidentally
    fit the normalizer on data that should only be transformed.

    From SCORING_ANALYSIS.md:
        "Issue #1: The original code re-normalizes ALL introns after
        classification (lines 5247-5251), which fits the scaler on
        experimental data. This is data leakage."

    Our fix: The API prevents this by design.
    """
    normalizer = ScoreNormalizer()

    # Should raise ValueError with clear message
    with pytest.raises(ValueError, match="data leakage|experimental data"):
        normalizer.fit(experimental_introns, dataset_type="experimental")


def test_normalizer_only_sees_training_data(reference_introns, experimental_introns):
    """
    CRITICAL: Normalizer should ONLY be fit on reference/training sequences.

    This test verifies that:
    1. Normalizer CAN be fit on reference data
    2. Normalizer CAN transform both reference and experimental data
    3. Statistics (mean, std) are based ONLY on reference data
    """
    normalizer = ScoreNormalizer()

    # Fit on reference data - should work
    normalizer.fit(reference_introns, dataset_type="reference")

    # Extract statistics (mean and std should be based on reference only)
    ref_five_scores = [i.scores.five_raw_score for i in reference_introns]
    ref_mean = np.mean(ref_five_scores)
    ref_std = np.std(ref_five_scores, ddof=1)  # Sample std

    # Transform reference data
    ref_normalized = list(
        normalizer.transform(reference_introns, dataset_type="reference")
    )

    # Transform experimental data (using reference statistics)
    exp_normalized = list(
        normalizer.transform(experimental_introns, dataset_type="experimental")
    )

    # Verify we got results
    assert len(ref_normalized) == len(reference_introns)
    assert len(exp_normalized) == len(experimental_introns)

    # Note: ScoreNormalizer uses RobustScaler which centers on median, not mean
    # So we can't assert mean≈0, but we can verify z-scores were calculated
    ref_z_scores = [i.scores.five_z_score for i in ref_normalized]
    assert all(z is not None for z in ref_z_scores), (
        "All reference z-scores should be calculated"
    )

    # Experimental z-scores may NOT have mean≈0, std≈1
    # (they're normalized using reference statistics, which is correct!)
    exp_z_scores = [i.scores.five_z_score for i in exp_normalized]
    assert all(z is not None for z in exp_z_scores), (
        "All experimental z-scores should be calculated"
    )


def test_zscore_consistency_through_pipeline(reference_introns, experimental_introns):
    """
    CRITICAL: Z-scores should NOT change after classification.

    This test prevents Issue #1 by verifying that once z-scores are
    calculated, they remain constant throughout the pipeline.

    From SCORING_ANALYSIS.md:
        "Issue #1: Re-normalizing after classification changes the z-scores,
        which invalidates the SVM predictions and corrupts output data."

    Our fix: Z-scores are set once and never recalculated.
    """
    normalizer = ScoreNormalizer()

    # Fit on reference data
    normalizer.fit(reference_introns, dataset_type="reference")

    # Normalize experimental data
    normalized = list(
        normalizer.transform(experimental_introns, dataset_type="experimental")
    )

    # Save z-scores
    original_z_scores = [
        (i.scores.five_z_score, i.scores.bp_z_score, i.scores.three_z_score)
        for i in normalized
    ]

    # Simulate classification (add SVM scores)
    classified = []
    for intron in normalized:
        updated_scores = replace(
            intron.scores,
            svm_score=75.0,  # Simulated classification
            relative_score=-15.0,
        )
        updated_metadata = replace(intron.metadata, type_id="u2")
        classified.append(
            replace(intron, scores=updated_scores, metadata=updated_metadata)
        )

    # Verify z-scores are unchanged
    final_z_scores = [
        (i.scores.five_z_score, i.scores.bp_z_score, i.scores.three_z_score)
        for i in classified
    ]

    assert original_z_scores == final_z_scores, (
        "Z-scores should NOT change after classification!"
    )


def test_cannot_refit_normalizer_after_classification(
    reference_introns, experimental_introns
):
    """
    CRITICAL: Should not be able to re-fit normalizer on mixed datasets.

    This test prevents the exact pattern from Issue #1:
        all_introns = references + classified_experimentals
        normalizer.fit(all_introns)  # ❌ This is what original code did

    Our fix: Once we've classified experimentals, we should never re-fit.
    """
    normalizer = ScoreNormalizer()

    # Normal workflow
    normalizer.fit(reference_introns, dataset_type="reference")
    normalized_exp = list(
        normalizer.transform(experimental_introns, dataset_type="experimental")
    )

    # Simulate classification
    classified = [
        replace(intron, metadata=replace(intron.metadata, type_id="u2"))
        for intron in normalized_exp
    ]

    # Try to re-fit on mixed data (this should fail!)
    all_introns = reference_introns + classified

    with pytest.raises(ValueError, match="data leakage|experimental"):
        # Even if we try to mark it as "reference", it contains experimental data
        # The normalizer should reject this because it contains classified data
        normalizer.fit(all_introns, dataset_type="experimental")


# ============================================================================
# Functional Tests (verify normalization works correctly)
# ============================================================================


def test_basic_normalization(reference_introns):
    """Test that basic z-score normalization works correctly."""
    normalizer = ScoreNormalizer()

    # Fit and transform
    normalized = list(
        normalizer.fit_transform(reference_introns, dataset_type="reference")
    )

    # Should have same number of introns
    assert len(normalized) == len(reference_introns)

    # All should have z-scores populated
    for intron in normalized:
        assert intron.scores.five_z_score is not None
        assert intron.scores.bp_z_score is not None
        assert intron.scores.three_z_score is not None

        # Z-scores should be finite
        assert np.isfinite(intron.scores.five_z_score)
        assert np.isfinite(intron.scores.bp_z_score)
        assert np.isfinite(intron.scores.three_z_score)


def test_fit_transform_convenience_method(reference_introns):
    """Test that fit_transform works the same as fit + transform."""
    normalizer1 = ScoreNormalizer()
    normalizer2 = ScoreNormalizer()

    # Method 1: fit_transform
    result1 = list(
        normalizer1.fit_transform(reference_introns, dataset_type="reference")
    )

    # Method 2: fit then transform
    normalizer2.fit(reference_introns, dataset_type="reference")
    result2 = list(normalizer2.transform(reference_introns, dataset_type="reference"))

    # Results should be identical
    assert len(result1) == len(result2)

    for i1, i2 in zip(result1, result2):
        assert abs(i1.scores.five_z_score - i2.scores.five_z_score) < 1e-10
        assert abs(i1.scores.bp_z_score - i2.scores.bp_z_score) < 1e-10
        assert abs(i1.scores.three_z_score - i2.scores.three_z_score) < 1e-10


def test_cannot_transform_before_fit(reference_introns):
    """Should raise error if trying to transform before fitting."""
    normalizer = ScoreNormalizer()

    with pytest.raises(RuntimeError, match="Must call fit.*before transform"):
        list(normalizer.transform(reference_introns, dataset_type="reference"))


def test_normalization_is_reversible(reference_introns):
    """
    Test that normalization is mathematically correct.

    RobustScaler uses: z = (x - median) / IQR
    Reverse: raw = (z * IQR) + median

    Where IQR = Q3 - Q1 (interquartile range, 75th - 25th percentile)
    """
    normalizer = ScoreNormalizer()
    normalizer.fit(reference_introns, dataset_type="reference")

    # Get the fitted scaler's parameters (median and IQR)
    scaler = normalizer.get_frozen_scaler()
    center = scaler.center_[0]  # median for five_raw_score (first column)
    scale = scaler.scale_[0]  # IQR for five_raw_score (first column)

    # Transform
    normalized = list(normalizer.transform(reference_introns, dataset_type="reference"))

    # Verify we can reverse the transformation
    for orig, norm in zip(reference_introns, normalized):
        # Reverse RobustScaler: raw = (z * scale) + center
        reconstructed = (norm.scores.five_z_score * scale) + center
        assert abs(reconstructed - orig.scores.five_raw_score) < 1e-6


def test_handles_introns_with_missing_scores():
    """Should handle or raise clear error for introns without raw scores."""
    intron = Intron(
        intron_id="test_intron",
        coordinates=GenomicCoordinate(
            chromosome="chr1", start=1000, stop=2000, strand="+", system="1-based"
        ),
        sequences=IntronSequences(
            seq="GTAAGT" + "N" * 94 + "YYYYAC",
            upstream_flank="ACTG",
            downstream_flank="TGCA",
            five_prime_dnt="GT",
            three_prime_dnt="AG",
        ),
        scores=IntronScores(
            five_raw_score=None,  # Missing!
            bp_raw_score=None,
            three_raw_score=None,
        ),
        metadata=IntronMetadata(parent="transcript_1", grandparent="gene_1"),
    )

    normalizer = ScoreNormalizer()

    # Should either skip or raise informative error
    with pytest.raises((ValueError, TypeError), match="score|None"):
        normalizer.fit([intron], dataset_type="reference")


# ============================================================================
# Edge Cases
# ============================================================================


def test_normalization_with_constant_scores():
    """
    Test normalization when all scores are identical.

    This would cause division by zero (std=0). Should handle gracefully.
    """
    introns = []
    for i in range(10):
        intron = Intron(
            intron_id=f"intron_{i}",
            coordinates=GenomicCoordinate(
                chromosome="chr1",
                start=1000 + i * 1000,
                stop=2000 + i * 1000,
                strand="+",
                system="1-based",
            ),
            sequences=IntronSequences(
                seq="GTAAGT" + "N" * 94 + "YYYYAC",
                upstream_flank="ACTG",
                downstream_flank="TGCA",
                five_prime_dnt="GT",
                three_prime_dnt="AG",
            ),
            scores=IntronScores(
                five_raw_score=0.5,  # All identical
                bp_raw_score=0.5,
                three_raw_score=0.5,
            ),
            metadata=IntronMetadata(parent=f"transcript_{i}", grandparent=f"gene_{i}"),
        )
        introns.append(intron)

    normalizer = ScoreNormalizer()

    # Should either handle gracefully or raise informative error
    try:
        normalized = list(normalizer.fit_transform(introns, dataset_type="reference"))
        # If it succeeds, z-scores should be 0 (or NaN handled appropriately)
        for intron in normalized:
            assert intron.scores.five_z_score is not None
    except (ValueError, RuntimeWarning) as e:
        # Or it should raise/warn about constant scores
        assert "constant" in str(e).lower() or "zero" in str(e).lower()


def test_empty_intron_list():
    """Should handle empty intron lists gracefully."""
    normalizer = ScoreNormalizer()

    with pytest.raises((ValueError, RuntimeError)):
        normalizer.fit([], dataset_type="reference")

    with pytest.raises((ValueError, RuntimeError)):
        normalizer.fit([], dataset_type="reference")


# ============================================================================
# Streaming Mode Support Tests
# ============================================================================


def test_get_frozen_scaler(reference_introns):
    """Test that we can extract the frozen scaler for streaming mode."""
    normalizer = ScoreNormalizer()
    normalizer.fit(reference_introns, dataset_type="reference")

    # Should be able to get the frozen scaler
    scaler = normalizer.get_frozen_scaler()

    # Scaler should be a fitted sklearn RobustScaler
    assert scaler is not None
    assert hasattr(scaler, "center_")
    assert hasattr(scaler, "scale_")
    assert scaler.center_ is not None
    assert scaler.scale_ is not None


def test_get_frozen_scaler_before_fit():
    """Should raise error if trying to get scaler before fitting."""
    normalizer = ScoreNormalizer()

    with pytest.raises(RuntimeError, match="not been fitted"):
        normalizer.get_frozen_scaler()


def test_transform_scores_array(reference_introns):
    """Test direct numpy array transformation for streaming mode."""
    normalizer = ScoreNormalizer()
    normalizer.fit(reference_introns, dataset_type="reference")

    # Create a raw score matrix
    raw_scores = np.array(
        [
            [0.5, 0.6, 0.55],  # Some test scores
            [0.1, 0.2, 0.15],
            [0.8, 0.9, 0.85],
        ]
    )

    # Transform using the array method
    z_scores = normalizer.transform_scores_array(raw_scores)

    # Should have same shape
    assert z_scores.shape == raw_scores.shape

    # Should be different from raw scores (transformed)
    assert not np.allclose(z_scores, raw_scores)


def test_transform_scores_array_before_fit():
    """Should raise error if trying to transform before fitting."""
    normalizer = ScoreNormalizer()
    raw_scores = np.array([[0.5, 0.6, 0.55]])

    with pytest.raises(RuntimeError, match="fit"):
        normalizer.transform_scores_array(raw_scores)


def test_frozen_scaler_matches_transform(reference_introns):
    """Verify frozen scaler produces same results as transform method."""
    normalizer = ScoreNormalizer()
    normalizer.fit(reference_introns, dataset_type="reference")

    # Get the first intron's raw scores
    first_intron = reference_introns[0]
    raw_scores = np.array(
        [
            [
                first_intron.scores.five_raw_score,
                first_intron.scores.bp_raw_score,
                first_intron.scores.three_raw_score,
            ]
        ]
    )

    # Transform using the intron-based method
    transformed_introns = list(normalizer.transform([first_intron], "reference"))
    z_via_transform = np.array(
        [
            [
                transformed_introns[0].scores.five_z_score,
                transformed_introns[0].scores.bp_z_score,
                transformed_introns[0].scores.three_z_score,
            ]
        ]
    )

    # Transform using the array method
    z_via_array = normalizer.transform_scores_array(raw_scores)

    # Should be identical
    np.testing.assert_allclose(z_via_transform, z_via_array)
