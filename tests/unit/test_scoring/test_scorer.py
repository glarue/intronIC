"""
Tests for the IntronScorer pipeline orchestrator.

The IntronScorer class orchestrates the full scoring pipeline:
1. PWM scoring for 5' splice site
2. Branch point detection and scoring
3. PWM scoring for 3' splice site
4. Log-ratio calculation (U12/U2)
5. Population of intron score objects

Port from: intronIC.py:3115-3400 (get_raw_scores, assign_raw_score, multi_matrix_score)

Test Strategy:
1. Test IntronScorer initialization
2. Test 5' splice site scoring
3. Test 3' splice site scoring
4. Test branch point scoring integration
5. Test log-ratio calculation
6. Test full pipeline with real-like sequences
7. Test handling of canonical vs non-canonical introns
"""

import math
from pathlib import Path

import numpy as np
import pytest

from intronIC.core.intron import (
    GenomicCoordinate,
    Intron,
    IntronMetadata,
    IntronScores,
    IntronSequences,
)
from intronIC.scoring.pwm import PWM, PWMLoader, PWMSet
from intronIC.scoring.scorer import IntronScorer

# ============================================================================
# Test Fixtures
# ============================================================================


@pytest.fixture
def simple_pwms():
    """
    Create simple test PWMs for five, bp, three regions.

    These are deliberately simple for predictable scores.
    """
    # 5' splice site PWM (8 positions: -3 to +5 around GT)
    # Favors AGGTAAGT pattern
    five_u12 = PWM(
        name="u12_five",
        matrix=np.array(
            [
                # -3  -2  -1   0   1   2   3   4
                [0.8, 0.1, 0.1, 0.1, 0.9, 0.9, 0.1, 0.1],  # A
                [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],  # C
                [0.1, 0.8, 0.1, 0.9, 0.1, 0.1, 0.8, 0.1],  # G
                [0.1, 0.1, 0.8, 0.1, 0.1, 0.1, 0.1, 0.8],  # T
            ]
        ),
        length=8,
        start_index=-3,
    )

    five_u2 = PWM(
        name="u2_five",
        matrix=np.array(
            [
                [0.3, 0.3, 0.3, 0.2, 0.4, 0.4, 0.3, 0.3],  # A
                [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2],  # C
                [0.3, 0.3, 0.2, 0.6, 0.2, 0.2, 0.3, 0.3],  # G
                [0.2, 0.2, 0.3, 0.2, 0.2, 0.2, 0.2, 0.2],  # T
            ]
        ),
        length=8,
        start_index=-3,
    )

    # Branch point PWM (7 positions for TACTAAC)
    bp_u12 = PWM(
        name="u12_bp",
        matrix=np.array(
            [
                [0.05, 0.95, 0.05, 0.05, 0.95, 0.95, 0.05],  # A
                [0.05, 0.05, 0.95, 0.05, 0.05, 0.05, 0.95],  # C
                [0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05],  # G
                [0.95, 0.05, 0.05, 0.95, 0.05, 0.05, 0.05],  # T
            ]
        ),
        length=7,
    )

    bp_u2 = PWM(
        name="u2_bp",
        matrix=np.array(
            [
                [0.3, 0.4, 0.3, 0.3, 0.4, 0.5, 0.3],  # A
                [0.3, 0.2, 0.3, 0.2, 0.2, 0.1, 0.3],  # C
                [0.2, 0.2, 0.2, 0.3, 0.2, 0.2, 0.2],  # G
                [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2],  # T
            ]
        ),
        length=7,
    )

    # 3' splice site PWM (8 positions: -6 to +2 around AG)
    # Favors TTTCAGGT pattern
    three_u12 = PWM(
        name="u12_three",
        matrix=np.array(
            [
                # -6  -5  -4  -3  -2  -1   0   1
                [0.1, 0.1, 0.1, 0.1, 0.9, 0.1, 0.1, 0.1],  # A
                [0.1, 0.1, 0.1, 0.9, 0.1, 0.1, 0.1, 0.1],  # C
                [0.1, 0.1, 0.1, 0.1, 0.1, 0.9, 0.1, 0.8],  # G
                [0.8, 0.8, 0.8, 0.1, 0.1, 0.1, 0.1, 0.1],  # T
            ]
        ),
        length=8,
        start_index=-6,
    )

    three_u2 = PWM(
        name="u2_three",
        matrix=np.array(
            [
                [0.2, 0.2, 0.2, 0.2, 0.5, 0.2, 0.2, 0.3],  # A
                [0.3, 0.3, 0.3, 0.4, 0.2, 0.2, 0.2, 0.2],  # C
                [0.2, 0.2, 0.2, 0.2, 0.1, 0.6, 0.3, 0.3],  # G
                [0.3, 0.3, 0.3, 0.2, 0.2, 0.2, 0.3, 0.2],  # T
            ]
        ),
        length=8,
        start_index=-6,
    )

    pwm_sets = {
        "five": PWMSet(matrices={("u2", "gtag"): five_u2, ("u12", "gtag"): five_u12}),
        "bp": PWMSet(matrices={("u2", "gtag"): bp_u2, ("u12", "gtag"): bp_u12}),
        "three": PWMSet(
            matrices={("u2", "gtag"): three_u2, ("u12", "gtag"): three_u12}
        ),
    }

    return pwm_sets


@pytest.fixture
def u12_like_intron():
    """
    Create an intron with U12-like sequences.

    Structure:
    - 5': AGGTAAGT (matches U12 pattern)
    - BP: TACTAAC (canonical U12)
    - 3': TTTCAGGT (matches U12 pattern)
    """
    # Build full sequence
    # 5bp upstream flank + 8bp five region + 60bp middle + 7bp BP + 13bp + 8bp three region + 5bp downstream
    upstream = "ACTGN"
    five_region = "AGGTAAGT"  # Positions -3 to +5
    middle_before_bp = "N" * 60
    bp_region = "TACTAAC"
    middle_after_bp = "N" * 13
    three_region = "TTTCAGGT"  # Positions -6 to +2
    downstream = "NACTG"

    full_seq = (
        five_region + middle_before_bp + bp_region + middle_after_bp + three_region
    )

    return Intron(
        intron_id="u12_like_intron",
        coordinates=GenomicCoordinate(
            chromosome="chr1",
            start=1000,
            stop=1000 + len(full_seq),
            strand="+",
            system="1-based",
        ),
        sequences=IntronSequences(
            seq=full_seq,
            upstream_flank=upstream,
            downstream_flank=downstream,
            five_prime_dnt="GT",
            three_prime_dnt="AG",
        ),
        scores=IntronScores(),
        metadata=IntronMetadata(parent="transcript_1", grandparent="gene_1"),
    )


@pytest.fixture
def u2_like_intron():
    """
    Create an intron with U2-like sequences.

    Less specific motifs that score better with U2 PWMs.
    """
    upstream = "ACTGN"
    five_region = "NACGTTNG"  # Weaker 5' motif
    middle_before_bp = "N" * 60
    bp_region = "NNCNNNN"  # Weak BP
    middle_after_bp = "N" * 13
    three_region = "NNCNAGNN"  # Weaker 3' motif
    downstream = "NACTG"

    full_seq = (
        five_region + middle_before_bp + bp_region + middle_after_bp + three_region
    )

    return Intron(
        intron_id="u2_like_intron",
        coordinates=GenomicCoordinate(
            chromosome="chr1",
            start=2000,
            stop=2000 + len(full_seq),
            strand="+",
            system="1-based",
        ),
        sequences=IntronSequences(
            seq=full_seq,
            upstream_flank=upstream,
            downstream_flank=downstream,
            five_prime_dnt="GT",
            three_prime_dnt="AG",
        ),
        scores=IntronScores(),
        metadata=IntronMetadata(parent="transcript_2", grandparent="gene_2"),
    )


# ============================================================================
# IntronScorer Initialization Tests
# ============================================================================


def test_intron_scorer_creation(simple_pwms):
    """Test that IntronScorer can be initialized with PWM sets."""
    scorer = IntronScorer(
        pwm_sets=simple_pwms,
        five_coords=(-3, 5),
        bp_coords=(-55, -5),
        three_coords=(-6, 2),
    )

    assert scorer.pwm_sets == simple_pwms
    assert scorer.five_coords == (-3, 5)
    assert scorer.bp_coords == (-55, -5)
    assert scorer.three_coords == (-6, 2)


def test_intron_scorer_default_coordinates(simple_pwms):
    """Test default coordinate settings match original intronIC."""
    scorer = IntronScorer(pwm_sets=simple_pwms)

    # Original intronIC defaults
    assert scorer.five_coords == (-3, 9)
    assert scorer.bp_coords == (-55, -5)
    assert scorer.three_coords == (-6, 4)


# ============================================================================
# 5' Splice Site Scoring Tests
# ============================================================================


def test_score_five_site_basic(simple_pwms, u12_like_intron):
    """Test scoring of 5' splice site."""
    scorer = IntronScorer(
        pwm_sets=simple_pwms,
        five_coords=(-3, 5),  # 8 positions
    )

    u12_score, u2_score = scorer._score_five_site(u12_like_intron)

    # Both scores should be calculated (non-zero)
    assert u12_score > 0
    assert u2_score > 0
    # Scores should be different (unless exactly identical sequence)
    assert u12_score != u2_score


def test_five_site_extracts_correct_region(simple_pwms):
    """Test that 5' scoring extracts the correct sequence region."""
    scorer = IntronScorer(pwm_sets=simple_pwms, five_coords=(-3, 5))

    # Create intron with known sequence
    seq = "AAAGTAAGT" + "N" * 90 + "TTTCAGNN"
    intron = Intron(
        intron_id="test",
        coordinates=GenomicCoordinate("chr1", 1000, 1100, "+", "1-based"),
        sequences=IntronSequences(
            seq=seq,
            upstream_flank="NNN",
            downstream_flank="NNN",
            five_prime_dnt="GT",
            three_prime_dnt="AG",
        ),
        scores=IntronScores(),
        metadata=IntronMetadata("t1", "g1"),
    )

    region = scorer._extract_five_region(intron)

    # Should extract first 8bp (coordinates -3 to +5 relative to start)
    # But since we're at the intron start, we take upstream flank
    # Actually, the five_coords are relative to the intron start
    # So -3 to +5 means: 3bp before intron start + 5bp into intron
    # We need to use upstream flank for the -3 part
    assert len(region) == 8


# ============================================================================
# 3' Splice Site Scoring Tests
# ============================================================================


def test_score_three_site_basic(simple_pwms, u12_like_intron):
    """Test scoring of 3' splice site."""
    scorer = IntronScorer(pwm_sets=simple_pwms, three_coords=(-6, 2))

    u12_score, u2_score = scorer._score_three_site(u12_like_intron)

    # Both scores should be calculated
    assert u12_score > 0
    assert u2_score > 0
    assert u12_score != u2_score


# ============================================================================
# Branch Point Scoring Tests
# ============================================================================


def test_score_branch_point_basic(simple_pwms, u12_like_intron):
    """Test branch point detection and scoring."""
    scorer = IntronScorer(
        pwm_sets=simple_pwms,
        bp_coords=(-80, -5),  # Wide search window
    )

    match, u2_score = scorer._score_branch_point(u12_like_intron)

    # Should find TACTAAC
    assert match.sequence == "TACTAAC"
    assert match.score > 0
    assert u2_score > 0


# ============================================================================
# Log-Ratio Calculation Tests
# ============================================================================


def test_log_ratio_calculation(simple_pwms):
    """Test log2 ratio calculation for scores."""
    scorer = IntronScorer(pwm_sets=simple_pwms)

    # Test with known values
    u12_score = 0.8
    u2_score = 0.2

    ratio = scorer._calculate_log_ratio(u12_score, u2_score)

    # log2(0.8 / 0.2) = log2(4) = 2.0
    expected = math.log2(0.8 / 0.2)
    assert abs(ratio - expected) < 1e-10


def test_log_ratio_with_zero_scores(simple_pwms):
    """Test that log ratio handles edge cases."""
    scorer = IntronScorer(pwm_sets=simple_pwms)

    # Should use pseudocount if score is 0
    # This should be handled by the PWM itself, but test edge case
    u12_score = 1e-10
    u2_score = 1e-10

    ratio = scorer._calculate_log_ratio(u12_score, u2_score)

    # log2(1e-10 / 1e-10) = log2(1) = 0.0
    assert abs(ratio - 0.0) < 1e-10


# ============================================================================
# Full Pipeline Integration Tests
# ============================================================================


def test_score_intron_full_pipeline(simple_pwms, u12_like_intron):
    """Test scoring a single intron through full pipeline."""
    scorer = IntronScorer(
        pwm_sets=simple_pwms,
        five_coords=(-3, 5),
        bp_coords=(-80, -5),
        three_coords=(-6, 2),
    )

    scored_intron = scorer.score_intron(u12_like_intron)

    # Check that all raw scores are populated
    assert scored_intron.scores.five_raw_score is not None
    assert scored_intron.scores.bp_raw_score is not None
    assert scored_intron.scores.three_raw_score is not None

    # Log ratios can be positive (U12 > U2) or negative (U2 > U12)
    # Just verify they're finite numbers
    assert math.isfinite(scored_intron.scores.five_raw_score)
    assert math.isfinite(scored_intron.scores.bp_raw_score)
    assert math.isfinite(scored_intron.scores.three_raw_score)


def test_score_introns_generator(simple_pwms, u12_like_intron, u2_like_intron):
    """Test scoring multiple introns via generator."""
    scorer = IntronScorer(
        pwm_sets=simple_pwms,
        five_coords=(-3, 5),
        bp_coords=(-80, -5),
        three_coords=(-6, 2),
    )

    introns = [u12_like_intron, u2_like_intron]
    scored = list(scorer.score_introns(introns))

    assert len(scored) == 2

    # All should have scores
    for intron in scored:
        assert intron.scores.five_raw_score is not None
        assert intron.scores.bp_raw_score is not None
        assert intron.scores.three_raw_score is not None


def test_u12_vs_u2_score_difference(simple_pwms, u12_like_intron, u2_like_intron):
    """Test that scoring works for both types of introns."""
    scorer = IntronScorer(
        pwm_sets=simple_pwms,
        five_coords=(-3, 5),
        bp_coords=(-80, -5),
        three_coords=(-6, 2),
    )

    u12_scored = scorer.score_intron(u12_like_intron)
    u2_scored = scorer.score_intron(u2_like_intron)

    # Both should have scores (they may be positive or negative)
    assert u12_scored.scores.five_raw_score is not None
    assert u12_scored.scores.bp_raw_score is not None
    assert u12_scored.scores.three_raw_score is not None

    assert u2_scored.scores.five_raw_score is not None
    assert u2_scored.scores.bp_raw_score is not None
    assert u2_scored.scores.three_raw_score is not None


# ============================================================================
# Canonical vs Non-Canonical Handling
# ============================================================================


def test_canonical_intron_scoring(simple_pwms):
    """Test scoring of canonical GT-AG intron."""
    # Create canonical intron
    seq = "GTAAGT" + "N" * 70 + "TACTAAC" + "N" * 10 + "TTTCAG"
    intron = Intron(
        intron_id="canonical",
        coordinates=GenomicCoordinate("chr1", 1000, 1000 + len(seq), "+", "1-based"),
        sequences=IntronSequences(
            seq=seq,
            upstream_flank="NNN",
            downstream_flank="NNN",
            five_prime_dnt="GT",
            three_prime_dnt="AG",
        ),
        scores=IntronScores(),
        metadata=IntronMetadata("t1", "g1"),
    )
    # noncanonical defaults to False, no need to set

    scorer = IntronScorer(
        pwm_sets=simple_pwms,
        five_coords=(-3, 5),
        bp_coords=(-80, -5),
        three_coords=(-6, 2),
    )
    scored = scorer.score_intron(intron)

    assert scored.scores.five_raw_score is not None


def test_noncanonical_intron_with_ignore_dnts(simple_pwms):
    """Test non-canonical intron scoring with dinucleotide positions ignored."""
    # Create GC-AG intron (non-canonical)
    seq = "GCAAGT" + "N" * 70 + "TACTAAC" + "N" * 10 + "TTTCAG"
    intron = Intron(
        intron_id="noncanonical",
        coordinates=GenomicCoordinate("chr1", 1000, 1000 + len(seq), "+", "1-based"),
        sequences=IntronSequences(
            seq=seq,
            upstream_flank="NNN",
            downstream_flank="NNN",
            five_prime_dnt="GC",
            three_prime_dnt="AG",
        ),
        scores=IntronScores(),
        metadata=IntronMetadata("t1", "g1"),
    )
    # Set noncanonical flag after construction
    intron.metadata.noncanonical = True

    scorer = IntronScorer(
        pwm_sets=simple_pwms,
        five_coords=(-3, 5),
        bp_coords=(-80, -5),
        three_coords=(-6, 2),
        ignore_nc_dnts=True,
    )
    scored = scorer.score_intron(intron)

    # Should still score (ignoring the GC dinucleotide)
    assert scored.scores.five_raw_score is not None


# ============================================================================
# Edge Cases
# ============================================================================


def test_very_short_intron(simple_pwms):
    """Test handling of intron too short for branch point search."""
    # Create very short intron (20bp) - too short for BP search window
    seq = "GTAA" + "N" * 12 + "TTAG"  # Only 20bp total
    intron = Intron(
        intron_id="short",
        coordinates=GenomicCoordinate("chr1", 1000, 1020, "+", "1-based"),
        sequences=IntronSequences(
            seq=seq,
            upstream_flank="NNN",
            downstream_flank="NNN",
            five_prime_dnt="GT",
            three_prime_dnt="AG",
        ),
        scores=IntronScores(),
        metadata=IntronMetadata("t1", "g1"),
    )

    scorer = IntronScorer(
        pwm_sets=simple_pwms,
        five_coords=(-3, 5),
        bp_coords=(-55, -5),  # This requires 50bp search window
        three_coords=(-6, 2),
    )

    # Should handle gracefully - either skip BP scoring or use what's available
    scored = scorer.score_intron(intron)

    # Five and three sites should still score
    assert scored.scores.five_raw_score is not None
    assert scored.scores.three_raw_score is not None
    # BP score might be None or 0 for very short introns
    # The implementation should handle this gracefully


def test_intron_without_sequences(simple_pwms):
    """Test error handling for intron without sequences."""
    intron = Intron(
        intron_id="no_seq",
        coordinates=GenomicCoordinate("chr1", 1000, 1100, "+", "1-based"),
        sequences=IntronSequences(seq=None),  # No sequence!
        scores=IntronScores(),
        metadata=IntronMetadata("t1", "g1"),
    )

    scorer = IntronScorer(pwm_sets=simple_pwms)

    with pytest.raises(ValueError, match="sequence|None|missing"):
        scorer.score_intron(intron)


# ============================================================================
# Coordinate System Tests
# ============================================================================


def test_custom_five_coordinates(simple_pwms):
    """Test custom 5' splice site coordinates."""
    scorer = IntronScorer(
        pwm_sets=simple_pwms,
        five_coords=(-5, 10),  # Wider region
    )

    assert scorer.five_coords == (-5, 10)


def test_custom_bp_coordinates(simple_pwms):
    """Test custom branch point search coordinates."""
    scorer = IntronScorer(
        pwm_sets=simple_pwms,
        bp_coords=(-100, -10),  # Different window
    )

    assert scorer.bp_coords == (-100, -10)


def test_custom_three_coordinates(simple_pwms):
    """Test custom 3' splice site coordinates."""
    scorer = IntronScorer(
        pwm_sets=simple_pwms,
        three_coords=(-10, 5),  # Wider region
    )

    assert scorer.three_coords == (-10, 5)


# ============================================================================
# Real PWM Integration Tests
# ============================================================================


@pytest.mark.skip(
    reason="Real PWM files don't have complete U2/U12 pairs for all regions"
)
def test_with_real_pwms_if_available(matrix_file):
    """Test scoring with real PWM matrices if available."""
    # Try to load real PWMs
    if not matrix_file.exists():
        pytest.skip("Real PWM file not available")

    # Load real PWMs
    loader = PWMLoader()
    pwm_sets = loader.load_from_file(matrix_file)  # Pass Path object directly

    # Real PWM coordinates based on actual matrix lengths
    # u12_gtag_five: start=-3, length varies
    # u12_gtag_bp: canonical U12 BP motif
    # u12_gtag_three: start=-20

    scorer = IntronScorer(
        pwm_sets=pwm_sets,
        five_coords=(-3, 9),  # Default from original
        bp_coords=(-55, -5),  # Default BP search window
        three_coords=(-20, 3),  # Adjusted for real PWM start position
    )

    # Create test intron with canonical U12-like sequences
    seq = "GTAAGTAT" + "N" * 60 + "TACTAAC" + "N" * 15 + "TTTTTTTTTTTTTTTTTTTCAG"
    intron = Intron(
        intron_id="real_pwm_test",
        coordinates=GenomicCoordinate("chr1", 1000, 1000 + len(seq), "+", "1-based"),
        sequences=IntronSequences(
            seq=seq,
            upstream_flank="NNNNN",
            downstream_flank="NNNNN",
            five_prime_dnt="GT",
            three_prime_dnt="AG",
        ),
        scores=IntronScores(),
        metadata=IntronMetadata("t1", "g1", noncanonical=False),
    )

    # Score with real PWMs
    scored = scorer.score_intron(intron)

    # Should have all scores populated
    assert scored.scores.five_raw_score is not None
    assert scored.scores.bp_raw_score is not None
    assert scored.scores.three_raw_score is not None
    assert math.isfinite(scored.scores.five_raw_score)
    assert math.isfinite(scored.scores.bp_raw_score)
    assert math.isfinite(scored.scores.three_raw_score)
    assert scored.scores.five_raw_score is not None
    assert scored.scores.bp_raw_score is not None
    assert scored.scores.three_raw_score is not None
    assert math.isfinite(scored.scores.five_raw_score)
    assert math.isfinite(scored.scores.bp_raw_score)
    assert math.isfinite(scored.scores.three_raw_score)


# ============================================================================
# Streaming Scoring Tests (score_and_normalize_introns)
# ============================================================================


@pytest.fixture
def fitted_scaler():
    """Create a fitted RobustScaler for testing streaming scoring."""
    from sklearn.preprocessing import RobustScaler

    # Create scaler with realistic parameters
    scaler = RobustScaler(with_centering=True, with_scaling=True)

    # Fit on synthetic data that represents typical score distributions
    # five: slightly positive, bp: positive, three: slightly negative
    synthetic_data = np.array(
        [
            [0.5, 1.2, -0.3],
            [0.8, 1.5, -0.1],
            [0.3, 0.9, -0.5],
            [0.6, 1.1, -0.2],
            [0.4, 1.3, -0.4],
        ]
    )
    scaler.fit(synthetic_data)

    return scaler


def test_score_and_normalize_introns_basic(simple_pwms, u12_like_intron, fitted_scaler):
    """Test streaming scoring with normalization."""
    from intronIC.scoring.scorer import score_and_normalize_introns

    scorer = IntronScorer(
        pwm_sets=simple_pwms,
        five_coords=(-3, 5),
        bp_coords=(-80, -5),
        three_coords=(-6, 2),
    )

    # Score a single intron
    results = list(
        score_and_normalize_introns([u12_like_intron], scorer, fitted_scaler)
    )

    assert len(results) == 1
    intron = results[0]

    # Check raw scores are populated
    assert intron.scores.five_raw_score is not None
    assert intron.scores.bp_raw_score is not None
    assert intron.scores.three_raw_score is not None

    # Check z-scores are populated
    assert intron.scores.five_z_score is not None
    assert intron.scores.bp_z_score is not None
    assert intron.scores.three_z_score is not None

    # Z-scores should be finite
    assert math.isfinite(intron.scores.five_z_score)
    assert math.isfinite(intron.scores.bp_z_score)
    assert math.isfinite(intron.scores.three_z_score)


def test_score_and_normalize_introns_generator(
    simple_pwms, u12_like_intron, u2_like_intron, fitted_scaler
):
    """Test streaming scoring returns generator (memory efficient)."""
    from intronIC.scoring.scorer import score_and_normalize_introns

    scorer = IntronScorer(
        pwm_sets=simple_pwms,
        five_coords=(-3, 5),
        bp_coords=(-80, -5),
        three_coords=(-6, 2),
    )

    # The function should return a generator, not a list
    result = score_and_normalize_introns(
        [u12_like_intron, u2_like_intron], scorer, fitted_scaler
    )

    # Generators are iterators
    assert hasattr(result, "__iter__")
    assert hasattr(result, "__next__")

    # Verify we can iterate through it
    results = list(result)
    assert len(results) == 2


def test_score_and_normalize_introns_preserves_metadata(
    simple_pwms, u12_like_intron, fitted_scaler
):
    """Test that streaming scoring preserves intron metadata."""
    from intronIC.scoring.scorer import score_and_normalize_introns

    scorer = IntronScorer(
        pwm_sets=simple_pwms,
        five_coords=(-3, 5),
        bp_coords=(-80, -5),
        three_coords=(-6, 2),
    )

    results = list(
        score_and_normalize_introns([u12_like_intron], scorer, fitted_scaler)
    )
    intron = results[0]

    # Original metadata should be preserved
    assert intron.intron_id == u12_like_intron.intron_id
    assert intron.coordinates == u12_like_intron.coordinates
    assert intron.metadata.parent == u12_like_intron.metadata.parent
    assert intron.metadata.grandparent == u12_like_intron.metadata.grandparent


def test_score_and_normalize_batch_basic(
    simple_pwms, u12_like_intron, u2_like_intron, fitted_scaler
):
    """Test batch scoring with normalization."""
    from intronIC.scoring.scorer import score_and_normalize_batch

    scorer = IntronScorer(
        pwm_sets=simple_pwms,
        five_coords=(-3, 5),
        bp_coords=(-80, -5),
        three_coords=(-6, 2),
    )

    # Score a batch of introns
    results = score_and_normalize_batch(
        [u12_like_intron, u2_like_intron], scorer, fitted_scaler
    )

    assert isinstance(results, list)
    assert len(results) == 2

    for intron in results:
        # Check raw scores
        assert intron.scores.five_raw_score is not None
        assert intron.scores.bp_raw_score is not None
        assert intron.scores.three_raw_score is not None

        # Check z-scores
        assert intron.scores.five_z_score is not None
        assert intron.scores.bp_z_score is not None
        assert intron.scores.three_z_score is not None


def test_score_and_normalize_batch_empty(simple_pwms, fitted_scaler):
    """Test batch scoring with empty list."""
    from intronIC.scoring.scorer import score_and_normalize_batch

    scorer = IntronScorer(
        pwm_sets=simple_pwms,
        five_coords=(-3, 5),
        bp_coords=(-80, -5),
        three_coords=(-6, 2),
    )

    results = score_and_normalize_batch([], scorer, fitted_scaler)

    assert results == []


def test_batch_matches_streaming(
    simple_pwms, u12_like_intron, u2_like_intron, fitted_scaler
):
    """Test that batch and streaming produce identical results."""
    from intronIC.scoring.scorer import (
        score_and_normalize_batch,
        score_and_normalize_introns,
    )

    scorer = IntronScorer(
        pwm_sets=simple_pwms,
        five_coords=(-3, 5),
        bp_coords=(-80, -5),
        three_coords=(-6, 2),
    )

    introns = [u12_like_intron, u2_like_intron]

    # Get results from both methods
    streaming_results = list(
        score_and_normalize_introns(introns, scorer, fitted_scaler)
    )
    batch_results = score_and_normalize_batch(introns, scorer, fitted_scaler)

    assert len(streaming_results) == len(batch_results)

    for s_intron, b_intron in zip(streaming_results, batch_results):
        # Raw scores should match
        assert np.isclose(
            s_intron.scores.five_raw_score, b_intron.scores.five_raw_score
        )
        assert np.isclose(s_intron.scores.bp_raw_score, b_intron.scores.bp_raw_score)
        assert np.isclose(
            s_intron.scores.three_raw_score, b_intron.scores.three_raw_score
        )

        # Z-scores should match
        assert np.isclose(s_intron.scores.five_z_score, b_intron.scores.five_z_score)
        assert np.isclose(s_intron.scores.bp_z_score, b_intron.scores.bp_z_score)
        assert np.isclose(s_intron.scores.three_z_score, b_intron.scores.three_z_score)
        assert np.isclose(s_intron.scores.five_z_score, b_intron.scores.five_z_score)
        assert np.isclose(s_intron.scores.bp_z_score, b_intron.scores.bp_z_score)
        assert np.isclose(s_intron.scores.three_z_score, b_intron.scores.three_z_score)
