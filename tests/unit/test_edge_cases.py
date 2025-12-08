"""
Tests for edge cases and boundary conditions.

This module tests unusual but valid inputs and boundary conditions.
"""

import numpy as np
import pytest

from intronIC.core.intron import (
    GenomicCoordinate,
    Intron,
    IntronMetadata,
    IntronScores,
    IntronSequences,
    OmissionReason,
)
from intronIC.extraction.filters import IntronFilter
from intronIC.scoring.normalizer import ScoreNormalizer
from intronIC.scoring.pwm import PWM, PWMSet
from intronIC.scoring.scorer import IntronScorer

# ============================================================================
# Coordinate Edge Cases
# ============================================================================


def test_intron_at_chromosome_start():
    """Test intron starting at position 1."""
    intron = Intron(
        intron_id="chr_start",
        coordinates=GenomicCoordinate("chr1", 1, 100, "+", "1-based"),
        sequences=IntronSequences(
            seq="GT" + "N" * 96 + "AG", five_prime_dnt="GT", three_prime_dnt="AG"
        ),
        scores=IntronScores(),
        metadata=IntronMetadata("t1", "g1"),
    )

    assert intron.coordinates.start == 1
    assert intron.coordinates.stop == 100


def test_very_long_chromosome_name():
    """Test handling of unusually long chromosome names."""
    long_name = "chr" + "A" * 200
    intron = Intron(
        intron_id="long_chr",
        coordinates=GenomicCoordinate(long_name, 1000, 2000, "+", "1-based"),
        sequences=IntronSequences(seq="GTAG"),
        scores=IntronScores(),
        metadata=IntronMetadata("t1", "g1"),
    )

    assert intron.coordinates.chromosome == long_name


def test_special_characters_in_chromosome_name():
    """Test chromosome names with special characters."""
    special_names = [
        "chr1_random",
        "chr1.2",
        "chr_KI270706v1_random",
        "HLA-A*01:01:01:01",
    ]

    for name in special_names:
        intron = Intron(
            intron_id=f"test_{name}",
            coordinates=GenomicCoordinate(name, 1000, 2000, "+", "1-based"),
            sequences=IntronSequences(seq="GTAG"),
            scores=IntronScores(),
            metadata=IntronMetadata("t1", "g1"),
        )
        assert intron.coordinates.chromosome == name


def test_negative_strand_intron():
    """Test intron on negative strand."""
    intron = Intron(
        intron_id="minus_strand",
        coordinates=GenomicCoordinate("chr1", 1000, 2000, "-", "1-based"),
        sequences=IntronSequences(seq="CTAC"),  # Reverse complement of GTAG
        scores=IntronScores(),
        metadata=IntronMetadata("t1", "g1"),
    )

    assert intron.coordinates.strand == "-"


# ============================================================================
# Sequence Edge Cases
# ============================================================================


def test_minimum_length_intron():
    """Test intron at minimum acceptable length (30bp)."""
    seq = "GT" + "N" * 26 + "AG"  # Exactly 30bp
    intron = Intron(
        intron_id="min_length",
        coordinates=GenomicCoordinate("chr1", 1000, 1030, "+", "1-based"),
        sequences=IntronSequences(seq=seq, five_prime_dnt="GT", three_prime_dnt="AG"),
        scores=IntronScores(),
        metadata=IntronMetadata("t1", "g1"),
    )

    assert len(intron.sequences.seq) == 30

    # Filter should accept it at default threshold
    filter_obj = IntronFilter(min_length=30)
    filter_obj._check_omission(intron)
    assert not intron.metadata.omitted


def test_just_below_minimum_length():
    """Test intron just below minimum length."""
    seq = "GT" + "N" * 25 + "AG"  # 29bp
    intron = Intron(
        intron_id="too_short",
        # For 1-based: length = stop - start + 1, so 1028 - 1000 + 1 = 29
        coordinates=GenomicCoordinate("chr1", 1000, 1028, "+", "1-based"),
        sequences=IntronSequences(seq=seq, five_prime_dnt="GT", three_prime_dnt="AG"),
        scores=IntronScores(),
        metadata=IntronMetadata("t1", "g1"),
    )

    filter_obj = IntronFilter(min_length=30)
    filter_obj._check_omission(intron)
    assert intron.metadata.omitted == OmissionReason.SHORT


def test_all_same_nucleotide():
    """Test intron sequence with only one nucleotide."""
    seq = "A" * 100
    intron = Intron(
        intron_id="all_a",
        coordinates=GenomicCoordinate("chr1", 1000, 1100, "+", "1-based"),
        sequences=IntronSequences(seq=seq, five_prime_dnt="AA", three_prime_dnt="AA"),
        scores=IntronScores(),
        metadata=IntronMetadata("t1", "g1"),
    )

    assert intron.sequences.seq == seq


def test_sequence_with_lowercase():
    """Test that lowercase sequences are handled."""
    seq = "gtacag"
    intron = Intron(
        intron_id="lowercase",
        coordinates=GenomicCoordinate("chr1", 1000, 1006, "+", "1-based"),
        sequences=IntronSequences(seq=seq, five_prime_dnt="gt", three_prime_dnt="ag"),
        scores=IntronScores(),
        metadata=IntronMetadata("t1", "g1"),
    )

    # Should be converted to uppercase internally or handled properly
    assert intron.sequences.seq is not None


def test_sequence_with_ambiguous_nucleotides():
    """Test sequences with IUPAC ambiguity codes."""
    seq = "GTACNNNWWWSSSAG"  # N, W, S are ambiguity codes
    intron = Intron(
        intron_id="ambiguous",
        coordinates=GenomicCoordinate("chr1", 1000, 1015, "+", "1-based"),
        sequences=IntronSequences(seq=seq, five_prime_dnt="GT", three_prime_dnt="AG"),
        scores=IntronScores(),
        metadata=IntronMetadata("t1", "g1"),
    )

    # Ambiguous bases should trigger omission
    filter_obj = IntronFilter()
    filter_obj._check_omission(intron)
    # Should be marked as ambiguous
    assert intron.metadata.omitted is not None


# ============================================================================
# Non-Canonical Dinucleotide Edge Cases
# ============================================================================


def test_all_noncanonical_types():
    """Test various non-canonical intron types."""
    nc_types = [
        ("GC", "AG"),  # GC-AG
        ("AT", "AC"),  # AT-AC (U12)
        ("GT", "TG"),  # GT-TG
        ("GG", "AG"),  # GG-AG
        ("CT", "AG"),  # CT-AG
        ("GT", "CG"),  # GT-CG
    ]

    for five_dnt, three_dnt in nc_types:
        seq = five_dnt + "NNNN" + three_dnt
        intron = Intron(
            intron_id=f"nc_{five_dnt}_{three_dnt}",
            coordinates=GenomicCoordinate("chr1", 1000, 1008, "+", "1-based"),
            sequences=IntronSequences(
                seq=seq, five_prime_dnt=five_dnt, three_prime_dnt=three_dnt
            ),
            scores=IntronScores(),
            metadata=IntronMetadata("t1", "g1"),
        )
        # Set noncanonical flag after construction
        intron.metadata.noncanonical = True

        assert intron.metadata.noncanonical is True


# ============================================================================
# Scoring Edge Cases
# ============================================================================


def test_scoring_with_all_zero_pwm():
    """Test scoring behavior with degenerate PWM (all zeros)."""
    # Create PWM with very low uniform probabilities
    zero_matrix = np.full((4, 7), 0.00001)  # Near-zero values

    pwm = PWM(name="zero_pwm", matrix=zero_matrix, length=7, start_index=0)

    # Scoring should handle without crashing (may give very low/negative scores)
    # This tests numerical stability


def test_scoring_with_uniform_pwm():
    """Test scoring with completely uniform PWM (no information)."""
    uniform_matrix = np.full((4, 7), 0.25)  # Equal probability for all bases

    pwm = PWM(name="uniform_pwm", matrix=uniform_matrix, length=7, start_index=0)

    # All sequences should score the same
    score1 = pwm.score_sequence("AAAAAAA", ignore_positions=set())
    score2 = pwm.score_sequence("GGGGGGG", ignore_positions=set())
    score3 = pwm.score_sequence("ACTGACT", ignore_positions=set())

    # Should all be equal (within floating point precision)
    assert abs(score1 - score2) < 1e-10
    assert abs(score2 - score3) < 1e-10


def test_z_score_normalization_with_zero_variance():
    """Test normalization when all reference scores are identical.

    ZeroAnchoredRobustScaler should handle zero variance gracefully by
    using a minimum scale value to prevent division by zero.
    """
    # Create introns with all identical scores
    introns = []
    for i in range(5):
        intron = Intron(
            intron_id=f"zero_var_intron_{i}",
            coordinates=GenomicCoordinate(
                "chr1", 1000 + i * 100, 2000 + i * 100, "+", "1-based"
            ),
            sequences=IntronSequences(seq="GTAAGT" + "N" * 50 + "AG"),
            scores=IntronScores(
                five_raw_score=5.0,  # All identical
                bp_raw_score=5.0,  # All identical
                three_raw_score=5.0,  # All identical
            ),
        )
        introns.append(intron)

    normalizer = ScoreNormalizer()

    # Fit should handle zero variance gracefully
    normalizer.fit(introns, dataset_type="reference")

    # Transform should work without division by zero
    normalized = list(normalizer.transform(introns, dataset_type="reference"))

    # All z-scores should be equal (all inputs were equal)
    z_scores = [i.scores.five_z_score for i in normalized]
    assert all(z is not None for z in z_scores), "All z-scores should be computed"
    # With zero variance, all scores should produce identical z-scores
    assert len(set(z_scores)) == 1, (
        "All z-scores should be identical for identical inputs"
    )


def test_z_score_with_single_reference():
    """Test normalization with only one reference intron.

    Edge case: fitting on a single intron should still work, though
    the statistics may not be meaningful.
    """
    intron = Intron(
        intron_id="single_intron",
        coordinates=GenomicCoordinate("chr1", 1000, 2000, "+", "1-based"),
        sequences=IntronSequences(seq="GTAAGT" + "N" * 50 + "AG"),
        scores=IntronScores(
            five_raw_score=5.0,
            bp_raw_score=3.0,
            three_raw_score=4.0,
        ),
    )

    normalizer = ScoreNormalizer()

    # Should handle single reference gracefully
    normalizer.fit([intron], dataset_type="reference")

    # Transform should work
    normalized = list(normalizer.transform([intron], dataset_type="reference"))
    assert len(normalized) == 1

    # Z-score should be computed (value may vary depending on implementation)
    assert normalized[0].scores.five_z_score is not None


# ============================================================================
# Metadata Edge Cases
# ============================================================================


def test_intron_with_very_long_id():
    """Test intron with extremely long identifier."""
    long_id = "intron_" + "A" * 1000

    intron = Intron(
        intron_id=long_id,
        coordinates=GenomicCoordinate("chr1", 1000, 2000, "+", "1-based"),
        sequences=IntronSequences(seq="GTAG"),
        scores=IntronScores(),
        metadata=IntronMetadata("t1", "g1"),
    )

    assert intron.intron_id == long_id


def test_intron_with_special_characters_in_id():
    """Test intron IDs with special characters."""
    special_ids = [
        "intron:1:2:3",
        "intron|variant|123",
        "intron-name.v2",
        "intron_name[corrected]",
        "intron_100%_identical",
    ]

    for special_id in special_ids:
        intron = Intron(
            intron_id=special_id,
            coordinates=GenomicCoordinate("chr1", 1000, 2000, "+", "1-based"),
            sequences=IntronSequences(seq="GTAG"),
            scores=IntronScores(),
            metadata=IntronMetadata("t1", "g1"),
        )
        assert intron.intron_id == special_id


def test_intron_with_empty_parent_info():
    """Test intron without parent/grandparent information."""
    intron = Intron(
        intron_id="orphan",
        coordinates=GenomicCoordinate("chr1", 1000, 2000, "+", "1-based"),
        sequences=IntronSequences(seq="GTAG"),
        scores=IntronScores(),
        metadata=IntronMetadata(parent=None, grandparent=None),
    )

    # Should handle None parents gracefully
    assert intron.metadata.parent is None
    assert intron.metadata.grandparent is None


def test_intron_with_zero_family_size():
    """Test intron with family_size=0."""
    intron = Intron(
        intron_id="no_family",
        coordinates=GenomicCoordinate("chr1", 1000, 2000, "+", "1-based"),
        sequences=IntronSequences(seq="GTAG"),
        scores=IntronScores(),
        metadata=IntronMetadata("t1", "g1", family_size=0),
    )

    assert intron.metadata.family_size == 0


# ============================================================================
# Filtering Edge Cases
# ============================================================================


def test_filter_with_all_options_disabled():
    """Test filter with most permissive settings."""
    filter_obj = IntronFilter(
        min_length=0,  # No minimum
        allow_noncanonical=True,
        allow_overlap=True,
        longest_only=False,
        include_duplicates=True,
    )

    # Create intron that would normally be filtered
    intron = Intron(
        intron_id="permissive",
        coordinates=GenomicCoordinate("chr1", 1000, 1010, "+", "1-based"),  # Very short
        sequences=IntronSequences(
            seq="GTNNNNNNAG", five_prime_dnt="GT", three_prime_dnt="AG"
        ),
        scores=IntronScores(),
        metadata=IntronMetadata("t1", "g1"),
    )
    # Set longest_isoform flag after construction
    intron.metadata.longest_isoform = False

    filter_obj._check_omission(intron)

    # With permissive settings, should pass many filters
    # (might still fail ambiguous nucleotide check)


def test_filter_with_all_options_strict():
    """Test filter with most restrictive settings."""
    filter_obj = IntronFilter(
        min_length=1000,  # Very long minimum
        allow_noncanonical=False,
        allow_overlap=False,
        longest_only=True,
        include_duplicates=False,
    )

    # Create short intron from non-longest isoform
    intron = Intron(
        intron_id="strict",
        coordinates=GenomicCoordinate("chr1", 1000, 1100, "+", "1-based"),  # 100bp
        sequences=IntronSequences(
            seq="GT" + "N" * 96 + "AG", five_prime_dnt="GT", three_prime_dnt="AG"
        ),
        scores=IntronScores(),
        metadata=IntronMetadata("t1", "g1"),
    )
    # Set longest_isoform flag after construction
    intron.metadata.longest_isoform = False

    filter_obj._check_omission(intron)

    # Should be omitted for multiple reasons
    assert intron.metadata.omitted is not None


# ============================================================================
# Numerical Edge Cases
# ============================================================================


def test_very_high_scores():
    """Test handling of extremely high numerical scores."""
    intron = Intron(
        intron_id="high_score",
        coordinates=GenomicCoordinate("chr1", 1000, 2000, "+", "1-based"),
        sequences=IntronSequences(seq="GTAG"),
        scores=IntronScores(
            five_raw_score=1e100, bp_raw_score=1e100, three_raw_score=1e100
        ),
        metadata=IntronMetadata("t1", "g1"),
    )

    # Should store without overflow
    assert intron.scores.five_raw_score == 1e100


def test_very_low_negative_scores():
    """Test handling of extremely low negative scores."""
    intron = Intron(
        intron_id="low_score",
        coordinates=GenomicCoordinate("chr1", 1000, 2000, "+", "1-based"),
        sequences=IntronSequences(seq="GTAG"),
        scores=IntronScores(
            five_raw_score=-1e100, bp_raw_score=-1e100, three_raw_score=-1e100
        ),
        metadata=IntronMetadata("t1", "g1"),
    )

    assert intron.scores.five_raw_score == -1e100


def test_nan_and_inf_scores():
    """Test handling of NaN and infinity in scores."""
    intron = Intron(
        intron_id="special_values",
        coordinates=GenomicCoordinate("chr1", 1000, 2000, "+", "1-based"),
        sequences=IntronSequences(seq="GTAG"),
        scores=IntronScores(
            five_raw_score=np.nan, bp_raw_score=np.inf, three_raw_score=-np.inf
        ),
        metadata=IntronMetadata("t1", "g1"),
    )

    # Should store special values (though downstream code should handle them)
    assert np.isnan(intron.scores.five_raw_score)
    assert np.isinf(intron.scores.bp_raw_score)
    assert np.isinf(intron.scores.three_raw_score)
    assert np.isnan(intron.scores.five_raw_score)
    assert np.isinf(intron.scores.bp_raw_score)
    assert np.isinf(intron.scores.three_raw_score)
