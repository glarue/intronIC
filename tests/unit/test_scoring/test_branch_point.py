"""
Tests for branch point detection.

Branch point detection uses a sliding window approach to find the optimal
branch point sequence within a search region upstream of the 3' splice site.

Port from: intronIC.py:2143-2178 (bp_score)

Test Strategy:
1. Test BranchPointMatch data structure
2. Test sliding window search algorithm
3. Test integration with Intron objects
4. Test edge cases (short sequences, no clear winner, etc.)
"""

from dataclasses import replace
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
from intronIC.scoring.branch_point import BranchPointMatch, BranchPointScorer
from intronIC.scoring.pwm import PWM

# ============================================================================
# Test Fixtures
# ============================================================================


@pytest.fixture
def u12_bp_pwm() -> PWM:
    """
    U12 branch point PWM (TACTAAC motif).

    This is a simplified version focusing on the canonical TACTAAC sequence.
    """
    # 7-position PWM for TACTAAC
    # High frequency for T-A-C-T-A-A-C pattern
    matrix = np.array(
        [
            # Position: 0    1    2    3    4    5    6
            [0.05, 0.95, 0.05, 0.05, 0.95, 0.95, 0.05],  # A
            [0.05, 0.05, 0.95, 0.05, 0.05, 0.05, 0.95],  # C
            [0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05],  # G
            [0.95, 0.05, 0.05, 0.95, 0.05, 0.05, 0.05],  # T
        ]
    )

    return PWM(name="u12_bp_test", matrix=matrix, length=7, pseudocount=0.0001)


@pytest.fixture
def u2_bp_pwm() -> PWM:
    """
    U2 branch point PWM (degenerate YNYURAY motif).

    This is more degenerate than U12, representing the variable U2 BP.
    """
    # 7-position PWM with more variability
    matrix = np.array(
        [
            # Position: 0    1    2    3    4    5    6
            [0.30, 0.20, 0.30, 0.40, 0.30, 0.50, 0.30],  # A
            [0.30, 0.30, 0.20, 0.20, 0.20, 0.10, 0.30],  # C
            [0.20, 0.20, 0.30, 0.20, 0.30, 0.20, 0.20],  # G
            [0.20, 0.30, 0.20, 0.20, 0.20, 0.20, 0.20],  # T
        ]
    )

    return PWM(name="u2_bp_test", matrix=matrix, length=7, pseudocount=0.0001)


@pytest.fixture
def simple_intron() -> Intron:
    """
    Create a simple intron with known sequence for testing.

    Structure:
    - 100bp total
    - Contains TACTAAC at position 70-76 (30bp before end)
    """
    # Build sequence: 70bp filler + TACTAAC + 23bp filler
    seq = "N" * 70 + "TACTAAC" + "N" * 23

    return Intron(
        intron_id="test_intron_1",
        coordinates=GenomicCoordinate(
            chromosome="chr1", start=1000, stop=1100, strand="+", system="1-based"
        ),
        sequences=IntronSequences(
            seq=seq,
            upstream_flank="ACTG",
            downstream_flank="TGCA",
            five_prime_dnt="GT",
            three_prime_dnt="AG",
        ),
        scores=IntronScores(),
        metadata=IntronMetadata(parent="transcript_1", grandparent="gene_1"),
    )


@pytest.fixture
def complex_intron() -> Intron:
    """
    Create intron with multiple possible branch points.

    Contains:
    - Weak TACGAAC at position 60-66
    - Strong TACTAAC at position 80-86
    Total length: 100bp
    """
    # Build sequence with two potential BPs
    seq = "N" * 60 + "TACGAAC" + "N" * 13 + "TACTAAC" + "N" * 13

    return Intron(
        intron_id="test_intron_2",
        coordinates=GenomicCoordinate(
            chromosome="chr1", start=2000, stop=2100, strand="+", system="1-based"
        ),
        sequences=IntronSequences(
            seq=seq,
            upstream_flank="ACTG",
            downstream_flank="TGCA",
            five_prime_dnt="GT",
            three_prime_dnt="AG",
        ),
        scores=IntronScores(),
        metadata=IntronMetadata(parent="transcript_2", grandparent="gene_2"),
    )


# ============================================================================
# BranchPointMatch Data Structure Tests
# ============================================================================


def test_branch_point_match_creation():
    """Test that BranchPointMatch can be created with required fields."""
    match = BranchPointMatch(
        sequence="TACTAAC",
        score=0.85,
        position=-30,  # 30bp before 3' end
        start_in_region=10,
        stop_in_region=17,
    )

    assert match.sequence == "TACTAAC"
    assert match.score == 0.85
    assert match.position == -30
    assert match.start_in_region == 10
    assert match.stop_in_region == 17


def test_branch_point_match_immutable():
    """BranchPointMatch should be immutable (frozen)."""
    match = BranchPointMatch(
        sequence="TACTAAC",
        score=0.85,
        position=-30,
        start_in_region=10,
        stop_in_region=17,
    )

    with pytest.raises((AttributeError, Exception)):
        match.score = 0.95


# ============================================================================
# BranchPointScorer Basic Tests
# ============================================================================


def test_branch_point_scorer_creation(u12_bp_pwm, u2_bp_pwm):
    """Test that BranchPointScorer can be initialized."""
    scorer = BranchPointScorer(u12_bp_pwm, u2_bp_pwm)

    assert scorer.u12_pwm == u12_bp_pwm
    assert scorer.u2_pwm == u2_bp_pwm


def test_sliding_window_search_simple_sequence(u12_bp_pwm, u2_bp_pwm):
    """Test sliding window search on a simple sequence."""
    scorer = BranchPointScorer(u12_bp_pwm, u2_bp_pwm)

    # Sequence: 5N + TACTAAC + 3N = 15bp total
    sequence = "NNNNN" + "TACTAAC" + "NNN"

    match = scorer._find_best_in_sequence(sequence, u12_bp_pwm, search_window_start=-55)

    # Should find TACTAAC
    assert match.sequence == "TACTAAC"
    assert match.start_in_region == 5
    assert match.stop_in_region == 12
    # Score should be high (all bases match high-frequency positions)
    assert match.score > 0.5


def test_sliding_window_finds_best_match(u12_bp_pwm, u2_bp_pwm):
    """Test that sliding window finds the highest-scoring match."""
    scorer = BranchPointScorer(u12_bp_pwm, u2_bp_pwm)

    # Sequence with two variants: weak + strong
    # TACGAAC (weak: G instead of T at position 3)
    # TACTAAC (strong: perfect match)
    sequence = "TACGAAC" + "NNN" + "TACTAAC"

    match = scorer._find_best_in_sequence(sequence, u12_bp_pwm, search_window_start=-55)

    # Should find the second (stronger) TACTAAC
    assert match.sequence == "TACTAAC"
    assert match.start_in_region == 10  # After weak variant + 3N


def test_sequence_too_short_for_pwm(u12_bp_pwm, u2_bp_pwm):
    """Should handle sequences shorter than PWM length."""
    scorer = BranchPointScorer(u12_bp_pwm, u2_bp_pwm)

    # Sequence shorter than 7bp PWM
    short_seq = "NNNNN"

    # Should either raise error or return None/sentinel
    with pytest.raises(ValueError, match="too short|shorter than"):
        scorer._find_best_in_sequence(short_seq, u12_bp_pwm, search_window_start=-55)


def test_sequence_exactly_pwm_length(u12_bp_pwm, u2_bp_pwm):
    """Test sequence exactly equal to PWM length."""
    scorer = BranchPointScorer(u12_bp_pwm, u2_bp_pwm)

    # Exactly 7bp
    sequence = "TACTAAC"

    match = scorer._find_best_in_sequence(sequence, u12_bp_pwm, search_window_start=-55)

    assert match.sequence == "TACTAAC"
    assert match.start_in_region == 0
    assert match.stop_in_region == 7


# ============================================================================
# Integration Tests with Intron Objects
# ============================================================================


def test_find_best_match_in_intron(u12_bp_pwm, u2_bp_pwm, simple_intron):
    """Test finding branch point in a full intron."""
    scorer = BranchPointScorer(u12_bp_pwm, u2_bp_pwm)

    # Search window: -50 to -20 relative to 3' end
    # Intron is 100bp, so this is positions 50-80
    match = scorer.find_best_match(simple_intron, search_window=(-50, -20))

    # Should find TACTAAC at position 70-76
    assert match.sequence == "TACTAAC"
    # Position should be relative to 3' end
    # TACTAAC starts at index 70, intron is 100bp, so -30bp from end
    assert match.position == -30


def test_search_window_extraction(u12_bp_pwm, u2_bp_pwm, simple_intron):
    """Test that search window is correctly extracted from intron."""
    scorer = BranchPointScorer(u12_bp_pwm, u2_bp_pwm)

    # Extract region (with default five_coords and three_coords)
    region, start_pos = scorer._extract_search_region(
        simple_intron,
        search_window=(-80, -60),
        five_coords=(-3, 9),
        three_coords=(-6, 4),
    )

    # Intron is 100bp, so:
    # -80 to -60 from end = positions 20-40
    # five_end = 9, three_start = 100 + (-6) = 94
    # start_pos = max(20, 9) = 20
    # stop_pos = min(40, 94) = 40
    assert len(region) == 20
    assert start_pos == 20
    assert region == simple_intron.sequences.seq[20:40]


def test_default_search_window(u12_bp_pwm, u2_bp_pwm, simple_intron):
    """Test using default search window (-55, -5)."""
    scorer = BranchPointScorer(u12_bp_pwm, u2_bp_pwm)

    # Should use default window if not specified
    match = scorer.find_best_match(simple_intron)

    # Should still find TACTAAC (it's at -30, within default -55 to -5)
    assert match.sequence == "TACTAAC"


def test_multiple_candidates_chooses_best(u12_bp_pwm, u2_bp_pwm, complex_intron):
    """Test that scorer chooses highest-scoring match when multiple present."""
    scorer = BranchPointScorer(u12_bp_pwm, u2_bp_pwm)

    # Complex intron has TACGAAC (weak) and TACTAAC (strong)
    match = scorer.find_best_match(complex_intron, search_window=(-50, -5))

    # Should find the stronger TACTAAC, not the weaker TACGAAC
    assert match.sequence == "TACTAAC"

    # Verify it found the second occurrence (at position 80-86)
    # Intron is 100bp, TACTAAC at 80-86
    # Search window (-50, -5) starts at position 50
    # TACTAAC is at offset 30 in the search region
    # Position: -50 + 30 = -20
    assert match.position == -20


# ============================================================================
# Edge Cases
# ============================================================================


def test_no_clear_winner_all_equal_scores(u12_bp_pwm, u2_bp_pwm):
    """Test behavior when all positions score equally (e.g., all Ns)."""
    scorer = BranchPointScorer(u12_bp_pwm, u2_bp_pwm)

    # All Ns - should score equally everywhere (all get pseudocount)
    sequence = "N" * 50

    match = scorer._find_best_in_sequence(sequence, u12_bp_pwm, search_window_start=-55)

    # Should return first occurrence when tied
    assert match.start_in_region == 0
    # Score should be low (all pseudocounts)
    assert match.score < 0.01


def test_search_window_at_intron_boundaries(u12_bp_pwm, u2_bp_pwm, simple_intron):
    """Test search window at the edges of the intron."""
    scorer = BranchPointScorer(u12_bp_pwm, u2_bp_pwm)

    # Search very close to 3' end
    match = scorer.find_best_match(simple_intron, search_window=(-20, -5))

    # Should still work, just searching a smaller region
    assert match is not None
    assert match.sequence is not None


def test_search_window_too_small_for_pwm(u12_bp_pwm, u2_bp_pwm, simple_intron):
    """Test search window smaller than PWM length returns None.

    Port from: intronIC.py:2944 - returns None for short sequences,
    which is then handled with pseudocount by the scorer.
    """
    scorer = BranchPointScorer(u12_bp_pwm, u2_bp_pwm)

    # Window of only 5bp, but PWM is 7bp - should return None
    match = scorer.find_best_match(simple_intron, search_window=(-10, -5))
    assert match is None, "Should return None for windows too small for PWM"


def test_short_intron_excludes_five_prime_region(u12_bp_pwm, u2_bp_pwm):
    """Test that short introns properly exclude 5' scoring region from BP search.

    This tests the fix for the issue where very short introns had BPS overlapping
    with the 5' splice site scoring region.

    Note: The original intronIC v1.5.1 ONLY excludes the 5' region. BP is allowed
    to overlap with the 3' region (typically 1 position overlap in short introns).
    """
    scorer = BranchPointScorer(u12_bp_pwm, u2_bp_pwm)

    # Create a 47bp intron (like the examples in the bug report)
    # Structure: 9bp (5' region) + 33bp (searchable) + 6bp (3' region) = 48bp total
    # Note: 1bp overlap at position 41 (searchable goes to 41, 3' starts at 41)
    # Put TACTAAC in the middle of the searchable region (position 20-26)
    seq = "GTAAGTNNN" + "NNNNNNNNN" + "NN" + "TACTAAC" + "N" * 14 + "TTTCAG"
    #     ^^^^^^^^^ 5' region (0-8, 9bp)
    #     0-8         9-17 (9bp)    18-19  20-26 (7bp)  27-40 (14bp)  41-46 (6bp 3' region)
    #                               TACTAAC at absolute position 20-26
    #     Total length: 9+9+2+7+14+6 = 47bp

    intron = Intron(
        intron_id="short_intron_47bp",
        coordinates=GenomicCoordinate(
            chromosome="chr1", start=1000, stop=1047, strand="+", system="1-based"
        ),
        sequences=IntronSequences(
            seq=seq,
            upstream_flank="CAG",
            downstream_flank="CCT",
            five_prime_dnt="GT",
            three_prime_dnt="AG",
        ),
        scores=IntronScores(),
        metadata=IntronMetadata(parent="transcript_short", grandparent="gene_short"),
    )

    # Search with standard window (-55, -5)
    # For 47bp intron: start would be 47+(-55)=-8, stop would be 47+(-5)=42
    # After clamping: start=max(-8, 9)=9, stop=min(42, 47)=42
    match = scorer.find_best_match(
        intron, search_window=(-55, -5), five_coords=(-3, 9), three_coords=(-6, 4)
    )

    # Should find TACTAAC
    assert match is not None, "Should find a match in the valid search region"
    assert match.sequence == "TACTAAC"

    # Verify that the search region excludes the 5' region
    # Position should be relative to 3' end
    # TACTAAC at position 20, intron length 47
    # Position: 20 - 47 = -27
    assert match.position == -27

    # Extract the actual search region to verify it excludes 5' region
    search_region, start_pos = scorer._extract_search_region(
        intron, search_window=(-55, -5), five_coords=(-3, 9), three_coords=(-6, 4)
    )

    # Search region should be positions [9, 42) = 33bp
    assert len(search_region) == 33, f"Expected 33bp, got {len(search_region)}"
    assert start_pos == 9
    assert search_region == seq[9:42]

    # Verify no overlap with 5' region (first 9bp)
    assert search_region != seq[0:33]  # Would include 5' region at start

    # Note: 1bp overlap with 3' region at position 41 is EXPECTED and ALLOWED
    # (matches original intronIC v1.5.1 behavior)


def test_negative_strand_intron(u12_bp_pwm, u2_bp_pwm):
    """Test branch point detection on negative-strand intron."""
    # Create negative-strand intron
    seq = "N" * 70 + "TACTAAC" + "N" * 23

    intron = Intron(
        intron_id="test_intron_minus",
        coordinates=GenomicCoordinate(
            chromosome="chr1", start=1000, stop=1100, strand="-", system="1-based"
        ),
        sequences=IntronSequences(
            seq=seq,  # Already reverse-complemented by extraction
            upstream_flank="ACTG",
            downstream_flank="TGCA",
            five_prime_dnt="GT",
            three_prime_dnt="AG",
        ),
        scores=IntronScores(),
        metadata=IntronMetadata(parent="transcript_3", grandparent="gene_3"),
    )

    scorer = BranchPointScorer(u12_bp_pwm, u2_bp_pwm)
    match = scorer.find_best_match(intron, search_window=(-50, -20))

    # Should still find TACTAAC (sequences already oriented correctly)
    assert match.sequence == "TACTAAC"


def test_ambiguous_bases_in_match(u12_bp_pwm, u2_bp_pwm):
    """Test handling of ambiguous bases (N, Y, R, etc.) in best match."""
    scorer = BranchPointScorer(u12_bp_pwm, u2_bp_pwm)

    # Sequence with some Ns in the match
    sequence = "NNNNN" + "TACTNNC" + "NNN"

    match = scorer._find_best_in_sequence(sequence, u12_bp_pwm, search_window_start=-55)

    # Should still find it, but with lower score (Ns get pseudocount)
    assert "TACTNNC" in match.sequence or match.sequence == "TACTNNC"
    assert match.score > 0  # Still scores > 0


def test_empty_intron_sequence(u12_bp_pwm, u2_bp_pwm):
    """Test handling of intron with no sequence."""
    intron = Intron(
        intron_id="test_intron_empty",
        coordinates=GenomicCoordinate(
            chromosome="chr1", start=1000, stop=1100, strand="+", system="1-based"
        ),
        sequences=IntronSequences(
            seq=None,  # No sequence!
            upstream_flank="ACTG",
            downstream_flank="TGCA",
        ),
        scores=IntronScores(),
        metadata=IntronMetadata(parent="transcript_4", grandparent="gene_4"),
    )

    scorer = BranchPointScorer(u12_bp_pwm, u2_bp_pwm)

    # Should raise informative error
    with pytest.raises(ValueError, match="sequence|None|missing"):
        scorer.find_best_match(intron)


# ============================================================================
# Algorithm Verification Tests
# ============================================================================


def test_algorithm_matches_original_logic(u12_bp_pwm, u2_bp_pwm):
    """
    Verify that our implementation matches original bp_score logic.

    Original algorithm (intronIC.py:2143-2178):
    1. Use sliding window of matrix length
    2. Score each window with seq_score
    3. Track best score, coords, and sequence
    4. Return all three
    """
    scorer = BranchPointScorer(u12_bp_pwm, u2_bp_pwm)

    # Test sequence with known structure
    sequence = "AAA" + "TACTAAC" + "TTT"

    match = scorer._find_best_in_sequence(sequence, u12_bp_pwm, search_window_start=-55)

    # Verify structure matches original return signature
    assert hasattr(match, "score")
    assert hasattr(match, "sequence")
    assert hasattr(match, "start_in_region")
    assert hasattr(match, "stop_in_region")

    # Verify it found the TACTAAC
    assert match.sequence == "TACTAAC"
    assert match.start_in_region == 3  # After "AAA"


def test_coordinates_are_relative_to_search_region(
    u12_bp_pwm, u2_bp_pwm, simple_intron
):
    """Test that returned coordinates are relative to the search region start."""
    scorer = BranchPointScorer(u12_bp_pwm, u2_bp_pwm)

    match = scorer.find_best_match(simple_intron, search_window=(-50, -20))

    # Coordinates should be relative to 3' end
    # TACTAAC is at position 70 in 100bp intron
    # Relative to 3' end: 70 - 100 = -30
    assert match.position == -30

    # start_in_region and stop_in_region are relative to the search region
    # Search region is positions 50-80 (30bp window)
    # TACTAAC at position 70-76 is at offset 20-26 in the search region
    assert match.start_in_region == 20
    assert match.stop_in_region == 27


# ============================================================================
# Performance Tests (Optional)
# ============================================================================


def test_long_search_region_performance(u12_bp_pwm, u2_bp_pwm):
    """Test that long search regions complete in reasonable time."""
    scorer = BranchPointScorer(u12_bp_pwm, u2_bp_pwm)

    # Very long sequence (1000bp)
    sequence = "N" * 1000

    # Should complete quickly (sliding window is efficient)
    match = scorer._find_best_in_sequence(sequence, u12_bp_pwm, search_window_start=-55)

    assert match is not None


# ============================================================================
# Integration with PWM Scoring
# ============================================================================


def test_uses_u12_pwm_for_scoring(u12_bp_pwm, u2_bp_pwm):
    """Test that scorer uses U12 PWM (not U2) for finding best match."""
    scorer = BranchPointScorer(u12_bp_pwm, u2_bp_pwm)

    # Sequence with perfect U12 BP
    sequence = "NNNNN" + "TACTAAC" + "NNN"

    match = scorer._find_best_in_sequence(sequence, u12_bp_pwm, search_window_start=-55)

    # Score should be high (matches U12 PWM well)
    # With our test PWM, perfect TACTAAC should score ~0.95^4 * 0.95^3 ≈ 0.7
    assert match.score > 0.5


def test_score_is_product_of_base_frequencies(u12_bp_pwm, u2_bp_pwm):
    """Verify score calculation matches PWM.score_sequence behavior."""
    scorer = BranchPointScorer(u12_bp_pwm, u2_bp_pwm)

    # Score a known sequence
    test_seq = "TACTAAC"

    # Score via BranchPointScorer
    match = scorer._find_best_in_sequence(test_seq, u12_bp_pwm, search_window_start=-55)

    # Score via PWM directly
    direct_score = u12_bp_pwm.score_sequence(test_seq)

    # Should match
    assert abs(match.score - direct_score) < 1e-10


# ============================================================================
# Streaming vs Standard Mode Consistency Tests
# ============================================================================


@pytest.mark.parametrize("length", [30, 40, 47, 55, 60, 70, 100, 150])
def test_bp_extraction_streaming_standard_consistency(u12_bp_pwm, u2_bp_pwm, length):
    """
    Test that BP extraction is identical between streaming and standard modes.

    Streaming mode uses Intron.extract_scoring_motifs() before scoring.
    Standard mode uses BranchPointScorer._extract_search_region() during scoring.

    Both should extract the same BP search region for identical introns.
    """
    five_coords = (-3, 9)
    bp_coords = (-55, -5)
    three_coords = (-6, 4)

    # Create test sequence with known pattern
    test_seq = "".join([chr(65 + (i % 26)) for i in range(length)])

    coord = GenomicCoordinate("test", 1, length, "+", "1-based")

    intron = Intron(
        intron_id=f"test_{length}bp",
        coordinates=coord,
        sequences=IntronSequences(
            seq=test_seq,
            upstream_flank="TTT",
            downstream_flank="CCC",
            five_prime_dnt="GT",
            three_prime_dnt="AG",
        ),
    )

    scorer = BranchPointScorer(u12_bp_pwm, u2_bp_pwm)

    # Method 1: extract_scoring_motifs (streaming mode)
    result_intron = intron.extract_scoring_motifs(five_coords, bp_coords, three_coords)
    streaming_bp = result_intron.motifs.bp_region

    # Method 2: BranchPointScorer._extract_search_region (standard mode)
    standard_bp, _ = scorer._extract_search_region(
        intron, bp_coords, five_coords, three_coords
    )

    # Both methods should produce identical BP regions
    assert streaming_bp == standard_bp, (
        f"BP extraction mismatch for {length}bp intron: "
        f"streaming={repr(streaming_bp[:30])}... vs standard={repr(standard_bp[:30])}..."
    )

    # Both methods should produce identical BP regions
    assert streaming_bp == standard_bp, (
        f"BP extraction mismatch for {length}bp intron: "
        f"streaming={repr(streaming_bp[:30])}... vs standard={repr(standard_bp[:30])}..."
    )
