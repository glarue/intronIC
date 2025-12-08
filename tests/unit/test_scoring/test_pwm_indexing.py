"""
Critical tests for PWM position indexing.

This module tests the critical fix for PWM position indexing that was causing
incorrect scores. The bug was fixed in commit c4244b1.

Key Issue:
    PWM.score_sequence() must correctly map sequence positions to matrix indices.
    The calculation is: matrix_index = (seq_start_pos + i) - start_index

    For example:
    - PWM has start_index=-20 (matrix position 0 represents logical position -20)
    - Sequence starts at logical position -3
    - First base (i=0) should map to matrix index 17:
      matrix_index = (-3 + 0) - (-20) = 17

Test Strategy:
    1. Test basic position mapping with explicit calculations
    2. Test edge cases (start_index=0, negative start_index)
    3. Test that sequences score consistently regardless of how they're split
    4. Test against known correct scores from original intronIC

Port from: commit c4244b1 (Fix critical PWM scoring bugs)
"""

import pytest
import numpy as np
import math
from pathlib import Path

from intronIC.scoring.pwm import PWM, PWMLoader


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def simple_pwm():
    """
    Create simple PWM for testing position mapping.

    Matrix starts at position -3 and has length 8 (positions -3 to +4).
    Each position has a unique pattern to verify correct indexing.
    """
    # Create distinctive pattern: each position favors a different base
    matrix = np.array([
        [0.7, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],  # A: high at pos -3
        [0.1, 0.7, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],  # C: high at pos -2
        [0.1, 0.1, 0.7, 0.1, 0.1, 0.1, 0.1, 0.1],  # G: high at pos -1
        [0.1, 0.1, 0.1, 0.7, 0.1, 0.1, 0.1, 0.1],  # T: high at pos 0
    ])

    return PWM(
        name="test_pwm",
        matrix=matrix,
        length=8,
        start_index=-3
    )


@pytest.fixture
def zero_start_pwm():
    """PWM starting at position 0."""
    matrix = np.array([
        [0.8, 0.1, 0.1, 0.1],  # A: high at pos 0
        [0.1, 0.8, 0.1, 0.1],  # C: high at pos 1
        [0.1, 0.1, 0.8, 0.1],  # G: high at pos 2
        [0.1, 0.1, 0.1, 0.8],  # T: high at pos 3
    ])

    return PWM(
        name="zero_start",
        matrix=matrix,
        length=4,
        start_index=0
    )


# ============================================================================
# Basic Position Mapping Tests
# ============================================================================

def test_position_mapping_with_negative_start(simple_pwm):
    """Test position mapping when PWM has negative start_index."""
    # Sequence starting at position -3 (same as PWM start)
    # Sequence: A at -3, C at -2, G at -1, T at 0
    seq = "ACGT"
    score = simple_pwm.score_sequence(seq, seq_start_position=-3)

    # Each base should hit its high-probability position
    expected_score = 0.7 * 0.7 * 0.7 * 0.7
    assert abs(score - expected_score) < 1e-10, \
        f"Score {score} doesn't match expected {expected_score}"


def test_position_mapping_with_offset(simple_pwm):
    """Test position mapping when sequence starts at different position."""
    # Sequence starting at position -2 (one position after PWM start)
    # Sequence: C at -2, G at -1, T at 0, A at +1
    seq = "CGTA"
    score = simple_pwm.score_sequence(seq, seq_start_position=-2)

    # C at -2: matrix[1, 1] = 0.7
    # G at -1: matrix[2, 2] = 0.7
    # T at 0: matrix[3, 3] = 0.7
    # A at +1: matrix[0, 4] = 0.1 (any base at position 4 has 0.1)
    expected_score = 0.7 * 0.7 * 0.7 * 0.1
    assert abs(score - expected_score) < 1e-10


def test_position_mapping_with_zero_start(zero_start_pwm):
    """Test position mapping when PWM starts at 0."""
    # Sequence starting at position 0
    seq = "ACGT"
    score = zero_start_pwm.score_sequence(seq, seq_start_position=0)

    # Each base should hit its high-probability position
    expected_score = 0.8 ** 4
    assert abs(score - expected_score) < 1e-10


def test_explicit_matrix_index_calculation(simple_pwm):
    """
    Test explicit matrix index calculation for each position.

    This verifies the formula: matrix_index = (seq_start_pos + i) - start_index
    """
    seq = "ACGT"
    seq_start_pos = -3
    start_index = simple_pwm.start_index  # -3

    # Manually calculate expected indices
    expected_indices = []
    for i in range(len(seq)):
        logical_pos = seq_start_pos + i  # -3, -2, -1, 0
        matrix_index = logical_pos - start_index  # 0, 1, 2, 3
        expected_indices.append(matrix_index)

    assert expected_indices == [0, 1, 2, 3], \
        "Matrix index calculation is incorrect"

    # Verify score uses these indices correctly
    score = simple_pwm.score_sequence(seq, seq_start_position=seq_start_pos)
    expected_score = (
        simple_pwm.matrix[0, 0] *  # A at index 0
        simple_pwm.matrix[1, 1] *  # C at index 1
        simple_pwm.matrix[2, 2] *  # G at index 2
        simple_pwm.matrix[3, 3]    # T at index 3
    )
    assert abs(score - expected_score) < 1e-10


# ============================================================================
# Edge Case Tests
# ============================================================================

def test_sequence_beyond_pwm_bounds(simple_pwm):
    """Test that bases outside PWM length are handled correctly."""
    # Sequence longer than PWM, starting before PWM
    seq = "AAACGT"  # Sequence starting before PWM starts
    score = simple_pwm.score_sequence(seq, seq_start_position=-5)

    # Positions -5, -4 are before PWM (start_index=-3), so skipped
    # Position -3 maps to seq[2]='A' at matrix index 0 (expects A, gets 0.7)
    # Position -2 maps to seq[3]='C' at matrix index 1 (expects C, gets 0.7)
    # Position -1 maps to seq[4]='G' at matrix index 2 (expects G, gets 0.7)
    # Position 0 maps to seq[5]='T' at matrix index 3 (expects T, gets 0.7)
    # All four positions match perfectly
    expected_score = 0.7 ** 4
    assert abs(score - expected_score) < 1e-10


def test_sequence_partially_overlapping_pwm(simple_pwm):
    """Test sequence that only partially overlaps PWM."""
    # Sequence starting at position +2 (near end of PWM)
    # PWM covers -3 to +4, so position +2 to +5 partially overlaps
    seq = "ACGT"
    score = simple_pwm.score_sequence(seq, seq_start_position=2)

    # Position +2 maps to matrix index 5: matrix[:, 5] all have 0.1
    # Position +3 maps to matrix index 6: matrix[:, 6] all have 0.1
    # Position +4 maps to matrix index 7: matrix[:, 7] all have 0.1
    # Position +5 is beyond PWM, should be skipped
    expected_score = 0.1 * 0.1 * 0.1  # Three positions at 0.1 each
    assert abs(score - expected_score) < 1e-10


# ============================================================================
# Consistency Tests
# ============================================================================

def test_scoring_split_sequences(simple_pwm):
    """Test that splitting a sequence doesn't change total score."""
    # Score full sequence
    full_seq = "ACGTACGT"
    full_score = simple_pwm.score_sequence(full_seq, seq_start_position=-3)

    # Score in two parts
    part1 = "ACGT"
    part2 = "ACGT"
    score1 = simple_pwm.score_sequence(part1, seq_start_position=-3)
    score2 = simple_pwm.score_sequence(part2, seq_start_position=1)

    # Product of parts should equal full sequence score
    combined_score = score1 * score2
    assert abs(full_score - combined_score) < 1e-10, \
        "Split scoring doesn't match full sequence scoring"


def test_subsequence_consistency(simple_pwm):
    """Test that subsequence scores are consistent with full sequence."""
    full_seq = "ACGTAA"  # Use valid bases instead of N
    full_score = simple_pwm.score_sequence(full_seq, seq_start_position=-3)

    # Score just the part that contributes
    subseq = "ACGT"
    sub_score = simple_pwm.score_sequence(subseq, seq_start_position=-3)

    # The AA part at the end should contribute matrix values at positions 4,5
    # Both A's at positions beyond the first 4 have probability 0.1
    aa_contribution = 0.1 * 0.1  # Both positions have 0.1 for A
    expected_full = sub_score * aa_contribution

    assert abs(full_score - expected_full) < 1e-10


# ============================================================================
# Real PWM Tests
# ============================================================================

@pytest.mark.slow
def test_five_prime_scoring_with_real_pwm(matrix_file):
    """
    Test 5' scoring with real PWM against known correct value.

    This is the critical test from commit c4244b1 that revealed the bug.
    The sequence TCAGTATCCTTC at position -3 should produce log ratio 18.22.
    """
    # Load real PWMs
    if not matrix_file.exists():
        pytest.skip("Real PWM file not available")

    loader = PWMLoader()
    pwm_sets = loader.load_from_file(matrix_file)

    # Get U12 GT-AG 5' PWM (most U12 introns are GT-AG)
    u12_pwm = pwm_sets['five'].select_best('u12', 'gtag')
    u2_pwm = pwm_sets['five'].select_best('u2', 'gtag')

    # Test sequence from ENST00000359435_5 (known U12 intron)
    test_seq = "TCAGTATCCTTC"
    seq_start_position = -3  # 5' region starts 3bp before intron

    # Score with both PWMs
    u12_score = u12_pwm.score_sequence(test_seq, seq_start_position=seq_start_position)
    u2_score = u2_pwm.score_sequence(test_seq, seq_start_position=seq_start_position)

    # Calculate log ratio
    log_ratio = math.log2(u12_score / u2_score)

    # Expected value from original intronIC
    expected_log_ratio = 18.216262900453128

    # Should match within reasonable tolerance
    assert abs(log_ratio - expected_log_ratio) < 0.01, \
        f"Log ratio {log_ratio} doesn't match expected {expected_log_ratio}"


# ============================================================================
# Ignore Positions Tests
# ============================================================================

def test_ignore_positions_uses_logical_positions(simple_pwm):
    """Test that ignore_positions uses logical positions, not matrix indices."""
    seq = "ACGT"
    seq_start_pos = -3

    # Ignore logical position -2 (which maps to matrix index 1)
    score_with_ignore = simple_pwm.score_sequence(
        seq,
        seq_start_position=seq_start_pos,
        ignore_positions={-2}
    )

    # Score without ignore
    score_without = simple_pwm.score_sequence(
        seq,
        seq_start_position=seq_start_pos
    )

    # With ignore: A(0.7) * 1.0 * G(0.7) * T(0.7) = 0.343
    # Without: A(0.7) * C(0.7) * G(0.7) * T(0.7) = 0.2401
    # Ratio should be 1/0.7 ≈ 1.43
    ratio = score_with_ignore / score_without
    expected_ratio = 1.0 / 0.7
    assert abs(ratio - expected_ratio) < 1e-10, \
        "Ignore positions not working correctly"


def test_ignore_multiple_positions(simple_pwm):
    """Test ignoring multiple positions."""
    seq = "ACGT"
    score = simple_pwm.score_sequence(
        seq,
        seq_start_position=-3,
        ignore_positions={-3, -1}  # Ignore first and third positions
    )

    # Only C at -2 and T at 0 contribute
    expected_score = 0.7 * 0.7
    assert abs(score - expected_score) < 1e-10


# ============================================================================
# Regression Tests
# ============================================================================

def test_position_calculation_example_from_bug_fix():
    """
    Test the specific example from the bug fix documentation.

    Example: PWM with start_index=-20, sequence starting at -3
    First base (i=0) should map to matrix index 17.
    """
    # Create PWM with start_index=-20 and length 40
    matrix = np.full((4, 40), 0.25)  # Uniform for simplicity
    matrix[0, 17] = 0.9  # Make position 17 special

    pwm = PWM(
        name="test",
        matrix=matrix,
        length=40,
        start_index=-20
    )

    # Sequence starting at logical position -3
    # First base should hit matrix index 17
    seq = "A"  # Just one base to test
    score = pwm.score_sequence(seq, seq_start_position=-3)

    # Should get the special value at matrix[0, 17]
    assert abs(score - 0.9) < 1e-10, \
        f"Matrix index calculation incorrect: got score {score}, expected 0.9"


def test_bug_would_have_used_wrong_index():
    """
    Test that demonstrates what the bug would have done wrong.

    The bug was using: pwm_position = i + self.start_index
    Correct is: matrix_index = (seq_start_pos + i) - self.start_index
    """
    matrix = np.full((4, 40), 0.25)
    matrix[0, 0] = 0.9   # Correct index for logical position -20
    matrix[0, 17] = 0.8  # Correct index for logical position -3

    pwm = PWM(name="test", matrix=matrix, length=40, start_index=-20)

    # Sequence at position -3, first base is A
    # Bug would calculate: pwm_position = 0 + (-20) = -20 → matrix index -20 (crash or wrap)
    # Correct calculates: matrix_index = (-3 + 0) - (-20) = 17
    seq = "A"
    score = pwm.score_sequence(seq, seq_start_position=-3)

    # Should use matrix index 17, not 0 or -20
    assert abs(score - 0.8) < 1e-10, \
        "Appears to be using incorrect index calculation (like the original bug)"
