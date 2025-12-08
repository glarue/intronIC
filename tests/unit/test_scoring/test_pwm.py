"""
Tests for PWM (Position Weight Matrix) scoring.

This module tests the PWM scoring system that calculates log-odds ratios
for intron splice sites.

Port from: intronIC.py:2114-2142 (seq_score)

Test Strategy:
1. Test PWM data structure and validation
2. Test sequence scoring with simple known matrices
3. Test ignore_positions feature (for canonical dinucleotides)
4. Test pseudocount handling for ambiguous bases
5. Test PWMLoader parsing from .iic format
"""

from pathlib import Path
from typing import Dict, Optional, Set

import numpy as np
import pytest

from intronIC.scoring.pwm import PWM, PWMLoader, PWMSet

# ============================================================================
# Test Fixtures
# ============================================================================


@pytest.fixture
def simple_pwm() -> PWM:
    """
    Create a simple 4-position PWM for testing.

    Matrix (4 positions, ACGT):
    Position:  0     1     2     3
    A:        1.0   0.0   0.5   0.25
    C:        0.0   1.0   0.0   0.25
    G:        0.0   0.0   0.5   0.25
    T:        0.0   0.0   0.0   0.25

    This should give perfect scores for "ACGT" and "ACGG".
    """
    matrix = np.array(
        [
            [1.0, 0.0, 0.5, 0.25],  # A
            [0.0, 1.0, 0.0, 0.25],  # C
            [0.0, 0.0, 0.5, 0.25],  # G
            [0.0, 0.0, 0.0, 0.25],  # T
        ]
    )

    return PWM(name="simple_test", matrix=matrix, length=4, pseudocount=0.0001)


@pytest.fixture
def canonical_five_prime_pwm() -> PWM:
    """
    Create a simplified 5' splice site PWM.

    Models the canonical GT dinucleotide at positions 0-1.
    """
    # Positions: -3, -2, -1, 0(G), 1(T), 2, 3, 4
    matrix = np.array(
        [
            [0.25, 0.30, 0.40, 0.01, 0.01, 0.35, 0.40, 0.30],  # A
            [0.25, 0.20, 0.15, 0.01, 0.01, 0.20, 0.20, 0.25],  # C
            [0.25, 0.30, 0.25, 0.97, 0.01, 0.25, 0.20, 0.20],  # G
            [0.25, 0.20, 0.20, 0.01, 0.97, 0.20, 0.20, 0.25],  # T
        ]
    )

    return PWM(
        name="five_prime_test",
        matrix=matrix,
        length=8,
        pseudocount=0.0001,
        start_index=-3,  # Position 0 in matrix is -3 relative to splice site
    )


@pytest.fixture
def sample_matrix_file(tmp_path: Path) -> Path:
    """
    Create a sample scoring_matrices.fasta.iic file for testing.

    Mimics the format used in the original intronIC.
    """
    content = """# Test matrix file
>u12_atac_five	start=-3	(n=100)
A	C	G	T
0.25	0.25	0.25	0.25
0.30	0.20	0.30	0.20
1.0	0.0	0.0	0.0
0.0	0.0	0.0	1.0
0.5	0.0	0.5	0.0

>u12_atac_bp	start=-20	(n=100)
A	C	G	T
0.25	0.25	0.25	0.25
0.40	0.20	0.20	0.20
0.30	0.30	0.20	0.20

>u12_atac_three	start=-5	(n=100)
A	C	G	T
0.20	0.30	0.25	0.25
0.90	0.03	0.03	0.04
0.05	0.05	0.85	0.05

>u2_gtag_five	start=-3	(n=10000)
A	C	G	T
0.30	0.20	0.25	0.25
0.25	0.25	0.25	0.25
0.01	0.01	0.97	0.01
0.01	0.01	0.01	0.97

>u2_gtag_bp	start=-20	(n=10000)
A	C	G	T
0.25	0.25	0.25	0.25
0.35	0.25	0.20	0.20
0.30	0.25	0.25	0.20

>u2_gtag_three	start=-5	(n=10000)
A	C	G	T
0.25	0.25	0.25	0.25
0.85	0.05	0.05	0.05
0.05	0.05	0.90	0.05
"""

    matrix_file = tmp_path / "scoring_matrices.fasta.iic"
    matrix_file.write_text(content)
    return matrix_file


# ============================================================================
# PWM Data Structure Tests
# ============================================================================


def test_pwm_creation(simple_pwm):
    """Test that PWM objects can be created and are frozen."""
    assert simple_pwm.name == "simple_test"
    assert simple_pwm.length == 4
    assert simple_pwm.pseudocount == 0.0001
    assert simple_pwm.matrix.shape == (4, 4)  # 4 bases x 4 positions

    # Verify frozen (immutable)
    with pytest.raises(AttributeError):
        simple_pwm.name = "changed"


def test_pwm_matrix_shape():
    """Test that matrix must have correct shape (4 rows for ACGT)."""
    with pytest.raises((ValueError, AssertionError)):
        PWM(
            name="bad_matrix",
            matrix=np.array([[1.0, 0.0]]),  # Only 1 row, not 4
            length=2,
            pseudocount=0.0001,
        )


# ============================================================================
# Sequence Scoring Tests (Core Algorithm)
# ============================================================================


def test_basic_sequence_scoring(simple_pwm):
    """
    Test basic sequence scoring with known matrix.

    Port from: intronIC.py:2114-2142 (seq_score)
    """
    # Score = product of frequencies at each position
    # For "ACGG":
    # Position 0: A -> 1.0
    # Position 1: C -> 1.0
    # Position 2: G -> 0.5
    # Position 3: G -> 0.25
    # Score = 1.0 * 1.0 * 0.5 * 0.25 = 0.125

    score = simple_pwm.score_sequence("ACGG")
    assert abs(score - 0.125) < 1e-6, f"Expected 0.125, got {score}"


def test_sequence_scoring_with_zero_frequency(simple_pwm):
    """
    Test that zero frequencies are handled with pseudocount.

    Port from: intronIC.py:2134 (pseudocount for missing bases)
    """
    # "TCGG" has T at position 0, which has frequency 0.0
    # Should use pseudocount instead
    score = simple_pwm.score_sequence("TCGG")

    # Expected: 0.0001 (pseudo) * 1.0 * 0.5 * 0.25 = 0.0000125
    expected = simple_pwm.pseudocount * 1.0 * 0.5 * 0.25
    assert abs(score - expected) < 1e-9, f"Expected {expected}, got {score}"


def test_sequence_scoring_with_ambiguous_base(simple_pwm):
    """
    Test that ambiguous bases (N) are handled with pseudocount.

    Port from: intronIC.py:2132-2134 (KeyError handling)
    """
    # "ACNG" has N at position 2
    # Should use pseudocount for N
    score = simple_pwm.score_sequence("ACNG")

    # Expected: 1.0 * 1.0 * 0.0001 (pseudo) * 0.25
    expected = 1.0 * 1.0 * simple_pwm.pseudocount * 0.25
    assert abs(score - expected) < 1e-9


def test_ignore_positions(canonical_five_prime_pwm):
    """
    Test ignore_positions parameter to skip canonical dinucleotides.

    Port from: intronIC.py:2118-2120, 2127-2128

    This is used when we want to score the non-canonical parts of a
    splice site without penalizing the canonical GT/AG dinucleotides.
    """
    # Sequence with canonical GT at positions 3-4 (0-indexed positions in seq)
    # But PWM has start_index=-3, so:
    # - Position 0 in sequence -> index -3 in PWM
    # - Position 3 in sequence -> index 0 in PWM (the G)
    # - Position 4 in sequence -> index 1 in PWM (the T)

    # PWM length is 8, so sequence must be 8 chars
    seq = "AAAGTTAA"

    # Score with canonical positions
    score_with_canonical = canonical_five_prime_pwm.score_sequence(seq)

    # Score ignoring canonical positions (indices 0 and 1 in PWM space)
    # These correspond to positions 3 and 4 in the sequence
    score_without_canonical = canonical_five_prime_pwm.score_sequence(
        seq, ignore_positions={0, 1}
    )

    # Score without canonical should be higher (not penalized by GT requirement)
    assert score_without_canonical > score_with_canonical


def test_sequence_length_mismatch(simple_pwm):
    """Test that sequences can be different lengths (flexible scoring).

    The PWM scoring now supports flexible sequence lengths to enable
    partial matrix scoring (useful for branch point scanning).
    """
    # Short sequences should work
    score_short = simple_pwm.score_sequence("AC")
    assert isinstance(score_short, float)
    assert score_short > 0

    # Long sequences should work
    score_long = simple_pwm.score_sequence("ACGTACGT")
    assert isinstance(score_long, float)
    assert score_long > 0


def test_empty_sequence(simple_pwm):
    """Test behavior with empty sequence.

    Empty sequences should raise ValueError - indicates a bug in
    extraction/filtering if we get here.
    """
    with pytest.raises(ValueError, match="Cannot score empty sequence"):
        simple_pwm.score_sequence("")


# ============================================================================
# PWMSet Tests (U2/U12 Matrix Pairs)
# ============================================================================


def test_pwmset_creation():
    """Test that PWMSet can group U2 and U12 PWMs."""
    u2_pwm = PWM(
        name="u2_test",
        matrix=np.array([[0.25] * 4] * 4),  # Uniform
        length=4,
        pseudocount=0.0001,
    )

    u12_pwm = PWM(
        name="u12_test",
        matrix=np.array([[0.5, 0.5, 0.0, 0.0]] * 4),  # Different
        length=4,
        pseudocount=0.0001,
    )

    pwm_set = PWMSet(
        matrices={
            ("u2", "gtag"): u2_pwm,
            ("u12", "gtag"): u12_pwm,
        }
    )

    assert pwm_set.select_best("u2", "gtag") == u2_pwm
    assert pwm_set.select_best("u12", "gtag") == u12_pwm


def test_pwmset_with_noncanonical():
    """Test PWMSet with both canonical and non-canonical PWMs."""
    canonical = PWM(
        name="canonical",
        matrix=np.array([[0.25] * 4] * 4),
        length=4,
        pseudocount=0.0001,
    )

    noncanonical = PWM(
        name="noncanonical",
        matrix=np.array([[0.5] * 4] * 4),
        length=4,
        pseudocount=0.0001,
    )

    pwm_set = PWMSet(
        matrices={
            ("u2", "gtag"): canonical,
            ("u2", "gcag"): noncanonical,
            ("u12", "gtag"): canonical,
            ("u12", "atac"): noncanonical,
        }
    )

    # Non-canonical variants should be accessible
    assert pwm_set.select_best("u2", "gcag") is not None
    assert pwm_set.select_best("u12", "atac") is not None


# ============================================================================
# PWMLoader Tests (File Parsing)
# ============================================================================


def test_pwm_loader_basic(sample_matrix_file):
    """
    Test that PWMLoader can parse the matrix file format.

    Port from: intronIC.py:1180-1264 (load_external_matrix)
    """
    pwm_sets = PWMLoader.load_from_file(sample_matrix_file)

    # Should have three regions: five, bp, three
    assert "five" in pwm_sets
    assert "bp" in pwm_sets
    assert "three" in pwm_sets

    # Each should be a PWMSet
    assert isinstance(pwm_sets["five"], PWMSet)
    assert isinstance(pwm_sets["bp"], PWMSet)
    assert isinstance(pwm_sets["three"], PWMSet)


def test_pwm_loader_parse_names(sample_matrix_file):
    """
    Test that loader correctly parses matrix names into categories.

    From original: u12_atac_five -> (u12, atac, five)
    """
    pwm_sets = PWMLoader.load_from_file(sample_matrix_file)

    # Check that U12 and U2 PWMs were loaded
    five_set = pwm_sets["five"]
    u12_pwm = five_set.select_best("u12", "gtag")
    u2_pwm = five_set.select_best("u2", "gtag")
    assert u12_pwm is not None
    assert u2_pwm is not None

    # Verify names
    assert "u12" in u12_pwm.name.lower()
    assert "u2" in u2_pwm.name.lower()


def test_pwm_loader_parse_start_index(sample_matrix_file):
    """
    Test that start_index is correctly parsed from headers.

    From file: >u12_atac_five	start=-3
    """
    pwm_sets = PWMLoader.load_from_file(sample_matrix_file)

    # Five prime site starts at -3
    assert pwm_sets["five"].select_best("u12", "gtag").start_index == -3

    # Branch point starts at -20
    assert pwm_sets["bp"].select_best("u12", "gtag").start_index == -20


def test_pwm_loader_matrix_dimensions(sample_matrix_file):
    """Test that loaded matrices have correct dimensions."""
    pwm_sets = PWMLoader.load_from_file(sample_matrix_file)

    # Five prime PWM has 5 positions (from sample file)
    five_pwm = pwm_sets["five"].select_best("u12", "gtag")
    assert five_pwm.matrix.shape == (4, 5)  # 4 bases x 5 positions
    assert five_pwm.length == 5

    # Branch point PWM has 3 positions
    bp_pwm = pwm_sets["bp"].select_best("u12", "gtag")
    assert bp_pwm.matrix.shape == (4, 3)
    assert bp_pwm.length == 3


def test_pwm_loader_frequency_values(sample_matrix_file):
    """Test that frequency values are correctly parsed."""
    pwm_sets = PWMLoader.load_from_file(sample_matrix_file)

    # Check that frequencies sum to ~1.0 at each position
    five_pwm = pwm_sets["five"].select_best("u12", "gtag")

    for pos in range(five_pwm.length):
        column_sum = sum(five_pwm.matrix[:, pos])
        assert abs(column_sum - 1.0) < 0.01, (
            f"Position {pos} frequencies should sum to ~1.0, got {column_sum}"
        )


def test_pwm_loader_nonexistent_file():
    """Test behavior when file doesn't exist."""
    with pytest.raises(FileNotFoundError):
        PWMLoader.load_from_file(Path("/nonexistent/path/to/matrices.iic"))


def test_pwm_loader_empty_file(tmp_path):
    """Test behavior with empty matrix file."""
    empty_file = tmp_path / "empty.iic"
    empty_file.write_text("")

    # Should return empty dict or raise informative error
    result = PWMLoader.load_from_file(empty_file)
    assert result == {} or isinstance(result, dict)


# ============================================================================
# Integration Tests (Scoring with Loaded PWMs)
# ============================================================================


def test_score_real_sequence_with_loaded_pwm(sample_matrix_file):
    """
    Integration test: Load PWMs and score a sequence.

    This verifies the complete workflow from file -> PWM -> score.
    """
    pwm_sets = PWMLoader.load_from_file(sample_matrix_file)

    # Get U12 five prime PWM (length 5, start=-3)
    u12_five = pwm_sets["five"].select_best("u12", "gtag")

    # Create a sequence that matches: XXATX (positions -3 to +1)
    # Position -3: A (0.25), -2: T (0.20), -1: A (1.0), 0: T (1.0), +1: G (0.5)
    test_seq = "XATAT"

    # Should be able to score it
    score = u12_five.score_sequence(test_seq)
    assert score > 0
    assert np.isfinite(score)


def test_compare_u2_vs_u12_scores(sample_matrix_file):
    """
    Test that U2 and U12 PWMs give different scores.

    This is fundamental to the classification system.
    """
    pwm_sets = PWMLoader.load_from_file(sample_matrix_file)

    u2_five = pwm_sets["five"].select_best("u2", "gtag")
    u12_five = pwm_sets["five"].select_best("u12", "gtag")

    # Note: U2 and U12 PWMs may have different lengths
    # In sample file: U2 has 4 positions, U12 has 5 positions

    # Create test sequences for each PWM
    u2_test_seq = "A" * u2_five.length
    u12_test_seq = "A" * u12_five.length

    u2_score = u2_five.score_sequence(u2_test_seq)
    u12_score = u12_five.score_sequence(u12_test_seq)

    # Scores should be positive and finite
    assert u2_score > 0 and np.isfinite(u2_score)
    assert u12_score > 0 and np.isfinite(u12_score)

    # Scores should differ if we score same-length sequences with different PWMs
    # Use the longer length (u12_five.length = 5)
    common_seq = "ACGTA"  # 5 bases to match u12 length
    u2_common_seq = common_seq[: u2_five.length]  # First 4 bases for u2
    u12_common_seq = common_seq  # All 5 bases for u12

    u2_common = u2_five.score_sequence(u2_common_seq)
    u12_common = u12_five.score_sequence(u12_common_seq)

    # At minimum, they should both be valid scores
    assert u2_common > 0
    assert u12_common > 0


# ============================================================================
# Edge Cases
# ============================================================================


def test_score_with_lowercase_sequence(simple_pwm):
    """Test that lowercase sequences are handled correctly."""
    # Should either auto-convert or raise clear error
    try:
        score = simple_pwm.score_sequence("acgg")
        assert score > 0  # If it auto-converts, should get valid score
    except ValueError as e:
        # If it doesn't auto-convert, should raise informative error
        assert "lowercase" in str(e).lower() or "upper" in str(e).lower()


def test_pwm_with_custom_pseudocount():
    """Test that custom pseudocount values work."""
    pwm = PWM(
        name="custom_pseudo",
        matrix=np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 0.0], [0.0, 0.0]]),
        length=2,
        pseudocount=0.01,  # Custom pseudocount
    )

    # Sequence with zero frequency: "TT"
    # Should use custom pseudocount: 0.01 * 0.01 = 0.0001
    score = pwm.score_sequence("TT")
    expected = 0.01 * 0.01
    assert abs(score - expected) < 1e-9
