"""
Deterministic tests for log-ratio score calculations.

This module provides mathematically straightforward tests with predefined PWMs
and sequences to verify that log-ratio calculations are correct. Any changes to
the scoring logic should cause these tests to fail if the behavior changes.

Log-ratio formula: log2(u12_score / u2_score)

Test Strategy:
    1. Create simple PWMs with known probabilities
    2. Score sequences with known outcomes
    3. Calculate expected log ratios mathematically
    4. Verify implementation matches mathematical expectation

This is critical for regression testing - these tests will catch any accidental
changes to the scoring pipeline.
"""

import pytest
import numpy as np
import math

from intronIC.scoring.pwm import PWM


# ============================================================================
# Simple Deterministic Tests
# ============================================================================

def test_log_ratio_with_identical_scores():
    """Test log ratio when U12 and U2 scores are identical."""
    # Create identical PWMs
    matrix = np.array([[0.25] * 4] * 4)  # Uniform distribution

    u12_pwm = PWM("u12", matrix, length=4)
    u2_pwm = PWM("u2", matrix, length=4)

    seq = "ACGT"
    u12_score = u12_pwm.score_sequence(seq)
    u2_score = u2_pwm.score_sequence(seq)

    # Scores should be identical
    assert abs(u12_score - u2_score) < 1e-10

    # Log ratio should be zero
    log_ratio = math.log2(u12_score / u2_score)
    assert abs(log_ratio - 0.0) < 1e-10, \
        f"Log ratio should be 0.0 for identical scores, got {log_ratio}"


def test_log_ratio_factor_of_two():
    """Test log ratio when U12 score is exactly 2x U2 score."""
    # U12 matrix: all positions favor A (0.5)
    u12_matrix = np.array([
        [0.5, 0.5, 0.5, 0.5],  # A
        [0.2, 0.2, 0.2, 0.2],  # C
        [0.2, 0.2, 0.2, 0.2],  # G
        [0.1, 0.1, 0.1, 0.1],  # T
    ])

    # U2 matrix: all positions have A at 0.25 (half of U12)
    u2_matrix = np.array([
        [0.25, 0.25, 0.25, 0.25],  # A
        [0.25, 0.25, 0.25, 0.25],  # C
        [0.25, 0.25, 0.25, 0.25],  # G
        [0.25, 0.25, 0.25, 0.25],  # T
    ])

    u12_pwm = PWM("u12", u12_matrix, length=4)
    u2_pwm = PWM("u2", u2_matrix, length=4)

    # Sequence of all A's
    seq = "AAAA"
    u12_score = u12_pwm.score_sequence(seq)  # 0.5^4 = 0.0625
    u2_score = u2_pwm.score_sequence(seq)     # 0.25^4 = 0.00390625

    # Ratio should be exactly 16 (0.5/0.25 = 2, raised to 4th power = 16)
    ratio = u12_score / u2_score
    expected_ratio = (0.5 / 0.25) ** 4
    assert abs(ratio - expected_ratio) < 1e-10

    # Log2(16) = 4
    log_ratio = math.log2(u12_score / u2_score)
    expected_log_ratio = 4.0
    assert abs(log_ratio - expected_log_ratio) < 1e-10, \
        f"Expected log ratio 4.0, got {log_ratio}"


def test_log_ratio_factor_of_four():
    """Test log ratio when U12 score is 4x U2 score."""
    # U12: all positions favor G (0.4)
    u12_matrix = np.array([
        [0.2, 0.2, 0.2, 0.2],  # A
        [0.2, 0.2, 0.2, 0.2],  # C
        [0.4, 0.4, 0.4, 0.4],  # G
        [0.2, 0.2, 0.2, 0.2],  # T
    ])

    # U2: all positions have G at 0.1 (1/4 of U12)
    u2_matrix = np.array([
        [0.3, 0.3, 0.3, 0.3],  # A
        [0.3, 0.3, 0.3, 0.3],  # C
        [0.1, 0.1, 0.1, 0.1],  # G
        [0.3, 0.3, 0.3, 0.3],  # T
    ])

    u12_pwm = PWM("u12", u12_matrix, length=4)
    u2_pwm = PWM("u2", u2_matrix, length=4)

    seq = "GGGG"
    u12_score = u12_pwm.score_sequence(seq)  # 0.4^4
    u2_score = u2_pwm.score_sequence(seq)     # 0.1^4

    # Ratio should be (0.4/0.1)^4 = 4^4 = 256
    ratio = u12_score / u2_score
    expected_ratio = 256.0
    assert abs(ratio - expected_ratio) < 1e-6

    # Log2(256) = 8
    log_ratio = math.log2(u12_score / u2_score)
    expected_log_ratio = 8.0
    assert abs(log_ratio - expected_log_ratio) < 1e-10, \
        f"Expected log ratio 8.0, got {log_ratio}"


def test_log_ratio_with_negative_value():
    """Test log ratio when U2 score is higher than U12 score (negative log ratio)."""
    # U12: all positions disfavor T (0.1)
    u12_matrix = np.array([
        [0.3, 0.3, 0.3, 0.3],  # A
        [0.3, 0.3, 0.3, 0.3],  # C
        [0.3, 0.3, 0.3, 0.3],  # G
        [0.1, 0.1, 0.1, 0.1],  # T
    ])

    # U2: all positions favor T (0.4)
    u2_matrix = np.array([
        [0.2, 0.2, 0.2, 0.2],  # A
        [0.2, 0.2, 0.2, 0.2],  # C
        [0.2, 0.2, 0.2, 0.2],  # G
        [0.4, 0.4, 0.4, 0.4],  # T
    ])

    u12_pwm = PWM("u12", u12_matrix, length=4)
    u2_pwm = PWM("u2", u2_matrix, length=4)

    seq = "TTTT"
    u12_score = u12_pwm.score_sequence(seq)  # 0.1^4 = 0.0001
    u2_score = u2_pwm.score_sequence(seq)     # 0.4^4 = 0.0256

    # Ratio should be 0.1/0.4 = 0.25, raised to 4th = (1/256)
    ratio = u12_score / u2_score
    expected_ratio = (0.1 / 0.4) ** 4
    assert abs(ratio - expected_ratio) < 1e-10

    # Log2(1/256) = -8
    log_ratio = math.log2(u12_score / u2_score)
    expected_log_ratio = -8.0
    assert abs(log_ratio - expected_log_ratio) < 1e-10, \
        f"Expected log ratio -8.0, got {log_ratio}"


# ============================================================================
# Position-Specific Tests
# ============================================================================

def test_log_ratio_single_position_difference():
    """Test log ratio when only one position differs between PWMs."""
    # Both PWMs uniform except position 2
    u12_matrix = np.array([
        [0.25, 0.25, 0.8, 0.25],  # A: high at position 2
        [0.25, 0.25, 0.1, 0.25],  # C
        [0.25, 0.25, 0.05, 0.25],  # G
        [0.25, 0.25, 0.05, 0.25],  # T
    ])

    u2_matrix = np.array([
        [0.25, 0.25, 0.1, 0.25],  # A: low at position 2
        [0.25, 0.25, 0.3, 0.25],  # C
        [0.25, 0.25, 0.3, 0.25],  # G
        [0.25, 0.25, 0.3, 0.25],  # T
    ])

    u12_pwm = PWM("u12", u12_matrix, length=4)
    u2_pwm = PWM("u2", u2_matrix, length=4)

    # Sequence with A at position 2
    seq = "NNAG"  # Position 2 is 'A'
    u12_score = u12_pwm.score_sequence(seq)
    u2_score = u2_pwm.score_sequence(seq)

    # Calculate expected ratio
    # Positions 0, 1, 3 contribute equally (0.25 each for both PWMs)
    # Position 2: U12 has 0.8, U2 has 0.1
    # Net ratio is 0.8/0.1 = 8
    expected_ratio = 0.8 / 0.1
    actual_ratio = u12_score / u2_score
    assert abs(actual_ratio - expected_ratio) < 1e-10

    log_ratio = math.log2(actual_ratio)
    expected_log_ratio = math.log2(8)  # = 3.0
    assert abs(log_ratio - expected_log_ratio) < 1e-10


# ============================================================================
# Multi-Position Tests
# ============================================================================

def test_log_ratio_all_positions_contribute():
    """Test log ratio where all positions contribute differently."""
    # Create PWMs where each position has a 2x difference
    u12_matrix = np.array([
        [0.4, 0.4, 0.4, 0.4],  # A
        [0.2, 0.2, 0.2, 0.2],  # C
        [0.2, 0.2, 0.2, 0.2],  # G
        [0.2, 0.2, 0.2, 0.2],  # T
    ])

    u2_matrix = np.array([
        [0.2, 0.2, 0.2, 0.2],  # A: half of U12
        [0.3, 0.3, 0.3, 0.3],  # C
        [0.3, 0.3, 0.3, 0.3],  # G
        [0.2, 0.2, 0.2, 0.2],  # T
    ])

    u12_pwm = PWM("u12", u12_matrix, length=4)
    u2_pwm = PWM("u2", u2_matrix, length=4)

    seq = "AAAA"
    u12_score = u12_pwm.score_sequence(seq)  # 0.4^4
    u2_score = u2_pwm.score_sequence(seq)     # 0.2^4

    # Each position contributes a factor of 2
    # Total ratio is 2^4 = 16
    expected_ratio = 16.0
    actual_ratio = u12_score / u2_score
    assert abs(actual_ratio - expected_ratio) < 1e-10

    log_ratio = math.log2(actual_ratio)
    expected_log_ratio = 4.0  # log2(16) = 4
    assert abs(log_ratio - expected_log_ratio) < 1e-10


# ============================================================================
# Start Position Tests
# ============================================================================

def test_log_ratio_with_start_position():
    """Test log ratio calculation with non-zero start position."""
    # PWMs with start_index = -3
    u12_matrix = np.array([
        [0.8, 0.2, 0.2, 0.2],  # A: high at position -3 (matrix index 0)
        [0.1, 0.1, 0.1, 0.1],  # C
        [0.05, 0.05, 0.05, 0.05],  # G
        [0.05, 0.05, 0.05, 0.05],  # T
    ])

    u2_matrix = np.array([
        [0.2, 0.2, 0.2, 0.2],  # A: low everywhere
        [0.3, 0.3, 0.3, 0.3],  # C
        [0.3, 0.3, 0.3, 0.3],  # G
        [0.2, 0.2, 0.2, 0.2],  # T
    ])

    u12_pwm = PWM("u12", u12_matrix, length=4, start_index=-3)
    u2_pwm = PWM("u2", u2_matrix, length=4, start_index=-3)

    # Sequence starting at position -3
    seq = "ACGT"
    u12_score = u12_pwm.score_sequence(seq, seq_start_position=-3)
    u2_score = u2_pwm.score_sequence(seq, seq_start_position=-3)

    # Calculate expected scores based on actual implementation:
    # Position -3 (A, matrix idx 0): u12=0.8, u2=0.2
    # Position -2 (C, matrix idx 1): u12=0.1, u2=0.3
    # Position -1 (G, matrix idx 2): u12=0.05, u2=0.3
    # Position 0 (T, matrix idx 3): u12=0.05, u2=0.2
    expected_u12 = 0.8 * 0.1 * 0.05 * 0.05
    expected_u2 = 0.2 * 0.3 * 0.3 * 0.2

    assert abs(u12_score - expected_u12) < 1e-10
    assert abs(u2_score - expected_u2) < 1e-10

    actual_ratio = u12_score / u2_score
    expected_ratio = expected_u12 / expected_u2
    assert abs(actual_ratio - expected_ratio) < 1e-10

    log_ratio = math.log2(actual_ratio)
    expected_log_ratio = math.log2(expected_ratio)
    assert abs(log_ratio - expected_log_ratio) < 1e-10


# ============================================================================
# Regression Test with Known Values
# ============================================================================

def test_log_ratio_exact_calculation():
    """
    Test exact log ratio calculation with specific numerical values.

    This serves as a regression test - any changes to the scoring logic
    that alter these specific values will be caught.
    """
    # Create PWMs with specific, easy-to-calculate values
    u12_matrix = np.array([
        [0.5, 0.125, 0.5, 0.25],  # A
        [0.25, 0.5, 0.25, 0.25],  # C
        [0.125, 0.25, 0.125, 0.25],  # G
        [0.125, 0.125, 0.125, 0.25],  # T
    ])

    u2_matrix = np.array([
        [0.25, 0.25, 0.25, 0.25],  # A
        [0.25, 0.25, 0.25, 0.25],  # C
        [0.25, 0.25, 0.25, 0.25],  # G
        [0.25, 0.25, 0.25, 0.25],  # T: uniform
    ])

    u12_pwm = PWM("u12", u12_matrix, length=4)
    u2_pwm = PWM("u2", u2_matrix, length=4)

    # Specific sequence
    seq = "ACAG"

    # Calculate expected scores manually:
    # U12: A(0.5) * C(0.5) * A(0.5) * G(0.25) = 0.03125
    # U2: 0.25 * 0.25 * 0.25 * 0.25 = 0.00390625
    expected_u12 = 0.5 * 0.5 * 0.5 * 0.25
    expected_u2 = 0.25 ** 4

    u12_score = u12_pwm.score_sequence(seq)
    u2_score = u2_pwm.score_sequence(seq)

    # Verify scores are exactly as calculated
    assert abs(u12_score - expected_u12) < 1e-10, \
        f"U12 score {u12_score} doesn't match expected {expected_u12}"
    assert abs(u2_score - expected_u2) < 1e-10, \
        f"U2 score {u2_score} doesn't match expected {expected_u2}"

    # Calculate log ratio
    expected_ratio = expected_u12 / expected_u2  # = 8.0
    expected_log_ratio = math.log2(expected_ratio)  # = 3.0

    log_ratio = math.log2(u12_score / u2_score)
    assert abs(log_ratio - expected_log_ratio) < 1e-10, \
        f"Log ratio {log_ratio} doesn't match expected {expected_log_ratio}"

    # Explicitly verify the expected value
    assert abs(log_ratio - 3.0) < 1e-10, \
        f"Log ratio should be exactly 3.0, got {log_ratio}"


@pytest.mark.parametrize("factor,expected_log", [
    (2.0, 1.0),    # 2x difference → log2(2) = 1
    (4.0, 2.0),    # 4x difference → log2(4) = 2
    (8.0, 3.0),    # 8x difference → log2(8) = 3
    (16.0, 4.0),   # 16x difference → log2(16) = 4
    (0.5, -1.0),   # 1/2x difference → log2(0.5) = -1
    (0.25, -2.0),  # 1/4x difference → log2(0.25) = -2
])
def test_log_ratio_parametrized(factor, expected_log):
    """
    Parametrized test for various factor differences.

    Tests that log2(factor) = expected_log for a range of factors.
    """
    # Create simple PWMs with the specified factor difference
    u12_val = 0.5 * factor if factor <= 1 else 0.5
    u2_val = 0.5 / factor if factor > 1 else 0.5

    # Normalize to ensure they sum correctly (rough approximation)
    u12_other = (1.0 - u12_val) / 3
    u2_other = (1.0 - u2_val) / 3

    u12_matrix = np.array([[u12_val] * 4, [u12_other] * 4, [u12_other] * 4, [u12_other] * 4])
    u2_matrix = np.array([[u2_val] * 4, [u2_other] * 4, [u2_other] * 4, [u2_other] * 4])

    u12_pwm = PWM("u12", u12_matrix, length=4)
    u2_pwm = PWM("u2", u2_matrix, length=4)

    seq = "AAAA"
    u12_score = u12_pwm.score_sequence(seq)
    u2_score = u2_pwm.score_sequence(seq)

    log_ratio = math.log2(u12_score / u2_score)

    # Each position contributes the factor, raised to 4th power
    expected_total_log = expected_log * 4
    assert abs(log_ratio - expected_total_log) < 0.1, \
        f"For factor {factor}, expected log ratio {expected_total_log}, got {log_ratio}"
