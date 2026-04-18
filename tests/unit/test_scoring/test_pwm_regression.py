"""
Regression tests for the PWM scoring pipeline.

Tests the full chain: JSON loading -> array indexing -> score_sequence() ->
log2(u12/u2) -> raw score, using fixed test PWMs with hand-computed expected
values. These tests are independent of the installed PWM data.

Separate structural smoke tests verify that the installed PWMs are well-formed
(correct dimensions, position 0 populated, etc.) without asserting specific
score values.

Design:
    Test PWMs (tests/data/test_scoring_pwms.json) use a cycling pattern where
    each position has one dominant base (freq=0.7, cycling A->C->G->T) and
    three weak bases (freq=0.1). U2 PWMs are uniform (0.25). This means:

    - Correct indexing: sequence matching the dominant cycle scores 0.7
    - Off-by-one: the wrong position scores 0.1 instead of 0.7
    - Detectable difference: log2(0.7/0.1) = 2.807 bits per position

    Special positions (splice junction, branch A) use 0.91 instead of 0.7
    to verify these biologically critical positions are correctly handled.

See also:
    docs/pwm_scoring_verification.md — historical verification notes
    scripts/verify_pwm_scoring.py — diagnostic/audit tool for installed PWMs
"""

import math
from pathlib import Path

import pytest

from intronIC.scoring.pwm import PWMLoader
from intronIC.scoring.scorer import IntronScorer


# ============================================================================
# Fixtures
# ============================================================================

TEST_PWM_JSON = Path(__file__).parents[2] / "data" / "test_scoring_pwms.json"


@pytest.fixture(scope="module")
def test_pwm_sets():
    """Load the fixed test PWMs."""
    return PWMLoader.load_from_file(TEST_PWM_JSON)


@pytest.fixture(scope="session")
def installed_pwm_sets(matrix_file):
    """Load the installed (real) PWMs for structural smoke tests."""
    return PWMLoader.load_from_file(matrix_file)


# ============================================================================
# Known-good test data
#
# Pre-computed by hand from the cycling pattern in test_scoring_pwms.json.
# U12 PWMs: dominant base at each position has freq 0.7 (cycling ACGT),
#           except splice junction / branch A positions which use 0.91.
# U2 PWMs:  uniform 0.25 at all positions.
#
# log2 LR = sum_i(log2(u12_freq_i / 0.25))
#         = n_dominant * log2(0.7/0.25) + n_forced * log2(0.91/0.25)
# ============================================================================

# 5'SS GT-AG: 12 positions (-3 to +8), start_index=-3, seq_start_position=-3
# Sequence matches dominant bases at all positions, pos 0 forced to G=0.91
FIVE_GTAG_SEQ = "ACGGACGTACGT"
FIVE_GTAG_LR = 18.2036335493  # 11 * log2(0.7/0.25) + 1 * log2(0.91/0.25)

# 3'SS GT-AG: 10 positions (-6 to +3), start_index=-6, seq_start_position=-6
# Matches dominant bases, pos -2 forced A=0.91, pos -1 forced G=0.91
THREE_GTAG_SEQ = "ACGTAGGTAC"
THREE_GTAG_LR = 15.6112915182  # 8 * log2(0.7/0.25) + 2 * log2(0.91/0.25)

# BPS GT-AG: 12 positions, start_index=0, seq_start_position=0 (default)
# JSON keys are -9 to +2 but scorer uses array indices 0-11
# Position 9 (bio 0, branch A) forced to A=0.91
BP_GTAG_SEQ = "ACGTACGTAAGT"
BP_GTAG_LR = 18.2036335493  # 11 * log2(0.7/0.25) + 1 * log2(0.91/0.25)

# 5'SS with deliberate mismatch at position +1 (array idx 4):
# Dominant base is A, sequence has T → freq 0.1 instead of 0.7
# Difference from matching score: log2(0.7/0.1) = log2(7) ≈ 2.807 bits
FIVE_MISMATCH_SEQ = "ACGGTCGTACGT"
FIVE_MISMATCH_LR = 15.3962786272  # 10 * log2(0.7/0.25) + 1 * log2(0.91/0.25) + 1 * log2(0.1/0.25)


# ============================================================================
# Part A: Scoring Logic Tests (fixed test PWMs)
# ============================================================================


class TestScoringLogic:
    """Tests using fixed test PWMs to verify scoring algorithm correctness."""

    def test_five_prime_gtag_score(self, test_pwm_sets):
        """5'SS GT-AG: full chain from JSON → score matches hand-computed value."""
        u12 = test_pwm_sets["five"].matrices[("u12", "gtag")]
        u2 = test_pwm_sets["five"].matrices[("u2", "gtag")]

        u12_score = u12.score_sequence(FIVE_GTAG_SEQ, seq_start_position=-3)
        u2_score = u2.score_sequence(FIVE_GTAG_SEQ, seq_start_position=-3)
        lr = math.log2(u12_score / u2_score)

        assert abs(lr - FIVE_GTAG_LR) < 0.0001, (
            f"5'SS GT-AG log2 LR = {lr:.6f}, expected {FIVE_GTAG_LR:.6f}"
        )

    def test_three_prime_gtag_score(self, test_pwm_sets):
        """3'SS GT-AG: full chain from JSON → score matches hand-computed value."""
        u12 = test_pwm_sets["three"].matrices[("u12", "gtag")]
        u2 = test_pwm_sets["three"].matrices[("u2", "gtag")]

        u12_score = u12.score_sequence(THREE_GTAG_SEQ, seq_start_position=-6)
        u2_score = u2.score_sequence(THREE_GTAG_SEQ, seq_start_position=-6)
        lr = math.log2(u12_score / u2_score)

        assert abs(lr - THREE_GTAG_LR) < 0.0001, (
            f"3'SS GT-AG log2 LR = {lr:.6f}, expected {THREE_GTAG_LR:.6f}"
        )

    def test_branch_point_gtag_score(self, test_pwm_sets):
        """BPS GT-AG: scored with default seq_start_position=0, matches expected."""
        u12 = test_pwm_sets["bp"].matrices[("u12", "gtag")]
        u2 = test_pwm_sets["bp"].matrices[("u2", "gtag")]

        # BPS scorer always uses seq_start_position=0 (the default)
        u12_score = u12.score_sequence(BP_GTAG_SEQ)
        u2_score = u2.score_sequence(BP_GTAG_SEQ)
        lr = math.log2(u12_score / u2_score)

        assert abs(lr - BP_GTAG_LR) < 0.0001, (
            f"BPS GT-AG log2 LR = {lr:.6f}, expected {BP_GTAG_LR:.6f}"
        )

    def test_mismatch_detectable(self, test_pwm_sets):
        """A single-position mismatch produces a detectably different score."""
        u12 = test_pwm_sets["five"].matrices[("u12", "gtag")]
        u2 = test_pwm_sets["five"].matrices[("u2", "gtag")]

        u12_score = u12.score_sequence(FIVE_MISMATCH_SEQ, seq_start_position=-3)
        u2_score = u2.score_sequence(FIVE_MISMATCH_SEQ, seq_start_position=-3)
        lr = math.log2(u12_score / u2_score)

        assert abs(lr - FIVE_MISMATCH_LR) < 0.0001, (
            f"Mismatch log2 LR = {lr:.6f}, expected {FIVE_MISMATCH_LR:.6f}"
        )
        # The difference from the matching score should be ~2.807 (log2(7))
        assert abs((FIVE_GTAG_LR - lr) - math.log2(7)) < 0.01

    def test_log_base_is_log2(self):
        """_calculate_log_ratio uses log2, not natural log."""
        lr = IntronScorer._calculate_log_ratio(0.8, 0.2)
        expected = math.log2(0.8 / 0.2)  # log2(4) = 2.0
        assert abs(lr - expected) < 1e-10
        # Explicitly check it's NOT natural log
        assert abs(lr - 2.0) < 1e-10
        assert abs(lr - math.log(4.0)) > 0.5  # ln(4) ≈ 1.386, very different

    def test_extraction_window_clipping(self, test_pwm_sets):
        """Extra bases beyond PWM length are silently clipped."""
        u12 = test_pwm_sets["five"].matrices[("u12", "gtag")]

        # Score the standard 12-char sequence
        score_12 = u12.score_sequence(FIVE_GTAG_SEQ, seq_start_position=-3)

        # Score a 13-char sequence (one extra base beyond PWM length)
        score_13 = u12.score_sequence(FIVE_GTAG_SEQ + "A", seq_start_position=-3)

        # Extra base at position +9 is beyond PWM (length=12, positions -3 to +8)
        # so it should be silently skipped
        assert abs(score_12 - score_13) < 1e-15, (
            "Extra base beyond PWM length should not affect score"
        )

    def test_bps_reference_offset(self, test_pwm_sets):
        """BPS reference_offset correctly identifies the branch A position."""
        u12_bp = test_pwm_sets["bp"].matrices[("u12", "gtag")]

        # For a 12-position BPS with JSON keys -9 to +2,
        # position 0 (branch A) is at array index 9
        assert u12_bp.reference_offset == 9

        # Verify the frequency at that position is the forced branch-A value
        from intronIC.scoring.pwm import BASE_TO_INDEX
        a_freq = u12_bp.matrix[BASE_TO_INDEX["A"], u12_bp.reference_offset]
        assert a_freq > 0.90, (
            f"Branch A frequency at reference_offset should be ~0.91, got {a_freq}"
        )

    def test_position_key_format_preserved(self, test_pwm_sets):
        """JSON keys with '+' prefix load correctly to the right array positions."""
        u12 = test_pwm_sets["five"].matrices[("u12", "gtag")]

        # Position +1 (JSON key "+1") should be at array index 4 (start=-3, so +1-(-3)=4)
        # In our cycling pattern, idx 4 has A as dominant (idx 4 % 4 = 0 → A)
        from intronIC.scoring.pwm import BASE_TO_INDEX
        a_freq = u12.matrix[BASE_TO_INDEX["A"], 4]
        assert a_freq == 0.7, (
            f"Position +1 (array idx 4) should have A=0.7, got {a_freq}. "
            "JSON key '+1' may not be loading to the correct array position."
        )

    def test_atac_five_prime_position_zero(self, test_pwm_sets):
        """AT-AC 5'SS position 0 has A as the dominant base (not G)."""
        u12_atac = test_pwm_sets["five"].matrices[("u12", "atac")]
        from intronIC.scoring.pwm import BASE_TO_INDEX

        # Position 0 is at array index 3 (start_index=-3)
        a_freq = u12_atac.matrix[BASE_TO_INDEX["A"], 3]
        g_freq = u12_atac.matrix[BASE_TO_INDEX["G"], 3]
        assert a_freq > 0.90, f"AT-AC 5'SS position 0 should favor A, got A={a_freq}"
        assert g_freq < 0.10, f"AT-AC 5'SS position 0 should not favor G, got G={g_freq}"


# ============================================================================
# Part B: Structural Smoke Tests (installed PWMs)
# ============================================================================


class TestInstalledPWMStructure:
    """Verify the installed PWM file is well-formed (no specific score assertions)."""

    def test_splice_site_dimensions(self, installed_pwm_sets):
        """All splice site PWMs have 40 positions."""
        for region in ["five", "three"]:
            for key, pwm in installed_pwm_sets[region].matrices.items():
                assert pwm.length == 40, (
                    f"{pwm.name}: expected 40 positions, got {pwm.length}"
                )

    def test_bps_dimensions(self, installed_pwm_sets):
        """All BPS PWMs have 12 positions."""
        for key, pwm in installed_pwm_sets["bp"].matrices.items():
            assert pwm.length == 12, (
                f"{pwm.name}: expected 12 positions, got {pwm.length}"
            )

    def test_gtag_five_position_zero_is_g(self, installed_pwm_sets):
        """GT-AG 5'SS position 0 should have G frequency near 1.0."""
        from intronIC.scoring.pwm import BASE_TO_INDEX
        u12 = installed_pwm_sets["five"].matrices[("u12", "gtag")]
        # Position 0 is at array index abs(start_index) = 20
        g_freq = u12.matrix[BASE_TO_INDEX["G"], u12.reference_offset]
        assert g_freq > 0.95, (
            f"GT-AG 5'SS position 0: G frequency should be >0.95, got {g_freq:.4f}"
        )

    def test_atac_five_position_zero_is_a(self, installed_pwm_sets):
        """AT-AC 5'SS position 0 should have A frequency near 1.0."""
        from intronIC.scoring.pwm import BASE_TO_INDEX
        u12 = installed_pwm_sets["five"].matrices[("u12", "atac")]
        a_freq = u12.matrix[BASE_TO_INDEX["A"], u12.reference_offset]
        assert a_freq > 0.95, (
            f"AT-AC 5'SS position 0: A frequency should be >0.95, got {a_freq:.4f}"
        )

    def test_bps_branch_a_is_adenine(self, installed_pwm_sets):
        """BPS position 0 (branch A) should have A frequency near 1.0."""
        from intronIC.scoring.pwm import BASE_TO_INDEX
        for key in [("u12", "gtag"), ("u12", "atac")]:
            if key not in installed_pwm_sets["bp"].matrices:
                continue
            u12 = installed_pwm_sets["bp"].matrices[key]
            a_freq = u12.matrix[BASE_TO_INDEX["A"], u12.reference_offset]
            assert a_freq > 0.95, (
                f"{u12.name} branch A: frequency should be >0.95, got {a_freq:.4f}"
            )

    def test_bps_reference_offset_consistency(self, installed_pwm_sets):
        """BPS reference_offset should be 9 for a 12-position PWM spanning -9 to +2."""
        for key, pwm in installed_pwm_sets["bp"].matrices.items():
            assert pwm.reference_offset == 9, (
                f"{pwm.name}: reference_offset should be 9, got {pwm.reference_offset}"
            )

    def test_splice_site_position_zero_not_uniform(self, installed_pwm_sets):
        """Position 0 in splice site PWMs should have non-uniform frequencies."""
        for region in ["five", "three"]:
            for key, pwm in installed_pwm_sets[region].matrices.items():
                col = pwm.matrix[:, pwm.reference_offset]
                max_freq = col.max()
                min_freq = col.min()
                assert max_freq - min_freq > 0.1, (
                    f"{pwm.name} position 0: frequencies look uniform "
                    f"(max={max_freq:.4f}, min={min_freq:.4f}). "
                    "Position 0 should reflect the splice junction base."
                )
