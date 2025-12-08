"""
Tests for PWM format equivalence between legacy .iic and JSON formats.

This module verifies that:
1. Both formats load the same matrix data
2. Scoring is identical for the same sequences across formats
3. All 12 matrices produce equivalent scores

Test Strategy:
- Load matrices from both legacy .iic and JSON formats
- Score identical sequences with both versions
- Assert scores are identical (within floating-point tolerance)
- Test all matrix types (five, bp, three) and intron types (u12, u2)
"""

from pathlib import Path
from typing import Dict

import numpy as np
import pytest

from intronIC.scoring.pwm import PWMLoader, PWMSet

# ============================================================================
# Path Configuration
# ============================================================================

# Get the actual data directory from the package
DATA_DIR = Path(__file__).parent.parent.parent.parent / "src" / "intronIC" / "data"
LEGACY_PWM_FILE = DATA_DIR / "archive" / "scoring_matrices.fasta.iic"
JSON_PWM_FILE = DATA_DIR / "intronIC_scoring_PWMs.json"


# ============================================================================
# Test Fixtures
# ============================================================================


@pytest.fixture
def legacy_pwm_file() -> Path:
    """Path to legacy .iic format PWM file."""
    return LEGACY_PWM_FILE


@pytest.fixture
def json_pwm_file() -> Path:
    """Path to JSON format PWM file."""
    return JSON_PWM_FILE


@pytest.fixture
def legacy_pwms(legacy_pwm_file) -> Dict[str, PWMSet]:
    """Load PWMs from legacy format."""
    if not legacy_pwm_file.exists():
        pytest.skip(f"Legacy PWM file not found: {legacy_pwm_file}")
    return PWMLoader.load_from_file(legacy_pwm_file)


@pytest.fixture
def json_pwms(json_pwm_file) -> Dict[str, PWMSet]:
    """Load PWMs from JSON format."""
    if not json_pwm_file.exists():
        pytest.skip(f"JSON PWM file not found: {json_pwm_file}")
    return PWMLoader.load_from_file(json_pwm_file)


@pytest.fixture
def test_sequences() -> Dict[str, str]:
    """
    Test sequences for scoring.

    These are realistic sequences that should score differently
    on U12 vs U2 matrices.
    """
    return {
        # U12-type AT-AC five prime (strong AT at positions 0-1)
        "u12_atac_five": "CAGTATCCTTC",  # positions -3 to +7
        # U2-type GT-AG five prime (strong GT at positions 0-1)
        "u2_gtag_five": "AAGGTAAGTAT",  # positions -3 to +7
        # U12 branch point (strong TCCTTAAC motif)
        "u12_bp": "ACTCCTTAACCGTA",  # typical U12 BP region
        # U2 branch point (weaker consensus)
        "u2_bp": "ACTGACTAACCGTA",  # typical U2 BP region
        # U12 three prime (strong YAG)
        "u12_three": "TTTTTTTTTTTCAG",  # poly-T tract + YAG
        # U2 three prime
        "u2_three": "CTCTCTCTCTCTAG",  # typical U2 3' region
    }


# ============================================================================
# Format Detection Tests
# ============================================================================


def test_format_detection_legacy(legacy_pwm_file):
    """Test that legacy .iic format is detected and loaded correctly."""
    if not legacy_pwm_file.exists():
        pytest.skip(f"Legacy PWM file not found: {legacy_pwm_file}")
    pwms = PWMLoader.load_from_file(legacy_pwm_file)

    # Should load all three regions
    assert "five" in pwms
    assert "bp" in pwms
    assert "three" in pwms


def test_format_detection_json(json_pwm_file):
    """Test that JSON format is detected and loaded correctly."""
    if not json_pwm_file.exists():
        pytest.skip(f"JSON PWM file not found: {json_pwm_file}")
    pwms = PWMLoader.load_from_file(json_pwm_file)

    # Should load all three regions
    assert "five" in pwms
    assert "bp" in pwms
    assert "three" in pwms


# ============================================================================
# Matrix Structure Equivalence Tests
# ============================================================================


def test_same_regions_loaded(legacy_pwms, json_pwms):
    """Test that both formats load the same regions."""
    assert set(legacy_pwms.keys()) == set(json_pwms.keys())


def test_same_matrix_names(legacy_pwms, json_pwms):
    """
    Test that both formats load the same matrix names within each region.

    Note: The JSON format contains additional U2 branch point matrices that
    are not present in the legacy .iic format. For five and three regions,
    the matrices should be identical. For bp, we verify the legacy is a
    subset of JSON.
    """
    # five and three regions should have identical keys
    for region in ["five", "three"]:
        legacy_set = legacy_pwms[region]
        json_set = json_pwms[region]

        assert set(legacy_set.matrices.keys()) == set(json_set.matrices.keys()), (
            f"Matrix keys differ for region {region}"
        )

    # bp region: legacy is a subset of JSON (JSON has u2-gtag-bp)
    legacy_bp_keys = set(legacy_pwms["bp"].matrices.keys())
    json_bp_keys = set(json_pwms["bp"].matrices.keys())

    assert legacy_bp_keys <= json_bp_keys, (
        f"Legacy bp matrices should be subset of JSON bp matrices.\n"
        f"Missing in JSON: {legacy_bp_keys - json_bp_keys}"
    )


def test_same_matrix_dimensions(legacy_pwms, json_pwms):
    """Test that matrices have identical dimensions across formats."""
    for region in ["five", "bp", "three"]:
        legacy_set = legacy_pwms[region]
        json_set = json_pwms[region]

        # Compare all matrices
        for key in legacy_set.matrices.keys():
            legacy_pwm = legacy_set.matrices[key]
            json_pwm = json_set.matrices[key]

            assert legacy_pwm.length == json_pwm.length, (
                f"{region} {key}: Legacy length {legacy_pwm.length} != JSON length {json_pwm.length}"
            )


def test_same_start_indices(legacy_pwms, json_pwms):
    """Test that matrices have identical start indices across formats."""
    for region in ["five", "bp", "three"]:
        legacy_set = legacy_pwms[region]
        json_set = json_pwms[region]

        for key in legacy_set.matrices.keys():
            legacy_pwm = legacy_set.matrices[key]
            json_pwm = json_set.matrices[key]

            assert legacy_pwm.start_index == json_pwm.start_index, (
                f"{region} {key}: Legacy start_index {legacy_pwm.start_index} != JSON start_index {json_pwm.start_index}"
            )


def test_same_matrix_values(legacy_pwms, json_pwms):
    """Test that matrix values are numerically equivalent."""
    for region in ["five", "bp", "three"]:
        legacy_set = legacy_pwms[region]
        json_set = json_pwms[region]

        for key in legacy_set.matrices.keys():
            legacy_pwm = legacy_set.matrices[key]
            json_pwm = json_set.matrices[key]

            # Compare actual matrix values
            np.testing.assert_allclose(
                legacy_pwm.matrix,
                json_pwm.matrix,
                rtol=1e-7,
                atol=1e-10,
                err_msg=f"{region} {key}: Matrix values differ",
            )


# ============================================================================
# Scoring Equivalence Tests
# ============================================================================


def test_five_prime_scoring_equivalence(legacy_pwms, json_pwms, test_sequences):
    """Test that 5' splice site scoring is equivalent across formats."""
    legacy_five = legacy_pwms["five"]
    json_five = json_pwms["five"]

    # Test each 5' matrix
    for key in legacy_five.matrices.keys():
        legacy_pwm = legacy_five.matrices[key]
        json_pwm = json_five.matrices[key]

        # Score the same sequence with both
        for seq_name, seq in test_sequences.items():
            if len(seq) >= legacy_pwm.length:
                test_seq = seq[: legacy_pwm.length]
                legacy_score = legacy_pwm.score_sequence(
                    test_seq, ignore_positions=set()
                )
                json_score = json_pwm.score_sequence(test_seq, ignore_positions=set())

                assert abs(legacy_score - json_score) < 1e-10, (
                    f"Score mismatch for {key} on {seq_name}: legacy={legacy_score}, json={json_score}"
                )


def test_branch_point_scoring_equivalence(legacy_pwms, json_pwms, test_sequences):
    """Test that branch point scoring is equivalent across formats."""
    legacy_bp = legacy_pwms["bp"]
    json_bp = json_pwms["bp"]

    for key in legacy_bp.matrices.keys():
        legacy_pwm = legacy_bp.matrices[key]
        json_pwm = json_bp.matrices[key]

        # Score BP-specific sequences
        for seq_name in ["u12_bp", "u2_bp"]:
            seq = test_sequences[seq_name]
            if len(seq) >= legacy_pwm.length:
                test_seq = seq[: legacy_pwm.length]
                legacy_score = legacy_pwm.score_sequence(
                    test_seq, ignore_positions=set()
                )
                json_score = json_pwm.score_sequence(test_seq, ignore_positions=set())

                assert abs(legacy_score - json_score) < 1e-10, (
                    f"Score mismatch for {key} on {seq_name}: legacy={legacy_score}, json={json_score}"
                )


def test_three_prime_scoring_equivalence(legacy_pwms, json_pwms, test_sequences):
    """Test that 3' splice site scoring is equivalent across formats."""
    legacy_three = legacy_pwms["three"]
    json_three = json_pwms["three"]

    for key in legacy_three.matrices.keys():
        legacy_pwm = legacy_three.matrices[key]
        json_pwm = json_three.matrices[key]

        # Score 3' sequences
        for seq_name in ["u12_three", "u2_three"]:
            seq = test_sequences[seq_name]
            if len(seq) >= legacy_pwm.length:
                test_seq = seq[: legacy_pwm.length]
                legacy_score = legacy_pwm.score_sequence(
                    test_seq, ignore_positions=set()
                )
                json_score = json_pwm.score_sequence(test_seq, ignore_positions=set())

                assert abs(legacy_score - json_score) < 1e-10, (
                    f"Score mismatch for {key} on {seq_name}: legacy={legacy_score}, json={json_score}"
                )


# ============================================================================
# Comprehensive Matrix Coverage Tests
# ============================================================================


def test_all_matrices_present(legacy_pwms, json_pwms):
    """
    Verify all expected matrices are present in both formats.

    Matrix keys are tuples:
    - five/three: (intron_class, splice_type) e.g., ('u12', 'atac'), ('u2', 'gtag')
    - bp: (intron_class, splice_type, adenosine_pos) e.g., ('u12', 'atac', 'A10')

    Note: JSON format is a superset of legacy format - it includes U2 branch
    point matrices that don't exist in the legacy .iic file.
    """
    # Check that both formats have the same regions
    assert set(legacy_pwms.keys()) == set(json_pwms.keys()), (
        f"Region mismatch: legacy={set(legacy_pwms.keys())}, json={set(json_pwms.keys())}"
    )

    # For five and three regions, verify identical matrix keys
    for region in ["five", "three"]:
        legacy_keys = set(legacy_pwms[region].matrices.keys())
        json_keys = set(json_pwms[region].matrices.keys())

        assert legacy_keys == json_keys, (
            f"Matrix key mismatch for {region}:\n  Legacy: {legacy_keys}\n  JSON: {json_keys}"
        )

    # For bp region, legacy should be subset of JSON (JSON has additional u2-gtag-bp)
    legacy_bp_keys = set(legacy_pwms["bp"].matrices.keys())
    json_bp_keys = set(json_pwms["bp"].matrices.keys())

    assert legacy_bp_keys <= json_bp_keys, (
        f"Legacy bp matrices should be subset of JSON:\n"
        f"  Missing in JSON: {legacy_bp_keys - json_bp_keys}"
    )

    # Verify expected number of matrices per region
    # five and three have 4 matrices each (u12/u2 x atac/gtag + gcag)
    assert len(legacy_pwms["five"].matrices) == 4, (
        f"Expected 4 five' matrices, got {len(legacy_pwms['five'].matrices)}"
    )
    assert len(legacy_pwms["three"].matrices) == 4, (
        f"Expected 4 three' matrices, got {len(legacy_pwms['three'].matrices)}"
    )

    # bp has 4 legacy matrices (u12 atac/gtag x A9/A10), JSON has 5 (adds u2-gtag-bp)
    assert len(legacy_pwms["bp"].matrices) == 4, (
        f"Expected 4 legacy bp matrices, got {len(legacy_pwms['bp'].matrices)}"
    )
    assert len(json_pwms["bp"].matrices) == 5, (
        f"Expected 5 JSON bp matrices (legacy + u2-gtag-bp), got {len(json_pwms['bp'].matrices)}"
    )


def test_random_sequence_scoring_equivalence(legacy_pwms, json_pwms):
    """Test scoring equivalence on random DNA sequences."""
    import random

    random.seed(42)  # Reproducible

    bases = "ACGT"

    for region in ["five", "bp", "three"]:
        legacy_set = legacy_pwms[region]
        json_set = json_pwms[region]

        for key in legacy_set.matrices.keys():
            legacy_pwm = legacy_set.matrices[key]
            json_pwm = json_set.matrices[key]

            # Generate random sequences of appropriate length
            for _ in range(10):
                seq = "".join(random.choice(bases) for _ in range(legacy_pwm.length))
                legacy_score = legacy_pwm.score_sequence(seq, ignore_positions=set())
                json_score = json_pwm.score_sequence(seq, ignore_positions=set())

                assert abs(legacy_score - json_score) < 1e-10, (
                    f"Score mismatch for {key} on random seq: legacy={legacy_score}, json={json_score}"
                )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
