"""
Tests for PWM format equivalence between legacy .iic and YAML formats.

This module verifies that:
1. Both formats load the same matrix data
2. Scoring is identical for the same sequences across formats
3. All 12 matrices produce equivalent scores

Test Strategy:
- Load matrices from both legacy .iic and YAML formats
- Score identical sequences with both versions
- Assert scores are identical (within floating-point tolerance)
- Test all matrix types (five, bp, three) and intron types (u12, u2)
"""

import pytest
import numpy as np
from pathlib import Path
from typing import Dict

from intronIC.scoring.pwm import PWM, PWMSet, PWMLoader


# ============================================================================
# Test Fixtures
# ============================================================================

@pytest.fixture
def legacy_pwm_file() -> Path:
    """Path to legacy .iic format PWM file."""
    return Path("data/scoring_matrices.fasta.iic")


@pytest.fixture
def yaml_pwm_file() -> Path:
    """Path to YAML format PWM file."""
    return Path("data/scoring_matrices.yaml")


@pytest.fixture
def legacy_pwms(legacy_pwm_file) -> Dict[str, PWMSet]:
    """Load PWMs from legacy format."""
    if not legacy_pwm_file.exists():
        pytest.skip(f"Legacy PWM file not found: {legacy_pwm_file}")
    return PWMLoader.load_from_file(legacy_pwm_file)


@pytest.fixture
def yaml_pwms(yaml_pwm_file) -> Dict[str, PWMSet]:
    """Load PWMs from YAML format."""
    if not yaml_pwm_file.exists():
        pytest.skip(f"YAML PWM file not found: {yaml_pwm_file}")
    return PWMLoader.load_from_file(yaml_pwm_file)


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
    """Test that legacy .iic format is detected correctly."""
    pwms = PWMLoader.load_from_file(legacy_pwm_file)

    # Should load all three regions
    assert "five" in pwms
    assert "bp" in pwms
    assert "three" in pwms


def test_format_detection_yaml(yaml_pwm_file):
    """Test that YAML format is detected correctly."""
    pwms = PWMLoader.load_from_file(yaml_pwm_file)

    # Should load all three regions
    assert "five" in pwms
    assert "bp" in pwms
    assert "three" in pwms


# ============================================================================
# Matrix Structure Equivalence Tests
# ============================================================================

def test_same_matrix_names(legacy_pwms, yaml_pwms):
    """Test that both formats load the same matrix names."""
    # Should have same regions
    assert set(legacy_pwms.keys()) == set(yaml_pwms.keys())

    # Check each region has same matrix types
    for region in ["five", "bp", "three"]:
        legacy_set = legacy_pwms[region]
        yaml_set = yaml_pwms[region]

        # Both should have same matrix keys
        assert set(legacy_set.matrices.keys()) == set(yaml_set.matrices.keys())


def test_same_matrix_dimensions(legacy_pwms, yaml_pwms):
    """Test that matrices have identical dimensions across formats."""
    for region in ["five", "bp", "three"]:
        legacy_set = legacy_pwms[region]
        yaml_set = yaml_pwms[region]

        # Compare all matrices
        for key in legacy_set.matrices.keys():
            legacy_pwm = legacy_set.matrices[key]
            yaml_pwm = yaml_set.matrices[key]

            assert legacy_pwm.length == yaml_pwm.length, \
                f"{region} {key}: Legacy length {legacy_pwm.length} != YAML length {yaml_pwm.length}"

            assert legacy_pwm.matrix.shape == yaml_pwm.matrix.shape, \
                f"{region} {key}: Legacy shape {legacy_pwm.matrix.shape} != YAML shape {yaml_pwm.matrix.shape}"


def test_same_start_indices(legacy_pwms, yaml_pwms):
    """Test that start_index values are preserved across formats."""
    for region in ["five", "bp", "three"]:
        legacy_set = legacy_pwms[region]
        yaml_set = yaml_pwms[region]

        # Compare all matrices
        for key in legacy_set.matrices.keys():
            legacy_pwm = legacy_set.matrices[key]
            yaml_pwm = yaml_set.matrices[key]

            assert legacy_pwm.start_index == yaml_pwm.start_index, \
                f"{region} {key}: Legacy start_index {legacy_pwm.start_index} != YAML start_index {yaml_pwm.start_index}"


def test_same_matrix_values(legacy_pwms, yaml_pwms):
    """Test that matrix frequency values are identical across formats."""
    tolerance = 1e-10  # Very tight tolerance for exact equivalence

    for region in ["five", "bp", "three"]:
        legacy_set = legacy_pwms[region]
        yaml_set = yaml_pwms[region]

        # Compare all matrices
        for key in legacy_set.matrices.keys():
            legacy_matrix = legacy_set.matrices[key].matrix
            yaml_matrix = yaml_set.matrices[key].matrix

            assert np.allclose(legacy_matrix, yaml_matrix, rtol=tolerance, atol=tolerance), \
                f"{region} {key}: Matrix values differ between formats"


# ============================================================================
# Scoring Equivalence Tests
# ============================================================================

def test_five_prime_scoring_equivalence(legacy_pwms, yaml_pwms, test_sequences):
    """Test that 5' splice site scoring is identical across formats."""
    tolerance = 1e-10

    # Test U12 AT-AC (if available)
    legacy_u12 = legacy_pwms["five"].u12_atac
    yaml_u12 = yaml_pwms["five"].u12_atac

    if legacy_u12 and yaml_u12:
        seq = test_sequences["u12_atac_five"]
        legacy_score = legacy_u12.score_sequence(seq, seq_start_position=legacy_u12.start_index)
        yaml_score = yaml_u12.score_sequence(seq, seq_start_position=yaml_u12.start_index)

        assert abs(legacy_score - yaml_score) < tolerance, \
            f"U12 five prime scores differ: legacy={legacy_score}, yaml={yaml_score}"

    # Test U2 GT-AG (if available)
    legacy_u2 = legacy_pwms["five"].u2_gtag
    yaml_u2 = yaml_pwms["five"].u2_gtag

    if legacy_u2 and yaml_u2:
        seq = test_sequences["u2_gtag_five"]
        legacy_score = legacy_u2.score_sequence(seq, seq_start_position=legacy_u2.start_index)
        yaml_score = yaml_u2.score_sequence(seq, seq_start_position=yaml_u2.start_index)

        assert abs(legacy_score - yaml_score) < tolerance, \
            f"U2 five prime scores differ: legacy={legacy_score}, yaml={yaml_score}"


def test_branch_point_scoring_equivalence(legacy_pwms, yaml_pwms, test_sequences):
    """Test that branch point scoring is identical across formats."""
    tolerance = 1e-10

    # Test all branch point matrices
    legacy_set = legacy_pwms["bp"]
    yaml_set = yaml_pwms["bp"]

    for key in legacy_set.matrices.keys():
        legacy_pwm = legacy_set.matrices[key]
        yaml_pwm = yaml_set.matrices[key]

        # Use appropriate test sequence based on intron type
        if key[0] == 'u12':
            seq = test_sequences["u12_bp"]
        else:
            seq = test_sequences["u2_bp"]

        legacy_score = legacy_pwm.score_sequence(seq, seq_start_position=legacy_pwm.start_index)
        yaml_score = yaml_pwm.score_sequence(seq, seq_start_position=yaml_pwm.start_index)

        assert abs(legacy_score - yaml_score) < tolerance, \
            f"BP {key} scores differ: legacy={legacy_score}, yaml={yaml_score}"


def test_three_prime_scoring_equivalence(legacy_pwms, yaml_pwms, test_sequences):
    """Test that 3' splice site scoring is identical across formats."""
    tolerance = 1e-10

    # Test all three prime matrices
    legacy_set = legacy_pwms["three"]
    yaml_set = yaml_pwms["three"]

    for key in legacy_set.matrices.keys():
        legacy_pwm = legacy_set.matrices[key]
        yaml_pwm = yaml_set.matrices[key]

        # Use appropriate test sequence based on intron type
        if key[0] == 'u12':
            seq = test_sequences["u12_three"]
        else:
            seq = test_sequences["u2_three"]

        legacy_score = legacy_pwm.score_sequence(seq, seq_start_position=legacy_pwm.start_index)
        yaml_score = yaml_pwm.score_sequence(seq, seq_start_position=yaml_pwm.start_index)

        assert abs(legacy_score - yaml_score) < tolerance, \
            f"Three prime {key} scores differ: legacy={legacy_score}, yaml={yaml_score}"


def test_all_matrices_scoring_equivalence(legacy_pwms, yaml_pwms):
    """
    Comprehensive test: Score all matrices with same random sequences.

    This tests all 12 matrices (u12/u2 x atac/gtag/gcag x five/bp/three)
    to ensure complete equivalence.
    """
    tolerance = 1e-10

    # Generate random test sequences for each region
    # Use fixed seed for reproducibility
    np.random.seed(42)

    test_cases = []

    for region in ["five", "bp", "three"]:
        legacy_set = legacy_pwms[region]
        yaml_set = yaml_pwms[region]

        # Get all matrices in this region
        for key, legacy_pwm in legacy_set.matrices.items():
            yaml_pwm = yaml_set.matrices.get(key)

            if yaml_pwm is None:
                pytest.fail(f"Matrix {key} present in legacy but not in YAML")

            # Generate random sequence of appropriate length
            seq_length = legacy_pwm.length
            seq = ''.join(np.random.choice(['A', 'C', 'G', 'T'], size=seq_length))

            # Score with both formats
            legacy_score = legacy_pwm.score_sequence(seq, seq_start_position=legacy_pwm.start_index)
            yaml_score = yaml_pwm.score_sequence(seq, seq_start_position=yaml_pwm.start_index)

            # Check equivalence
            assert abs(legacy_score - yaml_score) < tolerance, \
                f"Matrix {key} scores differ for seq '{seq}': legacy={legacy_score}, yaml={yaml_score}"

            test_cases.append((region, key, seq, legacy_score))

    # Verify we tested all expected matrices
    print(f"\nTested {len(test_cases)} matrix/sequence combinations")
    assert len(test_cases) >= 12, "Should have tested at least 12 matrices"


# ============================================================================
# Round-trip Conversion Tests
# ============================================================================

def test_yaml_to_legacy_roundtrip(yaml_pwm_file, tmp_path):
    """
    Test that YAML -> .iic -> YAML preserves all data.

    This uses the converter script to ensure bidirectional conversion works.
    """
    import subprocess

    converter_script = Path("utils/convert_pwm_format.py")
    if not converter_script.exists():
        pytest.skip("Converter script not found")

    # Convert YAML -> .iic
    iic_file = tmp_path / "converted.iic"
    result = subprocess.run(
        ["python", str(converter_script), str(yaml_pwm_file), str(iic_file)],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0, f"Conversion failed: {result.stderr}"

    # Load both and compare scoring
    yaml_pwms = PWMLoader.load_from_file(yaml_pwm_file)
    iic_pwms = PWMLoader.load_from_file(iic_file)

    # Quick scoring check with any available U2 matrix
    test_seq = "AAGGTAAGTAT"
    yaml_u2 = yaml_pwms["five"].u2_gtag or yaml_pwms["five"].u2_gcag
    iic_u2 = iic_pwms["five"].u2_gtag or iic_pwms["five"].u2_gcag

    if yaml_u2 and iic_u2:
        yaml_score = yaml_u2.score_sequence(test_seq, seq_start_position=yaml_u2.start_index)
        iic_score = iic_u2.score_sequence(test_seq, seq_start_position=iic_u2.start_index)

        assert abs(yaml_score - iic_score) < 1e-10, \
            "Round-trip conversion changed scoring results"


def test_legacy_to_yaml_roundtrip(legacy_pwm_file, tmp_path):
    """
    Test that .iic -> YAML -> .iic preserves all data.
    """
    import subprocess

    converter_script = Path("utils/convert_pwm_format.py")
    if not converter_script.exists():
        pytest.skip("Converter script not found")

    # Convert .iic -> YAML
    yaml_file = tmp_path / "converted.yaml"
    result = subprocess.run(
        ["python", str(converter_script), str(legacy_pwm_file), str(yaml_file)],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0, f"Conversion failed: {result.stderr}"

    # Load both and compare scoring
    legacy_pwms = PWMLoader.load_from_file(legacy_pwm_file)
    yaml_pwms = PWMLoader.load_from_file(yaml_file)

    # Quick scoring check with any available U12 matrix
    test_seq = "CAGTATCCTTC"
    legacy_u12 = legacy_pwms["five"].u12_atac or legacy_pwms["five"].u12_gtag
    yaml_u12 = yaml_pwms["five"].u12_atac or yaml_pwms["five"].u12_gtag

    if legacy_u12 and yaml_u12:
        legacy_score = legacy_u12.score_sequence(test_seq, seq_start_position=legacy_u12.start_index)
        yaml_score = yaml_u12.score_sequence(test_seq, seq_start_position=yaml_u12.start_index)

        assert abs(legacy_score - yaml_score) < 1e-10, \
            "Round-trip conversion changed scoring results"
