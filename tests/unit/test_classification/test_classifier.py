"""
Tests for IntronClassifier - initialization and validation tests.

This module tests only fast initialization and validation.
Tests that perform actual SVM training have been moved to:
tests/integration/test_classification_pipeline.py

Port from: intronIC.py:5038-5900
"""

import pytest

from intronIC.classification.classifier import IntronClassifier, ClassificationResult
from intronIC.core.intron import Intron, IntronScores, IntronSequences, GenomicCoordinate


# Test fixtures for validation tests

@pytest.fixture
def u12_reference():
    """Create minimal reference U12 introns with z-scores for validation tests."""
    introns = []
    for i in range(5):  # Just a few for validation
        intron = Intron(
            intron_id=f"ref_u12_{i}",
            coordinates=GenomicCoordinate(
                chromosome="chr1",
                start=1000 + i * 100,
                stop=1100 + i * 100,
                strand="+",
                system="1-based"
            ),
            sequences=IntronSequences(
                seq="GTATGT" + "N" * 50 + "TCCTTAAC",
                five_seq="GTATGT",
                three_seq="TCCTTAAC",
                bp_seq="TCCTTAAC"
            ),
            scores=IntronScores(
                five_z_score=2.0,
                bp_z_score=2.5,
                three_z_score=2.0,
            )
        )
        introns.append(intron)
    return introns


@pytest.fixture
def u2_reference():
    """Create minimal reference U2 introns with z-scores for validation tests."""
    introns = []
    for i in range(5):  # Just a few for validation
        intron = Intron(
            intron_id=f"ref_u2_{i}",
            coordinates=GenomicCoordinate(
                chromosome="chr1",
                start=10000 + i * 100,
                stop=10100 + i * 100,
                strand="+",
                system="1-based"
            ),
            sequences=IntronSequences(
                seq="GTAAGT" + "N" * 50 + "TTTCAG",
                five_seq="GTAAGT",
                three_seq="TTTCAG",
                bp_seq="CTAAC"
            ),
            scores=IntronScores(
                five_z_score=-1.0,
                bp_z_score=-1.5,
                three_z_score=-1.0,
            )
        )
        introns.append(intron)
    return introns


@pytest.fixture
def experimental_mixed():
    """Create minimal experimental introns for validation tests."""
    introns = []
    for i in range(2):
        intron = Intron(
            intron_id=f"exp_{i}",
            coordinates=GenomicCoordinate(
                chromosome="chr2",
                start=1000 + i * 100,
                stop=1100 + i * 100,
                strand="+",
                system="1-based"
            ),
            sequences=IntronSequences(
                seq="GTATGT" + "N" * 50 + "TCCTTAAC",
                five_seq="GTATGT",
                three_seq="TCCTTAAC"
            ),
            scores=IntronScores(
                five_z_score=2.0,
                bp_z_score=2.5,
                three_z_score=2.0,
            )
        )
        introns.append(intron)
    return introns


# Test IntronClassifier initialization

def test_classifier_initialization():
    """Test IntronClassifier initialization with default parameters."""
    classifier = IntronClassifier()
    assert classifier.n_optimization_rounds == 3
    assert classifier.n_ensemble_models == 3
    assert classifier.classification_threshold == 90.0
    assert classifier.subsample_u2 is True
    assert classifier.optimize_c is True


def test_classifier_custom_parameters():
    """Test IntronClassifier with custom parameters."""
    classifier = IntronClassifier(
        n_optimization_rounds=3,
        n_ensemble_models=5,
        classification_threshold=85.0,
        subsample_u2=False,
        random_state=123
    )
    assert classifier.n_optimization_rounds == 3
    assert classifier.n_ensemble_models == 5
    assert classifier.classification_threshold == 85.0
    assert classifier.subsample_u2 is False
    assert classifier.random_state == 123


def test_classifier_invalid_threshold():
    """Test IntronClassifier rejects invalid threshold."""
    with pytest.raises(ValueError, match="must be 0-100"):
        IntronClassifier(classification_threshold=-10)

    with pytest.raises(ValueError, match="must be 0-100"):
        IntronClassifier(classification_threshold=150)


def test_classifier_fixed_c_without_value():
    """Test that fixed_c is required when optimize_c=False."""
    with pytest.raises(ValueError, match="Must provide fixed_c"):
        IntronClassifier(optimize_c=False, fixed_c=None)


def test_classifier_fixed_c_with_value():
    """Test IntronClassifier with fixed C parameter."""
    classifier = IntronClassifier(optimize_c=False, fixed_c=1.0)
    assert classifier.optimize_c is False
    assert classifier.fixed_c == 1.0


# =============================================================================
# NOTE: Tests that perform actual SVM training have been moved to:
#       tests/integration/test_classification_pipeline.py
#
# This unit test file now contains only fast initialization and validation tests.
# =============================================================================


# Test validation

def test_classify_validates_reference_u12_z_scores(u2_reference, experimental_mixed):
    """Test that classify validates U12 reference introns have z-scores."""
    # Create U12 reference without z-scores
    bad_u12 = [
        Intron(
            intron_id="bad_u12",
            coordinates=GenomicCoordinate(
                chromosome="chr1",
                start=1000,
                stop=1100,
                strand="+",
                system="1-based"
            ),
            sequences=IntronSequences(seq="ATCG", five_seq="AT", three_seq="CG"),
            scores=IntronScores(
                five_z_score=None,  # Missing!
                bp_z_score=1.0,
                three_z_score=1.0
            )
        )
    ]

    classifier = IntronClassifier()

    with pytest.raises(ValueError, match="u12_reference.*Missing z-scores"):
        classifier.classify(bad_u12, u2_reference, experimental_mixed)


def test_classify_validates_reference_u2_z_scores(u12_reference, experimental_mixed):
    """Test that classify validates U2 reference introns have z-scores."""
    # Create U2 reference without z-scores
    bad_u2 = [
        Intron(
            intron_id="bad_u2",
            coordinates=GenomicCoordinate(
                chromosome="chr1",
                start=1000,
                stop=1100,
                strand="+",
                system="1-based"
            ),
            sequences=IntronSequences(seq="ATCG", five_seq="AT", three_seq="CG"),
            scores=None  # No scores at all!
        )
    ]

    classifier = IntronClassifier()

    with pytest.raises(ValueError, match="u2_reference.*No scores object"):
        classifier.classify(u12_reference, bad_u2, experimental_mixed)


def test_classify_validates_experimental_z_scores(u12_reference, u2_reference):
    """Test that classify validates experimental introns have z-scores."""
    # Create experimental introns without z-scores
    bad_exp = [
        Intron(
            intron_id="bad_exp",
            coordinates=GenomicCoordinate(
                chromosome="chr1",
                start=1000,
                stop=1100,
                strand="+",
                system="1-based"
            ),
            sequences=IntronSequences(seq="ATCG", five_seq="AT", three_seq="CG"),
            scores=IntronScores(
                five_z_score=1.0,
                bp_z_score=None,  # Missing!
                three_z_score=1.0
            )
        )
    ]

    classifier = IntronClassifier()

    with pytest.raises(ValueError, match="experimental.*Missing z-scores"):
        classifier.classify(u12_reference, u2_reference, bad_exp)


# Test edge cases
# Note: All tests that perform actual training have been moved to integration tests
