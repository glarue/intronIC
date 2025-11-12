"""
Tests for IntronClassifier - high-level classification pipeline orchestrator.

Tests the complete U2/U12 classification pipeline including optimization,
training, and prediction.

Port from: intronIC.py:5038-5900
"""

import pytest
import numpy as np
from dataclasses import replace

from classification.classifier import IntronClassifier, ClassificationResult
from core.intron import Intron, IntronScores, IntronSequences, GenomicCoordinate, IntronMetadata


# Test fixtures

@pytest.fixture
def u12_reference():
    """Create reference U12 introns with z-scores."""
    introns = []
    for i in range(50):  # Larger reference set for training
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
                five_raw_score=12.5,
                bp_raw_score=10.2,
                three_raw_score=15.3,
                five_z_score=2.0 + np.random.randn() * 0.3,
                bp_z_score=2.5 + np.random.randn() * 0.3,
                three_z_score=2.0 + np.random.randn() * 0.3,
            )
        )
        introns.append(intron)
    return introns


@pytest.fixture
def u2_reference():
    """Create reference U2 introns with z-scores."""
    introns = []
    for i in range(50):  # Larger reference set for training
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
                five_raw_score=5.2,
                bp_raw_score=3.8,
                three_raw_score=6.1,
                five_z_score=-1.0 + np.random.randn() * 0.3,
                bp_z_score=-1.5 + np.random.randn() * 0.3,
                three_z_score=-1.0 + np.random.randn() * 0.3,
            )
        )
        introns.append(intron)
    return introns


@pytest.fixture
def experimental_mixed():
    """Create experimental introns with mixed U12/U2-like features."""
    introns = []

    # 5 U12-like introns
    for i in range(5):
        intron = Intron(
            intron_id=f"exp_u12_like_{i}",
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
                five_z_score=2.2 + np.random.randn() * 0.2,
                bp_z_score=2.7 + np.random.randn() * 0.2,
                three_z_score=2.1 + np.random.randn() * 0.2,
            )
        )
        introns.append(intron)

    # 15 U2-like introns
    for i in range(15):
        intron = Intron(
            intron_id=f"exp_u2_like_{i}",
            coordinates=GenomicCoordinate(
                chromosome="chr2",
                start=2000 + i * 100,
                stop=2100 + i * 100,
                strand="+",
                system="1-based"
            ),
            sequences=IntronSequences(
                seq="GTAAGT" + "N" * 50 + "TTTCAG",
                five_seq="GTAAGT",
                three_seq="TTTCAG"
            ),
            scores=IntronScores(
                five_z_score=-0.9 + np.random.randn() * 0.2,
                bp_z_score=-1.3 + np.random.randn() * 0.2,
                three_z_score=-0.8 + np.random.randn() * 0.2,
            )
        )
        introns.append(intron)

    return introns


# Test IntronClassifier initialization

def test_classifier_initialization():
    """Test IntronClassifier initialization with default parameters."""
    classifier = IntronClassifier()
    assert classifier.n_optimization_rounds == 5
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


# Test complete classification pipeline

def test_classify_complete_pipeline(u12_reference, u2_reference, experimental_mixed):
    """Test complete classification pipeline."""
    classifier = IntronClassifier(
        n_optimization_rounds=2,  # Faster for testing
        n_ensemble_models=2,
        classification_threshold=50.0,  # Lower for easier testing
        random_state=42
    )

    result = classifier.classify(
        u12_reference=u12_reference,
        u2_reference=u2_reference,
        experimental=experimental_mixed
    )

    # Check result structure
    assert isinstance(result, ClassificationResult)
    assert len(result.classified_introns) == len(experimental_mixed)
    assert result.ensemble is not None
    assert len(result.ensemble.models) == 2
    assert result.parameters is not None
    assert result.n_u12_reference == len(u12_reference)
    assert result.n_u2_reference == len(u2_reference)

    # Check that all introns have been classified
    for intron in result.classified_introns:
        assert intron.scores is not None
        assert intron.scores.svm_score is not None
        assert 0 <= intron.scores.svm_score <= 100
        assert intron.metadata is not None
        assert intron.metadata.type_id in ['u2', 'u12']


def test_classify_with_fixed_c(u12_reference, u2_reference, experimental_mixed):
    """Test classification with fixed C parameter (no optimization)."""
    classifier = IntronClassifier(
        optimize_c=False,
        fixed_c=1.0,
        n_ensemble_models=2,
        classification_threshold=50.0,
        random_state=42
    )

    result = classifier.classify(
        u12_reference=u12_reference,
        u2_reference=u2_reference,
        experimental=experimental_mixed
    )

    assert result.parameters.C == 1.0
    assert result.parameters.round_found == 0  # Fixed, not optimized
    assert len(result.classified_introns) == len(experimental_mixed)


def test_classify_assigns_u12_and_u2(u12_reference, u2_reference, experimental_mixed):
    """Test that classification assigns both U12 and U2 types."""
    classifier = IntronClassifier(
        n_optimization_rounds=2,
        n_ensemble_models=2,
        classification_threshold=50.0,
        random_state=42
    )

    result = classifier.classify(
        u12_reference=u12_reference,
        u2_reference=u2_reference,
        experimental=experimental_mixed
    )

    # Count classifications
    u12_count = sum(
        1 for i in result.classified_introns
        if i.metadata and i.metadata.type_id == 'u12'
    )
    u2_count = sum(
        1 for i in result.classified_introns
        if i.metadata and i.metadata.type_id == 'u2'
    )

    # Should have both types
    assert u12_count > 0
    assert u2_count > 0
    assert u12_count + u2_count == len(experimental_mixed)


def test_classify_preserves_z_scores(u12_reference, u2_reference, experimental_mixed):
    """
    CRITICAL TEST: Verify z-scores are NOT re-normalized during classification.
    This is Issue #1 fix - prevents data leakage.
    """
    # Store original z-scores
    original_z_scores = {
        intron.intron_id: (
            intron.scores.five_z_score,
            intron.scores.bp_z_score,
            intron.scores.three_z_score
        )
        for intron in experimental_mixed
    }

    classifier = IntronClassifier(
        n_optimization_rounds=2,
        n_ensemble_models=2,
        random_state=42
    )

    result = classifier.classify(
        u12_reference=u12_reference,
        u2_reference=u2_reference,
        experimental=experimental_mixed
    )

    # Check that z-scores are EXACTLY the same
    for intron in result.classified_introns:
        original = original_z_scores[intron.intron_id]
        current = (
            intron.scores.five_z_score,
            intron.scores.bp_z_score,
            intron.scores.three_z_score
        )
        assert original == current, f"Z-scores changed for {intron.intron_id}!"


# Test batch classification

def test_classify_batch(u12_reference, u2_reference, experimental_mixed):
    """Test batch classification produces same results as regular."""
    classifier = IntronClassifier(
        n_optimization_rounds=2,
        n_ensemble_models=2,
        classification_threshold=50.0,
        random_state=42
    )

    # Regular classification
    result_regular = classifier.classify(
        u12_reference=u12_reference,
        u2_reference=u2_reference,
        experimental=experimental_mixed
    )

    # Batch classification with small batch size
    result_batch = classifier.classify_batch(
        u12_reference=u12_reference,
        u2_reference=u2_reference,
        experimental=experimental_mixed,
        batch_size=5
    )

    # Results should be identical
    assert len(result_regular.classified_introns) == len(result_batch.classified_introns)

    for reg, batch in zip(result_regular.classified_introns, result_batch.classified_introns):
        assert abs(reg.scores.svm_score - batch.scores.svm_score) < 1e-6
        assert reg.metadata.type_id == batch.metadata.type_id


# Test ClassificationResult methods

def test_classification_result_get_u12_predictions(u12_reference, u2_reference, experimental_mixed):
    """Test ClassificationResult.get_u12_predictions()."""
    classifier = IntronClassifier(
        n_optimization_rounds=2,
        n_ensemble_models=2,
        classification_threshold=50.0,
        random_state=42
    )

    result = classifier.classify(
        u12_reference=u12_reference,
        u2_reference=u2_reference,
        experimental=experimental_mixed
    )

    u12_predictions = result.get_u12_predictions(threshold=50.0)

    # All should be U12 with score >= threshold
    for intron in u12_predictions:
        assert intron.metadata.type_id == 'u12'
        assert intron.scores.svm_score >= 50.0


def test_classification_result_get_u2_predictions(u12_reference, u2_reference, experimental_mixed):
    """Test ClassificationResult.get_u2_predictions()."""
    classifier = IntronClassifier(
        n_optimization_rounds=2,
        n_ensemble_models=2,
        classification_threshold=50.0,
        random_state=42
    )

    result = classifier.classify(
        u12_reference=u12_reference,
        u2_reference=u2_reference,
        experimental=experimental_mixed
    )

    u2_predictions = result.get_u2_predictions(threshold=50.0)

    # All should be U2 with score < threshold
    for intron in u2_predictions:
        assert intron.metadata.type_id == 'u2'
        assert intron.scores.svm_score < 50.0


def test_classification_result_threshold_affects_filtering(u12_reference, u2_reference, experimental_mixed):
    """Test that threshold parameter affects get_u12_predictions filtering."""
    classifier = IntronClassifier(
        n_optimization_rounds=2,
        n_ensemble_models=2,
        classification_threshold=50.0,
        random_state=42
    )

    result = classifier.classify(
        u12_reference=u12_reference,
        u2_reference=u2_reference,
        experimental=experimental_mixed
    )

    # Lower threshold should give more U12 predictions
    u12_low = result.get_u12_predictions(threshold=40.0)
    u12_high = result.get_u12_predictions(threshold=60.0)

    assert len(u12_low) >= len(u12_high)


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

def test_classify_small_datasets(u12_reference, u2_reference):
    """Test classification with minimal experimental data."""
    # Single experimental intron
    single_exp = [
        Intron(
            intron_id="single",
            coordinates=GenomicCoordinate(
                chromosome="chr1",
                start=1000,
                stop=1100,
                strand="+",
                system="1-based"
            ),
            sequences=IntronSequences(seq="ATCG", five_seq="AT", three_seq="CG"),
            scores=IntronScores(
                five_z_score=2.0,
                bp_z_score=2.5,
                three_z_score=2.0
            )
        )
    ]

    classifier = IntronClassifier(
        n_optimization_rounds=2,
        n_ensemble_models=2,
        random_state=42
    )

    result = classifier.classify(u12_reference, u2_reference, single_exp)

    assert len(result.classified_introns) == 1
    assert result.classified_introns[0].scores.svm_score is not None


def test_classify_reproducibility(u12_reference, u2_reference, experimental_mixed):
    """Test that classification is reproducible with same random_state."""
    classifier1 = IntronClassifier(
        n_optimization_rounds=2,
        n_ensemble_models=2,
        random_state=42
    )
    result1 = classifier1.classify(u12_reference, u2_reference, experimental_mixed)

    classifier2 = IntronClassifier(
        n_optimization_rounds=2,
        n_ensemble_models=2,
        random_state=42
    )
    result2 = classifier2.classify(u12_reference, u2_reference, experimental_mixed)

    # Scores should be identical
    for i1, i2 in zip(result1.classified_introns, result2.classified_introns):
        assert abs(i1.scores.svm_score - i2.scores.svm_score) < 1e-6
        assert i1.metadata.type_id == i2.metadata.type_id
