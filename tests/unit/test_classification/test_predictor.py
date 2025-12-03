"""
Tests for SVMPredictor - ensemble classification.

Tests the classification algorithm that applies trained ensemble models
to classify introns as U2 or U12 type.

Port from: intronIC.py:5651-5900
"""

from dataclasses import replace

import numpy as np
import pytest
from sklearn.calibration import CalibratedClassifierCV
from sklearn.svm import LinearSVC

from intronIC.classification.optimizer import SVMParameters
from intronIC.classification.predictor import SVMPredictor
from intronIC.classification.trainer import SVMEnsemble, SVMModel
from intronIC.core.intron import (
    GenomicCoordinate,
    Intron,
    IntronScores,
    IntronSequences,
)

# Test fixtures


@pytest.fixture
def sample_u12_introns():
    """Create sample U12-type introns with raw scores."""
    introns = []
    for i in range(10):
        intron = Intron(
            intron_id=f"u12_{i}",
            coordinates=GenomicCoordinate(
                chromosome="chr1",
                start=1000 + i * 100,
                stop=1100 + i * 100,
                strand="+",
                system="1-based",
            ),
            sequences=IntronSequences(
                seq="GTATGT" + "N" * 50 + "TCCTTAAC",
                five_seq="GTATGT",
                three_seq="TCCTTAAC",
                bp_seq="TCCTTAAC",
            ),
            scores=IntronScores(
                five_raw_score=2.0 + np.random.randn() * 0.5,
                bp_raw_score=2.5 + np.random.randn() * 0.5,
                three_raw_score=2.0 + np.random.randn() * 0.5,
                five_z_score=2.0 + np.random.randn() * 0.5,
                bp_z_score=2.5 + np.random.randn() * 0.5,
                three_z_score=2.0 + np.random.randn() * 0.5,
            ),
        )
        introns.append(intron)
    return introns


@pytest.fixture
def sample_u2_introns():
    """Create sample U2-type introns with raw scores."""
    introns = []
    for i in range(10):
        intron = Intron(
            intron_id=f"u2_{i}",
            coordinates=GenomicCoordinate(
                chromosome="chr1",
                start=2000 + i * 100,
                stop=2100 + i * 100,
                strand="+",
                system="1-based",
            ),
            sequences=IntronSequences(
                seq="GTAAGT" + "N" * 50 + "TTTCAG",
                five_seq="GTAAGT",
                three_seq="TTTCAG",
                bp_seq="CTAAC",
            ),
            scores=IntronScores(
                five_raw_score=-1.0 + np.random.randn() * 0.5,
                bp_raw_score=-1.5 + np.random.randn() * 0.5,
                three_raw_score=-1.0 + np.random.randn() * 0.5,
                five_z_score=-1.0 + np.random.randn() * 0.5,
                bp_z_score=-1.5 + np.random.randn() * 0.5,
                three_z_score=-1.0 + np.random.randn() * 0.5,
            ),
        )
        introns.append(intron)
    return introns


@pytest.fixture
def trained_ensemble(sample_u12_introns, sample_u2_introns):
    """Create a trained ensemble for testing."""
    # Prepare training data
    u12_features = np.array(
        [
            [i.scores.five_z_score, i.scores.bp_z_score, i.scores.three_z_score]
            for i in sample_u12_introns
        ]
    )
    u2_features = np.array(
        [
            [i.scores.five_z_score, i.scores.bp_z_score, i.scores.three_z_score]
            for i in sample_u2_introns
        ]
    )
    X = np.vstack([u2_features, u12_features])
    y = np.array([0] * len(sample_u2_introns) + [1] * len(sample_u12_introns))

    # Train 3 models with proper LinearSVC + calibration
    models = []
    for i in range(3):
        base_svm = LinearSVC(
            C=1.0,
            class_weight="balanced",
            dual="auto",
            random_state=42 + i,
            max_iter=20000,
        )
        # Wrap in CalibratedClassifierCV for probability estimates
        svm = CalibratedClassifierCV(base_svm, method="sigmoid", cv=3)
        svm.fit(X, y)

        model = SVMModel(
            model=svm,
            train_size=len(X),
            u12_count=len(sample_u12_introns),
            u2_count=len(sample_u2_introns),
            parameters=SVMParameters(
                C=1.0,
                calibration_method="sigmoid",
                saturate_enabled=False,
                include_max=False,
                include_pairwise_mins=False,
                penalty="l2",
                class_weight_multiplier=1.0,
                loss="squared_hinge",
                dual=False,
                intercept_scaling=1.0,
                cv_score=0.85,
                round_found=1,
            ),
        )
        models.append(model)

    return SVMEnsemble(models=models)


# Test SVMPredictor initialization


def test_predictor_initialization():
    """Test SVMPredictor initialization with default threshold."""
    predictor = SVMPredictor()
    assert predictor.threshold == 90.0


def test_predictor_custom_threshold():
    """Test SVMPredictor with custom threshold."""
    predictor = SVMPredictor(threshold=85.0)
    assert predictor.threshold == 85.0


def test_predictor_invalid_threshold():
    """Test SVMPredictor rejects invalid thresholds."""
    with pytest.raises(ValueError, match="must be between 0 and 100"):
        SVMPredictor(threshold=-10)

    with pytest.raises(ValueError, match="must be between 0 and 100"):
        SVMPredictor(threshold=150)


# Test feature extraction


def test_prepare_features(sample_u12_introns):
    """Test feature matrix extraction."""
    predictor = SVMPredictor()
    X = predictor._prepare_features(sample_u12_introns)

    assert X.shape == (10, 3)  # 10 introns, 3 features
    assert np.all(np.isfinite(X))  # No NaNs or Infs


def test_prepare_features_no_scores():
    """Test feature extraction fails without scores."""
    intron = Intron(
        intron_id="test",
        coordinates=GenomicCoordinate(
            chromosome="chr1", start=1000, stop=1100, strand="+", system="1-based"
        ),
        sequences=IntronSequences(seq="ATCG", five_seq="AT", three_seq="CG"),
        scores=None,  # No scores
    )

    predictor = SVMPredictor()
    with pytest.raises(ValueError, match="has no scores"):
        predictor._prepare_features([intron])


def test_prepare_features_missing_z_scores():
    """Test feature extraction fails with incomplete raw scores."""
    intron = Intron(
        intron_id="test",
        coordinates=GenomicCoordinate(
            chromosome="chr1", start=1000, stop=1100, strand="+", system="1-based"
        ),
        sequences=IntronSequences(seq="ATCG", five_seq="AT", three_seq="CG"),
        scores=IntronScores(
            five_raw_score=1.0,
            bp_raw_score=None,  # Missing
            three_raw_score=1.0,
        ),
    )

    predictor = SVMPredictor()
    with pytest.raises(ValueError, match="missing raw scores"):
        predictor._prepare_features([intron])


# Test prediction


def test_predict_u12_introns(trained_ensemble, sample_u12_introns):
    """Test prediction on U12-type introns (should classify as U12)."""
    predictor = SVMPredictor(threshold=50.0)  # Lower threshold for testing
    classified = predictor.predict(trained_ensemble, sample_u12_introns)

    assert len(classified) == len(sample_u12_introns)

    # Check that U12 introns get high scores
    for intron in classified:
        assert intron.scores.svm_score is not None
        assert 0 <= intron.scores.svm_score <= 100
        # Should be classified as U12 with high score
        assert intron.scores.svm_score > 50.0


def test_predict_u2_introns(trained_ensemble, sample_u2_introns):
    """Test prediction on U2-type introns (should classify as U2)."""
    predictor = SVMPredictor(threshold=50.0)
    classified = predictor.predict(trained_ensemble, sample_u2_introns)

    assert len(classified) == len(sample_u2_introns)

    # Check that U2 introns get low scores
    for intron in classified:
        assert intron.scores.svm_score is not None
        assert 0 <= intron.scores.svm_score <= 100
        # Should be classified as U2 with low score
        assert intron.scores.svm_score < 50.0


def test_predict_assigns_type_id(
    trained_ensemble, sample_u12_introns, sample_u2_introns
):
    """Test that type_id is assigned based on threshold."""
    predictor = SVMPredictor(threshold=50.0)

    # Classify both types
    u12_classified = predictor.predict(trained_ensemble, sample_u12_introns)
    u2_classified = predictor.predict(trained_ensemble, sample_u2_introns)

    # U12 introns should have type_id = 'u12'
    for intron in u12_classified:
        assert intron.metadata is not None
        assert intron.metadata.type_id == "u12"

    # U2 introns should have type_id = 'u2'
    for intron in u2_classified:
        assert intron.metadata is not None
        assert intron.metadata.type_id == "u2"


def test_predict_threshold_affects_classification(trained_ensemble, sample_u12_introns):
    """Test that changing threshold affects type assignment."""
    # With very low threshold, all should be U12
    predictor_low = SVMPredictor(threshold=10.0)
    classified_low = predictor_low.predict(trained_ensemble, sample_u12_introns)
    u12_count_low = sum(
        1 for i in classified_low if i.metadata and i.metadata.type_id == "u12"
    )

    # With very high threshold, fewer should be U12
    predictor_high = SVMPredictor(threshold=99.0)
    classified_high = predictor_high.predict(trained_ensemble, sample_u12_introns)
    u12_count_high = sum(
        1 for i in classified_high if i.metadata and i.metadata.type_id == "u12"
    )

    # Lower threshold should result in more U12 classifications
    assert u12_count_low >= u12_count_high


def test_predict_calculates_relative_score(trained_ensemble, sample_u12_introns):
    """Test that relative_score (decision function) is calculated."""
    predictor = SVMPredictor()
    classified = predictor.predict(trained_ensemble, sample_u12_introns)

    for intron in classified:
        assert intron.scores.relative_score is not None
        assert np.isfinite(intron.scores.relative_score)


def test_predict_ensemble_averaging(sample_u12_introns, sample_u2_introns):
    """Test that predictions use ensemble averaging."""
    # Train two models with enough data for calibration
    train_X = np.array(
        [
            [2.0, 2.5, 2.0],  # U12
            [2.1, 2.6, 2.1],  # U12
            [2.2, 2.7, 2.2],  # U12
            [-1.0, -1.5, -1.0],  # U2
            [-1.1, -1.6, -1.1],  # U2
            [-1.2, -1.7, -1.2],  # U2
        ]
    )
    train_y = np.array([1, 1, 1, 0, 0, 0])

    # Create two different models
    models = []
    for seed in [42, 43]:
        base_svm = LinearSVC(
            C=1.0,
            class_weight="balanced",
            dual="auto",
            random_state=seed,
            max_iter=20000,
        )
        svm = CalibratedClassifierCV(base_svm, method="sigmoid", cv=2)
        svm.fit(train_X, train_y)

        model = SVMModel(
            model=svm,
            train_size=len(train_X),
            u12_count=3,
            u2_count=3,
            parameters=SVMParameters(
                C=1.0,
                calibration_method="sigmoid",
                saturate_enabled=False,
                include_max=False,
                include_pairwise_mins=False,
                penalty="l2",
                class_weight_multiplier=1.0,
                loss="squared_hinge",
                dual=False,
                intercept_scaling=1.0,
                cv_score=0.9,
                round_found=1,
            ),
        )
        models.append(model)

    ensemble = SVMEnsemble(models=models)

    # Test intron
    test_intron = sample_u12_introns[0]

    predictor = SVMPredictor()
    classified = predictor.predict(ensemble, [test_intron])

    # Should have predictions
    assert len(classified) == 1
    assert classified[0].scores.svm_score is not None


def test_predict_empty_ensemble():
    """Test prediction with empty ensemble."""
    empty_ensemble = SVMEnsemble(models=[])

    predictor = SVMPredictor()
    intron = Intron(
        intron_id="test",
        coordinates=GenomicCoordinate(
            chromosome="chr1", start=1000, stop=1100, strand="+", system="1-based"
        ),
        sequences=IntronSequences(seq="ATCG", five_seq="AT", three_seq="CG"),
        scores=IntronScores(five_z_score=1.0, bp_z_score=1.0, three_z_score=1.0),
    )

    with pytest.raises(ValueError, match="Ensemble has no models"):
        predictor.predict(empty_ensemble, [intron])


# Test batch prediction


def test_predict_batch(trained_ensemble, sample_u12_introns, sample_u2_introns):
    """Test batch prediction produces same results as regular prediction."""
    all_introns = sample_u12_introns + sample_u2_introns

    predictor = SVMPredictor()

    # Regular prediction
    regular = predictor.predict(trained_ensemble, all_introns)

    # Batch prediction with small batch size
    batch = predictor.predict_batch(trained_ensemble, all_introns, batch_size=5)

    assert len(regular) == len(batch)

    # Check scores match
    for r, b in zip(regular, batch):
        assert abs(r.scores.svm_score - b.scores.svm_score) < 1e-6
        assert r.metadata.type_id == b.metadata.type_id


def test_predict_batch_single_batch(trained_ensemble, sample_u12_introns):
    """Test batch prediction with batch_size >= n_introns."""
    predictor = SVMPredictor()

    # Batch size larger than dataset
    batch = predictor.predict_batch(
        trained_ensemble, sample_u12_introns, batch_size=100
    )

    assert len(batch) == len(sample_u12_introns)


# Test preservation of original data


def test_predict_preserves_intron_data(trained_ensemble, sample_u12_introns):
    """Test that prediction preserves original intron data."""
    original = sample_u12_introns[0]

    predictor = SVMPredictor()
    classified = predictor.predict(trained_ensemble, sample_u12_introns)

    result = classified[0]

    # Check that basic data is preserved
    assert result.intron_id == original.intron_id
    assert result.coordinates == original.coordinates
    assert result.sequences == original.sequences

    # Check that original z-scores are preserved
    assert result.scores.five_z_score == original.scores.five_z_score
    assert result.scores.bp_z_score == original.scores.bp_z_score
    assert result.scores.three_z_score == original.scores.three_z_score


def test_predict_does_not_renormalize(trained_ensemble, sample_u12_introns):
    """
    CRITICAL TEST: Verify z-scores are NOT re-normalized after classification.
    This is Issue #1 fix - prevents data leakage.
    """
    original_z_scores = [
        (i.scores.five_z_score, i.scores.bp_z_score, i.scores.three_z_score)
        for i in sample_u12_introns
    ]

    predictor = SVMPredictor()
    classified = predictor.predict(trained_ensemble, sample_u12_introns)

    classified_z_scores = [
        (i.scores.five_z_score, i.scores.bp_z_score, i.scores.three_z_score)
        for i in classified
    ]

    # Z-scores should be EXACTLY the same (no re-normalization)
    for orig, classif in zip(original_z_scores, classified_z_scores):
        assert orig == classif, "Z-scores were modified during classification!"


# Test edge cases


def test_predict_single_intron(trained_ensemble, sample_u12_introns):
    """Test prediction on single intron."""
    predictor = SVMPredictor()
    classified = predictor.predict(trained_ensemble, sample_u12_introns[:1])

    assert len(classified) == 1
    assert classified[0].scores.svm_score is not None
    assert classified[0].metadata is not None
    assert classified[0].metadata.type_id in ["u2", "u12"]


def test_predict_boundary_scores(trained_ensemble):
    """Test introns with scores very close to threshold."""
    # Create intron with features that should give ~50% probability
    intron = Intron(
        intron_id="boundary",
        coordinates=GenomicCoordinate(
            chromosome="chr1", start=1000, stop=1100, strand="+", system="1-based"
        ),
        sequences=IntronSequences(
            seq="GTAAGT" + "N" * 50 + "CAG",
            five_seq="GTAAGT",
            three_seq="CAG",
            bp_seq="CTAAC",
        ),
        scores=IntronScores(
            five_z_score=0.5,  # Medium score
            bp_z_score=0.5,
            three_z_score=0.5,
        ),
    )

    predictor = SVMPredictor(threshold=90.0)
    classified = predictor.predict(trained_ensemble, [intron])

    assert len(classified) == 1
    # Score should be between 0 and 100
    assert 0 <= classified[0].scores.svm_score <= 100


def test_predict_extreme_z_scores(trained_ensemble):
    """Test introns with very high/low z-scores."""
    # Very high z-scores (strong U12)
    high_intron = Intron(
        intron_id="high",
        coordinates=GenomicCoordinate(
            chromosome="chr1", start=1000, stop=1100, strand="+", system="1-based"
        ),
        sequences=IntronSequences(seq="ATCG", five_seq="AT", three_seq="CG"),
        scores=IntronScores(five_z_score=10.0, bp_z_score=10.0, three_z_score=10.0),
    )

    # Very low z-scores (strong U2)
    low_intron = Intron(
        intron_id="low",
        coordinates=GenomicCoordinate(
            chromosome="chr1", start=2000, stop=2100, strand="+", system="1-based"
        ),
        sequences=IntronSequences(seq="ATCG", five_seq="AT", three_seq="CG"),
        scores=IntronScores(five_z_score=-10.0, bp_z_score=-10.0, three_z_score=-10.0),
    )

    predictor = SVMPredictor(threshold=50.0)
    classified = predictor.predict(trained_ensemble, [high_intron, low_intron])

    # High z-scores should give high SVM score
    assert classified[0].scores.svm_score > 50.0
    assert classified[0].metadata is not None
    assert classified[0].metadata.type_id == "u12"

    # Low z-scores should give low SVM score
    assert classified[1].scores.svm_score < 50.0
    assert classified[1].metadata is not None
    assert classified[1].metadata.type_id == "u2"


# ============================================================================
# Streaming Classification Tests
# ============================================================================


def test_classify_introns_streaming_basic(trained_ensemble):
    """Test streaming classification produces same results as batch."""
    from intronIC.classification.predictor import classify_introns_streaming

    # Create test introns with z-scores
    introns = [
        Intron(
            intron_id="high",
            coordinates=GenomicCoordinate(
                chromosome="chr1", start=1000, stop=1100, strand="+", system="1-based"
            ),
            sequences=IntronSequences(seq="ATCG", five_seq="AT", three_seq="CG"),
            scores=IntronScores(five_z_score=2.0, bp_z_score=2.5, three_z_score=2.0),
        ),
        Intron(
            intron_id="low",
            coordinates=GenomicCoordinate(
                chromosome="chr1", start=2000, stop=2100, strand="+", system="1-based"
            ),
            sequences=IntronSequences(seq="ATCG", five_seq="AT", three_seq="CG"),
            scores=IntronScores(five_z_score=-2.0, bp_z_score=-2.5, three_z_score=-2.0),
        ),
    ]

    # Classify using streaming
    results = list(
        classify_introns_streaming(iter(introns), trained_ensemble, threshold=90.0)
    )

    assert len(results) == 2

    # Check classification results
    for result in results:
        assert result.scores.svm_score is not None
        assert result.scores.relative_score is not None
        assert result.scores.decision_distance is not None
        assert result.metadata.type_id in ("u12", "u2")


def test_classify_introns_streaming_is_generator(trained_ensemble):
    """Test that streaming classification returns a generator."""
    from intronIC.classification.predictor import classify_introns_streaming

    introns = [
        Intron(
            intron_id="test",
            coordinates=GenomicCoordinate(
                chromosome="chr1", start=1000, stop=1100, strand="+", system="1-based"
            ),
            sequences=IntronSequences(seq="ATCG"),
            scores=IntronScores(five_z_score=1.0, bp_z_score=1.0, three_z_score=1.0),
        ),
    ]

    result = classify_introns_streaming(iter(introns), trained_ensemble)

    # Should be a generator
    assert hasattr(result, "__iter__")
    assert hasattr(result, "__next__")


def test_classify_introns_streaming_preserves_metadata(trained_ensemble):
    """Test that streaming classification preserves original metadata."""
    from intronIC.classification.predictor import classify_introns_streaming
    from intronIC.core.intron import IntronMetadata

    intron = Intron(
        intron_id="test_meta",
        coordinates=GenomicCoordinate(
            chromosome="chr1", start=1000, stop=1100, strand="+", system="1-based"
        ),
        sequences=IntronSequences(seq="ATCG"),
        scores=IntronScores(five_z_score=1.0, bp_z_score=1.0, three_z_score=1.0),
        metadata=IntronMetadata(parent="transcript_1", grandparent="gene_1"),
    )

    results = list(classify_introns_streaming([intron], trained_ensemble))
    result = results[0]

    # Original metadata should be preserved
    assert result.metadata.parent == "transcript_1"
    assert result.metadata.grandparent == "gene_1"
    # type_id should be added
    assert result.metadata.type_id in ("u12", "u2")


def test_classify_introns_batch_basic(trained_ensemble):
    """Test batch classification function."""
    from intronIC.classification.predictor import classify_introns_batch

    introns = [
        Intron(
            intron_id="batch_1",
            coordinates=GenomicCoordinate(
                chromosome="chr1", start=1000, stop=1100, strand="+", system="1-based"
            ),
            sequences=IntronSequences(seq="ATCG"),
            scores=IntronScores(five_z_score=2.0, bp_z_score=2.5, three_z_score=2.0),
        ),
        Intron(
            intron_id="batch_2",
            coordinates=GenomicCoordinate(
                chromosome="chr1", start=2000, stop=2100, strand="+", system="1-based"
            ),
            sequences=IntronSequences(seq="ATCG"),
            scores=IntronScores(five_z_score=-2.0, bp_z_score=-2.5, three_z_score=-2.0),
        ),
    ]

    results = classify_introns_batch(introns, trained_ensemble, threshold=90.0)

    assert isinstance(results, list)
    assert len(results) == 2
    for result in results:
        assert result.scores.svm_score is not None
        assert result.metadata.type_id in ("u12", "u2")


def test_classify_introns_batch_empty(trained_ensemble):
    """Test batch classification with empty list."""
    from intronIC.classification.predictor import classify_introns_batch

    results = classify_introns_batch([], trained_ensemble)

    assert results == []


def test_streaming_matches_batch(trained_ensemble):
    """Test that streaming and batch produce identical results."""
    from intronIC.classification.predictor import (
        classify_introns_batch,
        classify_introns_streaming,
    )

    introns = [
        Intron(
            intron_id=f"test_{i}",
            coordinates=GenomicCoordinate(
                chromosome="chr1",
                start=1000 + i * 100,
                stop=1100 + i * 100,
                strand="+",
                system="1-based",
            ),
            sequences=IntronSequences(seq="ATCG"),
            scores=IntronScores(
                five_z_score=float(i - 2),
                bp_z_score=float(i - 1),
                three_z_score=float(i - 2),
            ),
        )
        for i in range(5)
    ]

    streaming_results = list(
        classify_introns_streaming(iter(introns), trained_ensemble)
    )
    batch_results = classify_introns_batch(introns, trained_ensemble)

    assert len(streaming_results) == len(batch_results)

    for s_intron, b_intron in zip(streaming_results, batch_results):
        # SVM scores should match
        assert np.isclose(
            s_intron.scores.svm_score, b_intron.scores.svm_score, rtol=1e-5
        )
        assert np.isclose(
            s_intron.scores.relative_score, b_intron.scores.relative_score, rtol=1e-5
        )
        assert np.isclose(
            s_intron.scores.decision_distance,
            b_intron.scores.decision_distance,
            rtol=1e-5,
        )
        # Type should match
        assert s_intron.metadata.type_id == b_intron.metadata.type_id
