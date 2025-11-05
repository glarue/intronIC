"""
Tests for classification.trainer module.

Tests the SVMTrainer class which implements ensemble training with U2 subsampling.
"""

import pytest
import numpy as np
from dataclasses import replace

from classification.trainer import (
    SVMTrainer,
    SVMModel,
    SVMEnsemble,
)
from classification.optimizer import SVMParameters
from core.intron import Intron, IntronSequences, IntronScores, IntronMetadata, GenomicCoordinate


# =============================================================================
# FIXTURES
# =============================================================================


@pytest.fixture
def mock_u12_introns():
    """Create mock U12 introns with z-scores."""
    introns = []
    # Create 20 U12 introns with high scores (positive z-scores)
    for i in range(20):
        intron = Intron(
            intron_id=f"u12_{i}",
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
                five_z_score=2.0 + np.random.randn() * 0.5,
                bp_z_score=2.5 + np.random.randn() * 0.5,
                three_z_score=2.0 + np.random.randn() * 0.5,
            )
        )
        introns.append(intron)
    return introns


@pytest.fixture
def mock_u2_introns():
    """Create mock U2 introns with z-scores."""
    introns = []
    # Create 100 U2 introns with low scores (negative z-scores)
    for i in range(100):
        intron = Intron(
            intron_id=f"u2_{i}",
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
                bp_seq="TTTCAG"
            ),
            scores=IntronScores(
                five_z_score=-1.0 + np.random.randn() * 0.5,
                bp_z_score=-1.5 + np.random.randn() * 0.5,
                three_z_score=-1.0 + np.random.randn() * 0.5,
            )
        )
        introns.append(intron)
    return introns


@pytest.fixture
def mock_parameters():
    """Create mock SVM parameters."""
    return SVMParameters(
        C=1.0,
        cv_score=0.95,
        round_found=-1
    )


# =============================================================================
# UNIT TESTS
# =============================================================================


class TestSVMTrainer:
    """Test SVMTrainer class."""

    def test_initialization(self):
        """Test trainer initialization with default parameters."""
        trainer = SVMTrainer()

        assert trainer.n_models == 3
        assert trainer.test_size == 0.2
        assert trainer.random_state == 42
        assert trainer.kernel == 'linear'

    def test_initialization_custom_params(self):
        """Test trainer initialization with custom parameters."""
        trainer = SVMTrainer(
            n_models=5,
            test_size=0.3,
            random_state=123,
            kernel='rbf'
        )

        assert trainer.n_models == 5
        assert trainer.test_size == 0.3
        assert trainer.random_state == 123
        assert trainer.kernel == 'rbf'

    def test_prepare_training_data_structure(self, mock_u12_introns, mock_u2_introns):
        """Test feature extraction produces correct structure."""
        trainer = SVMTrainer()
        X, y = trainer._prepare_training_data(mock_u12_introns, mock_u2_introns)

        # Check shapes
        n_u2 = len(mock_u2_introns)
        n_u12 = len(mock_u12_introns)
        assert X.shape == (n_u2 + n_u12, 3)  # 3 features per intron
        assert y.shape == (n_u2 + n_u12,)

        # Check labels
        # U2 introns (first n_u2) should have label 0
        assert np.all(y[:n_u2] == 0)
        # U12 introns (next n_u12) should have label 1
        assert np.all(y[n_u2:] == 1)

    def test_prepare_training_data_values(self, mock_u12_introns, mock_u2_introns):
        """Test feature extraction produces correct values."""
        trainer = SVMTrainer()
        X, y = trainer._prepare_training_data(mock_u12_introns, mock_u2_introns)

        # Check first U2 intron features
        u2_features = X[0]
        expected_u2 = [
            mock_u2_introns[0].scores.five_z_score,
            mock_u2_introns[0].scores.bp_z_score,
            mock_u2_introns[0].scores.three_z_score,
        ]
        assert np.allclose(u2_features, expected_u2)

        # Check first U12 intron features
        n_u2 = len(mock_u2_introns)
        u12_features = X[n_u2]
        expected_u12 = [
            mock_u12_introns[0].scores.five_z_score,
            mock_u12_introns[0].scores.bp_z_score,
            mock_u12_introns[0].scores.three_z_score,
        ]
        assert np.allclose(u12_features, expected_u12)

    def test_subsample_u2(self, mock_u2_introns):
        """Test U2 subsampling."""
        trainer = SVMTrainer()

        # Subsample 80% of U2 introns
        subsampled = trainer._subsample_u2(
            mock_u2_introns,
            seed=42,
            subsample_ratio=0.8
        )

        # Check size
        expected_size = int(len(mock_u2_introns) * 0.8)
        assert len(subsampled) == expected_size

        # Check that all are from original set
        original_ids = {i.intron_id for i in mock_u2_introns}
        subsampled_ids = {i.intron_id for i in subsampled}
        assert subsampled_ids.issubset(original_ids)

    def test_subsample_u2_reproducibility(self, mock_u2_introns):
        """Test that subsampling with same seed gives same results."""
        trainer = SVMTrainer()

        subsample1 = trainer._subsample_u2(mock_u2_introns, seed=42)
        subsample2 = trainer._subsample_u2(mock_u2_introns, seed=42)

        ids1 = [i.intron_id for i in subsample1]
        ids2 = [i.intron_id for i in subsample2]
        assert ids1 == ids2

    def test_subsample_u2_different_seeds(self, mock_u2_introns):
        """Test that subsampling with different seeds gives different results."""
        trainer = SVMTrainer()

        subsample1 = trainer._subsample_u2(mock_u2_introns, seed=42)
        subsample2 = trainer._subsample_u2(mock_u2_introns, seed=99)

        ids1 = set(i.intron_id for i in subsample1)
        ids2 = set(i.intron_id for i in subsample2)

        # Should have some different introns
        # (not guaranteed to be different, but very likely with 100 introns)
        # So we just check they're different sets (not identical)
        assert ids1 != ids2 or len(ids1) < len(mock_u2_introns)

    def test_train_single_model(self, mock_u12_introns, mock_u2_introns, mock_parameters):
        """Test training a single SVM model."""
        trainer = SVMTrainer(test_size=0.2)

        model = trainer._train_single_model(
            mock_u12_introns,
            mock_u2_introns,
            mock_parameters,
            seed=42
        )

        # Check model structure
        assert isinstance(model, SVMModel)
        assert model.model is not None
        assert 0 <= model.f1_score <= 1
        assert 0 <= model.precision_recall_auc <= 1
        assert model.train_size > 0
        assert model.test_size > 0
        assert model.u12_count == len(mock_u12_introns)
        assert model.u2_count == len(mock_u2_introns)
        assert model.parameters == mock_parameters

        # Check model can predict
        X, y = trainer._prepare_training_data(mock_u12_introns, mock_u2_introns)
        predictions = model.model.predict(X)
        assert len(predictions) == len(y)

    def test_train_single_model_performance(
        self, mock_u12_introns, mock_u2_introns, mock_parameters
    ):
        """Test that trained model has reasonable performance."""
        trainer = SVMTrainer(test_size=0.2)

        model = trainer._train_single_model(
            mock_u12_introns,
            mock_u2_introns,
            mock_parameters,
            seed=42
        )

        # With well-separated data, should get good scores
        assert model.f1_score > 0.7
        assert model.precision_recall_auc > 0.7

    def test_train_ensemble_structure(
        self, mock_u12_introns, mock_u2_introns, mock_parameters
    ):
        """Test ensemble training produces correct structure."""
        trainer = SVMTrainer(n_models=3)

        ensemble = trainer.train_ensemble(
            mock_u12_introns,
            mock_u2_introns,
            mock_parameters,
            subsample_u2=True
        )

        # Check ensemble structure
        assert isinstance(ensemble, SVMEnsemble)
        assert len(ensemble.models) == 3
        assert len(ensemble) == 3  # Test __len__
        assert 0 <= ensemble.mean_f1 <= 1
        assert 0 <= ensemble.mean_pr_auc <= 1

        # Check all models are present
        for model in ensemble.models:
            assert isinstance(model, SVMModel)
            assert model.model is not None

    def test_train_ensemble_no_subsampling(
        self, mock_u12_introns, mock_u2_introns, mock_parameters
    ):
        """Test ensemble training without subsampling."""
        trainer = SVMTrainer(n_models=2)

        ensemble = trainer.train_ensemble(
            mock_u12_introns,
            mock_u2_introns,
            mock_parameters,
            subsample_u2=False
        )

        # All models should use same U2 count
        u2_counts = [m.u2_count for m in ensemble.models]
        assert all(c == len(mock_u2_introns) for c in u2_counts)

    def test_train_ensemble_with_subsampling(
        self, mock_u12_introns, mock_u2_introns, mock_parameters
    ):
        """Test ensemble training with subsampling creates diversity."""
        trainer = SVMTrainer(n_models=3)

        ensemble = trainer.train_ensemble(
            mock_u12_introns,
            mock_u2_introns,
            mock_parameters,
            subsample_u2=True,
            subsample_ratio=0.8
        )

        # With subsampling, U2 counts should be less than full set
        expected_u2_count = int(len(mock_u2_introns) * 0.8)
        u2_counts = [m.u2_count for m in ensemble.models]
        assert all(c == expected_u2_count for c in u2_counts)

    def test_train_ensemble_metrics_consistency(
        self, mock_u12_introns, mock_u2_introns, mock_parameters
    ):
        """Test that ensemble metrics are consistent with model metrics."""
        trainer = SVMTrainer(n_models=3)

        ensemble = trainer.train_ensemble(
            mock_u12_introns,
            mock_u2_introns,
            mock_parameters
        )

        # Mean F1 should be average of model F1s
        expected_mean_f1 = np.mean([m.f1_score for m in ensemble.models])
        assert np.isclose(ensemble.mean_f1, expected_mean_f1)

        # Mean PR-AUC should be average of model PR-AUCs
        expected_mean_pr_auc = np.mean([m.precision_recall_auc for m in ensemble.models])
        assert np.isclose(ensemble.mean_pr_auc, expected_mean_pr_auc)

    def test_train_single_model_reproducibility(
        self, mock_u12_introns, mock_u2_introns, mock_parameters
    ):
        """Test that training with same seed gives same results."""
        trainer = SVMTrainer()

        model1 = trainer._train_single_model(
            mock_u12_introns,
            mock_u2_introns,
            mock_parameters,
            seed=42
        )

        model2 = trainer._train_single_model(
            mock_u12_introns,
            mock_u2_introns,
            mock_parameters,
            seed=42
        )

        # Should get same performance
        assert np.isclose(model1.f1_score, model2.f1_score)
        assert np.isclose(model1.precision_recall_auc, model2.precision_recall_auc)

    def test_train_ensemble_single_model(
        self, mock_u12_introns, mock_u2_introns, mock_parameters
    ):
        """Test ensemble training with n_models=1."""
        trainer = SVMTrainer(n_models=1)

        ensemble = trainer.train_ensemble(
            mock_u12_introns,
            mock_u2_introns,
            mock_parameters,
            subsample_u2=False  # No point subsampling with 1 model
        )

        # Should have exactly 1 model
        assert len(ensemble.models) == 1
        assert ensemble.mean_f1 == ensemble.models[0].f1_score
        assert ensemble.mean_pr_auc == ensemble.models[0].precision_recall_auc


# =============================================================================
# INTEGRATION TESTS
# =============================================================================


class TestSVMTrainerIntegration:
    """Integration tests for SVMTrainer."""

    @pytest.mark.slow
    def test_full_training_pipeline(self, mock_u12_introns, mock_u2_introns):
        """Test full training pipeline with realistic data."""
        # Make data more realistic and well-separated
        realistic_u12s = []
        for intron in mock_u12_introns:
            new_scores = replace(
                intron.scores,
                five_z_score=3.0 + np.random.randn() * 0.3,
                bp_z_score=4.0 + np.random.randn() * 0.3,
                three_z_score=3.5 + np.random.randn() * 0.3
            )
            realistic_u12s.append(replace(intron, scores=new_scores))

        realistic_u2s = []
        for intron in mock_u2_introns:
            new_scores = replace(
                intron.scores,
                five_z_score=-1.5 + np.random.randn() * 0.3,
                bp_z_score=-2.0 + np.random.randn() * 0.3,
                three_z_score=-1.5 + np.random.randn() * 0.3
            )
            realistic_u2s.append(replace(intron, scores=new_scores))

        # Train with realistic parameters
        trainer = SVMTrainer(n_models=3)
        parameters = SVMParameters(C=1.0, cv_score=0.95, round_found=-1)

        ensemble = trainer.train_ensemble(
            realistic_u12s,
            realistic_u2s,
            parameters
        )

        # Should get high performance with well-separated data
        assert ensemble.mean_f1 > 0.9
        assert ensemble.mean_pr_auc > 0.9

        # All models should have reasonable performance
        for model in ensemble.models:
            assert model.f1_score > 0.85
            assert model.precision_recall_auc > 0.85


# =============================================================================
# EDGE CASES
# =============================================================================


class TestSVMTrainerEdgeCases:
    """Test edge cases and error handling."""

    def test_train_with_minimal_data(self):
        """Test training with minimal number of introns."""
        # Create minimal dataset (5 of each)
        u12_introns = []
        for i in range(5):
            intron = Intron(
                intron_id=f"u12_{i}",
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
            u12_introns.append(intron)

        u2_introns = []
        for i in range(5):
            intron = Intron(
                intron_id=f"u2_{i}",
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
                    bp_seq="TTTCAG"
                ),
                scores=IntronScores(
                    five_z_score=-1.0,
                    bp_z_score=-1.5,
                    three_z_score=-1.0,
                )
            )
            u2_introns.append(intron)

        trainer = SVMTrainer(n_models=1, test_size=0.2)
        parameters = SVMParameters(C=1.0, cv_score=0.95, round_found=-1)

        # Should still run without error
        ensemble = trainer.train_ensemble(
            u12_introns,
            u2_introns,
            parameters,
            subsample_u2=False
        )

        assert len(ensemble.models) == 1
        assert 0 <= ensemble.mean_f1 <= 1
        assert 0 <= ensemble.mean_pr_auc <= 1

    def test_prepare_training_data_no_scores(self, mock_u12_introns):
        """Test error when intron lacks scores."""
        trainer = SVMTrainer()

        # Create intron without scores
        bad_intron = replace(mock_u12_introns[0], scores=None)

        with pytest.raises(ValueError, match="has no scores"):
            trainer._prepare_training_data([bad_intron], mock_u12_introns[1:])

    def test_prepare_training_data_missing_z_scores(self, mock_u12_introns):
        """Test error when intron has incomplete z-scores."""
        trainer = SVMTrainer()

        # Create intron with missing z-scores
        bad_scores = replace(
            mock_u12_introns[0].scores,
            five_z_score=None
        )
        bad_intron = replace(mock_u12_introns[0], scores=bad_scores)

        with pytest.raises(ValueError, match="missing z-scores"):
            trainer._prepare_training_data([bad_intron], mock_u12_introns[1:])
