"""
Tests for classification.optimizer module.

Tests the SVMOptimizer class which implements geometric grid search for
SVM hyperparameter optimization.
"""

from dataclasses import replace

import numpy as np
import pytest
from scipy.stats import gmean

from intronIC.classification.optimizer import (
    OptimizationRound,
    SVMOptimizer,
    SVMParameters,
)
from intronIC.core.intron import (
    GenomicCoordinate,
    Intron,
    IntronMetadata,
    IntronScores,
    IntronSequences,
)

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
                system="1-based",
            ),
            sequences=IntronSequences(
                seq="GTATGT" + "N" * 50 + "TCCTTAAC",
                five_seq="GTATGT",
                three_seq="TCCTTAAC",
                bp_seq="TCCTTAAC",
            ),
            scores=IntronScores(
                five_z_score=2.0 + np.random.randn() * 0.5,
                bp_z_score=2.5 + np.random.randn() * 0.5,
                three_z_score=2.0 + np.random.randn() * 0.5,
            ),
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
                system="1-based",
            ),
            sequences=IntronSequences(
                seq="GTAAGT" + "N" * 50 + "TTTCAG",
                five_seq="GTAAGT",
                three_seq="TTTCAG",
                bp_seq="TTTCAG",
            ),
            scores=IntronScores(
                five_z_score=-1.0 + np.random.randn() * 0.5,
                bp_z_score=-1.5 + np.random.randn() * 0.5,
                three_z_score=-1.0 + np.random.randn() * 0.5,
            ),
        )
        introns.append(intron)
    return introns


# =============================================================================
# UNIT TESTS
# =============================================================================


class TestSVMOptimizer:
    """Test SVMOptimizer class."""

    def test_initialization(self):
        """Test optimizer initialization with default parameters."""
        optimizer = SVMOptimizer()

        assert optimizer.n_rounds == 3
        assert optimizer.n_points_initial == 13
        assert optimizer.n_points_refine == 100
        assert optimizer.cv_folds == 7  # Updated default
        assert optimizer.random_state == 42
        assert optimizer.rounds_ == []

    def test_initialization_custom_params(self):
        """Test optimizer initialization with custom parameters."""
        optimizer = SVMOptimizer(
            n_rounds=3,
            n_points_initial=10,
            n_points_refine=50,
            cv_folds=3,
            random_state=123,
        )

        assert optimizer.n_rounds == 3
        assert optimizer.n_points_initial == 10
        assert optimizer.n_points_refine == 50
        assert optimizer.cv_folds == 3
        assert optimizer.random_state == 123

    def test_create_initial_grid_default_range(self):
        """Test initial grid creation with default range."""
        optimizer = SVMOptimizer()
        grid = optimizer._create_initial_grid((1e-6, 1e6))

        # Should have 13 points (default)
        assert len(grid) == 13

        # Should span 10^-6 to 10^6
        assert np.isclose(grid[0], 1e-6)
        assert np.isclose(grid[-1], 1e6)

        # Should be logarithmically spaced
        log_grid = np.log10(grid)
        log_diffs = np.diff(log_grid)
        assert np.allclose(log_diffs, log_diffs[0], rtol=1e-10)

    def test_create_initial_grid_custom_range(self):
        """Test initial grid creation with custom range."""
        optimizer = SVMOptimizer(n_points_initial=7)
        grid = optimizer._create_initial_grid((1e-3, 1e3))

        assert len(grid) == 7
        assert np.isclose(grid[0], 1e-3)
        assert np.isclose(grid[-1], 1e3)

    def test_index_of_nearest(self):
        """Test finding index of nearest value."""
        optimizer = SVMOptimizer()
        array = np.array([1, 5, 10, 20, 50, 100])

        # Exact matches
        assert optimizer._index_of_nearest(array, 10) == 2
        assert optimizer._index_of_nearest(array, 50) == 4

        # Nearest values
        assert optimizer._index_of_nearest(array, 11) == 2  # closer to 10
        assert optimizer._index_of_nearest(array, 16) == 3  # closer to 20
        assert optimizer._index_of_nearest(array, 30) == 3  # closer to 20
        assert optimizer._index_of_nearest(array, 60) == 4  # closer to 50

    def test_prepare_training_data_structure(self, mock_u12_introns, mock_u2_introns):
        """Test feature extraction produces correct structure."""
        optimizer = SVMOptimizer()
        X, y = optimizer._prepare_training_data(mock_u12_introns, mock_u2_introns)

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
        optimizer = SVMOptimizer()
        X, y = optimizer._prepare_training_data(mock_u12_introns, mock_u2_introns)

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

    def test_prepare_training_data_no_scores(self, mock_u12_introns):
        """Test error when intron lacks scores."""
        optimizer = SVMOptimizer()

        # Create intron without scores
        bad_intron = replace(mock_u12_introns[0], scores=None)

        with pytest.raises(ValueError, match="has no scores"):
            optimizer._prepare_training_data([bad_intron], mock_u12_introns[1:])

    def test_prepare_training_data_missing_z_scores(self, mock_u12_introns):
        """Test error when intron has incomplete z-scores."""
        optimizer = SVMOptimizer()

        # Create intron with missing z-scores
        bad_scores = replace(mock_u12_introns[0].scores, five_z_score=None)
        bad_intron = replace(mock_u12_introns[0], scores=bad_scores)

        with pytest.raises(ValueError, match="missing z-scores"):
            optimizer._prepare_training_data([bad_intron], mock_u12_introns[1:])

    def test_refine_grid_middle_value(self):
        """Test grid refinement around middle value."""
        optimizer = SVMOptimizer(n_points_refine=50)
        current_grid = np.array([1, 10, 100, 1000, 10000])
        best_C = 100  # Middle value

        refined_grid = optimizer._refine_grid(current_grid, best_C)

        # Should have 50 points
        assert len(refined_grid) == 50

        # Should span neighbors (10 to 1000)
        assert np.isclose(refined_grid[0], 10)
        assert np.isclose(refined_grid[-1], 1000)

    def test_refine_grid_edge_values(self):
        """Test grid refinement at edges."""
        optimizer = SVMOptimizer(n_points_refine=50)
        current_grid = np.array([1, 10, 100, 1000, 10000])

        # Test at low edge - new implementation extends below for edge cases
        refined_low = optimizer._refine_grid(current_grid, 1)
        # Should extend below 1 and go up to 10
        assert refined_low[0] < 1  # Extended below
        assert np.isclose(refined_low[-1], 10)

        # Test at high edge - new implementation extends above for edge cases
        refined_high = optimizer._refine_grid(current_grid, 10000)
        # Should start at 1000 and extend above 10000
        assert np.isclose(refined_high[0], 1000)
        assert refined_high[-1] > 10000  # Extended above

    @pytest.mark.slow
    def test_grid_search_returns_valid_scores(self, mock_u12_introns, mock_u2_introns):
        """Test that grid search produces valid scores."""
        optimizer = SVMOptimizer(cv_folds=3)
        X, y = optimizer._prepare_training_data(mock_u12_introns, mock_u2_introns)

        C_grid = np.logspace(-2, 2, num=5)
        round_result = optimizer._grid_search_round(X, y, C_grid, round_idx=0)

        # All scores should be between 0 and 1
        for score in round_result.scores:
            assert 0.0 <= score <= 1.0

        # With our well-separated mock data, best score should be high
        assert round_result.best_score > 0.8

    @pytest.mark.slow
    def test_grid_search_round(self, mock_u12_introns, mock_u2_introns):
        """Test single grid search round."""
        optimizer = SVMOptimizer(cv_folds=3)
        X, y = optimizer._prepare_training_data(mock_u12_introns, mock_u2_introns)

        C_grid = np.logspace(-2, 2, num=5)  # Small grid for speed
        round_result = optimizer._grid_search_round(X, y, C_grid, round_idx=0)

        # Check result structure
        assert isinstance(round_result, OptimizationRound)
        assert len(round_result.grid_points) == 5
        # New optimizer tests multiple parameter combinations per C value
        assert len(round_result.scores) >= 5  # At least one score per C
        assert round_result.best_C > 0
        assert 0 <= round_result.best_score <= 1
        assert len(round_result.rank_one_Cs) >= 1

        # Best C should be approximately in the grid range
        assert round_result.best_C >= C_grid[0] * 0.5  # Allow some tolerance
        assert round_result.best_C <= C_grid[-1] * 2.0

    @pytest.mark.slow
    def test_optimize_full_run(self, mock_u12_introns, mock_u2_introns):
        """Test full optimization with 2 rounds (faster for testing)."""
        optimizer = SVMOptimizer(
            n_rounds=2, n_points_initial=7, n_points_refine=10, cv_folds=3
        )

        params = optimizer.optimize(mock_u12_introns, mock_u2_introns)

        # Check result structure
        assert isinstance(params, SVMParameters)
        assert params.C > 0
        assert 0 <= params.cv_score <= 1
        assert params.round_found == -1  # Indicates averaged result

        # Should have completed 2 rounds
        assert len(optimizer.rounds_) == 2

        # Each round should have results
        for i, round_result in enumerate(optimizer.rounds_):
            assert len(round_result.grid_points) > 0
            assert len(round_result.scores) > 0
            assert round_result.best_C > 0

        # Rounds should refine (later rounds have tighter ranges)
        range_round1 = np.ptp(optimizer.rounds_[0].grid_points)
        range_round2 = np.ptp(optimizer.rounds_[1].grid_points)
        assert range_round2 < range_round1

    @pytest.mark.slow
    def test_optimize_convergence(self, mock_u12_introns, mock_u2_introns):
        """Test that optimization converges to similar values."""
        optimizer1 = SVMOptimizer(
            n_rounds=2,
            n_points_initial=7,
            n_points_refine=10,
            cv_folds=3,
            random_state=42,
        )

        optimizer2 = SVMOptimizer(
            n_rounds=2,
            n_points_initial=7,
            n_points_refine=10,
            cv_folds=3,
            random_state=42,
        )

        params1 = optimizer1.optimize(mock_u12_introns, mock_u2_introns)
        params2 = optimizer2.optimize(mock_u12_introns, mock_u2_introns)

        # With same random seed, should get identical results
        assert np.isclose(params1.C, params2.C)
        assert np.isclose(params1.cv_score, params2.cv_score)

    @pytest.mark.slow
    def test_optimize_reproducibility_with_seed(
        self, mock_u12_introns, mock_u2_introns
    ):
        """Test reproducibility with different random seeds."""
        optimizer1 = SVMOptimizer(
            n_rounds=1, n_points_initial=7, cv_folds=3, random_state=42
        )

        optimizer2 = SVMOptimizer(
            n_rounds=1, n_points_initial=7, cv_folds=3, random_state=99
        )

        params1 = optimizer1.optimize(mock_u12_introns, mock_u2_introns)
        params2 = optimizer2.optimize(mock_u12_introns, mock_u2_introns)

        # Different seeds might give slightly different results
        # but should be in same ballpark
        assert params1.C > 0 and params2.C > 0
        # Allow factor of 10 difference
        assert 0.1 < params1.C / params2.C < 10

    @pytest.mark.slow
    def test_optimize_uses_default_range(self, mock_u12_introns, mock_u2_introns):
        """Test optimization uses appropriate default range."""
        optimizer = SVMOptimizer(n_rounds=1, n_points_initial=5, cv_folds=3)

        params = optimizer.optimize(mock_u12_introns, mock_u2_introns)

        # Should find a valid C value
        assert params.C > 0
        assert 0 <= params.cv_score <= 1

        # Should have grid points
        initial_grid = optimizer.rounds_[0].grid_points
        assert len(initial_grid) == 5


# =============================================================================
# INTEGRATION TESTS
# =============================================================================


class TestSVMOptimizerIntegration:
    """Integration tests for SVMOptimizer."""

    @pytest.mark.slow
    def test_optimizer_with_realistic_data(self, mock_u12_introns, mock_u2_introns):
        """Test optimizer with realistic intron data."""
        # Adjust mock data to be more realistic
        # U12s should have higher raw scores, U2s lower
        # Since scores are frozen, we need to replace them
        realistic_u12s = []
        for intron in mock_u12_introns:
            new_scores = replace(
                intron.scores,
                five_raw_score=3.0 + np.random.randn() * 0.3,
                bp_raw_score=4.0 + np.random.randn() * 0.3,
                three_raw_score=3.5 + np.random.randn() * 0.3,
                five_z_score=3.0 + np.random.randn() * 0.3,
                bp_z_score=4.0 + np.random.randn() * 0.3,
                three_z_score=3.5 + np.random.randn() * 0.3,
            )
            realistic_u12s.append(replace(intron, scores=new_scores))

        realistic_u2s = []
        for intron in mock_u2_introns:
            new_scores = replace(
                intron.scores,
                five_raw_score=-1.5 + np.random.randn() * 0.3,
                bp_raw_score=-2.0 + np.random.randn() * 0.3,
                three_raw_score=-1.5 + np.random.randn() * 0.3,
                five_z_score=-1.5 + np.random.randn() * 0.3,
                bp_z_score=-2.0 + np.random.randn() * 0.3,
                three_z_score=-1.5 + np.random.randn() * 0.3,
            )
            realistic_u2s.append(replace(intron, scores=new_scores))

        optimizer = SVMOptimizer(
            n_rounds=3, n_points_initial=7, n_points_refine=20, cv_folds=3
        )

        params = optimizer.optimize(realistic_u12s, realistic_u2s)

        # With well-separated data, should get good calibration (low log-loss)
        # cv_score now represents log-loss (lower is better)
        assert params.cv_score < 0.5  # Log-loss should be low for well-separated data
        assert params.C > 0

        # Should have completed 3 rounds
        assert len(optimizer.rounds_) == 3

        # Each round should have good performance
        scores = [r.best_score for r in optimizer.rounds_]
        # Scores are balanced_accuracy, should be high
        assert all(s > 0.8 for s in scores)


# =============================================================================
# EDGE CASES
# =============================================================================


class TestSVMOptimizerEdgeCases:
    """Test edge cases and error handling."""

    @pytest.mark.slow
    def test_optimizer_with_minimal_data(self):
        """Test optimizer with minimal number of introns."""
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
                    system="1-based",
                ),
                sequences=IntronSequences(
                    seq="GTATGT" + "N" * 50 + "TCCTTAAC",
                    five_seq="GTATGT",
                    three_seq="TCCTTAAC",
                    bp_seq="TCCTTAAC",
                ),
                scores=IntronScores(
                    five_z_score=2.0,
                    bp_z_score=2.5,
                    three_z_score=2.0,
                ),
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
                    system="1-based",
                ),
                sequences=IntronSequences(
                    seq="GTAAGT" + "N" * 50 + "TTTCAG",
                    five_seq="GTAAGT",
                    three_seq="TTTCAG",
                    bp_seq="TTTCAG",
                ),
                scores=IntronScores(
                    five_z_score=-1.0,
                    bp_z_score=-1.5,
                    three_z_score=-1.0,
                ),
            )
            u2_introns.append(intron)

        optimizer = SVMOptimizer(
            n_rounds=1,
            n_points_initial=3,
            cv_folds=2,  # Small folds for small dataset
        )

        # Should still run without error
        params = optimizer.optimize(u12_introns, u2_introns)
        assert params.C > 0
        assert 0 <= params.cv_score <= 1
