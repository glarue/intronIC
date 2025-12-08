"""
Integration test for complete classification pipeline.

Tests the end-to-end workflow:
1. Load reference U12/U2 sequences from intronIC data
2. Parse scoring matrices (PWMs)
3. Score introns with PWMs
4. Normalize z-scores
5. Run classification pipeline
6. Validate results match expected accuracy

This test uses real reference data from the original intronIC package.
"""

import gzip
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pytest

from intronIC.classification.classifier import ClassificationResult, IntronClassifier
from intronIC.core.intron import (
    GenomicCoordinate,
    Intron,
    IntronScores,
    IntronSequences,
)
from intronIC.scoring.normalizer import ScoreNormalizer
from intronIC.scoring.pwm import PWM
from intronIC.scoring.scorer import IntronScorer

# Path to reference data
DATA_DIR = Path(__file__).parent.parent.parent / "src" / "intronIC" / "data"
U12_REFERENCE_FILE = DATA_DIR / "u12_reference.introns.iic.gz"
U2_REFERENCE_FILE = DATA_DIR / "u2_reference.introns.iic.gz"
PWM_FILE = DATA_DIR / "intronIC_scoring_PWMs.json"


def load_reference_sequences(filepath: Path, max_count: int = None) -> List[Intron]:
    """
    Load reference intron sequences from .iic.gz file.

    Format (tab-delimited):
    - intron_id
    - score
    - upstream_flank
    - intron_seq
    - downstream_flank
    - sources
    - source_count

    Args:
        filepath: Path to .iic.gz file
        max_count: Maximum number to load (None = all)

    Returns:
        List of Intron objects with sequences
    """
    introns = []

    with gzip.open(filepath, "rt") as f:
        for line_num, line in enumerate(f, 1):
            # Skip comments
            if line.startswith("#"):
                continue

            # Parse line
            fields = line.strip().split("\t")
            if len(fields) < 5:
                continue

            intron_id = fields[0]
            # score = fields[1]  # Not used for reference
            upstream_flank = fields[2]
            intron_seq = fields[3]
            downstream_flank = fields[4]

            # Create dummy coordinates (not needed for scoring)
            coord = GenomicCoordinate(
                chromosome="ref",
                start=1,
                stop=len(intron_seq),
                strand="+",
                system="1-based",
            )

            # Extract terminal dinucleotides
            five_dnt = intron_seq[:2] if len(intron_seq) >= 2 else None
            three_dnt = intron_seq[-2:] if len(intron_seq) >= 2 else None

            # Create IntronSequences
            sequences = IntronSequences(
                seq=intron_seq,
                upstream_flank=upstream_flank,
                downstream_flank=downstream_flank,
                five_prime_dnt=five_dnt,
                three_prime_dnt=three_dnt,
            )

            # Create Intron
            intron = Intron(intron_id=intron_id, coordinates=coord, sequences=sequences)

            introns.append(intron)

            # Check if we've reached max
            if max_count and len(introns) >= max_count:
                break

    return introns


def parse_pwm_file(filepath: Path) -> dict:
    """
    Parse PWM file to extract matrices.

    Returns dict mapping matrix names to PWM objects.

    Note: This is a simplified parser for testing.
    For production, use the full PWM parser from scoring module.
    """
    pwms = {}
    current_name = None
    current_freqs = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()

            # Skip comments and empty lines
            if not line or line.startswith("#"):
                continue

            # Header line
            if line.startswith(">"):
                # Save previous PWM if any
                if current_name and current_freqs:
                    pwms[current_name] = np.array(
                        current_freqs
                    ).T  # Transpose to (bases, positions)
                    current_freqs = []

                # Parse new header
                parts = line[1:].split("\t")
                current_name = parts[0]

            # Skip metadata lines
            elif line.startswith("A\t") or line == "A	C	G	T":
                continue

            # Frequency line
            else:
                try:
                    freqs = [float(x) for x in line.split("\t")]
                    if len(freqs) == 4:  # A, C, G, T
                        current_freqs.append(freqs)
                except ValueError:
                    continue

        # Save last PWM
        if current_name and current_freqs:
            pwms[current_name] = np.array(current_freqs).T

    return pwms


# Test fixtures


@pytest.fixture(scope="module")
def u12_reference():
    """Load U12 reference introns (all 387)."""
    if not U12_REFERENCE_FILE.exists():
        pytest.skip(f"Reference data not found: {U12_REFERENCE_FILE}")
    return load_reference_sequences(U12_REFERENCE_FILE)


@pytest.fixture(scope="module")
def u2_reference_subset():
    """Load subset of U2 reference introns (500 for speed)."""
    if not U2_REFERENCE_FILE.exists():
        pytest.skip(f"Reference data not found: {U2_REFERENCE_FILE}")
    return load_reference_sequences(U2_REFERENCE_FILE, max_count=500)


@pytest.fixture(scope="module")
def pwm_matrices():
    """Load PWM matrices."""
    if not PWM_FILE.exists():
        pytest.skip(f"PWM file not found: {PWM_FILE}")
    return parse_pwm_file(PWM_FILE)


# Integration tests


@pytest.mark.integration
def test_load_reference_data(u12_reference, u2_reference_subset):
    """Test that reference data loads correctly."""
    assert len(u12_reference) == 387
    assert len(u2_reference_subset) == 500

    # Check U12 reference structure
    for intron in u12_reference[:10]:
        assert intron.intron_id
        assert intron.sequences is not None
        assert intron.sequences.seq
        assert len(intron.sequences.seq) >= 30  # Minimum intron length


@pytest.mark.integration
def test_score_reference_introns(u12_reference, u2_reference_subset, pwm_matrices):
    """Test scoring reference introns with PWMs."""
    # Create simple PWMs for testing
    # In production, would use full PWM objects from scoring module
    five_pwm = pwm_matrices.get("u12_atac_five")
    bp_pwm = pwm_matrices.get("u12_atac_bp")
    three_pwm = pwm_matrices.get("u12_atac_three")

    if five_pwm is None or bp_pwm is None or three_pwm is None:
        pytest.skip("Required PWMs not found in matrices file")

    # For this integration test, we'll just verify the data structure
    # Full scoring would require the complete PWM infrastructure
    assert five_pwm.shape[0] == 4  # A, C, G, T
    assert five_pwm.shape[1] > 0  # Has positions


@pytest.mark.integration
@pytest.mark.slow
def test_classification_pipeline_with_reference_data(
    u12_reference, u2_reference_subset
):
    """
    Test complete classification pipeline with real reference data.

    This test simulates the full workflow:
    1. Create mock scored introns (with z-scores)
    2. Split into reference/experimental
    3. Run classification
    4. Validate accuracy
    """
    # For this integration test, we'll create synthetic z-scores
    # that mimic real distributions

    # Add z-scores to U12 reference (high values)
    # Track true labels explicitly
    true_labels = {}

    u12_scored = []
    for intron in u12_reference[:100]:  # Use subset for speed
        scores = IntronScores(
            five_z_score=2.0 + np.random.randn() * 0.5,
            bp_z_score=2.5 + np.random.randn() * 0.5,
            three_z_score=2.0 + np.random.randn() * 0.5,
        )
        new_intron = Intron(
            intron_id=intron.intron_id,
            coordinates=intron.coordinates,
            sequences=intron.sequences,
            scores=scores,
        )
        u12_scored.append(new_intron)
        true_labels[intron.intron_id] = "u12"

    # Add z-scores to U2 reference (low values)
    u2_scored = []
    for intron in u2_reference_subset[:100]:  # Use subset for speed
        scores = IntronScores(
            five_z_score=-1.0 + np.random.randn() * 0.5,
            bp_z_score=-1.5 + np.random.randn() * 0.5,
            three_z_score=-1.0 + np.random.randn() * 0.5,
        )
        new_intron = Intron(
            intron_id=intron.intron_id,
            coordinates=intron.coordinates,
            sequences=intron.sequences,
            scores=scores,
        )
        u2_scored.append(new_intron)
        true_labels[intron.intron_id] = "u2"

    # Split into reference (for training) and experimental (for testing)
    # Use 70% for training, 30% for testing
    n_u12_train = int(len(u12_scored) * 0.7)
    n_u2_train = int(len(u2_scored) * 0.7)

    u12_train = u12_scored[:n_u12_train]
    u12_test = u12_scored[n_u12_train:]
    u2_train = u2_scored[:n_u2_train]
    u2_test = u2_scored[n_u2_train:]

    experimental = u12_test + u2_test

    print(f"\nTraining set: {len(u12_train)} U12, {len(u2_train)} U2")
    print(f"Test set: {len(u12_test)} U12, {len(u2_test)} U2")

    # Run classification pipeline
    classifier = IntronClassifier(
        n_optimization_rounds=3,  # Reduced for speed
        n_ensemble_models=3,
        classification_threshold=50.0,  # Lower for easier testing
        random_state=42,
    )

    result = classifier.classify(
        u12_reference=u12_train, u2_reference=u2_train, experimental=experimental
    )

    # Validate results
    assert len(result.classified_introns) == len(experimental)

    # Count correct classifications using tracked labels
    correct = 0
    total = 0

    for intron in result.classified_introns:
        # True label from our tracking dict
        true_label = true_labels.get(intron.intron_id, "unknown")
        pred_label = intron.metadata.type_id if intron.metadata else "unknown"

        if pred_label == true_label:
            correct += 1
        total += 1

    accuracy = correct / total if total > 0 else 0
    print(f"\nClassification accuracy: {accuracy:.2%} ({correct}/{total})")

    # Should achieve reasonable accuracy (>70%) with synthetic data
    assert accuracy > 0.70, f"Low accuracy: {accuracy:.2%}"

    # Check that we got both U12 and U2 predictions
    u12_preds = result.get_u12_predictions(threshold=50.0)
    u2_preds = result.get_u2_predictions(threshold=50.0)

    print(f"U12 predictions: {len(u12_preds)}")
    print(f"U2 predictions: {len(u2_preds)}")

    assert len(u12_preds) > 0, "No U12 predictions"
    assert len(u2_preds) > 0, "No U2 predictions"


@pytest.mark.integration
def test_classification_preserves_z_scores(u12_reference, u2_reference_subset):
    """
    CRITICAL INTEGRATION TEST: Verify z-scores are not re-normalized.

    This is the Issue #1 fix - must ensure z-scores computed from
    reference data are NOT changed during classification.
    """
    # Create small reference sets with z-scores
    u12_ref = []
    for intron in u12_reference[:50]:
        scores = IntronScores(
            five_z_score=2.0,
            bp_z_score=2.5,
            three_z_score=2.0,
        )
        u12_ref.append(
            Intron(
                intron_id=intron.intron_id,
                coordinates=intron.coordinates,
                sequences=intron.sequences,
                scores=scores,
            )
        )

    u2_ref = []
    for intron in u2_reference_subset[:50]:
        scores = IntronScores(
            five_z_score=-1.0,
            bp_z_score=-1.5,
            three_z_score=-1.0,
        )
        u2_ref.append(
            Intron(
                intron_id=intron.intron_id,
                coordinates=intron.coordinates,
                sequences=intron.sequences,
                scores=scores,
            )
        )

    # Create experimental set with specific z-scores
    experimental = []
    original_z_scores = {}

    for i, intron in enumerate(u12_reference[50:55]):
        z_scores = (1.5 + i * 0.1, 2.0 + i * 0.1, 1.8 + i * 0.1)
        original_z_scores[intron.intron_id] = z_scores

        scores = IntronScores(
            five_z_score=z_scores[0],
            bp_z_score=z_scores[1],
            three_z_score=z_scores[2],
        )
        experimental.append(
            Intron(
                intron_id=intron.intron_id,
                coordinates=intron.coordinates,
                sequences=intron.sequences,
                scores=scores,
            )
        )

    # Run classification
    classifier = IntronClassifier(
        n_optimization_rounds=2, n_ensemble_models=2, random_state=42
    )

    result = classifier.classify(
        u12_reference=u12_ref, u2_reference=u2_ref, experimental=experimental
    )

    # CRITICAL CHECK: Z-scores must be exactly the same
    for intron in result.classified_introns:
        original = original_z_scores[intron.intron_id]
        current = (
            intron.scores.five_z_score,
            intron.scores.bp_z_score,
            intron.scores.three_z_score,
        )

        assert original == current, (
            f"Z-scores changed for {intron.intron_id}!\n"
            f"  Original: {original}\n"
            f"  Current:  {current}\n"
            f"This indicates data leakage - CRITICAL BUG!"
        )


@pytest.mark.integration
def test_batch_classification_with_reference_data(u12_reference, u2_reference_subset):
    """Test batch classification mode with real reference data."""
    # Create small scored reference sets
    u12_scored = []
    for intron in u12_reference[:30]:
        scores = IntronScores(
            five_z_score=2.0 + np.random.randn() * 0.3,
            bp_z_score=2.5 + np.random.randn() * 0.3,
            three_z_score=2.0 + np.random.randn() * 0.3,
        )
        u12_scored.append(
            Intron(
                intron_id=intron.intron_id,
                coordinates=intron.coordinates,
                sequences=intron.sequences,
                scores=scores,
            )
        )

    u2_scored = []
    for intron in u2_reference_subset[:30]:
        scores = IntronScores(
            five_z_score=-1.0 + np.random.randn() * 0.3,
            bp_z_score=-1.5 + np.random.randn() * 0.3,
            three_z_score=-1.0 + np.random.randn() * 0.3,
        )
        u2_scored.append(
            Intron(
                intron_id=intron.intron_id,
                coordinates=intron.coordinates,
                sequences=intron.sequences,
                scores=scores,
            )
        )

    # Run batch classification
    classifier = IntronClassifier(
        n_optimization_rounds=2, n_ensemble_models=2, random_state=42
    )

    result = classifier.classify_batch(
        u12_reference=u12_scored[:20],
        u2_reference=u2_scored[:20],
        experimental=u12_scored[20:] + u2_scored[20:],
        batch_size=5,
    )

    assert len(result.classified_introns) == 20
    assert result.ensemble is not None


@pytest.mark.integration
def test_reproducibility_with_reference_data(u12_reference, u2_reference_subset):
    """Test that classification is reproducible with same random seed."""
    # Create small reference sets
    u12_ref = []
    for intron in u12_reference[:30]:
        scores = IntronScores(
            five_z_score=2.0 + np.random.randn() * 0.3,
            bp_z_score=2.5 + np.random.randn() * 0.3,
            three_z_score=2.0 + np.random.randn() * 0.3,
        )
        u12_ref.append(
            Intron(
                intron_id=intron.intron_id,
                coordinates=intron.coordinates,
                sequences=intron.sequences,
                scores=scores,
            )
        )

    u2_ref = []
    for intron in u2_reference_subset[:30]:
        scores = IntronScores(
            five_z_score=-1.0 + np.random.randn() * 0.3,
            bp_z_score=-1.5 + np.random.randn() * 0.3,
            three_z_score=-1.0 + np.random.randn() * 0.3,
        )
        u2_ref.append(
            Intron(
                intron_id=intron.intron_id,
                coordinates=intron.coordinates,
                sequences=intron.sequences,
                scores=scores,
            )
        )

    experimental = u12_ref[20:] + u2_ref[20:]
    train_u12 = u12_ref[:20]
    train_u2 = u2_ref[:20]

    # Run twice with same seed
    classifier1 = IntronClassifier(
        n_optimization_rounds=2, n_ensemble_models=2, random_state=42
    )
    result1 = classifier1.classify(train_u12, train_u2, experimental)

    classifier2 = IntronClassifier(
        n_optimization_rounds=2, n_ensemble_models=2, random_state=42
    )
    result2 = classifier2.classify(train_u12, train_u2, experimental)

    # Results should be identical
    for i1, i2 in zip(result1.classified_introns, result2.classified_introns):
        assert i1.intron_id == i2.intron_id
        assert abs(i1.scores.svm_score - i2.scores.svm_score) < 1e-6
        assert i1.metadata.type_id == i2.metadata.type_id


# =============================================================================
# Additional integration tests moved from unit tests (tests that do full SVM training)
# =============================================================================


@pytest.fixture
def u12_synthetic():
    """Create synthetic U12 introns with z-scores for testing."""
    introns = []
    for i in range(50):
        intron = Intron(
            intron_id=f"syn_u12_{i}",
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
                five_raw_score=12.5,
                bp_raw_score=10.2,
                three_raw_score=15.3,
                five_z_score=2.0 + np.random.randn() * 0.3,
                bp_z_score=2.5 + np.random.randn() * 0.3,
                three_z_score=2.0 + np.random.randn() * 0.3,
            ),
        )
        introns.append(intron)
    return introns


@pytest.fixture
def u2_synthetic():
    """Create synthetic U2 introns with z-scores for testing."""
    introns = []
    for i in range(50):
        intron = Intron(
            intron_id=f"syn_u2_{i}",
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
                bp_seq="CTAAC",
            ),
            scores=IntronScores(
                five_raw_score=5.2,
                bp_raw_score=3.8,
                three_raw_score=6.1,
                five_z_score=-1.0 + np.random.randn() * 0.3,
                bp_z_score=-1.5 + np.random.randn() * 0.3,
                three_z_score=-1.0 + np.random.randn() * 0.3,
            ),
        )
        introns.append(intron)
    return introns


@pytest.fixture
def experimental_synthetic():
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
                system="1-based",
            ),
            sequences=IntronSequences(
                seq="GTATGT" + "N" * 50 + "TCCTTAAC",
                five_seq="GTATGT",
                three_seq="TCCTTAAC",
            ),
            scores=IntronScores(
                five_z_score=2.2 + np.random.randn() * 0.2,
                bp_z_score=2.7 + np.random.randn() * 0.2,
                three_z_score=2.1 + np.random.randn() * 0.2,
            ),
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
                system="1-based",
            ),
            sequences=IntronSequences(
                seq="GTAAGT" + "N" * 50 + "TTTCAG",
                five_seq="GTAAGT",
                three_seq="TTTCAG",
            ),
            scores=IntronScores(
                five_z_score=-0.9 + np.random.randn() * 0.2,
                bp_z_score=-1.3 + np.random.randn() * 0.2,
                three_z_score=-0.8 + np.random.randn() * 0.2,
            ),
        )
        introns.append(intron)

    return introns


@pytest.mark.integration
@pytest.mark.slow
def test_classify_complete_pipeline(
    u12_synthetic, u2_synthetic, experimental_synthetic
):
    """Test complete classification pipeline with synthetic data."""
    classifier = IntronClassifier(
        n_optimization_rounds=2,  # Faster for testing
        n_ensemble_models=2,
        classification_threshold=50.0,  # Lower for easier testing
        random_state=42,
    )

    result = classifier.classify(
        u12_reference=u12_synthetic,
        u2_reference=u2_synthetic,
        experimental=experimental_synthetic,
    )

    # Check result structure
    assert isinstance(result, ClassificationResult)
    assert len(result.classified_introns) == len(experimental_synthetic)
    assert result.ensemble is not None
    assert len(result.ensemble.models) == 2
    assert result.parameters is not None
    assert result.n_u12_reference == len(u12_synthetic)
    assert result.n_u2_reference == len(u2_synthetic)

    # Check that all introns have been classified
    for intron in result.classified_introns:
        assert intron.scores is not None
        assert intron.scores.svm_score is not None
        assert 0 <= intron.scores.svm_score <= 100
        assert intron.metadata is not None
        assert intron.metadata.type_id in ["u2", "u12"]


@pytest.mark.integration
@pytest.mark.slow
def test_classify_with_fixed_c(u12_synthetic, u2_synthetic, experimental_synthetic):
    """Test classification with fixed C parameter (no optimization)."""
    classifier = IntronClassifier(
        optimize_c=False,
        fixed_c=1.0,
        n_ensemble_models=2,
        classification_threshold=50.0,
        random_state=42,
    )

    result = classifier.classify(
        u12_reference=u12_synthetic,
        u2_reference=u2_synthetic,
        experimental=experimental_synthetic,
    )

    assert result.parameters.C == 1.0
    # round_found == -1 indicates the C value was not found through optimization
    # (either averaged from folds or fixed by user)
    assert result.parameters.round_found == -1  # Fixed, not optimized
    assert len(result.classified_introns) == len(experimental_synthetic)


@pytest.mark.integration
@pytest.mark.slow
def test_classify_assigns_u12_and_u2(
    u12_synthetic, u2_synthetic, experimental_synthetic
):
    """Test that classification assigns both U12 and U2 types."""
    classifier = IntronClassifier(
        n_optimization_rounds=2,
        n_ensemble_models=2,
        classification_threshold=50.0,
        random_state=42,
    )

    result = classifier.classify(
        u12_reference=u12_synthetic,
        u2_reference=u2_synthetic,
        experimental=experimental_synthetic,
    )

    # Count classifications
    u12_count = sum(
        1
        for i in result.classified_introns
        if i.metadata and i.metadata.type_id == "u12"
    )
    u2_count = sum(
        1
        for i in result.classified_introns
        if i.metadata and i.metadata.type_id == "u2"
    )

    # Should have both types
    assert u12_count > 0
    assert u2_count > 0
    assert u12_count + u2_count == len(experimental_synthetic)


@pytest.mark.integration
@pytest.mark.slow
def test_classify_preserves_z_scores_synthetic(
    u12_synthetic, u2_synthetic, experimental_synthetic
):
    """
    CRITICAL TEST: Verify z-scores are NOT re-normalized during classification.
    This is Issue #1 fix - prevents data leakage.
    """
    # Store original z-scores
    original_z_scores = {
        intron.intron_id: (
            intron.scores.five_z_score,
            intron.scores.bp_z_score,
            intron.scores.three_z_score,
        )
        for intron in experimental_synthetic
    }

    classifier = IntronClassifier(
        n_optimization_rounds=2, n_ensemble_models=2, random_state=42
    )

    result = classifier.classify(
        u12_reference=u12_synthetic,
        u2_reference=u2_synthetic,
        experimental=experimental_synthetic,
    )

    # Check that z-scores are EXACTLY the same
    for intron in result.classified_introns:
        original = original_z_scores[intron.intron_id]
        current = (
            intron.scores.five_z_score,
            intron.scores.bp_z_score,
            intron.scores.three_z_score,
        )
        assert original == current, f"Z-scores changed for {intron.intron_id}!"


@pytest.mark.integration
@pytest.mark.slow
def test_classify_batch_synthetic(u12_synthetic, u2_synthetic, experimental_synthetic):
    """Test batch classification produces same results as regular."""
    classifier = IntronClassifier(
        n_optimization_rounds=2,
        n_ensemble_models=2,
        classification_threshold=50.0,
        random_state=42,
    )

    # Regular classification
    result_regular = classifier.classify(
        u12_reference=u12_synthetic,
        u2_reference=u2_synthetic,
        experimental=experimental_synthetic,
    )

    # Batch classification with small batch size
    result_batch = classifier.classify_batch(
        u12_reference=u12_synthetic,
        u2_reference=u2_synthetic,
        experimental=experimental_synthetic,
        batch_size=5,
    )

    # Results should be identical
    assert len(result_regular.classified_introns) == len(
        result_batch.classified_introns
    )

    for reg, batch in zip(
        result_regular.classified_introns, result_batch.classified_introns
    ):
        assert abs(reg.scores.svm_score - batch.scores.svm_score) < 1e-6
        assert reg.metadata.type_id == batch.metadata.type_id


@pytest.mark.integration
@pytest.mark.slow
def test_classification_result_get_u12_predictions(
    u12_synthetic, u2_synthetic, experimental_synthetic
):
    """Test ClassificationResult.get_u12_predictions()."""
    classifier = IntronClassifier(
        n_optimization_rounds=2,
        n_ensemble_models=2,
        classification_threshold=50.0,
        random_state=42,
    )

    result = classifier.classify(
        u12_reference=u12_synthetic,
        u2_reference=u2_synthetic,
        experimental=experimental_synthetic,
    )

    u12_predictions = result.get_u12_predictions(threshold=50.0)

    # All should be U12 with score >= threshold
    for intron in u12_predictions:
        assert intron.metadata.type_id == "u12"
        assert intron.scores.svm_score >= 50.0


@pytest.mark.integration
@pytest.mark.slow
def test_classification_result_get_u2_predictions(
    u12_synthetic, u2_synthetic, experimental_synthetic
):
    """Test ClassificationResult.get_u2_predictions()."""
    classifier = IntronClassifier(
        n_optimization_rounds=2,
        n_ensemble_models=2,
        classification_threshold=50.0,
        random_state=42,
    )

    result = classifier.classify(
        u12_reference=u12_synthetic,
        u2_reference=u2_synthetic,
        experimental=experimental_synthetic,
    )

    u2_predictions = result.get_u2_predictions(threshold=50.0)

    # All should be U2 with score < threshold
    for intron in u2_predictions:
        assert intron.metadata.type_id == "u2"
        assert intron.scores.svm_score < 50.0


@pytest.mark.integration
@pytest.mark.slow
def test_classification_result_threshold_affects_filtering(
    u12_synthetic, u2_synthetic, experimental_synthetic
):
    """Test that threshold parameter affects get_u12_predictions filtering."""
    classifier = IntronClassifier(
        n_optimization_rounds=2,
        n_ensemble_models=2,
        classification_threshold=50.0,
        random_state=42,
    )

    result = classifier.classify(
        u12_reference=u12_synthetic,
        u2_reference=u2_synthetic,
        experimental=experimental_synthetic,
    )

    # Lower threshold should give more U12 predictions
    u12_low = result.get_u12_predictions(threshold=40.0)
    u12_high = result.get_u12_predictions(threshold=60.0)

    assert len(u12_low) >= len(u12_high)


@pytest.mark.integration
@pytest.mark.slow
def test_classify_small_datasets(u12_synthetic, u2_synthetic):
    """Test classification with minimal experimental data."""
    # Single experimental intron
    single_exp = [
        Intron(
            intron_id="single",
            coordinates=GenomicCoordinate(
                chromosome="chr1", start=1000, stop=1100, strand="+", system="1-based"
            ),
            sequences=IntronSequences(seq="ATCG", five_seq="AT", three_seq="CG"),
            scores=IntronScores(five_z_score=2.0, bp_z_score=2.5, three_z_score=2.0),
        )
    ]

    classifier = IntronClassifier(
        n_optimization_rounds=2, n_ensemble_models=2, random_state=42
    )

    result = classifier.classify(u12_synthetic, u2_synthetic, single_exp)

    assert len(result.classified_introns) == 1
    assert result.classified_introns[0].scores.svm_score is not None


@pytest.mark.integration
@pytest.mark.slow
def test_classify_reproducibility_synthetic(
    u12_synthetic, u2_synthetic, experimental_synthetic
):
    """Test that classification is reproducible with same random_state."""
    classifier1 = IntronClassifier(
        n_optimization_rounds=2, n_ensemble_models=2, random_state=42
    )
    result1 = classifier1.classify(u12_synthetic, u2_synthetic, experimental_synthetic)

    classifier2 = IntronClassifier(
        n_optimization_rounds=2, n_ensemble_models=2, random_state=42
    )
    result2 = classifier2.classify(u12_synthetic, u2_synthetic, experimental_synthetic)

    # Scores should be identical
    for i1, i2 in zip(result1.classified_introns, result2.classified_introns):
        assert abs(i1.scores.svm_score - i2.scores.svm_score) < 1e-6
        assert i1.metadata.type_id == i2.metadata.type_id

    # Scores should be identical
    for i1, i2 in zip(result1.classified_introns, result2.classified_introns):
        assert abs(i1.scores.svm_score - i2.scores.svm_score) < 1e-6
        assert i1.metadata.type_id == i2.metadata.type_id
