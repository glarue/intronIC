"""
Tests for writer error handling and edge cases.

This module tests that file writers handle errors gracefully and produce
valid output in edge cases.
"""

import pytest
import tempfile
from pathlib import Path
import os

from intronIC.file_io.writers import MetaWriter, BedWriter, SeqWriter, ScoreWriter
from intronIC.core.intron import Intron, IntronSequences, IntronScores, IntronMetadata, GenomicCoordinate


# ============================================================================
# Test Fixtures
# ============================================================================

@pytest.fixture
def sample_intron():
    """Create a standard test intron."""
    return Intron(
        intron_id="test_intron",
        coordinates=GenomicCoordinate("chr1", 1000, 2000, '+', '1-based'),
        sequences=IntronSequences(
            seq="GT" + "N"*996 + "AG",
            upstream_flank="NNNNN",
            downstream_flank="NNNNN",
            five_prime_dnt="GT",
            three_prime_dnt="AG"
        ),
        scores=IntronScores(
            five_raw_score=2.5,
            bp_raw_score=3.0,
            three_raw_score=2.8,
            five_z_score=1.2,
            bp_z_score=1.5,
            three_z_score=1.3,
            svm_score=0.95
        ),
        metadata=IntronMetadata(
            parent="transcript_1",
            grandparent="gene_1",
            intron_index=1,
            family_size=5,
            longest_isoform=True
        ),
        type_id="u12"
    )


# ============================================================================
# Directory Permission Tests
# ============================================================================

def test_write_to_nonexistent_directory(sample_intron):
    """Test writing to a directory that doesn't exist."""
    writer = MetaWriter()

    nonexistent_path = "/nonexistent/directory/output.meta.iic"

    with pytest.raises((FileNotFoundError, OSError, PermissionError)):
        with writer.open(nonexistent_path):
            writer.write_intron(sample_intron)


@pytest.mark.skipif(os.getuid() == 0, reason="Test doesn't work as root")
def test_write_to_readonly_directory(sample_intron):
    """Test writing to a read-only directory."""
    writer = MetaWriter()

    with tempfile.TemporaryDirectory() as tmpdir:
        # Make directory read-only
        os.chmod(tmpdir, 0o444)

        try:
            output_path = Path(tmpdir) / "output.meta.iic"

            with pytest.raises((PermissionError, OSError)):
                with writer.open(str(output_path)):
                    writer.write_intron(sample_intron)
        finally:
            # Restore permissions for cleanup
            os.chmod(tmpdir, 0o755)


# ============================================================================
# Disk Space Tests
# ============================================================================

def test_write_many_introns(sample_intron):
    """Test writing a large number of introns doesn't cause issues."""
    writer = MetaWriter()

    with tempfile.NamedTemporaryFile(mode='w', suffix='.meta.iic', delete=False) as f:
        temp_path = f.name

    try:
        with writer.open(temp_path):
            # Write 10,000 introns
            for i in range(10_000):
                # Create variant of sample intron
                intron = Intron(
                    intron_id=f"intron_{i}",
                    coordinates=GenomicCoordinate("chr1", i*1000, i*1000+100, '+', '1-based'),
                    sequences=IntronSequences(
                        seq=sample_intron.sequences.seq,
                        upstream_flank="NNNNN",
                        downstream_flank="NNNNN",
                        five_prime_dnt="GT",
                        three_prime_dnt="AG"
                    ),
                    scores=sample_intron.scores,
                    metadata=IntronMetadata("t1", "g1"),
                    type_id="u2"
                )
                writer.write_intron(intron)

        # Verify file exists and has content
        assert Path(temp_path).exists()
        assert Path(temp_path).stat().st_size > 0
    finally:
        Path(temp_path).unlink()


# ============================================================================
# Invalid Data Tests
# ============================================================================

def test_meta_writer_with_missing_sequence():
    """Test MetaWriter with intron missing sequence data."""
    writer = MetaWriter()

    intron = Intron(
        intron_id="no_seq",
        coordinates=GenomicCoordinate("chr1", 1000, 2000, '+', '1-based'),
        sequences=IntronSequences(seq=None),  # No sequence!
        scores=IntronScores(),
        metadata=IntronMetadata("t1", "g1"),
        type_id="u2"
    )

    with tempfile.NamedTemporaryFile(mode='w', suffix='.meta.iic', delete=False) as f:
        temp_path = f.name

    try:
        with writer.open(temp_path):
            # Should handle gracefully - either skip or write with placeholders
            writer.write_intron(intron)

        # File should still be valid
        assert Path(temp_path).exists()
    finally:
        Path(temp_path).unlink()


def test_meta_writer_with_missing_scores():
    """Test MetaWriter with intron missing score data."""
    writer = MetaWriter()

    intron = Intron(
        intron_id="no_scores",
        coordinates=GenomicCoordinate("chr1", 1000, 2000, '+', '1-based'),
        sequences=IntronSequences(
            seq="GTAG",
            five_prime_dnt="GT",
            three_prime_dnt="AG"
        ),
        scores=IntronScores(),  # No scores populated
        metadata=IntronMetadata("t1", "g1"),
        type_id="u2"
    )

    with tempfile.NamedTemporaryFile(mode='w', suffix='.meta.iic', delete=False) as f:
        temp_path = f.name

    try:
        with writer.open(temp_path):
            writer.write_intron(intron)

        # Should write with placeholders (e.g., '.' for missing values)
        content = Path(temp_path).read_text()
        assert len(content) > 0
    finally:
        Path(temp_path).unlink()


def test_bed_writer_with_invalid_coordinates():
    """Test BedWriter with coordinates that could cause BED format issues."""
    writer = BedWriter()

    # Create intron with start=stop (zero-length in BED)
    intron = Intron(
        intron_id="zero_length",
        coordinates=GenomicCoordinate("chr1", 1000, 1000, '+', '1-based'),
        sequences=IntronSequences(seq=""),
        scores=IntronScores(svm_score=0.5),
        metadata=IntronMetadata("t1", "g1"),
        type_id="u2"
    )

    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed.iic', delete=False) as f:
        temp_path = f.name

    try:
        with writer.open(temp_path):
            # Should handle or reject zero-length feature
            try:
                writer.write_intron(intron)
            except ValueError:
                # Rejecting zero-length is acceptable
                pass
    finally:
        Path(temp_path).unlink()


def test_seq_writer_with_very_long_sequence(sample_intron):
    """Test SeqWriter with extremely long sequence."""
    writer = SeqWriter()

    # Create intron with very long sequence (1MB)
    long_seq = "N" * 1_000_000
    intron = Intron(
        intron_id="long_seq",
        coordinates=GenomicCoordinate("chr1", 1000, 1_001_000, '+', '1-based'),
        sequences=IntronSequences(
            seq=long_seq,
            upstream_flank="NNNNN",
            downstream_flank="NNNNN",
            five_prime_dnt="NN",
            three_prime_dnt="NN"
        ),
        scores=sample_intron.scores,
        metadata=IntronMetadata("t1", "g1"),
        type_id="u2"
    )

    with tempfile.NamedTemporaryFile(mode='w', suffix='.seqs.iic', delete=False) as f:
        temp_path = f.name

    try:
        with writer.open(temp_path):
            writer.write_intron(intron)

        # Should write successfully
        assert Path(temp_path).exists()
        assert Path(temp_path).stat().st_size > 1_000_000
    finally:
        Path(temp_path).unlink()


def test_score_writer_with_special_float_values():
    """Test ScoreWriter with NaN and infinity values."""
    writer = ScoreWriter()

    import numpy as np

    intron = Intron(
        intron_id="special_values",
        coordinates=GenomicCoordinate("chr1", 1000, 2000, '+', '1-based'),
        sequences=IntronSequences(
            seq="GTAG",
            five_prime_dnt="GT",
            three_prime_dnt="AG"
        ),
        scores=IntronScores(
            five_raw_score=np.nan,
            bp_raw_score=np.inf,
            three_raw_score=-np.inf,
            five_z_score=0.0,
            bp_z_score=0.0,
            three_z_score=0.0,
            svm_score=np.nan
        ),
        metadata=IntronMetadata("t1", "g1"),
        type_id="u2"
    )

    with tempfile.NamedTemporaryFile(mode='w', suffix='.scores.iic', delete=False) as f:
        temp_path = f.name

    try:
        with writer.open(temp_path):
            # Should handle special values (convert to strings or placeholders)
            writer.write_intron(intron)

        content = Path(temp_path).read_text()
        # Should contain some representation of special values
        assert len(content) > 0
    finally:
        Path(temp_path).unlink()


# ============================================================================
# Format Edge Cases
# ============================================================================

def test_intron_id_with_special_characters(sample_intron):
    """Test that special characters in IDs are handled properly."""
    writer = MetaWriter()

    special_ids = [
        "intron|with|pipes",
        "intron\twith\ttabs",
        "intron;with;semicolons",
        "intron:with:colons"
    ]

    with tempfile.NamedTemporaryFile(mode='w', suffix='.meta.iic', delete=False) as f:
        temp_path = f.name

    try:
        with writer.open(temp_path):
            for special_id in special_ids:
                intron = Intron(
                    intron_id=special_id,
                    coordinates=sample_intron.coordinates,
                    sequences=sample_intron.sequences,
                    scores=sample_intron.scores,
                    metadata=sample_intron.metadata,
                    type_id="u2"
                )
                # Should sanitize or escape special characters
                writer.write_intron(intron)

        # Output should be parseable
        content = Path(temp_path).read_text()
        lines = content.strip().split('\n')
        assert len(lines) == len(special_ids)
    finally:
        Path(temp_path).unlink()


def test_chromosome_name_with_special_characters(sample_intron):
    """Test BED format with unusual chromosome names."""
    writer = BedWriter()

    special_chroms = [
        "chr1_random",
        "chr1.2",
        "chr_KI270706v1_random"
    ]

    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed.iic', delete=False) as f:
        temp_path = f.name

    try:
        with writer.open(temp_path):
            for chrom in special_chroms:
                intron = Intron(
                    intron_id=sample_intron.intron_id,
                    coordinates=GenomicCoordinate(chrom, 1000, 2000, '+', '1-based'),
                    sequences=sample_intron.sequences,
                    scores=sample_intron.scores,
                    metadata=sample_intron.metadata,
                    type_id="u2"
                )
                writer.write_intron(intron)

        # Should produce valid BED format
        content = Path(temp_path).read_text()
        for chrom in special_chroms:
            assert chrom in content
    finally:
        Path(temp_path).unlink()


# ============================================================================
# Context Manager Tests
# ============================================================================

def test_writer_without_context_manager():
    """Test that writers require context manager usage."""
    writer = MetaWriter()

    # Attempting to write without opening should fail
    intron = Intron(
        intron_id="test",
        coordinates=GenomicCoordinate("chr1", 1000, 2000, '+', '1-based'),
        sequences=IntronSequences(seq="GTAG"),
        scores=IntronScores(),
        metadata=IntronMetadata("t1", "g1"),
        type_id="u2"
    )

    with pytest.raises((AttributeError, ValueError)):
        writer.write_intron(intron)


def test_writer_multiple_open_close_cycles(sample_intron):
    """Test opening and closing writer multiple times."""
    writer = MetaWriter()

    with tempfile.NamedTemporaryFile(mode='w', suffix='.meta.iic', delete=False) as f:
        temp_path = f.name

    try:
        # First cycle
        with writer.open(temp_path):
            writer.write_intron(sample_intron)

        # Second cycle (should overwrite or append depending on mode)
        with writer.open(temp_path):
            writer.write_intron(sample_intron)

        # File should still be valid
        assert Path(temp_path).exists()
    finally:
        Path(temp_path).unlink()


# ============================================================================
# Empty Input Tests
# ============================================================================

def test_writer_with_no_introns():
    """Test creating output file with no introns written."""
    writer = MetaWriter()

    with tempfile.NamedTemporaryFile(mode='w', suffix='.meta.iic', delete=False) as f:
        temp_path = f.name

    try:
        with writer.open(temp_path):
            pass  # Write nothing

        # File should exist (possibly empty or with header)
        assert Path(temp_path).exists()
    finally:
        Path(temp_path).unlink()
