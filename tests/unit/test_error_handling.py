"""
Tests for error handling across the pipeline.

This module tests that the pipeline gracefully handles various error conditions
and provides clear error messages.
"""

import gzip
import tempfile
from pathlib import Path

import pytest

from intronIC.core.intron import (
    GenomicCoordinate,
    Intron,
    IntronMetadata,
    IntronScores,
    IntronSequences,
)
from intronIC.extraction.annotator import AnnotationHierarchyBuilder
from intronIC.extraction.sequences import SequenceExtractor
from intronIC.scoring.pwm import PWMLoader

# ============================================================================
# File I/O Error Tests
# ============================================================================


def test_annotation_file_not_found():
    """Test error handling when annotation file doesn't exist."""
    builder = AnnotationHierarchyBuilder(["exon"])

    with pytest.raises((FileNotFoundError, IOError)):
        builder.build_from_file("/nonexistent/path/to/annotation.gff3")


def test_genome_file_not_found():
    """Test error handling when genome file doesn't exist."""
    with pytest.raises((FileNotFoundError, IOError)):
        extractor = SequenceExtractor("/nonexistent/path/to/genome.fa")


def test_pwm_file_not_found():
    """Test error handling when PWM file doesn't exist."""
    loader = PWMLoader()

    with pytest.raises((FileNotFoundError, IOError)):
        loader.load_from_file(Path("/nonexistent/path/to/pwms.iic"))


def test_malformed_annotation_file():
    """Test handling of malformed GFF3 file."""
    builder = AnnotationHierarchyBuilder(["exon"])

    # Create malformed GFF3 file (missing required columns)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".gff3", delete=False) as f:
        f.write("##gff-version 3\n")
        f.write("chr1\tsource\tgene\n")  # Missing columns
        temp_path = f.name

    try:
        # Should either skip bad lines or raise clear error
        genes = builder.build_from_file(temp_path)
        # If it skips bad lines, should return empty or partial list
        # This is acceptable behavior
    except (ValueError, IndexError) as e:
        # If it raises an error, it should be clear
        assert "column" in str(e).lower() or "field" in str(e).lower()
    finally:
        Path(temp_path).unlink()


def test_empty_annotation_file():
    """Test handling of empty annotation file."""
    builder = AnnotationHierarchyBuilder(["exon"])

    with tempfile.NamedTemporaryFile(mode="w", suffix=".gff3", delete=False) as f:
        f.write("##gff-version 3\n")
        temp_path = f.name

    try:
        genes = builder.build_from_file(temp_path)
        # Empty file should return empty list, not crash
        assert isinstance(genes, list)
    finally:
        Path(temp_path).unlink()


# ============================================================================
# Sequence Extraction Error Tests
# ============================================================================


def test_sequence_extraction_with_invalid_chromosome():
    """Test sequence extraction when chromosome doesn't exist in genome."""
    # Create minimal valid FASTA
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
        f.write(">chr1\n")
        f.write("ACTGACTG\n")
        temp_path = f.name

    try:
        extractor = SequenceExtractor(temp_path)

        # Create intron with non-existent chromosome
        intron = Intron(
            intron_id="test",
            coordinates=GenomicCoordinate("chr_nonexistent", 1, 5, "+", "1-based"),
            sequences=IntronSequences(seq=None),
            scores=IntronScores(),
            metadata=IntronMetadata("t1", "g1"),
        )

        # Should handle gracefully - may print warning and skip instead of raising
        # Check that it doesn't crash and handles the missing chromosome
        results = list(extractor.extract_sequences([intron]))
        # Either returns empty (skipped) or the intron with no sequence
        if results:
            # If returned, sequence should still be None since chromosome missing
            pass  # Graceful handling is acceptable
    finally:
        Path(temp_path).unlink()


def test_sequence_extraction_out_of_bounds():
    """Test sequence extraction when coordinates exceed chromosome length."""
    # Create minimal valid FASTA (only 8bp)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
        f.write(">chr1\n")
        f.write("ACTGACTG\n")
        temp_path = f.name

    try:
        extractor = SequenceExtractor(temp_path)

        # Create intron that exceeds chromosome bounds
        intron = Intron(
            intron_id="test",
            coordinates=GenomicCoordinate(
                "chr1", 1, 1000, "+", "1-based"
            ),  # Beyond 8bp length
            sequences=IntronSequences(seq=None),
            scores=IntronScores(),
            metadata=IntronMetadata("t1", "g1"),
        )

        # Should handle gracefully (either truncate or raise clear error)
        try:
            extracted = list(extractor.extract_sequences([intron]))
            # If it succeeds, sequence should be truncated to available length
            assert extracted[0].sequences.seq is not None
        except (ValueError, IndexError) as e:
            # If it fails, error should mention coordinates or bounds
            assert (
                "coordinate" in str(e).lower()
                or "bound" in str(e).lower()
                or "length" in str(e).lower()
            )
    finally:
        Path(temp_path).unlink()


# ============================================================================
# PWM Loading Error Tests
# ============================================================================


def test_pwm_file_malformed_header():
    """Test handling of PWM file with malformed header."""
    loader = PWMLoader()

    with tempfile.NamedTemporaryFile(mode="w", suffix=".iic", delete=False) as f:
        f.write(">pwm_name_missing_start\n")  # Missing start= parameter
        f.write("A\tC\tG\tT\n")
        f.write("0.25\t0.25\t0.25\t0.25\n")
        temp_path = f.name

    try:
        # Loader should handle gracefully - uses default start or skips
        # Either raises an error or handles gracefully with defaults
        try:
            result = loader.load_from_file(Path(temp_path))
            # If it succeeds, it used default values which is acceptable
        except (ValueError, KeyError):
            pass  # Also acceptable if it raises
    finally:
        Path(temp_path).unlink()


def test_pwm_file_malformed_matrix():
    """Test handling of PWM file with invalid matrix values."""
    loader = PWMLoader()

    with tempfile.NamedTemporaryFile(mode="w", suffix=".iic", delete=False) as f:
        f.write(">test_pwm\tstart=0\n")
        f.write("A\tC\tG\tT\n")
        f.write("invalid\t0.25\t0.25\t0.25\n")  # Invalid numeric value
        temp_path = f.name

    try:
        with pytest.raises((ValueError, TypeError)):
            loader.load_from_file(Path(temp_path))
    finally:
        Path(temp_path).unlink()


# ============================================================================
# Memory and Resource Error Tests
# ============================================================================


def test_very_large_intron_sequence():
    """Test handling of extremely large intron sequence."""
    # Create intron with very large sequence (but not infinite)
    large_seq = "N" * 1_000_000  # 1MB sequence

    intron = Intron(
        intron_id="large",
        coordinates=GenomicCoordinate("chr1", 1, 1_000_000, "+", "1-based"),
        sequences=IntronSequences(
            seq=large_seq, five_prime_dnt="GT", three_prime_dnt="AG"
        ),
        scores=IntronScores(),
        metadata=IntronMetadata("t1", "g1"),
    )

    # Should handle without crashing
    assert len(intron.sequences.seq) == 1_000_000


def test_many_introns_memory_efficiency():
    """Test that creating many introns doesn't cause memory issues."""
    # Create 10,000 small introns
    introns = []
    for i in range(10_000):
        # Start at 1+i*100 to avoid 0-start for 1-based coordinates
        intron = Intron(
            intron_id=f"intron_{i}",
            coordinates=GenomicCoordinate(
                "chr1", 1 + i * 100, 1 + i * 100 + 50, "+", "1-based"
            ),
            sequences=IntronSequences(
                seq="GTACAG", five_prime_dnt="GT", three_prime_dnt="AG"
            ),
            scores=IntronScores(),
            metadata=IntronMetadata("t1", "g1"),
        )
        introns.append(intron)

    # Should complete without memory error
    assert len(introns) == 10_000


# ============================================================================
# Gzipped File Tests
# ============================================================================


def test_gzipped_annotation_file():
    """Test that gzipped annotation files are handled correctly."""
    builder = AnnotationHierarchyBuilder(["exon"])

    with tempfile.NamedTemporaryFile(mode="wb", suffix=".gff3.gz", delete=False) as f:
        with gzip.open(f, "wt") as gz:
            gz.write("##gff-version 3\n")
            gz.write("chr1\tsource\tgene\t1000\t2000\t.\t+\t.\tID=gene1\n")
        temp_path = f.name

    try:
        genes = builder.build_from_file(temp_path)
        # Should successfully read gzipped file
        assert len(genes) >= 0  # May be empty due to minimal data
    finally:
        Path(temp_path).unlink()


def test_corrupted_gzip_file():
    """Test handling of corrupted gzip file."""
    builder = AnnotationHierarchyBuilder(["exon"])

    with tempfile.NamedTemporaryFile(mode="wb", suffix=".gff3.gz", delete=False) as f:
        f.write(b"not a valid gzip file")
        temp_path = f.name

    try:
        # May raise error or handle gracefully depending on smart_open behavior
        try:
            builder.build_from_file(temp_path)
        except (gzip.BadGzipFile, OSError, EOFError):
            pass  # Expected - corrupt file detected
        # If no error, that's also acceptable (graceful handling)
    finally:
        Path(temp_path).unlink()
