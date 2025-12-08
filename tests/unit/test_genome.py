"""
Unit tests for genome reading functionality.

Tests:
- parse_fasta: Streaming FASTA parser
- GenomeReader: Genome file reading (streaming and cached modes)
- Subsequence extraction with coordinates

Author: intronIC refactoring project
Date: 2025-11-02
"""

import pytest
from pathlib import Path
from intronIC.file_io.genome import parse_fasta, GenomeReader
from intronIC.utils.coordinates import GenomicCoordinate


# Test data paths - use fixtures in tests


class TestParseFasta:
    """Test parse_fasta streaming parser."""

    def test_parse_simple_fasta(self, tmp_path):
        """Test parsing a simple FASTA file."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(
            ">seq1\n"
            "ATCGATCG\n"
            ">seq2\n"
            "GGCCGGCC\n"
        )

        sequences = list(parse_fasta(fasta_file))
        assert len(sequences) == 2
        assert sequences[0] == ("seq1", "ATCGATCG")
        assert sequences[1] == ("seq2", "GGCCGGCC")

    def test_parse_multiline_sequences(self, tmp_path):
        """Test parsing FASTA with multi-line sequences."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(
            ">seq1\n"
            "ATCG\n"
            "ATCG\n"
            ">seq2\n"
            "GGCC\n"
            "GGCC\n"
        )

        sequences = list(parse_fasta(fasta_file))
        assert sequences[0] == ("seq1", "ATCGATCG")
        assert sequences[1] == ("seq2", "GGCCGGCC")

    def test_parse_with_empty_lines(self, tmp_path):
        """Test parsing FASTA with empty lines."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(
            ">seq1\n"
            "ATCG\n"
            "\n"
            "ATCG\n"
            "\n"
            ">seq2\n"
            "GGCC\n"
        )

        sequences = list(parse_fasta(fasta_file))
        assert sequences[0] == ("seq1", "ATCGATCG")
        assert sequences[1] == ("seq2", "GGCC")

    def test_parse_uppercases_sequences(self, tmp_path):
        """Test that sequences are converted to uppercase."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(">seq1\natcg\n")

        sequences = list(parse_fasta(fasta_file))
        assert sequences[0] == ("seq1", "ATCG")

    def test_parse_header_first_word_only(self, tmp_path):
        """Test that only first word of header is used as name."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(">seq1 additional info here\nATCG\n")

        sequences = list(parse_fasta(fasta_file))
        assert sequences[0][0] == "seq1"


class TestGenomeReaderStreaming:
    """Test GenomeReader in streaming mode."""

    def test_create_reader_streaming(self, tmp_path):
        """Test creating reader in streaming mode."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(">chr1\nATCG\n")

        reader = GenomeReader(fasta_file, cached=False)
        assert reader.is_cached == False
        assert reader.cache is None

    def test_stream_sequences(self, tmp_path):
        """Test streaming sequences from file."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(
            ">chr1\nATCG\n"
            ">chr2\nGGCC\n"
        )

        reader = GenomeReader(fasta_file, cached=False)
        sequences = list(reader.stream())

        assert len(sequences) == 2
        assert sequences[0] == ("chr1", "ATCG")
        assert sequences[1] == ("chr2", "GGCC")

    def test_get_sequence_fails_in_streaming_mode(self, tmp_path):
        """Test that get_sequence raises error in streaming mode."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(">chr1\nATCG\n")

        reader = GenomeReader(fasta_file, cached=False)

        with pytest.raises(RuntimeError, match="requires cached mode"):
            reader.get_sequence("chr1")


class TestGenomeReaderCached:
    """Test GenomeReader in cached mode."""

    def test_create_reader_cached(self, tmp_path):
        """Test creating reader in cached mode."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(">chr1\nATCG\n")

        reader = GenomeReader(fasta_file, cached=True)
        assert reader.is_cached == True
        assert reader.cache is not None
        assert "chr1" in reader.cache

    def test_get_sequence(self, tmp_path):
        """Test getting sequence by chromosome name."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(
            ">chr1\nATCG\n"
            ">chr2\nGGCC\n"
        )

        reader = GenomeReader(fasta_file, cached=True)
        assert reader.get_sequence("chr1") == "ATCG"
        assert reader.get_sequence("chr2") == "GGCC"

    def test_get_sequence_not_found(self, tmp_path):
        """Test KeyError when chromosome not found."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(">chr1\nATCG\n")

        reader = GenomeReader(fasta_file, cached=True)

        with pytest.raises(KeyError, match="chr999"):
            reader.get_sequence("chr999")

    def test_has_chromosome(self, tmp_path):
        """Test checking if chromosome exists."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(">chr1\nATCG\n")

        reader = GenomeReader(fasta_file, cached=True)
        assert reader.has_chromosome("chr1") == True
        assert reader.has_chromosome("chr999") == False

    def test_get_chromosome_names(self, tmp_path):
        """Test getting list of chromosome names."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(
            ">chr1\nATCG\n"
            ">chr2\nGGCC\n"
        )

        reader = GenomeReader(fasta_file, cached=True)
        names = reader.get_chromosome_names()
        assert "chr1" in names
        assert "chr2" in names
        assert len(names) == 2

    def test_get_chromosome_length(self, tmp_path):
        """Test getting chromosome length."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(">chr1\nATCGATCG\n")

        reader = GenomeReader(fasta_file, cached=True)
        assert reader.get_chromosome_length("chr1") == 8


class TestSubsequenceExtraction:
    """Test extracting subsequences using coordinates."""

    def test_extract_simple_subsequence(self, tmp_path):
        """Test extracting a simple subsequence."""
        fasta_file = tmp_path / "test.fa"
        # Create a 20bp sequence
        fasta_file.write_text(">chr1\nATCGATCGATCGATCGATCG\n")

        reader = GenomeReader(fasta_file, cached=True)
        coord = GenomicCoordinate("chr1", 1, 4, '+', '1-based')

        # Extract positions 1-4 (1-based, inclusive) = ATCG
        subseq = reader.extract_subsequence(coord)
        assert subseq == "ATCG"

    def test_extract_with_upstream_flank(self, tmp_path):
        """Test extracting with upstream flanking sequence."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(">chr1\nATCGATCGATCGATCGATCG\n")

        reader = GenomeReader(fasta_file, cached=True)
        coord = GenomicCoordinate("chr1", 5, 8, '+', '1-based')

        # Extract 5-8 with 2bp upstream flank
        # Positions 3-8 (1-based) = CGATCG
        subseq = reader.extract_subsequence(coord, upstream_flank=2)
        assert subseq == "CGATCG"

    def test_extract_with_downstream_flank(self, tmp_path):
        """Test extracting with downstream flanking sequence."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(">chr1\nATCGATCGATCGATCGATCG\n")

        reader = GenomeReader(fasta_file, cached=True)
        coord = GenomicCoordinate("chr1", 1, 4, '+', '1-based')

        # Extract 1-4 with 2bp downstream flank
        # Positions 1-6 (1-based) = ATCGAT
        subseq = reader.extract_subsequence(coord, downstream_flank=2)
        assert subseq == "ATCGAT"

    def test_extract_negative_strand_reverse_complements(self, tmp_path):
        """Test that negative strand sequences are reverse complemented."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(">chr1\nATCGATCGATCGATCGATCG\n")

        reader = GenomeReader(fasta_file, cached=True)
        coord = GenomicCoordinate("chr1", 1, 4, '-', '1-based')

        # Extract 1-4 on negative strand
        # ATCG -> reverse complement -> CGAT
        subseq = reader.extract_subsequence(coord)
        assert subseq == "CGAT"

    def test_extract_out_of_bounds_raises_error(self, tmp_path):
        """Test that out-of-bounds coordinates raise ValueError."""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(">chr1\nATCGATCG\n")

        reader = GenomeReader(fasta_file, cached=True)
        coord = GenomicCoordinate("chr1", 1, 100, '+', '1-based')

        with pytest.raises(ValueError, match="out of bounds"):
            reader.extract_subsequence(coord)


class TestWithRealData:
    """Test with real chr19 genome data."""

    def test_load_chr19_cached(self, test_data_dir):
        """Test loading chr19 genome in cached mode."""
        chr19_genome = test_data_dir / "Homo_sapiens.Chr19.Ensembl_91.fa.gz"
        if not chr19_genome.exists():
            pytest.skip("chr19 test data not available")

        reader = GenomeReader(chr19_genome, cached=True)

        assert reader.is_cached
        assert "19" in reader.get_chromosome_names()

    def test_chr19_sequence_length(self, test_data_dir):
        """Test chr19 sequence length is reasonable."""
        chr19_genome = test_data_dir / "Homo_sapiens.Chr19.Ensembl_91.fa.gz"
        if not chr19_genome.exists():
            pytest.skip("chr19 test data not available")

        reader = GenomeReader(chr19_genome, cached=True)

        chr19_length = reader.get_chromosome_length("19")
        # Chr19 is ~58M bp, allow some variation
        assert 50_000_000 < chr19_length < 70_000_000

    def test_extract_from_chr19(self, test_data_dir):
        """Test extracting a subsequence from chr19."""
        chr19_genome = test_data_dir / "Homo_sapiens.Chr19.Ensembl_91.fa.gz"
        if not chr19_genome.exists():
            pytest.skip("chr19 test data not available")

        reader = GenomeReader(chr19_genome, cached=True)

        # Extract first 100bp
        coord = GenomicCoordinate("19", 1, 100, '+', '1-based')
        subseq = reader.extract_subsequence(coord)

        assert len(subseq) == 100
        assert all(base in 'ATCGN' for base in subseq)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
