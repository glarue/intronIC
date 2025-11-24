"""
Unit tests for coordinate handling and conversion.

This is CRITICAL for intronIC because coordinate bugs are the #1 risk area.

Coordinate Systems:
- BED format: 0-based, half-open [start, stop)
- GFF3/GTF format: 1-based, closed [start, stop]
- intronIC internal: 1-based, closed [start, stop]

Example:
  First 10 bases of a chromosome:
  - BED:  [0, 10)  means bases 1-10 (0-based start, exclusive stop)
  - GFF3: [1, 10]  means bases 1-10 (1-based start, inclusive stop)
  - Both represent the same genomic region!

Test Coverage:
1. Basic coordinate conversions (0-based ↔ 1-based)
2. BED format parsing and conversion
3. BED output generation
4. Strand handling (especially reverse strand)
5. Coordinate validation (start < stop, valid strand)
6. Edge cases (chromosome boundaries, single-base features)
"""

import pytest

# Import production implementations
from intronIC.utils.coordinates import (
    GenomicCoordinate,
    bed_to_internal,
    internal_to_bed,
    gff_to_internal,
)


# =============================================================================
# Test Suite
# =============================================================================

class TestCoordinateConversion:
    """Test basic coordinate system conversions."""

    def test_bed_to_internal_simple(self):
        """Test BED → internal conversion for a simple case."""
        # BED: [0, 10) represents bases 1-10
        coord = bed_to_internal("chr1", 0, 10, '+')
        assert coord.chromosome == "chr1"
        assert coord.start == 1
        assert coord.stop == 10
        assert coord.strand == '+'
        assert coord.system == '1-based'

    def test_internal_to_bed_simple(self):
        """Test internal → BED conversion for a simple case."""
        # Internal: [1, 10] represents bases 1-10
        coord = GenomicCoordinate("chr1", 1, 10, '+', '1-based')
        chrom, start, stop, strand = internal_to_bed(coord)
        assert chrom == "chr1"
        assert start == 0
        assert stop == 10
        assert strand == '+'

    def test_round_trip_bed_conversion(self):
        """Test BED → internal → BED preserves coordinates."""
        # Start with BED coordinates
        original_bed = ("chr1", 1000, 2000, '+')

        # Convert to internal
        coord = bed_to_internal(*original_bed)

        # Convert back to BED
        result_bed = internal_to_bed(coord)

        assert result_bed == original_bed

    def test_gff_to_internal_is_noop(self):
        """Test that GFF → internal is a no-op (both 1-based)."""
        coord = gff_to_internal("chr1", 1, 10, '+')
        assert coord.start == 1
        assert coord.stop == 10
        assert coord.system == '1-based'

    def test_intron_coordinate_example(self):
        """
        Test realistic intron coordinates.

        Example: Intron between exon 1 (bases 100-200) and exon 2 (bases 300-400)
        Intron spans bases 201-299 (1-based, 99 bases)
        """
        # In BED format: [200, 299)
        coord = bed_to_internal("chr1", 200, 299, '+')
        assert coord.start == 201  # First base of intron
        assert coord.stop == 299   # Last base of intron
        assert coord.stop - coord.start + 1 == 99  # Length


class TestBEDParsing:
    """Test BED file parsing with coordinate conversion."""

    def test_parse_bed_line_positive_strand(self):
        """Test parsing a BED line on positive strand."""
        bed_line = "chr1\t1000\t2000\tintron_1\t50\t+"
        fields = bed_line.split('\t')

        coord = bed_to_internal(fields[0], int(fields[1]), int(fields[2]), fields[5])

        assert coord.chromosome == "chr1"
        assert coord.start == 1001  # BED 1000 → internal 1001
        assert coord.stop == 2000
        assert coord.strand == '+'

    def test_parse_bed_line_negative_strand(self):
        """Test parsing a BED line on negative strand."""
        bed_line = "chr2\t5000\t6000\tintron_2\t75\t-"
        fields = bed_line.split('\t')

        coord = bed_to_internal(fields[0], int(fields[1]), int(fields[2]), fields[5])

        assert coord.chromosome == "chr2"
        assert coord.start == 5001
        assert coord.stop == 6000
        assert coord.strand == '-'

    def test_bed_single_base_feature(self):
        """
        Test single-base feature in BED format.

        BED: [100, 101) represents a single base at position 101 (1-based)

        Note: Single-base features (start==stop) are now allowed to match
        original intronIC behavior, which extracts very short (<30bp) introns
        and marks them as omitted during filtering.
        """
        # BED [100, 101) converts to internal [101, 101]
        # This is now accepted (1bp feature)
        coord = bed_to_internal("chr1", 100, 101, '+')
        assert coord.start == 101
        assert coord.stop == 101
        assert coord.length == 1


class TestBEDOutput:
    """Test BED output generation with coordinate conversion."""

    def test_generate_bed_output_positive_strand(self):
        """Test generating BED output from internal coordinates."""
        coord = GenomicCoordinate("chr1", 1001, 2000, '+', '1-based')
        chrom, start, stop, strand = internal_to_bed(coord)

        bed_line = f"{chrom}\t{start}\t{stop}\tintron_1\t50\t{strand}"
        assert bed_line == "chr1\t1000\t2000\tintron_1\t50\t+"

    def test_generate_bed_output_negative_strand(self):
        """Test generating BED output on negative strand."""
        coord = GenomicCoordinate("chr2", 5001, 6000, '-', '1-based')
        chrom, start, stop, strand = internal_to_bed(coord)

        bed_line = f"{chrom}\t{start}\t{stop}\tintron_2\t75\t{strand}"
        assert bed_line == "chr2\t5000\t6000\tintron_2\t75\t-"


class TestStrandHandling:
    """Test coordinate handling with strand considerations."""

    def test_positive_strand_coordinate_order(self):
        """On positive strand, start < stop always."""
        coord = GenomicCoordinate("chr1", 100, 200, '+')
        assert coord.start < coord.stop

    def test_negative_strand_coordinate_order(self):
        """On negative strand, start < stop still (genomic coordinates)."""
        coord = GenomicCoordinate("chr1", 100, 200, '-')
        assert coord.start < coord.stop
        # Note: Even on negative strand, genomic coordinates have start < stop
        # The "biological" direction is reversed, but coordinates are not

    def test_strand_validation_positive(self):
        """Test that '+' strand is accepted."""
        coord = GenomicCoordinate("chr1", 100, 200, '+')
        assert coord.strand == '+'

    def test_strand_validation_negative(self):
        """Test that '-' strand is accepted."""
        coord = GenomicCoordinate("chr1", 100, 200, '-')
        assert coord.strand == '-'

    def test_invalid_strand_rejected(self):
        """Test that invalid strand symbols are rejected."""
        with pytest.raises(ValueError, match="strand must be"):
            GenomicCoordinate("chr1", 100, 200, '.')

        with pytest.raises(ValueError, match="strand must be"):
            GenomicCoordinate("chr1", 100, 200, '?')


class TestCoordinateValidation:
    """Test coordinate validation rules."""

    def test_start_less_than_stop_valid(self):
        """Test that start < stop is accepted."""
        coord = GenomicCoordinate("chr1", 100, 200, '+')
        assert coord.start == 100
        assert coord.stop == 200

    def test_start_equals_stop_valid(self):
        """Test that start == stop is now allowed (1bp features)."""
        coord = GenomicCoordinate("chr1", 100, 100, '+')
        assert coord.start == 100
        assert coord.stop == 100
        assert coord.length == 1

    def test_start_greater_than_stop_invalid(self):
        """Test that start > stop is rejected."""
        with pytest.raises(ValueError, match="start .* must be <= stop"):
            GenomicCoordinate("chr1", 200, 100, '+')

    def test_negative_start_invalid(self):
        """Test that negative start coordinate is rejected for 1-based system."""
        with pytest.raises(ValueError, match="start coordinate must be >= 1"):
            GenomicCoordinate("chr1", -1, 100, '+')

    def test_zero_start_valid_for_bed(self):
        """Test that 0-based coordinates can have start=0."""
        coord = GenomicCoordinate("chr1", 0, 100, '+', '0-based')
        assert coord.start == 0

    def test_chromosome_name_formats(self):
        """Test various chromosome name formats are accepted."""
        # Standard formats
        coord1 = GenomicCoordinate("chr1", 1, 100, '+')
        assert coord1.chromosome == "chr1"

        coord2 = GenomicCoordinate("1", 1, 100, '+')
        assert coord2.chromosome == "1"

        coord3 = GenomicCoordinate("chrX", 1, 100, '+')
        assert coord3.chromosome == "chrX"

        coord4 = GenomicCoordinate("scaffold_123", 1, 100, '+')
        assert coord4.chromosome == "scaffold_123"


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_chromosome_start_bed(self):
        """Test coordinates at the very start of a chromosome (BED)."""
        # BED: [0, 100) is the first 100 bases
        coord = bed_to_internal("chr1", 0, 100, '+')
        assert coord.start == 1  # First base in 1-based system
        assert coord.stop == 100

    def test_chromosome_start_gff(self):
        """Test coordinates at the very start of a chromosome (GFF)."""
        # GFF: [1, 100] is the first 100 bases
        coord = gff_to_internal("chr1", 1, 100, '+')
        assert coord.start == 1
        assert coord.stop == 100

    def test_large_coordinates(self):
        """Test very large coordinate values (e.g., end of human chr1)."""
        # Human chr1 is ~249 million bases
        coord = GenomicCoordinate("chr1", 249_000_000, 249_250_621, '+')
        assert coord.start == 249_000_000
        assert coord.stop == 249_250_621

    def test_minimum_feature_size(self):
        """
        Test minimum feature size (1 base).

        Single-base features are now allowed to match original intronIC behavior.
        They will be marked as omitted during filtering (< 30bp minimum).
        """
        # Single base [100, 100] is now accepted (1bp feature)
        coord = GenomicCoordinate("chr1", 100, 100, '+')
        assert coord.length == 1

    def test_typical_intron_sizes(self):
        """Test coordinate handling for typical intron sizes."""
        # Small intron (30 bp minimum in intronIC)
        small = GenomicCoordinate("chr1", 1000, 1029, '+')
        assert small.stop - small.start + 1 == 30

        # Medium intron (~1kb)
        medium = GenomicCoordinate("chr1", 1000, 1999, '+')
        assert medium.stop - medium.start + 1 == 1000

        # Large intron (~100kb)
        large = GenomicCoordinate("chr1", 1000, 100999, '+')
        assert large.stop - large.start + 1 == 100000


class TestRealWorldExamples:
    """Test with real examples from intronIC gold standard."""

    def test_chr19_first_intron_bed_parsing(self):
        """
        Test parsing the first intron from chr19 gold standard.

        From homo_sapiens_gold.bed.iic:
        19  58345  58861  ENST00000587541:i1/1  0.00  +
        """
        coord = bed_to_internal("19", 58345, 58861, '+')

        assert coord.chromosome == "19"
        assert coord.start == 58346  # BED 58345 → internal 58346
        assert coord.stop == 58861
        assert coord.strand == '+'

        # Length should be 516 bases
        length = coord.stop - coord.start + 1
        assert length == 516

    def test_chr19_intron_bed_output(self):
        """Test generating BED output matching gold standard format."""
        # Internal representation of intron
        coord = GenomicCoordinate("19", 58346, 58861, '+', '1-based')

        # Convert to BED
        chrom, start, stop, strand = internal_to_bed(coord)

        # Should match original BED line coordinates
        assert chrom == "19"
        assert start == 58345
        assert stop == 58861
        assert strand == '+'

    def test_negative_strand_intron(self):
        """Test negative strand intron handling."""
        # Example negative strand intron
        coord = bed_to_internal("19", 100000, 101000, '-')

        assert coord.strand == '-'
        assert coord.start == 100001
        assert coord.stop == 101000

        # Length calculation should be strand-independent
        length = coord.stop - coord.start + 1
        assert length == 1000


class TestCoordinateSystem:
    """Test coordinate system tracking and validation."""

    def test_default_system_is_1based(self):
        """Test that default coordinate system is 1-based."""
        coord = GenomicCoordinate("chr1", 1, 100, '+')
        assert coord.system == '1-based'

    def test_explicit_0based_system(self):
        """Test creating coordinates with explicit 0-based system."""
        coord = GenomicCoordinate("chr1", 0, 100, '+', '0-based')
        assert coord.system == '0-based'
        assert coord.start == 0

    def test_explicit_1based_system(self):
        """Test creating coordinates with explicit 1-based system."""
        coord = GenomicCoordinate("chr1", 1, 100, '+', '1-based')
        assert coord.system == '1-based'
        assert coord.start == 1

    def test_invalid_system_rejected(self):
        """Test that invalid coordinate systems are rejected."""
        with pytest.raises(ValueError, match="system must be"):
            GenomicCoordinate("chr1", 1, 100, '+', 'invalid')

    def test_cannot_convert_0based_to_bed(self):
        """Test that we can only convert 1-based coords to BED."""
        coord = GenomicCoordinate("chr1", 0, 100, '+', '0-based')

        # This should fail because we expect internal coords to be 1-based
        with pytest.raises(ValueError, match="Can only convert 1-based"):
            internal_to_bed(coord)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
