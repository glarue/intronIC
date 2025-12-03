"""
Integration tests for parser → writer pipeline.

Tests the complete flow:
1. Parse annotation files (GFF3/GTF)
2. Parse genome files (FASTA)
3. Create intron objects
4. Write output in various formats

Uses real chr19 test data.

Author: intronIC refactoring project
Date: 2025-11-02
"""

from pathlib import Path

import pytest

from intronIC.core.intron import Intron, IntronFlags, IntronMetadata, IntronSequences
from intronIC.file_io.genome import GenomeReader
from intronIC.file_io.parsers import BEDParser, BioGLAnnotationParser, SequenceParser
from intronIC.file_io.writers import BEDWriter, MetaWriter, ScoreWriter, SequenceWriter
from intronIC.utils.coordinates import GenomicCoordinate

# Test data paths - use src layout path
TEST_DATA_DIR = (
    Path(__file__).parent.parent.parent / "src" / "intronIC" / "data" / "test_data"
)
CHR19_GENOME = TEST_DATA_DIR / "Homo_sapiens.Chr19.Ensembl_91.fa.gz"
CHR19_ANNOTATION = TEST_DATA_DIR / "Homo_sapiens.Chr19.Ensembl_91.gff3.gz"


# ============================================================================
# Annotation Parser Integration Tests
# ============================================================================


class TestAnnotationParserIntegration:
    """Test annotation parsing with real data."""

    @pytest.mark.skipif(
        not CHR19_ANNOTATION.exists(), reason="Chr19 annotation not available"
    )
    def test_parse_chr19_annotation(self):
        """Test parsing full chr19 annotation."""
        parser = BioGLAnnotationParser()
        features = list(parser.parse_file(CHR19_ANNOTATION))

        # Should have many features
        assert len(features) > 1000

        # Check feature type distribution
        feature_types = [f.feat_type for f in features]
        assert "gene" in feature_types
        assert "exon" in feature_types or "cds" in feature_types

    @pytest.mark.skipif(
        not CHR19_ANNOTATION.exists(), reason="Chr19 annotation not available"
    )
    def test_parse_chr19_structure(self):
        """Test that parsed chr19 data has proper structure."""
        parser = BioGLAnnotationParser()
        features = list(parser.parse_file(CHR19_ANNOTATION))

        # All features should have valid coordinates
        for feature in features[:100]:  # Check first 100
            assert feature.region is not None
            assert feature.start > 0
            assert feature.stop >= feature.start
            assert feature.strand in ["+", "-"]

    @pytest.mark.skipif(
        not CHR19_ANNOTATION.exists(), reason="Chr19 annotation not available"
    )
    def test_parse_chr19_parent_relationships(self):
        """Test that parent relationships are captured."""
        parser = BioGLAnnotationParser()
        features = list(parser.parse_file(CHR19_ANNOTATION))

        # Find exons/CDS features - they should have parents
        exons = [f for f in features if f.feat_type in ["exon", "cds"]]
        assert len(exons) > 0

        # Most exons should have parents
        with_parents = [e for e in exons if e.parent and e.parent != [None]]
        assert len(with_parents) > len(exons) * 0.5  # At least 50% should have parents


# ============================================================================
# Genome Reader Integration Tests
# ============================================================================


class TestGenomeReaderIntegration:
    """Test genome reading with real data."""

    @pytest.mark.skipif(not CHR19_GENOME.exists(), reason="Chr19 genome not available")
    def test_read_chr19_genome(self):
        """Test reading full chr19 genome."""
        reader = GenomeReader(CHR19_GENOME, cached=True)

        # Should have at least chr19
        chrom_names = reader.get_chromosome_names()
        assert len(chrom_names) > 0

        # Check one chromosome
        if "19" in chrom_names:
            chr19 = reader.get_sequence("19")
            assert len(chr19) > 1_000_000  # Chr19 should be > 1Mb
            assert all(c in "ACGTN" for c in chr19[:100])  # Should be uppercase DNA

    @pytest.mark.skipif(not CHR19_GENOME.exists(), reason="Chr19 genome not available")
    def test_extract_subsequence_from_chr19(self):
        """Test extracting subsequences from chr19."""
        reader = GenomeReader(CHR19_GENOME, cached=True)

        # Try to extract a small region
        coord = GenomicCoordinate("19", 1000, 1100, "+", "1-based")
        subseq = reader.extract_subsequence(coord)

        assert subseq is not None
        assert len(subseq) == 101  # 1000-1100 inclusive = 101 bases
        assert all(c in "ACGTN" for c in subseq)  # Should be uppercase DNA


# ============================================================================
# Parser → Intron → Writer Pipeline Tests
# ============================================================================


class TestParserWriterPipeline:
    """Test complete parser→intron→writer pipeline."""

    def test_bed_roundtrip(self, tmp_path):
        """Test BED file: parse → write → parse again."""
        # Create test BED file
        bed_file = tmp_path / "test.bed"
        bed_file.write_text(
            "chr1\t1000\t2000\tintron1\t95.5\t+\nchr1\t3000\t4000\tintron2\t12.3\t-\n"
        )

        # Parse BED file
        parser = BEDParser()
        bed_lines = list(parser.parse_file(bed_file))

        assert len(bed_lines) == 2
        assert bed_lines[0].chrom == "chr1"
        assert bed_lines[0].start == 1000  # 0-based
        assert bed_lines[0].stop == 2000

        # Create introns from BED data
        introns = []
        for bed in bed_lines:
            # Convert 0-based BED to 1-based coordinates
            coord = GenomicCoordinate(
                bed.chrom,
                bed.start + 1,  # Convert to 1-based
                bed.stop,
                bed.strand,
                "1-based",
            )
            intron = Intron(bed.name, coord)
            introns.append(intron)

        # Write back to BED format
        output_bed = tmp_path / "output.bed"
        with BEDWriter(output_bed) as writer:
            writer.write_introns(introns)

        # Parse output and verify
        output_lines = list(parser.parse_file(output_bed))
        assert len(output_lines) == 2
        assert output_lines[0].chrom == "chr1"
        assert output_lines[0].start == 1000  # Should be back to 0-based

    def test_sequence_roundtrip(self, tmp_path):
        """Test sequence file: parse → write → parse again."""
        # Create test sequence file
        # Format: name  score  upstream_flank  sequence  downstream_flank
        seq_file = tmp_path / "test.iic"
        seq_file.write_text(
            "intron1\t95.5\tAGGCT\tGTAAGTNNNNNTTTAG\tCATGG\n"
            "intron2\t12.3\tTTGCA\tGTAAGANNNNNTTTAG\tGGTAC\n"
        )

        # Parse sequence file
        parser = SequenceParser()
        seq_lines = list(parser.parse_file(seq_file))

        assert len(seq_lines) == 2
        assert seq_lines[0].name == "intron1"
        assert seq_lines[0].score == 95.5

        # Create introns from sequence data
        introns = []
        for seq_line in seq_lines:
            coord = GenomicCoordinate("chr1", 1000, 2000, "+", "1-based")
            sequences = IntronSequences(
                seq=seq_line.sequence,
                upstream_flank=seq_line.upstream_flank,
                downstream_flank=seq_line.downstream_flank,
            )
            intron = Intron(seq_line.name, coord, sequences=sequences)
            introns.append(intron)

        # Write back to sequence format
        output_seq = tmp_path / "output.iic"
        with SequenceWriter(output_seq) as writer:
            writer.write_introns(introns, include_score=False)

        # Parse output and verify
        output_lines = list(parser.parse_file(output_seq))
        assert len(output_lines) == 2
        assert output_lines[0].sequence == "GTAAGTNNNNNTTTAG"

    def test_multiple_format_output(self, tmp_path):
        """Test writing the same intron to multiple formats."""
        # Create a fully populated intron
        coord = GenomicCoordinate("chr1", 1001, 2000, "+", "1-based")
        sequences = IntronSequences(
            seq="GTAAGT" + "N" * 988 + "TTTAG",
            five_seq="AGGCTGTAAGT",
            three_seq="TTTAGCATGG",
            bp_seq="TACTAAC",
            upstream_flank="AGGCT",
            downstream_flank="CATGG",
        )
        metadata = IntronMetadata(
            parent="TRANS001",
            grandparent="GENE001",
            index=3,
            family_size=5,
            type_id="u2",
            flags=IntronFlags.LONGEST_ISOFORM,
        )
        intron = Intron("test_intron", coord, sequences=sequences, metadata=metadata)

        # Write to all formats
        bed_file = tmp_path / "test.bed.iic"
        meta_file = tmp_path / "test.meta.iic"
        seq_file = tmp_path / "test.introns.iic"
        score_file = tmp_path / "test.score_info.iic"

        with BEDWriter(bed_file) as writer:
            writer.write_intron(intron, species_name="test_species")

        with MetaWriter(meta_file) as writer:
            writer.write_header()
            writer.write_intron(intron, species_name="test_species")

        with SequenceWriter(seq_file) as writer:
            writer.write_intron(
                intron, species_name="test_species", include_score=False
            )

        with ScoreWriter(score_file) as writer:
            writer.write_header()
            writer.write_intron(intron, species_name="test_species")

        # Verify all files exist and have content
        assert bed_file.exists()
        assert meta_file.exists()
        assert seq_file.exists()
        assert score_file.exists()

        # Check that each has at least one line of data
        assert len(bed_file.read_text().strip().split("\n")) == 1
        assert len(meta_file.read_text().strip().split("\n")) == 2  # Header + data
        assert len(seq_file.read_text().strip().split("\n")) == 1
        assert len(score_file.read_text().strip().split("\n")) == 2  # Header + data

        # Verify names are consistent across files
        bed_name = bed_file.read_text().split("\t")[3].split(";")[0]
        meta_name = meta_file.read_text().split("\n")[1].split("\t")[0]
        seq_name = seq_file.read_text().split("\t")[0]
        score_name = score_file.read_text().split("\n")[1].split("\t")[0]

        # The expected name format is: {species_abbrev}-{grandparent}@{parent}_{index}({family_size})
        # For species "test_species" -> abbreviated to "TesSpe"
        expected_name = "TesSpe-GENE001@TRANS001_3(5)"
        assert bed_name == expected_name
        assert meta_name == expected_name
        assert seq_name == expected_name
        assert score_name == expected_name


# ============================================================================
# Real Data Integration Test
# ============================================================================


class TestRealDataPipeline:
    """Test pipeline with real chr19 data."""

    @pytest.mark.skipif(
        not (CHR19_ANNOTATION.exists() and CHR19_GENOME.exists()),
        reason="Chr19 test data not available",
    )
    def test_parse_chr19_and_write_outputs(self, tmp_path):
        """Test parsing chr19 annotation and writing sample outputs."""
        # Parse annotation
        parser = BioGLAnnotationParser()
        features = list(parser.parse_file(CHR19_ANNOTATION))

        # Find some exons to create test introns
        exons = [f for f in features if f.feat_type == "exon"][
            :10
        ]  # Use first 10 exons

        # Create sample introns (simplified - no actual intron extraction)
        introns = []
        for i, exon in enumerate(exons[:5]):  # Just 5 for testing
            coord = GenomicCoordinate(
                exon.region, exon.start, exon.stop, exon.strand, "1-based"
            )
            metadata = IntronMetadata(
                parent=exon.parent[0] if exon.parent else None,
                index=i + 1,
                family_size=5,
            )
            intron = Intron(f"test_intron_{i}", coord, metadata=metadata)
            introns.append(intron)

        # Write to various formats
        bed_file = tmp_path / "chr19_sample.bed.iic"
        meta_file = tmp_path / "chr19_sample.meta.iic"

        with BEDWriter(bed_file) as writer:
            count = writer.write_introns(introns)

        assert count == 5
        assert bed_file.exists()

        with MetaWriter(meta_file) as writer:
            writer.write_header()
            count = writer.write_introns(introns)

        assert count == 5
        assert meta_file.exists()

        # Verify output format
        bed_lines = bed_file.read_text().strip().split("\n")
        assert len(bed_lines) == 5

        meta_lines = meta_file.read_text().strip().split("\n")
        assert len(meta_lines) == 6  # Header + 5 introns

        # Each BED line should have 7 fields (including attributes)
        for line in bed_lines:
            fields = line.split("\t")
            assert len(fields) == 7

        # Each metadata line should have 15 fields (including attributes)
        for line in meta_lines[1:]:  # Skip header
            fields = line.split("\t")
            assert len(fields) == 15
            assert len(fields) == 15
