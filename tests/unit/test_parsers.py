"""
Unit tests for file format parsers.

Tests:
- BioGLAnnotationParser: GFF3/GTF annotation parsing
- BEDParser: BED format parsing
- SequenceParser: .iic sequence file parsing

Author: intronIC refactoring project
Date: 2025-11-02
"""

import gzip
from pathlib import Path

import pytest

from intronIC.file_io.parsers import (
    AnnotationLine,
    BEDLine,
    BEDParser,
    BioGLAnnotationParser,
    SequenceLine,
    SequenceParser,
)

# Test data paths - use fixtures in tests


# ============================================================================
# BioGLAnnotationParser Tests
# ============================================================================


class TestBioGLAnnotationParser:
    """Test annotation parser using biogl backend."""

    def test_parse_gff3_gene_line(self):
        """Test parsing a GFF3 gene line."""
        parser = BioGLAnnotationParser()
        line = "chr1\tENSEMBL\tgene\t1000\t2000\t.\t+\t.\tID=ENSG001;Name=GENE1"

        result = parser.parse_line(line, 1)

        assert result is not None
        assert result.name == "ENSG001"
        assert result.feat_type == "gene"
        assert result.region == "chr1"
        assert result.strand == "+"
        assert result.start == 1000
        assert result.stop == 2000
        assert result.line_number == 1
        assert (
            result.parent == []
        )  # biogl returns empty list for features without parents
        assert result.phase is None

    def test_parse_gff3_transcript_line(self):
        """Test parsing a GFF3 transcript/mRNA line."""
        parser = BioGLAnnotationParser()
        line = "chr1\tENSEMBL\tmRNA\t1000\t2000\t.\t+\t.\tID=ENST001;Parent=ENSG001"

        result = parser.parse_line(line, 2)

        assert result is not None
        assert result.name == "ENST001"
        assert result.feat_type == "transcript"  # biogl normalizes mRNA to transcript
        assert result.parent == ["ENSG001"]
        assert result.grandparent is None  # biogl may or may not populate this
        assert result.line_number == 2

    def test_parse_gff3_exon_line(self):
        """Test parsing a GFF3 exon line."""
        parser = BioGLAnnotationParser()
        line = "chr1\tENSEMBL\texon\t1000\t1200\t.\t+\t.\tID=ENSE001;Parent=ENST001"

        result = parser.parse_line(line, 3)

        assert result is not None
        assert result.name == "ENSE001"
        assert result.feat_type == "exon"
        assert result.parent == ["ENST001"]
        assert result.start == 1000
        assert result.stop == 1200

    def test_parse_gff3_cds_with_phase(self):
        """Test parsing a GFF3 CDS line with phase."""
        parser = BioGLAnnotationParser()
        line = "chr1\tENSEMBL\tCDS\t1000\t1200\t.\t+\t0\tID=CDS001;Parent=ENST001"

        result = parser.parse_line(line, 4)

        assert result is not None
        assert result.feat_type == "cds"  # normalized
        assert result.phase == 0

    def test_parse_negative_strand(self):
        """Test parsing negative strand feature."""
        parser = BioGLAnnotationParser()
        line = "chr1\tENSEMBL\texon\t1000\t1200\t.\t-\t.\tID=ENSE001;Parent=ENST001"

        result = parser.parse_line(line, 5)

        assert result is not None
        assert result.strand == "-"

    def test_parse_comment_line(self):
        """Test that comment lines return None."""
        parser = BioGLAnnotationParser()
        line = "#This is a comment"

        result = parser.parse_line(line, 1)

        assert result is None

    def test_parse_header_line(self):
        """Test that GFF3 header lines return None."""
        parser = BioGLAnnotationParser()
        line = "##gff-version 3"

        result = parser.parse_line(line, 1)

        assert result is None

    def test_parse_invalid_line(self):
        """Test that invalid lines return None."""
        parser = BioGLAnnotationParser()
        line = "not a valid gff3 line"

        result = parser.parse_line(line, 1)

        assert result is None

    def test_parse_phase_dot(self):
        """Test that phase '.' is converted to None."""
        parser = BioGLAnnotationParser()
        line = "chr1\tENSEMBL\texon\t1000\t1200\t.\t+\t.\tID=ENSE001;Parent=ENST001"

        result = parser.parse_line(line, 1)

        assert result is not None
        assert result.phase is None

    def test_parse_phase_values(self):
        """Test parsing all valid phase values (0, 1, 2)."""
        parser = BioGLAnnotationParser()

        for phase_val in [0, 1, 2]:
            line = f"chr1\tENSEMBL\tCDS\t1000\t1200\t.\t+\t{phase_val}\tID=CDS001;Parent=ENST001"
            result = parser.parse_line(line, 1)
            assert result is not None
            assert result.phase == phase_val

    def test_parse_file_basic(self, tmp_path):
        """Test parsing a complete GFF3 file."""
        gff3_file = tmp_path / "test.gff3"
        gff3_file.write_text(
            "##gff-version 3\n"
            "chr1\tENSEMBL\tgene\t1000\t2000\t.\t+\t.\tID=ENSG001\n"
            "chr1\tENSEMBL\tmRNA\t1000\t2000\t.\t+\t.\tID=ENST001;Parent=ENSG001\n"
            "chr1\tENSEMBL\texon\t1000\t1200\t.\t+\t.\tID=ENSE001;Parent=ENST001\n"
            "chr1\tENSEMBL\texon\t1500\t2000\t.\t+\t.\tID=ENSE002;Parent=ENST001\n"
        )

        parser = BioGLAnnotationParser()
        results = list(parser.parse_file(gff3_file))

        # Should have 4 features (header is skipped)
        assert len(results) == 4
        assert results[0].feat_type == "gene"
        assert results[1].feat_type == "transcript"  # biogl normalizes mRNA
        assert results[2].feat_type == "exon"
        assert results[3].feat_type == "exon"

    def test_parse_gzipped_file(self, tmp_path):
        """Test parsing a gzipped GFF3 file."""
        gff3_file = tmp_path / "test.gff3.gz"
        content = ("chr1\tENSEMBL\tgene\t1000\t2000\t.\t+\t.\tID=ENSG001\n").encode(
            "utf-8"
        )

        with gzip.open(gff3_file, "wb") as f:
            f.write(content)

        parser = BioGLAnnotationParser()
        results = list(parser.parse_file(gff3_file))

        assert len(results) == 1
        assert results[0].feat_type == "gene"

    def test_parse_real_chr19_annotation(self, test_data_dir):
        """Test parsing real chr19 GFF3 data."""
        chr19_annotation = test_data_dir / "Homo_sapiens.Chr19.Ensembl_91.gff3.gz"
        if not chr19_annotation.exists():
            pytest.skip("Chr19 test data not available")

        parser = BioGLAnnotationParser()
        results = list(parser.parse_file(chr19_annotation))

        # Should have many features
        assert len(results) > 1000

        # Check we have various feature types
        feature_types = {r.feat_type for r in results}
        assert "gene" in feature_types
        assert "exon" in feature_types or "cds" in feature_types

    def test_parse_chr19_gene_structure(self, test_data_dir):
        """Test that chr19 parsing captures gene structure correctly."""
        chr19_annotation = test_data_dir / "Homo_sapiens.Chr19.Ensembl_91.gff3.gz"
        if not chr19_annotation.exists():
            pytest.skip("Chr19 test data not available")

        parser = BioGLAnnotationParser()
        results = list(parser.parse_file(chr19_annotation))

        # Find a gene and its children
        genes = [r for r in results if r.feat_type == "gene"]
        assert len(genes) > 0

        # Check that genes have proper coordinates
        for gene in genes[:10]:  # Check first 10
            assert gene.start > 0
            assert gene.stop > gene.start
            assert gene.strand in ["+", "-"]


# ============================================================================
# BEDParser Tests
# ============================================================================


class TestBEDParser:
    """Test BED format parser."""

    def test_parse_bed3(self):
        """Test parsing minimal BED3 format (chrom, start, stop)."""
        parser = BEDParser()
        line = "chr1\t1000\t2000"

        result = parser.parse_line(line)

        assert result is not None
        assert result.chrom == "chr1"
        assert result.start == 1000
        assert result.stop == 2000
        assert result.name == "."
        assert result.score == "."
        assert result.strand == "."

    def test_parse_bed6(self):
        """Test parsing standard BED6 format."""
        parser = BEDParser()
        line = "chr1\t1000\t2000\tfeature1\t100\t+"

        result = parser.parse_line(line)

        assert result is not None
        assert result.chrom == "chr1"
        assert result.start == 1000
        assert result.stop == 2000
        assert result.name == "feature1"
        assert result.score == "100"
        assert result.strand == "+"

    def test_parse_bed_negative_strand(self):
        """Test parsing BED with negative strand."""
        parser = BEDParser()
        line = "chr1\t1000\t2000\tfeature1\t100\t-"

        result = parser.parse_line(line)

        assert result is not None
        assert result.strand == "-"

    def test_parse_bed_with_extra_columns(self):
        """Test parsing BED with extra columns beyond standard 6."""
        parser = BEDParser()
        line = "chr1\t1000\t2000\tfeature1\t100\t+\t1000\t2000\t255,0,0"

        result = parser.parse_line(line)

        assert result is not None
        assert len(result.extra_fields) == 3
        assert result.extra_fields[0] == "1000"
        assert result.extra_fields[1] == "2000"
        assert result.extra_fields[2] == "255,0,0"

    def test_parse_bed_0based_coordinates(self):
        """Test that BED coordinates are correctly parsed as 0-based."""
        parser = BEDParser()
        line = "chr1\t0\t100\tfeature1\t100\t+"

        result = parser.parse_line(line)

        # Start should be 0 (BED is 0-based)
        assert result.start == 0
        assert result.stop == 100

    def test_parse_comment_hash(self):
        """Test that comment lines starting with # are skipped."""
        parser = BEDParser()
        line = "#This is a comment"

        result = parser.parse_line(line)

        assert result is None

    def test_parse_track_line(self):
        """Test that track lines are skipped."""
        parser = BEDParser()
        line = 'track name=myTrack description="My Track"'

        result = parser.parse_line(line)

        assert result is None

    def test_parse_browser_line(self):
        """Test that browser lines are skipped."""
        parser = BEDParser()
        line = "browser position chr1:1000-2000"

        result = parser.parse_line(line)

        assert result is None

    def test_parse_empty_line(self):
        """Test that empty lines are skipped."""
        parser = BEDParser()
        line = ""

        result = parser.parse_line(line)

        assert result is None

    def test_parse_whitespace_only_line(self):
        """Test that whitespace-only lines are skipped."""
        parser = BEDParser()
        line = "   \t  \n"

        result = parser.parse_line(line)

        assert result is None

    def test_parse_invalid_coordinates(self):
        """Test that lines with non-numeric coordinates return None."""
        parser = BEDParser()
        line = "chr1\tabc\t2000\tfeature1\t100\t+"

        result = parser.parse_line(line)

        assert result is None

    def test_parse_too_few_columns(self):
        """Test that lines with fewer than 3 columns return None."""
        parser = BEDParser()
        line = "chr1\t1000"

        result = parser.parse_line(line)

        assert result is None

    def test_parse_file_basic(self, tmp_path):
        """Test parsing a complete BED file."""
        bed_file = tmp_path / "test.bed"
        bed_file.write_text(
            "#Comment line\n"
            "chr1\t1000\t2000\tfeature1\t100\t+\n"
            "chr1\t3000\t4000\tfeature2\t200\t-\n"
            "chr2\t5000\t6000\tfeature3\t150\t+\n"
        )

        parser = BEDParser()
        results = list(parser.parse_file(bed_file))

        assert len(results) == 3
        assert results[0].name == "feature1"
        assert results[1].name == "feature2"
        assert results[2].name == "feature3"

    def test_parse_gzipped_bed_file(self, tmp_path):
        """Test parsing a gzipped BED file."""
        bed_file = tmp_path / "test.bed.gz"
        content = "chr1\t1000\t2000\tfeature1\t100\t+\n".encode("utf-8")

        with gzip.open(bed_file, "wb") as f:
            f.write(content)

        parser = BEDParser()
        results = list(parser.parse_file(bed_file))

        assert len(results) == 1
        assert results[0].name == "feature1"


# ============================================================================
# SequenceParser Tests
# ============================================================================


class TestSequenceParser:
    """Test .iic sequence file parser."""

    def test_parse_basic_sequence_line(self):
        """Test parsing a basic 4-column sequence line."""
        parser = SequenceParser()
        line = "intron1\tAGGCT\tGTAAGTTTTTTTTTTTTTTTTTTTTTTTCAG\tCATGG"

        result = parser.parse_line(line)

        assert result is not None
        assert result.name == "intron1"
        assert result.upstream_flank == "AGGCT"
        assert (
            result.sequence == "GTAAGTTTTTTTTTTTTTTTTTTTTTTTCAG"
        )  # Match the input exactly
        assert result.downstream_flank == "CATGG"
        assert result.score is None

    def test_parse_sequence_with_score(self):
        """Test parsing a 5-column sequence line with score."""
        parser = SequenceParser()
        # Format: name, score, upstream, sequence, downstream
        line = "intron1\t95.5\tAGGCT\tGTAAGT\tCATGG"

        result = parser.parse_line(line)

        assert result is not None
        assert result.name == "intron1"
        assert result.score == 95.5

    def test_parse_sequence_with_zero_score(self):
        """Test parsing sequence with score of 0."""
        parser = SequenceParser()
        # Format: name, score, upstream, sequence, downstream
        line = "intron1\t0.0\tAGGCT\tGTAAGT\tCATGG"

        result = parser.parse_line(line)

        assert result is not None
        assert result.score == 0.0

    def test_parse_sequence_with_100_score(self):
        """Test parsing sequence with maximum score."""
        parser = SequenceParser()
        # Format: name, score, upstream, sequence, downstream
        line = "intron1\t100.0\tAGGCT\tGTAAGT\tCATGG"

        result = parser.parse_line(line)

        assert result is not None
        assert result.score == 100.0

    def test_parse_long_sequence(self):
        """Test parsing a long intron sequence."""
        parser = SequenceParser()
        long_seq = "A" * 10000
        line = f"intron1\tAGGCT\t{long_seq}\tCATGG"

        result = parser.parse_line(line)

        assert result is not None
        assert len(result.sequence) == 10000

    def test_parse_empty_flanks(self):
        """Test parsing sequence with empty flanking sequences."""
        parser = SequenceParser()
        # With 5 fields and empty flanks on the edges, trailing tabs get stripped
        # So we test with empty upstream in the middle which is preserved
        # Format: name, score, upstream (empty), sequence, downstream (x to preserve)
        line = "intron1\t95.0\t\tGTAAGT\tx"  # 5 fields, empty upstream

        result = parser.parse_line(line)

        assert result is not None
        assert result.upstream_flank == ""
        assert result.sequence == "GTAAGT"
        assert result.downstream_flank == "x"
        assert result.score == 95.0

    def test_parse_empty_line(self):
        """Test that empty lines return None."""
        parser = SequenceParser()
        line = ""

        result = parser.parse_line(line)

        assert result is None

    def test_parse_whitespace_only_line(self):
        """Test that whitespace-only lines return None."""
        parser = SequenceParser()
        line = "   \t  \n"

        result = parser.parse_line(line)

        assert result is None

    def test_parse_too_few_columns(self):
        """Test that lines with fewer than 4 columns return None."""
        parser = SequenceParser()
        line = "intron1\tAGGCT\tGTAAGT"

        result = parser.parse_line(line)

        assert result is None

    def test_parse_invalid_score(self):
        """Test that lines with non-numeric score return None."""
        parser = SequenceParser()
        line = "intron1\tAGGCT\tGTAAGT\tCATGG\tabc"

        result = parser.parse_line(line)

        assert result is None

    def test_parse_file_basic(self, tmp_path):
        """Test parsing a complete .iic sequence file."""
        iic_file = tmp_path / "test.iic"
        # Format: name, score, upstream, sequence, downstream
        iic_file.write_text(
            "intron1\t95.5\tAGGCT\tGTAAGT\tCATGG\n"
            "intron2\t87.3\tTTGCA\tGTAAGA\tGGTAC\n"
            "intron3\t12.1\tCCAGT\tGTAAGC\tATTGG\n"
        )

        parser = SequenceParser()
        results = list(parser.parse_file(iic_file))

        assert len(results) == 3
        assert results[0].name == "intron1"
        assert results[0].score == 95.5
        assert results[1].name == "intron2"
        assert results[1].score == 87.3
        assert results[2].name == "intron3"
        assert results[2].score == 12.1

    def test_parse_gzipped_iic_file(self, tmp_path):
        """Test parsing a gzipped .iic file."""
        iic_file = tmp_path / "test.iic.gz"
        # Format: name, score, upstream, sequence, downstream
        content = "intron1\t95.5\tAGGCT\tGTAAGT\tCATGG\n".encode("utf-8")

        with gzip.open(iic_file, "wb") as f:
            f.write(content)

        parser = SequenceParser()
        results = list(parser.parse_file(iic_file))

        assert len(results) == 1
        assert results[0].name == "intron1"
        assert results[0].score == 95.5

    def test_parse_file_without_scores(self, tmp_path):
        """Test parsing .iic file without score column."""
        iic_file = tmp_path / "test.iic"
        iic_file.write_text(
            "intron1\tAGGCT\tGTAAGT\tCATGG\nintron2\tTTGCA\tGTAAGA\tGGTAC\n"
        )

        parser = SequenceParser()
        results = list(parser.parse_file(iic_file))

        assert len(results) == 2
        assert results[0].score is None
        assert results[1].score is None


# ============================================================================
# Data Structure Tests
# ============================================================================


class TestDataStructures:
    """Test the dataclass structures themselves."""

    def test_annotation_line_creation(self):
        """Test creating an AnnotationLine object."""
        line = AnnotationLine(
            name="ENSG001",
            feat_type="gene",
            parent=[],
            grandparent=None,
            region="chr1",
            strand="+",
            start=1000,
            stop=2000,
            line_number=1,
            phase=None,
        )

        assert line.name == "ENSG001"
        assert line.feat_type == "gene"
        assert line.start == 1000

    def test_bed_line_creation(self):
        """Test creating a BEDLine object."""
        bed = BEDLine(
            chrom="chr1",
            start=1000,
            stop=2000,
            name="feature1",
            score="100",
            strand="+",
        )

        assert bed.chrom == "chr1"
        assert bed.start == 1000
        assert bed.strand == "+"

    def test_sequence_line_creation(self):
        """Test creating a SequenceLine object."""
        seq = SequenceLine(
            name="intron1",
            upstream_flank="AGGCT",
            sequence="GTAAGT",
            downstream_flank="CATGG",
            score=95.5,
        )

        assert seq.name == "intron1"
        assert seq.score == 95.5


# ============================================================================
# Graph Cycle Prevention Tests
# ============================================================================


class TestCyclePrevention:
    """Test that ID cleaning prevents graph cycles."""

    def test_refseq_style_no_cycles(self):
        """
        Test that RefSeq-style prefixes (gene-, rna-) are NOT cleaned,
        preventing ID collisions and graph cycles.

        This simulates the tRNA annotation pattern that caused 671 self-loops
        in Basidiobolus when gene- and rna- prefixes were incorrectly removed.
        """
        from networkx import is_directed_acyclic_graph

        from intronIC.extraction.annotator import AnnotationHierarchyBuilder

        # RefSeq-style annotation: gene and tRNA share same base ID
        # gene-K493DRAFT_t33 and rna-K493DRAFT_t33
        # If we incorrectly clean to just "K493DRAFT_t33", we get a self-loop!
        gff_lines = [
            "chr1\tGenbank\tgene\t1000\t2000\t.\t+\t.\tID=gene-TEST_t1;Name=TEST_t1",
            "chr1\tGenbank\ttRNA\t1000\t2000\t.\t+\t.\tID=rna-TEST_t1;Parent=gene-TEST_t1",
            "chr1\tGenbank\texon\t1000\t2000\t.\t+\t.\tID=exon-TEST_t1-1;Parent=rna-TEST_t1",
        ]

        builder = AnnotationHierarchyBuilder(
            child_features=["exon"],
            clean_names=True,  # Use default cleaning
        )

        # Parse lines
        from intronIC.file_io.parsers import BioGLAnnotationParser

        parser = BioGLAnnotationParser(clean_names=True)
        annotations = list(parser.parse_lines(gff_lines))

        # Build graph
        from networkx import DiGraph

        feat_graph = DiGraph()

        for ann in annotations:
            if ann.start > ann.stop:
                continue

            features = builder._create_features_from_annotation(ann)

            for feat in features:
                feat_type = feat.attributes.get("_orig_feat_type", "")
                if not feat_type and hasattr(feat, "feature_type"):
                    feat_type = feat.feature_type
                feat_type = feat_type.lower() if feat_type else "unknown"

                parent_name = getattr(feat, "parent_id", None) or feat.attributes.get(
                    "_parent_name"
                )
                grandparent_name = feat.attributes.get("_grandparent_name")

                if feat_type in builder.child_features:
                    name = f"{feat_type}_Parent={parent_name}:{feat.start}_{feat.stop}"
                else:
                    name = feat.feature_id

                # Build graph edges
                if parent_name is not None:
                    feat_graph.add_edge(parent_name, name)

                if grandparent_name is not None:
                    feat_graph.add_edge(grandparent_name, parent_name)

        # CRITICAL: Graph must be acyclic
        # If gene- and rna- prefixes were removed, we'd have:
        # - gene-TEST_t1 -> TEST_t1 (after cleaning)
        # - rna-TEST_t1 -> TEST_t1 (after cleaning)
        # - Edge: TEST_t1 -> TEST_t1 (SELF-LOOP!)
        assert is_directed_acyclic_graph(feat_graph), (
            "Graph contains cycles! RefSeq-style prefixes (gene-, rna-) must NOT be cleaned"
        )

    def test_ensembl_style_cleaning(self):
        """
        Test that Ensembl-style prefixes (gene:, transcript:) ARE cleaned
        without causing cycles.

        Ensembl IDs are globally unique even without prefixes (ENSG*, ENST*),
        so cleaning is safe and produces cleaner output.
        """
        from networkx import is_directed_acyclic_graph

        from intronIC.extraction.annotator import AnnotationHierarchyBuilder

        # Ensembl-style annotation with redundant prefixes
        gff_lines = [
            "chr1\tENSEMBL\tgene\t1000\t2000\t.\t+\t.\tID=gene:ENSG00000001;Name=GENE1",
            "chr1\tENSEMBL\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:ENST00000001;Parent=gene:ENSG00000001",
            "chr1\tENSEMBL\texon\t1000\t1500\t.\t+\t.\tID=exon:ENSE00000001;Parent=transcript:ENST00000001",
        ]

        builder = AnnotationHierarchyBuilder(child_features=["exon"], clean_names=True)

        from intronIC.file_io.parsers import BioGLAnnotationParser

        parser = BioGLAnnotationParser(clean_names=True)
        annotations = list(parser.parse_lines(gff_lines))

        # Check that prefixes were cleaned
        gene_ann = [a for a in annotations if a.feat_type == "gene"][0]
        # biogl converts mRNA feature type to 'transcript'
        transcript_ann = [a for a in annotations if a.feat_type == "transcript"][0]

        assert gene_ann.name == "ENSG00000001", (
            f"Expected ENSG00000001, got {gene_ann.name}"
        )
        assert transcript_ann.name == "ENST00000001", (
            f"Expected ENST00000001, got {transcript_ann.name}"
        )

        # Build graph
        from networkx import DiGraph

        feat_graph = DiGraph()

        for ann in annotations:
            if ann.start > ann.stop:
                continue

            features = builder._create_features_from_annotation(ann)

            for feat in features:
                feat_type = feat.attributes.get("_orig_feat_type", "")
                if not feat_type and hasattr(feat, "feature_type"):
                    feat_type = feat.feature_type
                feat_type = feat_type.lower() if feat_type else "unknown"

                parent_name = getattr(feat, "parent_id", None) or feat.attributes.get(
                    "_parent_name"
                )
                grandparent_name = feat.attributes.get("_grandparent_name")

                if feat_type in builder.child_features:
                    name = f"{feat_type}_Parent={parent_name}:{feat.start}_{feat.stop}"
                else:
                    name = feat.feature_id

                if parent_name is not None:
                    feat_graph.add_edge(parent_name, name)

                if grandparent_name is not None:
                    feat_graph.add_edge(grandparent_name, parent_name)

        # Graph must still be acyclic even with cleaning
        assert is_directed_acyclic_graph(feat_graph), (
            "Graph contains cycles! Ensembl-style cleaning should be safe"
        )

    def test_mirna_bidirectional_cycle_prevention(self):
        """
        Test that miRNA annotations don't create bidirectional cycles.

        This simulates the human genome pattern where:
        - gene-MIR6859-1 (gene)
        - rna-NR_106918.1 (transcript, parent=gene-MIR6859-1)
        - rna-MIR6859-1 (miRNA, parent=rna-NR_106918.1)

        If gene- and rna- are cleaned, we get:
        - MIR6859-1 -> NR_106918.1
        - NR_106918.1 -> MIR6859-1
        Creating a bidirectional cycle!
        """
        from networkx import is_directed_acyclic_graph

        from intronIC.extraction.annotator import AnnotationHierarchyBuilder

        gff_lines = [
            "chr1\tRefSeq\tgene\t1000\t2000\t.\t+\t.\tID=gene-MIR6859-1;Name=MIR6859-1",
            "chr1\tRefSeq\ttranscript\t1000\t2000\t.\t+\t.\tID=rna-NR_106918.1;Parent=gene-MIR6859-1",
            "chr1\tRefSeq\texon\t1000\t2000\t.\t+\t.\tID=exon-1;Parent=rna-NR_106918.1",
            "chr1\tRefSeq\tmiRNA\t1000\t1500\t.\t+\t.\tID=rna-MIR6859-1;Parent=rna-NR_106918.1",
        ]

        builder = AnnotationHierarchyBuilder(child_features=["exon"], clean_names=True)

        from intronIC.file_io.parsers import BioGLAnnotationParser

        parser = BioGLAnnotationParser(clean_names=True)
        annotations = list(parser.parse_lines(gff_lines))

        # Build graph
        from networkx import DiGraph

        feat_graph = DiGraph()

        for ann in annotations:
            if ann.start > ann.stop:
                continue

            features = builder._create_features_from_annotation(ann)

            for feat in features:
                feat_type = feat.attributes.get("_orig_feat_type", "")
                if not feat_type and hasattr(feat, "feature_type"):
                    feat_type = feat.feature_type
                feat_type = feat_type.lower() if feat_type else "unknown"

                parent_name = getattr(feat, "parent_id", None) or feat.attributes.get(
                    "_parent_name"
                )
                grandparent_name = feat.attributes.get("_grandparent_name")

                if feat_type in builder.child_features:
                    name = f"{feat_type}_Parent={parent_name}:{feat.start}_{feat.stop}"
                else:
                    name = feat.feature_id

                if parent_name is not None:
                    feat_graph.add_edge(parent_name, name)

                if grandparent_name is not None:
                    feat_graph.add_edge(grandparent_name, parent_name)

        # CRITICAL: Must not have cycles
        # With improper cleaning:
        # - gene-MIR6859-1 cleans to MIR6859-1
        # - rna-MIR6859-1 cleans to MIR6859-1 (collision!)
        # - Edge: MIR6859-1 -> NR_106918.1 -> MIR6859-1 (cycle!)
        assert is_directed_acyclic_graph(feat_graph), (
            "Graph contains cycles! miRNA gene/rna- ID collision detected"
        )
        # - rna-MIR6859-1 cleans to MIR6859-1 (collision!)
        # - Edge: MIR6859-1 -> NR_106918.1 -> MIR6859-1 (cycle!)
        assert is_directed_acyclic_graph(feat_graph), (
            "Graph contains cycles! miRNA gene/rna- ID collision detected"
        )
