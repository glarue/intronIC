"""
Unit tests for core data models (GenomeFeature hierarchy).

Tests:
- GenomeFeature: Basic genomic feature representation
- Gene: Gene with multiple transcripts
- Transcript: Transcript with multiple exons
- Exon: Exon representation with coding information

Author: intronIC refactoring project
Date: 2025-11-02
"""

import pytest
from intronIC.utils.coordinates import GenomicCoordinate
from intronIC.core.models import GenomeFeature, Gene, Transcript, Exon


class TestGenomeFeature:
    """Test GenomeFeature base class."""

    def test_create_feature(self):
        """Test creating a basic genome feature."""
        coord = GenomicCoordinate("chr1", 1000, 2000, '+', '1-based')
        feature = GenomeFeature("feat1", coord, "exon")

        assert feature.feature_id == "feat1"
        assert feature.feature_type == "exon"
        assert feature.chromosome == "chr1"
        assert feature.start == 1000
        assert feature.stop == 2000
        assert feature.strand == '+'
        assert feature.length == 1001

    def test_feature_with_attributes(self):
        """Test feature with additional attributes."""
        coord = GenomicCoordinate("chr1", 1000, 2000, '+', '1-based')
        attrs = {"Name": "MyExon", "score": 100}
        feature = GenomeFeature("feat1", coord, "exon", attrs)

        assert feature.attributes["Name"] == "MyExon"
        assert feature.attributes["score"] == 100

    def test_empty_feature_id_rejected(self):
        """Test that empty feature_id is rejected."""
        coord = GenomicCoordinate("chr1", 1000, 2000, '+', '1-based')
        with pytest.raises(ValueError, match="feature_id cannot be empty"):
            GenomeFeature("", coord, "exon")

    def test_empty_feature_type_rejected(self):
        """Test that empty feature_type is rejected."""
        coord = GenomicCoordinate("chr1", 1000, 2000, '+', '1-based')
        with pytest.raises(ValueError, match="feature_type cannot be empty"):
            GenomeFeature("feat1", coord, "")

    def test_requires_1based_coordinates(self):
        """Test that only 1-based coordinates are accepted."""
        coord = GenomicCoordinate("chr1", 1000, 2000, '+', '0-based')
        with pytest.raises(ValueError, match="requires 1-based coordinates"):
            GenomeFeature("feat1", coord, "exon")

    def test_str_representation(self):
        """Test string representation."""
        coord = GenomicCoordinate("chr1", 1000, 2000, '+', '1-based')
        feature = GenomeFeature("feat1", coord, "exon")
        str_repr = str(feature)

        assert "exon" in str_repr
        assert "feat1" in str_repr
        assert "chr1" in str_repr


class TestGene:
    """Test Gene class."""

    def test_create_gene(self):
        """Test creating a gene."""
        coord = GenomicCoordinate("chr1", 1000, 5000, '+', '1-based')
        gene = Gene("GENE001", coord)

        assert gene.feature_id == "GENE001"
        assert gene.feature_type == "gene"
        assert gene.num_children == 0
        assert len(gene.children) == 0

    def test_gene_with_name(self):
        """Test gene with gene_name attribute."""
        coord = GenomicCoordinate("chr1", 1000, 5000, '+', '1-based')
        attrs = {"gene_name": "BRCA1"}
        gene = Gene("ENSG001", coord, attrs)

        assert gene.gene_name == "BRCA1"

    def test_add_child_transcript(self):
        """Test adding transcripts to gene."""
        coord = GenomicCoordinate("chr1", 1000, 5000, '+', '1-based')
        gene = Gene("GENE001", coord)

        gene.add_child("TRANS001")
        gene.add_child("TRANS002")

        assert gene.num_children == 2
        assert "TRANS001" in gene.children
        assert "TRANS002" in gene.children
        assert "TRANS001" in gene.transcript_ids

    def test_remove_child_transcript(self):
        """Test removing transcripts from gene."""
        coord = GenomicCoordinate("chr1", 1000, 5000, '+', '1-based')
        gene = Gene("GENE001", coord)

        gene.add_child("TRANS001")
        gene.add_child("TRANS002")
        gene.remove_child("TRANS001")

        assert gene.num_children == 1
        assert "TRANS001" not in gene.children
        assert "TRANS002" in gene.children

    def test_gene_coordinates(self):
        """Test gene coordinate properties."""
        coord = GenomicCoordinate("chr2", 10000, 20000, '-', '1-based')
        gene = Gene("GENE002", coord)

        assert gene.chromosome == "chr2"
        assert gene.start == 10000
        assert gene.stop == 20000
        assert gene.strand == '-'
        assert gene.length == 10001

    def test_gene_str_representation(self):
        """Test gene string representation."""
        coord = GenomicCoordinate("chr1", 1000, 5000, '+', '1-based')
        gene = Gene("ENSG001", coord, {"gene_name": "TP53"})
        gene.add_child("TRANS001")

        str_repr = str(gene)
        assert "TP53" in str_repr
        assert "1 transcript" in str_repr


class TestTranscript:
    """Test Transcript class."""

    def test_create_transcript(self):
        """Test creating a transcript."""
        coord = GenomicCoordinate("chr1", 1000, 5000, '+', '1-based')
        trans = Transcript("TRANS001", coord)

        assert trans.feature_id == "TRANS001"
        assert trans.feature_type == "transcript"
        assert trans.num_children == 0
        assert trans.parent_id is None

    def test_transcript_with_parent(self):
        """Test transcript with parent gene."""
        coord = GenomicCoordinate("chr1", 1000, 5000, '+', '1-based')
        trans = Transcript("TRANS001", coord, parent_id="GENE001")

        assert trans.parent_id == "GENE001"
        assert trans.gene_id == "GENE001"

    def test_add_child_exon(self):
        """Test adding exons to transcript."""
        coord = GenomicCoordinate("chr1", 1000, 5000, '+', '1-based')
        trans = Transcript("TRANS001", coord)

        trans.add_child("EXON001")
        trans.add_child("EXON002")
        trans.add_child("EXON003")

        assert trans.num_children == 3
        assert "EXON001" in trans.exon_ids
        assert "EXON002" in trans.exon_ids

    def test_transcript_with_name(self):
        """Test transcript with transcript_name attribute."""
        coord = GenomicCoordinate("chr1", 1000, 5000, '+', '1-based')
        attrs = {"transcript_name": "BRCA1-001"}
        trans = Transcript("ENST001", coord, attrs)

        assert trans.transcript_name == "BRCA1-001"

    def test_transcript_str_representation(self):
        """Test transcript string representation."""
        coord = GenomicCoordinate("chr1", 1000, 5000, '+', '1-based')
        trans = Transcript("TRANS001", coord, parent_id="GENE001")
        trans.add_child("EXON001")
        trans.add_child("EXON002")

        str_repr = str(trans)
        assert "TRANS001" in str_repr
        assert "2 exons" in str_repr
        assert "GENE001" in str_repr


class TestExon:
    """Test Exon class."""

    def test_create_exon(self):
        """Test creating an exon."""
        coord = GenomicCoordinate("chr1", 1000, 1200, '+', '1-based')
        exon = Exon("EXON001", coord)

        assert exon.feature_id == "EXON001"
        assert exon.chromosome == "chr1"
        assert exon.start == 1000
        assert exon.stop == 1200
        assert exon.length == 201
        assert exon.is_coding == False
        assert exon.phase is None

    def test_coding_exon_with_phase(self):
        """Test creating a coding exon with phase."""
        coord = GenomicCoordinate("chr1", 1000, 1200, '+', '1-based')
        exon = Exon("EXON001", coord, phase=0, is_coding=True)

        assert exon.is_coding == True
        assert exon.phase == 0

    def test_exon_with_parent(self):
        """Test exon with parent transcript."""
        coord = GenomicCoordinate("chr1", 1000, 1200, '+', '1-based')
        exon = Exon("EXON001", coord, parent_id="TRANS001")

        assert exon.parent_id == "TRANS001"
        assert exon.transcript_id == "TRANS001"

    def test_exon_with_number(self):
        """Test exon with exon_number attribute."""
        coord = GenomicCoordinate("chr1", 1000, 1200, '+', '1-based')
        attrs = {"exon_number": "1"}
        exon = Exon("EXON001", coord, attributes=attrs)

        assert exon.exon_number == 1

    def test_invalid_phase_rejected(self):
        """Test that invalid phase values are rejected."""
        coord = GenomicCoordinate("chr1", 1000, 1200, '+', '1-based')
        with pytest.raises(ValueError, match="phase must be"):
            Exon("EXON001", coord, phase=3)

        with pytest.raises(ValueError, match="phase must be"):
            Exon("EXON001", coord, phase=-1)

    def test_exon_phases_valid(self):
        """Test that all valid phases (0, 1, 2) are accepted."""
        coord = GenomicCoordinate("chr1", 1000, 1200, '+', '1-based')

        exon0 = Exon("EXON0", coord, phase=0)
        exon1 = Exon("EXON1", coord, phase=1)
        exon2 = Exon("EXON2", coord, phase=2)

        assert exon0.phase == 0
        assert exon1.phase == 1
        assert exon2.phase == 2

    def test_exon_immutable(self):
        """Test that exon is immutable (frozen)."""
        coord = GenomicCoordinate("chr1", 1000, 1200, '+', '1-based')
        exon = Exon("EXON001", coord)

        with pytest.raises(Exception):  # FrozenInstanceError
            exon.phase = 1

    def test_exon_str_representation(self):
        """Test exon string representation."""
        coord = GenomicCoordinate("chr1", 1000, 1200, '+', '1-based')
        exon = Exon("EXON001", coord, parent_id="TRANS001", phase=0, is_coding=True)

        str_repr = str(exon)
        assert "EXON001" in str_repr
        assert "CDS" in str_repr
        assert "TRANS001" in str_repr


class TestHierarchyIntegration:
    """Test integration of Gene -> Transcript -> Exon hierarchy."""

    def test_gene_transcript_exon_hierarchy(self):
        """Test creating a complete gene hierarchy."""
        # Create gene
        gene_coord = GenomicCoordinate("chr1", 1000, 5000, '+', '1-based')
        gene = Gene("GENE001", gene_coord, {"gene_name": "MyGene"})

        # Create transcript
        trans_coord = GenomicCoordinate("chr1", 1000, 5000, '+', '1-based')
        trans = Transcript("TRANS001", trans_coord, parent_id="GENE001")

        # Create exons
        exon1_coord = GenomicCoordinate("chr1", 1000, 1200, '+', '1-based')
        exon2_coord = GenomicCoordinate("chr1", 2000, 2300, '+', '1-based')
        exon3_coord = GenomicCoordinate("chr1", 4500, 5000, '+', '1-based')

        exon1 = Exon("EXON001", exon1_coord, parent_id="TRANS001", phase=0, is_coding=True)
        exon2 = Exon("EXON002", exon2_coord, parent_id="TRANS001", phase=0, is_coding=True)
        exon3 = Exon("EXON003", exon3_coord, parent_id="TRANS001", phase=0, is_coding=True)

        # Build hierarchy
        gene.add_child("TRANS001")
        trans.add_child("EXON001")
        trans.add_child("EXON002")
        trans.add_child("EXON003")

        # Verify hierarchy
        assert gene.num_children == 1
        assert "TRANS001" in gene.children
        assert trans.num_children == 3
        assert trans.parent_id == gene.feature_id
        assert exon1.parent_id == trans.feature_id
        assert exon2.parent_id == trans.feature_id
        assert exon3.parent_id == trans.feature_id

    def test_gene_with_multiple_transcripts(self):
        """Test gene with multiple transcript isoforms."""
        gene_coord = GenomicCoordinate("chr1", 1000, 5000, '+', '1-based')
        gene = Gene("GENE001", gene_coord)

        # Create two transcripts
        trans1 = Transcript("TRANS001", gene_coord, parent_id="GENE001")
        trans2 = Transcript("TRANS002", gene_coord, parent_id="GENE001")

        # Add to gene
        gene.add_child("TRANS001")
        gene.add_child("TRANS002")

        assert gene.num_children == 2
        assert trans1.gene_id == gene.feature_id
        assert trans2.gene_id == gene.feature_id


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
