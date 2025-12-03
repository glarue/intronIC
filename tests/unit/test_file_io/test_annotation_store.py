"""
Tests for StreamingAnnotationStore.

Tests the SQLite-backed annotation storage for memory-efficient streaming mode.
"""

import pytest

from intronIC.file_io.annotation_store import StreamingAnnotationStore


class TestStreamingAnnotationStore:
    """Tests for StreamingAnnotationStore class."""

    @pytest.fixture
    def sample_gff3(self, tmp_path):
        """Create a sample GFF3 file for testing."""
        gff3_content = """\
##gff-version 3
chr1\t.\tgene\t1000\t5000\t.\t+\t.\tID=gene1;Name=TestGene1
chr1\t.\tmRNA\t1000\t5000\t.\t+\t.\tID=trans1;Parent=gene1
chr1\t.\texon\t1000\t1200\t.\t+\t.\tID=exon1;Parent=trans1
chr1\t.\texon\t2000\t2500\t.\t+\t.\tID=exon2;Parent=trans1
chr1\t.\texon\t4500\t5000\t.\t+\t.\tID=exon3;Parent=trans1
chr2\t.\tgene\t500\t2000\t.\t-\t.\tID=gene2;Name=TestGene2
chr2\t.\tmRNA\t500\t2000\t.\t-\t.\tID=trans2;Parent=gene2
chr2\t.\texon\t500\t800\t.\t-\t.\tID=exon4;Parent=trans2
chr2\t.\texon\t1500\t2000\t.\t-\t.\tID=exon5;Parent=trans2
chr3\t.\tgene\t100\t300\t.\t+\t.\tID=gene3;Name=TestGene3
"""
        gff3_path = tmp_path / "test.gff3"
        gff3_path.write_text(gff3_content)
        return gff3_path

    @pytest.fixture
    def db_path(self, tmp_path):
        """Get a path for the test database."""
        return tmp_path / "test_annotations.db"

    def test_create_from_file(self, sample_gff3, db_path):
        """Test creating annotation store from GFF3 file."""
        store = StreamingAnnotationStore.create_from_file(
            annotation_path=sample_gff3,
            db_path=db_path,
        )

        assert db_path.exists()
        assert store.get_contig_count() == 3
        assert store.get_total_annotations() == 10

        store.cleanup()
        assert not db_path.exists()

    def test_get_contigs(self, sample_gff3, db_path):
        """Test retrieving list of contigs."""
        store = StreamingAnnotationStore.create_from_file(
            annotation_path=sample_gff3,
            db_path=db_path,
        )

        contigs = store.get_contigs()
        assert contigs == ["chr1", "chr2", "chr3"]

        store.cleanup()

    def test_get_contigs_with_counts(self, sample_gff3, db_path):
        """Test retrieving contigs with their annotation counts."""
        store = StreamingAnnotationStore.create_from_file(
            annotation_path=sample_gff3,
            db_path=db_path,
        )

        contigs_with_counts = store.get_contigs_with_counts()
        # chr1: gene, mRNA, 3 exons = 5
        # chr2: gene, mRNA, 2 exons = 4
        # chr3: gene = 1
        assert contigs_with_counts == [("chr1", 5), ("chr2", 4), ("chr3", 1)]

        # Verify counts sum to total
        total = sum(count for _, count in contigs_with_counts)
        assert total == store.get_total_annotations()

        store.cleanup()

    def test_get_annotations_for_contig(self, sample_gff3, db_path):
        """Test retrieving annotations for a specific contig."""
        store = StreamingAnnotationStore.create_from_file(
            annotation_path=sample_gff3,
            db_path=db_path,
        )

        # Get chr1 annotations
        chr1_anns = store.get_annotations_for_contig("chr1")
        assert len(chr1_anns) == 5  # gene, mRNA, 3 exons

        # Verify they're in order by line number
        line_numbers = [a.line_number for a in chr1_anns]
        assert line_numbers == sorted(line_numbers)

        # Verify annotation data
        gene = chr1_anns[0]
        assert gene.name == "gene1"
        assert gene.feat_type == "gene"
        assert gene.start == 1000
        assert gene.stop == 5000
        assert gene.strand == "+"

        store.cleanup()

    def test_get_annotations_for_nonexistent_contig(self, sample_gff3, db_path):
        """Test retrieving annotations for a contig that doesn't exist."""
        store = StreamingAnnotationStore.create_from_file(
            annotation_path=sample_gff3,
            db_path=db_path,
        )

        # Should return empty list, not error
        chr99_anns = store.get_annotations_for_contig("chr99")
        assert chr99_anns == []

        store.cleanup()

    def test_iter_annotations_for_contig(self, sample_gff3, db_path):
        """Test iterating over annotations for a contig."""
        store = StreamingAnnotationStore.create_from_file(
            annotation_path=sample_gff3,
            db_path=db_path,
        )

        # Use iterator
        chr2_anns = list(store.iter_annotations_for_contig("chr2"))
        assert len(chr2_anns) == 4  # gene, mRNA, 2 exons

        # Verify it's the same as get_annotations_for_contig
        chr2_anns_list = store.get_annotations_for_contig("chr2")
        assert len(chr2_anns) == len(chr2_anns_list)

        store.cleanup()

    def test_context_manager(self, sample_gff3, db_path):
        """Test using store as context manager."""
        with StreamingAnnotationStore.create_from_file(
            annotation_path=sample_gff3,
            db_path=db_path,
        ) as store:
            assert store.get_contig_count() == 3

        # Connection should be closed but file still exists
        assert db_path.exists()

        # Cleanup manually
        store.cleanup()

    def test_file_exists_error(self, sample_gff3, db_path):
        """Test that creating store fails if database already exists."""
        # Create first store
        store1 = StreamingAnnotationStore.create_from_file(
            annotation_path=sample_gff3,
            db_path=db_path,
        )
        store1.close()

        # Try to create another at same path - should fail
        with pytest.raises(FileExistsError):
            StreamingAnnotationStore.create_from_file(
                annotation_path=sample_gff3,
                db_path=db_path,
            )

        store1.cleanup()

    def test_annotation_not_found_error(self, db_path):
        """Test that creating store fails if annotation file doesn't exist."""
        with pytest.raises(FileNotFoundError):
            StreamingAnnotationStore.create_from_file(
                annotation_path="/nonexistent/file.gff3",
                db_path=db_path,
            )

    def test_parent_serialization(self, tmp_path, db_path):
        """Test that parent lists are correctly serialized/deserialized."""
        # Create GFF3 with feature that has parent
        gff3_content = """\
##gff-version 3
chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1
chr1\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=trans1;Parent=gene1
chr1\t.\texon\t1000\t1500\t.\t+\t.\tID=exon1;Parent=trans1
"""
        gff3_path = tmp_path / "parent_test.gff3"
        gff3_path.write_text(gff3_content)

        store = StreamingAnnotationStore.create_from_file(
            annotation_path=gff3_path,
            db_path=db_path,
        )

        anns = store.get_annotations_for_contig("chr1")

        # Check parent of mRNA (biogl normalizes mRNA to "transcript")
        transcript = [a for a in anns if a.feat_type == "transcript"][0]
        assert transcript.parent == ["gene1"]

        # Check parent of exon
        exon = [a for a in anns if a.feat_type == "exon"][0]
        assert exon.parent == ["trans1"]

        store.cleanup()
