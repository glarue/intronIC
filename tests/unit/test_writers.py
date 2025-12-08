"""
Unit tests for file format writers.

Tests:
- BEDWriter: BED format output
- MetaWriter: Metadata output
- SequenceWriter: Sequence file output
- ScoreWriter: Detailed scoring output
- MappingWriter: Mapping file output

Author: intronIC refactoring project
Date: 2025-11-02
"""

from pathlib import Path

import pytest

from intronIC.core.intron import (
    Intron,
    IntronMetadata,
    IntronScores,
    IntronSequences,
    OmissionReason,
)
from intronIC.file_io.writers import (
    BEDWriter,
    MappingWriter,
    MetaWriter,
    ScoreWriter,
    SequenceWriter,
)
from intronIC.utils.coordinates import GenomicCoordinate

# ============================================================================
# Test Fixtures
# ============================================================================


@pytest.fixture
def basic_intron():
    """Create a basic intron with coordinates only."""
    coord = GenomicCoordinate("chr1", 1001, 2000, "+", "1-based")
    return Intron("intron_001", coord)


@pytest.fixture
def full_intron():
    """Create a fully populated intron with all data."""
    coord = GenomicCoordinate("chr1", 1001, 2000, "+", "1-based")

    scores = IntronScores(
        svm_score=95.5,
        relative_score=5.5,
        decision_distance=2.3,
        five_raw_score=12.345678,
        five_z_score=2.456,
        bp_raw_score=8.765432,
        bp_z_score=1.234,
        three_raw_score=15.987654,
        three_z_score=3.789,
    )

    sequences = IntronSequences(
        seq="GTAAGT" + "N" * 988 + "TTTAG",
        five_seq="AGGCTGTAAGT",
        three_seq="TTTAGCATGG",
        bp_seq="TACTAAC",
        bp_region_seq="TACTAACTATGCTATG",  # Branch point region
        upstream_flank="AGGCT",
        downstream_flank="CATGG",
    )

    metadata = IntronMetadata(
        parent="TRANS001",
        grandparent="GENE001",
        index=3,
        family_size=5,
        type_id="u12",
        phase=0,
        fractional_position=0.5,  # Set explicitly for testing
    )
    # Set flag properties after creation
    metadata.noncanonical = False
    metadata.longest_isoform = True

    return Intron("intron_full", coord, scores, sequences, metadata)


@pytest.fixture
def minimal_scored_intron():
    """Create intron with minimal score data."""
    coord = GenomicCoordinate("chr2", 5000, 6000, "-", "1-based")
    scores = IntronScores(svm_score=12.3)
    metadata = IntronMetadata(parent="TRANS002", index=1, family_size=2, type_id="u2")

    return Intron("intron_minimal", coord, scores=scores, metadata=metadata)


@pytest.fixture
def intron_with_tags():
    """Create intron with various tags (noncanonical, duplicate, etc.)."""
    coord = GenomicCoordinate("chr3", 10000, 11000, "+", "1-based")

    metadata = IntronMetadata(
        parent="TRANS003",
        index=2,
        family_size=4,
        type_id="u2",
        duplicate="intron_rep",
        omitted=OmissionReason.SHORT,
    )
    # Set flag properties after creation
    metadata.noncanonical = True
    metadata.longest_isoform = False
    metadata.corrected = True

    sequences = IntronSequences(
        seq="ATAACT" + "N" * 988 + "TTCAG",
        five_seq="AGGCTATAACT",
        three_seq="TTCAGCATGG",
    )

    return Intron("intron_tagged", coord, sequences=sequences, metadata=metadata)


# ============================================================================
# BEDWriter Tests
# ============================================================================


class TestBEDWriter:
    """Test BED format writer."""

    def test_create_writer(self, tmp_path):
        """Test creating a BED writer."""
        bed_file = tmp_path / "test.bed"
        writer = BEDWriter(bed_file)
        assert writer.file_path == bed_file
        assert writer.introns_written == 0

    def test_context_manager(self, tmp_path, basic_intron):
        """Test writer as context manager."""
        bed_file = tmp_path / "test.bed"

        with BEDWriter(bed_file) as writer:
            writer.write_intron(basic_intron)

        assert bed_file.exists()

    def test_write_basic_intron(self, tmp_path, basic_intron):
        """Test writing a basic intron."""
        bed_file = tmp_path / "test.bed"

        with BEDWriter(bed_file) as writer:
            writer.write_intron(basic_intron)

        lines = bed_file.read_text().strip().split("\n")
        assert len(lines) == 1

        fields = lines[0].split("\t")
        assert len(fields) == 7  # Now includes attributes column
        assert fields[0] == "chr1"
        assert fields[1] == "1000"  # 0-based start
        assert fields[2] == "2000"  # 1-based stop
        assert fields[4] == "NA"  # no score
        assert fields[5] == "+"
        assert fields[6] == "NA"  # no attributes for basic intron

    def test_write_scored_intron(self, tmp_path, full_intron):
        """Test writing an intron with SVM score."""
        bed_file = tmp_path / "test.bed"

        with BEDWriter(bed_file) as writer:
            writer.write_intron(full_intron)

        lines = bed_file.read_text().strip().split("\n")
        fields = lines[0].split("\t")

        assert fields[4] == "95.5"  # SVM score

    def test_write_with_species_name(self, tmp_path, full_intron):
        """Test writing with species name prefix."""
        bed_file = tmp_path / "test.bed"

        with BEDWriter(bed_file) as writer:
            writer.write_intron(full_intron, species_name="homo_sapiens")

        lines = bed_file.read_text().strip().split("\n")
        fields = lines[0].split("\t")

        # Name should include species prefix - check name field contains expected parts
        name_field = fields[3]
        assert "homo_sapiens" in name_field or "TRANS001" in name_field

    def test_write_simple_name(self, tmp_path, full_intron):
        """Test writing with simple naming (no species)."""
        bed_file = tmp_path / "test.bed"

        with BEDWriter(bed_file) as writer:
            writer.write_intron(full_intron, simple_name=True)

        lines = bed_file.read_text().strip().split("\n")
        fields = lines[0].split("\t")

        # With simple_name=True and no intron_number, it uses intron_id
        # Format: {species_prefix}-i{intron_id}
        name_field = fields[3]
        assert "-i" in name_field  # Simple naming uses -i prefix
        assert "intron_full" in name_field or "XXXXXX" in name_field

    def test_write_intron_with_tags(self, tmp_path, intron_with_tags):
        """Test writing intron with tags and attributes."""
        bed_file = tmp_path / "test.bed"

        with BEDWriter(bed_file) as writer:
            writer.write_intron(intron_with_tags)

        lines = bed_file.read_text().strip().split("\n")
        fields = lines[0].split("\t")

        # Omit tag should be present in name (only [o:s] is in name)
        name_field = fields[3]
        assert "[o:s]" in name_field  # omitted:short

        # Verbose attributes should be in attributes column
        attrs_field = fields[6]
        assert "noncanonical" in attrs_field
        assert "not_longest_isoform" in attrs_field
        assert "duplicate" in attrs_field
        assert "corrected" in attrs_field
        assert "omitted_short" in attrs_field  # omitted='s' becomes omitted_short

    def test_write_multiple_introns(self, tmp_path, basic_intron, full_intron):
        """Test writing multiple introns."""
        bed_file = tmp_path / "test.bed"

        with BEDWriter(bed_file) as writer:
            count = writer.write_introns([basic_intron, full_intron])

        assert count == 2
        lines = bed_file.read_text().strip().split("\n")
        assert len(lines) == 2

    def test_write_negative_strand(self, tmp_path, minimal_scored_intron):
        """Test writing negative strand intron."""
        bed_file = tmp_path / "test.bed"

        with BEDWriter(bed_file) as writer:
            writer.write_intron(minimal_scored_intron)

        lines = bed_file.read_text().strip().split("\n")
        fields = lines[0].split("\t")

        assert fields[5] == "-"

    def test_0based_coordinates(self, tmp_path, basic_intron):
        """Test that start coordinate is 0-based."""
        bed_file = tmp_path / "test.bed"

        with BEDWriter(bed_file) as writer:
            writer.write_intron(basic_intron)

        lines = bed_file.read_text().strip().split("\n")
        fields = lines[0].split("\t")

        # Start should be 1001-1 = 1000 (0-based)
        assert fields[1] == "1000"
        # Stop should remain 2000 (1-based inclusive)
        assert fields[2] == "2000"


# ============================================================================
# MetaWriter Tests
# ============================================================================


class TestMetaWriter:
    """Test metadata format writer."""

    def test_create_writer(self, tmp_path):
        """Test creating a metadata writer."""
        meta_file = tmp_path / "test.meta.iic"
        writer = MetaWriter(meta_file)
        assert writer.file_path == meta_file
        assert writer.introns_written == 0

    def test_write_header(self, tmp_path):
        """Test writing metadata header."""
        meta_file = tmp_path / "test.meta.iic"

        with MetaWriter(meta_file) as writer:
            writer.write_header()

        lines = meta_file.read_text().strip().split("\n")
        assert len(lines) == 1

        header = lines[0].split("\t")
        assert "name" in header
        assert "rel_score" in header
        assert "dnts" in header
        assert "length" in header
        assert "parent" in header
        assert "type_id" in header
        assert "attributes" in header  # New column

    def test_write_basic_intron(self, tmp_path, basic_intron):
        """Test writing a basic intron metadata."""
        meta_file = tmp_path / "test.meta.iic"

        with MetaWriter(meta_file) as writer:
            writer.write_header()
            writer.write_intron(basic_intron)

        lines = meta_file.read_text().strip().split("\n")
        assert len(lines) == 2  # header + 1 intron

        fields = lines[1].split("\t")
        assert len(fields) == 15  # 15 metadata fields (including attributes)

    def test_write_full_intron(self, tmp_path, full_intron):
        """Test writing fully populated intron."""
        meta_file = tmp_path / "test.meta.iic"

        with MetaWriter(meta_file) as writer:
            writer.write_header()
            writer.write_intron(full_intron)

        lines = meta_file.read_text().strip().split("\n")
        fields = lines[1].split("\t")

        # Check key fields
        assert "TRANS001_3" in fields[0]  # name
        assert fields[1] == "5.5"  # relative score
        assert "GT-AG" in fields[2]  # terminal dinucleotides
        assert fields[5] == "1000"  # length
        assert fields[6] == "TRANS001"  # parent
        assert fields[7] == "GENE001"  # grandparent
        assert fields[8] == "3"  # index
        assert fields[9] == "5"  # family_size
        assert fields[12] == "u12"  # type_id

    def test_null_values(self, tmp_path, basic_intron):
        """Test that missing values are replaced with null placeholder."""
        meta_file = tmp_path / "test.meta.iic"

        with MetaWriter(meta_file) as writer:
            writer.write_intron(basic_intron)

        lines = meta_file.read_text().strip().split("\n")
        fields = lines[0].split("\t")

        # Many fields should be 'NA' for basic intron
        assert fields[1] == "NA"  # rel_score (no scores)
        assert fields[6] == "NA"  # parent (no metadata)

    def test_fractional_position(self, tmp_path, full_intron):
        """Test fractional position output in meta.iic."""
        meta_file = tmp_path / "test.meta.iic"

        with MetaWriter(meta_file) as writer:
            writer.write_intron(full_intron)

        lines = meta_file.read_text().strip().split("\n")
        fields = lines[0].split("\t")

        # Fractional position is stored field, set to 0.5 in fixture
        assert fields[10] == "0.5"

    def test_write_multiple_introns(self, tmp_path, basic_intron, full_intron):
        """Test writing multiple introns."""
        meta_file = tmp_path / "test.meta.iic"

        with MetaWriter(meta_file) as writer:
            writer.write_header()
            count = writer.write_introns([basic_intron, full_intron])

        assert count == 2
        lines = meta_file.read_text().strip().split("\n")
        assert len(lines) == 3  # header + 2 introns


# ============================================================================
# SequenceWriter Tests
# ============================================================================


class TestSequenceWriter:
    """Test sequence format writer."""

    def test_create_writer(self, tmp_path):
        """Test creating a sequence writer."""
        seq_file = tmp_path / "test.introns.iic"
        writer = SequenceWriter(seq_file)
        assert writer.file_path == seq_file

    def test_write_intron_with_sequences(self, tmp_path, full_intron):
        """Test writing intron with sequences."""
        seq_file = tmp_path / "test.introns.iic"

        with SequenceWriter(seq_file) as writer:
            writer.write_intron(full_intron, include_score=True)

        lines = seq_file.read_text().strip().split("\n")
        assert len(lines) == 1

        fields = lines[0].split("\t")
        assert len(fields) == 5  # name, score, upstream, seq, downstream
        assert fields[1] == "95.5"  # score
        assert fields[2] == "AGGCT"  # upstream flank
        assert "GTAAGT" in fields[3]  # sequence starts with donor
        assert fields[4] == "CATGG"  # downstream flank

    def test_write_without_score(self, tmp_path, full_intron):
        """Test writing sequences without score column."""
        seq_file = tmp_path / "test.introns.iic"

        with SequenceWriter(seq_file) as writer:
            writer.write_intron(full_intron, include_score=False)

        lines = seq_file.read_text().strip().split("\n")
        fields = lines[0].split("\t")

        assert len(fields) == 4  # name, upstream, seq, downstream (no score)
        assert fields[1] == "AGGCT"  # upstream flank is now field 1

    def test_write_intron_without_sequences_fails(self, tmp_path, basic_intron):
        """Test that writing intron without sequences raises error."""
        seq_file = tmp_path / "test.introns.iic"

        with pytest.raises(ValueError, match="no sequence data"):
            with SequenceWriter(seq_file) as writer:
                writer.write_intron(basic_intron)

    def test_empty_flanks(self, tmp_path):
        """Test writing sequences with empty flanks."""
        seq_file = tmp_path / "test.introns.iic"

        # Create intron with sequence but no flanks (None values become empty strings)
        coord = GenomicCoordinate("chr1", 1000, 2000, "+", "1-based")
        sequences = IntronSequences(seq="GTAAGTNNNNTTTAG")
        intron = Intron("test", coord, sequences=sequences)

        with SequenceWriter(seq_file) as writer:
            writer.write_intron(intron, include_score=False)

        # Read raw line
        lines = seq_file.read_text().strip().split("\n")
        line = lines[0]

        fields = line.split("\t")
        assert (
            len(fields) >= 3
        )  # At minimum: name, upstream_flank, seq (downstream might be empty)
        assert fields[1] == ""  # upstream flank should be empty
        assert fields[2] == "GTAAGTNNNNTTTAG"  # sequence
        # fields[3] may or may not exist depending on trailing tab handling

    def test_write_multiple_introns(self, tmp_path, full_intron):
        """Test writing multiple introns."""
        seq_file = tmp_path / "test.introns.iic"

        # Create second intron
        coord2 = GenomicCoordinate("chr2", 5000, 6000, "-", "1-based")
        seqs2 = IntronSequences(
            seq="GTAAGT" + "N" * 988 + "TTTAG",
            upstream_flank="AGGCT",
            downstream_flank="CATGG",
        )
        intron2 = Intron("intron2", coord2, sequences=seqs2)

        with SequenceWriter(seq_file) as writer:
            count = writer.write_introns([full_intron, intron2])

        assert count == 2
        lines = seq_file.read_text().strip().split("\n")
        assert len(lines) == 2

    def test_null_score_when_unavailable(self, tmp_path):
        """Test that score is 'NA' when unavailable."""
        seq_file = tmp_path / "test.introns.iic"

        coord = GenomicCoordinate("chr1", 1000, 2000, "+", "1-based")
        sequences = IntronSequences(seq="GTAAGTNNNNTTTAG")
        intron = Intron("test", coord, sequences=sequences)

        with SequenceWriter(seq_file) as writer:
            writer.write_intron(intron, include_score=True)

        lines = seq_file.read_text().strip().split("\n")
        fields = lines[0].split("\t")

        assert fields[1] == "NA"  # score is null


# ============================================================================
# ScoreWriter Tests
# ============================================================================


class TestScoreWriter:
    """Test detailed score format writer."""

    def test_create_writer(self, tmp_path):
        """Test creating a score writer."""
        score_file = tmp_path / "test.score_info.iic"
        writer = ScoreWriter(score_file)
        assert writer.file_path == score_file

    def test_write_header(self, tmp_path):
        """Test writing score file header."""
        score_file = tmp_path / "test.score_info.iic"

        with ScoreWriter(score_file) as writer:
            writer.write_header()

        lines = score_file.read_text().strip().split("\n")
        header = lines[0].split("\t")

        assert "name" in header
        assert "svm_score" in header
        assert "5'_raw" in header
        assert "bp_z" in header
        assert "3'_seq" in header

    def test_write_full_scores(self, tmp_path, full_intron):
        """Test writing intron with all scores."""
        score_file = tmp_path / "test.score_info.iic"

        with ScoreWriter(score_file) as writer:
            writer.write_header()
            writer.write_intron(full_intron)

        lines = score_file.read_text().strip().split("\n")
        fields = lines[1].split("\t")

        assert len(fields) == 18
        # Header: name, rel_score, svm_score, 5'_seq, 5'_raw, 5'_z,
        #         bp_seq, bp_seq_u2, bp_raw, bp_z, 3'_seq, 3'_raw, 3'_z,
        #         min(5,bp), min(5,3), max(5,bp), max(5,3), decision_dist
        # Check scores are rounded correctly
        assert fields[1] == "5.5"  # rel_score
        assert fields[2] == "95.5"  # svm_score
        assert fields[4] == "12.345678"  # five_raw (index 4, after 5'_seq)
        assert fields[5] == "2.456"  # five_z (index 5)

    def test_write_partial_scores(self, tmp_path):
        """Test writing intron with minimal scores."""
        score_file = tmp_path / "test.score_info.iic"

        # Create intron with only SVM score
        coord = GenomicCoordinate("chr2", 5000, 6000, "-", "1-based")
        scores = IntronScores(svm_score=12.3, relative_score=2.3)
        metadata = IntronMetadata(
            parent="TRANS002", index=1, family_size=2, type_id="u2"
        )
        intron = Intron("intron_minimal", coord, scores=scores, metadata=metadata)

        with ScoreWriter(score_file) as writer:
            writer.write_intron(intron)

        lines = score_file.read_text().strip().split("\n")
        fields = lines[0].split("\t")

        # Check that we have all 18 fields
        assert len(fields) == 18
        # Relative and SVM scores should be available
        assert fields[1] == "2.3"  # rel_score available
        assert fields[2] == "12.3"  # svm_score available
        # PWM scores should be null
        assert fields[5] == "NA"  # five_raw not available
        assert fields[6] == "NA"  # five_z not available

    def test_null_values_for_no_scores(self, tmp_path, basic_intron):
        """Test that all score fields are null when no scores."""
        score_file = tmp_path / "test.score_info.iic"

        with ScoreWriter(score_file) as writer:
            writer.write_intron(basic_intron)

        lines = score_file.read_text().strip().split("\n")
        fields = lines[0].split("\t")

        # All score fields should be 'NA'
        for i in [1, 2, 3, 5, 6, 9, 10, 12, 13]:  # score field indices
            assert fields[i] == "NA"

    def test_write_sequences_in_scores(self, tmp_path, full_intron):
        """Test that scoring sequences are included."""
        score_file = tmp_path / "test.score_info.iic"

        with ScoreWriter(score_file) as writer:
            writer.write_intron(full_intron)

        lines = score_file.read_text().strip().split("\n")
        fields = lines[0].split("\t")

        # Header: name, rel_score, svm_score, 5'_seq, 5'_raw, 5'_z,
        #         bp_seq, bp_seq_u2, bp_raw, bp_z, 3'_seq, 3'_raw, 3'_z, ...
        assert fields[3] == "AGGCTGTAAGT"  # five_seq (index 3)
        assert fields[6] == "TACTAAC"  # bp_seq (index 6)
        # Note: bp_region_seq is now bp_seq_u2 field at index 7
        assert fields[10] == "TTTAGCATGG"  # three_seq (index 10)


# ============================================================================
# MappingWriter Tests
# ============================================================================


class TestMappingWriter:
    """Test mapping file writer."""

    def test_create_writer(self, tmp_path):
        """Test creating a mapping writer."""
        map_file = tmp_path / "test.dupe_map.iic"
        writer = MappingWriter(map_file)
        assert writer.file_path == map_file

    def test_write_single_mapping(self, tmp_path):
        """Test writing a single mapping."""
        map_file = tmp_path / "test.dupe_map.iic"

        with MappingWriter(map_file) as writer:
            writer.write_mapping("representative1", "duplicate1")

        lines = map_file.read_text().strip().split("\n")
        assert len(lines) == 1

        fields = lines[0].split("\t")
        assert len(fields) == 2
        assert fields[0] == "representative1"
        assert fields[1] == "duplicate1"

    def test_write_multiple_mappings_dict(self, tmp_path):
        """Test writing mappings from dictionary."""
        map_file = tmp_path / "test.dupe_map.iic"

        mappings = {"rep1": ["dup1", "dup2", "dup3"], "rep2": ["dup4", "dup5"]}

        with MappingWriter(map_file) as writer:
            count = writer.write_mappings(mappings)

        assert count == 5  # 3 + 2 mappings
        lines = map_file.read_text().strip().split("\n")
        assert len(lines) == 5

    def test_write_overlap_map(self, tmp_path):
        """Test writing overlap mappings."""
        map_file = tmp_path / "test.overlap_map.iic"

        overlaps = {"kept1": ["excluded1", "excluded2"]}

        with MappingWriter(map_file) as writer:
            writer.write_mappings(overlaps)

        lines = map_file.read_text().strip().split("\n")
        assert len(lines) == 2

        # Both lines should reference kept1
        assert all("kept1" in line for line in lines)

    def test_context_manager(self, tmp_path):
        """Test mapping writer as context manager."""
        map_file = tmp_path / "test.map"

        with MappingWriter(map_file) as writer:
            writer.write_mapping("a", "b")

        assert map_file.exists()
        assert writer.mappings_written == 1


# ============================================================================
# Integration Tests
# ============================================================================


class TestWriterIntegration:
    """Test interactions between multiple writers."""

    def test_all_writers_same_intron(self, tmp_path, full_intron):
        """Test that all writers can write the same intron."""
        bed_file = tmp_path / "test.bed.iic"
        meta_file = tmp_path / "test.meta.iic"
        seq_file = tmp_path / "test.introns.iic"
        score_file = tmp_path / "test.score_info.iic"

        # Write to all formats
        with BEDWriter(bed_file) as writer:
            writer.write_intron(full_intron)

        with MetaWriter(meta_file) as writer:
            writer.write_header()
            writer.write_intron(full_intron)

        with SequenceWriter(seq_file) as writer:
            writer.write_intron(full_intron)

        with ScoreWriter(score_file) as writer:
            writer.write_header()
            writer.write_intron(full_intron)

        # All files should exist
        assert bed_file.exists()
        assert meta_file.exists()
        assert seq_file.exists()
        assert score_file.exists()

    def test_consistent_naming_across_writers(self, tmp_path, full_intron):
        """Test that intron names are consistent across all writers."""
        bed_file = tmp_path / "test.bed.iic"
        meta_file = tmp_path / "test.meta.iic"
        seq_file = tmp_path / "test.introns.iic"

        species = "test_species"

        with BEDWriter(bed_file) as writer:
            writer.write_intron(full_intron, species_name=species)

        with MetaWriter(meta_file) as writer:
            writer.write_intron(full_intron, species_name=species)

        with SequenceWriter(seq_file) as writer:
            writer.write_intron(full_intron, species_name=species)

        # Extract names from each file
        bed_name = bed_file.read_text().strip().split("\t")[3].split(";")[0]
        meta_name = meta_file.read_text().strip().split("\t")[0]
        seq_name = seq_file.read_text().strip().split("\t")[0]

        # All names should be consistent with each other
        # Format is now abbreviated: TesSpe-GENE001@TRANS001_3(5)
        assert bed_name == meta_name
        assert meta_name == seq_name
        # Should contain parent info
        assert "TRANS001" in bed_name
        assert "GENE001" in bed_name
