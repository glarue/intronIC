"""
Unit tests for Intron composition classes.

Tests:
- IntronScores: Scoring data
- IntronSequences: Sequence data
- IntronMetadata: Metadata and tags
- Intron: Main class with composition

Author: intronIC refactoring project
Date: 2025-11-02
"""

import pytest
from intronIC.utils.coordinates import GenomicCoordinate
from intronIC.core.models import Exon
from intronIC.core.intron import (
    Intron,
    IntronScores,
    IntronSequences,
    IntronMetadata,
)


class TestIntronScores:
    """Test IntronScores dataclass."""

    def test_create_empty_scores(self):
        """Test creating empty scores object."""
        scores = IntronScores()

        assert scores.svm_score is None
        assert scores.five_z_score is None
        assert scores.has_all_scores() == False

    def test_create_complete_scores(self):
        """Test creating complete scores object."""
        scores = IntronScores(
            five_raw_score=12.5,
            bp_raw_score=8.3,
            three_raw_score=15.2,
            five_z_score=2.1,
            bp_z_score=1.8,
            three_z_score=2.5,
            svm_score=95.5,
            relative_score=5.5,
            decision_distance=1.2
        )

        assert scores.svm_score == 95.5
        assert scores.five_z_score == 2.1
        assert scores.has_all_scores() == True

    def test_is_high_confidence(self):
        """Test high confidence checking."""
        scores_high = IntronScores(svm_score=95.0)
        scores_low = IntronScores(svm_score=85.0)
        scores_none = IntronScores()

        assert scores_high.is_high_confidence(threshold=90.0) == True
        assert scores_low.is_high_confidence(threshold=90.0) == False
        assert scores_none.is_high_confidence(threshold=90.0) == False

    def test_has_all_scores(self):
        """Test has_all_scores method."""
        # Complete scores
        scores_complete = IntronScores(
            five_raw_score=1.0,
            bp_raw_score=2.0,
            three_raw_score=3.0,
            five_z_score=1.0,
            bp_z_score=2.0,
            three_z_score=3.0,
            svm_score=90.0
        )
        assert scores_complete.has_all_scores() == True

        # Missing some scores
        scores_partial = IntronScores(
            five_raw_score=1.0,
            svm_score=90.0
        )
        assert scores_partial.has_all_scores() == False

    def test_scores_immutable(self):
        """Test that IntronScores is immutable (frozen)."""
        scores = IntronScores(svm_score=95.0)

        with pytest.raises(Exception):  # FrozenInstanceError
            scores.svm_score = 80.0


class TestIntronSequences:
    """Test IntronSequences dataclass."""

    def test_create_empty_sequences(self):
        """Test creating empty sequences object."""
        seqs = IntronSequences()

        assert seqs.seq is None
        assert seqs.five_seq is None
        assert seqs.has_sequences() == False

    def test_create_complete_sequences(self):
        """Test creating complete sequences object."""
        seqs = IntronSequences(
            seq="GTAAGT" + "A" * 100 + "TTTAG",
            five_seq="GTAAGT",
            three_seq="TTTAG",
            bp_seq="TACTAAC",
            bp_region_seq="ATACTAACTA",
            upstream_flank="AGGCT",
            downstream_flank="CATGG"
        )

        assert seqs.seq is not None
        assert seqs.five_seq == "GTAAGT"
        assert seqs.three_seq == "TTTAG"
        assert seqs.has_sequences() == True
        assert seqs.has_flanks() == True

    def test_terminal_dinucleotides_canonical(self):
        """Test terminal dinucleotides extraction for canonical introns."""
        seqs = IntronSequences(seq="GTAAGTCCCCTTTAG")

        assert seqs.terminal_dinucleotides == "GT-AG"

    def test_terminal_dinucleotides_noncanonical(self):
        """Test terminal dinucleotides for non-canonical introns."""
        seqs_atac = IntronSequences(seq="ATATCCCCCCCAC")
        seqs_gcag = IntronSequences(seq="GCAAGTCCCCCAG")

        assert seqs_atac.terminal_dinucleotides == "AT-AC"
        assert seqs_gcag.terminal_dinucleotides == "GC-AG"

    def test_terminal_dinucleotides_none_when_no_seq(self):
        """Test terminal dinucleotides returns None when sequence absent."""
        seqs = IntronSequences()

        assert seqs.terminal_dinucleotides is None

    def test_terminal_dinucleotides_short_sequence(self):
        """Test terminal dinucleotides returns None for very short sequences."""
        seqs = IntronSequences(seq="ATG")  # Only 3 bases

        assert seqs.terminal_dinucleotides is None

    def test_sequences_immutable(self):
        """Test that IntronSequences is immutable (frozen)."""
        seqs = IntronSequences(seq="GTAAGT")

        with pytest.raises(Exception):  # FrozenInstanceError
            seqs.seq = "ATCGAT"


class TestIntronMetadata:
    """Test IntronMetadata dataclass."""

    def test_create_empty_metadata(self):
        """Test creating empty metadata object."""
        meta = IntronMetadata()

        assert meta.parent is None
        assert meta.grandparent is None
        assert meta.type_id == 'unknown'
        assert meta.is_omitted() == False
        assert meta.is_duplicate() == False

    def test_create_complete_metadata(self):
        """Test creating complete metadata object."""
        meta = IntronMetadata(
            parent="TRANS001",
            grandparent="GENE001",
            index=1,
            family_size=5,
            type_id='u12'
        )
        # Set flag properties after creation
        meta.noncanonical = False
        meta.longest_isoform = True

        assert meta.parent == "TRANS001"
        assert meta.grandparent == "GENE001"
        assert meta.index == 1
        assert meta.family_size == 5
        assert meta.type_id == 'u12'
        assert meta.longest_isoform == True

    def test_is_omitted(self):
        """Test omission status checking."""
        meta_omitted = IntronMetadata(omitted='s')  # short
        meta_not_omitted = IntronMetadata()

        assert meta_omitted.is_omitted() == True
        assert meta_not_omitted.is_omitted() == False

    def test_is_duplicate(self):
        """Test duplicate status checking."""
        meta_dup = IntronMetadata(duplicate="intron_1")
        meta_not_dup = IntronMetadata()

        assert meta_dup.is_duplicate() == True
        assert meta_not_dup.is_duplicate() == False

    def test_is_canonical(self):
        """Test canonical status checking."""
        meta_canonical = IntronMetadata()
        meta_canonical.noncanonical = False

        meta_noncanonical = IntronMetadata()
        meta_noncanonical.noncanonical = True

        assert meta_canonical.is_canonical() == True
        assert meta_noncanonical.is_canonical() == False

    def test_fractional_position(self):
        """Test fractional position storage."""
        # fractional_position is now a stored field, not computed
        # It's calculated during intron generation based on cumulative exon lengths

        # First intron (would be 0.0 if evenly spaced)
        meta1 = IntronMetadata(index=1, family_size=5, fractional_position=0.0)
        assert meta1.fractional_position == 0.0

        # Middle intron (could be any value based on actual exon lengths)
        meta3 = IntronMetadata(index=3, family_size=5, fractional_position=0.52)
        assert meta3.fractional_position == 0.52

        # Last intron (would be 1.0 for last intron... but we use cumulative before last exon)
        # So last intron is typically less than 1.0
        meta5 = IntronMetadata(index=5, family_size=5, fractional_position=0.95)
        assert meta5.fractional_position == 0.95

        # Unset position
        meta_unknown = IntronMetadata()
        assert meta_unknown.fractional_position is None

    def test_metadata_mutable(self):
        """Test that IntronMetadata is mutable (not frozen)."""
        meta = IntronMetadata()

        # Should be able to modify
        meta.type_id = 'u12'
        meta.omitted = 's'
        meta.duplicate = "other_intron"

        assert meta.type_id == 'u12'
        assert meta.omitted == 's'
        assert meta.duplicate == "other_intron"


class TestIntron:
    """Test Intron main class."""

    def test_create_basic_intron(self):
        """Test creating a basic intron with just coordinates."""
        coord = GenomicCoordinate("chr1", 1001, 2000, '+', '1-based')
        intron = Intron("intron_1", coord)

        assert intron.intron_id == "intron_1"
        assert intron.chromosome == "chr1"
        assert intron.start == 1001
        assert intron.stop == 2000
        assert intron.strand == '+'
        assert intron.length == 1000
        assert intron.has_scores == False
        assert intron.has_sequences == False

    def test_intron_with_scores(self):
        """Test intron with scoring data."""
        coord = GenomicCoordinate("chr1", 1001, 2000, '+', '1-based')
        scores = IntronScores(
            five_z_score=2.0,
            bp_z_score=1.5,
            three_z_score=2.2,
            five_raw_score=10.0,
            bp_raw_score=8.0,
            three_raw_score=12.0,
            svm_score=95.0
        )
        intron = Intron("intron_1", coord, scores=scores)

        assert intron.has_scores == True
        assert intron.svm_score == 95.0

    def test_intron_with_sequences(self):
        """Test intron with sequence data."""
        coord = GenomicCoordinate("chr1", 1001, 2000, '+', '1-based')
        seqs = IntronSequences(
            seq="GTAAGT" + "A" * 988 + "TTTAG",
            five_seq="GTAAGT",
            three_seq="TTTAG"
        )
        intron = Intron("intron_1", coord, sequences=seqs)

        assert intron.has_sequences == True
        assert intron.terminal_dinucleotides == "GT-AG"

    def test_intron_with_metadata(self):
        """Test intron with metadata."""
        coord = GenomicCoordinate("chr1", 1001, 2000, '+', '1-based')
        meta = IntronMetadata(
            parent="TRANS001",
            grandparent="GENE001",
            type_id='u12'
        )
        intron = Intron("intron_1", coord, metadata=meta)

        assert intron.has_metadata == True
        assert intron.type_id == 'u12'

    def test_intron_from_exon_pair(self):
        """Test creating intron from adjacent exons."""
        coord1 = GenomicCoordinate("chr1", 1000, 1200, '+', '1-based')
        coord2 = GenomicCoordinate("chr1", 1500, 1700, '+', '1-based')
        exon1 = Exon("exon1", coord1, parent_id="trans1")
        exon2 = Exon("exon2", coord2, parent_id="trans1")

        intron = Intron.from_exon_pair(exon1, exon2, "intron_1")

        assert intron.start == 1201  # exon1.stop + 1
        assert intron.stop == 1499    # exon2.start - 1
        assert intron.length == 299
        assert intron.metadata.parent == "trans1"

    def test_intron_from_exon_pair_different_chromosomes(self):
        """Test that from_exon_pair rejects exons on different chromosomes."""
        coord1 = GenomicCoordinate("chr1", 1000, 1200, '+', '1-based')
        coord2 = GenomicCoordinate("chr2", 1500, 1700, '+', '1-based')
        exon1 = Exon("exon1", coord1)
        exon2 = Exon("exon2", coord2)

        with pytest.raises(ValueError, match="different chromosomes"):
            Intron.from_exon_pair(exon1, exon2)

    def test_intron_from_exon_pair_different_strands(self):
        """Test that from_exon_pair rejects exons on different strands."""
        coord1 = GenomicCoordinate("chr1", 1000, 1200, '+', '1-based')
        coord2 = GenomicCoordinate("chr1", 1500, 1700, '-', '1-based')
        exon1 = Exon("exon1", coord1)
        exon2 = Exon("exon2", coord2)

        with pytest.raises(ValueError, match="different strands"):
            Intron.from_exon_pair(exon1, exon2)

    def test_intron_from_exon_pair_overlapping(self):
        """Test that from_exon_pair rejects overlapping exons."""
        coord1 = GenomicCoordinate("chr1", 1000, 1500, '+', '1-based')
        coord2 = GenomicCoordinate("chr1", 1400, 1700, '+', '1-based')
        exon1 = Exon("exon1", coord1)
        exon2 = Exon("exon2", coord2)

        with pytest.raises(ValueError, match="overlap"):
            Intron.from_exon_pair(exon1, exon2)

    def test_with_scores(self):
        """Test with_scores method creates new intron."""
        coord = GenomicCoordinate("chr1", 1001, 2000, '+', '1-based')
        intron1 = Intron("intron_1", coord)

        scores = IntronScores(svm_score=95.0)
        intron2 = intron1.with_scores(scores)

        # Original unchanged
        assert intron1.scores is None
        assert intron1.svm_score is None

        # New intron has scores
        assert intron2.scores is not None
        assert intron2.svm_score == 95.0

    def test_with_sequences(self):
        """Test with_sequences method creates new intron."""
        coord = GenomicCoordinate("chr1", 1001, 2000, '+', '1-based')
        intron1 = Intron("intron_1", coord)

        seqs = IntronSequences(seq="GTAAGT")
        intron2 = intron1.with_sequences(seqs)

        # Original unchanged
        assert intron1.sequences is None
        # New intron has sequences
        assert intron2.sequences is not None

    def test_with_metadata(self):
        """Test with_metadata method creates new intron."""
        coord = GenomicCoordinate("chr1", 1001, 2000, '+', '1-based')
        intron1 = Intron("intron_1", coord)

        meta = IntronMetadata(parent="TRANS001")
        intron2 = intron1.with_metadata(meta)

        # Original unchanged
        assert intron1.metadata is None
        # New intron has metadata
        assert intron2.metadata is not None
        assert intron2.metadata.parent == "TRANS001"

    def test_intron_immutable(self):
        """Test that Intron is immutable (frozen)."""
        coord = GenomicCoordinate("chr1", 1001, 2000, '+', '1-based')
        intron = Intron("intron_1", coord)

        with pytest.raises(Exception):  # FrozenInstanceError
            intron.intron_id = "new_id"

    def test_empty_intron_id_rejected(self):
        """Test that empty intron_id is rejected."""
        coord = GenomicCoordinate("chr1", 1001, 2000, '+', '1-based')

        with pytest.raises(ValueError, match="intron_id cannot be empty"):
            Intron("", coord)

    def test_requires_1based_coordinates(self):
        """Test that only 1-based coordinates are accepted."""
        coord = GenomicCoordinate("chr1", 1000, 2000, '+', '0-based')

        with pytest.raises(ValueError, match="requires 1-based coordinates"):
            Intron("intron_1", coord)

    def test_clear_sequences(self):
        """Test clear_sequences clears large fields but preserves scoring sequences."""
        coord = GenomicCoordinate("chr1", 1001, 2000, '+', '1-based')
        seqs = IntronSequences(
            seq="GTAAGT" + "A" * 988 + "TTTAG",
            five_seq="GTAAGT",
            three_seq="TTTAG",
            bp_seq="TACTAAC",
            bp_seq_u2="ATCCTTTT",
            bp_region_seq="ATACTAACTA",
            upstream_flank="AGGCT",
            downstream_flank="CATGG",
            five_display_seq="GTAAGTAAAA",
            three_display_seq="ATTTAG",
            five_prime_dnt="GT",
            three_prime_dnt="AG",
            bp_relative_coords=(5, 12)
        )
        intron = Intron("intron_1", coord, sequences=seqs)

        cleared = intron.clear_sequences()

        # Large sequences should be cleared
        assert cleared.sequences.seq is None
        assert cleared.sequences.upstream_flank is None
        assert cleared.sequences.downstream_flank is None
        assert cleared.sequences.bp_region_seq is None
        assert cleared.sequences.five_display_seq is None
        assert cleared.sequences.three_display_seq is None

        # Small scoring sequences should be preserved
        assert cleared.sequences.five_seq == "GTAAGT"
        assert cleared.sequences.three_seq == "TTTAG"
        assert cleared.sequences.bp_seq == "TACTAAC"
        assert cleared.sequences.bp_seq_u2 == "ATCCTTTT"
        assert cleared.sequences.bp_relative_coords == (5, 12)

        # Dinucleotides should be preserved
        assert cleared.sequences.five_prime_dnt == "GT"
        assert cleared.sequences.three_prime_dnt == "AG"

        # Original should be unchanged (immutability)
        assert intron.sequences.seq is not None
        assert intron.sequences.upstream_flank == "AGGCT"

    def test_clear_sequences_when_no_sequences(self):
        """Test clear_sequences returns same intron when no sequences present."""
        coord = GenomicCoordinate("chr1", 1001, 2000, '+', '1-based')
        intron = Intron("intron_1", coord)

        cleared = intron.clear_sequences()

        assert cleared is intron  # Should return same object
        assert cleared.sequences is None

    def test_clear_all_sequences(self):
        """Test clear_all_sequences clears all sequence fields."""
        coord = GenomicCoordinate("chr1", 1001, 2000, '+', '1-based')
        seqs = IntronSequences(
            seq="GTAAGT" + "A" * 988 + "TTTAG",
            five_seq="GTAAGT",
            three_seq="TTTAG",
            bp_seq="TACTAAC",
            bp_region_seq="ATACTAACTA",
            upstream_flank="AGGCT",
            downstream_flank="CATGG"
        )
        intron = Intron("intron_1", coord, sequences=seqs)

        cleared = intron.clear_all_sequences()

        # All sequences should be cleared
        assert cleared.sequences is None
        assert cleared.has_sequences == False

        # Original should be unchanged (immutability)
        assert intron.sequences is not None
        assert intron.sequences.seq is not None

    def test_clear_all_sequences_when_no_sequences(self):
        """Test clear_all_sequences returns same intron when no sequences present."""
        coord = GenomicCoordinate("chr1", 1001, 2000, '+', '1-based')
        intron = Intron("intron_1", coord)

        cleared = intron.clear_all_sequences()

        assert cleared is intron  # Should return same object
        assert cleared.sequences is None

    def test_clear_sequences_preserves_other_components(self):
        """Test that clearing sequences preserves scores and metadata."""
        coord = GenomicCoordinate("chr1", 1001, 2000, '+', '1-based')
        scores = IntronScores(svm_score=95.0)
        seqs = IntronSequences(seq="GTAAGTTTAG", five_seq="GTAAGT", three_seq="TTTAG")
        meta = IntronMetadata(parent="TRANS001", type_id='u12')

        intron = Intron("intron_1", coord, scores, seqs, meta)
        cleared = intron.clear_sequences()

        # Sequences should be partially cleared
        assert cleared.sequences.seq is None
        assert cleared.sequences.five_seq == "GTAAGT"

        # Other components should be unchanged
        assert cleared.scores is not None
        assert cleared.svm_score == 95.0
        assert cleared.metadata is not None
        assert cleared.metadata.parent == "TRANS001"
        assert cleared.type_id == 'u12'


class TestIntronIntegration:
    """Test integration scenarios with complete Intron objects."""

    def test_complete_intron_object(self):
        """Test creating a fully populated intron object."""
        coord = GenomicCoordinate("chr1", 1001, 2000, '+', '1-based')

        scores = IntronScores(
            five_raw_score=12.5,
            bp_raw_score=8.3,
            three_raw_score=15.2,
            five_z_score=2.1,
            bp_z_score=1.8,
            three_z_score=2.5,
            svm_score=95.5
        )

        seqs = IntronSequences(
            seq="GTAAGT" + "A" * 988 + "TTTAG",
            five_seq="GTAAGT",
            three_seq="TTTAG",
            bp_seq="TACTAAC",
            upstream_flank="AGGCT",
            downstream_flank="CATGG"
        )

        meta = IntronMetadata(
            parent="ENST001",
            grandparent="ENSG001",
            index=2,
            family_size=5,
            type_id='u12'
        )
        meta.longest_isoform = True

        intron = Intron("intron_chr1_1001_2000", coord, scores, seqs, meta)

        # Verify all components
        assert intron.has_scores == True
        assert intron.has_sequences == True
        assert intron.has_metadata == True

        assert intron.svm_score == 95.5
        assert intron.terminal_dinucleotides == "GT-AG"
        assert intron.type_id == 'u12'
        assert intron.length == 1000

    def test_progressive_intron_building(self):
        """Test building intron progressively through pipeline stages."""
        # Stage 1: Extraction (coordinates only)
        coord = GenomicCoordinate("chr1", 1001, 2000, '+', '1-based')
        intron = Intron("intron_1", coord)
        assert intron.has_sequences == False
        assert intron.has_scores == False

        # Stage 2: Add sequences
        seqs = IntronSequences(
            seq="GTAAGT" + "A" * 988 + "TTTAG",
            five_seq="GTAAGT",
            three_seq="TTTAG"
        )
        intron = intron.with_sequences(seqs)
        assert intron.has_sequences == True
        assert intron.has_scores == False

        # Stage 3: Add scores
        scores = IntronScores(svm_score=95.0)
        intron = intron.with_scores(scores)
        assert intron.has_sequences == True
        assert intron.has_scores == False  # Not all scores present

        # Stage 4: Complete scores
        complete_scores = IntronScores(
            five_raw_score=12.5,
            bp_raw_score=8.3,
            three_raw_score=15.2,
            five_z_score=2.1,
            bp_z_score=1.8,
            three_z_score=2.5,
            svm_score=95.0
        )
        intron = intron.with_scores(complete_scores)
        assert intron.has_scores == True


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
