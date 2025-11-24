"""
Unit tests for extraction/filters.py pre-filtering functions.

Tests the pre-filtering logic that determines which introns need
sequence extraction before actually extracting sequences.
"""

import pytest
from typing import Set, Tuple

from intronIC.core.intron import Intron, IntronMetadata
from intronIC.utils.coordinates import GenomicCoordinate
from intronIC.extraction.filters import (
    should_extract_sequences_for,
    prefilter_introns,
    PrefilterResult
)


class TestShouldExtractSequencesFor:
    """Tests for should_extract_sequences_for() function."""

    def create_test_intron(
        self,
        start: int = 100,
        stop: int = 200,
        grandparent: str = "gene1",
        parent: str = "transcript1"
    ) -> Intron:
        """Helper to create a test intron without sequences."""
        coords = GenomicCoordinate(
            chromosome="chr1",
            start=start,
            stop=stop,
            strand="+",
            system="1-based"
        )
        metadata = IntronMetadata(
            grandparent=grandparent,
            parent=parent
        )
        return Intron(
            intron_id=f"intron_{start}_{stop}",
            coordinates=coords,
            metadata=metadata
        )

    def test_too_short_intron_skipped(self):
        """Test that introns shorter than min_length are skipped."""
        intron = self.create_test_intron(start=100, stop=120)  # Length 20

        result = should_extract_sequences_for(
            intron=intron,
            min_length=30,
            longest_only=False,
            longest_isoforms={},
            seen_coordinates=set(),
            include_duplicates=False
        )

        assert result is False

    def test_long_enough_intron_extracted(self):
        """Test that introns meeting min_length are extracted."""
        intron = self.create_test_intron(start=100, stop=200)  # Length 100

        result = should_extract_sequences_for(
            intron=intron,
            min_length=30,
            longest_only=False,
            longest_isoforms={},
            seen_coordinates=set(),
            include_duplicates=False
        )

        assert result is True

    def test_longest_only_skips_non_longest_isoform(self):
        """Test that non-longest isoforms are skipped when longest_only=True."""
        intron = self.create_test_intron(
            grandparent="gene1",
            parent="transcript2"  # Not the longest
        )
        longest_isoforms = {"gene1": "transcript1"}  # transcript1 is longest

        result = should_extract_sequences_for(
            intron=intron,
            min_length=30,
            longest_only=True,
            longest_isoforms=longest_isoforms,
            seen_coordinates=set(),
            include_duplicates=False
        )

        assert result is False

    def test_longest_only_extracts_longest_isoform(self):
        """Test that longest isoforms are extracted when longest_only=True."""
        intron = self.create_test_intron(
            grandparent="gene1",
            parent="transcript1"  # The longest
        )
        longest_isoforms = {"gene1": "transcript1"}

        result = should_extract_sequences_for(
            intron=intron,
            min_length=30,
            longest_only=True,
            longest_isoforms=longest_isoforms,
            seen_coordinates=set(),
            include_duplicates=False
        )

        assert result is True

    def test_duplicate_skipped_when_not_including_duplicates(self):
        """Test that duplicate coordinates are skipped when include_duplicates=False."""
        intron = self.create_test_intron(start=100, stop=200)
        seen_coords: Set[Tuple[int, int]] = {(100, 200)}  # Already seen

        result = should_extract_sequences_for(
            intron=intron,
            min_length=30,
            longest_only=False,
            longest_isoforms={},
            seen_coordinates=seen_coords,
            include_duplicates=False
        )

        assert result is False

    def test_duplicate_always_skipped_for_extraction(self):
        """
        Test that duplicate coordinates are always skipped for extraction.

        The include_duplicates flag affects output, not extraction.
        Sequences are always reused for duplicates to save memory.
        """
        intron = self.create_test_intron(start=100, stop=200)
        seen_coords: Set[Tuple[int, int]] = {(100, 200)}  # Already seen

        # With include_duplicates=False
        result_without = should_extract_sequences_for(
            intron=intron,
            min_length=30,
            longest_only=False,
            longest_isoforms={},
            seen_coordinates=seen_coords,
            include_duplicates=False
        )

        # With include_duplicates=True
        result_with = should_extract_sequences_for(
            intron=intron,
            min_length=30,
            longest_only=False,
            longest_isoforms={},
            seen_coordinates=seen_coords,
            include_duplicates=True
        )

        # Both should skip extraction (reuse sequences)
        assert result_without is False
        assert result_with is False

    def test_first_occurrence_of_coordinates_extracted(self):
        """Test that first occurrence of coordinates is always extracted."""
        intron = self.create_test_intron(start=100, stop=200)
        seen_coords: Set[Tuple[int, int]] = set()  # Not seen yet

        result = should_extract_sequences_for(
            intron=intron,
            min_length=30,
            longest_only=False,
            longest_isoforms={},
            seen_coordinates=seen_coords,
            include_duplicates=False
        )

        assert result is True


class TestPrefilterIntrons:
    """Tests for prefilter_introns() function."""

    def create_test_introns(self):
        """Create a set of test introns with various characteristics."""
        introns = []

        # Intron 1: Normal, longest isoform, chr1
        coords1 = GenomicCoordinate(chromosome="chr1", start=100, stop=200, strand="+", system="1-based")
        metadata1 = IntronMetadata(
            grandparent="gene1",
            parent="transcript1",
            parent_length=1000
        )
        introns.append(Intron(
            intron_id="intron1",
            coordinates=coords1,
            metadata=metadata1
        ))

        # Intron 2: Duplicate coordinates of intron1
        coords2 = GenomicCoordinate(chromosome="chr1", start=100, stop=200, strand="+", system="1-based")
        metadata2 = IntronMetadata(
            grandparent="gene1",
            parent="transcript1",
            parent_length=1000
        )
        introns.append(Intron(
            intron_id="intron2",
            coordinates=coords2,
            metadata=metadata2
        ))

        # Intron 3: Too short
        coords3 = GenomicCoordinate(chromosome="chr1", start=300, stop=320, strand="+", system="1-based")
        metadata3 = IntronMetadata(
            grandparent="gene2",
            parent="transcript3",
            parent_length=500
        )
        introns.append(Intron(
            intron_id="intron3",
            coordinates=coords3,
            metadata=metadata3
        ))

        # Intron 4: Not longest isoform (shorter transcript from gene1)
        coords4 = GenomicCoordinate(chromosome="chr1", start=400, stop=500, strand="+", system="1-based")
        metadata4 = IntronMetadata(
            grandparent="gene1",
            parent="transcript2",  # Different transcript from same gene
            parent_length=800  # Shorter than transcript1
        )
        introns.append(Intron(
            intron_id="intron4",
            coordinates=coords4,
            metadata=metadata4
        ))

        # Intron 5: Normal, different chromosome
        coords5 = GenomicCoordinate(chromosome="chr2", start=100, stop=300, strand="+", system="1-based")
        metadata5 = IntronMetadata(
            grandparent="gene3",
            parent="transcript5",
            parent_length=1200
        )
        introns.append(Intron(
            intron_id="intron5",
            coordinates=coords5,
            metadata=metadata5
        ))

        return introns

    def test_basic_prefiltering(self):
        """Test basic pre-filtering without strict options."""
        introns = self.create_test_introns()

        result = prefilter_introns(
            introns=introns,
            min_length=30,
            longest_only=False,
            include_duplicates=False
        )

        # Should extract: intron1, intron4, intron5
        # Should skip: intron2 (duplicate), intron3 (too short)
        assert len(result.extract_list) == 3
        assert len(result.skip_list) == 2
        assert result.stats['total'] == 5
        assert result.stats['extract'] == 3
        assert result.stats['skip'] == 2
        assert result.stats['too_short'] == 1
        assert result.stats['duplicate'] == 1

    def test_longest_only_prefiltering(self):
        """Test pre-filtering with longest_only=True."""
        introns = self.create_test_introns()

        result = prefilter_introns(
            introns=introns,
            min_length=30,
            longest_only=True,
            include_duplicates=False
        )

        # Should extract: intron1, intron5 (longest isoforms only)
        # Should skip: intron2 (duplicate), intron3 (too short), intron4 (not longest)
        assert len(result.extract_list) == 2
        assert len(result.skip_list) == 3
        assert result.stats['not_longest_isoform'] == 1

    def test_include_duplicates_prefiltering(self):
        """
        Test pre-filtering with include_duplicates=True.

        Note: include_duplicates affects output, not extraction.
        Duplicates are always skipped during extraction to save memory.
        """
        introns = self.create_test_introns()

        result = prefilter_introns(
            introns=introns,
            min_length=30,
            longest_only=False,
            include_duplicates=True
        )

        # Duplicates are ALWAYS skipped for extraction (sequences reused)
        # Should extract: intron1, intron4, intron5
        # Should skip: intron2 (duplicate), intron3 (too short)
        assert len(result.extract_list) == 3
        assert len(result.skip_list) == 2

    def test_empty_input(self):
        """Test pre-filtering with empty intron list."""
        result = prefilter_introns(
            introns=[],
            min_length=30,
            longest_only=False,
            include_duplicates=False
        )

        assert len(result.extract_list) == 0
        assert len(result.skip_list) == 0
        assert result.stats['total'] == 0

    def test_all_introns_extracted(self):
        """Test case where all introns pass pre-filtering."""
        introns = []
        for i in range(5):
            coords = GenomicCoordinate(
                chromosome="chr1",
                start=(i + 1) * 1000,
                stop=(i + 1) * 1000 + 500,
                strand="+",
                system="1-based"
            )
            metadata = IntronMetadata(
                grandparent=f"gene{i}",
                parent=f"transcript{i}",
                parent_length=1000
            )
            introns.append(Intron(
                intron_id=f"intron{i}",
                coordinates=coords,
                metadata=metadata
            ))

        result = prefilter_introns(
            introns=introns,
            min_length=30,
            longest_only=False,
            include_duplicates=False
        )

        assert len(result.extract_list) == 5
        assert len(result.skip_list) == 0
        assert result.stats['extract'] == 5
        assert result.stats['skip'] == 0

    def test_all_introns_skipped(self):
        """Test case where all introns are too short."""
        introns = []
        for i in range(5):
            coords = GenomicCoordinate(
                chromosome="chr1",
                start=(i + 1) * 100,
                stop=(i + 1) * 100 + 10,  # Only 10bp long
                strand="+",
                system="1-based"
            )
            metadata = IntronMetadata(
                grandparent=f"gene{i}",
                parent=f"transcript{i}"
            )
            introns.append(Intron(
                intron_id=f"intron{i}",
                coordinates=coords,
                metadata=metadata
            ))

        result = prefilter_introns(
            introns=introns,
            min_length=30,
            longest_only=False,
            include_duplicates=False
        )

        assert len(result.extract_list) == 0
        assert len(result.skip_list) == 5
        assert result.stats['too_short'] == 5

    def test_prefilter_result_structure(self):
        """Test that PrefilterResult has correct structure."""
        introns = self.create_test_introns()
        result = prefilter_introns(introns)

        assert isinstance(result, PrefilterResult)
        assert isinstance(result.extract_list, list)
        assert isinstance(result.skip_list, list)
        assert isinstance(result.stats, dict)
        assert 'total' in result.stats
        assert 'extract' in result.stats
        assert 'skip' in result.stats
        assert 'too_short' in result.stats
        assert 'duplicate' in result.stats
        assert 'not_longest_isoform' in result.stats


class TestPrefilteringIntegration:
    """Integration tests for pre-filtering with realistic scenarios."""

    def test_human_genome_like_scenario(self):
        """
        Test pre-filtering with a scenario similar to human genome.

        Simulates:
        - 2.1M total introns
        - ~10% duplicates
        - ~85% from non-longest isoforms
        - ~5% too short

        Expected: ~5-15% extraction rate with longest_only=True
        """
        introns = []
        num_introns = 1000  # Scaled down from 2.1M for testing

        # 85% from non-longest isoforms
        for i in range(int(num_introns * 0.85)):
            coords = GenomicCoordinate(
                chromosome="chr1",
                start=(i + 1) * 1000,
                stop=(i + 1) * 1000 + 100,
                strand="+",
                system="1-based"
            )
            metadata = IntronMetadata(
                grandparent=f"gene{i % 100}",
                parent=f"transcript_short_{i}",
                parent_length=500  # Shorter
            )
            introns.append(Intron(
                intron_id=f"intron{i}",
                coordinates=coords,
                metadata=metadata
            ))

        # 10% longest isoforms
        for i in range(int(num_introns * 0.10)):
            coords = GenomicCoordinate(
                chromosome="chr1",
                start=(i + 1000) * 1000,
                stop=(i + 1000) * 1000 + 100,
                strand="+",
                system="1-based"
            )
            metadata = IntronMetadata(
                grandparent=f"gene{i % 100}",
                parent=f"transcript_long_{i}",
                parent_length=1000  # Longer
            )
            introns.append(Intron(
                intron_id=f"intron{i + 850}",
                coordinates=coords,
                metadata=metadata
            ))

        # 5% too short
        for i in range(int(num_introns * 0.05)):
            coords = GenomicCoordinate(
                chromosome="chr1",
                start=(i + 2000) * 1000,
                stop=(i + 2000) * 1000 + 10,  # Too short
                strand="+",
                system="1-based"
            )
            metadata = IntronMetadata(
                grandparent=f"gene_short_{i}",
                parent=f"transcript{i}"
            )
            introns.append(Intron(
                intron_id=f"intron{i + 950}",
                coordinates=coords,
                metadata=metadata
            ))

        result = prefilter_introns(
            introns=introns,
            min_length=30,
            longest_only=True,
            include_duplicates=False
        )

        # Should extract roughly 10-15% (longest isoforms only, minus too-short)
        extraction_rate = len(result.extract_list) / len(introns)
        assert 0.05 <= extraction_rate <= 0.20  # Reasonable range

        # Should skip majority
        skip_rate = len(result.skip_list) / len(introns)
        assert skip_rate >= 0.80

        print(f"\nHuman-like scenario results:")
        print(f"  Total introns: {result.stats['total']}")
        print(f"  Extract: {result.stats['extract']} ({extraction_rate:.1%})")
        print(f"  Skip: {result.stats['skip']} ({skip_rate:.1%})")
        print(f"  Memory savings: ~{skip_rate:.0%}")
