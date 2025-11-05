"""
Filter and deduplicate introns.

This module handles quality control filtering, duplicate detection,
longest isoform identification, and overlap detection for introns.
"""

from typing import List, Dict, Set, Tuple, Optional
from collections import defaultdict
from dataclasses import dataclass

from core.intron import Intron
from utils.sequences import is_valid_dna, has_ambiguous_bases


@dataclass
class FilterStats:
    """Statistics from filtering operations."""
    omitted_short: int = 0
    omitted_ambiguous: int = 0
    omitted_noncanonical: int = 0
    omitted_overlap: int = 0
    omitted_isoform: int = 0
    duplicates: int = 0
    total_introns: int = 0
    kept_introns: int = 0


class IntronFilter:
    """
    Filters and deduplicates introns based on quality criteria.

    This class implements the filtering logic from the original intronIC,
    including:
    - Length filtering
    - Sequence quality filtering (ambiguous bases)
    - Duplicate detection
    - Longest isoform identification
    - Coordinate overlap detection

    Examples:
        >>> filter_obj = IntronFilter(
        ...     min_length=30,
        ...     allow_noncanonical=False,
        ...     allow_overlap=False,
        ...     longest_only=False
        ... )
        >>> filtered = list(filter_obj.filter_introns(introns))
    """

    def __init__(
        self,
        min_length: int = 30,
        bp_matrix_length: int = 7,
        scoring_regions: List[str] = None,
        allow_noncanonical: bool = False,
        allow_overlap: bool = False,
        longest_only: bool = False,
        include_duplicates: bool = False
    ):
        """
        Initialize the intron filter.

        Args:
            min_length: Minimum intron length (default: 30)
            bp_matrix_length: Minimum valid bp region length (default: 7)
            scoring_regions: Regions to check for ambiguous bases (default: ['five', 'three'])
            allow_noncanonical: Include non-canonical introns (default: False)
            allow_overlap: Include overlapping introns (default: False)
            longest_only: Only keep longest isoform (default: False)
            include_duplicates: Include duplicate introns (default: False)
        """
        self.min_length = min_length
        self.bp_matrix_length = bp_matrix_length
        self.scoring_regions = scoring_regions or ['five', 'three']
        self.allow_noncanonical = allow_noncanonical
        self.allow_overlap = allow_overlap
        self.longest_only = longest_only
        self.include_duplicates = include_duplicates

        # Tracking structures
        self.intron_index: Dict[Tuple, Dict] = defaultdict(lambda: defaultdict(dict))
        self.longest_isoforms: Dict[str, str] = {}
        self.duplicate_map: Dict[int, Set[str]] = defaultdict(set)
        self.overlap_map: Dict[int, Set[str]] = defaultdict(set)
        self.stats = FilterStats()

    def filter_introns(self, introns: List[Intron]) -> List[Intron]:
        """
        Filter a list of introns based on quality criteria.

        Args:
            introns: List of Intron objects

        Returns:
            List of filtered Intron objects

        Examples:
            >>> filter_obj = IntronFilter(min_length=50)
            >>> filtered = filter_obj.filter_introns(intron_list)
            >>> print(f"Kept {len(filtered)}/{len(intron_list)} introns")
        """
        filtered_introns = []
        self.stats.total_introns = len(introns)

        # First pass: Identify longest transcripts per gene
        self._identify_longest_isoforms(introns)

        # Second pass: Filter introns
        for intron in introns:
            # Step 1: Check omission criteria
            self._check_omission(intron)

            # Step 2: Tag duplicates and longest isoforms
            self._tag_intron(intron)

            # Step 3: Re-check omission (may change based on tags)
            self._check_omission(intron)

            # Step 4: Update statistics
            self._update_stats(intron)

            # Step 5: Decide if intron should be kept
            if self._should_keep(intron):
                filtered_introns.append(intron)
                self.stats.kept_introns += 1

        return filtered_introns

    def _identify_longest_isoforms(self, introns: List[Intron]) -> None:
        """
        First pass: identify longest transcript per gene.

        Populates self.longest_isoforms dictionary with the longest transcript
        for each gene based on transcript length.

        Args:
            introns: List of all introns
        """
        for intron in introns:
            grandparent = intron.metadata.grandparent
            parent = intron.metadata.parent
            parent_length = intron.metadata.parent_length or 0

            if grandparent:
                if grandparent not in self.longest_isoforms:
                    # First transcript for this gene
                    self.longest_isoforms[grandparent] = {
                        'transcript': parent,
                        'length': parent_length
                    }
                else:
                    # Update if this transcript is longer
                    current = self.longest_isoforms[grandparent]
                    if parent_length > current['length']:
                        current['transcript'] = parent
                        current['length'] = parent_length

    def _check_omission(self, intron: Intron) -> None:
        """
        Check if intron meets omission criteria.

        Updates intron.metadata.omitted with reason code:
        - 's' = short
        - 'a' = ambiguous sequence
        - 'n' = noncanonical
        - 'v' = coordinate overlap
        - 'i' = not in longest isoform

        Args:
            intron: Intron object to check
        """
        # Check length
        if intron.length < self.min_length:
            intron.metadata.omitted = 's'
            return

        # Check for ambiguous bases in scoring regions
        if intron.sequences:
            for region in self.scoring_regions:
                if region == 'five' and intron.sequences.five_seq:
                    if has_ambiguous_bases(intron.sequences.five_seq):
                        intron.metadata.omitted = 'a'
                        return
                elif region == 'three' and intron.sequences.three_seq:
                    if has_ambiguous_bases(intron.sequences.three_seq):
                        intron.metadata.omitted = 'a'
                        return

            # Check bp region length and quality
            if intron.sequences.bp_region_seq:
                bp_seq = intron.sequences.bp_region_seq
                # Count longest stretch of valid bases
                valid_length = self._longest_valid_stretch(bp_seq)

                if valid_length < self.bp_matrix_length:
                    if len(bp_seq) < self.bp_matrix_length:
                        intron.metadata.omitted = 's'
                    else:
                        intron.metadata.omitted = 'a'
                    return

        # Check noncanonical
        if not self.allow_noncanonical and intron.metadata.noncanonical:
            intron.metadata.omitted = 'n'
            return

        # Check longest isoform
        if self.longest_only and not intron.metadata.longest_isoform:
            intron.metadata.omitted = 'i'
            return

        # Check overlap
        if not self.allow_overlap and intron.metadata.overlap:
            intron.metadata.omitted = 'v'
            return

    @staticmethod
    def _longest_valid_stretch(seq: str) -> int:
        """
        Find longest stretch of valid DNA bases (A, C, T, G).

        Args:
            seq: DNA sequence

        Returns:
            Length of longest valid stretch
        """
        import re
        matches = re.findall(r'[ATCG]+', seq.upper())
        if matches:
            return len(max(matches, key=len))
        return 0

    def _tag_intron(self, intron: Intron) -> None:
        """
        Tag intron as duplicate, overlapping, or from longest isoform.

        Args:
            intron: Intron object to tag
        """
        # Create region key (separate indices for omitted vs not omitted)
        region_id = (
            intron.coordinates.chromosome,
            intron.coordinates.strand,
            not bool(intron.metadata.omitted)
        )

        # Create coordinate key
        coord_key = (intron.coordinates.start, intron.coordinates.stop)

        region_idx = self.intron_index[region_id]

        # Check for duplicate
        if coord_key not in region_idx:
            # First occurrence of these coordinates
            intron.metadata.duplicate = None
            region_idx[coord_key] = {
                'parent_length': intron.metadata.parent_length or 0,
                'family_size': intron.metadata.family_size or 0,
                'unique_id': id(intron)  # Unique identifier for this intron
            }
        else:
            # Duplicate found - reference the original
            intron.metadata.duplicate = region_idx[coord_key]['unique_id']
            intron.metadata.overlap = intron.metadata.duplicate

        # Check for longest isoform
        # longest_isoforms dictionary already populated by _identify_longest_isoforms
        parent = intron.metadata.parent
        grandparent = intron.metadata.grandparent

        if grandparent and grandparent in self.longest_isoforms:
            # Check if this intron's transcript is the longest for its gene
            longest_transcript = self.longest_isoforms[grandparent]['transcript']
            intron.metadata.longest_isoform = (parent == longest_transcript)
        else:
            # No grandparent info, assume longest
            intron.metadata.longest_isoform = True

        # Check for coordinate overlap (only if not duplicate and not omitted)
        if not intron.metadata.duplicate and not intron.metadata.omitted:
            if (not intron.metadata.longest_isoform and
                not self.allow_overlap and
                not self.longest_only):
                # Check if coordinates overlap with any existing intron
                seen_coords = list(region_idx.keys())
                overlap = self._check_coord_overlap(coord_key, seen_coords)
                if overlap:
                    intron.metadata.overlap = region_idx[overlap]['unique_id']
                else:
                    intron.metadata.overlap = None

    @staticmethod
    def _check_coord_overlap(
        coord: Tuple[int, int],
        seen_coords: List[Tuple[int, int]]
    ) -> Optional[Tuple[int, int]]:
        """
        Check if coordinates overlap with any seen coordinates.

        Args:
            coord: (start, stop) tuple
            seen_coords: List of (start, stop) tuples

        Returns:
            Overlapping coordinates or None
        """
        for seen in seen_coords:
            # Skip self
            if coord == seen:
                continue

            # Check overlap using the elegant formula
            # (a.start - b.stop) * (a.stop - b.start) < 0 indicates overlap
            val = (coord[0] - seen[1]) * (coord[1] - seen[0])
            if val < 0:
                return seen

        return None

    def _update_stats(self, intron: Intron) -> None:
        """
        Update filtering statistics.

        Args:
            intron: Intron object
        """
        if intron.metadata.duplicate:
            self.stats.duplicates += 1

        omitted = intron.metadata.omitted
        if omitted == 's':
            self.stats.omitted_short += 1
        elif omitted == 'a':
            self.stats.omitted_ambiguous += 1
        elif omitted == 'n':
            self.stats.omitted_noncanonical += 1
        elif omitted == 'v':
            self.stats.omitted_overlap += 1
        elif omitted == 'i':
            self.stats.omitted_isoform += 1

    def _should_keep(self, intron: Intron) -> bool:
        """
        Determine if intron should be kept in filtered set.

        Args:
            intron: Intron object

        Returns:
            True if intron should be kept, False otherwise
        """
        # Omitted introns are not kept
        if intron.metadata.omitted:
            return False

        # Duplicates only kept if explicitly allowed
        if intron.metadata.duplicate and not self.include_duplicates:
            return False

        return True

    def get_duplicate_map(self) -> Dict[str, Set[str]]:
        """
        Get mapping of representative introns to their duplicates.

        Returns:
            Dictionary mapping intron names to sets of duplicate names
        """
        return dict(self.duplicate_map)

    def get_overlap_map(self) -> Dict[str, Set[str]]:
        """
        Get mapping of introns to their overlapping introns.

        Returns:
            Dictionary mapping intron names to sets of overlapping names
        """
        return dict(self.overlap_map)

    def get_stats(self) -> FilterStats:
        """
        Get filtering statistics.

        Returns:
            FilterStats object with counts
        """
        return self.stats


def filter_introns(
    introns: List[Intron],
    min_length: int = 30,
    allow_noncanonical: bool = False,
    allow_overlap: bool = False,
    longest_only: bool = False
) -> List[Intron]:
    """
    Convenience function to filter introns.

    This is a functional wrapper for backwards compatibility with
    the original intronIC API.

    Args:
        introns: List of Intron objects
        min_length: Minimum intron length
        allow_noncanonical: Include non-canonical introns
        allow_overlap: Include overlapping introns
        longest_only: Only keep longest isoform

    Returns:
        List of filtered Intron objects

    Examples:
        >>> filtered = filter_introns(intron_list, min_length=50)
    """
    filter_obj = IntronFilter(
        min_length=min_length,
        allow_noncanonical=allow_noncanonical,
        allow_overlap=allow_overlap,
        longest_only=longest_only
    )
    return filter_obj.filter_introns(introns)
