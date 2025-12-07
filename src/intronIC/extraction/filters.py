"""
Filter and deduplicate introns.

This module handles quality control filtering, duplicate detection,
longest isoform identification, and overlap detection for introns.
"""

from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

from intronIC.core.intron import Intron, OmissionReason
from intronIC.utils.sequences import has_ambiguous_bases, is_valid_dna


@dataclass
class FilterStats:
    """Statistics from filtering operations.

    Tracks counts for each exclusion category. Categories are mutually exclusive
    with duplicates taking priority (most specific reason for exclusion).
    """

    # Counts by category (mutually exclusive)
    duplicates: int = 0
    short: int = 0
    ambiguous: int = 0
    noncanonical: int = 0
    overlap: int = 0
    isoform: int = 0

    # Totals
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
        include_duplicates: bool = False,
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
        self.scoring_regions = scoring_regions or ["five", "three"]
        self.allow_noncanonical = allow_noncanonical
        self.allow_overlap = allow_overlap
        self.longest_only = longest_only
        self.include_duplicates = include_duplicates

        # Tracking structures
        self.intron_index: Dict[Tuple, Dict] = defaultdict(lambda: defaultdict(dict))
        self.longest_isoforms: Dict[str, str] = {}
        self.duplicate_map: Dict[str, Set[str]] = defaultdict(set)
        self.overlap_map: Dict[str, Set[str]] = defaultdict(set)
        self.stats = FilterStats()

    def filter_introns(self, introns: List[Intron]) -> List[Intron]:
        """
        Filter a list of introns based on quality criteria.

        All introns are processed and returned with their metadata updated:
        - metadata.omitted set to appropriate OmissionReason for filtered introns
        - metadata.duplicate set for coordinate duplicates
        - metadata.longest_isoform set for isoform filtering

        This allows callers to include omitted introns in output files (matching
        original intronIC behavior where .meta.iic includes both scored and omitted).

        Args:
            introns: List of Intron objects

        Returns:
            List of ALL Intron objects (with metadata updated)

        Examples:
            >>> filter_obj = IntronFilter(min_length=50)
            >>> all_introns = filter_obj.filter_introns(intron_list)
            >>> kept = [i for i in all_introns if filter_obj._should_keep(i)]
            >>> print(f"Kept {len(kept)}/{len(all_introns)} introns")
        """
        self.stats.total_introns = len(introns)

        # Sort introns using hierarchical sort (matching original's get_sub_seqs sorting)
        # This is CRITICAL: longest isoform identification depends on processing order!
        sorted_introns = self._sort_introns(introns)

        # First pass: Identify longest transcripts per gene (using sorted order)
        self._identify_longest_isoforms(sorted_introns)

        # Second pass: Filter introns (using sorted order)
        for intron in sorted_introns:
            # Step 1: Check omission criteria
            self._check_omission(intron)

            # Step 2: Tag duplicates and longest isoforms
            self._tag_intron(intron)

            # Step 3: Re-check omission (may change based on tags)
            self._check_omission(intron)

            # Step 4: Update statistics
            self._update_stats(intron)

            # Step 5: Track kept count
            if self._should_keep(intron):
                self.stats.kept_introns += 1

        # Return ALL introns (callers can filter by _should_keep if needed)
        return sorted_introns

    @staticmethod
    def _sort_introns(introns: List[Intron]) -> List[Intron]:
        """
        Sort introns using hierarchical_sort_attrs logic from original.

        Matches original intronIC's sorting in get_sub_seqs (line 2674):
        - defined_by (CDS before exon)
        - parent_length descending (KEY for longest isoform!)
        - parent (transcript ID)
        - family_size descending
        - index (intron position)
        - line_number (final tiebreaker for reproducibility)

        Args:
            introns: List of introns to sort

        Returns:
            Sorted list of introns
        """
        return sorted(
            introns,
            key=lambda i: (
                i.metadata.defined_by or "",  # CDS before exon
                -(i.metadata.parent_length or 0),  # Descending by parent length
                i.metadata.parent or "",  # Transcript ID
                -(i.metadata.family_size or 0),  # Descending by family size
                i.metadata.index or 0,  # Intron index
                i.metadata.line_number or 0,  # Final tiebreaker
            ),
        )

    def _identify_longest_isoforms(self, introns: List[Intron]) -> None:
        """
        First pass: identify "longest" transcript per gene.

        Expects introns to already be sorted by _sort_introns() (hierarchical order).
        Since introns are sorted by parent_length descending, the first transcript
        seen for each gene is the longest.

        Args:
            introns: List of introns (must already be sorted!)
        """
        for intron in introns:
            grandparent = intron.metadata.grandparent
            parent = intron.metadata.parent

            if grandparent:
                if grandparent not in self.longest_isoforms:
                    # First transcript for this gene (after hierarchical sorting)
                    # Since sorted by parent_length descending, this IS the longest
                    self.longest_isoforms[grandparent] = parent

    def _check_omission(self, intron: Intron) -> None:
        """
        Check if intron meets omission criteria.

        Updates intron.metadata.omitted with OmissionReason enum:
        - OmissionReason.SHORT = short
        - OmissionReason.AMBIGUOUS = ambiguous sequence
        - OmissionReason.NONCANONICAL = noncanonical
        - OmissionReason.OVERLAP = coordinate overlap
        - OmissionReason.ISOFORM = not in longest isoform

        Args:
            intron: Intron object to check
        """
        # Clear omitted flag to re-evaluate (important for second pass after tagging)
        intron.metadata.omitted = OmissionReason.NONE

        # Check length
        if intron.length < self.min_length:
            intron.metadata.omitted = OmissionReason.SHORT
            return

        # Check for ambiguous bases in scoring regions
        if intron.sequences:
            for region in self.scoring_regions:
                if region == "five" and intron.sequences.five_seq:
                    if has_ambiguous_bases(intron.sequences.five_seq):
                        intron.metadata.omitted = OmissionReason.AMBIGUOUS
                        return
                elif region == "three" and intron.sequences.three_seq:
                    if has_ambiguous_bases(intron.sequences.three_seq):
                        intron.metadata.omitted = OmissionReason.AMBIGUOUS
                        return

            # Check bp region length and quality
            if intron.sequences.bp_region_seq:
                bp_seq = intron.sequences.bp_region_seq
                # Count longest stretch of valid bases
                valid_length = self._longest_valid_stretch(bp_seq)

                if valid_length < self.bp_matrix_length:
                    if len(bp_seq) < self.bp_matrix_length:
                        intron.metadata.omitted = OmissionReason.SHORT
                    else:
                        intron.metadata.omitted = OmissionReason.AMBIGUOUS
                    return

        # Check noncanonical
        if not self.allow_noncanonical and intron.metadata.noncanonical:
            intron.metadata.omitted = OmissionReason.NONCANONICAL
            return

        # Check longest isoform
        # IMPORTANT: Use 'is False' not 'not' to match original behavior
        # Original only omits when explicitly False, not when None
        # (see intronIC.py line 718: "if longest_only and self.longest_isoform is False")
        if self.longest_only and intron.metadata.longest_isoform is False:
            intron.metadata.omitted = OmissionReason.ISOFORM
            return

        # Check overlap
        if not self.allow_overlap and intron.metadata.overlap:
            intron.metadata.omitted = OmissionReason.OVERLAP
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

        matches = re.findall(r"[ATCG]+", seq.upper())
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
            intron.metadata.omitted == OmissionReason.NONE,
        )

        # Create coordinate key
        coord_key = (intron.coordinates.start, intron.coordinates.stop)

        region_idx = self.intron_index[region_id]

        # Check for duplicate (skip for sequence-only introns without real coordinates)
        from intronIC.core.intron import IntronFlags

        is_sequence_only = IntronFlags.SEQUENCE_ONLY in intron.metadata.flags

        if is_sequence_only:
            # Skip duplicate detection for sequence-only introns (no real coordinates)
            intron.metadata.duplicate = None
        elif coord_key not in region_idx:
            # First occurrence of these coordinates
            intron.metadata.duplicate = None
            region_idx[coord_key] = {
                "parent_length": intron.metadata.parent_length or 0,
                "family_size": intron.metadata.family_size or 0,
                "unique_id": id(intron),  # Unique identifier for this intron
                "intron_id": intron.intron_id,  # Store ID for duplicate mapping
            }
        else:
            # Duplicate found - reference the original
            intron.metadata.duplicate = region_idx[coord_key]["unique_id"]
            intron.metadata.overlap = intron.metadata.duplicate
            # Add dynamic tag for duplicates
            intron.metadata.dynamic_tags.add("[d]")
            # Record the duplicate mapping (representative -> duplicate)
            representative_id = region_idx[coord_key]["intron_id"]
            self.duplicate_map[representative_id].add(intron.intron_id)

        # Check for longest isoform
        # longest_isoforms dictionary already populated by _identify_longest_isoforms
        parent = intron.metadata.parent
        grandparent = intron.metadata.grandparent

        if grandparent and grandparent in self.longest_isoforms:
            # Check if this intron's transcript is the "longest" (first seen) for its gene
            longest_transcript = self.longest_isoforms[grandparent]
            intron.metadata.longest_isoform = parent == longest_transcript
            # Add dynamic tag for alternative isoforms
            if not intron.metadata.longest_isoform:
                intron.metadata.dynamic_tags.add("[i]")
        else:
            # No grandparent info, assume longest
            intron.metadata.longest_isoform = True

        # Check for coordinate overlap (only if not duplicate and not omitted)
        if (
            not intron.metadata.duplicate
            and intron.metadata.omitted == OmissionReason.NONE
        ):
            if (
                not intron.metadata.longest_isoform
                and not self.allow_overlap
                and not self.longest_only
            ):
                # Check if coordinates overlap with any existing intron
                seen_coords = list(region_idx.keys())
                overlap = self._check_coord_overlap(coord_key, seen_coords)
                if overlap:
                    intron.metadata.overlap = region_idx[overlap]["unique_id"]
                    # Record the overlap mapping (representative -> overlapping)
                    representative_id = region_idx[overlap]["intron_id"]
                    self.overlap_map[representative_id].add(intron.intron_id)
                else:
                    intron.metadata.overlap = None

    @staticmethod
    def _check_coord_overlap(
        coord: Tuple[int, int], seen_coords: List[Tuple[int, int]]
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

        Tracks how many introns fall into each category based on their properties
        (not just omission status). This allows reporting "found X, excluded Y"
        even when user options include a category.

        Categories are mutually exclusive with duplicates taking priority
        (most specific reason for exclusion).

        Args:
            intron: Intron object
        """
        # Count duplicates first (highest priority - most specific reason)
        if intron.metadata.duplicate:
            self.stats.duplicates += 1
            return

        # For non-duplicates, count by property/category
        # Check omission reason first (always accurate for short/ambiguous)
        omitted = intron.metadata.omitted
        if omitted == OmissionReason.SHORT:
            self.stats.short += 1
        elif omitted == OmissionReason.AMBIGUOUS:
            self.stats.ambiguous += 1
        # For optional categories, check the property itself (not just omission)
        # This way we count them even when user includes them
        elif intron.metadata.noncanonical:
            self.stats.noncanonical += 1
        elif intron.metadata.overlap:
            self.stats.overlap += 1
        elif intron.metadata.longest_isoform is False:
            self.stats.isoform += 1

    def _should_keep(self, intron: Intron) -> bool:
        """
        Determine if intron should be kept in filtered set.

        Args:
            intron: Intron object

        Returns:
            True if intron should be kept, False otherwise
        """
        # Omitted introns are not kept
        if intron.metadata.omitted != OmissionReason.NONE:
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


def should_extract_sequences_for(
    intron: Intron,
    min_length: int,
    longest_only: bool,
    longest_isoforms: Dict[str, str],
    seen_coordinates: Set[Tuple[int, int]],
    include_duplicates: bool,
) -> bool:
    """
    Determine if sequences should be extracted for this intron based on metadata only.

    This function allows pre-filtering before expensive sequence extraction.
    It checks criteria that can be determined without sequences:
    - Length (too short to ever be useful)
    - Longest isoform (if longest_only=True and not longest isoform)
    - Duplicates (if not include_duplicates and is duplicate)

    Args:
        intron: Intron object (without sequences)
        min_length: Minimum intron length
        longest_only: Only extract for longest isoform
        longest_isoforms: Dict mapping gene ID to longest transcript ID
        seen_coordinates: Set of (start, stop) tuples already seen
        include_duplicates: Whether to extract for all duplicates

    Returns:
        True if sequences should be extracted, False otherwise
    """
    # Check length - too short, skip extraction
    if intron.length < min_length:
        return False

    # Check longest isoform - if longest_only=True, only extract for longest
    if longest_only:
        grandparent = intron.metadata.grandparent
        parent = intron.metadata.parent
        if grandparent and grandparent in longest_isoforms:
            is_longest = parent == longest_isoforms[grandparent]
            if not is_longest:
                return False

    # Check duplicates - only extract for first occurrence
    # We reuse sequences for all duplicates regardless of include_duplicates flag
    # The include_duplicates flag affects filtering/output, not extraction
    coord_key = (intron.coordinates.start, intron.coordinates.stop)
    if coord_key in seen_coordinates:
        # This is a duplicate - skip extraction (will reuse from first occurrence)
        return False

    # Extract by default
    return True


@dataclass
class PrefilterResult:
    """Result of pre-filtering introns before sequence extraction."""

    extract_list: List[Intron]  # Introns that need sequences extracted
    skip_list: List[Intron]  # Introns that can skip extraction
    stats: Dict[str, int]  # Statistics about filtering decisions


def prefilter_introns(
    introns: List[Intron],
    min_length: int = 30,
    longest_only: bool = False,
    include_duplicates: bool = False,
) -> PrefilterResult:
    """
    Pre-filter introns before sequence extraction based on metadata only.

    This function separates introns into two groups:
    1. extract_list: Introns that need sequences extracted
    2. skip_list: Introns that can skip extraction (too short, wrong isoform, duplicates)

    By filtering before extraction, we avoid the expensive operation of extracting
    sequences for ~85-90% of introns that will be omitted anyway.

    Args:
        introns: List of Intron objects (without sequences)
        min_length: Minimum intron length
        longest_only: Only extract for longest isoform per gene
        include_duplicates: Extract for all duplicate coordinates

    Returns:
        PrefilterResult containing extract_list, skip_list, and statistics

    Examples:
        >>> result = prefilter_introns(introns, min_length=50, longest_only=True)
        >>> print(f"Extracting {len(result.extract_list)}/{len(introns)} introns")
        >>> # Extract sequences only for result.extract_list
        >>> introns_with_seq = extract_sequences(result.extract_list, ...)
    """
    extract_list = []
    skip_list = []
    stats = {
        "total": len(introns),
        "too_short": 0,
        "not_longest_isoform": 0,
        "duplicate": 0,
        "extract": 0,
        "skip": 0,
    }

    # Sort introns using same hierarchical sort as filtering
    # This ensures longest isoform identification is correct
    sorted_introns = IntronFilter._sort_introns(introns)

    # Identify longest isoforms (same logic as IntronFilter)
    longest_isoforms: Dict[str, str] = {}
    for intron in sorted_introns:
        grandparent = intron.metadata.grandparent
        parent = intron.metadata.parent
        if grandparent and grandparent not in longest_isoforms:
            longest_isoforms[grandparent] = parent

    # Track seen coordinates per (chr, strand) to identify duplicates
    seen_by_region: Dict[Tuple, Set[Tuple[int, int]]] = defaultdict(set)

    # Process each intron
    for intron in sorted_introns:
        region_key = (intron.coordinates.chromosome, intron.coordinates.strand)
        coord_key = (intron.coordinates.start, intron.coordinates.stop)
        seen_coords = seen_by_region[region_key]

        # Decide if we should extract sequences
        should_extract = should_extract_sequences_for(
            intron=intron,
            min_length=min_length,
            longest_only=longest_only,
            longest_isoforms=longest_isoforms,
            seen_coordinates=seen_coords,
            include_duplicates=include_duplicates,
        )

        if should_extract:
            extract_list.append(intron)
            stats["extract"] += 1
            # Mark coordinates as seen AFTER deciding
            seen_by_region[region_key].add(coord_key)
        else:
            skip_list.append(intron)
            stats["skip"] += 1

            # Track why it was skipped (for statistics)
            if intron.length < min_length:
                stats["too_short"] += 1
            elif longest_only:
                grandparent = intron.metadata.grandparent
                parent = intron.metadata.parent
                if grandparent and grandparent in longest_isoforms:
                    if parent != longest_isoforms[grandparent]:
                        stats["not_longest_isoform"] += 1
            elif coord_key in seen_coords:
                stats["duplicate"] += 1

    return PrefilterResult(extract_list=extract_list, skip_list=skip_list, stats=stats)


def filter_introns(
    introns: List[Intron],
    min_length: int = 30,
    allow_noncanonical: bool = False,
    allow_overlap: bool = False,
    longest_only: bool = False,
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
        longest_only=longest_only,
    )
    return filter_obj.filter_introns(introns)
