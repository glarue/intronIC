"""
Generate introns from pairs of exons.

This module handles the core logic of creating Intron objects from consecutive
Exon objects in a transcript.
"""

from itertools import islice
from typing import TYPE_CHECKING, Dict, Iterator, List, Optional, Tuple, Union

import numpy as np

from intronIC.core.intron import Intron
from intronIC.core.models import Exon, Gene, Transcript

if TYPE_CHECKING:
    from intronIC.cli.messenger import Messenger


class IntronGenerator:
    """
    Generates Intron objects from exon pairs within transcripts.

    This class implements the "intronator" logic that creates introns
    from consecutive exons, handling strand direction and checking for
    overlapping exons.

    Examples:
        >>> generator = IntronGenerator()
        >>> exons = [exon1, exon2, exon3]  # From a transcript
        >>> introns = list(generator.generate_from_exons(exons))
        >>> print(f"Created {len(introns)} introns")
    """

    def __init__(self, debug: bool = False, messenger: Optional["Messenger"] = None):
        """
        Initialize the IntronGenerator.

        Args:
            debug: If True, print detailed messages about touching exons
            messenger: Optional Messenger instance for logging (if None, uses print)
        """
        self.debug = debug
        self.messenger = messenger
        self.touching_exons_count = 0  # Track total touching exons found

    @staticmethod
    def _sliding_window(
        seq: List[Exon], window_size: int = 2
    ) -> Iterator[Tuple[Exon, ...]]:
        """
        Generate sliding window over sequence.

        Args:
            seq: Sequence to window over
            window_size: Size of window (default: 2 for exon pairs)

        Yields:
            Tuples of consecutive items
        """
        it = iter(seq)
        result = tuple(islice(it, window_size))
        if len(result) == window_size:
            yield result
        for elem in it:
            result = result[1:] + (elem,)
            yield result

    @staticmethod
    def _check_overlap(exon1: Exon, exon2: Exon) -> bool:
        """
        Check if two exons overlap.

        Uses the elegant formula: (a.start - b.stop) * (a.stop - b.start) < 0
        indicates overlap.

        Args:
            exon1: First exon
            exon2: Second exon

        Returns:
            True if exons overlap, False otherwise

        Examples:
            >>> e1 = Exon(feature_id='e1', ...)
            >>> e2 = Exon(feature_id='e2', ...)
            >>> IntronGenerator._check_overlap(e1, e2)
            True
        """
        c1 = exon1.coordinates
        c2 = exon2.coordinates
        val = (c1.start - c2.stop) * (c1.stop - c2.start)
        return val < 0

    @staticmethod
    def _sort_exons_by_coding_direction(exons: List[Exon]) -> List[Exon]:
        """
        Sort exons in coding direction (5' to 3').

        For positive strand: ascending by start position
        For negative strand: descending by start position

        Args:
            exons: List of Exon objects

        Returns:
            Sorted list of exons
        """
        if not exons:
            return []

        # All exons in a transcript should have the same strand
        strand = exons[0].coordinates.strand

        if strand == "-":
            # Negative strand: sort descending (highest coordinate first)
            return sorted(exons, key=lambda e: e.coordinates.start, reverse=True)
        else:
            # Positive strand: sort ascending (lowest coordinate first)
            return sorted(exons, key=lambda e: e.coordinates.start)

    def generate_from_exons(self, exons: List[Exon]) -> Iterator[Intron]:
        """
        Generate introns from a list of exons.

        Args:
            exons: List of Exon objects from a single transcript

        Yields:
            Intron objects created from consecutive exon pairs

        Examples:
            >>> generator = IntronGenerator()
            >>> exons = [exon1, exon2, exon3]
            >>> for intron in generator.generate_from_exons(exons):
            ...     print(f"Intron: {intron.start}-{intron.stop}")
        """
        if len(exons) < 2:
            # Need at least 2 exons to make an intron
            return

        # Sort in coding direction
        sorted_exons = self._sort_exons_by_coding_direction(exons)

        # Generate introns from consecutive pairs
        for index, (upstream_exon, downstream_exon) in enumerate(
            self._sliding_window(sorted_exons), start=1
        ):
            # Check for overlapping exons (annotation error)
            if self._check_overlap(upstream_exon, downstream_exon):
                print(
                    f"[!] Warning: Overlapping exons in {upstream_exon.parent_id}: "
                    f"({upstream_exon.coordinates.start}, {upstream_exon.coordinates.stop}) and "
                    f"({downstream_exon.coordinates.start}, {downstream_exon.coordinates.stop}) - skipping"
                )
                continue

            # Create intron from exon pair
            try:
                intron = Intron.from_exon_pair(upstream_exon, downstream_exon)
            except ValueError as e:
                # Skip invalid introns (too short, overlapping, etc.)
                if "overlap or touch" in str(e):
                    continue  # Already warned about overlap
                else:
                    print(f"[!] Warning: Could not create intron: {e}")
                    continue

            # Store exon IDs for fractional position calculation
            intron.metadata.upstream_exon_id = upstream_exon.feature_id
            intron.metadata.downstream_exon_id = downstream_exon.feature_id

            # Set intron index (1-based position in this feature type's exon pairs)
            # Note: This will be recalculated in generate_from_transcript() to account
            # for touching exons and maintain ordinal numbering
            intron.metadata.index = index

            # Set line_number as average of both exons (matching original intronIC.py:573)
            # This enables tie-breaking in hierarchical sort for duplicate parent attributes
            upstream_line = upstream_exon.attributes.get("_line_number", 0)
            downstream_line = downstream_exon.attributes.get("_line_number", 0)
            intron.metadata.line_number = (upstream_line + downstream_line) / 2

            yield intron

    def generate_from_transcript(
        self,
        transcript: Transcript,
        feature_index: Dict[str, Union[Gene, Transcript, Exon]],
    ) -> Iterator[Intron]:
        """
        Generate introns from all exons in a transcript.

        Implements two-pass algorithm matching original intronIC:
        1. Generate introns from CDS features (is_coding=True)
        2. Generate introns from exon features (is_coding=False)
        3. Only add exon-derived introns that don't overlap CDS-derived ones

        This ensures CDS-defined intron boundaries take priority, and
        exon-only introns (e.g., in UTR regions) fill in the gaps.

        Design Decisions:
            **Touching/Zero-length Exon Pairs**: When adjacent exons have no gap
            between them (annotation errors where exon N ends at position X and
            exon N+1 starts at position X+1), these do NOT produce introns and
            are NOT counted in family_size. This matches the original intronIC
            v1.5.1 behavior. The intron index sequence remains contiguous (1, 2, 3...)
            with no gaps for these annotation artifacts.

            **Mixed CDS/Exon Transcripts**: For transcripts with both CDS and exon
            features, introns are generated from CDS first (prioritized for phase
            information), then exon-only introns (typically in UTR regions) are
            added if they don't overlap existing CDS introns. All introns are then
            sorted by genomic position and assigned sequential indices, ensuring
            proper ordering (e.g., 5' UTR introns come before CDS introns).

            **Family Size**: Represents the actual number of introns output for
            the transcript, not a theoretical count including annotation errors.

        Args:
            transcript: Transcript object with children (IDs)
            feature_index: Dictionary mapping feature IDs to objects

        Yields:
            Intron objects

        Examples:
            >>> generator = IntronGenerator()
            >>> for intron in generator.generate_from_transcript(transcript, feat_index):
            ...     print(intron.metadata.parent)
        """
        # Resolve exon IDs to Exon objects and separate by feature type
        cds_features = []
        exon_features = []

        for child_id in transcript.children:
            child = feature_index.get(child_id)
            if child and isinstance(child, Exon):
                # Verify exon belongs to this transcript
                if child.parent_id == transcript.feature_id:
                    if child.is_coding:
                        cds_features.append(child)
                    else:
                        exon_features.append(child)
                else:
                    print(
                        f"[!] Warning: Exon {child.feature_id} claims parent {child.parent_id} but is child of {transcript.feature_id}"
                    )

        # Calculate coding length (sum of CDS/exon feature lengths, not genomic span)
        # Matches original intronIC.py logic: prefer CDS length, fall back to exon length
        # This is used for longest isoform determination (parent_length sort key)
        if cds_features:
            coding_length = sum(cds.length for cds in cds_features)
        elif exon_features:
            coding_length = sum(exon.length for exon in exon_features)
        else:
            coding_length = 0  # No exons/CDS (shouldn't happen for valid transcripts)

        # Two-pass algorithm: CDS first, then exon (with overlap checking)
        # This ensures CDS-defined intron boundaries take priority, and
        # exon-only introns (e.g., in UTR regions) fill in the gaps.
        non_redundant_introns = []

        # Pass 1: Generate introns from CDS features
        if cds_features:
            for intron in self.generate_from_exons(cds_features):
                intron.metadata.defined_by = "cds"
                non_redundant_introns.append(intron)

        # Pass 2: Generate introns from exon features, excluding overlaps
        if exon_features:
            # Get coordinates of existing (CDS) introns for overlap checking
            existing_coords = [
                (i.coordinates.start, i.coordinates.stop) for i in non_redundant_introns
            ]

            for intron in self.generate_from_exons(exon_features):
                # Check if this exon-derived intron overlaps any CDS-derived intron
                intron_coords = (intron.coordinates.start, intron.coordinates.stop)
                if not self._check_intron_overlap(intron_coords, existing_coords):
                    # This is an exon-only intron (e.g., in UTR region)
                    intron.metadata.defined_by = "exon"
                    non_redundant_introns.append(intron)

        if not non_redundant_introns:
            return

        # Sort introns by genomic position (coding direction)
        strand = (
            non_redundant_introns[0].coordinates.strand
            if non_redundant_introns
            else "+"
        )
        if strand == "-":
            # Negative strand: sort descending (highest coord = first intron)
            non_redundant_introns.sort(key=lambda i: i.coordinates.start, reverse=True)
        else:
            # Positive strand: sort ascending (lowest coord = first intron)
            non_redundant_introns.sort(key=lambda i: i.coordinates.start)

        # Family size = number of actual introns output
        # Design decision: touching/zero-length exon pairs (annotation errors) are
        # silently skipped and NOT included in family_size. This matches the original
        # intronIC behavior (v1.5.1 line 421: family_size = len(non_redundant)).
        family_size = len(non_redundant_introns)

        # Assign sequential ordinal indices after sorting
        # This ensures proper ordering even when mixing CDS and exon-only introns
        # (e.g., 5' UTR introns come before CDS introns in coding direction)
        for index, intron in enumerate(non_redundant_introns, start=1):
            intron.metadata.index = index

        # Calculate fractional positions based on cumulative exon lengths
        # (matching original intronIC.py:429-434)
        exon_lengths = []
        for intron in non_redundant_introns:
            # Get upstream exon length
            upstream_exon_id = intron.metadata.upstream_exon_id
            if upstream_exon_id and upstream_exon_id in feature_index:
                upstream_exon = feature_index[upstream_exon_id]
                exon_lengths.append(upstream_exon.length)
            else:
                # Fallback if exon ID not found (shouldn't happen)
                exon_lengths.append(0)

        # Add last exon (downstream of last intron)
        last_intron = non_redundant_introns[-1]
        if (
            last_intron.metadata.downstream_exon_id
            and last_intron.metadata.downstream_exon_id in feature_index
        ):
            last_exon = feature_index[last_intron.metadata.downstream_exon_id]
            exon_lengths.append(last_exon.length)
        else:
            exon_lengths.append(0)

        # Calculate cumulative sum and fractional positions
        exon_cumsum = np.array(exon_lengths)[:-1].cumsum()
        aggregate_length = sum(exon_lengths)

        # Original multiplied by 100, but we store as 0.0-1.0 for clarity
        if aggregate_length > 0:
            frac_positions = (exon_cumsum / aggregate_length).round(4)
        else:
            frac_positions = np.zeros(len(non_redundant_introns))

        for array_index, intron in enumerate(non_redundant_introns):
            # Set metadata (index already assigned above via sequential enumeration)
            intron.metadata.family_size = family_size
            intron.metadata.parent = transcript.feature_id
            intron.metadata.parent_length = (
                coding_length  # Sum of CDS/exon lengths (not genomic span)
            )
            # Use array_index for fractional position array lookup (0-based)
            intron.metadata.fractional_position = float(frac_positions[array_index])

            # Set grandparent if available
            if transcript.parent_id and transcript.parent_id in feature_index:
                intron.metadata.grandparent = transcript.parent_id

            yield intron

    @staticmethod
    def _check_intron_overlap(
        intron_coords: Tuple[int, int], existing_coords: List[Tuple[int, int]]
    ) -> bool:
        """
        Check if intron coordinates overlap with or match any existing intron coordinates.

        Args:
            intron_coords: (start, stop) tuple for intron to check
            existing_coords: List of (start, stop) tuples for existing introns

        Returns:
            True if overlap or exact match found, False otherwise
        """
        for existing in existing_coords:
            # Using the elegant overlap formula: (a.start - b.stop) * (a.stop - b.start) < 0
            # Note: For exact matches (same coordinates), val = 0, so use <= 0 to detect them
            val = (intron_coords[0] - existing[1]) * (intron_coords[1] - existing[0])
            if val <= 0:
                return True
        return False

    def generate_from_gene(
        self, gene: Gene, feature_index: Dict[str, Union[Gene, Transcript, Exon]]
    ) -> Iterator[Intron]:
        """
        Generate introns from all transcripts in a gene.

        Args:
            gene: Gene object with transcript children (IDs)
            feature_index: Dictionary mapping feature IDs to objects

        Yields:
            Intron objects from all transcripts

        Examples:
            >>> generator = IntronGenerator()
            >>> introns = list(generator.generate_from_gene(gene, feat_index))
        """
        # Resolve transcript IDs to Transcript objects
        transcripts = []
        for child_id in gene.children:
            child = feature_index.get(child_id)
            if child and isinstance(child, Transcript):
                transcripts.append(child)

        for transcript in transcripts:
            yield from self.generate_from_transcript(transcript, feature_index)

    def generate_from_genes(
        self, genes: List[Gene], feature_index: Dict[str, Union[Gene, Transcript, Exon]]
    ) -> Iterator[Intron]:
        """
        Generate introns from a list of genes.

        Args:
            genes: List of Gene objects
            feature_index: Dictionary mapping feature IDs to objects

        Yields:
            All introns from all genes

        Examples:
            >>> generator = IntronGenerator()
            >>> introns = list(generator.generate_from_genes(gene_list, feat_index))
            >>> print(f"Total: {len(introns)} introns")
        """
        for gene in genes:
            yield from self.generate_from_gene(gene, feature_index)


def generate_introns_from_exons(exons: List[Exon]) -> Iterator[Intron]:
    """
    Convenience function to generate introns from exons.

    This is a functional wrapper for backwards compatibility with
    the original intronIC API.

    Args:
        exons: List of Exon objects from a transcript

    Yields:
        Intron objects

    Examples:
        >>> introns = list(generate_introns_from_exons(exon_list))
    """
    generator = IntronGenerator()
    return generator.generate_from_exons(exons)
