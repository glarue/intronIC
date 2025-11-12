"""
Generate introns from pairs of exons.

This module handles the core logic of creating Intron objects from consecutive
Exon objects in a transcript.
"""

from itertools import islice
from typing import List, Iterator, Tuple, Dict, Union

from core.models import Exon, Transcript, Gene
from core.intron import Intron


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

    @staticmethod
    def _sliding_window(seq: List[Exon], window_size: int = 2) -> Iterator[Tuple[Exon, ...]]:
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

        if strand == '-':
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

            # TODO: Store exon IDs for potential annotation updates
            # (would need to add upstream_exon/downstream_exon to IntronMetadata)
            # intron.metadata.upstream_exon = upstream_exon.feature_id
            # intron.metadata.downstream_exon = downstream_exon.feature_id

            # Set intron index (1-based position in transcript, matching original)
            intron.metadata.index = index

            # Set line_number as average of both exons (matching original intronIC.py:573)
            # This enables tie-breaking in hierarchical sort for duplicate parent attributes
            upstream_line = upstream_exon.attributes.get('_line_number', 0)
            downstream_line = downstream_exon.attributes.get('_line_number', 0)
            intron.metadata.line_number = (upstream_line + downstream_line) / 2

            yield intron

    def generate_from_transcript(
        self,
        transcript: Transcript,
        feature_index: Dict[str, Union[Gene, Transcript, Exon]]
    ) -> Iterator[Intron]:
        """
        Generate introns from all exons in a transcript.

        Implements two-pass algorithm matching original intronIC:
        1. Generate introns from CDS features (is_coding=True)
        2. Generate introns from exon features (is_coding=False)
        3. Only add exon-derived introns that don't overlap CDS-derived ones

        This ensures CDS-defined intron boundaries take priority, and
        exon-only introns (e.g., in UTR regions) fill in the gaps.

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
                    print(f"[!] Warning: Exon {child.feature_id} claims parent {child.parent_id} but is child of {transcript.feature_id}")

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
        non_redundant_introns = []

        # Pass 1: Generate introns from CDS features
        if cds_features:
            for intron in self.generate_from_exons(cds_features):
                intron.metadata.defined_by = 'cds'
                non_redundant_introns.append(intron)

        # Pass 2: Generate introns from exon features, excluding overlaps
        if exon_features:
            # Get coordinates of existing (CDS) introns for overlap checking
            existing_coords = [(i.coordinates.start, i.coordinates.stop)
                             for i in non_redundant_introns]

            for intron in self.generate_from_exons(exon_features):
                # Check if this exon-derived intron overlaps any CDS-derived intron
                intron_coords = (intron.coordinates.start, intron.coordinates.stop)
                if not self._check_intron_overlap(intron_coords, existing_coords):
                    # This is an exon-only intron (e.g., in UTR region)
                    intron.metadata.defined_by = 'exon'
                    non_redundant_introns.append(intron)

        if not non_redundant_introns:
            return

        # Set family size (total number of introns in transcript)
        family_size = len(non_redundant_introns)

        # Re-index introns sequentially (since they came from two separate generate calls)
        # Introns need to be sorted by genomic position first
        strand = non_redundant_introns[0].coordinates.strand if non_redundant_introns else '+'
        if strand == '-':
            # Negative strand: sort descending
            non_redundant_introns.sort(key=lambda i: i.coordinates.start, reverse=True)
        else:
            # Positive strand: sort ascending
            non_redundant_introns.sort(key=lambda i: i.coordinates.start)

        for index, intron in enumerate(non_redundant_introns, start=1):
            # Update intron metadata (1-based indexing to match original)
            intron.metadata.index = index  # Re-index sequentially
            intron.metadata.family_size = family_size
            intron.metadata.parent = transcript.feature_id
            intron.metadata.parent_length = coding_length  # Sum of CDS/exon lengths (not genomic span)

            # Set grandparent if available
            if transcript.parent_id and transcript.parent_id in feature_index:
                intron.metadata.grandparent = transcript.parent_id

            # Note: fractional_position is a computed property based on index and family_size
            # No need to set it manually

            yield intron

    @staticmethod
    def _check_intron_overlap(
        intron_coords: Tuple[int, int],
        existing_coords: List[Tuple[int, int]]
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
        self,
        gene: Gene,
        feature_index: Dict[str, Union[Gene, Transcript, Exon]]
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
        self,
        genes: List[Gene],
        feature_index: Dict[str, Union[Gene, Transcript, Exon]]
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
