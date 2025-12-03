"""
Extract sequences for introns from genome files.

This module handles extracting intron sequences, flanking exonic sequences,
and scoring region sequences (5'ss, BP, 3'ss) from genome FASTA files.
"""

from collections import defaultdict
from typing import Dict, Iterator, List, Optional, Tuple

from intronIC.core.intron import Intron, IntronSequences
from intronIC.file_io.genome import GenomeReader
from intronIC.utils.sequences import reverse_complement

# Canonical splice site dinucleotides
CANONICAL_FIVE_PRIME = {"GT", "GC"}
CANONICAL_THREE_PRIME = {"AG"}
CANONICAL_PAIRS = {
    ("GT", "AG"),
    ("GC", "AG"),  # Major spliceosome
    ("AT", "AC"),  # U12-type
}


class SequenceExtractor:
    """
    Extracts sequences for introns from genome data.

    This class handles:
    - Extracting full intron sequences
    - Extracting flanking exonic sequences
    - Extracting scoring region sequences (5'ss, BP, 3'ss)
    - Identifying terminal dinucleotides
    - Checking for canonical splice sites

    Examples:
        >>> extractor = SequenceExtractor('genome.fa')
        >>> for intron in extractor.extract_sequences(introns):
        ...     print(f"{intron.name}: {len(intron.sequences.seq)} bp")
    """

    def __init__(self, genome_file: str, use_cache: bool = True):
        """
        Initialize the sequence extractor.

        Args:
            genome_file: Path to genome FASTA file
            use_cache: Whether to cache genome in memory (default: True)
        """
        self.genome_file = genome_file
        self.genome_reader = GenomeReader(genome_file, cached=use_cache)
        self.use_cache = use_cache

    @classmethod
    def from_indexed_reader(
        cls, genome_file: str, indexed_reader
    ) -> "SequenceExtractor":
        """
        Create a SequenceExtractor using an IndexedGenomeReader.

        This is useful for streaming mode where we want to use indexed
        random access instead of loading entire chromosomes.

        Args:
            genome_file: Path to genome FASTA file (for reference)
            indexed_reader: An IndexedGenomeReader instance with fetch() method

        Returns:
            SequenceExtractor configured for indexed access
        """
        instance = cls.__new__(cls)
        instance.genome_file = genome_file
        instance.genome_reader = indexed_reader
        instance.use_cache = False  # Use fetch() mode
        return instance

    def extract_sequences(
        self,
        introns: List[Intron],
        flank_size: int | Tuple[int, int] = 200,
        five_score_coords: Tuple[int, int] = (-3, 9),
        three_score_coords: Tuple[int, int] = (-6, 4),
        bp_coords: Tuple[int, int] = (-55, -5),
    ) -> Iterator[Intron]:
        """
        Extract sequences for a list of introns.

        Args:
            introns: List of Intron objects
            flank_size: Size of flanking regions (int or (upstream, downstream) tuple)
            five_score_coords: Coordinates for 5' splice site scoring region
            three_score_coords: Coordinates for 3' splice site scoring region
            bp_coords: Coordinates for branch point scoring region

        Yields:
            Intron objects with sequences populated

        Examples:
            >>> extractor = SequenceExtractor('genome.fa')
            >>> introns_with_seqs = list(
            ...     extractor.extract_sequences(intron_list, flank_size=100)
            ... )
        """
        # Group introns by region for efficient genome access
        introns_by_region = self._group_by_region(introns)

        # Process each chromosome/region
        for region_name, region_introns in introns_by_region.items():
            # For cached mode: load entire chromosome once
            # For indexed mode: we'll use fetch() per intron instead
            region_seq = None
            if self.use_cache:
                try:
                    region_seq = self.genome_reader.get_sequence(region_name).upper()
                except KeyError:
                    print(
                        f"[!] Warning: Region '{region_name}' not found in genome, skipping"
                    )
                    continue

            # Process each intron in this region
            for intron in region_introns:
                # Extract sequences and populate intron
                intron = self._extract_intron_sequences(
                    intron,
                    region_seq,
                    flank_size,
                    five_score_coords,
                    three_score_coords,
                    bp_coords,
                )

                yield intron

    def extract_sequences_with_deduplication(
        self,
        introns: List[Intron],
        flank_size: int | Tuple[int, int] = 200,
        five_score_coords: Tuple[int, int] = (-3, 9),
        three_score_coords: Tuple[int, int] = (-6, 4),
        bp_coords: Tuple[int, int] = (-55, -5),
    ) -> Iterator[Intron]:
        """
        Extract sequences with deduplication - reuse sequences for duplicate coordinates.

        This method provides significant performance and memory benefits by extracting
        sequences only once per unique coordinate set and reusing them for all introns
        with the same coordinates.

        Memory savings: For datasets with ~10% duplicates, saves ~10% of extraction time
        and memory during extraction phase.

        Args:
            introns: List of Intron objects
            flank_size: Size of flanking regions (int or (upstream, downstream) tuple)
            five_score_coords: Coordinates for 5' splice site scoring region
            three_score_coords: Coordinates for 3' splice site scoring region
            bp_coords: Coordinates for branch point scoring region

        Yields:
            Intron objects with sequences populated

        Examples:
            >>> extractor = SequenceExtractor('genome.fa')
            >>> introns_with_seqs = list(
            ...     extractor.extract_sequences_with_deduplication(intron_list)
            ... )
        """
        # Group introns by region for efficient genome access
        introns_by_region = self._group_by_region(introns)

        # Process each chromosome/region
        for region_name, region_introns in introns_by_region.items():
            # Get region sequence
            try:
                region_seq = self.genome_reader.get_sequence(region_name).upper()
            except KeyError:
                print(
                    f"[!] Warning: Region '{region_name}' not found in genome, skipping"
                )
                continue

            # Group introns by coordinates within this region
            # Key: (start, stop, strand)
            # Value: List of introns with same coordinates
            coord_groups = defaultdict(list)
            for intron in region_introns:
                coord_key = (
                    intron.coordinates.start,
                    intron.coordinates.stop,
                    intron.coordinates.strand,
                )
                coord_groups[coord_key].append(intron)

            # Process each coordinate group
            for coord_key, duplicate_introns in coord_groups.items():
                # Extract sequences for first intron in group
                first_intron = duplicate_introns[0]
                first_with_seqs = self._extract_intron_sequences(
                    first_intron,
                    region_seq,
                    flank_size,
                    five_score_coords,
                    three_score_coords,
                    bp_coords,
                )

                # Yield first intron with extracted sequences
                yield first_with_seqs

                # For remaining introns (duplicates), reuse the sequences
                if len(duplicate_introns) > 1:
                    shared_sequences = first_with_seqs.sequences
                    # Get noncanonical status from first intron to propagate to duplicates
                    is_noncanonical = (
                        first_with_seqs.metadata.noncanonical
                        if first_with_seqs.metadata
                        else False
                    )
                    for dup_intron in duplicate_introns[1:]:
                        # Attach same sequences object to duplicate intron
                        dup_with_seqs = dup_intron.with_sequences(shared_sequences)
                        # Propagate noncanonical flag to duplicate's metadata
                        if dup_with_seqs.metadata:
                            dup_with_seqs.metadata.noncanonical = is_noncanonical
                            # Also propagate dynamic tag for non-canonical
                            if is_noncanonical:
                                dup_with_seqs.metadata.dynamic_tags.add("[n]")
                            else:
                                dup_with_seqs.metadata.dynamic_tags.discard("[n]")
                        yield dup_with_seqs

    def _group_by_region(self, introns: List[Intron]) -> Dict[str, List[Intron]]:
        """
        Group introns by chromosome/region.

        Args:
            introns: List of Intron objects

        Returns:
            Dictionary mapping region names to introns
        """
        grouped = defaultdict(list)
        for intron in introns:
            grouped[intron.coordinates.chromosome].append(intron)
        return grouped

    def _extract_intron_sequences(
        self,
        intron: Intron,
        region_seq: Optional[str],
        flank_size: int | Tuple[int, int],
        five_score_coords: Tuple[int, int],
        three_score_coords: Tuple[int, int],
        bp_coords: Tuple[int, int],
    ) -> Intron:
        """
        Extract all sequences for a single intron.

        Args:
            intron: Intron object
            region_seq: Full chromosome/region sequence (None for indexed mode)
            flank_size: Flanking sequence size(s)
            five_score_coords: 5'ss scoring coordinates
            three_score_coords: 3'ss scoring coordinates
            bp_coords: Branch point scoring coordinates

        Returns:
            Intron with sequences populated
        """
        coord = intron.coordinates

        # Handle flank sizes
        if isinstance(flank_size, int):
            upstream_flank_size = downstream_flank_size = flank_size
        else:
            upstream_flank_size, downstream_flank_size = flank_size

        # Extract intron sequence with flanks
        # For cached mode: slice from region_seq
        # For indexed mode: use fetch() to get just what we need
        if region_seq is not None:
            # Cached mode: use provided region sequence
            # Coordinates are 1-based, Python slicing is 0-based
            start_idx = coord.start - 1  # Convert to 0-based
            stop_idx = coord.stop  # Stop is inclusive in 1-based, exclusive in 0-based

            # Calculate flank boundaries
            upstream_start = max(0, start_idx - upstream_flank_size)
            downstream_end = min(len(region_seq), stop_idx + downstream_flank_size)

            # Extract full sequence with flanks
            full_seq = region_seq[upstream_start:downstream_end]

            # Calculate actual boundaries of intron within full_seq
            intron_start_in_full = start_idx - upstream_start
            intron_stop_in_full = intron_start_in_full + (stop_idx - start_idx)
        else:
            # Indexed mode: fetch just the region we need (with flanks)
            # Calculate boundaries in 1-based coordinates
            fetch_start = max(1, coord.start - upstream_flank_size)
            fetch_stop = coord.stop + downstream_flank_size

            # Fetch from indexed genome (handles coordinate conversion internally)
            full_seq = self.genome_reader.fetch(
                coord.chromosome, fetch_start, fetch_stop
            ).upper()

            # For indexed mode, calculate offset differently
            # fetch_start is the actual start position in 1-based coords
            # full_seq starts at fetch_start
            actual_upstream_size = coord.start - fetch_start
            intron_length = (
                coord.stop - coord.start + 1
            )  # +1 because both are inclusive

            intron_start_in_full = actual_upstream_size
            intron_stop_in_full = intron_start_in_full + intron_length

        # Split into components
        upstream_flank = full_seq[:intron_start_in_full]
        intron_seq = full_seq[intron_start_in_full:intron_stop_in_full]
        downstream_flank = full_seq[intron_stop_in_full:]

        # Handle reverse complement for negative strand
        if coord.strand == "-":
            upstream_flank = reverse_complement(upstream_flank)
            intron_seq = reverse_complement(intron_seq)
            downstream_flank = reverse_complement(downstream_flank)
            # Swap flanks for negative strand
            upstream_flank, downstream_flank = downstream_flank, upstream_flank

        # Extract scoring sequences
        scoring_seq = upstream_flank + intron_seq + downstream_flank
        upstream_len = len(upstream_flank)
        downstream_len = len(downstream_flank)

        # Adjust coordinates for scoring regions (relative to scoring_seq)
        five_start = five_score_coords[0] + upstream_len
        five_end = five_score_coords[1] + upstream_len
        three_start = three_score_coords[0] - downstream_len
        three_end = three_score_coords[1] - downstream_len

        # Handle zero indices (Python slicing None = end/start)
        five_start = five_start if five_start != 0 else None
        five_end = five_end if five_end != 0 else None
        three_start = three_start if three_start != 0 else None
        three_end = three_end if three_end != 0 else None

        five_seq = scoring_seq[five_start:five_end]
        three_seq = scoring_seq[three_start:three_end]

        # Extract branch point region
        # BP coords are relative to 3' splice site (acceptor)
        # Adjust to be 0-based and relative to intron end
        bp_start_rel = bp_coords[0]  # Negative value (e.g., -55)
        bp_end_rel = bp_coords[1]  # Negative value (e.g., -5)

        if coord.strand == "+":
            # For + strand: 3' splice site is at coord.stop
            # bp_start_rel is more negative (e.g., -55), so further upstream
            # bp_end_rel is less negative (e.g., -5), so closer to 3'ss
            bp_start_abs = coord.stop + bp_start_rel - 1  # -1 for 0-based
            bp_end_abs = coord.stop + bp_end_rel - 1
        else:
            # For - strand: 3' splice site is at coord.start
            # On - strand, the coordinates are reversed: more negative means further
            # in genomic coords (higher position), so we swap start/end
            # This matches v1.5.1: bpc = (bp_stop, bp_start) for reverse strand
            bp_end_abs = coord.start - bp_start_rel - 1  # Note: assigned to end
            bp_start_abs = coord.start - bp_end_rel - 1  # Note: assigned to start

        # Extract BP region (different handling for cached vs indexed mode)
        if region_seq is not None:
            # Cached mode: slice from region_seq
            bp_start_abs = max(0, bp_start_abs)
            bp_end_abs = min(len(region_seq), bp_end_abs)
            bp_region_seq = region_seq[bp_start_abs:bp_end_abs]
        else:
            # Indexed mode: fetch from genome
            # Convert back to 1-based for fetch
            bp_fetch_start = max(1, bp_start_abs + 1)
            bp_fetch_stop = bp_end_abs + 1

            bp_region_seq = self.genome_reader.fetch(
                coord.chromosome, bp_fetch_start, bp_fetch_stop
            ).upper()

        if coord.strand == "-":
            bp_region_seq = reverse_complement(bp_region_seq)

        # Identify terminal dinucleotides
        five_prime_dnt = intron_seq[:2] if len(intron_seq) >= 2 else intron_seq
        three_prime_dnt = intron_seq[-2:] if len(intron_seq) >= 2 else intron_seq

        # Check if canonical
        is_canonical = self._is_canonical(five_prime_dnt, three_prime_dnt)

        # Display sequences for motif schematic (following original implementation)
        # five_display_seq: first 10bp of intron for 5' boundary display
        five_display_length = 10
        five_display_seq = (
            intron_seq[:five_display_length]
            if len(intron_seq) >= five_display_length
            else intron_seq
        )

        # three_display_seq: from bp search region end to intron end
        # bp_coords[1] is relative to 3' end (e.g., -5 means 5bp from end)
        # This will be used for 3' boundary display in motif schematic
        three_display_seq = (
            intron_seq[bp_coords[1] :]
            if len(intron_seq) >= abs(bp_coords[1])
            else intron_seq
        )

        # Create sequences object
        sequences = IntronSequences(
            seq=intron_seq,
            upstream_flank=upstream_flank,
            downstream_flank=downstream_flank,
            five_seq=five_seq,
            three_seq=three_seq,
            bp_region_seq=bp_region_seq,
            five_prime_dnt=five_prime_dnt,
            three_prime_dnt=three_prime_dnt,
            five_display_seq=five_display_seq,
            three_display_seq=three_display_seq,
            # bp_seq_u2 and bp_relative_coords populated during PWM scoring
        )

        # Update intron metadata
        intron.metadata.noncanonical = not is_canonical

        # Update dynamic tag for non-canonical introns
        if intron.metadata.noncanonical:
            intron.metadata.dynamic_tags.add("[n]")
        else:
            # Remove [n] tag if intron is canonical (e.g., after correction)
            intron.metadata.dynamic_tags.discard("[n]")

        # Return intron with sequences
        return intron.with_sequences(sequences)

    @staticmethod
    def _is_canonical(five_dnt: str, three_dnt: str) -> bool:
        """
        Check if dinucleotides form canonical splice sites.

        Args:
            five_dnt: 5' splice site dinucleotide
            three_dnt: 3' splice site dinucleotide

        Returns:
            True if canonical, False otherwise

        Examples:
            >>> SequenceExtractor._is_canonical('GT', 'AG')
            True
            >>> SequenceExtractor._is_canonical('GC', 'AG')
            True
            >>> SequenceExtractor._is_canonical('AT', 'AC')
            True
            >>> SequenceExtractor._is_canonical('GG', 'AG')
            False
        """
        return (five_dnt, three_dnt) in CANONICAL_PAIRS

    def close(self):
        """Close the genome reader."""
        self.genome_reader.close()

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()


def extract_sequences_for_introns(
    introns: List[Intron],
    genome_file: str,
    flank_size: int = 200,
    five_score_coords: Tuple[int, int] = (-3, 9),
    three_score_coords: Tuple[int, int] = (-6, 4),
    bp_coords: Tuple[int, int] = (-55, -5),
) -> Iterator[Intron]:
    """
    Convenience function to extract sequences for introns.

    This is a functional wrapper for backwards compatibility with
    the original intronIC API.

    Args:
        introns: List of Intron objects
        genome_file: Path to genome FASTA file
        flank_size: Size of flanking regions
        five_score_coords: 5' splice site scoring coordinates
        three_score_coords: 3' splice site scoring coordinates
        bp_coords: Branch point scoring coordinates

    Yields:
        Intron objects with sequences populated

    Examples:
        >>> introns_with_seqs = list(
        ...     extract_sequences_for_introns(introns, 'genome.fa')
        ... )
    """
    with SequenceExtractor(genome_file) as extractor:
        yield from extractor.extract_sequences(
            introns, flank_size, five_score_coords, three_score_coords, bp_coords
        )
