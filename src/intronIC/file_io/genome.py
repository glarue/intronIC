"""
Genome file reading and sequence extraction.

This module provides efficient FASTA file parsing with support for:
- Streaming mode (memory efficient, one chromosome at a time)
- Cached mode (load all sequences into memory for fast repeated access)
- Gzip compression support
- Subsequence extraction by genomic coordinates

Design considerations:
- Memory efficiency: Don't load entire genome unless requested
- Performance: Cache sequences when beneficial
- Flexibility: Support both streaming and random access patterns

Author: intronIC refactoring project
Date: 2025-11-02
"""

from pathlib import Path
from typing import Iterator, Tuple, Dict, Optional, Union
from smart_open import open as smart_open
from intronIC.utils.coordinates import GenomicCoordinate


def parse_fasta(file_path: Union[str, Path]) -> Iterator[Tuple[str, str]]:
    """
    Parse a FASTA file (optionally gzipped) and yield (name, sequence) tuples.

    This is a streaming parser that yields one sequence at a time to minimize
    memory usage. Suitable for processing large genome files.

    Args:
        file_path: Path to FASTA file (.fa, .fasta, .fa.gz, .fasta.gz)

    Yields:
        Tuple of (sequence_name, sequence_string)
        - sequence_name: The header line without '>'
        - sequence_string: The assembled sequence (stripped of whitespace)

    Examples:
        >>> # This would work with a real FASTA file
        >>> # for name, seq in parse_fasta("genome.fa"):
        >>> #     print(f"{name}: {len(seq)} bp")

    Note:
        - Automatically detects and handles gzip compression
        - Strips whitespace from sequence lines
        - Converts sequences to uppercase
    """
    file_path = Path(file_path)

    with smart_open(file_path, 'rt') as f:
        name = None
        seq_lines = []

        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith('>'):
                # Yield previous sequence if exists
                if name is not None:
                    yield name, ''.join(seq_lines).upper()

                # Start new sequence
                name = line[1:].split()[0]  # Take first word after '>'
                seq_lines = []
            else:
                # Accumulate sequence lines
                seq_lines.append(line)

        # Yield final sequence
        if name is not None:
            yield name, ''.join(seq_lines).upper()


class GenomeReader:
    """
    Efficient genome sequence reader with optional caching.

    Provides two modes of operation:
    1. Streaming mode (default): Memory efficient, processes one chromosome at a time
    2. Cached mode: Loads all sequences into memory for fast repeated access

    Attributes:
        file_path: Path to genome FASTA file
        cache: Dictionary of chromosome_name -> sequence (if cached=True)
        is_cached: Whether sequences are loaded into memory

    Examples:
        >>> # Streaming mode (memory efficient)
        >>> reader = GenomeReader("genome.fa")
        >>> for chrom, seq in reader.stream():
        ...     print(f"{chrom}: {len(seq)} bp")

        >>> # Cached mode (faster repeated access)
        >>> reader = GenomeReader("genome.fa", cached=True)
        >>> seq = reader.get_sequence("chr1")
        >>> print(f"chr1: {len(seq)} bp")

        >>> # Extract subsequence using coordinates
        >>> from utils.coordinates import GenomicCoordinate
        >>> coord = GenomicCoordinate("chr1", 1000, 2000, '+', '1-based')
        >>> subseq = reader.extract_subsequence(coord)
    """

    def __init__(self, file_path: Union[str, Path], cached: bool = False):
        """
        Initialize GenomeReader.

        Args:
            file_path: Path to FASTA file
            cached: If True, load all sequences into memory on init

        Raises:
            FileNotFoundError: If file doesn't exist
        """
        self.file_path = Path(file_path)

        if not self.file_path.exists():
            raise FileNotFoundError(f"Genome file not found: {self.file_path}")

        self.cache: Optional[Dict[str, str]] = None
        self.is_cached = cached

        if cached:
            self._load_cache()

    def _load_cache(self) -> None:
        """Load all sequences into memory."""
        self.cache = {}
        for name, seq in parse_fasta(self.file_path):
            self.cache[name] = seq
        self.is_cached = True

    def load_sequences(self, sequence_names: list[str]) -> None:
        """
        Load only specific sequences into cache (selective caching).

        This is much more memory-efficient than loading the entire genome
        when you only need a subset of contigs/chromosomes.

        Args:
            sequence_names: List of sequence names to load

        Examples:
            >>> reader = GenomeReader("genome.fa", cached=False)
            >>> reader.load_sequences(["chr1", "chr2"])  # Only load chr1 and chr2
            >>> seq = reader.get_sequence("chr1")  # Now works!

        Note:
            - Stops reading the file once all requested sequences are found
            - Overwrites any existing cache
            - Sets is_cached=True after loading
        """
        if not sequence_names:
            return

        self.cache = {}
        target_names = set(sequence_names)

        for name, seq in parse_fasta(self.file_path):
            if name in target_names:
                self.cache[name] = seq
                # Stop early if we've found all requested sequences
                if len(self.cache) == len(target_names):
                    break

        self.is_cached = True

    def stream(self) -> Iterator[Tuple[str, str]]:
        """
        Stream sequences from genome file (memory efficient).

        Yields:
            Tuple of (chromosome_name, sequence)

        Note:
            If the reader is cached, yields from cache instead of file.
            This maintains consistent API regardless of caching mode.
        """
        if self.is_cached and self.cache is not None:
            yield from self.cache.items()
        else:
            yield from parse_fasta(self.file_path)

    def get_sequence(self, chromosome: str) -> str:
        """
        Get full sequence for a chromosome.

        Args:
            chromosome: Chromosome/contig name

        Returns:
            Complete chromosome sequence

        Raises:
            KeyError: If chromosome not found
            RuntimeError: If not in cached mode

        Note:
            Requires cached mode. For streaming access, use stream() instead.
        """
        if not self.is_cached or self.cache is None:
            raise RuntimeError(
                "get_sequence() requires cached mode. "
                "Use GenomeReader(file_path, cached=True) or use stream() instead."
            )

        if chromosome not in self.cache:
            raise KeyError(f"Chromosome '{chromosome}' not found in genome")

        return self.cache[chromosome]

    def has_chromosome(self, chromosome: str) -> bool:
        """
        Check if chromosome exists in genome.

        Args:
            chromosome: Chromosome/contig name

        Returns:
            True if chromosome exists

        Note:
            Requires cached mode for O(1) lookup. In streaming mode,
            would require full scan of file.
        """
        if not self.is_cached or self.cache is None:
            raise RuntimeError(
                "has_chromosome() requires cached mode. "
                "Use GenomeReader(file_path, cached=True)."
            )

        return chromosome in self.cache

    def extract_subsequence(
        self,
        coordinate: GenomicCoordinate,
        upstream_flank: int = 0,
        downstream_flank: int = 0
    ) -> str:
        """
        Extract a subsequence from the genome using genomic coordinates.

        Args:
            coordinate: GenomicCoordinate specifying location
            upstream_flank: Number of bases to include upstream
            downstream_flank: Number of bases to include downstream

        Returns:
            Extracted sequence (reverse complemented if on negative strand)

        Raises:
            KeyError: If chromosome not found
            RuntimeError: If not in cached mode
            ValueError: If coordinates are out of bounds

        Examples:
            >>> reader = GenomeReader("genome.fa", cached=True)
            >>> coord = GenomicCoordinate("chr1", 1000, 2000, '+', '1-based')
            >>> seq = reader.extract_subsequence(coord)
            >>> len(seq)
            1001

            >>> # With flanking sequence
            >>> seq_with_flanks = reader.extract_subsequence(coord, upstream_flank=100, downstream_flank=100)
            >>> len(seq_with_flanks)
            1201

        Note:
            - Uses 1-based coordinates (GenomicCoordinate default)
            - Automatically reverse complements for negative strand
            - Flanking regions are relative to strand direction
        """
        from intronIC.utils.sequences import reverse_complement

        if not self.is_cached or self.cache is None:
            raise RuntimeError(
                "extract_subsequence() requires cached mode. "
                "Use GenomeReader(file_path, cached=True)."
            )

        # Get chromosome sequence
        chrom_seq = self.get_sequence(coordinate.chromosome)

        # Convert to 0-based for string slicing
        # 1-based [start, stop] → 0-based [start-1, stop)
        start_0based = coordinate.start - 1
        stop_0based = coordinate.stop  # stop is inclusive in 1-based, becomes exclusive in 0-based

        # Add flanking regions
        flank_start = max(0, start_0based - upstream_flank)
        flank_stop = min(len(chrom_seq), stop_0based + downstream_flank)

        # Validate coordinates
        if start_0based < 0 or stop_0based > len(chrom_seq):
            raise ValueError(
                f"Coordinates out of bounds: {coordinate.chromosome}:{coordinate.start}-{coordinate.stop} "
                f"(chromosome length: {len(chrom_seq)})"
            )

        # Extract sequence
        seq = chrom_seq[flank_start:flank_stop]

        # Reverse complement if on negative strand
        if coordinate.strand == '-':
            seq = reverse_complement(seq)

        return seq

    def get_chromosome_names(self) -> list:
        """
        Get list of all chromosome names in genome.

        Returns:
            List of chromosome names

        Note:
            Requires cached mode.
        """
        if not self.is_cached or self.cache is None:
            raise RuntimeError(
                "get_chromosome_names() requires cached mode. "
                "Use GenomeReader(file_path, cached=True)."
            )

        return list(self.cache.keys())

    def get_chromosome_length(self, chromosome: str) -> int:
        """
        Get length of a chromosome.

        Args:
            chromosome: Chromosome name

        Returns:
            Length in bases

        Note:
            Requires cached mode.
        """
        return len(self.get_sequence(chromosome))

    def __repr__(self) -> str:
        """Developer-friendly representation."""
        mode = "cached" if self.is_cached else "streaming"
        n_chroms = len(self.cache) if self.cache else "unknown"
        return f"GenomeReader('{self.file_path}', mode={mode}, chromosomes={n_chroms})"

    def __str__(self) -> str:
        """Human-readable representation."""
        mode = "cached" if self.is_cached else "streaming"
        return f"GenomeReader[{mode}]: {self.file_path.name}"


if __name__ == "__main__":
    import doctest
    doctest.testmod()
