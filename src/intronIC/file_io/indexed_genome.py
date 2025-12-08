"""
Indexed FASTA genome reader for memory-efficient parallel processing.

This module provides random-access to genome sequences using indexed FASTA files,
avoiding the need to load the entire genome into memory. Perfect for parallel
processing where multiple workers need genome access.

Uses pyfastx for efficient indexed access (supports standard gzip files).

Author: intronIC refactoring project
Date: 2025-11-21
"""

from pathlib import Path
from typing import Optional, Union


class IndexedGenomeReader:
    """
    Memory-efficient genome reader using indexed FASTA.

    This class provides random access to genome sequences without loading
    the entire genome into memory. Each instance uses minimal memory (~few MB)
    making it perfect for parallel processing.

    The genome file is automatically indexed on first access.
    Supports both plain and gzipped FASTA files (standard gzip or bgzip).

    Attributes:
        file_path: Path to FASTA file
        fasta: pyfastx.Fasta handle (lazy-loaded)

    Examples:
        >>> reader = IndexedGenomeReader("genome.fa")
        >>> seq = reader.fetch("chr1", 1000, 2000)  # 1-based coordinates
        >>> len(seq)
        1000

    Memory Profile:
        - Without cache: ~5-10 MB per instance (just file handle + index)
        - With cache: ~250 MB (full genome in RAM, like old approach)
        - Recommended: No cache for parallel processing

    Performance:
        - Random access: ~0.1-1ms per fetch (depends on disk speed)
        - OS page cache: Frequently accessed regions are cached by kernel
        - For ~200K introns: Total fetch time ~20-200 seconds (trivial)
    """

    def __init__(self, file_path: Union[str, Path], use_cache: bool = False):
        """
        Initialize indexed genome reader.

        Args:
            file_path: Path to FASTA file (.fa, .fasta, optionally .gz)
            use_cache: If True, load entire genome into memory (not recommended)

        Note:
            The index file will be created automatically if missing.
            Supports both standard gzip (.gz) and plain FASTA files.
        """
        self.file_path = Path(file_path)
        self._fasta = None  # Lazy-loaded
        self._cache = {} if use_cache else None
        self._index_status_logged = False  # Track if we've logged index status

        if not self.file_path.exists():
            raise FileNotFoundError(f"Genome file not found: {self.file_path}")

    def _ensure_index(self) -> None:
        """Check if index exists (logging handled by main process)."""
        # Index status is logged by the main process before workers start
        # Workers silently reuse or create the index as needed
        pass

    @property
    def fasta(self):
        """Lazy-load pyfastx.Fasta handle."""
        if self._fasta is None:
            try:
                import pyfastx
            except ImportError:
                raise ImportError(
                    "pyfastx is required for indexed genome access. "
                    "Install with: pixi add --pypi pyfastx"
                )

            # Log index status on first access
            if not self._index_status_logged:
                self._ensure_index()
                self._index_status_logged = True

            # pyfastx automatically creates index if missing
            try:
                self._fasta = pyfastx.Fasta(str(self.file_path))
            except Exception as e:
                raise RuntimeError(
                    f"Failed to load genome file '{self.file_path}' with pyfastx. "
                    f"Error: {e}. "
                    f"This may happen if the file is corrupted, not a valid FASTA, "
                    f"or if there are permission issues."
                ) from e
        return self._fasta

    def fetch(self, chromosome: str, start: int, stop: int) -> str:
        """
        Fetch sequence from genome using 1-based coordinates.

        Args:
            chromosome: Chromosome/contig name
            start: Start position (1-based, inclusive)
            stop: Stop position (1-based, inclusive)

        Returns:
            Uppercase sequence string

        Note:
            Input uses 1-based coordinates (intronIC internal format).
            pyfastx.fetch() also uses 1-based coordinates.
        """
        # Check cache first
        if self._cache is not None:
            key = (chromosome, start, stop)
            if key in self._cache:
                return self._cache[key]

        # pyfastx fetch() uses 1-based coordinates
        # Returns string directly
        try:
            seq = self.fasta.fetch(chromosome, (start, stop))

            # Handle case where fetch returns None
            if seq is None:
                raise ValueError(
                    f"pyfastx returned None for {chromosome}:{start}-{stop}. "
                    f"This may indicate chromosome not found in genome or invalid coordinates."
                )

            seq = seq.upper()
        except Exception as e:
            raise RuntimeError(
                f"Failed to fetch sequence from {chromosome}:{start}-{stop}. "
                f"Error: {e}. "
                f"Check that chromosome name exists in genome and coordinates are valid."
            ) from e

        # Cache if enabled
        if self._cache is not None:
            self._cache[(chromosome, start, stop)] = seq

        return seq

    def get_sequence(self, chromosome: str) -> str:
        """
        Get entire chromosome sequence.

        Args:
            chromosome: Chromosome/contig name

        Returns:
            Full chromosome sequence (uppercase)

        Warning:
            This loads the entire chromosome into memory.
            Use fetch() for subsequences instead.
        """
        # pyfastx returns sequence object, use .seq to get string
        return self.fasta[chromosome].seq.upper()

    def close(self) -> None:
        """Close FASTA file handle (pyfastx.Fasta auto-manages resources)."""
        if self._fasta is not None:
            # pyfastx.Fasta doesn't have a close() method - it auto-manages resources
            # Just clear the reference to allow garbage collection
            self._fasta = None

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()

    def __del__(self):
        """Cleanup on deletion."""
        # Only cleanup if not already done
        if hasattr(self, '_fasta'):
            self.close()


# Module-level global for worker processes
_WORKER_GENOME = None


def init_worker_genome(genome_path: str) -> None:
    """
    Initialize genome reader in worker process.

    This is called once per worker via Pool(initializer=...).
    Creates a lightweight indexed reader that workers can use.

    Args:
        genome_path: Path to indexed FASTA file
    """
    global _WORKER_GENOME
    _WORKER_GENOME = IndexedGenomeReader(genome_path, use_cache=False)


def get_worker_genome() -> IndexedGenomeReader:
    """
    Get the worker's genome reader.

    Returns:
        IndexedGenomeReader instance for this worker

    Raises:
        RuntimeError: If worker genome not initialized
    """
    global _WORKER_GENOME
    if _WORKER_GENOME is None:
        raise RuntimeError(
            "Worker genome not initialized. "
            "Use Pool(initializer=init_worker_genome, initargs=(genome_path,))"
        )
    return _WORKER_GENOME


def get_contig_lengths(genome_path: Union[str, Path]) -> dict[str, int]:
    """
    Get contig lengths from pyfastx index (no extra file scan).

    This reads metadata from the .fxi index file that was already created,
    so there's no additional overhead beyond opening the index.

    Args:
        genome_path: Path to FASTA file (index must already exist)

    Returns:
        Dictionary mapping contig names to their lengths in base pairs

    Examples:
        >>> lengths = get_contig_lengths("genome.fa")
        >>> print(f"chr1: {lengths['chr1']:,} bp")
        chr1: 248,956,422 bp

    Note:
        Requires pyfastx index (.fxi file) to exist.
        This is typically created automatically on first access.
    """
    try:
        import pyfastx
    except ImportError:
        raise ImportError(
            "pyfastx is required for indexed genome access. "
            "Install with: pixi add --pypi pyfastx"
        )

    genome_path = Path(genome_path)
    fa = pyfastx.Fasta(str(genome_path))

    # Get lengths for all contigs in the genome
    return {name: len(fa[name]) for name in fa.keys()}
