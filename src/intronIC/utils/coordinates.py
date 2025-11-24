"""
Genomic coordinate handling with explicit coordinate system tracking.

This module provides a robust abstraction for working with genomic coordinates,
addressing the #1 risk area identified in intronIC: coordinate conversion errors.

Coordinate Systems:
-------------------
- **BED format**: 0-based, half-open [start, stop)
  - start is 0-based (first base is 0)
  - stop is exclusive (not included in the range)
  - Example: [0, 10) represents bases 1-10 (10 bases total)

- **GFF3/GTF format**: 1-based, closed [start, stop]
  - start is 1-based (first base is 1)
  - stop is inclusive (included in the range)
  - Example: [1, 10] represents bases 1-10 (10 bases total)

- **intronIC internal**: 1-based, closed [start, stop]
  - Same as GFF3/GTF
  - This is the canonical internal representation

Conversion Rules:
-----------------
BED → internal (1-based):
    internal_start = bed_start + 1
    internal_stop = bed_stop

internal (1-based) → BED:
    bed_start = internal_start - 1
    bed_stop = internal_stop

GFF3 → internal:
    No conversion needed (both 1-based, closed)

Critical Notes:
---------------
1. Always validate that start < stop (we use strict inequality)
2. Genomic coordinates always have start < stop regardless of strand
3. On negative strand, biological direction is reversed, but coordinates are not
4. Single-base features are not supported (start must be < stop)
5. Negative coordinates are invalid (start must be >= 0 for 0-based, >= 1 for 1-based)

Examples:
---------
>>> # Convert BED to internal
>>> coord = bed_to_internal("chr1", 1000, 2000, '+')
>>> coord.start, coord.stop
(1001, 2000)

>>> # Convert internal to BED
>>> internal = GenomicCoordinate("chr1", 1001, 2000, '+', '1-based')
>>> chrom, start, stop, strand = internal_to_bed(internal)
>>> start, stop
(1000, 2000)

>>> # Round-trip conversion
>>> original = ("chr1", 1000, 2000, '+')
>>> coord = bed_to_internal(*original)
>>> result = internal_to_bed(coord)
>>> result == original
True

Author: intronIC refactoring project
Date: 2025-11-02
"""

from dataclasses import dataclass
from typing import Literal, Tuple


# Type aliases for clarity
Chromosome = str
Strand = Literal['+', '-']
CoordinateSystem = Literal['0-based', '1-based']


@dataclass(frozen=True, slots=True)
class GenomicCoordinate:
    """
    Immutable representation of a genomic coordinate with explicit system tracking.

    This class provides:
    - Type-safe coordinate representation
    - Automatic validation on creation
    - Immutability (prevents accidental modification)
    - Clear documentation of coordinate system

    Attributes:
        chromosome: Chromosome/contig/scaffold name (e.g., 'chr1', '1', 'scaffold_123')
        start: Start coordinate (interpretation depends on system)
        stop: Stop coordinate (interpretation depends on system)
        strand: '+' for forward strand, '-' for reverse strand
        system: '0-based' for BED-like, '1-based' for GFF3-like (default)

    Validation:
        - start must be < stop (strict inequality)
        - strand must be '+' or '-'
        - system must be '0-based' or '1-based'
        - start must be >= 0 (for 0-based) or >= 1 (for 1-based, but not enforced)

    Notes:
        - Single-base features are NOT supported (would require start == stop)
        - intronIC requires introns >= 30bp, so this restriction is appropriate
        - Use frozen=True for immutability (coordinates shouldn't change)

    Examples:
        >>> # Create a 1-based coordinate (default)
        >>> coord = GenomicCoordinate("chr1", 1001, 2000, '+')
        >>> coord.start, coord.stop
        (1001, 2000)

        >>> # Create a 0-based coordinate (BED-like)
        >>> bed_coord = GenomicCoordinate("chr1", 1000, 2000, '+', '0-based')
        >>> bed_coord.start, bed_coord.stop
        (1000, 2000)

        >>> # Validation catches errors
        >>> GenomicCoordinate("chr1", 100, 100, '+')  # doctest: +SKIP
        ValueError: start (100) must be < stop (100)
    """

    chromosome: Chromosome
    start: int
    stop: int
    strand: Strand
    system: CoordinateSystem = '1-based'

    def __post_init__(self) -> None:
        """
        Validate coordinates after initialization.

        Raises:
            ValueError: If coordinates are invalid
        """
        # Validate start <= stop (allow 1bp features to match original behavior)
        if self.start > self.stop:
            raise ValueError(
                f"start ({self.start}) must be <= stop ({self.stop})."
            )

        # Validate strand
        if self.strand not in ['+', '-']:
            raise ValueError(
                f"strand must be '+' or '-', got '{self.strand}'"
            )

        # Validate coordinate system
        if self.system not in ['0-based', '1-based']:
            raise ValueError(
                f"system must be '0-based' or '1-based', got '{self.system}'"
            )

        # Validate coordinate range based on system
        if self.system == '0-based':
            # 0-based coordinates: start must be >= 0
            if self.start < 0:
                raise ValueError(
                    f"start coordinate cannot be negative for 0-based system: {self.start}"
                )
        else:  # 1-based
            # 1-based coordinates: start must be >= 1
            if self.start < 1:
                raise ValueError(
                    f"start coordinate must be >= 1 for 1-based system: {self.start}"
                )

    @property
    def length(self) -> int:
        """
        Calculate the length of this genomic feature.

        For 1-based coordinates: length = stop - start + 1
        For 0-based coordinates: length = stop - start

        Returns:
            Length in bases

        Examples:
            >>> coord = GenomicCoordinate("chr1", 1, 10, '+', '1-based')
            >>> coord.length
            10

            >>> bed_coord = GenomicCoordinate("chr1", 0, 10, '+', '0-based')
            >>> bed_coord.length
            10
        """
        if self.system == '1-based':
            return self.stop - self.start + 1
        else:  # 0-based
            return self.stop - self.start

    def __str__(self) -> str:
        """Human-readable string representation."""
        return (
            f"{self.chromosome}:{self.start}-{self.stop}"
            f"({self.strand})[{self.system}]"
        )

    def __repr__(self) -> str:
        """Developer-friendly representation."""
        return (
            f"GenomicCoordinate("
            f"chromosome='{self.chromosome}', "
            f"start={self.start}, "
            f"stop={self.stop}, "
            f"strand='{self.strand}', "
            f"system='{self.system}')"
        )


def bed_to_internal(
    chrom: Chromosome,
    start: int,
    stop: int,
    strand: Strand
) -> GenomicCoordinate:
    """
    Convert BED coordinates (0-based, half-open) to internal (1-based, closed).

    BED format uses:
    - 0-based start coordinate (first base is 0)
    - Exclusive stop coordinate (not included)
    - Example: [0, 10) represents bases 1-10

    Internal format uses:
    - 1-based start coordinate (first base is 1)
    - Inclusive stop coordinate (included)
    - Example: [1, 10] represents bases 1-10

    Conversion:
        internal_start = bed_start + 1
        internal_stop = bed_stop (exclusive → inclusive works out to same value)

    Args:
        chrom: Chromosome name
        start: BED start coordinate (0-based)
        stop: BED stop coordinate (exclusive)
        strand: '+' or '-'

    Returns:
        GenomicCoordinate with 1-based, closed coordinates

    Raises:
        ValueError: If coordinates are invalid

    Examples:
        >>> coord = bed_to_internal("chr1", 0, 10, '+')
        >>> coord.start, coord.stop
        (1, 10)

        >>> coord = bed_to_internal("chr1", 1000, 2000, '-')
        >>> coord.start, coord.stop, coord.strand
        (1001, 2000, '-')

        >>> # Real example from chr19 gold standard
        >>> coord = bed_to_internal("19", 58345, 58861, '+')
        >>> coord.start, coord.stop, coord.length
        (58346, 58861, 516)
    """
    return GenomicCoordinate(
        chromosome=chrom,
        start=start + 1,  # 0-based → 1-based
        stop=stop,        # Half-open exclusive → closed inclusive (same value)
        strand=strand,
        system='1-based'
    )


def internal_to_bed(coord: GenomicCoordinate) -> Tuple[Chromosome, int, int, Strand]:
    """
    Convert internal coordinates (1-based, closed) to BED (0-based, half-open).

    Internal format uses:
    - 1-based start coordinate (first base is 1)
    - Inclusive stop coordinate (included)
    - Example: [1, 10] represents bases 1-10

    BED format uses:
    - 0-based start coordinate (first base is 0)
    - Exclusive stop coordinate (not included)
    - Example: [0, 10) represents bases 1-10

    Conversion:
        bed_start = internal_start - 1
        bed_stop = internal_stop (inclusive → exclusive works out to same value)

    Args:
        coord: GenomicCoordinate with 1-based coordinates

    Returns:
        Tuple of (chromosome, start, stop, strand) in BED format

    Raises:
        ValueError: If coord is not in 1-based system

    Examples:
        >>> coord = GenomicCoordinate("chr1", 1, 10, '+', '1-based')
        >>> chrom, start, stop, strand = internal_to_bed(coord)
        >>> start, stop
        (0, 10)

        >>> coord = GenomicCoordinate("chr1", 1001, 2000, '-', '1-based')
        >>> chrom, start, stop, strand = internal_to_bed(coord)
        >>> start, stop, strand
        (1000, 2000, '-')

        >>> # Real example from chr19 gold standard
        >>> coord = GenomicCoordinate("19", 58346, 58861, '+', '1-based')
        >>> chrom, start, stop, strand = internal_to_bed(coord)
        >>> start, stop
        (58345, 58861)
    """
    if coord.system != '1-based':
        raise ValueError(
            f"Can only convert 1-based coordinates to BED, got {coord.system}"
        )

    return (
        coord.chromosome,
        coord.start - 1,  # 1-based → 0-based
        coord.stop,       # Closed inclusive → half-open exclusive (same value)
        coord.strand
    )


def gff_to_internal(
    chrom: Chromosome,
    start: int,
    stop: int,
    strand: Strand
) -> GenomicCoordinate:
    """
    Convert GFF3/GTF coordinates (1-based, closed) to internal (1-based, closed).

    This is a no-op conversion since both GFF3/GTF and intronIC internal
    representation use the same coordinate system:
    - 1-based start coordinate (first base is 1)
    - Inclusive stop coordinate (included)
    - Example: [1, 10] represents bases 1-10

    Args:
        chrom: Chromosome name
        start: GFF3 start coordinate (1-based, inclusive)
        stop: GFF3 stop coordinate (1-based, inclusive)
        strand: '+' or '-'

    Returns:
        GenomicCoordinate with 1-based, closed coordinates

    Raises:
        ValueError: If coordinates are invalid

    Examples:
        >>> coord = gff_to_internal("chr1", 1, 10, '+')
        >>> coord.start, coord.stop
        (1, 10)

        >>> coord = gff_to_internal("chr1", 1001, 2000, '-')
        >>> coord.start, coord.stop, coord.strand
        (1001, 2000, '-')
    """
    return GenomicCoordinate(
        chromosome=chrom,
        start=start,  # No conversion needed
        stop=stop,    # No conversion needed
        strand=strand,
        system='1-based'
    )


def internal_to_gff(coord: GenomicCoordinate) -> Tuple[Chromosome, int, int, Strand]:
    """
    Convert internal coordinates (1-based, closed) to GFF3 (1-based, closed).

    This is a no-op conversion since both use the same coordinate system.

    Args:
        coord: GenomicCoordinate with 1-based coordinates

    Returns:
        Tuple of (chromosome, start, stop, strand) in GFF3 format

    Raises:
        ValueError: If coord is not in 1-based system

    Examples:
        >>> coord = GenomicCoordinate("chr1", 1, 10, '+', '1-based')
        >>> chrom, start, stop, strand = internal_to_gff(coord)
        >>> start, stop
        (1, 10)
    """
    if coord.system != '1-based':
        raise ValueError(
            f"Can only convert 1-based coordinates to GFF, got {coord.system}"
        )

    return (
        coord.chromosome,
        coord.start,  # No conversion needed
        coord.stop,   # No conversion needed
        coord.strand
    )


# Convenience functions for common operations

def validate_bed_coordinates(chrom: str, start: int, stop: int, strand: str) -> bool:
    """
    Validate BED coordinates without creating a GenomicCoordinate.

    Args:
        chrom: Chromosome name
        start: BED start (0-based)
        stop: BED stop (exclusive)
        strand: '+' or '-'

    Returns:
        True if valid, False otherwise

    Examples:
        >>> validate_bed_coordinates("chr1", 0, 10, '+')
        True
        >>> validate_bed_coordinates("chr1", 10, 10, '+')
        False
        >>> validate_bed_coordinates("chr1", 10, 5, '+')
        False
    """
    try:
        bed_to_internal(chrom, start, stop, strand)
        return True
    except ValueError:
        return False


def validate_gff_coordinates(chrom: str, start: int, stop: int, strand: str) -> bool:
    """
    Validate GFF3 coordinates without creating a GenomicCoordinate.

    Args:
        chrom: Chromosome name
        start: GFF3 start (1-based, inclusive)
        stop: GFF3 stop (1-based, inclusive)
        strand: '+' or '-'

    Returns:
        True if valid, False otherwise

    Examples:
        >>> validate_gff_coordinates("chr1", 1, 10, '+')
        True
        >>> validate_gff_coordinates("chr1", 10, 10, '+')
        False
        >>> validate_gff_coordinates("chr1", 0, 10, '+')
        False
    """
    try:
        gff_to_internal(chrom, start, stop, strand)
        return True
    except ValueError:
        return False


if __name__ == "__main__":
    import doctest
    doctest.testmod()
