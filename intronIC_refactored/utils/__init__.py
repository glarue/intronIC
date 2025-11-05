"""
Utility modules for intronIC.

This package contains utility functions and classes used throughout intronIC.
"""

from .coordinates import (
    GenomicCoordinate,
    bed_to_internal,
    internal_to_bed,
    gff_to_internal,
    internal_to_gff,
    validate_bed_coordinates,
    validate_gff_coordinates,
)

from .sequences import (
    reverse_complement,
    is_valid_dna,
    has_ambiguous_bases,
    gc_content,
    count_bases,
    extract_subsequence,
    normalize_sequence,
    sliding_window,
)

__all__ = [
    # Coordinates
    'GenomicCoordinate',
    'bed_to_internal',
    'internal_to_bed',
    'gff_to_internal',
    'internal_to_gff',
    'validate_bed_coordinates',
    'validate_gff_coordinates',
    # Sequences
    'reverse_complement',
    'is_valid_dna',
    'has_ambiguous_bases',
    'gc_content',
    'count_bases',
    'extract_subsequence',
    'normalize_sequence',
    'sliding_window',
]
