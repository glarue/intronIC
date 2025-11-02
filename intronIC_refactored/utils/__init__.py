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

__all__ = [
    'GenomicCoordinate',
    'bed_to_internal',
    'internal_to_bed',
    'gff_to_internal',
    'internal_to_gff',
    'validate_bed_coordinates',
    'validate_gff_coordinates',
]
