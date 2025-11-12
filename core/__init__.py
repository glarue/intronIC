"""
Core data models for intronIC.

This package contains the core data structures for representing genomic features.
"""

from .models import (
    GenomeFeature,
    Parent,
    Gene,
    Transcript,
    Exon,
)

from .intron import (
    Intron,
    IntronScores,
    IntronSequences,
    IntronMetadata,
)

__all__ = [
    # Feature hierarchy
    'GenomeFeature',
    'Parent',
    'Gene',
    'Transcript',
    'Exon',
    # Intron composition
    'Intron',
    'IntronScores',
    'IntronSequences',
    'IntronMetadata',
]
