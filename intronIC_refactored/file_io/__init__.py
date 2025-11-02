"""
I/O modules for intronIC.

This package contains modules for reading and writing various file formats.
"""

from .genome import GenomeReader, parse_fasta
from .parsers import (
    AnnotationLine,
    AnnotationParser,
    BioGLAnnotationParser,
    BEDLine,
    BEDParser,
    SequenceLine,
    SequenceParser,
)
from .writers import (
    BEDWriter,
    MetaWriter,
    SequenceWriter,
    ScoreWriter,
    MappingWriter,
)

__all__ = [
    # Genome reading
    'GenomeReader',
    'parse_fasta',
    # Annotation parsing
    'AnnotationLine',
    'AnnotationParser',
    'BioGLAnnotationParser',
    # BED parsing
    'BEDLine',
    'BEDParser',
    # Sequence file parsing
    'SequenceLine',
    'SequenceParser',
    # Writers
    'BEDWriter',
    'MetaWriter',
    'SequenceWriter',
    'ScoreWriter',
    'MappingWriter',
]
