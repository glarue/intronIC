"""
Extraction pipeline for intronIC.

This package contains modules for extracting introns from genomic data:
- annotator: Build gene/transcript/exon hierarchies from annotations
- intronator: Generate introns from exon pairs
- sequences: Extract intron sequences from genome
- filters: Filter and deduplicate introns
"""

from intronIC.extraction.annotator import AnnotationHierarchyBuilder
from intronIC.extraction.intronator import IntronGenerator
from intronIC.extraction.sequences import SequenceExtractor
from intronIC.extraction.filters import IntronFilter

__all__ = [
    'AnnotationHierarchyBuilder',
    'IntronGenerator',
    'SequenceExtractor',
    'IntronFilter',
]
