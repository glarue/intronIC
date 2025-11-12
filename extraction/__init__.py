"""
Extraction pipeline for intronIC.

This package contains modules for extracting introns from genomic data:
- annotator: Build gene/transcript/exon hierarchies from annotations
- intronator: Generate introns from exon pairs
- sequences: Extract intron sequences from genome
- filters: Filter and deduplicate introns
"""

from extraction.annotator import AnnotationHierarchyBuilder
from extraction.intronator import IntronGenerator
from extraction.sequences import SequenceExtractor
from extraction.filters import IntronFilter

__all__ = [
    'AnnotationHierarchyBuilder',
    'IntronGenerator',
    'SequenceExtractor',
    'IntronFilter',
]
