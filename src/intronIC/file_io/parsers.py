"""
Parsers for various genomic file formats.

This module provides parsers for:
- GFF3/GTF annotations (modular, swappable backend)
- BED files
- .iic sequence files

Design principle: Use Protocol/ABC to define interfaces so parser
implementations can be swapped without changing downstream code.

Author: intronIC refactoring project
Date: 2025-11-02
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Protocol, Iterator, List, Optional, Union, Tuple
from abc import ABC, abstractmethod
from smart_open import open as smart_open


# ============================================================================
# Annotation Parsing (GFF3/GTF) - Modular interface
# ============================================================================

@dataclass(slots=True)
class AnnotationLine:
    """
    Standardized representation of a parsed annotation line.

    This is the common interface that all annotation parsers must return.
    Decouples downstream code from specific parser implementations.

    Attributes:
        name: Feature ID
        feat_type: Feature type (gene, transcript, exon, cds, etc.)
        parent: List of parent IDs (can be multiple in GFF3)
        grandparent: Grandparent ID (gene for exons/CDS)
        region: Chromosome/contig name
        strand: Strand ('+' or '-')
        start: Start coordinate (format depends on parser, caller must convert)
        stop: Stop coordinate (format depends on parser, caller must convert)
        line_number: Original line number in file
        phase: Coding phase (0, 1, 2) or None
        attributes: Dictionary of additional attributes

    Examples:
        >>> line = AnnotationLine(
        ...     name="EXON001",
        ...     feat_type="exon",
        ...     parent=["TRANS001"],
        ...     grandparent="GENE001",
        ...     region="chr1",
        ...     strand="+",
        ...     start=1000,
        ...     stop=1200,
        ...     line_number=42,
        ...     phase=0
        ... )
        >>> line.feat_type
        'exon'
    """

    name: str
    feat_type: str
    parent: List[str]
    grandparent: Optional[str]
    region: str
    strand: str
    start: int
    stop: int
    line_number: int
    phase: Optional[int] = None
    attributes: dict = field(default_factory=dict)


class AnnotationParser(Protocol):
    """
    Protocol defining the interface for annotation parsers.

    Any annotation parser (biogl, GTFtools, custom, etc.) must implement
    this interface to be usable by intronIC.

    This allows us to swap parser implementations without changing
    downstream code.
    """

    def parse_line(self, line: str, line_number: int) -> Optional[AnnotationLine]:
        """
        Parse a single annotation line.

        Args:
            line: Raw line from annotation file
            line_number: Line number (for error reporting)

        Returns:
            AnnotationLine object, or None if line should be skipped
            (e.g., comment, header, invalid)

        Note:
            Implementations should handle:
            - Comment lines (starting with #)
            - Header lines
            - Malformed lines (return None)
            - Multiple parents (GFF3 feature)
        """
        ...

    def parse_file(self, file_path: Union[str, Path]) -> Iterator[AnnotationLine]:
        """
        Parse an entire annotation file.

        Args:
            file_path: Path to GFF3/GTF file

        Yields:
            AnnotationLine objects

        Note:
            Should handle gzipped files automatically.
        """
        ...


# ============================================================================
# BioGL-based implementation (current default)
# ============================================================================

class BioGLAnnotationParser:
    """
    Annotation parser using biogl.GxfParse backend.

    This is a thin adapter that wraps biogl's GxfParse to provide
    the standardized AnnotationLine interface.

    Advantages of biogl:
    - Already a dependency
    - Permissive parsing (handles edge cases)
    - Extracts all needed metadata
    - Battle-tested in original intronIC

    Args:
        clean_names: If True, remove 'transcript:' and 'gene:' prefixes from IDs

    Examples:
        >>> parser = BioGLAnnotationParser()
        >>> # line = parser.parse_line("chr1\\tENSEMBL\\texon\\t1000\\t2000\\t.\\t+\\t.\\tID=exon1;Parent=trans1", 1)
        >>> # line.feat_type
        >>> # 'exon'
    """

    def __init__(self, clean_names: bool = True):
        """Initialize parser.

        Args:
            clean_names: If True, remove Ensembl-style 'gene:' and 'transcript:' prefixes.
                         This only removes colon-separated prefixes (Ensembl convention).
                         RefSeq-style hyphen/underscore prefixes (gene-, rna_) are preserved
                         to avoid ID collisions.
        """
        from biogl import GxfParse
        self._gxf_parse = GxfParse
        self._clean_names = clean_names

    @staticmethod
    def _clean_id(id_str: Optional[str]) -> Optional[str]:
        """
        Remove Ensembl-style GFF3 ID prefixes (only colon-separated).

        This conservatively removes only "gene:" and "transcript:" prefixes
        that use colon separators (Ensembl style), which are truly redundant.

        Does NOT remove hyphen/underscore prefixes (e.g., "gene-", "rna_")
        which are used in RefSeq/NCBI annotations for disambiguation.

        Args:
            id_str: Feature ID string

        Returns:
            Cleaned ID string

        Examples:
            >>> _clean_id("transcript:ENST00000123456")
            'ENST00000123456'
            >>> _clean_id("gene:ENSG00000139618")
            'ENSG00000139618'
            >>> _clean_id("gene-MIR6859-1")
            'gene-MIR6859-1'
            >>> _clean_id("rna-NM_001234")
            'rna-NM_001234'
        """
        if not id_str:
            return id_str

        # Only remove Ensembl-style "gene:" and "transcript:" prefixes
        # These are redundant because Ensembl IDs already indicate type
        # (ENSG* = gene, ENST* = transcript, ENSP* = protein)
        if id_str.startswith('gene:'):
            return id_str[5:]
        elif id_str.startswith('transcript:'):
            return id_str[11:]

        return id_str

    def parse_line(self, line: str, line_number: int) -> Optional[AnnotationLine]:
        """
        Parse a single annotation line using biogl.

        Args:
            line: Raw GFF3/GTF line
            line_number: Line number

        Returns:
            AnnotationLine or None if line is invalid/comment
        """
        try:
            line_info = self._gxf_parse(line, line_number)
        except TypeError:
            # biogl raises TypeError for non-annotation lines (comments, headers, etc.)
            return None

        # Convert biogl output to AnnotationLine
        # Note: biogl may return multiple parents in some cases
        parent_list = line_info.parent if isinstance(line_info.parent, list) else [line_info.parent]

        # Handle phase - biogl may return '.' or None
        phase = None
        if hasattr(line_info, 'phase') and line_info.phase is not None:
            if line_info.phase != '.':
                try:
                    phase = int(line_info.phase)
                except (ValueError, TypeError):
                    phase = None

        # Generate unique ID for features without IDs (common for exons in Ensembl GFF3)
        # biogl returns 'exon_None' string (not Python None) for exons without IDs
        # Use location-based identifier: feattype_parent_chrom_start_stop_strand
        feature_name = line_info.name
        if feature_name is None or feature_name == 'exon_None' or 'None' in str(feature_name):
            parent_str = parent_list[0] if parent_list else "none"
            feature_name = (
                f"{line_info.feat_type}_{parent_str}_"
                f"{line_info.region}_{line_info.start}_{line_info.stop}_{line_info.strand}"
            )

        # Clean ID prefixes if requested
        if self._clean_names:
            feature_name = self._clean_id(feature_name)
            parent_list = [self._clean_id(p) for p in parent_list if p]
            grandparent = self._clean_id(line_info.grandparent)
        else:
            grandparent = line_info.grandparent

        return AnnotationLine(
            name=feature_name,
            feat_type=line_info.feat_type.lower(),  # Normalize to lowercase
            parent=parent_list,
            grandparent=grandparent,
            region=line_info.region,
            strand=line_info.strand,
            start=line_info.start,
            stop=line_info.stop,
            line_number=line_info.line_number,
            phase=phase,
            attributes={}  # biogl doesn't expose raw attributes dict, could extend if needed
        )

    def parse_file(self, file_path: Union[str, Path]) -> Iterator[AnnotationLine]:
        """
        Parse entire annotation file.

        Args:
            file_path: Path to GFF3/GTF file

        Yields:
            AnnotationLine objects
        """
        from biogl import flex_open

        file_path = Path(file_path)

        with flex_open(str(file_path)) as f:
            for line_num, line in enumerate(f, start=1):
                parsed = self.parse_line(line, line_num)
                if parsed is not None:
                    yield parsed

    def parse_lines(self, lines: List[str]) -> Iterator[AnnotationLine]:
        """
        Parse a list of annotation lines (streaming-friendly).

        This method enables memory-efficient streaming processing by accepting
        pre-extracted lines instead of requiring full file access.

        Args:
            lines: List of annotation lines (GFF3/GTF format)

        Yields:
            AnnotationLine objects

        Examples:
            >>> parser = BioGLAnnotationParser()
            >>> lines = ["chr1\\tENSEMBL\\texon\\t1000\\t2000\\t.\\t+\\t.\\tID=exon1;Parent=trans1"]
            >>> for ann in parser.parse_lines(lines):
            ...     print(ann.feat_type)
            exon

        Note:
            Line numbers in the returned AnnotationLine objects will be based on
            the position within the provided list (1-indexed), not the original
            file line numbers.
        """
        for line_num, line in enumerate(lines, start=1):
            parsed = self.parse_line(line, line_num)
            if parsed is not None:
                yield parsed


# ============================================================================
# BED Format Parsing
# ============================================================================

@dataclass
class BEDLine:
    """
    Standardized representation of a BED line.

    BED format is simple:
    chrom  start  stop  name  score  strand  [additional columns...]

    Coordinates are 0-based, half-open: [start, stop)

    Attributes:
        chrom: Chromosome name
        start: Start position (0-based)
        stop: Stop position (0-based, exclusive)
        name: Feature name
        score: Feature score (or '.')
        strand: Strand ('+', '-', or '.')
        extra_fields: Any additional columns beyond standard 6

    Examples:
        >>> bed = BEDLine("chr1", 1000, 2000, "feature1", "100", "+")
        >>> bed.chrom
        'chr1'
        >>> bed.start
        1000
    """

    chrom: str
    start: int
    stop: int
    name: str = "."
    score: str = "."
    strand: str = "."
    extra_fields: List[str] = field(default_factory=list)


class BEDParser:
    """
    Parser for BED format files.

    BED format specification:
    - Tab-delimited
    - Minimum 3 columns: chrom, start (0-based), stop (0-based, exclusive)
    - Optional columns: name, score, strand, thickStart, thickStop, itemRgb, etc.
    - Comments start with '#' or 'track' or 'browser'

    Examples:
        >>> parser = BEDParser()
        >>> # bed = parser.parse_line("chr1\\t1000\\t2000\\tfeature1\\t100\\t+")
        >>> # bed.start
        >>> # 1000
    """

    def parse_line(self, line: str) -> Optional[BEDLine]:
        """
        Parse a single BED line.

        Args:
            line: Raw BED line

        Returns:
            BEDLine object or None if comment/invalid
        """
        line = line.strip()

        # Skip empty lines
        if not line:
            return None

        # Skip comment lines
        if line.startswith('#') or line.startswith('track') or line.startswith('browser'):
            return None

        fields = line.split('\t')

        # Need at least 3 fields (chrom, start, stop)
        if len(fields) < 3:
            return None

        try:
            chrom = fields[0]
            start = int(fields[1])
            stop = int(fields[2])
            name = fields[3] if len(fields) > 3 else "."
            score = fields[4] if len(fields) > 4 else "."
            strand = fields[5] if len(fields) > 5 else "."
            extra = fields[6:] if len(fields) > 6 else []

            return BEDLine(
                chrom=chrom,
                start=start,
                stop=stop,
                name=name,
                score=score,
                strand=strand,
                extra_fields=extra
            )
        except (ValueError, IndexError):
            return None

    def parse_file(self, file_path: Union[str, Path]) -> Iterator[BEDLine]:
        """
        Parse entire BED file.

        Args:
            file_path: Path to BED file

        Yields:
            BEDLine objects
        """
        file_path = Path(file_path)

        with smart_open(file_path, 'rt') as f:
            for line in f:
                parsed = self.parse_line(line)
                if parsed is not None:
                    yield parsed


# ============================================================================
# .iic Sequence File Parsing
# ============================================================================

@dataclass
class SequenceLine:
    """
    Representation of a line from .iic sequence file.

    intronIC .iic format (tab-delimited):
    name  upstream_flank  sequence  downstream_flank  [optional: score]

    Attributes:
        name: Intron/feature name
        upstream_flank: 5' flanking exonic sequence
        sequence: Intron sequence
        downstream_flank: 3' flanking exonic sequence
        score: Optional SVM score

    Examples:
        >>> seq = SequenceLine("intron1", "AGGCT", "GTAAGT...TTTAG", "CATGG")
        >>> seq.name
        'intron1'
    """

    name: str
    upstream_flank: str
    sequence: str
    downstream_flank: str
    score: Optional[float] = None


class SequenceParser:
    """
    Parser for intronIC .iic sequence files.

    Format:
    - Tab-delimited
    - 4-5 columns: name, [score], upstream_flank, sequence, downstream_flank
    - No header
    - Score can be a number or "NA" (parsed as None)

    Examples:
        >>> parser = SequenceParser()
        >>> # seq = parser.parse_line("intron1\\t95.5\\tAGGCT\\tGTAAGT...TTTAG\\tCATGG")
        >>> # seq.score
        >>> # 95.5
    """

    def parse_line(self, line: str) -> Optional[SequenceLine]:
        """
        Parse a single sequence line.

        Args:
            line: Raw .iic line

        Returns:
            SequenceLine object or None if invalid
        """
        line = line.strip()

        if not line:
            return None

        fields = line.split('\t')

        # Need at least 4 fields (name + 3 sequence fields)
        # Can have 5 fields if score is included
        if len(fields) < 4:
            return None

        try:
            name = fields[0]

            # Determine if we have score or not
            if len(fields) == 5:
                # Format: name, score, upstream, sequence, downstream
                score_str = fields[1]
                score = None if score_str == "NA" else float(score_str)
                upstream_flank = fields[2]
                sequence = fields[3]
                downstream_flank = fields[4]
            else:
                # Format: name, upstream, sequence, downstream (no score)
                score = None
                upstream_flank = fields[1]
                sequence = fields[2]
                downstream_flank = fields[3]

            return SequenceLine(
                name=name,
                upstream_flank=upstream_flank,
                sequence=sequence,
                downstream_flank=downstream_flank,
                score=score
            )
        except (ValueError, IndexError):
            return None

    def parse_file(self, file_path: Union[str, Path]) -> Iterator[SequenceLine]:
        """
        Parse entire .iic sequence file.

        Args:
            file_path: Path to .iic file

        Yields:
            SequenceLine objects
        """
        file_path = Path(file_path)

        with smart_open(file_path, 'rt') as f:
            for line in f:
                parsed = self.parse_line(line)
                if parsed is not None:
                    yield parsed


# ============================================================================
# Chromosome extraction utilities for per-chromosome processing
# ============================================================================

def get_chromosomes_from_annotation(annotation_file: Union[str, Path]) -> List[str]:
    """
    Extract list of unique chromosomes from annotation file.

    This is used for per-chromosome processing to reduce memory usage.

    Args:
        annotation_file: Path to GFF3/GTF annotation file

    Returns:
        Sorted list of unique chromosome names

    Examples:
        >>> chroms = get_chromosomes_from_annotation('genome.gff3')
        >>> print(chroms)
        ['chr1', 'chr2', 'chr3', ...]

    Note:
        This does a single fast pass through the annotation file,
        extracting only the chromosome column.
    """
    from biogl import flex_open

    chromosomes = set()
    with flex_open(str(annotation_file)) as f:
        for line in f:
            # Skip comments and headers
            if line.startswith('#'):
                continue

            # Parse chromosome from first column
            fields = line.strip().split('\t')
            if len(fields) >= 9:  # Valid GFF3/GTF line
                chromosomes.add(fields[0])

    return sorted(chromosomes)


def filter_annotation_lines_by_chromosome(
    annotation_file: Union[str, Path],
    chromosome: str
) -> Iterator[AnnotationLine]:
    """
    Parse annotation file but only return lines for a specific contig.

    This enables per-contig processing without loading the entire
    annotation into memory at once.

    Args:
        annotation_file: Path to GFF3/GTF annotation file
        chromosome: Contig/chromosome name to filter for

    Yields:
        AnnotationLine objects for the specified contig

    Examples:
        >>> parser = BioGLAnnotationParser()
        >>> chr1_lines = filter_annotation_lines_by_chromosome('genome.gff3', 'chr1')
        >>> chr1_genes = builder.build_from_annotations(chr1_lines)

    Note:
        This is memory-efficient because it yields lines one at a time
        instead of loading all lines into memory.
    """
    parser = BioGLAnnotationParser()

    for ann_line in parser.parse_file(annotation_file):
        if ann_line.region == chromosome:
            yield ann_line


if __name__ == "__main__":
    import doctest
    doctest.testmod()
