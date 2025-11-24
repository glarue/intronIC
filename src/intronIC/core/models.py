"""
Core data models for genomic features.

This module defines the hierarchical data model for representing genomic annotations:
- GenomeFeature: Base class for all genomic features
- Parent: Abstract parent class for hierarchical features (Gene, Transcript)
- Gene: Represents a gene with multiple transcripts
- Transcript: Represents a transcript with multiple exons
- Exon: Represents an exon (coding or non-coding)

Design principles:
- Immutable where possible (using frozen dataclasses)
- Comprehensive type hints
- Validation at creation time
- 1-based coordinate system (intronIC internal format)

Author: intronIC refactoring project
Date: 2025-11-02
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Optional, List, Set, Dict, Any
from intronIC.utils.coordinates import GenomicCoordinate, Strand


@dataclass(frozen=True, slots=True)
class GenomeFeature:
    """
    Base class for all genomic features.

    Represents a feature with genomic coordinates and basic attributes.
    All coordinates are in 1-based, closed interval format (intronIC internal).

    Attributes:
        feature_id: Unique identifier for this feature
        coordinates: GenomicCoordinate object (chromosome, start, stop, strand)
        feature_type: Type of feature (gene, transcript, exon, etc.)
        attributes: Additional GFF3/GTF attributes (e.g., Name, gene_id, etc.)

    Examples:
        >>> coord = GenomicCoordinate("chr1", 1000, 2000, '+', '1-based')
        >>> feature = GenomeFeature("exon001", coord, "exon", {})
        >>> feature.chromosome
        'chr1'
        >>> feature.length
        1001
    """

    feature_id: str
    coordinates: GenomicCoordinate
    feature_type: str
    attributes: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate feature after initialization."""
        if not self.feature_id:
            raise ValueError("feature_id cannot be empty")

        if not self.feature_type:
            raise ValueError("feature_type cannot be empty")

        # Ensure coordinates are in 1-based system
        if self.coordinates.system != '1-based':
            raise ValueError(
                f"GenomeFeature requires 1-based coordinates, got {self.coordinates.system}"
            )

    # Convenience properties for accessing coordinate information

    @property
    def chromosome(self) -> str:
        """Chromosome/contig name."""
        return self.coordinates.chromosome

    @property
    def start(self) -> int:
        """Start coordinate (1-based, inclusive)."""
        return self.coordinates.start

    @property
    def stop(self) -> int:
        """Stop coordinate (1-based, inclusive)."""
        return self.coordinates.stop

    @property
    def strand(self) -> Strand:
        """Strand (+ or -)."""
        return self.coordinates.strand

    @property
    def length(self) -> int:
        """Feature length in bases."""
        return self.coordinates.length

    def __str__(self) -> str:
        """Human-readable string representation."""
        return (
            f"{self.feature_type}:{self.feature_id} "
            f"{self.coordinates}"
        )

    def __repr__(self) -> str:
        """Developer-friendly representation."""
        return (
            f"GenomeFeature("
            f"feature_id='{self.feature_id}', "
            f"coordinates={self.coordinates!r}, "
            f"feature_type='{self.feature_type}', "
            f"attributes={self.attributes!r})"
        )


@dataclass(slots=True)
class Parent(ABC):
    """
    Abstract base class for hierarchical parent features (Gene, Transcript).

    Parent features can have children (e.g., Gene has Transcripts, Transcript has Exons).
    Unlike leaf features (GenomeFeature), Parents are mutable to allow adding children.

    Attributes:
        feature_id: Unique identifier
        coordinates: GenomicCoordinate
        attributes: Additional attributes
        children: Set of child feature IDs

    Note:
        Not frozen (mutable) to allow adding children after creation.
        Use add_child() method to maintain consistency.
        Each subclass should define its own feature_type property.
    """

    feature_id: str
    coordinates: GenomicCoordinate
    attributes: Dict[str, Any] = field(default_factory=dict)
    children: Set[str] = field(default_factory=set)

    def __post_init__(self) -> None:
        """Validate parent feature after initialization."""
        if not self.feature_id:
            raise ValueError("feature_id cannot be empty")

        # Ensure coordinates are in 1-based system
        if self.coordinates.system != '1-based':
            raise ValueError(
                f"Parent feature requires 1-based coordinates, got {self.coordinates.system}"
            )

    # Convenience properties (same as GenomeFeature)

    @property
    def chromosome(self) -> str:
        """Chromosome/contig name."""
        return self.coordinates.chromosome

    @property
    def start(self) -> int:
        """Start coordinate (1-based, inclusive)."""
        return self.coordinates.start

    @property
    def stop(self) -> int:
        """Stop coordinate (1-based, inclusive)."""
        return self.coordinates.stop

    @property
    def strand(self) -> Strand:
        """Strand (+ or -)."""
        return self.coordinates.strand

    @property
    def length(self) -> int:
        """Feature length in bases."""
        return self.coordinates.length

    # Child management

    def add_child(self, child_id: str) -> None:
        """
        Add a child feature ID to this parent.

        Args:
            child_id: ID of the child feature

        Examples:
            >>> coord = GenomicCoordinate("chr1", 1000, 5000, '+', '1-based')
            >>> gene = Gene("GENE001", coord, {})
            >>> gene.add_child("TRANS001")
            >>> "TRANS001" in gene.children
            True
        """
        self.children.add(child_id)

    def remove_child(self, child_id: str) -> None:
        """
        Remove a child feature ID from this parent.

        Args:
            child_id: ID of the child feature to remove
        """
        self.children.discard(child_id)

    @property
    def num_children(self) -> int:
        """Number of child features."""
        return len(self.children)

    @abstractmethod
    def __str__(self) -> str:
        """String representation (must be implemented by subclasses)."""
        pass


@dataclass(slots=True)
class Gene(Parent):
    """
    Represents a gene with one or more transcripts.

    A gene is a collection of transcripts (isoforms) at a genomic locus.

    Attributes:
        feature_id: Gene ID (e.g., ENSG00000012345)
        coordinates: Gene coordinates (union of all transcript coordinates)
        attributes: Additional attributes (e.g., gene_name, gene_biotype)
        children: Set of transcript IDs belonging to this gene

    Examples:
        >>> coord = GenomicCoordinate("chr1", 1000, 5000, '+', '1-based')
        >>> gene = Gene("ENSG001", coord, {"gene_name": "BRCA1"})
        >>> gene.add_child("ENST001")
        >>> gene.num_children
        1
    """

    @property
    def feature_type(self) -> str:
        """Feature type (always 'gene')."""
        return 'gene'

    @property
    def gene_name(self) -> Optional[str]:
        """Gene name (e.g., BRCA1), if available."""
        return self.attributes.get('gene_name') or self.attributes.get('Name')

    @property
    def transcript_ids(self) -> Set[str]:
        """Set of transcript IDs (alias for children)."""
        return self.children

    def __str__(self) -> str:
        """Human-readable string representation."""
        name = self.gene_name or self.feature_id
        return (
            f"Gene:{name} "
            f"{self.chromosome}:{self.start}-{self.stop}({self.strand}) "
            f"[{self.num_children} transcripts]"
        )

    def __repr__(self) -> str:
        """Developer-friendly representation."""
        return (
            f"Gene("
            f"feature_id='{self.feature_id}', "
            f"coordinates={self.coordinates!r}, "
            f"num_children={self.num_children})"
        )


@dataclass(slots=True)
class Transcript(Parent):
    """
    Represents a transcript (mRNA isoform) with one or more exons.

    A transcript is a specific isoform of a gene, containing exons.

    Attributes:
        feature_id: Transcript ID (e.g., ENST00000012345)
        coordinates: Transcript coordinates (union of all exon coordinates)
        attributes: Additional attributes (e.g., transcript_name, transcript_biotype)
        children: Set of exon IDs belonging to this transcript
        parent_id: Gene ID this transcript belongs to

    Examples:
        >>> coord = GenomicCoordinate("chr1", 1000, 5000, '+', '1-based')
        >>> trans = Transcript("ENST001", coord, parent_id="ENSG001")
        >>> trans.add_child("EXON001")
        >>> trans.parent_id
        'ENSG001'
    """

    parent_id: Optional[str] = None

    @property
    def feature_type(self) -> str:
        """Feature type (always 'transcript')."""
        return 'transcript'

    @property
    def transcript_name(self) -> Optional[str]:
        """Transcript name, if available."""
        return self.attributes.get('transcript_name') or self.attributes.get('Name')

    @property
    def exon_ids(self) -> Set[str]:
        """Set of exon IDs (alias for children)."""
        return self.children

    @property
    def gene_id(self) -> Optional[str]:
        """Parent gene ID (alias for parent_id)."""
        return self.parent_id

    def __str__(self) -> str:
        """Human-readable string representation."""
        name = self.transcript_name or self.feature_id
        parent_info = f" (gene:{self.parent_id})" if self.parent_id else ""
        return (
            f"Transcript:{name} "
            f"{self.chromosome}:{self.start}-{self.stop}({self.strand}) "
            f"[{self.num_children} exons]{parent_info}"
        )

    def __repr__(self) -> str:
        """Developer-friendly representation."""
        return (
            f"Transcript("
            f"feature_id='{self.feature_id}', "
            f"coordinates={self.coordinates!r}, "
            f"parent_id='{self.parent_id}', "
            f"num_children={self.num_children})"
        )


@dataclass(frozen=True, slots=True)
class Exon:
    """
    Represents an exon (coding or non-coding region).

    Exons are the building blocks of transcripts. Introns are derived from
    the gaps between consecutive exons.

    Attributes:
        feature_id: Exon ID
        coordinates: Exon coordinates
        parent_id: Transcript ID this exon belongs to
        attributes: Additional attributes (e.g., exon_number)
        phase: Coding phase (0, 1, 2, or None for non-coding)
        is_coding: Whether this exon is coding (CDS) or non-coding (UTR)

    Examples:
        >>> coord = GenomicCoordinate("chr1", 1000, 1200, '+', '1-based')
        >>> exon = Exon("EXON001", coord, "ENST001", {}, 0, True)
        >>> exon.length
        201
        >>> exon.is_coding
        True
    """

    feature_id: str
    coordinates: GenomicCoordinate
    parent_id: Optional[str] = None
    attributes: Dict[str, Any] = field(default_factory=dict)
    phase: Optional[int] = None
    is_coding: bool = False

    def __post_init__(self) -> None:
        """Validate exon after initialization."""
        if not self.feature_id:
            raise ValueError("feature_id cannot be empty")

        # Ensure coordinates are in 1-based system
        if self.coordinates.system != '1-based':
            raise ValueError(
                f"Exon requires 1-based coordinates, got {self.coordinates.system}"
            )

        # Validate phase
        if self.phase is not None and self.phase not in [0, 1, 2]:
            raise ValueError(f"phase must be 0, 1, 2, or None, got {self.phase}")

    # Convenience properties

    @property
    def chromosome(self) -> str:
        """Chromosome/contig name."""
        return self.coordinates.chromosome

    @property
    def start(self) -> int:
        """Start coordinate (1-based, inclusive)."""
        return self.coordinates.start

    @property
    def stop(self) -> int:
        """Stop coordinate (1-based, inclusive)."""
        return self.coordinates.stop

    @property
    def strand(self) -> Strand:
        """Strand (+ or -)."""
        return self.coordinates.strand

    @property
    def length(self) -> int:
        """Exon length in bases."""
        return self.coordinates.length

    @property
    def transcript_id(self) -> Optional[str]:
        """Parent transcript ID (alias for parent_id)."""
        return self.parent_id

    @property
    def exon_number(self) -> Optional[int]:
        """Exon number within transcript, if available."""
        num = self.attributes.get('exon_number')
        return int(num) if num is not None else None

    def __str__(self) -> str:
        """Human-readable string representation."""
        coding_info = " [CDS]" if self.is_coding else " [UTR]"
        parent_info = f" (transcript:{self.parent_id})" if self.parent_id else ""
        return (
            f"Exon:{self.feature_id} "
            f"{self.coordinates}"
            f"{coding_info}{parent_info}"
        )

    def __repr__(self) -> str:
        """Developer-friendly representation."""
        return (
            f"Exon("
            f"feature_id='{self.feature_id}', "
            f"coordinates={self.coordinates!r}, "
            f"parent_id='{self.parent_id}', "
            f"phase={self.phase}, "
            f"is_coding={self.is_coding})"
        )


if __name__ == "__main__":
    import doctest
    doctest.testmod()
