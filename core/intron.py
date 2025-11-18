"""
Intron data models using composition pattern.

This module defines the core Intron class and its composed components:
- IntronScores: All scoring-related attributes
- IntronSequences: All sequence-related attributes
- IntronMetadata: All metadata attributes (parent, tags, etc.)
- Intron: Main class composing all components

Design rationale:
- Composition over inheritance reduces coupling
- Clear separation of concerns (scores vs sequences vs metadata)
- Each component can be None if data not yet available
- Easier to test individual components
- More maintainable than 41-attribute monolithic class

Author: intronIC refactoring project
Date: 2025-11-02
"""

from dataclasses import dataclass, field
from typing import Literal, Optional

from utils.coordinates import GenomicCoordinate

# Type aliases for clarity
IntronType = Literal["u2", "u12", "unknown"]
OmissionCode = Optional[
    Literal["s", "a", "n", "i", "v", "d"]
]  # short, ambiguous, noncanonical, isoform, overlap, duplicate


@dataclass(frozen=True, slots=True)
class IntronScores:
    """
    Scoring data for an intron.

    Contains all PWM scores, z-scores, and classification results.

    Attributes:
        five_raw_score: Raw PWM score for 5' splice site
        bp_raw_score: Raw PWM score for branch point
        three_raw_score: Raw PWM score for 3' splice site
        five_z_score: Normalized z-score for 5' splice site
        bp_z_score: Normalized z-score for branch point
        three_z_score: Normalized z-score for 3' splice site
        min_5_bp: BothEndsStrong feature: min(5', BP) = "both must be strong"
        min_5_3: BothEndsStrong feature: min(5', 3') = "both must be strong"
        max_5_bp: BothEndsStrong feature: max(5', BP) = "at least one strong" (optional)
        max_5_3: BothEndsStrong feature: max(5', 3') = "at least one strong" (optional)
        svm_score: SVM probability score (0-100%)
        relative_score: Score relative to threshold
        decision_distance: Distance from SVM decision boundary

    Examples:
        >>> scores = IntronScores(
        ...     five_raw_score=12.5, bp_raw_score=8.3, three_raw_score=15.2,
        ...     five_z_score=2.1, bp_z_score=1.8, three_z_score=2.5,
        ...     svm_score=95.5, relative_score=5.5, decision_distance=1.2
        ... )
        >>> scores.svm_score
        95.5
        >>> scores.is_high_confidence(threshold=90.0)
        True
    """

    five_raw_score: Optional[float] = None
    bp_raw_score: Optional[float] = None
    three_raw_score: Optional[float] = None
    five_z_score: Optional[float] = None
    bp_z_score: Optional[float] = None
    three_z_score: Optional[float] = None
    min_5_bp: Optional[float] = None
    min_5_3: Optional[float] = None
    max_5_bp: Optional[float] = None
    max_5_3: Optional[float] = None
    svm_score: Optional[float] = None
    relative_score: Optional[float] = None
    decision_distance: Optional[float] = None

    def is_high_confidence(self, threshold: float = 90.0) -> bool:
        """
        Check if this intron has high-confidence U12 classification.

        Args:
            threshold: SVM score threshold (default: 90.0)

        Returns:
            True if svm_score >= threshold, False otherwise
        """
        if self.svm_score is None:
            return False
        return self.svm_score >= threshold

    def has_all_scores(self) -> bool:
        """Check if all scoring fields are populated."""
        return all(
            [
                self.five_raw_score is not None,
                self.bp_raw_score is not None,
                self.three_raw_score is not None,
                self.five_z_score is not None,
                self.bp_z_score is not None,
                self.three_z_score is not None,
                self.svm_score is not None,
            ]
        )

    def __str__(self) -> str:
        """Human-readable string representation."""
        if self.svm_score is not None:
            return f"IntronScores(svm={self.svm_score:.1f}%, relative={self.relative_score:+.1f})"
        return "IntronScores(not scored)"


@dataclass(frozen=True, slots=True)
class IntronSequences:
    """
    Sequence data for an intron.

    Contains the intron sequence and all extracted sub-sequences.

    Attributes:
        seq: Full intron sequence
        five_seq: 5' splice site sequence
        three_seq: 3' splice site sequence
        bp_seq: Branch point sequence (best match)
        bp_region_seq: Full branch point search region
        upstream_flank: Exonic sequence upstream of intron
        downstream_flank: Exonic sequence downstream of intron
        five_prime_dnt: Terminal 5' dinucleotide (e.g., 'GT')
        three_prime_dnt: Terminal 3' dinucleotide (e.g., 'AG')
        five_display_seq: First 10bp of intron for display (motif schematic)
        three_display_seq: From BP end to intron end for display
        bp_seq_u2: U2-type branch point sequence
        bp_relative_coords: (start, stop) position of BP within bp_region_seq

    Examples:
        >>> seqs = IntronSequences(
        ...     seq="GTAAGT...TTTAG",
        ...     five_seq="GTAAGT",
        ...     three_seq="TTTAG",
        ...     bp_seq="TACTAAC",
        ...     upstream_flank="AGGCT",
        ...     downstream_flank="CATGG"
        ... )
        >>> seqs.five_seq
        'GTAAGT'
        >>> seqs.has_sequences()
        True
    """

    seq: Optional[str] = None
    five_seq: Optional[str] = None
    three_seq: Optional[str] = None
    bp_seq: Optional[str] = None
    bp_region_seq: Optional[str] = None
    upstream_flank: Optional[str] = None
    downstream_flank: Optional[str] = None
    five_prime_dnt: Optional[str] = None
    three_prime_dnt: Optional[str] = None
    # Display sequences for motif schematic generation
    five_display_seq: Optional[str] = None
    three_display_seq: Optional[str] = None
    # U2 branch point information
    bp_seq_u2: Optional[str] = None
    # Branch point coordinates within bp_region_seq
    bp_relative_coords: Optional[tuple[int, int]] = None

    def has_sequences(self) -> bool:
        """Check if core sequences are populated."""
        return all(
            [
                self.seq is not None,
                self.five_seq is not None,
                self.three_seq is not None,
            ]
        )

    def has_flanks(self) -> bool:
        """Check if flanking sequences are populated."""
        return all(
            [
                self.upstream_flank is not None,
                self.downstream_flank is not None,
            ]
        )

    @property
    def terminal_dinucleotides(self) -> Optional[str]:
        """
        Extract terminal dinucleotides (e.g., 'GT-AG', 'AT-AC').

        Uses stored five_prime_dnt and three_prime_dnt fields if available
        (memory-efficient), falls back to extracting from seq if needed.

        Returns:
            String like 'GT-AG' or None if sequence not available
        """
        # Prefer stored dnts (available after memory optimization clears seq)
        if self.five_prime_dnt and self.three_prime_dnt:
            return f"{self.five_prime_dnt}-{self.three_prime_dnt}"
        # Fall back to extracting from seq (for backwards compatibility)
        if self.seq is None or len(self.seq) < 4:
            return None
        return f"{self.seq[:2]}-{self.seq[-2:]}"

    def __str__(self) -> str:
        """Human-readable string representation."""
        if self.seq:
            dnts = self.terminal_dinucleotides or "??"
            return f"IntronSequences({dnts}, {len(self.seq)}bp)"
        return "IntronSequences(no sequence)"


@dataclass(slots=True)
class IntronMetadata:
    """
    Metadata for an intron.

    Contains parent relationships, tags, and classification info.
    Unlike Scores and Sequences, this is mutable to allow tagging
    during filtering and classification.

    Attributes:
        parent: Transcript ID this intron belongs to
        grandparent: Gene ID this intron belongs to
        index: Ordinal position in transcript (1-based)
        family_size: Total number of introns in this transcript
        parent_length: Length of parent transcript (for tiebreaking)
        line_number: Annotation line number (for final tiebreaking)
        type_id: Classification ('u2', 'u12', 'unknown')
        noncanonical: Whether intron has non-standard boundaries
        omitted: Omission reason code (None if not omitted)
        duplicate: Reference to representative if duplicate
        overlap: Overlapping coordinate record
        longest_isoform: Whether from longest transcript
        corrected: Whether boundaries were adjusted
        phase: Coding phase information
        dynamic_tags: Set of dynamic tags ([c:N], [d], [e], [n], [i], etc.)
        correction_distance: Distance boundaries were shifted (for [c:N] tag)

    Note:
        Mutable (not frozen) to allow updating tags during pipeline.
    """

    parent: Optional[str] = None
    grandparent: Optional[str] = None
    index: Optional[int] = None
    family_size: Optional[int] = None
    parent_length: Optional[int] = None
    line_number: Optional[int] = None  # Annotation line number for hierarchical sort tiebreaker
    type_id: IntronType = "unknown"
    noncanonical: bool = False
    omitted: OmissionCode = None
    duplicate: Optional[str] = None
    overlap: Optional[str] = None
    longest_isoform: bool = False
    corrected: bool = False
    phase: Optional[int] = None
    defined_by: Optional[str] = None  # 'cds' or 'exon' - which feature type defined this intron
    # Dynamic tagging system for various intron states
    dynamic_tags: set[str] = field(default_factory=set)
    correction_distance: Optional[int] = None

    def is_omitted(self) -> bool:
        """Check if this intron should be omitted."""
        return self.omitted is not None

    def is_duplicate(self) -> bool:
        """Check if this intron is a duplicate."""
        return self.duplicate is not None

    def is_canonical(self) -> bool:
        """Check if this intron has canonical boundaries."""
        return not self.noncanonical

    @property
    def fractional_position(self) -> Optional[float]:
        """
        Calculate fractional position in transcript (0.0-1.0).

        Returns:
            Fraction (0.0 = first intron, 1.0 = last), or None if index/family_size unknown
        """
        if self.index is None or self.family_size is None or self.family_size <= 1:
            return None
        return (self.index - 1) / (self.family_size - 1)

    def __str__(self) -> str:
        """Human-readable string representation."""
        parts = []
        if self.parent:
            parts.append(f"parent:{self.parent}")
        if self.index is not None and self.family_size is not None:
            parts.append(f"i{self.index}/{self.family_size}")
        if self.type_id != "unknown":
            parts.append(self.type_id)
        if self.omitted:
            parts.append(f"omitted:{self.omitted}")

        return f"IntronMetadata({', '.join(parts)})" if parts else "IntronMetadata()"


@dataclass(frozen=True, slots=True)
class Intron:
    """
    Main Intron class using composition pattern.

    Composes three main components:
    - coordinates: GenomicCoordinate (location)
    - scores: IntronScores (classification data)
    - sequences: IntronSequences (sequence data)
    - metadata: IntronMetadata (parent, tags, etc.)

    This design separates concerns and makes it easy to:
    - Create introns with partial data (e.g., coordinates only)
    - Update scores/sequences as they become available
    - Test components independently
    - Maintain clear boundaries between pipeline stages

    Attributes:
        intron_id: Unique identifier
        coordinates: Genomic location
        scores: Scoring data (optional)
        sequences: Sequence data (optional)
        metadata: Metadata (optional)

    Examples:
        >>> from utils.coordinates import GenomicCoordinate
        >>> coord = GenomicCoordinate("chr1", 1001, 2000, '+', '1-based')
        >>> intron = Intron("intron_1", coord)
        >>> intron.length
        1000
        >>> intron.chromosome
        'chr1'

        >>> # With sequences
        >>> seqs = IntronSequences(seq="GTAAGT...TTTAG", five_seq="GTAAGT", three_seq="TTTAG")
        >>> intron2 = Intron("intron_2", coord, sequences=seqs)
        >>> intron2.has_sequences
        True
    """

    intron_id: str
    coordinates: GenomicCoordinate
    scores: Optional[IntronScores] = None
    sequences: Optional[IntronSequences] = None
    metadata: Optional[IntronMetadata] = None

    def __post_init__(self) -> None:
        """Validate intron after initialization."""
        if not self.intron_id:
            raise ValueError("intron_id cannot be empty")

        # Ensure coordinates are in 1-based system
        if self.coordinates.system != "1-based":
            raise ValueError(
                f"Intron requires 1-based coordinates, got {self.coordinates.system}"
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
    def strand(self) -> str:
        """Strand (+ or -)."""
        return self.coordinates.strand

    @property
    def length(self) -> int:
        """Intron length in bases."""
        return self.coordinates.length

    # Status properties

    @property
    def has_scores(self) -> bool:
        """Check if scoring data is available."""
        return self.scores is not None and self.scores.has_all_scores()

    @property
    def has_sequences(self) -> bool:
        """Check if sequence data is available."""
        return self.sequences is not None and self.sequences.has_sequences()

    @property
    def has_metadata(self) -> bool:
        """Check if metadata is available."""
        return self.metadata is not None

    # Convenience properties for common attributes

    @property
    def svm_score(self) -> Optional[float]:
        """SVM score (0-100%), or None if not scored."""
        return self.scores.svm_score if self.scores else None

    @property
    def type_id(self) -> IntronType:
        """Intron type ('u2', 'u12', 'unknown')."""
        return self.metadata.type_id if self.metadata else "unknown"

    @property
    def terminal_dinucleotides(self) -> Optional[str]:
        """Terminal dinucleotides (e.g., 'GT-AG'), or None if sequence unavailable."""
        return self.sequences.terminal_dinucleotides if self.sequences else None

    # Factory methods

    @classmethod
    def from_exon_pair(
        cls,
        exon1,  # Type would be Exon but avoiding circular import
        exon2,  # Type would be Exon
        intron_id: Optional[str] = None,
    ) -> "Intron":
        """
        Create an Intron from a pair of adjacent exons.

        Args:
            exon1: First exon (upstream in genomic coordinates)
            exon2: Second exon (downstream in genomic coordinates)
            intron_id: Optional ID (will be auto-generated if not provided)

        Returns:
            Intron object with coordinates derived from exon gap

        Raises:
            ValueError: If exons are not adjacent or on different strands/chromosomes

        Examples:
            >>> from core.models import Exon
            >>> from utils.coordinates import GenomicCoordinate
            >>> coord1 = GenomicCoordinate("chr1", 1000, 1200, '+', '1-based')
            >>> coord2 = GenomicCoordinate("chr1", 1500, 1700, '+', '1-based')
            >>> exon1 = Exon("exon1", coord1, parent_id="trans1")
            >>> exon2 = Exon("exon2", coord2, parent_id="trans1")
            >>> intron = Intron.from_exon_pair(exon1, exon2, "intron_1")
            >>> intron.start
            1201
            >>> intron.stop
            1499
            >>> intron.length
            299
        """
        # Validate exon compatibility
        if exon1.chromosome != exon2.chromosome:
            raise ValueError(
                f"Exons on different chromosomes: {exon1.chromosome} vs {exon2.chromosome}"
            )

        if exon1.strand != exon2.strand:
            raise ValueError(
                f"Exons on different strands: {exon1.strand} vs {exon2.strand}"
            )

        # Determine intron coordinates (gap between exons)
        # Use min/max to handle exons in any order (works for both strands)
        intron_start = min(exon1.stop, exon2.stop) + 1
        intron_stop = max(exon1.start, exon2.start) - 1

        if intron_start > intron_stop:
            raise ValueError(
                f"Exons overlap or touch: min_stop={min(exon1.stop, exon2.stop)}, "
                f"max_start={max(exon1.start, exon2.start)}"
            )

        # Create coordinate
        coord = GenomicCoordinate(
            exon1.chromosome, intron_start, intron_stop, exon1.strand, "1-based"
        )

        # Auto-generate ID if not provided
        if intron_id is None:
            parent_id = exon1.parent_id or "unknown"
            intron_id = f"{parent_id}:intron_{intron_start}_{intron_stop}"

        # Derive phase from upstream exon (CDS) phase annotation (if available)
        # Formula: phase = (exon_length - exon_phase) % 3
        # This calculates the reading frame at the start of the next exon
        phase = None
        if exon1.phase is not None:
            phase = (exon1.length - exon1.phase) % 3

        # Create metadata with parent info
        metadata = IntronMetadata(
            parent=exon1.parent_id,
            grandparent=None,  # Will be filled in later if available
            phase=phase,
        )

        return cls(intron_id=intron_id, coordinates=coord, metadata=metadata)

    # Convenience methods for updating mutable metadata

    def with_scores(self, scores: IntronScores) -> "Intron":
        """
        Create a new Intron with updated scores.

        Since Intron is frozen, this returns a new instance.

        Args:
            scores: IntronScores object

        Returns:
            New Intron with updated scores
        """
        return Intron(
            intron_id=self.intron_id,
            coordinates=self.coordinates,
            scores=scores,
            sequences=self.sequences,
            metadata=self.metadata,
        )

    def with_sequences(self, sequences: IntronSequences) -> "Intron":
        """
        Create a new Intron with updated sequences.

        Args:
            sequences: IntronSequences object

        Returns:
            New Intron with updated sequences
        """
        return Intron(
            intron_id=self.intron_id,
            coordinates=self.coordinates,
            scores=self.scores,
            sequences=sequences,
            metadata=self.metadata,
        )

    def with_metadata(self, metadata: IntronMetadata) -> "Intron":
        """
        Create a new Intron with updated metadata.

        Args:
            metadata: IntronMetadata object

        Returns:
            New Intron with updated metadata
        """
        return Intron(
            intron_id=self.intron_id,
            coordinates=self.coordinates,
            scores=self.scores,
            sequences=self.sequences,
            metadata=metadata,
        )

    def __str__(self) -> str:
        """Human-readable string representation."""
        parts = [f"Intron:{self.intron_id}"]
        parts.append(str(self.coordinates))

        if self.sequences:
            parts.append(str(self.sequences))

        if self.scores:
            parts.append(str(self.scores))

        return " | ".join(parts)

    def __repr__(self) -> str:
        """Developer-friendly representation."""
        return (
            f"Intron("
            f"intron_id='{self.intron_id}', "
            f"coordinates={self.coordinates!r}, "
            f"has_scores={self.has_scores}, "
            f"has_sequences={self.has_sequences})"
        )


if __name__ == "__main__":
    import doctest

    doctest.testmod()
