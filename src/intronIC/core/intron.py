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
from enum import IntEnum, IntFlag
from typing import Literal, Optional

from intronIC.utils.coordinates import GenomicCoordinate

# Type aliases for clarity
IntronType = Literal["u2", "u12", "unknown"]


class OmissionReason(IntEnum):
    """
    Omission reason codes (compact integer storage).

    These match the original intronIC single-character codes but store as integers.
    Conversion to human-readable strings happens only at write time.

    Examples:
        >>> reason = OmissionReason.SHORT
        >>> reason.value  # 1
        >>> reason.code   # 's'
        >>> reason.name   # 'SHORT'
        >>> reason.verbose  # 'omitted_short'
    """

    NONE = 0  # Not omitted (default)
    SHORT = 1  # 's' - Too short
    AMBIGUOUS = 2  # 'a' - Contains ambiguous bases
    NONCANONICAL = 3  # 'n' - Non-canonical splice sites
    ISOFORM = 4  # 'i' - Not from longest isoform
    OVERLAP = 5  # 'v' - Overlapping coordinates
    DUPLICATE = 6  # 'd' - Duplicate coordinates

    @property
    def code(self) -> str:
        """Get single-character code for backward compatibility."""
        _CODE_MAP = {
            OmissionReason.NONE: "",
            OmissionReason.SHORT: "s",
            OmissionReason.AMBIGUOUS: "a",
            OmissionReason.NONCANONICAL: "n",
            OmissionReason.ISOFORM: "i",
            OmissionReason.OVERLAP: "v",
            OmissionReason.DUPLICATE: "d",
        }
        return _CODE_MAP[self]

    @property
    def verbose(self) -> str:
        """Get verbose name for output files."""
        _VERBOSE_MAP = {
            OmissionReason.NONE: "",
            OmissionReason.SHORT: "omitted_short",
            OmissionReason.AMBIGUOUS: "omitted_ambiguous",
            OmissionReason.NONCANONICAL: "omitted_noncanonical",
            OmissionReason.ISOFORM: "omitted_not_longest_isoform",
            OmissionReason.OVERLAP: "omitted_overlap",
            OmissionReason.DUPLICATE: "duplicate",
        }
        return _VERBOSE_MAP[self]

    @classmethod
    def from_code(cls, code: str) -> "OmissionReason":
        """Create from single-character code."""
        _REVERSE_MAP = {
            "s": cls.SHORT,
            "a": cls.AMBIGUOUS,
            "n": cls.NONCANONICAL,
            "i": cls.ISOFORM,
            "v": cls.OVERLAP,
            "d": cls.DUPLICATE,
        }
        return _REVERSE_MAP.get(code, cls.NONE)


class IntronFlags(IntFlag):
    """
    Bit flags for intron properties (ultra-compact storage).

    Multiple flags can be combined using bitwise OR (|).
    Storage: single integer (28 bytes in Python, could be 4 bytes in numpy).

    Examples:
        >>> flags = IntronFlags.NONCANONICAL | IntronFlags.CORRECTED
        >>> IntronFlags.NONCANONICAL in flags  # True
        >>> IntronFlags.LONGEST_ISOFORM in flags  # False
    """

    NONE = 0
    NONCANONICAL = 1 << 0  # 0x001 - Non-standard splice sites
    LONGEST_ISOFORM = 1 << 1  # 0x002 - From longest transcript
    CORRECTED = 1 << 2  # 0x004 - Boundaries were adjusted
    DUPLICATE = 1 << 3  # 0x008 - Duplicate coordinates
    EDGE_CASE = 1 << 4  # 0x010 - Edge case marker
    SEQUENCE_ONLY = 1 << 5  # 0x020 - Loaded from sequence file (no real coordinates)


# Type alias for backward compatibility
OmissionCode = Optional[OmissionReason]


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


@dataclass(frozen=True, slots=True)
class ScoringMotifs:
    """
    Minimal sequence data needed for PWM scoring.

    This dataclass stores only the short sequence regions required for scoring,
    enabling memory-efficient streaming mode where full sequences are written
    to disk and only these motifs are kept in RAM.

    Memory comparison (per intron):
    - Full IntronSequences: ~1-10 KB (depends on intron length)
    - ScoringMotifs: ~100-150 bytes (fixed size)

    For human genome (~250k introns):
    - Full sequences in memory: ~500 MB - 2 GB
    - Motifs only: ~25-40 MB

    Attributes:
        five_region: 5' splice site region for PWM scoring
                     Default coords (-3, 9): 12bp (3bp exon + 9bp intron)
        three_region: 3' splice site region for PWM scoring
                      Default coords (-6, 4): 10bp (6bp intron + 4bp exon)
        bp_region: Branch point search region for PWM scoring
                   Default coords (-55, -5): 50bp from intron 3' end
        terminal_dnts: Terminal dinucleotides (e.g., 'GT-AG') for matrix selection
        upstream_flank: Upstream exonic flank (kept for final output formatting)
        downstream_flank: Downstream exonic flank (kept for final output formatting)

    Usage:
        # Extract motifs from full sequences before clearing them
        motifs = intron.extract_scoring_motifs(five_coords, bp_coords, three_coords)

        # In IntronScorer, use motifs if available
        if intron.motifs:
            five_region = intron.motifs.five_region
        else:
            five_region = extract_from_full_sequence(intron.sequences)

    Note:
        This is used in streaming mode (--streaming flag) where sequences are
        written to a temporary SQLite database during extraction, then merged
        with scores after classification to produce the final .introns.iic file.
    """

    five_region: str
    three_region: str
    bp_region: str
    terminal_dnts: str
    upstream_flank: Optional[str] = None
    downstream_flank: Optional[str] = None

    def __str__(self) -> str:
        """Human-readable string representation."""
        return (
            f"ScoringMotifs(5'={len(self.five_region)}bp, "
            f"BP={len(self.bp_region)}bp, "
            f"3'={len(self.three_region)}bp, "
            f"dnts={self.terminal_dnts})"
        )


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
        omitted: Omission reason (integer enum, 0=not omitted)
        duplicate: Reference to representative intron name if duplicate
        overlap: Overlapping coordinate record
        flags: Compact bit flags (noncanonical, longest_isoform, corrected, etc.)
        phase: Coding phase information
        defined_by: Which feature type defined this intron ('cds' or 'exon')
        upstream_exon_id: ID of exon/CDS feature upstream of this intron
        downstream_exon_id: ID of exon/CDS feature downstream of this intron
        fractional_position: Actual position in transcript based on cumulative exon lengths (0.0-1.0)
        dynamic_tags: Set of dynamic tags for special cases ([c:N], etc.)
        correction_distance: Distance boundaries were shifted (for [c:N] tag)

    Note:
        Mutable (not frozen) to allow updating tags during pipeline.
        Using integer codes and bit flags for memory efficiency.
    """

    parent: Optional[str] = None
    grandparent: Optional[str] = None
    index: Optional[int] = None
    family_size: Optional[int] = None
    parent_length: Optional[int] = None
    line_number: Optional[int] = (
        None  # Annotation line number for hierarchical sort tiebreaker
    )
    type_id: IntronType = "unknown"
    omitted: OmissionReason = (
        OmissionReason.NONE
    )  # Integer enum (4 bytes vs ~50 bytes for str)
    duplicate: Optional[str] = None
    overlap: Optional[str] = None
    flags: IntronFlags = IntronFlags.NONE  # Bit flags for boolean properties
    phase: Optional[int] = None
    defined_by: Optional[str] = (
        None  # 'cds' or 'exon' - which feature type defined this intron
    )
    upstream_exon_id: Optional[str] = None  # ID of exon/CDS upstream of this intron
    downstream_exon_id: Optional[str] = None  # ID of exon/CDS downstream of this intron
    fractional_position: Optional[float] = (
        None  # Actual position in transcript (0.0-1.0)
    )
    # Dynamic tagging system for special cases (rarely used, so overhead acceptable)
    dynamic_tags: set[str] = field(default_factory=set)
    correction_distance: Optional[int] = None

    # Convenience properties for backward compatibility
    @property
    def noncanonical(self) -> bool:
        """Check if intron has non-canonical boundaries."""
        return IntronFlags.NONCANONICAL in self.flags

    @noncanonical.setter
    def noncanonical(self, value: bool):
        """Set non-canonical flag."""
        if value:
            self.flags |= IntronFlags.NONCANONICAL
        else:
            self.flags &= ~IntronFlags.NONCANONICAL

    @property
    def longest_isoform(self) -> bool:
        """Check if from longest transcript."""
        return IntronFlags.LONGEST_ISOFORM in self.flags

    @longest_isoform.setter
    def longest_isoform(self, value: bool):
        """Set longest isoform flag."""
        if value:
            self.flags |= IntronFlags.LONGEST_ISOFORM
        else:
            self.flags &= ~IntronFlags.LONGEST_ISOFORM

    @property
    def corrected(self) -> bool:
        """Check if boundaries were adjusted."""
        return IntronFlags.CORRECTED in self.flags

    @corrected.setter
    def corrected(self, value: bool):
        """Set corrected flag."""
        if value:
            self.flags |= IntronFlags.CORRECTED
        else:
            self.flags &= ~IntronFlags.CORRECTED

    def is_omitted(self) -> bool:
        """Check if this intron should be omitted."""
        return self.omitted != OmissionReason.NONE

    def is_duplicate(self) -> bool:
        """Check if this intron is a duplicate."""
        return self.duplicate is not None

    def is_canonical(self) -> bool:
        """Check if this intron has canonical boundaries."""
        return not self.noncanonical

    def __str__(self) -> str:
        """Human-readable string representation."""
        parts = []
        if self.parent:
            parts.append(f"parent:{self.parent}")
        if self.index is not None and self.family_size is not None:
            parts.append(f"i{self.index}/{self.family_size}")
        if self.type_id != "unknown":
            parts.append(self.type_id)
        if self.is_omitted():
            parts.append(f"omitted:{self.omitted.code}")

        # Show flags
        flag_parts = []
        if self.noncanonical:
            flag_parts.append("NC")
        if self.longest_isoform:
            flag_parts.append("longest")
        if self.corrected:
            flag_parts.append("corrected")
        if flag_parts:
            parts.append(f"[{','.join(flag_parts)}]")

        return f"IntronMetadata({', '.join(parts)})" if parts else "IntronMetadata()"


@dataclass(frozen=True, slots=True)
class Intron:
    """
    Main Intron class using composition pattern.

    Composes four main components:
    - coordinates: GenomicCoordinate (location)
    - scores: IntronScores (classification data)
    - sequences: IntronSequences (full sequence data)
    - metadata: IntronMetadata (parent, tags, etc.)
    - motifs: ScoringMotifs (minimal sequence data for scoring, streaming mode)

    This design separates concerns and makes it easy to:
    - Create introns with partial data (e.g., coordinates only)
    - Update scores/sequences as they become available
    - Test components independently
    - Maintain clear boundaries between pipeline stages

    Streaming mode optimization:
        In streaming mode (--streaming flag), full sequences are written to
        a temporary SQLite database during extraction. Only the minimal
        scoring motifs (ScoringMotifs) are kept in memory, reducing RAM
        usage by ~85% for large genomes.

    Attributes:
        intron_id: Unique identifier
        coordinates: Genomic location
        scores: Scoring data (optional)
        sequences: Full sequence data (optional, cleared in streaming mode)
        metadata: Metadata (optional)
        motifs: Minimal scoring motifs (optional, used in streaming mode)

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

        >>> # In streaming mode, extract motifs before clearing sequences
        >>> motifs = intron2.extract_scoring_motifs((-3, 9), (-55, -5), (-6, 4))
        >>> intron3 = replace(motifs, sequences=None)  # Clear full sequences
        >>> intron3.motifs.five_region  # Still available for scoring
        'CAGGTAAGT...'
    """

    intron_id: str
    coordinates: GenomicCoordinate
    scores: Optional[IntronScores] = None
    sequences: Optional[IntronSequences] = None
    metadata: Optional[IntronMetadata] = None
    motifs: Optional[ScoringMotifs] = None

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
            motifs=self.motifs,
        )

    def with_sequences(self, sequences: Optional[IntronSequences]) -> "Intron":
        """
        Create a new Intron with updated sequences.

        Args:
            sequences: IntronSequences object (or None to clear)

        Returns:
            New Intron with updated sequences
        """
        return Intron(
            intron_id=self.intron_id,
            coordinates=self.coordinates,
            scores=self.scores,
            sequences=sequences,
            metadata=self.metadata,
            motifs=self.motifs,
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
            motifs=self.motifs,
        )

    def clear_sequences(self) -> "Intron":
        """
        Create a new Intron with large sequence fields cleared.

        This clears the memory-heavy sequence fields (seq, upstream_flank,
        downstream_flank, bp_region_seq) while preserving the small scoring
        sequences (five_seq, three_seq, bp_seq, bp_seq_u12) needed for
        scoring and classification.

        Memory savings per intron:
        - seq: ~500 bytes
        - upstream_flank: ~200 bytes
        - downstream_flank: ~200 bytes
        - bp_region_seq: ~50 bytes
        Total: ~950 bytes saved per intron

        Returns:
            New Intron with large sequences cleared (None)

        Examples:
            >>> intron = Intron("test", coord, sequences=IntronSequences(
            ...     seq="ATCG"*100, five_seq="GTAAGT", three_seq="TTTAG"
            ... ))
            >>> cleared = intron.clear_sequences()
            >>> cleared.sequences.seq is None
            True
            >>> cleared.sequences.five_seq  # Scoring sequences preserved
            'GTAAGT'
        """
        if self.sequences is None:
            return self  # No sequences to clear

        cleared_sequences = IntronSequences(
            # Clear large sequences
            seq=None,
            upstream_flank=None,
            downstream_flank=None,
            bp_region_seq=None,
            five_display_seq=None,
            three_display_seq=None,
            # Preserve small scoring sequences
            five_seq=self.sequences.five_seq,
            three_seq=self.sequences.three_seq,
            bp_seq=self.sequences.bp_seq,
            bp_seq_u2=self.sequences.bp_seq_u2,
            bp_relative_coords=self.sequences.bp_relative_coords,
            # Preserve stored dinucleotides (memory-efficient fields)
            five_prime_dnt=self.sequences.five_prime_dnt,
            three_prime_dnt=self.sequences.three_prime_dnt,
        )

        return self.with_sequences(cleared_sequences)

    def clear_all_sequences(self) -> "Intron":
        """
        Create a new Intron with ALL sequence fields cleared.

        This clears all sequence-related fields, including the small scoring
        sequences. Use this after scoring is complete and sequences are no
        longer needed.

        Memory savings per intron: ~13 KB (sequences + Python object overhead)

        Returns:
            New Intron with all sequences cleared (None)

        Examples:
            >>> intron = Intron("test", coord, sequences=IntronSequences(
            ...     seq="ATCG"*100, five_seq="GTAAGT", three_seq="TTTAG"
            ... ))
            >>> cleared = intron.clear_all_sequences()
            >>> cleared.sequences is None
            True
        """
        if self.sequences is None:
            return self  # No sequences to clear

        return self.with_sequences(None)

    def with_motifs(self, motifs: ScoringMotifs) -> "Intron":
        """
        Create a new Intron with scoring motifs.

        Args:
            motifs: ScoringMotifs object

        Returns:
            New Intron with motifs set
        """
        return Intron(
            intron_id=self.intron_id,
            coordinates=self.coordinates,
            scores=self.scores,
            sequences=self.sequences,
            metadata=self.metadata,
            motifs=motifs,
        )

    def extract_scoring_motifs(
        self,
        five_coords: tuple[int, int],
        bp_coords: tuple[int, int],
        three_coords: tuple[int, int],
    ) -> "Intron":
        """
        Extract minimal scoring motifs from full sequences.

        This extracts only the short sequence regions needed for PWM scoring,
        enabling streaming mode where full sequences are written to disk and
        only these motifs are kept in RAM.

        The extraction logic matches IntronScorer._extract_*_region methods
        to ensure scoring produces identical results.

        Args:
            five_coords: (start, stop) for 5' region relative to intron start
                         e.g., (-3, 9) means 3bp upstream + 9bp into intron
            bp_coords: (start, stop) for BP region relative to intron 3' end
                       e.g., (-55, -5) means positions -55 to -5 from end
            three_coords: (start, stop) for 3' region relative to intron end
                          e.g., (-6, 4) means 6bp from intron + 4bp downstream

        Returns:
            New Intron with motifs populated (sequences unchanged)

        Raises:
            ValueError: If sequences are not available

        Example:
            >>> intron_with_motifs = intron.extract_scoring_motifs(
            ...     five_coords=(-3, 9),
            ...     bp_coords=(-55, -5),
            ...     three_coords=(-6, 4)
            ... )
            >>> # Now safe to clear full sequences
            >>> lightweight = intron_with_motifs.with_sequences(None)
            >>> lightweight.motifs.five_region  # Still available for scoring
        """
        if self.sequences is None or self.sequences.seq is None:
            raise ValueError(
                f"Cannot extract scoring motifs from intron {self.intron_id}: "
                "no sequences available"
            )

        seq = self.sequences.seq
        upstream_flank = self.sequences.upstream_flank or ""
        downstream_flank = self.sequences.downstream_flank or ""

        # Extract 5' region (matching IntronScorer._extract_five_region)
        five_start, five_stop = five_coords
        if five_start < 0:
            # Need upstream flank
            upstream_needed = abs(five_start)
            if len(upstream_flank) >= upstream_needed:
                upstream_part = upstream_flank[-upstream_needed:]
            else:
                upstream_part = "N" * upstream_needed
            intron_part = seq[:five_stop] if five_stop > 0 else ""
            five_region = upstream_part + intron_part
        else:
            five_region = seq[five_start:five_stop]

        # Extract 3' region (matching IntronScorer._extract_three_region)
        three_start, three_stop = three_coords
        intron_length = len(seq)

        if three_start < 0:
            intron_start_idx = intron_length + three_start
        else:
            intron_start_idx = intron_length

        if three_stop <= 0:
            intron_stop_idx = intron_length + three_stop
        else:
            intron_stop_idx = intron_length

        if intron_start_idx < intron_stop_idx:
            intron_part = seq[intron_start_idx:intron_stop_idx]
        else:
            intron_part = ""

        if three_stop > 0:
            if len(downstream_flank) >= three_stop:
                downstream_part = downstream_flank[:three_stop]
            else:
                downstream_part = "N" * three_stop
            three_region = intron_part + downstream_part
        else:
            three_region = intron_part

        # Extract BP region (matching IntronScorer/BranchPointScorer._extract_search_region)
        bp_start, bp_stop = bp_coords
        # BP coords are relative to intron end (negative values)
        bp_start_idx = intron_length + bp_start
        bp_stop_idx = intron_length + bp_stop

        # Clamp to valid range, EXCLUDING the 5' splice site scoring region
        # This matches BranchPointScorer._extract_search_region() which does:
        #   start_pos = max(start_pos, five_end)
        # where five_end is the end position of the 5' scoring region.
        # This prevents the BP search from overlapping with the 5' splice site,
        # which is critical for short introns.
        five_end = five_coords[1]  # e.g., 9 for (-3, 9) - end of 5' region
        bp_start_idx = max(bp_start_idx, five_end)
        bp_stop_idx = max(0, min(intron_length, bp_stop_idx))

        bp_region = seq[bp_start_idx:bp_stop_idx] if bp_start_idx < bp_stop_idx else ""

        # Get terminal dinucleotides
        terminal_dnts = self.sequences.terminal_dinucleotides or ""

        motifs = ScoringMotifs(
            five_region=five_region,
            three_region=three_region,
            bp_region=bp_region,
            terminal_dnts=terminal_dnts,
            upstream_flank=upstream_flank,
            downstream_flank=downstream_flank,
        )

        return self.with_motifs(motifs)

    @property
    def has_motifs(self) -> bool:
        """Check if scoring motifs are available."""
        return self.motifs is not None

    def __str__(self) -> str:
        """Human-readable string representation."""
        parts = [f"Intron:{self.intron_id}"]
        parts.append(str(self.coordinates))

        if self.sequences:
            parts.append(str(self.sequences))
        elif self.motifs:
            parts.append(str(self.motifs))

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
