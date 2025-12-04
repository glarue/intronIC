"""
Writers for various intronIC output formats.

This module provides writers for:
- BED format (.bed.iic)
- Metadata format (.meta.iic)
- Sequence format (.introns.iic)
- Score details format (.score_info.iic)
- Mapping files (.dupe_map.iic, .overlap_map.iic)

Design principle: Generator-friendly - accept iterables and write one
item at a time to minimize memory usage.

Author: intronIC refactoring project
Date: 2025-11-02
"""

from pathlib import Path
from typing import Dict, Iterable, Optional, TextIO, Tuple, Union

from intronIC.core.intron import Intron, IntronFlags, OmissionReason


# ============================================================================
# Base Writer Class
# ============================================================================


class BaseWriter:
    """Base class for all intronIC output writers with shared file handling."""

    def __init__(self, file_path: Union[str, Path], mode: str = "w"):
        """
        Initialize writer.

        Args:
            file_path: Path to output file
            mode: File open mode ('w' for write, 'a' for append)
        """
        self.file_path = Path(file_path)
        self.file: Optional[TextIO] = None
        self.items_written = 0
        self._mode = mode

    def open(self) -> None:
        """Open output file for writing."""
        self.file = open(self.file_path, self._mode)

    def close(self) -> None:
        """Close output file."""
        if self.file:
            self.file.close()
            self.file = None

    def __enter__(self):
        """Context manager entry."""
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()

    def _require_open(self):
        """Raise error if file is not open."""
        if not self.file:
            raise ValueError("File not open. Call open() first or use context manager.")


# ============================================================================
# Attribute Tag Mapping
# ============================================================================


# NOTE: This mapping is now primarily for backward compatibility and dynamic tags.
# For omission reasons, we use OmissionReason.verbose property.
# For flags, we check IntronFlags directly.
TAG_TO_ATTRIBUTE: Dict[str, str] = {
    # Dynamic tags (still used for special cases)
    "c": "corrected",  # [c] - Splice site boundaries were adjusted
    "d": "duplicate",  # [d] - Duplicate coordinates (excluded from analysis)
    "e": "edge_case",  # [e] - Edge case marker
}


def generate_attributes(intron: Intron) -> str:
    """
    Generate verbose comma-separated attributes string from intron metadata.

    Uses the new enum-based system for compact and auditable storage.

    Args:
        intron: Intron object

    Returns:
        Comma-separated attributes (e.g., "noncanonical,not_longest_isoform")
        or '.' if no attributes

    Examples:
        >>> coord = GenomicCoordinate("chr1", 1000, 2000, '+', '1-based')
        >>> metadata = IntronMetadata(noncanonical=True, longest_isoform=False)
        >>> intron = Intron("test", coord, metadata=metadata)
        >>> generate_attributes(intron)
        'noncanonical,not_longest_isoform'
    """
    attrs = []

    if intron.metadata:
        # Check flags (uses properties that access IntronFlags)
        if intron.metadata.noncanonical:
            attrs.append("noncanonical")
        if not intron.metadata.longest_isoform:
            attrs.append("not_longest_isoform")
        if intron.metadata.corrected:
            attrs.append("corrected")
        if intron.metadata.duplicate:
            attrs.append("duplicate")

        # Omission reason (uses OmissionReason enum)
        if intron.metadata.is_omitted():
            attrs.append(intron.metadata.omitted.verbose)

    return ",".join(attrs) if attrs else "NA"


# ============================================================================
# Formatting Helper Functions
# ============================================================================


def generate_species_abbreviation(species_name: str) -> str:
    """
    Generate 3+3 character abbreviation from species name.

    Follows original intronIC convention: first 3 chars of genus +
    first 3 chars of species epithet, with title case.

    Port from: intronIC.py (inferred from output format)

    Args:
        species_name: Full species name (e.g., "homo_sapiens" or "Homo sapiens")

    Returns:
        6-character abbreviation (e.g., "HomSap")

    Examples:
        >>> generate_species_abbreviation("homo_sapiens")
        'HomSap'
        >>> generate_species_abbreviation("drosophila_melanogaster")
        'DroMel'
        >>> generate_species_abbreviation("c_elegans")
        'CEle'
        >>> generate_species_abbreviation("Arabidopsis thaliana")
        'AraTha'
    """
    # Normalize: replace underscores with spaces, split on whitespace
    parts = species_name.replace("_", " ").strip().split()

    if len(parts) >= 2:
        # Standard binomial nomenclature: Genus species
        genus = parts[0]
        species = parts[1]

        # First 3 chars of each, with title case (capitalize first letter)
        genus_abbrev = genus[:3].capitalize()
        species_abbrev = species[:3].capitalize()

        return genus_abbrev + species_abbrev

    elif len(parts) == 1:
        # Single name (unusual case)
        name = parts[0]
        if len(name) >= 6:
            # Take first 6 chars
            return name[:6].capitalize()
        else:
            # Pad short names with 'X'
            return (name[:1].upper() + name[1:].lower()).ljust(6, "X")

    else:
        # Empty or invalid - return placeholder
        return "XXXXXX"


def format_omission_tag(omitted: OmissionReason) -> str:
    """
    Format omission tag for intron name.

    Port from: intronIC.py:629-632
    Updated to use OmissionReason enum.

    Args:
        omitted: OmissionReason enum value

    Returns:
        Formatted tag like ';[o:s]' or empty string if not omitted

    Examples:
        >>> format_omission_tag(OmissionReason.SHORT)
        ';[o:s]'
        >>> format_omission_tag(OmissionReason.NONCANONICAL)
        ';[o:n]'
        >>> format_omission_tag(OmissionReason.NONE)
        ''
    """
    if omitted != OmissionReason.NONE:
        return f";[o:{omitted.code}]"
    return ""


def format_dynamic_tags(tags: set[str]) -> str:
    """
    Format dynamic tags for intron name.

    Port from: intronIC.py:633-636

    Dynamic tags track various intron states:
    - [c:N]: Boundary corrected by N bases
    - [d]: Duplicate marker
    - [e]: Edge case
    - [n]: Non-canonical
    - [i]: Isoform-related
    - Terminal dinucleotides (e.g., GC-AG)

    Args:
        tags: Set of tag strings (may or may not have brackets)

    Returns:
        Formatted tag string like ';[c:5];[d]' or empty string if no tags

    Examples:
        >>> format_dynamic_tags({'[c:5]', '[d]'})
        ';[c:5];[d]'
        >>> format_dynamic_tags({'c:5', 'd'})
        ';[c:5];[d]'
        >>> format_dynamic_tags(set())
        ''
        >>> format_dynamic_tags({'[n]', 'GC-AG'})
        ';[GC-AG];[n]'
    """
    if not tags:
        return ""

    # Ensure all tags have brackets
    formatted_tags = []
    for tag in sorted(tags):  # Sort for deterministic output
        if not tag.startswith("["):
            tag = f"[{tag}]"
        formatted_tags.append(tag)

    return ";" + ";".join(formatted_tags)


def generate_intron_name(
    intron: Intron,
    species_name: Optional[str] = None,
    simple_name: bool = False,
    no_abbreviate: bool = False,
    intron_number: Optional[int] = None,
) -> str:
    """
    Generate standardized intron name for all output formats.

    Port from: intronIC.py:622-646 (get_name)

    Format: {species_abbrev}-{grandparent}@{parent}-intron_{index}({family_size}){omit_tag}{dynamic_tags}

    Example: HomSap-gene:ENSG00000196218@transcript:ENST00000355481-intron_69(104);[o:i];[i]

    Args:
        intron: Intron object
        species_name: Full species name (e.g., "homo_sapiens")
        simple_name: Use simplified format (species-i_{number})
        no_abbreviate: Use full species name instead of abbreviation
        intron_number: Simple incrementing number (used with simple_name)

    Returns:
        Formatted intron name string

    Examples:
        >>> generate_intron_name(intron, "homo_sapiens", False)
        'HomSap-gene:ENSG00000196218@transcript:ENST00000355481_69(104)'
        >>> generate_intron_name(intron, "homo_sapiens", False, True)
        'homo_sapiens-gene:ENSG00000196218@transcript:ENST00000355481_69(104)'
        >>> generate_intron_name(omitted_intron, "homo_sapiens", False)
        'HomSap-gene:ENSG00000196218@transcript:ENST00000355481_69(104);[o:i]'
        >>> generate_intron_name(intron, "homo_sapiens", True, False, 123)
        'HomSap-i123'
    """
    if simple_name:
        # Simple format: species-i{number}{tags}
        if no_abbreviate:
            species_prefix = species_name if species_name else "unknown"
        else:
            species_prefix = (
                generate_species_abbreviation(species_name)
                if species_name
                else "XXXXXX"
            )
        omit_tag = (
            format_omission_tag(intron.metadata.omitted) if intron.metadata else ""
        )
        dyn_tag = (
            format_dynamic_tags(intron.metadata.dynamic_tags) if intron.metadata else ""
        )
        # Use intron_number if provided, otherwise fall back to intron_id
        identifier = (
            str(intron_number) if intron_number is not None else intron.intron_id
        )
        return f"{species_prefix}-i{identifier}{omit_tag}{dyn_tag}"

    # Full format
    if not intron.metadata:
        # Fallback if no metadata - return intron_id
        return intron.intron_id

    # If metadata exists but all key fields are None (loaded from sequence file),
    # fall back to intron_id to preserve original naming
    if (
        intron.metadata.grandparent is None
        and intron.metadata.parent is None
        and intron.metadata.index is None
    ):
        # This is likely loaded from a sequence file - preserve original name
        omit_tag = format_omission_tag(intron.metadata.omitted)
        dyn_tag = format_dynamic_tags(intron.metadata.dynamic_tags)
        return f"{intron.intron_id}{omit_tag}{dyn_tag}"

    # Species prefix (abbreviated or full depending on flag)
    if no_abbreviate:
        species_prefix = species_name if species_name else "unknown"
    else:
        species_prefix = (
            generate_species_abbreviation(species_name) if species_name else "XXXXXX"
        )

    # Gene ID (grandparent) - preserve "gene:" prefix if present
    grandparent = intron.metadata.grandparent if intron.metadata.grandparent else "?"

    # Transcript ID (parent) - preserve "transcript:" prefix if present
    parent = intron.metadata.parent if intron.metadata.parent else "?"

    # Index and family size
    index = intron.metadata.index if intron.metadata.index is not None else "?"
    family_size = (
        intron.metadata.family_size if intron.metadata.family_size is not None else "?"
    )

    # Format tags
    omit_tag = format_omission_tag(intron.metadata.omitted)
    dyn_tag = format_dynamic_tags(intron.metadata.dynamic_tags)

    # Build name: species-grandparent@parent_index(family_size)tags
    # Note: Removed "-intron" as it doesn't add information and takes up space
    name = f"{species_prefix}-{grandparent}@{parent}_{index}({family_size}){omit_tag}{dyn_tag}"

    return name


def generate_intron_label(
    intron: Intron,
    species_name: Optional[str] = None,
    simple_name: bool = False,
    no_abbreviate: bool = False,
    intron_number: Optional[int] = None,
) -> str:
    """
    Generate intron label for BED output (name + score).

    Port from: intronIC.py:649-659 (get_label)

    Format: {species-}parent_index(family_size);svm_score[;tags]
    Or with simple_name + intron_number: {species-}i_{number};svm_score[;tags]

    Args:
        intron: Intron object
        species_name: Species name (optional)
        simple_name: Exclude species prefix and use simple numbering
        no_abbreviate: Use full species name instead of abbreviation
        intron_number: Simple incrementing number (used with simple_name)

    Returns:
        Intron label string

    Examples:
        >>> generate_intron_label(intron, "homo_sapiens", False)
        'HomSap-ENST00000397910_1(83);0.00'
        >>> generate_intron_label(intron, "homo_sapiens", False, True)
        'homo_sapiens-ENST00000397910_1(83);0.00'
        >>> generate_intron_label(intron, "homo_sapiens", True, False, 123)
        'HomSap-i_123;0.00'
    """
    parts = []

    # Add species prefix (even with simple_name, we still show species)
    if species_name:
        if no_abbreviate:
            parts.append(species_name)
        else:
            species_abbrev = generate_species_abbreviation(species_name)
            parts.append(species_abbrev)

    # Use simple numbering if requested and number provided
    if simple_name and intron_number is not None:
        parts.append(f"i_{intron_number}")
    # Otherwise use parent_index(family_size) format
    elif intron.metadata and intron.metadata.parent:
        parent = intron.metadata.parent
        index = intron.metadata.index if intron.metadata.index else 1
        family_size = intron.metadata.family_size if intron.metadata.family_size else 1
        parts.append(f"{parent}_{index}({family_size})")
    else:
        parts.append(intron.intron_id)

    # Join with dash, not underscore (matches original format)
    name = "-".join(parts) if parts else intron.intron_id

    # Add SVM score if available
    if intron.svm_score is not None:
        name += f";{intron.svm_score:.2f}"

    # Add tags with semicolon separator
    tags = []
    if intron.metadata:
        if intron.metadata.noncanonical:
            tags.append("[n]")
        if not intron.metadata.longest_isoform:
            tags.append("[i]")
        if intron.metadata.corrected:
            if intron.metadata.correction_distance is not None:
                tags.append(f"[c:{intron.metadata.correction_distance}]")
            else:
                tags.append("[c]")
        if intron.metadata.duplicate:
            tags.append("[d]")
        if intron.metadata.omitted != OmissionReason.NONE:
            tags.append(f"[o:{intron.metadata.omitted.code}]")

    if tags:
        name += ";" + ";".join(tags)

    return name


def generate_motif_schematic(intron: Intron, exonic: int = 3) -> str:
    """
    Generate motif schematic string for .meta.iic output.

    Port from: intronIC.py:725-742 (motif_string)

    Format: {exon_3bp}|{5'_10bp}...{bp_u12}/{bp_u2}...{3'_display}|{exon_3bp}

    Example: AAG|GTCGGGGCTT...TACTAAC/CACAG...TTTAG|TCC

    Args:
        intron: Intron object with sequences
        exonic: Number of exonic bases to show (default: 3)

    Returns:
        Motif schematic string or 'NA' if sequences missing

    Examples:
        >>> # intron with all sequences populated
        >>> schematic = generate_motif_schematic(intron)
        >>> schematic
        'AAG|GTCGGGGCTT...TACTAAC/CACAG...TTTAG|TCC'
    """
    if not intron.sequences:
        return "NA"

    seqs = intron.sequences

    # Check for required sequences
    if not all(
        [
            seqs.upstream_flank,
            seqs.five_display_seq,
            seqs.three_display_seq,
            seqs.downstream_flank,
        ]
    ):
        return "NA"

    # Five prime boundary: {last 3bp of exon}|{first 10bp of intron}
    five_boundary = f"{seqs.upstream_flank[-exonic:]}|{seqs.five_display_seq}"

    # Three prime boundary: {last Nbp of intron}|{first 3bp of exon}
    three_boundary = f"{seqs.three_display_seq}|{seqs.downstream_flank[:exonic]}"

    # Branch point display: {U12_bp}/{U2_bp} or just {U12_bp} if U2 missing
    bps_display = None
    if seqs.bp_seq and seqs.bp_seq_u2:
        bps_display = f"{seqs.bp_seq}/{seqs.bp_seq_u2}"
    elif seqs.bp_seq:
        bps_display = seqs.bp_seq

    # Assemble schematic with '...' separators
    schematic_parts = [five_boundary]
    if bps_display:
        schematic_parts.append(bps_display)
    schematic_parts.append(three_boundary)

    return "...".join(schematic_parts)


def annotate_sequence(sequence: str, start: int, stop: int) -> str:
    """
    Add brackets around substring.

    Port from: intronIC.py:3208-3211

    Args:
        sequence: Full sequence
        start: Start position (0-based)
        stop: Stop position (0-based, exclusive)

    Returns:
        Annotated sequence with [brackets] around substring

    Examples:
        >>> annotate_sequence("ABCDEFGH", 2, 5)
        'AB[CDE]FGH'
        >>> annotate_sequence("TTGACAGGTACTAACGACTGA", 8, 15)
        'TTGACAGG[TACTAAC]GACTGA'
    """
    return sequence[:start] + "[" + sequence[start:stop] + "]" + sequence[stop:]


def generate_bp_context(intron: Intron) -> str:
    """
    Generate branch point context string for .meta.iic output.

    Port from: intronIC.py:744-752 (bps_context), 3223-3226

    Format: {bp_region with [brackets] around bp_seq} + {three_display_seq}

    Example: TTGACAGGCAGTGATAT[TACTAAC]GACTGAGTTTAG

    The BP sequence is wrapped in brackets within the bp_region_seq,
    then the three_display_seq is appended.

    Args:
        intron: Intron object with sequences

    Returns:
        BP context string or 'NA' if sequences missing

    Examples:
        >>> # intron with BP information populated
        >>> context = generate_bp_context(intron)
        >>> context
        'TTGACAGGCAGTGATAT[TACTAAC]GACTGAGTTTAG'
    """
    if not intron.sequences:
        return "NA"

    seqs = intron.sequences

    # Check for required sequences
    if not all([seqs.bp_region_seq, seqs.bp_relative_coords, seqs.three_display_seq]):
        return "NA"

    try:
        start, stop = seqs.bp_relative_coords
        # Annotate BP region with brackets around the BP sequence
        annotated_bp_region = annotate_sequence(seqs.bp_region_seq, start, stop)
        # Append three_display_seq
        context = annotated_bp_region + seqs.three_display_seq
        return context
    except Exception:
        # If any error (e.g., invalid coordinates), return placeholder
        return "NA"


# ============================================================================
# BED Format Writer
# ============================================================================


class BEDWriter(BaseWriter):
    """
    Writer for BED format output (.bed.iic).

    BED format (6 columns):
        chrom  start  stop  name  score  strand

    Notes:
        - Start is 0-based (BED convention), stop is 1-based
        - Score is SVM score (0-100%) or '.' if unavailable
        - Name includes intron label with tags

    Examples:
        >>> from pathlib import Path
        >>> writer = BEDWriter(Path("output.bed"))
        >>> writer.write_header()  # Optional, BED has no standard header
        >>> # writer.write_intron(intron)
        >>> writer.close()
    """

    @property
    def introns_written(self) -> int:
        """Backward-compatible alias for items_written."""
        return self.items_written

    def write_header(self) -> None:
        """Write BED header (optional). BED has no standard header."""
        self._require_open()

    def write_intron(
        self,
        intron: Intron,
        species_name: Optional[str] = None,
        simple_name: bool = False,
        no_abbreviate: bool = False,
        intron_number: Optional[int] = None,
    ) -> None:
        """
        Write a single intron in BED format.

        Args:
            intron: Intron object to write
            species_name: Species name for intron label (optional)
            simple_name: Use simple naming (no species prefix)
            no_abbreviate: Use full species name instead of abbreviation
            intron_number: Simple incrementing number (used with simple_name)

        Format:
            chrom  start(0-based)  stop  name  svm_score  strand  attributes
        """
        self._require_open()

        # Get BED start (0-based)
        start_0based = intron.start - 1

        # Get SVM score or 'NA' if unavailable
        score = "NA" if intron.svm_score is None else str(intron.svm_score)

        # Generate intron name using same format as meta.iic
        name = generate_intron_name(
            intron, species_name, simple_name, no_abbreviate, intron_number
        )

        # Generate verbose attributes
        attributes = generate_attributes(intron)

        # Write BED line
        fields = [
            intron.chromosome,
            str(start_0based),
            str(intron.stop),
            name,
            score,
            intron.strand,
            attributes,
        ]
        self.file.write("\t".join(fields) + "\n")
        self.items_written += 1

    def write_introns(
        self,
        introns: Iterable[Intron],
        species_name: Optional[str] = None,
        simple_name: bool = False,
        no_abbreviate: bool = False,
    ) -> int:
        """
        Write multiple introns.

        Args:
            introns: Iterable of Intron objects
            species_name: Species name for labels
            simple_name: Use simple naming
            no_abbreviate: Use full species name instead of abbreviation

        Returns:
            Number of introns written
        """
        count = 0
        for idx, intron in enumerate(introns, start=1):
            # Pass enumeration counter when using simple_name
            intron_number = idx if simple_name else None
            self.write_intron(
                intron, species_name, simple_name, no_abbreviate, intron_number
            )
            count += 1
        return count


# ============================================================================
# Metadata Format Writer
# ============================================================================


class MetaWriter(BaseWriter):
    """
    Writer for metadata format output (.meta.iic).

    Format (tab-delimited):
        name  rel_score  dnts  motif  bp_context  length  parent  grandparent
        index  family_size  frac_pos  phase  type_id  feature

    This comprehensive format includes all intron metadata for downstream analysis.

    Examples:
        >>> writer = MetaWriter(Path("output.meta.iic"))
        >>> # with writer:
        >>> #     writer.write_header()
        >>> #     writer.write_introns(introns)
    """

    @property
    def introns_written(self) -> int:
        """Backward-compatible alias for items_written."""
        return self.items_written

    def write_header(self) -> None:
        """Write metadata file header."""
        self._require_open()

        header_fields = [
            "name",
            "rel_score",
            "dnts",
            "motif_schematic",
            "bp_context",
            "length",
            "parent",
            "grandparent",
            "index",
            "family_size",
            "frac_pos",
            "phase",
            "type_id",
            "feature",
            "attributes",
        ]
        self.file.write("\t".join(header_fields) + "\n")

    def write_intron(
        self,
        intron: Intron,
        species_name: Optional[str] = None,
        simple_name: bool = False,
        no_abbreviate: bool = False,
        null: str = "NA",
        intron_number: Optional[int] = None,
    ) -> None:
        """
        Write a single intron's metadata.

        Args:
            intron: Intron object to write
            species_name: Species name for intron name
            simple_name: Use simple naming
            no_abbreviate: Use full species name instead of abbreviation
            null: Placeholder for missing values
            intron_number: Simple incrementing number (used with simple_name)

        Format:
            name  rel_score  dnts  motif  bp_context  length  parent  grandparent
            index  family_size  frac_pos  phase  type_id  feature
        """
        self._require_open()

        # Generate intron name using shared function
        name = generate_intron_name(
            intron, species_name, simple_name, no_abbreviate, intron_number
        )

        # Relative score (rounded to 4 decimal places)
        rel_score = null
        if intron.scores and intron.scores.relative_score is not None:
            rel_score = str(round(intron.scores.relative_score, 4))

        # Terminal dinucleotides (e.g., 'GT-AG')
        dnts = null
        if intron.sequences and intron.sequences.terminal_dinucleotides:
            dnts = intron.sequences.terminal_dinucleotides

        # Motif schematic (using new formatting function)
        motif = generate_motif_schematic(intron)

        # Branch point context (using new formatting function)
        bp_context = generate_bp_context(intron)

        # Length
        length = str(intron.length)

        # Parent/grandparent
        parent = null
        grandparent = null
        index = null
        family_size = null
        frac_pos = null
        phase = null

        if intron.metadata:
            parent = intron.metadata.parent if intron.metadata.parent else null
            grandparent = (
                intron.metadata.grandparent if intron.metadata.grandparent else null
            )
            index = (
                str(intron.metadata.index)
                if intron.metadata.index is not None
                else null
            )
            family_size = (
                str(intron.metadata.family_size)
                if intron.metadata.family_size
                else null
            )
            frac_pos_val = intron.metadata.fractional_position
            frac_pos = str(round(frac_pos_val, 4)) if frac_pos_val is not None else null
            phase = (
                str(intron.metadata.phase)
                if intron.metadata.phase is not None
                else null
            )

        # Type ID - write '.' if unknown (for omitted introns)
        type_id = null if intron.type_id == "unknown" else intron.type_id

        # Feature type (exon/cds) - defined_by field tracks which feature type defined this intron
        feature = null
        if intron.metadata and intron.metadata.defined_by:
            feature = intron.metadata.defined_by

        # Generate verbose attributes
        attributes = generate_attributes(intron)

        fields = [
            name,
            rel_score,
            dnts,
            motif,
            bp_context,
            length,
            parent,
            grandparent,
            index,
            family_size,
            frac_pos,
            phase,
            type_id,
            feature,
            attributes,
        ]

        self.file.write("\t".join(fields) + "\n")
        self.items_written += 1

    def write_introns(
        self,
        introns: Iterable[Intron],
        species_name: Optional[str] = None,
        simple_name: bool = False,
    ) -> int:
        """Write multiple introns. Returns number written."""
        count = 0
        for idx, intron in enumerate(introns, start=1):
            intron_number = idx if simple_name else None
            self.write_intron(
                intron, species_name, simple_name, intron_number=intron_number
            )
            count += 1
        return count


# ============================================================================
# Sequence Format Writer
# ============================================================================


class SequenceWriter(BaseWriter):
    """
    Writer for sequence format output (.introns.iic).

    Format (tab-delimited):
        name  [score]  upstream_flank  sequence  downstream_flank

    This format stores intron sequences with flanking exonic sequences
    for downstream analysis or re-scoring.

    Examples:
        >>> writer = SequenceWriter(Path("output.introns.iic"))
        >>> # with writer:
        >>> #     writer.write_introns(introns, include_score=True)
    """

    @property
    def introns_written(self) -> int:
        """Backward-compatible alias for items_written."""
        return self.items_written

    def write_intron(
        self,
        intron: Intron,
        species_name: Optional[str] = None,
        simple_name: bool = False,
        no_abbreviate: bool = False,
        include_score: bool = True,
        intron_number: Optional[int] = None,
    ) -> None:
        """
        Write a single intron's sequences.

        Args:
            intron: Intron object to write
            species_name: Species name for intron name
            simple_name: Use simple naming
            no_abbreviate: Use full species name instead of abbreviation
            include_score: Include SVM score in output
            intron_number: Simple incrementing number (used with simple_name)

        Format:
            name  [score]  upstream_flank  sequence  downstream_flank
        """
        self._require_open()

        if not intron.sequences or not intron.sequences.seq:
            raise ValueError(f"Intron {intron.intron_id} has no sequence data")

        # Generate intron name using shared function
        name = generate_intron_name(
            intron, species_name, simple_name, no_abbreviate, intron_number
        )

        # Get sequences (with defaults)
        upstream = intron.sequences.upstream_flank or ""
        sequence = intron.sequences.seq
        downstream = intron.sequences.downstream_flank or ""

        fields = [name]

        # Optionally include score
        if include_score:
            score = str(intron.svm_score) if intron.svm_score is not None else "NA"
            fields.append(score)

        fields.extend([upstream, sequence, downstream])

        self.file.write("\t".join(fields) + "\n")
        self.items_written += 1

    def write_from_row(
        self,
        intron_id: str,
        formatted_name: str,
        upstream_flank: Optional[str],
        seq: str,
        downstream_flank: Optional[str],
        terminal_dnts: Optional[str],
        svm_score: Optional[float],
    ) -> None:
        """
        Write intron sequence data from raw row values (for streaming mode).

        This method writes sequence data directly from SQLite SequenceRow values
        without requiring a full Intron object. Used when reading sequences back
        from the streaming SQLite store.

        Args:
            intron_id: Unique intron identifier (used for score lookup)
            formatted_name: Pre-computed formatted output name
            upstream_flank: Upstream flanking sequence
            seq: Intron sequence
            downstream_flank: Downstream flanking sequence
            terminal_dnts: Terminal dinucleotides (e.g., "GT-AG")
            svm_score: SVM classification score (or None if not scored)

        Format:
            name  score  upstream_flank  sequence  downstream_flank
        """
        self._require_open()

        fields = [formatted_name]

        # Include score
        score_str = f"{svm_score:.2f}" if svm_score is not None else "NA"
        fields.append(score_str)

        # Add sequences (with defaults for None)
        fields.extend([upstream_flank or "", seq, downstream_flank or ""])

        self.file.write("\t".join(fields) + "\n")
        self.items_written += 1

    def write_introns(
        self,
        introns: Iterable[Intron],
        species_name: Optional[str] = None,
        simple_name: bool = False,
        no_abbreviate: bool = False,
        include_score: bool = True,
    ) -> int:
        """
        Write multiple introns.

        Args:
            introns: Iterable of Intron objects
            species_name: Species name
            simple_name: Use simple naming
            no_abbreviate: Use full species name instead of abbreviation
            include_score: Include SVM score

        Returns:
            Number of introns written
        """
        count = 0
        for idx, intron in enumerate(introns, start=1):
            # Pass enumeration counter when using simple_name
            intron_number = idx if simple_name else None
            self.write_intron(
                intron,
                species_name,
                simple_name,
                no_abbreviate,
                include_score,
                intron_number,
            )
            count += 1
        return count


# ============================================================================
# Score Details Writer
# ============================================================================


class ScoreWriter(BaseWriter):
    """
    Writer for detailed scoring information (.score_info.iic).

    Format (tab-delimited):
        name  rel_score  svm_score  decision_dist  5'_seq  5'_raw  5'_z
        bp_seq  bp_region  bp_raw  bp_z  3'_seq  3'_raw  3'_z

    This comprehensive format includes all scoring details for
    in-depth analysis and debugging.

    Examples:
        >>> writer = ScoreWriter(Path("output.score_info.iic"))
        >>> # with writer:
        >>> #     writer.write_header()
        >>> #     writer.write_introns(introns)
    """

    @property
    def introns_written(self) -> int:
        """Backward-compatible alias for items_written."""
        return self.items_written

    def write_header(self) -> None:
        """Write score file header."""
        self._require_open()

        header_fields = [
            "name",
            "rel_score",
            "svm_score",
            "5'_seq",
            "5'_raw",
            "5'_z",
            "bp_seq",
            "bp_seq_u2",
            "bp_raw",
            "bp_z",
            "3'_seq",
            "3'_raw",
            "3'_z",
            "min(5,bp)",
            "min(5,3)",
            "max(5,bp)",
            "max(5,3)",
            "decision_dist",
        ]
        self.file.write("\t".join(header_fields) + "\n")

    def write_intron(
        self,
        intron: Intron,
        species_name: Optional[str] = None,
        simple_name: bool = False,
        no_abbreviate: bool = False,
        null: str = "NA",
    ) -> None:
        """
        Write a single intron's detailed scores.

        Args:
            intron: Intron object to write
            species_name: Species name for intron name
            simple_name: Use simple naming
            no_abbreviate: Use full species name instead of abbreviation
            null: Placeholder for missing values
        """
        self._require_open()

        # Generate intron name using shared function
        name = generate_intron_name(intron, species_name, simple_name, no_abbreviate)

        # Default all values to null
        rel_score = null
        svm_score = null
        decision_dist = null
        five_seq = null
        five_raw = null
        five_z = null
        bp_seq = null
        bp_seq_u2 = null
        bp_raw = null
        bp_z = null
        three_seq = null
        three_raw = null
        three_z = null
        min_5_bp = null
        min_5_3 = null
        max_5_bp = null
        max_5_3 = null

        # Fill in scores if available
        if intron.scores:
            if intron.scores.relative_score is not None:
                rel_score = str(round(intron.scores.relative_score, 4))
            if intron.scores.svm_score is not None:
                svm_score = str(round(intron.scores.svm_score, 2))
            if intron.scores.decision_distance is not None:
                decision_dist = str(round(intron.scores.decision_distance, 4))

            # Five prime scores
            if intron.scores.five_raw_score is not None:
                five_raw = str(round(intron.scores.five_raw_score, 6))
            if intron.scores.five_z_score is not None:
                five_z = str(round(intron.scores.five_z_score, 4))

            # Branch point scores
            if intron.scores.bp_raw_score is not None:
                bp_raw = str(round(intron.scores.bp_raw_score, 6))
            if intron.scores.bp_z_score is not None:
                bp_z = str(round(intron.scores.bp_z_score, 4))

            # Three prime scores
            if intron.scores.three_raw_score is not None:
                three_raw = str(round(intron.scores.three_raw_score, 6))
            if intron.scores.three_z_score is not None:
                three_z = str(round(intron.scores.three_z_score, 4))

            # BothEndsStrong augmented features
            if intron.scores.min_5_bp is not None:
                min_5_bp = str(round(intron.scores.min_5_bp, 4))
            if intron.scores.min_5_3 is not None:
                min_5_3 = str(round(intron.scores.min_5_3, 4))
            if intron.scores.max_5_bp is not None:
                max_5_bp = str(round(intron.scores.max_5_bp, 4))
            if intron.scores.max_5_3 is not None:
                max_5_3 = str(round(intron.scores.max_5_3, 4))

        # Fill in sequences if available
        if intron.sequences:
            if intron.sequences.five_seq:
                five_seq = intron.sequences.five_seq
            if intron.sequences.bp_seq:
                bp_seq = intron.sequences.bp_seq
            if intron.sequences.bp_seq_u2:
                bp_seq_u2 = intron.sequences.bp_seq_u2
            if intron.sequences.three_seq:
                three_seq = intron.sequences.three_seq

        fields = [
            name,
            rel_score,
            svm_score,
            five_seq,
            five_raw,
            five_z,
            bp_seq,
            bp_seq_u2,
            bp_raw,
            bp_z,
            three_seq,
            three_raw,
            three_z,
            min_5_bp,
            min_5_3,
            max_5_bp,
            max_5_3,
            decision_dist,
        ]

        self.file.write("\t".join(fields) + "\n")
        self.items_written += 1

    def write_introns(
        self,
        introns: Iterable[Intron],
        species_name: Optional[str] = None,
        simple_name: bool = False,
    ) -> int:
        """
        Write multiple introns.

        Args:
            introns: Iterable of Intron objects
            species_name: Species name
            simple_name: Use simple naming

        Returns:
            Number of introns written
        """
        count = 0
        for intron in introns:
            self.write_intron(intron, species_name, simple_name)
            count += 1
        return count


# ============================================================================
# Mapping File Writers
# ============================================================================


class MappingWriter(BaseWriter):
    """
    Writer for mapping files (duplicate and overlap maps).

    Format (tab-delimited):
        representative_name  duplicate/overlap_name

    These files map excluded introns back to their representative
    for downstream analysis.

    Examples:
        >>> writer = MappingWriter(Path("output.dupe_map.iic"))
        >>> # with writer:
        >>> #     writer.write_mapping("rep1", "dup1")
        >>> #     writer.write_mapping("rep1", "dup2")
    """

    @property
    def mappings_written(self) -> int:
        """Backward-compatible alias for items_written."""
        return self.items_written

    def write_mapping(self, representative: str, excluded: str) -> None:
        """
        Write a single mapping entry.

        Args:
            representative: Name of the representative intron
            excluded: Name of the excluded (duplicate/overlapping) intron
        """
        self._require_open()
        self.file.write(f"{representative}\t{excluded}\n")
        self.items_written += 1

    def write_mappings(self, mappings: Dict[str, Iterable[str]]) -> int:
        """
        Write multiple mappings.

        Args:
            mappings: Dictionary mapping representative -> set of excluded names

        Returns:
            Number of mappings written
        """
        count = 0
        for representative, excluded_set in mappings.items():
            for excluded in excluded_set:
                self.write_mapping(representative, excluded)
                count += 1
        return count


# ============================================================================
# Streaming Output Writer (for true streaming classification)
# ============================================================================


class StreamingOutputWriter:
    """
    Unified output writer for streaming classification.

    Manages all output files (BED, meta, sequences, scores) in a single
    context manager, allowing per-intron writes without accumulating
    introns in memory.

    Used in true streaming mode where introns are processed per-contig
    and written immediately to minimize memory usage.

    Example:
        >>> from pathlib import Path
        >>> config = StreamingOutputConfig(
        ...     output_dir=Path("output"),
        ...     base_name="test",
        ...     species_name="homo_sapiens",
        ... )
        >>> with StreamingOutputWriter(config) as writer:
        ...     for intron in classify_introns_streaming(scored, ensemble):
        ...         writer.write_intron(intron)
        ...     writer.write_summary()
    """

    def __init__(
        self,
        output_dir: Path,
        base_name: str,
        species_name: Optional[str] = None,
        simple_name: bool = False,
        no_abbreviate: bool = False,
        write_bed: bool = True,
        write_sequences: bool = True,
        write_scores: bool = True,
        no_headers: bool = False,
    ):
        """
        Initialize streaming output writer.

        Args:
            output_dir: Directory for output files
            base_name: Base name for output files
            species_name: Species name for intron naming
            simple_name: Use simple naming format
            no_abbreviate: Use full species name
            write_bed: Write BED format output
            write_sequences: Write sequence format output
            write_scores: Write score details output
            no_headers: Omit column headers from output files
        """
        self.output_dir = Path(output_dir)
        self.base_name = base_name
        self.species_name = species_name
        self.simple_name = simple_name
        self.no_abbreviate = no_abbreviate
        self.write_bed = write_bed
        self.write_sequences = write_sequences
        self.write_scores = write_scores
        self.no_headers = no_headers

        # Writers (initialized on __enter__)
        self._bed_writer: Optional[BEDWriter] = None
        self._meta_writer: Optional[MetaWriter] = None
        self._seq_writer: Optional[SequenceWriter] = None
        self._score_writer: Optional[ScoreWriter] = None

        # Counters for summary
        self.total_written = 0
        self.u12_count = 0
        self.u2_count = 0
        self.high_confidence_u12 = 0
        self.threshold = 90.0  # Default threshold for high-confidence

    def __enter__(self) -> "StreamingOutputWriter":
        """Open all output files."""
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Always write metadata
        meta_path = self.output_dir / f"{self.base_name}.meta.iic"
        self._meta_writer = MetaWriter(meta_path)
        self._meta_writer.open()
        if not self.no_headers:
            self._meta_writer.write_header()

        # Optionally write BED
        if self.write_bed:
            bed_path = self.output_dir / f"{self.base_name}.bed.iic"
            self._bed_writer = BEDWriter(bed_path)
            self._bed_writer.open()

        # Optionally write sequences
        if self.write_sequences:
            seq_path = self.output_dir / f"{self.base_name}.introns.iic"
            self._seq_writer = SequenceWriter(seq_path)
            self._seq_writer.open()

        # Optionally write scores
        if self.write_scores:
            score_path = self.output_dir / f"{self.base_name}.score_info.iic"
            self._score_writer = ScoreWriter(score_path)
            self._score_writer.open()
            if not self.no_headers:
                self._score_writer.write_header()

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Close all output files."""
        if self._meta_writer:
            self._meta_writer.close()
        if self._bed_writer:
            self._bed_writer.close()
        if self._seq_writer:
            self._seq_writer.close()
        if self._score_writer:
            self._score_writer.close()

    def write_intron(self, intron: Intron, intron_number: Optional[int] = None) -> None:
        """
        Write a single intron to all output files.

        Args:
            intron: Classified intron to write
            intron_number: Optional simple numbering for output
        """
        # Write to metadata (always)
        if self._meta_writer:
            self._meta_writer.write_intron(
                intron,
                species_name=self.species_name,
                simple_name=self.simple_name,
                no_abbreviate=self.no_abbreviate,
                intron_number=intron_number,
            )

        # Write to BED (if enabled)
        if self._bed_writer:
            self._bed_writer.write_intron(
                intron,
                species_name=self.species_name,
                simple_name=self.simple_name,
                no_abbreviate=self.no_abbreviate,
                intron_number=intron_number,
            )

        # Write to sequences (if enabled)
        if self._seq_writer:
            self._seq_writer.write_intron(
                intron,
                species_name=self.species_name,
                simple_name=self.simple_name,
                no_abbreviate=self.no_abbreviate,
                intron_number=intron_number,
            )

        # Write to scores (if enabled)
        if self._score_writer:
            self._score_writer.write_intron(
                intron,
                species_name=self.species_name,
                simple_name=self.simple_name,
                no_abbreviate=self.no_abbreviate,
            )

        # Update counters
        self.total_written += 1
        if intron.metadata and intron.metadata.type_id == "u12":
            self.u12_count += 1
            if intron.scores and intron.scores.svm_score is not None:
                if intron.scores.svm_score >= self.threshold:
                    self.high_confidence_u12 += 1
        else:
            self.u2_count += 1

    def get_summary(self) -> dict:
        """
        Get classification summary statistics.

        Returns:
            Dictionary with classification counts and percentages
        """
        u12_pct = (
            (self.u12_count / self.total_written * 100)
            if self.total_written > 0
            else 0.0
        )
        high_conf_pct = (
            (self.high_confidence_u12 / self.total_written * 100)
            if self.total_written > 0
            else 0.0
        )

        return {
            "total_introns": self.total_written,
            "u12_count": self.u12_count,
            "u2_count": self.u2_count,
            "u12_percentage": u12_pct,
            "high_confidence_u12": self.high_confidence_u12,
            "high_confidence_percentage": high_conf_pct,
            "threshold": self.threshold,
        }


if __name__ == "__main__":
    import doctest

    doctest.testmod()
