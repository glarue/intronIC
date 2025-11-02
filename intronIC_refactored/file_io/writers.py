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
from typing import Union, Iterable, Optional, TextIO, Dict, Tuple
from core.intron import Intron


# ============================================================================
# Attribute Tag Mapping
# ============================================================================

# Maps tag codes to verbose attribute names
# Used to generate both compact tags (e.g., [n][i]) and verbose attributes column
TAG_TO_ATTRIBUTE: Dict[str, str] = {
    # Boolean flags (appear as tags in intron name)
    'n': 'noncanonical',          # [n] - Non-canonical splice sites (not GT-AG, GC-AG, AT-AC)
    'i': 'not_longest_isoform',   # [i] - Not from the longest transcript isoform
    'c': 'corrected',             # [c] - Splice site boundaries were adjusted
    'd': 'duplicate',             # [d] - Duplicate coordinates (excluded from analysis)

    # Omission codes (appear as [o:X] tags in intron name)
    # These are INDEPENDENT of property tags - an intron can have both [n] and [o:s]
    'o:s': 'omitted_short',       # [o:s] - Too short (below minimum length threshold)
    'o:a': 'omitted_ambiguous',   # [o:a] - Contains ambiguous bases (N, etc.)
    'o:n': 'omitted_noncanonical',  # [o:n] - Omitted because non-canonical
    'o:v': 'omitted_overlap',     # [o:v] - Overlapping coordinates with another intron
    'o:i': 'omitted_not_longest_isoform',  # [o:i] - Omitted because not in longest isoform
}


def generate_attributes(intron: Intron) -> str:
    """
    Generate verbose comma-separated attributes string from intron metadata.

    Uses TAG_TO_ATTRIBUTE mapping to convert compact tags to
    human-readable attribute names.

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
        if intron.metadata.noncanonical:
            attrs.append(TAG_TO_ATTRIBUTE['n'])
        if not intron.metadata.longest_isoform:
            attrs.append(TAG_TO_ATTRIBUTE['i'])
        if intron.metadata.corrected:
            attrs.append(TAG_TO_ATTRIBUTE['c'])
        if intron.metadata.duplicate:
            attrs.append(TAG_TO_ATTRIBUTE['d'])
        if intron.metadata.omitted:
            # Map omission code to verbose name using o: prefix
            omit_code = f"o:{intron.metadata.omitted}"
            if omit_code in TAG_TO_ATTRIBUTE:
                attrs.append(TAG_TO_ATTRIBUTE[omit_code])
            else:
                # Fallback for unknown codes
                attrs.append(f"omitted_{intron.metadata.omitted}")

    return ','.join(attrs) if attrs else '.'


# ============================================================================
# BED Format Writer
# ============================================================================

class BEDWriter:
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

    def __init__(self, file_path: Union[str, Path]):
        """
        Initialize BED writer.

        Args:
            file_path: Path to output file
        """
        self.file_path = Path(file_path)
        self.file: Optional[TextIO] = None
        self.introns_written = 0

    def open(self) -> None:
        """Open output file for writing."""
        self.file = open(self.file_path, 'w')

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

    def write_header(self) -> None:
        """
        Write BED header (optional).

        BED format doesn't have a standard header, but we can add a track line.
        """
        if not self.file:
            raise ValueError("File not open. Call open() first or use context manager.")

    def write_intron(
        self,
        intron: Intron,
        species_name: Optional[str] = None,
        simple_name: bool = False
    ) -> None:
        """
        Write a single intron in BED format.

        Args:
            intron: Intron object to write
            species_name: Species name for intron label (optional)
            simple_name: Use simple naming (no species prefix)

        Format:
            chrom  start(0-based)  stop  label  svm_score  strand  attributes
        """
        if not self.file:
            raise ValueError("File not open. Call open() first or use context manager.")

        # Get BED start (0-based)
        start_0based = intron.start - 1

        # Get SVM score or '.' if unavailable
        score = '.' if intron.svm_score is None else str(intron.svm_score)

        # Generate intron label
        label = self._generate_label(intron, species_name, simple_name)

        # Generate verbose attributes
        attributes = generate_attributes(intron)

        # Write BED line
        fields = [
            intron.chromosome,
            str(start_0based),
            str(intron.stop),
            label,
            score,
            intron.strand,
            attributes
        ]
        self.file.write('\t'.join(fields) + '\n')
        self.introns_written += 1

    def write_introns(
        self,
        introns: Iterable[Intron],
        species_name: Optional[str] = None,
        simple_name: bool = False
    ) -> int:
        """
        Write multiple introns.

        Args:
            introns: Iterable of Intron objects
            species_name: Species name for labels
            simple_name: Use simple naming

        Returns:
            Number of introns written
        """
        count = 0
        for intron in introns:
            self.write_intron(intron, species_name, simple_name)
            count += 1
        return count

    def _generate_label(
        self,
        intron: Intron,
        species_name: Optional[str],
        simple_name: bool
    ) -> str:
        """
        Generate intron label for BED name field.

        Format: [species_]parent_index;svm_score[tags]

        Args:
            intron: Intron object
            species_name: Species name (optional)
            simple_name: Exclude species prefix

        Returns:
            Intron label string
        """
        parts = []

        # Add species prefix if requested
        if species_name and not simple_name:
            parts.append(species_name)

        # Add parent_index
        if intron.metadata and intron.metadata.parent:
            parent = intron.metadata.parent
            index = intron.metadata.index if intron.metadata.index else 1
            parts.append(f"{parent}_{index}")
        else:
            parts.append(intron.intron_id)

        name = '_'.join(parts) if parts else intron.intron_id

        # Add SVM score if available
        if intron.svm_score is not None:
            name += f";{intron.svm_score:.2f}"

        # Add tags
        tags = self._generate_tags(intron)
        if tags:
            name += tags

        return name

    def _generate_tags(self, intron: Intron) -> str:
        """
        Generate tag string for intron.

        Tags indicate special properties:
        - [n] = non-canonical
        - [i] = not longest isoform
        - [c] = corrected (or [c:N] if distance available)
        - [d] = duplicate

        Omission tags (independent, can appear with property tags):
        - [o:s] = omitted:short
        - [o:a] = omitted:ambiguous
        - [o:n] = omitted:noncanonical
        - [o:v] = omitted:overlap
        - [o:i] = omitted:not_longest_isoform

        Args:
            intron: Intron object

        Returns:
            Tag string (e.g., "[n][i]" or "[o:s]" or "[n][o:s]")
        """
        tags = []

        if intron.metadata:
            if intron.metadata.noncanonical:
                tags.append("[n]")
            if not intron.metadata.longest_isoform:
                tags.append("[i]")
            if intron.metadata.corrected:
                tags.append("[c]")  # Could add distance if available
            if intron.metadata.duplicate:
                tags.append("[d]")
            if intron.metadata.omitted:
                # Omission tags use o: prefix (e.g., [o:s] for omitted:short)
                tags.append(f"[o:{intron.metadata.omitted}]")

        return ''.join(tags)


# ============================================================================
# Metadata Format Writer
# ============================================================================

class MetaWriter:
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

    def __init__(self, file_path: Union[str, Path]):
        """
        Initialize metadata writer.

        Args:
            file_path: Path to output file
        """
        self.file_path = Path(file_path)
        self.file: Optional[TextIO] = None
        self.introns_written = 0

    def open(self) -> None:
        """Open output file for writing."""
        self.file = open(self.file_path, 'w')

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

    def write_header(self) -> None:
        """Write metadata file header."""
        if not self.file:
            raise ValueError("File not open. Call open() first or use context manager.")

        header_fields = [
            "name", "rel_score", "dnts", "motif_schematic", "bp_context",
            "length", "parent", "grandparent", "index", "family_size",
            "frac_pos", "phase", "type_id", "feature", "attributes"
        ]
        self.file.write('\t'.join(header_fields) + '\n')

    def write_intron(
        self,
        intron: Intron,
        species_name: Optional[str] = None,
        simple_name: bool = False,
        null: str = '.'
    ) -> None:
        """
        Write a single intron's metadata.

        Args:
            intron: Intron object to write
            species_name: Species name for intron name
            simple_name: Use simple naming
            null: Placeholder for missing values

        Format:
            name  rel_score  dnts  motif  bp_context  length  parent  grandparent
            index  family_size  frac_pos  phase  type_id  feature
        """
        if not self.file:
            raise ValueError("File not open. Call open() first or use context manager.")

        # Generate intron name
        name = self._generate_name(intron, species_name, simple_name)

        # Relative score (rounded to 4 decimal places)
        rel_score = null
        if intron.scores and intron.scores.relative_score is not None:
            rel_score = str(round(intron.scores.relative_score, 4))

        # Terminal dinucleotides (e.g., 'GT-AG')
        dnts = null
        if intron.sequences and intron.sequences.terminal_dinucleotides:
            dnts = intron.sequences.terminal_dinucleotides

        # Motif schematic (simplified - would need actual implementation)
        motif = null  # TODO: implement motif_string property

        # Branch point context
        bp_context = null
        if intron.sequences and intron.sequences.bp_seq:
            bp_context = intron.sequences.bp_seq

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
            grandparent = intron.metadata.grandparent if intron.metadata.grandparent else null
            index = str(intron.metadata.index) if intron.metadata.index is not None else null
            family_size = str(intron.metadata.family_size) if intron.metadata.family_size else null
            frac_pos_val = intron.metadata.fractional_position
            frac_pos = str(round(frac_pos_val, 4)) if frac_pos_val is not None else null
            phase = str(intron.metadata.exon_phase) if intron.metadata.exon_phase is not None else null

        # Type ID
        type_id = intron.type_id

        # Feature type (exon/cds)
        feature = null  # TODO: track feature type

        # Generate verbose attributes
        attributes = generate_attributes(intron)

        fields = [
            name, rel_score, dnts, motif, bp_context, length, parent, grandparent,
            index, family_size, frac_pos, phase, type_id, feature, attributes
        ]

        self.file.write('\t'.join(fields) + '\n')
        self.introns_written += 1

    def write_introns(
        self,
        introns: Iterable[Intron],
        species_name: Optional[str] = None,
        simple_name: bool = False
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

    def _generate_name(
        self,
        intron: Intron,
        species_name: Optional[str],
        simple_name: bool
    ) -> str:
        """
        Generate intron name for metadata.

        Format: [species_]parent_index

        Args:
            intron: Intron object
            species_name: Species name (optional)
            simple_name: Exclude species prefix

        Returns:
            Intron name string
        """
        parts = []

        if species_name and not simple_name:
            parts.append(species_name)

        if intron.metadata and intron.metadata.parent:
            parent = intron.metadata.parent
            index = intron.metadata.index if intron.metadata.index else 1
            parts.append(f"{parent}_{index}")
        else:
            parts.append(intron.intron_id)

        return '_'.join(parts) if parts else intron.intron_id


# ============================================================================
# Sequence Format Writer
# ============================================================================

class SequenceWriter:
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

    def __init__(self, file_path: Union[str, Path]):
        """
        Initialize sequence writer.

        Args:
            file_path: Path to output file
        """
        self.file_path = Path(file_path)
        self.file: Optional[TextIO] = None
        self.introns_written = 0

    def open(self) -> None:
        """Open output file for writing."""
        self.file = open(self.file_path, 'w')

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

    def write_intron(
        self,
        intron: Intron,
        species_name: Optional[str] = None,
        simple_name: bool = False,
        include_score: bool = True
    ) -> None:
        """
        Write a single intron's sequences.

        Args:
            intron: Intron object to write
            species_name: Species name for intron name
            simple_name: Use simple naming
            include_score: Include SVM score in output

        Format:
            name  [score]  upstream_flank  sequence  downstream_flank
        """
        if not self.file:
            raise ValueError("File not open. Call open() first or use context manager.")

        if not intron.sequences or not intron.sequences.seq:
            raise ValueError(f"Intron {intron.intron_id} has no sequence data")

        # Generate intron name
        name = self._generate_name(intron, species_name, simple_name)

        # Get sequences (with defaults)
        upstream = intron.sequences.upstream_flank or ""
        sequence = intron.sequences.seq
        downstream = intron.sequences.downstream_flank or ""

        fields = [name]

        # Optionally include score
        if include_score:
            score = str(intron.svm_score) if intron.svm_score is not None else "."
            fields.append(score)

        fields.extend([upstream, sequence, downstream])

        self.file.write('\t'.join(fields) + '\n')
        self.introns_written += 1

    def write_introns(
        self,
        introns: Iterable[Intron],
        species_name: Optional[str] = None,
        simple_name: bool = False,
        include_score: bool = True
    ) -> int:
        """
        Write multiple introns.

        Args:
            introns: Iterable of Intron objects
            species_name: Species name
            simple_name: Use simple naming
            include_score: Include SVM score

        Returns:
            Number of introns written
        """
        count = 0
        for intron in introns:
            self.write_intron(intron, species_name, simple_name, include_score)
            count += 1
        return count

    def _generate_name(
        self,
        intron: Intron,
        species_name: Optional[str],
        simple_name: bool
    ) -> str:
        """Generate intron name (same as MetaWriter)."""
        parts = []

        if species_name and not simple_name:
            parts.append(species_name)

        if intron.metadata and intron.metadata.parent:
            parent = intron.metadata.parent
            index = intron.metadata.index if intron.metadata.index else 1
            parts.append(f"{parent}_{index}")
        else:
            parts.append(intron.intron_id)

        return '_'.join(parts) if parts else intron.intron_id


# ============================================================================
# Score Details Writer
# ============================================================================

class ScoreWriter:
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

    def __init__(self, file_path: Union[str, Path]):
        """
        Initialize score writer.

        Args:
            file_path: Path to output file
        """
        self.file_path = Path(file_path)
        self.file: Optional[TextIO] = None
        self.introns_written = 0

    def open(self) -> None:
        """Open output file for writing."""
        self.file = open(self.file_path, 'w')

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

    def write_header(self) -> None:
        """Write score file header."""
        if not self.file:
            raise ValueError("File not open. Call open() first or use context manager.")

        header_fields = [
            "name", "rel_score", "svm_score", "decision_dist",
            "5'_seq", "5'_raw", "5'_z",
            "bp_seq", "bp_region", "bp_raw", "bp_z",
            "3'_seq", "3'_raw", "3'_z"
        ]
        self.file.write('\t'.join(header_fields) + '\n')

    def write_intron(
        self,
        intron: Intron,
        species_name: Optional[str] = None,
        simple_name: bool = False,
        null: str = '.'
    ) -> None:
        """
        Write a single intron's detailed scores.

        Args:
            intron: Intron object to write
            species_name: Species name for intron name
            simple_name: Use simple naming
            null: Placeholder for missing values
        """
        if not self.file:
            raise ValueError("File not open. Call open() first or use context manager.")

        # Generate intron name
        name = self._generate_name(intron, species_name, simple_name)

        # Default all values to null
        rel_score = null
        svm_score = null
        decision_dist = null
        five_seq = null
        five_raw = null
        five_z = null
        bp_seq = null
        bp_region = null
        bp_raw = null
        bp_z = null
        three_seq = null
        three_raw = null
        three_z = null

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

        # Fill in sequences if available
        if intron.sequences:
            if intron.sequences.five_seq:
                five_seq = intron.sequences.five_seq
            if intron.sequences.bp_seq:
                bp_seq = intron.sequences.bp_seq
            if intron.sequences.bp_region_seq:
                bp_region = intron.sequences.bp_region_seq
            if intron.sequences.three_seq:
                three_seq = intron.sequences.three_seq

        fields = [
            name, rel_score, svm_score, decision_dist,
            five_seq, five_raw, five_z,
            bp_seq, bp_region, bp_raw, bp_z,
            three_seq, three_raw, three_z
        ]

        self.file.write('\t'.join(fields) + '\n')
        self.introns_written += 1

    def write_introns(
        self,
        introns: Iterable[Intron],
        species_name: Optional[str] = None,
        simple_name: bool = False
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

    def _generate_name(
        self,
        intron: Intron,
        species_name: Optional[str],
        simple_name: bool
    ) -> str:
        """Generate intron name (same as MetaWriter)."""
        parts = []

        if species_name and not simple_name:
            parts.append(species_name)

        if intron.metadata and intron.metadata.parent:
            parent = intron.metadata.parent
            index = intron.metadata.index if intron.metadata.index else 1
            parts.append(f"{parent}_{index}")
        else:
            parts.append(intron.intron_id)

        return '_'.join(parts) if parts else intron.intron_id


# ============================================================================
# Mapping File Writers
# ============================================================================

class MappingWriter:
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

    def __init__(self, file_path: Union[str, Path]):
        """
        Initialize mapping writer.

        Args:
            file_path: Path to output file
        """
        self.file_path = Path(file_path)
        self.file: Optional[TextIO] = None
        self.mappings_written = 0

    def open(self) -> None:
        """Open output file for writing."""
        self.file = open(self.file_path, 'w')

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

    def write_mapping(self, representative: str, excluded: str) -> None:
        """
        Write a single mapping entry.

        Args:
            representative: Name of the representative intron
            excluded: Name of the excluded (duplicate/overlapping) intron
        """
        if not self.file:
            raise ValueError("File not open. Call open() first or use context manager.")

        self.file.write(f"{representative}\t{excluded}\n")
        self.mappings_written += 1

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


if __name__ == "__main__":
    import doctest
    doctest.testmod()
