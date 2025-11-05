"""
U12 boundary correction for misannotated non-canonical introns.

This module implements the U12 boundary correction algorithm that searches
for misannotated U12-type intron boundaries. U12 introns have a distinctive
5' splice site motif (ATATCCTT or similar), but annotations may miss this
by a few nucleotides, resulting in "non-canonical" boundaries.

The correction:
1. Searches for strong U12 5' SS motif in a 17bp window
2. Calculates the shift needed to align with the motif
3. Adjusts BOTH boundaries symmetrically (maintains intron length change)
4. Updates phase and adds correction tag
5. IMPORTANT: Requires sequence re-extraction after correction

Port from: intronIC.py:2298-2353 (u12_correction)

Algorithm details documented in: U12_CORRECTION_SPEC.md

Author: intronIC refactoring project
Date: 2025-11-05
"""

import re
from typing import Optional, Tuple
from collections import deque
from dataclasses import replace

from core.intron import Intron
from utils.coordinates import GenomicCoordinate


# U12-specific 5' splice site motifs
# Strict motif (default): [AG]TATCC([ACTG]T|T[ACTG])
# Matches patterns like: ATATCCAT, ATATCCTT, GTATCCAT, etc.
STRICT_U12_MOTIF = re.compile(r'[AG]TATCC([ACTG]T|T[ACTG])')

# Lax motif (alternative, not currently used): [AG]TATC[CT]
LAX_U12_MOTIF = re.compile(r'[AG]TATC[CT]')


def search_u12_boundary(
    intron: Intron,
    upstream_flank: str,
    downstream_intron: str,
    use_strict: bool = True
) -> Optional[int]:
    """
    Search for misannotated U12 boundary in non-canonical intron.

    Searches a 17bp window (last 5bp of exon + first 12bp of intron) for
    a strong U12-specific 5' splice site motif. If found at a different
    position than expected, returns the shift needed to correct.

    Port from: intronIC.py:2298-2353 (u12_correction)

    Args:
        intron: Intron to check (should be non-canonical)
        upstream_flank: Upstream exonic sequence (≥5bp needed)
        downstream_intron: Beginning of intron sequence (≥12bp needed)
        use_strict: Use strict motif (default) vs lax motif

    Returns:
        Shift amount (int) if correction needed, None if:
        - No motif found
        - Motif at expected position (shift=0)
        - Insufficient sequence length

    Example:
        >>> # Intron with GC-AG boundaries (non-canonical 5')
        >>> # Search finds ATATCCTT at position 3 (should be at 5)
        >>> shift = search_u12_boundary(intron, "...CGATC", "GCTATCCTTATA...")
        >>> # Returns -2 (move boundaries 2bp upstream)
    """
    # Port from: intronIC.py:2325-2326
    up_n = 5    # Bases of upstream flank to search
    down_n = 12  # Bases into intron to search

    # Validate sequence lengths
    if len(upstream_flank) < up_n or len(downstream_intron) < down_n:
        return None

    # Select motif pattern
    # Port from: intronIC.py:2327-2336
    motif = STRICT_U12_MOTIF if use_strict else LAX_U12_MOTIF

    # Build search region: last 5bp of exon + first 12bp of intron
    # Port from: intronIC.py:2337
    search_region = upstream_flank[-up_n:] + downstream_intron[:down_n]

    # Search for U12 motif
    # Port from: intronIC.py:2338-2340
    match = motif.search(search_region)
    if not match:
        return None  # No U12 motif found

    # Calculate shift needed
    # Port from: intronIC.py:2341-2344
    match_index = match.start()  # Position where motif starts (0-16)
    shift = match_index - up_n    # Shift relative to expected position

    if shift == 0:
        return None  # Motif at expected position, no correction needed

    return shift


def apply_u12_correction(
    intron: Intron,
    shift: int
) -> Intron:
    """
    Apply boundary correction to intron.

    Creates a new intron with adjusted coordinates, rotated phase,
    and correction tag. The shift is applied symmetrically to both
    boundaries (maintaining intron length modulo the shift).

    IMPORTANT: After correction, sequences MUST be re-extracted using
    the new coordinates. This function only updates metadata and coordinates.

    Port from: intronIC.py:2345-2353 (in-place modification)

    Args:
        intron: Original intron
        shift: Number of bases to shift (negative = upstream, positive = downstream)

    Returns:
        New Intron with corrected coordinates and metadata

    Example:
        >>> # Original: chr1:1000-2000:+ (GC-AG)
        >>> # Shift: -2 (move 2bp upstream)
        >>> corrected = apply_u12_correction(intron, -2)
        >>> # Result: chr1:998-1998:+ with [c:-2] tag
    """
    # Rotate phase by shift amount
    # Port from: intronIC.py:2315-2323, 2346
    new_phase = _rotate_phase(intron.metadata.phase if intron.metadata else None, shift)

    # Calculate new coordinates
    # Port from: intronIC.py:2348-2351
    # Shift direction depends on strand
    coord_shift = shift
    if intron.coordinates.strand == '-':
        coord_shift *= -1  # Reverse for negative strand

    new_start = intron.coordinates.start + coord_shift
    new_stop = intron.coordinates.stop + coord_shift

    # Create new coordinate object
    new_coordinates = GenomicCoordinate(
        chromosome=intron.coordinates.chromosome,
        start=new_start,
        stop=new_stop,
        strand=intron.coordinates.strand,
        system='1-based'
    )

    # Update metadata with correction info
    # Port from: intronIC.py:2345, 2347
    updated_metadata = replace(
        intron.metadata,
        corrected=shift,  # Store shift amount
        phase=new_phase
    )

    # Create corrected intron
    # Note: Sequences are NOT updated here - caller must re-extract
    corrected_intron = replace(
        intron,
        coordinates=new_coordinates,
        metadata=updated_metadata,
        sequences=None  # Clear sequences - must be re-extracted
    )

    return corrected_intron


def _rotate_phase(phase: Optional[int], shift: int) -> Optional[int]:
    """
    Rotate coding phase by shift amount.

    When intron boundaries shift, the reading frame changes. This function
    rotates the phase using a circular deque to maintain correct frame.

    Port from: intronIC.py:2315-2323 (_shift_phase helper)

    Args:
        phase: Original phase (0, 1, 2, or None)
        shift: Amount to shift

    Returns:
        New phase after rotation, or None if phase was None

    Examples:
        >>> _rotate_phase(0, -2)  # Shift -2: 0 → 1
        1
        >>> _rotate_phase(1, -2)  # Shift -2: 1 → 2
        2
        >>> _rotate_phase(2, -2)  # Shift -2: 2 → 0
        0
        >>> _rotate_phase(None, -2)  # Non-coding exon
        None
    """
    if phase is None:
        return None

    # Create circular deque of phases
    phases = deque([0, 1, 2])

    try:
        current_index = phases.index(int(phase))
    except (ValueError, TypeError):
        # Invalid phase value, return unchanged
        return phase

    # Rotate by negative shift amount
    # (negative rotation means earlier positions move forward)
    phases.rotate(-shift)

    return phases[current_index]


def correct_intron_if_needed(
    intron: Intron,
    correction_enabled: bool = True,
    use_strict_motif: bool = True
) -> Tuple[Intron, bool]:
    """
    Check intron and apply U12 correction if needed.

    This is the main entry point for U12 boundary correction.

    Only non-canonical introns are checked. If a strong U12 motif is found
    at a shifted position, the intron is corrected and returned with
    sequences cleared (caller must re-extract).

    Port from: intronIC.py:2692 (conditional call to u12_correction)

    Args:
        intron: Intron to check
        correction_enabled: Whether correction is enabled (--no_nc_ss_adjustment flag)
        use_strict_motif: Use strict vs lax U12 motif

    Returns:
        Tuple of (intron, was_corrected)
        - intron: Original or corrected intron
        - was_corrected: True if correction was applied

    Example:
        >>> corrected, changed = correct_intron_if_needed(nc_intron, True)
        >>> if changed:
        ...     # Re-extract sequences with new coordinates
        ...     corrected = extract_sequences(corrected, genome)
    """
    # Skip if correction disabled
    if not correction_enabled:
        return intron, False

    # Skip if not non-canonical
    if not (intron.metadata and intron.metadata.noncanonical):
        return intron, False

    # Skip if no sequences available
    if not intron.sequences or not intron.sequences.seq:
        return intron, False

    # Extract upstream flank and intron start
    upstream_flank = intron.sequences.upstream_flank or ""
    intron_seq = intron.sequences.seq or ""

    # Search for U12 motif
    shift = search_u12_boundary(
        intron,
        upstream_flank,
        intron_seq,
        use_strict=use_strict_motif
    )

    # No correction needed
    if shift is None:
        return intron, False

    # Apply correction
    corrected_intron = apply_u12_correction(intron, shift)

    return corrected_intron, True


if __name__ == "__main__":
    import doctest
    doctest.testmod()
