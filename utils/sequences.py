"""
Sequence utility functions for DNA/RNA manipulation.

This module provides functions for common sequence operations:
- Reverse complement
- Sequence validation
- GC content calculation
- Ambiguity checking

Author: intronIC refactoring project
Date: 2025-11-02
"""

from typing import Dict


# Reverse complement lookup table
COMPLEMENT_TABLE: Dict[str, str] = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',
    'N': 'N',
    'a': 't',
    't': 'a',
    'c': 'g',
    'g': 'c',
    'n': 'n',
    # Ambiguity codes (IUPAC)
    'R': 'Y',  # A or G -> T or C
    'Y': 'R',  # C or T -> G or A
    'S': 'S',  # G or C
    'W': 'W',  # A or T
    'K': 'M',  # G or T -> C or A
    'M': 'K',  # A or C -> T or G
    'B': 'V',  # C or G or T -> G or C or A
    'D': 'H',  # A or G or T -> T or C or A
    'H': 'D',  # A or C or T -> T or G or A
    'V': 'B',  # A or C or G -> T or G or C
}


def reverse_complement(seq: str, validate: bool = False) -> str:
    """
    Calculate the reverse complement of a DNA sequence.

    Args:
        seq: DNA sequence (case-insensitive)
        validate: If True, validate sequence before converting (default: False)

    Returns:
        Reverse complement sequence

    Raises:
        ValueError: If validate=True and sequence contains invalid characters

    Examples:
        >>> reverse_complement("ATCG")
        'CGAT'

        >>> reverse_complement("GTAAGT")
        'ACTTAC'

        >>> reverse_complement("atcg")
        'cgat'

        >>> reverse_complement("GTNNAG")
        'CTNNAC'

        >>> reverse_complement("invalid", validate=True)
        Traceback (most recent call last):
            ...
        ValueError: Invalid DNA base(s) in sequence: invalid characters found
    """
    if validate and not is_valid_dna(seq, allow_ambiguous=True):
        raise ValueError(
            f"Invalid DNA base(s) in sequence: invalid characters found"
        )

    # Build reverse complement
    try:
        return ''.join(COMPLEMENT_TABLE[base] for base in reversed(seq))
    except KeyError as e:
        raise ValueError(
            f"Invalid DNA base in sequence: {e.args[0]}"
        ) from e


def is_valid_dna(seq: str, allow_ambiguous: bool = False) -> bool:
    """
    Check if a sequence contains only valid DNA bases.

    Args:
        seq: DNA sequence to validate
        allow_ambiguous: If True, allow IUPAC ambiguity codes (default: False)

    Returns:
        True if sequence is valid, False otherwise

    Examples:
        >>> is_valid_dna("ATCG")
        True

        >>> is_valid_dna("ATCGN", allow_ambiguous=True)
        True

        >>> is_valid_dna("ATCGN", allow_ambiguous=False)
        False

        >>> is_valid_dna("ATCGX")
        False

        >>> is_valid_dna("")
        True

        >>> is_valid_dna("atcg")  # Case-insensitive
        True
    """
    if not seq:
        return True

    if allow_ambiguous:
        valid_bases = set('ATCGNRYSWKMBDHVatcgn ryswkmbdhv')
    else:
        valid_bases = set('ATCGatcg')

    return all(base in valid_bases for base in seq)


def has_ambiguous_bases(seq: str) -> bool:
    """
    Check if a sequence contains ambiguous bases (non-ATCG).

    Args:
        seq: DNA sequence

    Returns:
        True if sequence contains N or other ambiguous bases

    Examples:
        >>> has_ambiguous_bases("ATCG")
        False

        >>> has_ambiguous_bases("ATNGC")
        True

        >>> has_ambiguous_bases("ATWGC")
        True

        >>> has_ambiguous_bases("")
        False
    """
    if not seq:
        return False

    canonical = set('ATCGatcg')
    return any(base not in canonical for base in seq)


def gc_content(seq: str) -> float:
    """
    Calculate GC content as a percentage.

    Args:
        seq: DNA sequence

    Returns:
        GC content as percentage (0.0 - 100.0)

    Raises:
        ValueError: If sequence is empty

    Examples:
        >>> gc_content("ATCG")
        50.0

        >>> gc_content("AAAA")
        0.0

        >>> gc_content("GGCC")
        100.0

        >>> gc_content("ATCGATCG")
        50.0

        >>> gc_content("atcg")  # Case-insensitive
        50.0
    """
    if not seq:
        raise ValueError("Cannot calculate GC content of empty sequence")

    seq_upper = seq.upper()
    gc_count = seq_upper.count('G') + seq_upper.count('C')

    return (gc_count / len(seq)) * 100.0


def count_bases(seq: str) -> Dict[str, int]:
    """
    Count occurrences of each base in a sequence.

    Args:
        seq: DNA sequence

    Returns:
        Dictionary mapping base to count

    Examples:
        >>> count_bases("ATCG")
        {'A': 1, 'T': 1, 'C': 1, 'G': 1, 'N': 0}

        >>> count_bases("AAATTTCCCGGG")
        {'A': 3, 'T': 3, 'C': 3, 'G': 3, 'N': 0}

        >>> count_bases("ATNGC")
        {'A': 1, 'T': 1, 'C': 1, 'G': 1, 'N': 1}
    """
    seq_upper = seq.upper()
    return {
        'A': seq_upper.count('A'),
        'T': seq_upper.count('T'),
        'C': seq_upper.count('C'),
        'G': seq_upper.count('G'),
        'N': seq_upper.count('N'),
    }


def extract_subsequence(seq: str, start: int, stop: int, strand: str = '+') -> str:
    """
    Extract a subsequence with optional reverse complement.

    Args:
        seq: Full sequence
        start: Start position (0-based, inclusive)
        stop: Stop position (0-based, exclusive)
        strand: '+' for forward, '-' for reverse complement

    Returns:
        Extracted subsequence (reverse complemented if strand is '-')

    Raises:
        ValueError: If coordinates are invalid

    Examples:
        >>> extract_subsequence("ATCGATCG", 0, 4)
        'ATCG'

        >>> extract_subsequence("ATCGATCG", 4, 8)
        'ATCG'

        >>> extract_subsequence("ATCGATCG", 0, 4, strand='-')
        'CGAT'

        >>> extract_subsequence("ATCGATCG", -1, 4)
        Traceback (most recent call last):
            ...
        ValueError: Invalid coordinates: start=-1, stop=4 for sequence length 8
    """
    if start < 0 or stop > len(seq) or start >= stop:
        raise ValueError(
            f"Invalid coordinates: start={start}, stop={stop} for sequence length {len(seq)}"
        )

    subseq = seq[start:stop]

    if strand == '-':
        return reverse_complement(subseq)
    elif strand == '+':
        return subseq
    else:
        raise ValueError(f"Invalid strand: {strand} (must be '+' or '-')")


def normalize_sequence(seq: str) -> str:
    """
    Normalize a sequence to uppercase.

    Args:
        seq: DNA sequence

    Returns:
        Uppercase sequence

    Examples:
        >>> normalize_sequence("atcg")
        'ATCG'

        >>> normalize_sequence("AtCg")
        'ATCG'

        >>> normalize_sequence("ATCG")
        'ATCG'
    """
    return seq.upper()


def sliding_window(seq: str, window_size: int, step: int = 1):
    """
    Generate sliding windows over a sequence.

    Args:
        seq: DNA sequence
        window_size: Size of each window
        step: Step size between windows (default: 1)

    Yields:
        Tuples of (start_position, window_sequence)

    Examples:
        >>> list(sliding_window("ATCGATCG", 4))
        [(0, 'ATCG'), (1, 'TCGA'), (2, 'CGAT'), (3, 'GATC'), (4, 'ATCG')]

        >>> list(sliding_window("ATCGATCG", 4, step=2))
        [(0, 'ATCG'), (2, 'CGAT'), (4, 'ATCG')]

        >>> list(sliding_window("ATCG", 10))
        []
    """
    if window_size > len(seq):
        return

    for i in range(0, len(seq) - window_size + 1, step):
        yield (i, seq[i:i + window_size])


if __name__ == "__main__":
    import doctest
    doctest.testmod()
