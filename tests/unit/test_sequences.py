"""
Unit tests for sequence utility functions.

Tests:
- reverse_complement: DNA reverse complement
- is_valid_dna: Sequence validation
- has_ambiguous_bases: Ambiguity detection
- gc_content: GC percentage calculation
- count_bases: Base counting
- extract_subsequence: Subsequence extraction with strand
- normalize_sequence: Case normalization
- sliding_window: Window generation

Author: intronIC refactoring project
Date: 2025-11-02
"""

import pytest
from intronIC.utils.sequences import (
    reverse_complement,
    is_valid_dna,
    has_ambiguous_bases,
    gc_content,
    count_bases,
    extract_subsequence,
    normalize_sequence,
    sliding_window,
)


class TestReverseComplement:
    """Test reverse_complement function."""

    def test_simple_sequences(self):
        """Test reverse complement of simple sequences."""
        assert reverse_complement("A") == "T"
        assert reverse_complement("T") == "A"
        assert reverse_complement("C") == "G"
        assert reverse_complement("G") == "C"

    def test_longer_sequences(self):
        """Test reverse complement of longer sequences."""
        assert reverse_complement("ATCG") == "CGAT"
        assert reverse_complement("GTAAGT") == "ACTTAC"
        assert reverse_complement("AAAATTTT") == "AAAATTTT"

    def test_case_insensitive(self):
        """Test that function handles lowercase."""
        assert reverse_complement("atcg") == "cgat"
        assert reverse_complement("AtCg") == "cGaT"

    def test_with_n_bases(self):
        """Test reverse complement with N bases."""
        assert reverse_complement("ATNGC") == "GCNAT"
        assert reverse_complement("NNNNN") == "NNNNN"

    def test_with_ambiguous_bases(self):
        """Test reverse complement with IUPAC ambiguity codes."""
        assert reverse_complement("ATRYC") == "GRYAT"  # R->Y, Y->R
        assert reverse_complement("WSMK") == "MKSW"    # Reverse: KMSW, then K->M, M->K, S->S, W->W

    def test_empty_sequence(self):
        """Test reverse complement of empty sequence."""
        assert reverse_complement("") == ""

    def test_with_validation(self):
        """Test reverse complement with validation enabled."""
        # Valid sequence should work
        assert reverse_complement("ATCG", validate=True) == "CGAT"

        # Invalid sequence should raise error
        with pytest.raises(ValueError, match="Invalid DNA base"):
            reverse_complement("ATCGX", validate=True)

    def test_invalid_base_without_validation(self):
        """Test that invalid bases raise error even without validation."""
        with pytest.raises(ValueError, match="Invalid DNA base"):
            reverse_complement("ATCGX", validate=False)


class TestIsValidDNA:
    """Test is_valid_dna function."""

    def test_valid_sequences(self):
        """Test valid DNA sequences."""
        assert is_valid_dna("ATCG") == True
        assert is_valid_dna("AAAA") == True
        assert is_valid_dna("TTTTCCCCGGGG") == True

    def test_case_insensitive(self):
        """Test case insensitivity."""
        assert is_valid_dna("atcg") == True
        assert is_valid_dna("AtCg") == True

    def test_empty_sequence(self):
        """Test empty sequence is considered valid."""
        assert is_valid_dna("") == True

    def test_with_n_no_ambiguous(self):
        """Test N is invalid when ambiguous not allowed."""
        assert is_valid_dna("ATCGN", allow_ambiguous=False) == False

    def test_with_n_allow_ambiguous(self):
        """Test N is valid when ambiguous allowed."""
        assert is_valid_dna("ATCGN", allow_ambiguous=True) == True

    def test_with_other_ambiguous_codes(self):
        """Test IUPAC ambiguity codes."""
        assert is_valid_dna("ATCGRYSWKM", allow_ambiguous=True) == True
        assert is_valid_dna("ATCGRYSWKM", allow_ambiguous=False) == False

    def test_invalid_characters(self):
        """Test sequences with invalid characters."""
        assert is_valid_dna("ATCGX") == False
        assert is_valid_dna("ATCG123") == False
        assert is_valid_dna("ATCG-") == False


class TestHasAmbiguousBases:
    """Test has_ambiguous_bases function."""

    def test_no_ambiguous(self):
        """Test sequences with no ambiguous bases."""
        assert has_ambiguous_bases("ATCG") == False
        assert has_ambiguous_bases("AAATTTCCCGGG") == False

    def test_with_n(self):
        """Test sequences with N."""
        assert has_ambiguous_bases("ATNGC") == True
        assert has_ambiguous_bases("NNNNN") == True

    def test_with_other_ambiguous(self):
        """Test sequences with other ambiguity codes."""
        assert has_ambiguous_bases("ATRYC") == True
        assert has_ambiguous_bases("WSMK") == True

    def test_empty_sequence(self):
        """Test empty sequence has no ambiguous bases."""
        assert has_ambiguous_bases("") == False

    def test_case_insensitive(self):
        """Test case insensitivity."""
        assert has_ambiguous_bases("atngc") == True
        assert has_ambiguous_bases("atcg") == False


class TestGCContent:
    """Test gc_content function."""

    def test_balanced_gc(self):
        """Test sequence with 50% GC."""
        assert gc_content("ATCG") == 50.0
        assert gc_content("ATCGATCG") == 50.0

    def test_all_gc(self):
        """Test sequence with 100% GC."""
        assert gc_content("GGCC") == 100.0
        assert gc_content("CCCCGGGG") == 100.0

    def test_no_gc(self):
        """Test sequence with 0% GC."""
        assert gc_content("AAAA") == 0.0
        assert gc_content("TTTTAAAA") == 0.0

    def test_with_n(self):
        """Test GC content with N bases (counted in denominator)."""
        # ATCGN = 2 GC out of 5 bases = 40%
        assert gc_content("ATCGN") == 40.0

    def test_case_insensitive(self):
        """Test case insensitivity."""
        assert gc_content("atcg") == 50.0
        assert gc_content("AtCg") == 50.0

    def test_empty_sequence_raises(self):
        """Test that empty sequence raises error."""
        with pytest.raises(ValueError, match="empty sequence"):
            gc_content("")

    def test_typical_values(self):
        """Test with typical genomic sequences."""
        # Human genome average ~41%
        assert 40.0 < gc_content("ATGCATGCATGCATGC") < 60.0


class TestCountBases:
    """Test count_bases function."""

    def test_balanced_sequence(self):
        """Test sequence with equal bases."""
        counts = count_bases("ATCG")
        assert counts['A'] == 1
        assert counts['T'] == 1
        assert counts['C'] == 1
        assert counts['G'] == 1
        assert counts['N'] == 0

    def test_homopolymer(self):
        """Test homopolymer sequences."""
        counts = count_bases("AAAA")
        assert counts['A'] == 4
        assert counts['T'] == 0
        assert counts['C'] == 0
        assert counts['G'] == 0

    def test_with_n(self):
        """Test counting N bases."""
        counts = count_bases("ATNGC")
        assert counts['A'] == 1
        assert counts['T'] == 1
        assert counts['N'] == 1
        assert counts['G'] == 1
        assert counts['C'] == 1

    def test_case_insensitive(self):
        """Test case insensitivity."""
        counts = count_bases("atcg")
        assert counts['A'] == 1
        assert counts['T'] == 1
        assert counts['C'] == 1
        assert counts['G'] == 1

    def test_empty_sequence(self):
        """Test counting bases in empty sequence."""
        counts = count_bases("")
        assert counts['A'] == 0
        assert counts['T'] == 0
        assert counts['C'] == 0
        assert counts['G'] == 0
        assert counts['N'] == 0


class TestExtractSubsequence:
    """Test extract_subsequence function."""

    def test_forward_strand(self):
        """Test extraction on forward strand."""
        seq = "ATCGATCG"
        assert extract_subsequence(seq, 0, 4, '+') == "ATCG"
        assert extract_subsequence(seq, 4, 8, '+') == "ATCG"
        assert extract_subsequence(seq, 2, 6, '+') == "CGAT"

    def test_reverse_strand(self):
        """Test extraction with reverse complement on reverse strand."""
        seq = "ATCGATCG"
        # Extract ATCG (0:4), reverse complement = CGAT
        assert extract_subsequence(seq, 0, 4, '-') == "CGAT"

    def test_full_sequence(self):
        """Test extracting full sequence."""
        seq = "ATCGATCG"
        assert extract_subsequence(seq, 0, 8, '+') == "ATCGATCG"

    def test_single_base(self):
        """Test extracting single base."""
        seq = "ATCGATCG"
        assert extract_subsequence(seq, 0, 1, '+') == "A"
        assert extract_subsequence(seq, 0, 1, '-') == "T"

    def test_invalid_coordinates(self):
        """Test that invalid coordinates raise errors."""
        seq = "ATCGATCG"

        # Negative start
        with pytest.raises(ValueError, match="Invalid coordinates"):
            extract_subsequence(seq, -1, 4, '+')

        # Stop beyond length
        with pytest.raises(ValueError, match="Invalid coordinates"):
            extract_subsequence(seq, 0, 10, '+')

        # Start >= stop
        with pytest.raises(ValueError, match="Invalid coordinates"):
            extract_subsequence(seq, 4, 4, '+')

        with pytest.raises(ValueError, match="Invalid coordinates"):
            extract_subsequence(seq, 6, 4, '+')

    def test_invalid_strand(self):
        """Test that invalid strand raises error."""
        seq = "ATCGATCG"

        with pytest.raises(ValueError, match="Invalid strand"):
            extract_subsequence(seq, 0, 4, '.')


class TestNormalizeSequence:
    """Test normalize_sequence function."""

    def test_lowercase_to_uppercase(self):
        """Test converting lowercase to uppercase."""
        assert normalize_sequence("atcg") == "ATCG"
        assert normalize_sequence("atcgatcg") == "ATCGATCG"

    def test_already_uppercase(self):
        """Test sequence already uppercase."""
        assert normalize_sequence("ATCG") == "ATCG"

    def test_mixed_case(self):
        """Test mixed case sequence."""
        assert normalize_sequence("AtCg") == "ATCG"
        assert normalize_sequence("aTcG") == "ATCG"

    def test_empty_sequence(self):
        """Test empty sequence."""
        assert normalize_sequence("") == ""


class TestSlidingWindow:
    """Test sliding_window function."""

    def test_window_size_4_step_1(self):
        """Test sliding window with default step."""
        windows = list(sliding_window("ATCGATCG", 4))
        expected = [
            (0, "ATCG"),
            (1, "TCGA"),
            (2, "CGAT"),
            (3, "GATC"),
            (4, "ATCG")
        ]
        assert windows == expected

    def test_window_size_4_step_2(self):
        """Test sliding window with step=2."""
        windows = list(sliding_window("ATCGATCG", 4, step=2))
        expected = [
            (0, "ATCG"),
            (2, "CGAT"),
            (4, "ATCG")
        ]
        assert windows == expected

    def test_window_larger_than_sequence(self):
        """Test window larger than sequence returns empty."""
        windows = list(sliding_window("ATCG", 10))
        assert windows == []

    def test_window_equals_sequence(self):
        """Test window equal to sequence length."""
        windows = list(sliding_window("ATCG", 4))
        assert windows == [(0, "ATCG")]

    def test_window_size_1(self):
        """Test sliding window of size 1 (each base)."""
        windows = list(sliding_window("ATCG", 1))
        expected = [
            (0, "A"),
            (1, "T"),
            (2, "C"),
            (3, "G")
        ]
        assert windows == expected

    def test_step_larger_than_window(self):
        """Test step size larger than window."""
        windows = list(sliding_window("ATCGATCGATCG", 3, step=5))
        expected = [
            (0, "ATC"),
            (5, "TCG"),
            (10, "CG")  # Note: This is only 2 bases, wait this shouldn't be included
        ]
        # Actually, let me recalculate. The sequence is 12 bases long.
        # Window size is 3, step is 5.
        # i=0: seq[0:3] = "ATC"
        # i=5: seq[5:8] = "TCG"
        # i=10: seq[10:13] = "CG" - but that's only 2 bases, so it shouldn't be included
        # Wait, let me check the function logic. It goes up to len(seq) - window_size + 1
        # So for len=12, window=3: range(0, 12-3+1, 5) = range(0, 10, 5) = [0, 5]
        # So it should be:
        windows = list(sliding_window("ATCGATCGATCG", 3, step=5))
        expected = [
            (0, "ATC"),
            (5, "TCG")
        ]
        assert windows == expected


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
