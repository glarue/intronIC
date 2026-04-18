"""
Tests for U12 boundary correction.

Tests the require_canonical option and verifies that boundary correction
handles both canonical and non-canonical correction results correctly.
"""

import pytest
from intronIC.extraction.boundary_correction import (
    correct_intron_if_needed,
    search_u12_boundary,
    get_corrected_dinucleotides,
)
from intronIC.core.intron import (
    Intron,
    GenomicCoordinate,
    IntronSequences,
    IntronScores,
    IntronMetadata,
)


def _make_intron(upstream, seq, downstream, noncanonical=True):
    """Helper to build a test intron."""
    intron = Intron(
        intron_id="test",
        coordinates=GenomicCoordinate("chr1", 1000, 1000 + len(seq), "+", "1-based"),
        sequences=IntronSequences(
            seq=seq.upper(),
            upstream_flank=upstream.upper(),
            downstream_flank=downstream.upper(),
            five_prime_dnt=seq[:2].upper(),
            three_prime_dnt=seq[-2:].upper(),
        ),
        scores=IntronScores(),
        metadata=IntronMetadata("t", "g"),
    )
    intron.metadata.noncanonical = noncanonical
    return intron


class TestSearchU12Boundary:
    """Test the U12 motif search."""

    def test_finds_shifted_motif(self):
        """Finds ATATCCTT motif at shifted position."""
        # Upstream ends with ...AAGGT, intron starts with TATATCCTCC
        # Motif ATATCC at position 6 in search window → shift = 1
        upstream = "AGAACAAGGT"
        seq = "TATATCCTCCGCC" + "A" * 40 + "AG"
        intron = _make_intron(upstream, seq, "CCCCCCCCCC")
        shift = search_u12_boundary(intron, upstream.upper(), seq.upper())
        assert shift == 1

    def test_no_motif_returns_none(self):
        """Returns None when no U12 motif present."""
        upstream = "AAAAAAAAAA"
        seq = "GTCCCCCCCCCCCCCCCCCCCCCCCCCAG"
        intron = _make_intron(upstream, seq, "CCCCCCCCCC")
        shift = search_u12_boundary(intron, upstream.upper(), seq.upper())
        assert shift is None

    def test_motif_at_expected_position_returns_none(self):
        """Returns None when motif is already at the correct position (shift=0)."""
        upstream = "AAAAACCCAG"
        seq = "ATATCCTTTCCCCCCCCCCCCCCCCCCCAG"
        intron = _make_intron(upstream, seq, "CCCCCCCCCC")
        shift = search_u12_boundary(intron, upstream.upper(), seq.upper())
        assert shift is None  # shift=0 → no correction needed


class TestRequireCanonical:
    """Test the require_canonical option in correct_intron_if_needed."""

    def test_canonical_result_corrected_both_modes(self):
        """Correction producing canonical GT-AG is applied in both modes."""
        # shift=+2 → new 5' = GT, new 3' = AG (canonical)
        upstream = "CCCCCCCCAA"
        seq = "AAGTATCCTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTT"
        downstream = "AGCCCCCCCC"

        intron_strict = _make_intron(upstream, seq, downstream)
        _, corrected_strict = correct_intron_if_needed(
            intron_strict, require_canonical=True
        )
        assert corrected_strict

        intron_lax = _make_intron(upstream, seq, downstream)
        _, corrected_lax = correct_intron_if_needed(
            intron_lax, require_canonical=False
        )
        assert corrected_lax

    def test_noncanonical_result_rejected_when_required(self):
        """Correction producing non-canonical 3' is rejected when require_canonical=True."""
        # GAPDHP76-like: TA-GC → shift +1 → AT-CC (non-canonical 3')
        upstream = "AGAACAAGGT"
        seq = ("TATATCCTCCGCCTTTGTTTCCCTTGCCTTCTGGCCTGGCCTGTTGCC"
               "CTGCTGCGCAACAATCCCAACTGCTCCTGAGCCACCGTGCCTGGC")
        downstream = "CGTCTCCTCT"

        intron = _make_intron(upstream, seq, downstream)
        _, corrected = correct_intron_if_needed(intron, require_canonical=True)
        assert not corrected

    def test_noncanonical_result_accepted_when_not_required(self):
        """Correction producing non-canonical 3' is accepted when require_canonical=False."""
        upstream = "AGAACAAGGT"
        seq = ("TATATCCTCCGCCTTTGTTTCCCTTGCCTTCTGGCCTGGCCTGTTGCC"
               "CTGCTGCGCAACAATCCCAACTGCTCCTGAGCCACCGTGCCTGGC")
        downstream = "CGTCTCCTCT"

        intron = _make_intron(upstream, seq, downstream)
        _, corrected = correct_intron_if_needed(intron, require_canonical=False)
        assert corrected

    def test_canonical_intron_never_corrected(self):
        """Canonical GT-AG introns are never checked for correction."""
        seq = "GTATCCTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAG"
        intron = _make_intron("AAAAAAAAAA", seq, "TTTTTTTTTT", noncanonical=False)
        _, corrected = correct_intron_if_needed(intron, require_canonical=True)
        assert not corrected

    def test_correction_disabled(self):
        """No correction when correction_enabled=False."""
        upstream = "CCCCCCCCAA"
        seq = "AAGTATCCTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTT"
        downstream = "AGCCCCCCCC"
        intron = _make_intron(upstream, seq, downstream)
        _, corrected = correct_intron_if_needed(
            intron, correction_enabled=False, require_canonical=False
        )
        assert not corrected


class TestGetCorrectedDinucleotides:
    """Test dinucleotide prediction after shift."""

    def test_shift_produces_gt_ag(self):
        upstream = "CCCCCCCCAA"
        seq = "AAGTATCCTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTT"
        downstream = "AGCCCCCCCC"
        result = get_corrected_dinucleotides(upstream, seq, downstream, shift=2)
        assert result == ("GT", "AG")

    def test_shift_produces_at_cc(self):
        upstream = "AGAACAAGGT"
        seq = ("TATATCCTCCGCCTTTGTTTCCCTTGCCTTCTGGCCTGGCCTGTTGCC"
               "CTGCTGCGCAACAATCCCAACTGCTCCTGAGCCACCGTGCCTGGC")
        downstream = "CGTCTCCTCT"
        result = get_corrected_dinucleotides(upstream, seq, downstream, shift=1)
        assert result == ("AT", "CC")
