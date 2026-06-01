"""
Tests for species-specific U2 background correction.
"""

import numpy as np
import pytest

from intronIC.scoring.background import (
    BackgroundConfig,
    SpeciesBackground,
    _RegionAccumulator,
    _BPSAccumulator,
)
from intronIC.scoring.pwm import PWM, PWMSet


# ============================================================================
# Fixtures
# ============================================================================

def _make_uniform_pwm(name, length, start_index, ref_offset=0):
    """Create a uniform (0.25) PWM for testing."""
    return PWM(
        name=name,
        matrix=np.full((4, length), 0.25),
        length=length,
        start_index=start_index,
        reference_offset=ref_offset,
    )


@pytest.fixture
def human_u2_pwm_sets():
    """Minimal human U2 PWM sets for testing."""
    five_pwm = _make_uniform_pwm('u2_gtag_five', 40, -20, ref_offset=20)
    three_pwm = _make_uniform_pwm('u2_gtag_three', 40, -20, ref_offset=20)
    bp_pwm = _make_uniform_pwm('u2_gtag_bp', 12, 0, ref_offset=9)
    return {
        'five': PWMSet(matrices={('u2', 'gtag'): five_pwm}),
        'three': PWMSet(matrices={('u2', 'gtag'): three_pwm}),
        'bp': PWMSet(matrices={('u2', 'gtag'): bp_pwm}),
    }


@pytest.fixture
def default_config():
    return BackgroundConfig(
        enabled=True, n0=1000, trim_percentile=5,
        pseudocount_per_base=1.0, n_iterations=1, min_introns=10,
    )


# ============================================================================
# RegionAccumulator tests
# ============================================================================

class TestRegionAccumulator:

    def test_accumulate_counts(self):
        acc = _RegionAccumulator(n_positions=4)
        acc.add('GT', 'ACGT')
        acc.add('GT', 'ACGT')
        acc.add('GT', 'AAAA')

        freqs, n = acc.get_frequencies('GT', pseudocount=0)
        assert n == 3
        # Position 0: A=3, C=0, G=0, T=0
        assert freqs[0, 0] == 1.0  # all A at position 0
        # Position 1: C=2, A=1
        assert freqs[1, 1] == 2/3  # C at position 1

    def test_pseudocount_prevents_zeros(self):
        acc = _RegionAccumulator(n_positions=4)
        acc.add('GT', 'AAAA')

        freqs, n = acc.get_frequencies('GT', pseudocount=1.0)
        # Position 0: A=1+1=2, C=0+1=1, G=0+1=1, T=0+1=1 → total=5
        assert freqs[0, 0] == pytest.approx(2/5)
        assert freqs[0, 1] == pytest.approx(1/5)

    def test_missing_subtype_returns_none(self):
        acc = _RegionAccumulator(n_positions=4)
        assert acc.get_frequencies('GC') is None

    def test_short_sequence_ignored(self):
        acc = _RegionAccumulator(n_positions=4)
        acc.add('GT', 'AC')  # too short
        assert acc.get_frequencies('GT') is None

    def test_subtract_introns(self):
        acc = _RegionAccumulator(n_positions=4)
        acc.add('GT', 'ACGT')
        acc.add('GT', 'ACGT')
        acc.add('GT', 'TTTT')

        # Remove one ACGT
        acc.subtract_introns('GT', ['ACGT'])

        freqs, n = acc.get_frequencies('GT', pseudocount=0)
        assert n == 2
        # Position 0: A=1 (from remaining ACGT), T=1 (from TTTT)
        assert freqs[0, 0] == pytest.approx(0.5)
        assert freqs[0, 3] == pytest.approx(0.5)


# ============================================================================
# BPSAccumulator tests
# ============================================================================

class TestBPSAccumulator:

    def test_uniform_composition(self):
        acc = _BPSAccumulator()
        acc.add('GT', 'ACGTACGTACGT')  # equal base composition
        freqs, n = acc.get_frequencies('GT', n_positions=12, pseudocount=0)
        assert n == 1
        # Each base appears 3 times in 12 chars → 0.25 each
        assert freqs[0, 0] == pytest.approx(0.25)
        assert freqs.shape == (12, 4)

    def test_tiled_to_n_positions(self):
        acc = _BPSAccumulator()
        acc.add('GT', 'AAAAAAAAAAAA')
        freqs, n = acc.get_frequencies('GT', n_positions=12, pseudocount=0)
        # All A → freqs should be [1, 0, 0, 0] at every position
        assert np.all(freqs[:, 0] == 1.0)
        assert np.all(freqs[:, 1] == 0.0)


# ============================================================================
# SpeciesBackground integration tests
# ============================================================================

class TestSpeciesBackground:

    def test_min_introns_fallback(self, human_u2_pwm_sets, default_config):
        """Below min_introns, returns human U2 PWMs unchanged."""
        config = BackgroundConfig(
            enabled=True, n0=1000, trim_percentile=5,
            pseudocount_per_base=1.0, n_iterations=1,
            min_introns=100,  # require 100
        )
        bg = SpeciesBackground(human_u2_pwm_sets, config)

        # Add only 5 introns
        for i in range(5):
            bg.accumulate(f'intron_{i}', 'GT', 'AG',
                          'ACGTACGTACGT', 'ACGTACGTAC', 'A' * 50)

        result = bg.build_corrected_pwm_sets()
        # Should be unchanged human PWMs
        assert ('u2', 'gtag') in result['five'].matrices

    def test_blending_weight(self, human_u2_pwm_sets, default_config):
        """Blending weight w = n/(n+n0) produces correct mix."""
        config = BackgroundConfig(
            enabled=True, n0=100, trim_percentile=0,  # no trimming
            pseudocount_per_base=0.0001, n_iterations=0,
            min_introns=10,
        )
        bg = SpeciesBackground(human_u2_pwm_sets, config, five_len=4)

        # Add 100 introns all with 'AAAA' at 5'SS → empirical is A=1.0
        for i in range(100):
            bg.accumulate(f'intron_{i}', 'GT', 'AG',
                          'AAAA', 'ACGTACGTAC', 'A' * 50)

        result = bg.build_corrected_pwm_sets()
        u2_five = result['five'].matrices[('u2', 'gtag')]

        # w = 100/(100+100) = 0.5
        # At scored positions (-3 to +0, indices 17-20 in the 40-pos PWM):
        # blended_A = 0.5 * ~1.0 + 0.5 * 0.25 ≈ 0.625
        # Check position -3 (array index 17)
        a_freq = u2_five.matrix[0, 17]  # A at position -3
        assert a_freq > 0.5  # should be ~0.625

    def test_accumulate_and_count(self, human_u2_pwm_sets, default_config):
        bg = SpeciesBackground(human_u2_pwm_sets, default_config)
        assert bg.n_introns == 0
        bg.accumulate('i1', 'GT', 'AG', 'ACGTACGTACGT', 'ACGTACGTAC', 'A'*50)
        assert bg.n_introns == 1

    def test_trim_excludes_high_scorers(self, human_u2_pwm_sets):
        """Trimming removes introns with high initial scores."""
        config = BackgroundConfig(
            enabled=True, n0=10, trim_percentile=50,  # trim top 50%
            pseudocount_per_base=1.0, n_iterations=1,
            min_introns=2,
        )
        bg = SpeciesBackground(human_u2_pwm_sets, config, five_len=4)

        # Add 4 introns, 2 with high scores and 2 with low
        bg.accumulate('high1', 'GT', 'AG', 'GGGG', 'ACGTACGTAC', 'A'*50)
        bg.accumulate('high2', 'GT', 'AG', 'GGGG', 'ACGTACGTAC', 'A'*50)
        bg.accumulate('low1', 'GT', 'AG', 'AAAA', 'ACGTACGTAC', 'A'*50)
        bg.accumulate('low2', 'GT', 'AG', 'AAAA', 'ACGTACGTAC', 'A'*50)

        bg.set_initial_scores('high1', five_raw=20.0, bp_raw=0.0, three_raw=0.0)
        bg.set_initial_scores('high2', five_raw=18.0, bp_raw=0.0, three_raw=0.0)
        bg.set_initial_scores('low1', five_raw=-5.0, bp_raw=0.0, three_raw=0.0)
        bg.set_initial_scores('low2', five_raw=-3.0, bp_raw=0.0, three_raw=0.0)

        result = bg.build_corrected_pwm_sets()
        # After trimming top 50% (high1, high2), only low1+low2 remain
        # Both have 'AAAA', so empirical should be A-enriched
        u2_five = result['five'].matrices[('u2', 'gtag')]
        # At position -3 (idx 17): A should be dominant after trimming GGGG
        a_freq = u2_five.matrix[0, 17]
        g_freq = u2_five.matrix[2, 17]
        assert a_freq > g_freq

    def test_low_n_noncanonical_dnt_does_not_clobber(
        self, human_u2_pwm_sets, default_config
    ):
        """Regression (U2 subtype-clobber bug): a rare non-canonical 5' dnt that
        defaults to 'gtag' (e.g. TT) must NOT overwrite the high-n GT-derived gtag
        background. The previous code sorted dnts and did last-writer-wins, so a
        handful of TT introns (n=3, blend weight ~0 => ~human prior) clobbered the
        GT background (n=2000), silently disabling the species correction on any
        annotation carrying non-canonical dnts. Fix: highest-n dnt defines the slot."""
        bg = SpeciesBackground(
            human_u2_pwm_sets, default_config,
            five_len=12, three_len=10, bp_len=12,
        )
        # 2000 canonical GT-AG introns, 5' window all-A => strongly species-shifted
        # away from the uniform (0.25) human prior.
        for i in range(2000):
            bg.accumulate(f'gt_{i}', 'GT', 'AG', 'A' * 12, 'C' * 10, 'A' * 50)
        # 3 non-canonical TT introns (5' dnt 'TT' -> defaults to 'gtag'), 5' all-G.
        # 'TT' > 'GT' alphabetically, so last-writer-wins would let these win.
        for i in range(3):
            bg.accumulate(f'tt_{i}', 'TT', 'AG', 'G' * 12, 'C' * 10, 'A' * 50)

        five = bg.build_corrected_pwm_sets()['five'].matrices[('u2', 'gtag')].matrix
        # Scored window is cols 17-28, A is row 0. Under the fix the gtag background
        # reflects the 2000 GT introns (A ~0.75); under the clobber bug it reverts
        # to ~human uniform (~0.25).
        a_scored = float(five[0, 17:29].mean())
        assert a_scored > 0.5, (
            f"gtag 5' A-freq={a_scored:.3f} in scored window: expected species-shifted "
            f"(~0.75 from 2000 GT introns); ~0.25 means 3 TT introns clobbered it to human"
        )
