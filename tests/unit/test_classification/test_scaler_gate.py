"""Unit tests for the Path C' scaler gate logic.

Covers compute_asymmetry, decide_route, apply_route. Synthetic populations
exercise the three regimes (adaptive / frozen / strict) plus the asymmetric
n_total floor edge case (ChlRei-class species with low n_total + strict
pattern).
"""
import numpy as np
import pytest

from intronIC.classification.scaler_gate import (
    apply_route,
    compute_asymmetry,
    decide_route,
)


def _make_scores(n_both, n_only_a, n_only_f, n_neither, thr=90.0):
    """Construct adaptive/frozen score arrays with prescribed asymmetry."""
    a = np.concatenate([
        np.full(n_both, 95.0),     # both modes call
        np.full(n_only_a, 95.0),   # adaptive only
        np.full(n_only_f, 50.0),   # frozen only (a below thr)
        np.full(n_neither, 10.0),  # neither
    ])
    f = np.concatenate([
        np.full(n_both, 95.0),
        np.full(n_only_a, 50.0),
        np.full(n_only_f, 95.0),
        np.full(n_neither, 10.0),
    ])
    return a, f


class TestComputeAsymmetry:
    def test_simple_case(self):
        a, f = _make_scores(n_both=10, n_only_a=2, n_only_f=3, n_neither=100)
        s = compute_asymmetry(a, f, threshold=90.0)
        assert s["n_total"] == 15
        assert s["n_both"] == 10
        assert s["n_only_a"] == 2
        assert s["n_only_f"] == 3
        assert s["frac_only_a"] == pytest.approx(2 / 15)
        assert s["frac_only_f"] == pytest.approx(3 / 15)
        assert s["frac_both"] == pytest.approx(10 / 15)

    def test_empty(self):
        a = np.zeros(100)
        f = np.zeros(100)
        s = compute_asymmetry(a, f)
        assert s["n_total"] == 0
        # Should not divide by zero
        assert s["frac_only_a"] == 0.0
        assert s["frac_only_f"] == 0.0

    def test_shape_mismatch_raises(self):
        with pytest.raises(ValueError):
            compute_asymmetry(np.zeros(10), np.zeros(11))


class TestDecideRoute:
    """Test the three regimes + edge cases."""

    def test_adaptive_well_agreed_vertebrate(self):
        # HomSap-like: high both, modest only_a
        # HomSap actual: only_a%=15.4%, only_f%=0%, both=653, n_total=772
        s = compute_asymmetry(*_make_scores(653, 119, 0, 0))
        assert decide_route(s) == "adaptive"

    def test_frozen_apostasia_class(self):
        # ApoShe-like: only_f% >= 10%, n_total >= 5
        # Actual: only_f%=19.9%, n_total=141
        s = compute_asymmetry(*_make_scores(113, 0, 28, 0))
        assert decide_route(s) == "frozen"

    def test_strict_tetrahymena_class(self):
        # TetThe-like: only_a 100%, n_both 0
        s = compute_asymmetry(*_make_scores(0, 49, 0, 0))
        assert decide_route(s) == "strict"

    def test_strict_low_ntotal_chlrei(self):
        # ChlRei-class: only_a 100%, n_both=0, n_total=2.
        # Revised rule: strict has no n_total floor.
        s = compute_asymmetry(*_make_scores(0, 2, 0, 0))
        assert decide_route(s) == "strict", (
            "Low n_total with 100% only_a + zero both should route strict"
        )

    def test_adaptive_quiet_default(self):
        # Saccer-class: no high-confidence calls in any mode.
        s = compute_asymmetry(*_make_scores(0, 0, 0, 100))
        assert decide_route(s) == "adaptive"

    def test_frozen_not_fired_on_low_ntotal(self):
        # Edge case: only_f=1 of n_total=2 (50%). Frozen has n_total floor,
        # so this should NOT route frozen (avoid over-firing on noise).
        s = compute_asymmetry(*_make_scores(1, 0, 1, 0))
        assert decide_route(s) == "adaptive"

    def test_intermediate_only_a_routes_adaptive(self):
        # 20% only_a, big both base (HomSap-like): adaptive
        s = compute_asymmetry(*_make_scores(400, 100, 0, 0))
        assert s["frac_only_a"] == pytest.approx(0.20)
        assert decide_route(s) == "adaptive"

    def test_high_only_a_with_some_both_routes_adaptive(self):
        # 50% only_a but n_both > n_both_max: should route adaptive,
        # not strict. (Strict requires both: only_a high AND n_both small.)
        # Build: n_both=10, n_only_a=10 → only_a=50%, but n_both=10 > 5
        s = compute_asymmetry(*_make_scores(10, 10, 0, 0))
        assert s["frac_only_a"] == pytest.approx(0.50)
        assert s["n_both"] == 10
        assert decide_route(s) == "adaptive"


class TestApplyRoute:
    def test_adaptive_picks_a(self):
        a = np.array([95.0, 50.0, 30.0])
        f = np.array([20.0, 95.0, 30.0])
        final, called = apply_route(a, f, "adaptive")
        assert np.array_equal(final, a)
        assert np.array_equal(called, np.array([True, False, False]))

    def test_frozen_picks_f(self):
        a = np.array([95.0, 50.0, 30.0])
        f = np.array([20.0, 95.0, 30.0])
        final, called = apply_route(a, f, "frozen")
        assert np.array_equal(final, f)
        assert np.array_equal(called, np.array([False, True, False]))

    def test_strict_min_and_both_required(self):
        a = np.array([95.0, 95.0, 50.0])
        f = np.array([95.0, 50.0, 95.0])
        final, called = apply_route(a, f, "strict")
        # final = min(a, f)
        assert np.array_equal(final, np.array([95.0, 50.0, 50.0]))
        # only the first row has both >= 90
        assert np.array_equal(called, np.array([True, False, False]))

    def test_invalid_route_raises(self):
        a = np.array([95.0])
        f = np.array([95.0])
        with pytest.raises(ValueError):
            apply_route(a, f, "nonsense")  # type: ignore[arg-type]


class TestThresholdRobustness:
    """Confirm the gate is stable in the empirically-tested ranges."""

    def test_only_f_thr_range_stable_for_apostasia(self):
        # ApoShe-like: only_f% = 19.9%. Should stay frozen for any only_f_thr
        # in [0.01, 0.19].
        s = compute_asymmetry(*_make_scores(113, 0, 28, 0))
        for tau in [0.01, 0.05, 0.10, 0.15, 0.19]:
            assert decide_route(s, only_f_thr=tau) == "frozen", f"failed at tau={tau}"
        # And drops out at >= 0.20
        assert decide_route(s, only_f_thr=0.20) == "adaptive"

    def test_only_a_thr_range_stable_for_tetrahymena(self):
        # TetThe: only_a=100%. Strict for any only_a_thr in [0.01, 1.0].
        s = compute_asymmetry(*_make_scores(0, 49, 0, 0))
        for tau in [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]:
            assert decide_route(s, only_a_thr=tau) == "strict", f"failed at tau={tau}"

    def test_n_both_max_doesnt_affect_clean_strict_case(self):
        # TetThe has n_both=0, so any n_both_max >= 0 should route strict.
        s = compute_asymmetry(*_make_scores(0, 49, 0, 0))
        for n_max in [0, 1, 2, 5, 10, 50]:
            assert decide_route(s, n_both_max=n_max) == "strict"
