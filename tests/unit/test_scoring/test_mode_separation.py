"""Unit tests for `intronIC.scoring.mode_separation`.

Covers:
- candidate-weight sigmoid shape
- weighted-median quantile helper (zero-weight fallback, edge cases)
- mode-separation fit + transform: U2 → z=0, U12 → z=1 by construction
- compute_support2 consensus feature
- evaluate_gate: passes / n_floor / location-prior / degenerate-separation /
  no-valley fail modes

These are isolation tests; the heavier classify-pipeline integration
lives in tests/integration.
"""
import numpy as np
import pytest

from intronIC.scoring.mode_separation import (
    DEFAULT_MU_U12_TOLERANCE_5P,
    DEFAULT_UNIVERSAL_ANCHORS,
    ModeSeparationStats,
    apply_mode_separation_z,
    candidate_weight_from_svm,
    compute_support2,
    evaluate_gate,
    fit_mode_separation,
)


# ---------------------------------------------------------------------------
# candidate_weight_from_svm
# ---------------------------------------------------------------------------

def test_candidate_weight_sigmoid_shape():
    svm = np.array([60.0, 80.0, 90.0, 100.0, 120.0])
    w = candidate_weight_from_svm(svm)
    assert w[2] == pytest.approx(0.5, abs=1e-9)          # at center
    assert w[1] < 0.2                                     # below center
    assert w[3] > 0.8                                     # above center
    assert np.all(np.diff(w) > 0)                         # monotonic


def test_candidate_weight_default_center_steepness_match_module_defaults():
    # Defaults: center=90, steepness=5 — at svm=85 weight should be ~0.27
    w = candidate_weight_from_svm(np.array([85.0]))[0]
    assert w == pytest.approx(0.2689414, abs=1e-5)


# ---------------------------------------------------------------------------
# fit_mode_separation + apply_mode_separation_z
# ---------------------------------------------------------------------------

def test_mode_sep_collapses_to_unit_anchor():
    rng = np.random.default_rng(0)
    n_u2 = 5000
    n_u12 = 1000  # large enough to make weighted-median estimate stable
    u2 = rng.normal(-40.0, 8.0, n_u2)
    u12 = rng.normal(15.0, 4.0, n_u12)
    raw = np.concatenate([u2, u12])
    # Perfect candidate weights: U2 = 0, U12 = 1
    cw = np.concatenate([np.zeros(n_u2), np.ones(n_u12)])

    stats = fit_mode_separation(raw, cw)
    assert stats.mu_u2 == pytest.approx(-40.0, abs=1.0)
    assert stats.mu_u12 == pytest.approx(15.0, abs=1.0)
    assert stats.separation > 50.0

    z = apply_mode_separation_z(raw, stats)
    # After mode-sep the U2 region centers on 0 and U12 region on 1, modulo
    # within-population scatter (~1.5/sep). Use loose tolerance — the test
    # asserts the calibration shape, not numerical equality.
    assert z[:n_u2].mean() == pytest.approx(0.0, abs=0.1)
    assert z[n_u2:].mean() == pytest.approx(1.0, abs=0.1)


def test_mode_sep_degenerate_raises():
    stats = ModeSeparationStats(
        mu_u2=5.0, mu_u12=2.0, spread_u2=1.0, separation=-3.0,
        n_eff_candidates=10.0,
    )
    with pytest.raises(ValueError, match="degenerate"):
        apply_mode_separation_z(np.array([1.0]), stats)


def test_mode_sep_shape_mismatch_raises():
    raw = np.array([1.0, 2.0, 3.0])
    cw = np.array([0.5, 0.5])
    with pytest.raises(ValueError, match="same shape"):
        fit_mode_separation(raw, cw)


def test_mode_sep_n_eff_equals_weight_sum():
    cw = np.array([0.1, 0.3, 0.8, 0.2])
    raw = np.array([0.0, 1.0, 2.0, 3.0])
    stats = fit_mode_separation(raw, cw)
    assert stats.n_eff_candidates == pytest.approx(0.1 + 0.3 + 0.8 + 0.2)


def test_weighted_quantile_zero_weights_falls_back_to_plain():
    raw = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    cw = np.zeros(5)
    stats = fit_mode_separation(raw, cw)
    # Both modes degenerate to the plain median when all weights are zero.
    assert stats.mu_u2 == pytest.approx(3.0)
    assert stats.mu_u12 == pytest.approx(3.0)


# ---------------------------------------------------------------------------
# compute_support2
# ---------------------------------------------------------------------------

def test_support2_picks_second_largest_clipped_positive():
    # Per-intron z-triple, then clip negatives, then take the second-largest.
    z5 = np.array([0.8, -0.2,  0.0])
    zbp = np.array([0.5,  0.6,  0.4])
    z3 = np.array([0.3, -0.1,  0.7])
    # Row 0: clipped (0.8, 0.5, 0.3) → sorted (0.3, 0.5, 0.8) → 2nd-largest 0.5
    # Row 1: clipped (0.0, 0.6, 0.0) → sorted (0.0, 0.0, 0.6) → 2nd-largest 0.0
    # Row 2: clipped (0.0, 0.4, 0.7) → sorted (0.0, 0.4, 0.7) → 2nd-largest 0.4
    out = compute_support2(z5, zbp, z3)
    assert np.allclose(out, [0.5, 0.0, 0.4])


# ---------------------------------------------------------------------------
# evaluate_gate
# ---------------------------------------------------------------------------

def _good_points():
    rng = np.random.default_rng(0)
    u2 = rng.normal([-40.0, -5.0, 0.0], [4.0, 2.0, 2.0], (500, 3))
    u12 = rng.normal([15.0,  5.0, 1.0], [3.0, 2.0, 2.0], (60, 3))
    return u12, u2


def _strong_valley(u12_pts, u2_pts):
    return {"median_depth": 0.65}


def _weak_valley(u12_pts, u2_pts):
    return {"median_depth": 0.10}


def test_evaluate_gate_passes_clean_species():
    u12, u2 = _good_points()
    g = evaluate_gate(
        u12, u2,
        mu_u12_5p_raw=DEFAULT_UNIVERSAL_ANCHORS["five_raw"],
        n_eff_candidates=60.0,
        separation_5p=55.0,
        valley_depth_fn=_strong_valley,
    )
    assert g.passes
    assert g.reason == "ok"
    assert g.valley_depth == pytest.approx(0.65)


def test_evaluate_gate_below_n_floor_fails():
    u12, u2 = _good_points()
    g = evaluate_gate(
        u12, u2,
        mu_u12_5p_raw=DEFAULT_UNIVERSAL_ANCHORS["five_raw"],
        n_eff_candidates=2.0,  # below default floor 5
        separation_5p=55.0,
        valley_depth_fn=_strong_valley,
    )
    assert not g.passes
    assert g.reason == "below_n_floor"


def test_evaluate_gate_degenerate_separation_fails():
    u12, u2 = _good_points()
    g = evaluate_gate(
        u12, u2,
        mu_u12_5p_raw=DEFAULT_UNIVERSAL_ANCHORS["five_raw"],
        n_eff_candidates=60.0,
        separation_5p=-0.5,
        valley_depth_fn=_strong_valley,
    )
    assert not g.passes
    assert g.reason == "degenerate_separation"


def test_evaluate_gate_location_prior_catches_mode_drift():
    """v3-era-style failure: μ_U12 lands far from the cross-species anchor."""
    u12, u2 = _good_points()
    drifted_mu_u12 = DEFAULT_UNIVERSAL_ANCHORS["five_raw"] + DEFAULT_MU_U12_TOLERANCE_5P + 2.0
    g = evaluate_gate(
        u12, u2,
        mu_u12_5p_raw=drifted_mu_u12,
        n_eff_candidates=60.0,
        separation_5p=55.0,
        valley_depth_fn=_strong_valley,  # would pass on valley alone
    )
    assert not g.passes
    assert g.reason == "u12_mode_outside_prior_range"
    assert abs(g.mu_u12_offset) > DEFAULT_MU_U12_TOLERANCE_5P


def test_evaluate_gate_no_valley_fails():
    u12, u2 = _good_points()
    g = evaluate_gate(
        u12, u2,
        mu_u12_5p_raw=DEFAULT_UNIVERSAL_ANCHORS["five_raw"],
        n_eff_candidates=60.0,
        separation_5p=55.0,
        valley_depth_fn=_weak_valley,
    )
    assert not g.passes
    assert g.reason == "no_kde_valley"
    assert g.valley_depth == pytest.approx(0.10)
