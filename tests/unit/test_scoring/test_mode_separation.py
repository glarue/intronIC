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
    DEFAULT_DISCOUNT_K_OVERCALL,
    DEFAULT_DISCOUNT_K_WEAKMOT,
    DEFAULT_DISCOUNT_TAU_MOTIF,
    DEFAULT_DISCOUNT_TAU_OVERCALL,
    ModeSeparationStats,
    apply_continuous_discount,
    apply_mode_separation_z,
    candidate_weight_from_svm,
    compute_boundary_mass,
    compute_support2,
    continuous_discount_logit_shift,
    evaluate_gate,
    fit_mode_separation,
    voting_fraction,
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


# Mocks return the full compute_valley_depth contract used by the v3 gate (Steps 1-4):
# the decision keys are gap_fraction (#4) and centroid_sigma (#3); median_depth is diagnostic.
def _strong_valley(u12_pts, u2_pts):
    return {"median_depth": 0.65, "gap_fraction": 0.45,
            "gap_fraction_ucl": 0.47, "centroid_sigma": 4.8}


def _weak_valley(u12_pts, u2_pts):
    # No positive gap → fails check #4 (below_gap_fraction).
    return {"median_depth": 0.10, "gap_fraction": -0.20,
            "gap_fraction_ucl": -0.15, "centroid_sigma": 2.0}


def _low_csig_valley(u12_pts, u2_pts):
    # Positive gap but weak centroid separation → fails check #3 (low_centroid_sigma).
    return {"median_depth": 0.50, "gap_fraction": 0.30,
            "gap_fraction_ucl": 0.32, "centroid_sigma": 3.0}


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
    assert g.valley_depth == pytest.approx(0.65)        # diagnostic
    assert g.gap_fraction == pytest.approx(0.45)        # check #4 (decision)
    assert g.centroid_sigma == pytest.approx(4.8)       # check #3 replacement (decision)


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


def test_evaluate_gate_mode_drift_no_longer_overkills():
    """v3 Step 4: the absolute μ_U12 anchor check was REMOVED (it over-killed 24 snRNA-confirmed
    real bearers whose μ_U12 is shifted low by estimator bias / motif divergence). A drifted μ_U12
    with genuine species-internal separation now PASSES; the offset is retained as a diagnostic."""
    u12, u2 = _good_points()
    drifted_mu_u12 = DEFAULT_UNIVERSAL_ANCHORS["five_raw"] - 9.0  # far below the old band, like Anopheles
    g = evaluate_gate(
        u12, u2,
        mu_u12_5p_raw=drifted_mu_u12,
        n_eff_candidates=60.0,
        separation_5p=55.0,
        valley_depth_fn=_strong_valley,  # strong species-internal separation (gap + csig)
    )
    assert g.passes
    assert g.reason == "ok"
    assert abs(g.mu_u12_offset) > DEFAULT_MU_U12_TOLERANCE_5P  # offset still computed (diagnostic)


def test_evaluate_gate_low_centroid_sigma_fails():
    """v3 Step 4: the centroid_sigma floor (check #3 replacement) catches a weak candidate cluster
    that has a positive gap but insufficient U2-σ separation (z-anchor-corruption guard)."""
    u12, u2 = _good_points()
    g = evaluate_gate(
        u12, u2,
        mu_u12_5p_raw=DEFAULT_UNIVERSAL_ANCHORS["five_raw"],
        n_eff_candidates=60.0,
        separation_5p=55.0,
        valley_depth_fn=_low_csig_valley,
    )
    assert not g.passes
    assert g.reason == "low_centroid_sigma"
    assert g.centroid_sigma == pytest.approx(3.0)


def test_evaluate_gate_below_gap_fraction_fails():
    """v3 Step 2: no positive gap on the Fisher axis → below_gap_fraction (was 'no_kde_valley')."""
    u12, u2 = _good_points()
    g = evaluate_gate(
        u12, u2,
        mu_u12_5p_raw=DEFAULT_UNIVERSAL_ANCHORS["five_raw"],
        n_eff_candidates=60.0,
        separation_5p=55.0,
        valley_depth_fn=_weak_valley,
    )
    assert not g.passes
    assert g.reason == "below_gap_fraction"
    assert g.valley_depth == pytest.approx(0.10)  # median_depth still surfaced as diagnostic


# ---------------------------------------------------------------------------
# v2.7: continuous per-intron discount
# ---------------------------------------------------------------------------

def test_continuous_discount_zero_in_healthy_regime():
    """Strong-motif U12-call has no penalty (raw_sum > tau_motif, svm_vs_naive < 0)."""
    # logit(0.99) ≈ +4.6; raw_sum = +22 → svm_vs_naive ≈ -17.4 (well below 0)
    # raw_sum 22 > tau_motif 10 → no weak-motif penalty
    shift = continuous_discount_logit_shift(
        logit_svm=np.array([4.6]),
        raw_sum=np.array([22.0]),
    )
    assert shift[0] == pytest.approx(0.0)


def test_continuous_discount_overcall_penalty():
    """High svm + negative raw_sum → strong overcall penalty (default k_oc=2.0)."""
    # logit(0.99) ≈ 4.6, raw_sum = -3 → svm_vs_naive = 7.6
    # k_overcall=2.0 → overcall penalty = 15.2
    # k_weakmot=0.0 (default) → no weak-motif penalty
    # Total shift = 15.2
    shift = continuous_discount_logit_shift(
        logit_svm=np.array([4.6]),
        raw_sum=np.array([-3.0]),
    )
    assert shift[0] == pytest.approx(15.2, abs=0.01)


def test_continuous_discount_weakmot_disabled_by_default():
    """v2.7 default: k_weakmot=0 → no weak-motif penalty.

    Borderline TP regime (moderate svm, weak motif, no overcall) gets
    NO penalty under defaults, preserving IPA-validated TPs like SUDS3,
    ARPC5, ap4e1, mios that have raw_sum < 10.
    """
    # logit(0.7) ≈ 0.847, raw_sum = 5 → svm_vs_naive = -4.15
    # k_overcall=2.0, tau=0 → max(0, -4.15) = 0 → no overcall penalty
    # k_weakmot=0.0 (default) → no weak-motif penalty
    # Total shift = 0
    shift = continuous_discount_logit_shift(
        logit_svm=np.array([0.847]),
        raw_sum=np.array([5.0]),
    )
    assert shift[0] == pytest.approx(0.0, abs=0.01)


def test_continuous_discount_weakmot_when_enabled():
    """Explicit k_weakmot > 0 fires the weak-motif penalty (opt-in)."""
    shift = continuous_discount_logit_shift(
        logit_svm=np.array([0.847]),
        raw_sum=np.array([5.0]),
        k_weakmot=0.20, tau_motif=10.0,
    )
    # 0.20 * (10 - 5) = 1.0
    assert shift[0] == pytest.approx(1.0, abs=0.01)


def test_apply_continuous_discount_monotone_in_overcall():
    """Higher svm_vs_naive → larger overcall penalty → lower adjusted score.

    Defaults: k_overcall=2.0, k_weakmot=0.0. So only overcall fires.
    """
    svm = np.array([99.0, 99.0, 99.0])
    # All three at svm=99 → logit≈4.595. svm_vs_naive = 4.595 - raw_sum.
    # raw_sums chosen so:
    #   +22: svm_vs_naive=-17.4 → no penalty (well below tau_oc=0)
    #   -1:  svm_vs_naive=+5.6 → moderate overcall penalty = 2*5.6 = 11.2
    #   -8:  svm_vs_naive=+12.6 → strong overcall penalty = 2*12.6 = 25.2
    raw_sums = np.array([+22.0, -1.0, -8.0])
    adj = apply_continuous_discount(svm, raw_sums)
    assert adj[0] >= 98.0   # no penalty (svm_vs_naive < 0)
    assert adj[1] < adj[0]  # moderate overcall penalty
    assert adj[2] < adj[1]  # stronger overcall penalty


def test_apply_continuous_discount_no_increase():
    """Discount is always non-positive — adjusted ≤ original."""
    rng = np.random.default_rng(0)
    svm = rng.uniform(0, 100, 100)
    raw_sums = rng.uniform(-30, 30, 100)
    adj = apply_continuous_discount(svm, raw_sums)
    # Allow tiny numerical error
    assert np.all(adj <= svm + 1e-6)


def test_apply_continuous_discount_custom_params():
    """k=0 disables discount entirely."""
    svm = np.array([99.0, 95.0, 50.0])
    raw_sums = np.array([22.0, -3.0, 5.0])
    adj = apply_continuous_discount(
        svm, raw_sums,
        k_overcall=0.0, k_weakmot=0.0,
    )
    np.testing.assert_array_almost_equal(adj, svm, decimal=4)


# ---------------------------------------------------------------------------
# v2.7: voting fraction
# ---------------------------------------------------------------------------

def test_voting_fraction_unanimous():
    probas = np.full((10, 5), 0.9)  # all 10 models call all 5 introns U12
    v = voting_fraction(probas)
    np.testing.assert_array_equal(v, np.ones(5))


def test_voting_fraction_split():
    probas = np.array([[0.99] * 5, [0.6] * 5, [0.4] * 5, [0.01] * 5])
    # 2 of 4 models above 0.5 for each intron
    v = voting_fraction(probas)
    np.testing.assert_array_equal(v, np.full(5, 0.5))


def test_voting_fraction_custom_threshold():
    probas = np.array([[0.7] * 5, [0.3] * 5])
    v = voting_fraction(probas, threshold=0.6)
    np.testing.assert_array_equal(v, np.full(5, 0.5))  # only first model > 0.6


def test_voting_fraction_shape_error():
    with pytest.raises(ValueError, match="2-D"):
        voting_fraction(np.array([0.5, 0.6, 0.7]))  # 1-D input


# ---------------------------------------------------------------------------
# v2.7: boundary mass
# ---------------------------------------------------------------------------

def test_boundary_mass_clean_bimodal():
    """Bimodal distribution (all extremes) → BM ≈ 0."""
    mean_p = np.concatenate([np.full(900, 0.01), np.full(100, 0.99)])
    bm = compute_boundary_mass(mean_p)
    assert bm == pytest.approx(0.0)


def test_boundary_mass_diffuse():
    """All in middle → BM ≈ 1.0."""
    mean_p = np.full(1000, 0.5)
    bm = compute_boundary_mass(mean_p)
    assert bm == pytest.approx(1.0)


def test_boundary_mass_custom_bounds():
    mean_p = np.array([0.05, 0.15, 0.35, 0.55, 0.75, 0.95])
    bm = compute_boundary_mass(mean_p, lower=0.2, upper=0.8)
    # 0.35, 0.55, 0.75 in [0.2, 0.8] → 3/6 = 0.5
    assert bm == pytest.approx(0.5)


def test_boundary_mass_empty():
    bm = compute_boundary_mass(np.array([]))
    assert np.isnan(bm)
