"""Unit tests for the species-level U12 adjudicator (depth_tail q-driver + operational guards).

Covers the 2026-06-27 revised design (docs/raw_gated_scoring.md §0d):
  - PRIMARY q-driver is the size-invariant ``depth_tail`` (not the size-aware ``excess_z``);
  - ``excess_z``/``p_gumbel`` are retained as a labelled SECONDARY diagnostic;
  - every result carries an explicit ``AdjStatus`` (no silent NaN): ADJUDICATED / UNDETERMINED /
    LOW_N / DEGENERATE_TAIL / SCHEMA_FAIL;
  - low-N uncertainty surfaces as a WIDE CI -> UNDETERMINED, never as a confident "loss".
"""
import numpy as np
import pytest

from intronIC.scoring.species_adjudicator import (
    AdjudicatorParams,
    AdjStatus,
    adjudicate,
    q_from_depth_tail,
    ADJUDICATOR_PARAMS_VERSION,
)

P = AdjudicatorParams()


def _synth_genome(med_u2, call_margin, n_u2=5000, n_call=20, call_sd=0.2, seed=1):
    """A toy genome: a U2 bulk well below the call threshold + a tight cluster of motif-strong calls."""
    rng = np.random.RandomState(seed)
    u2 = rng.normal(med_u2, 1.0, n_u2)
    calls = rng.normal(call_margin, call_sd, n_call) if n_call else np.empty(0)
    margin = np.concatenate([u2, calls])
    p_motif = 1.0 / (1.0 + np.exp(-(P.platt_a * margin + P.platt_c)))
    return margin, p_motif


# --------------------------------------------------------------------------------------------------
# depth_tail q-driver
# --------------------------------------------------------------------------------------------------
def test_q_boundary_is_in_the_panel_gap():
    """q=0.5 must sit inside the empirical loss/bearer gap depth_tail in (2.91, 3.39)."""
    boundary = -P.q_b / P.q_a
    assert 2.91 < boundary < 3.39


@pytest.mark.parametrize("depth_tail,expect", [
    (1.30, "deep_loss"),    # Symbiodinium CONFLICT — must read clearly loss-side
    (1.01, "deep_loss"),    # paramecium loss
    (2.91, "loss"),         # chlamydomonas — the deepest loss, just below the boundary
    (3.39, "bearer"),       # neocallimastix — the shallowest bearer
    (4.53, "deep_bearer"),  # human
    (6.09, "deep_bearer"),  # drosophila
])
def test_q_from_depth_tail_panel_points(depth_tail, expect):
    """Every panel genome must land on the correct side of the q=0.5 boundary (robust to the exact fit)."""
    q = float(q_from_depth_tail(depth_tail, P))
    if expect == "deep_loss":
        assert q < 0.10
    elif expect == "loss":
        assert q < 0.50
    elif expect == "bearer":
        assert q > 0.50
    else:  # deep_bearer
        assert q > 0.90


def test_q_is_monotonic_increasing_in_depth_tail():
    xs = np.linspace(0, 8, 50)
    qs = q_from_depth_tail(xs, P)
    assert np.all(np.diff(qs) > 0)


# --------------------------------------------------------------------------------------------------
# end-to-end adjudication
# --------------------------------------------------------------------------------------------------
def test_deep_bearer_is_confident():
    margin, p_motif = _synth_genome(-2.0, 5.0)
    r = adjudicate(margin, p_motif, P)
    assert r.status is AdjStatus.ADJUDICATED
    assert r.q > 0.9 and r.q_lo > 0.5
    assert r.depth_tail > 3.39
    assert r.assessable


def test_shallow_loss_reads_loss_side():
    margin, p_motif = _synth_genome(-2.0, 3.0)
    r = adjudicate(margin, p_motif, P)
    assert r.status is AdjStatus.ADJUDICATED
    assert r.q < 0.5 and r.q_hi < 0.5


def test_few_borderline_calls_are_undetermined_not_loss():
    """A handful of borderline calls must widen the CI to straddle 0.5 (UNDETERMINED), never silent loss."""
    rng = np.random.RandomState(3)
    u2 = rng.normal(-2.0, 1.0, 5000)
    q999 = np.percentile(u2, 99.9)
    mad = 1.4826 * np.median(np.abs(u2 - np.median(u2)))
    calls = rng.normal(q999 + 3.0 * mad, 0.3, 4)
    margin = np.concatenate([u2, calls])
    p_motif = 1.0 / (1.0 + np.exp(-(P.platt_a * margin + P.platt_c)))
    r = adjudicate(margin, p_motif, P)
    assert r.status is AdjStatus.UNDETERMINED
    assert r.q_lo <= 0.5 <= r.q_hi


def test_depth_driven_not_count_driven():
    """A CLEAN low-N bearer (two genuinely deep calls, no U2-tail FP contamination) stays confident at N=2:
    depth, not count, drives certainty. (A low U2 background is used so the only calls are the real, deep
    ones — the case the down-sampling result in q_bootstrap.py demonstrated.)"""
    margin, p_motif = _synth_genome(-5.0, 6.0, n_call=2, call_sd=0.05)
    r = adjudicate(margin, p_motif, P)
    assert r.n_call == 2
    assert r.q > 0.9
    assert r.status is AdjStatus.ADJUDICATED


def test_few_deep_calls_buried_in_fps_are_not_silently_confident():
    """The conservative flip-side: 2 deep calls drowned by a shallow U2-tail FP majority must NOT read as a
    confident bearer — the robust median peels the deep minority -> low-q/wide-CI (UNDETERMINED), never a
    confident over-call. (Honest single-genome limit; cross-species conservation is what would rescue it.)"""
    margin, p_motif = _synth_genome(-2.0, 6.0, n_call=2, call_sd=0.05)
    r = adjudicate(margin, p_motif, P)
    confident_bearer = (r.status is AdjStatus.ADJUDICATED and r.q_lo > 0.5)
    assert not confident_bearer


# --------------------------------------------------------------------------------------------------
# operational guards / status codes
# --------------------------------------------------------------------------------------------------
def test_low_n_u2_is_not_assessable():
    """Too few species U2 (below min_u2) -> LOW_N; the file side then defaults to q_eff=1 (no suppression)."""
    margin, p_motif = _synth_genome(-2.0, 5.0, n_u2=100)
    r = adjudicate(margin, p_motif, P)
    assert r.status is AdjStatus.LOW_N
    assert not r.assessable
    assert r.q == 0.0


def test_lowering_min_u2_lets_a_small_genome_self_adjudicate():
    """Option B opt-in: lowering min_u2 lets a small genome assess against its OWN U2 tail (not LOW_N)."""
    margin, p_motif = _synth_genome(-2.0, 5.0, n_u2=100)
    r_default = adjudicate(margin, p_motif, P)
    r_optin = adjudicate(margin, p_motif, AdjudicatorParams(min_u2=50))
    assert r_default.status is AdjStatus.LOW_N            # 100 U2 < default min_u2=200
    assert r_optin.status in (AdjStatus.ADJUDICATED, AdjStatus.UNDETERMINED)   # 100 U2 >= 50 -> own tail
    assert np.isfinite(r_optin.depth_tail)


def test_degenerate_tail_when_u2_scale_collapses():
    u2 = np.full(5000, -2.0)              # MAD == 0
    calls = np.full(20, 5.0)
    margin = np.concatenate([u2, calls])
    p_motif = 1.0 / (1.0 + np.exp(-(P.platt_a * margin + P.platt_c)))
    r = adjudicate(margin, p_motif, P)
    assert r.status is AdjStatus.DEGENERATE_TAIL
    assert not r.assessable


def test_schema_fail_on_shape_mismatch():
    r = adjudicate(np.zeros(10), np.zeros(9), P)
    assert r.status is AdjStatus.SCHEMA_FAIL
    assert not r.assessable


def test_schema_fail_on_all_nan():
    r = adjudicate(np.full(500, np.nan), np.full(500, np.nan), P)
    assert r.status is AdjStatus.SCHEMA_FAIL


def test_nan_rows_are_dropped_not_fatal():
    """A few NaN rows must be filtered, not poison the whole genome."""
    margin, p_motif = _synth_genome(-2.0, 5.0)
    margin[10] = np.nan
    p_motif[20] = np.nan
    r = adjudicate(margin, p_motif, P)
    assert r.status is AdjStatus.ADJUDICATED
    assert r.q > 0.9


# --------------------------------------------------------------------------------------------------
# secondary diagnostic + provenance
# --------------------------------------------------------------------------------------------------
def test_secondary_excess_z_available_for_a_normal_genome():
    margin, p_motif = _synth_genome(-2.0, 5.0)
    r = adjudicate(margin, p_motif, P)
    assert r.secondary_available
    assert np.isfinite(r.excess_z)
    assert np.isfinite(r.z_expmax)
    # a deep bearer should also be a significance outlier (positive excess_z)
    assert r.excess_z > 0


def test_params_version_is_pinned():
    assert P.params_version == ADJUDICATOR_PARAMS_VERSION
    assert "zexcess" in ADJUDICATOR_PARAMS_VERSION  # PRIMARY = z_excess gap gate (depth_tail->q superseded)


def _synthetic(n_u2, n_u12, u12_margin, seed=0):
    from intronIC.scoring.species_adjudicator import p_motif_from_margin
    rng = np.random.RandomState(seed)
    m = np.concatenate([rng.normal(-2.0, 1.0, n_u2), rng.normal(u12_margin, 0.3, n_u12)])
    return m, p_motif_from_margin(m, P)


def test_zexcess_gate_calls_a_numerous_population_bearer():
    """The PRIMARY z_excess gate must call a real (numerous) U12 population DETECTED even when the
    count-blind depth_tail is shallow -- the exact failure depth_tail->q had on divergent bearers."""
    from intronIC.scoring.species_adjudicator import MotifCategory
    m, p = _synthetic(40000, 150, 3.0)
    r = adjudicate(m, p)
    assert np.isfinite(r.z_excess) and r.z_excess >= P.bearer_floor_z
    assert r.motif_category == MotifCategory.DETECTED


def test_zexcess_gate_calls_a_loss():
    from intronIC.scoring.species_adjudicator import MotifCategory
    m, p = _synthetic(40000, 3, 1.5)
    r = adjudicate(m, p)
    assert r.z_excess <= P.loss_ceiling_z
    assert r.motif_category == MotifCategory.NOT_DETECTED


def test_classify_motif_category_gap_is_inconclusive():
    """The empirical gap between the frozen anchors is the abstain (INCONCLUSIVE) zone; non-finite -> UNASSESSABLE.
    cs_p95=nan AND p_gumbel_p95=nan keep both strength paths inert so this exercises the z_excess gate alone."""
    from intronIC.scoring.species_adjudicator import classify_motif_category, MotifCategory
    nan = float("nan")
    assert classify_motif_category(7.0, nan, nan, 0, P) == MotifCategory.DETECTED     # >= bearer_floor 5.50
    assert classify_motif_category(1.0, nan, nan, 0, P) == MotifCategory.NOT_DETECTED
    assert classify_motif_category(3.2, nan, nan, 0, P) == MotifCategory.INCONCLUSIVE   # in the widened gap [2.60, 5.50)
    assert classify_motif_category(5.0, nan, nan, 0, P) == MotifCategory.INCONCLUSIVE   # near-floor uncertain zone (was DETECTED @4.00)
    assert classify_motif_category(nan, nan, nan, 0, P) == MotifCategory.UNASSESSABLE


def test_call_strength_gate_recovers_divergent_bearers():
    """Call-strength CO-FALLBACK (cs_p95): with p_gumbel_p95=nan (primary inert), a genuinely strong POINT
    cs_p95 (>= cs_point_threshold, enough calls) still promotes to DETECTED even when z_excess (count) misses
    it — the null-free fallback path — while a below-threshold cs_p95, or too-few-calls, is NOT promoted."""
    from intronIC.scoring.species_adjudicator import classify_motif_category, MotifCategory, AdjudicatorParams
    P = AdjudicatorParams()
    thr, mc = P.cs_point_threshold, P.cs_min_calls
    nan = float("nan")
    # divergent bearer: low z (would be NOT_DETECTED / gap) but strong cs_p95 with ample calls -> DETECTED
    assert classify_motif_category(1.07, thr + 1.0, nan, 11, P) == MotifCategory.DETECTED
    assert classify_motif_category(3.2, thr + 0.5, nan, 20, P) == MotifCategory.DETECTED   # cs resolves the z-gap
    # relic loss: cs_p95 BELOW threshold (the p95 discounts the lone outlier) -> stays NOT_DETECTED
    assert classify_motif_category(1.92, thr - 1.3, nan, 102, P) == MotifCategory.NOT_DETECTED
    # lone relic: cs above threshold but too FEW calls (p95 == max, outlier-prone) -> gate off -> NOT_DETECTED
    assert classify_motif_category(1.0, thr + 0.5, nan, mc - 1, P) == MotifCategory.NOT_DETECTED
    # disabling the gate reverts to z-only behavior
    Poff = AdjudicatorParams(cs_gate_enabled=False)
    assert classify_motif_category(1.07, thr + 1.0, nan, 11, Poff) == MotifCategory.NOT_DETECTED


def test_p_gumbel_gate_is_the_primary_strength_driver():
    """PRIMARY strength driver (2026-07-03b): a significant per-genome Gumbel tail-prob (p_gumbel_p95 <=
    p_gumbel_threshold, enough calls) promotes to DETECTED even when the null-free cs_p95 co-fallback is BELOW
    its threshold — the per-genome test recovers what the absolute constant would miss. An INsignificant
    p_gumbel_p95 does not promote; too-few-calls gates both paths off."""
    from intronIC.scoring.species_adjudicator import classify_motif_category, MotifCategory, AdjudicatorParams
    P = AdjudicatorParams()
    pt, mc = P.p_gumbel_threshold, P.cs_min_calls
    nan = float("nan")
    # significant per-genome outlier (p below threshold) but cs_p95 below its co-fallback -> DETECTED via p_gumbel
    assert classify_motif_category(3.2, P.cs_point_threshold - 1.0, pt / 2, 9, P) == MotifCategory.DETECTED
    # NON-significant p_gumbel AND cs_p95 below fallback -> not promoted (stays in the gap / loss)
    assert classify_motif_category(3.2, P.cs_point_threshold - 1.0, pt * 20, 9, P) == MotifCategory.INCONCLUSIVE
    assert classify_motif_category(1.0, P.cs_point_threshold - 1.0, pt * 20, 9, P) == MotifCategory.NOT_DETECTED
    # significant p_gumbel but too FEW calls -> gate off -> not promoted
    assert classify_motif_category(1.0, nan, pt / 2, mc - 1, P) == MotifCategory.NOT_DETECTED
    # disabling the strength gate ignores p_gumbel too
    Poff = AdjudicatorParams(cs_gate_enabled=False)
    assert classify_motif_category(1.07, nan, pt / 2, 11, Poff) == MotifCategory.NOT_DETECTED


def test_zero_calls_well_powered_is_not_detected():
    """A genome with an ample U2 background but ZERO strong calls is a definite NOT_DETECTED (the clearest
    loss), NOT UNASSESSABLE -- the plasmodium case (7069 U2, 0 calls). Regression for zexcess_gap_2026-06-30b."""
    from intronIC.scoring.species_adjudicator import MotifCategory, AdjStatus, p_motif_from_margin
    m = np.full(8000, -3.0)                    # all deep in U2 territory -> 0 strong calls, plenty of U2
    p = p_motif_from_margin(m, P)
    r = adjudicate(m, p)
    assert r.n_call == 0 and r.n_u2 >= P.min_u2
    assert r.motif_category == MotifCategory.NOT_DETECTED
    assert r.z_excess == 0.0                   # 0 by construction (no calls)
    assert r.status == AdjStatus.LOW_N         # secondary q-pipeline degenerate, but the gate is decisive


def test_zero_calls_underpowered_is_unassessable():
    """Too few U2 to reference the tail -> genuinely UNASSESSABLE (giardia: 8 U2, 0 calls)."""
    from intronIC.scoring.species_adjudicator import MotifCategory, p_motif_from_margin
    m = np.full(8, -3.0)                        # below min_u2 -> cannot reference its own tail
    p = p_motif_from_margin(m, P)
    r = adjudicate(m, p)
    assert r.motif_category == MotifCategory.UNASSESSABLE


def test_tied_calls_do_not_give_degenerate_ci():
    """Regression: a plain bootstrap-of-median gives a width-0 CI on tied/few calls, silently defeating the
    UNDETERMINED band. The count-aware smoothed bootstrap must give a NON-degenerate CI, widening as the
    call count shrinks (few calls -> wide)."""
    rng = np.random.RandomState(0)
    u2 = rng.normal(-4.8, 1.0, 5000)   # clean U2 (no tail-calls); borderline calls land near the boundary

    def ci_width(call_value, k):
        calls = np.full(k, call_value)
        margin = np.concatenate([u2, calls])
        p_motif = 1.0 / (1.0 + np.exp(-(P.platt_a * margin + P.platt_c)))
        r = adjudicate(margin, p_motif, P)
        return r.q_hi - r.q_lo

    # borderline (depth_tail ~3.4), tied calls -> CI must NOT be degenerate, and must widen for fewer calls
    w3 = ci_width(1.5, 3)
    w60 = ci_width(1.5, 60)
    assert w3 > 0.05                      # not the old width-0 degeneracy
    assert w3 > w60                       # count-aware: fewer calls -> wider CI

    # with the smoothing DISABLED (floor=0), tied calls reproduce the degenerate width-0 CI (documents the bug)
    calls = np.full(3, 1.5); margin = np.concatenate([u2, calls])
    p_motif = 1.0 / (1.0 + np.exp(-(P.platt_a * margin + P.platt_c)))
    r0 = adjudicate(margin, p_motif, AdjudicatorParams(ci_smooth_floor_frac=0.0))
    assert (r0.q_hi - r0.q_lo) < 1e-6     # the degeneracy the fix removes


def test_adjudication_is_order_invariant():
    """The whole adjudication — including the bootstrap q-CI — must be invariant to the input row order.
    The CI bootstrap resamples calls by index, so without sorting the calls the streaming (per-contig) and
    in-memory paths (same call values, different order) produced different P_adj_lo/P_adj_hi. Regression for
    the streaming==in-memory parity break under the raw-feature default."""
    rng = np.random.RandomState(7)
    u2 = rng.normal(-4.8, 1.0, 3000)
    calls = rng.normal(1.6, 0.4, 30)
    margin = np.concatenate([u2, calls])
    p_motif = 1.0 / (1.0 + np.exp(-(P.platt_a * margin + P.platt_c)))
    r1 = adjudicate(margin, p_motif, P)
    perm = rng.permutation(len(margin))
    r2 = adjudicate(margin[perm], p_motif[perm], P)
    assert r1.q == r2.q
    assert r1.q_lo == r2.q_lo and r1.q_hi == r2.q_hi    # the bootstrap CI must match exactly across orderings
    assert r1.depth_tail == r2.depth_tail


def test_bootstrap_is_reproducible():
    margin, p_motif = _synth_genome(-2.0, 5.0)
    r1 = adjudicate(margin, p_motif, P)
    r2 = adjudicate(margin, p_motif, P)
    assert r1.q_lo == r2.q_lo and r1.q_hi == r2.q_hi


def test_params_from_dict_ignores_unknown_keys():
    params = AdjudicatorParams.from_dict({"q_a": 3.0, "not_a_field": 99})
    assert params.q_a == 3.0
    assert not hasattr(params, "not_a_field")


# --------------------------------------------------------------------------------------------------
# file-side post-process (apply_pmotif_adjudication) — the inference hook
# --------------------------------------------------------------------------------------------------
def _tiny_ensemble():
    """A minimal real ensemble matching the bundle structure (m.model.calibrated_classifiers_[i].estimator).

    Trained on a toy 6-feature U12/U2 split so decision_function separates strong motifs from weak ones.
    """
    from types import SimpleNamespace
    from sklearn.svm import SVC
    from sklearn.calibration import CalibratedClassifierCV
    rng = np.random.RandomState(0)
    u2 = np.column_stack([rng.normal(-3, 2, 80), rng.normal(-3, 2, 80), rng.normal(-2, 1.5, 80),
                          rng.normal(-5, 3, 80), rng.normal(2, 2, 80), rng.normal(0, 1, 80)])
    u12 = np.column_stack([rng.normal(14, 1.5, 80), rng.normal(8, 1.5, 80), rng.normal(1.5, 1, 80),
                           rng.normal(-12, 3, 80), rng.normal(7, 1, 80), rng.normal(8, 1, 80)])
    X = np.vstack([u2, u12]); y = np.r_[np.zeros(80), np.ones(80)]
    clf = CalibratedClassifierCV(SVC(C=50, gamma=0.01), cv=3).fit(X, y)
    return [SimpleNamespace(model=clf)]


def _write_score_info(path, n_u2=400, n_u12=20, inject_nan=True):
    import pandas as pd
    rng = np.random.RandomState(1)
    u2 = np.column_stack([rng.normal(-3, 2, n_u2), rng.normal(-3, 2, n_u2), rng.normal(-2, 1.5, n_u2),
                          rng.normal(-5, 3, n_u2), rng.normal(2, 2, n_u2)])
    u12 = np.column_stack([rng.normal(14, 1.5, n_u12), rng.normal(8, 1.5, n_u12), rng.normal(1.5, 1, n_u12),
                           rng.normal(-12, 3, n_u12), rng.normal(7, 1, n_u12)])
    A = np.vstack([u2, u12])
    df = pd.DataFrame({
        "name": [f"i{i}" for i in range(len(A))],
        "5'_raw": A[:, 0], "bp_raw": A[:, 1], "3'_raw": A[:, 2],
        "bp_offset": A[:, 3], "bp_scan_confidence": A[:, 4],
        "type_id": ["u2"] * len(A),
    })
    if inject_nan:
        df.loc[2, "5'_raw"] = np.nan
    df.to_csv(path, sep="\t", index=False)
    return len(A), n_u2


def test_apply_pmotif_adjudication_writes_columns_and_calls(tmp_path):
    from intronIC.scoring.species_adjudicator import apply_pmotif_adjudication
    import pandas as pd
    from intronIC.scoring.species_adjudicator import MotifCategory
    p = tmp_path / "x.score_info.iic"
    # a real U12 *population* (150 strong calls) -> z_excess detects it -> DETECTED. (A handful of
    # calls is no longer a "bearer": the primary gate is a count-based population statistic, not depth alone.)
    n, n_u2 = _write_score_info(p, n_u12=150)
    res = apply_pmotif_adjudication(str(p), _tiny_ensemble(), params=P)
    out = pd.read_csv(p, sep="\t", keep_default_na=False)

    # two-number columns present + row alignment preserved
    for c in ("P_motif", "z_excess", "motif_category", "type_id", "rel_score"):
        assert c in out.columns
    assert len(out) == n
    assert res.motif_category is MotifCategory.DETECTED   # 400 U2 + a real U12 population -> motif population DETECTED
    assert np.isfinite(res.z_excess) and res.z_excess >= P.bearer_floor_z
    assert res.n_adj_called == int((out["type_id"] == "u12").sum())
    # ungated motif-call count: for a non-NOT_DETECTED genome the gate doesn't suppress, so it equals the
    # gated call count == #(P_motif>=0.5). (For NOT_DETECTED genomes it would EXCEED n_adj_called==0.)
    pm_all = pd.to_numeric(out["P_motif"], errors="coerce").to_numpy()
    assert res.n_motif_called == int(np.nansum(pm_all >= 0.5)) == res.n_adj_called
    # motif_category is constant within the run (one genome)
    assert set(out["motif_category"]) == {MotifCategory.DETECTED.value}

    pm = pd.to_numeric(out["P_motif"], errors="coerce").to_numpy()
    padj = pd.to_numeric(out["P_adj"], errors="coerce").to_numpy()
    fin = np.isfinite(pm) & np.isfinite(padj)
    # legacy P_adj == q_eff * P_motif with the BINARY gate (q_eff==1 for a non-loss species)
    assert res.q == 1.0 and np.allclose(padj[fin], pm[fin], atol=1e-6)
    # the U12 population (last n=150) should be the called set; U2 bulk mostly not called
    assert (out["type_id"].to_numpy()[-150:] == "u12").sum() >= 140
    assert (out["type_id"].to_numpy()[:n_u2] == "u12").mean() < 0.1
    # injected NaN row stays unscored / uncalled
    assert out["P_motif"][2] == "nan" and out["type_id"][2] == "u2"


def test_apply_pmotif_adjudication_low_n_falls_back_to_pmotif(tmp_path):
    """Below min_u2 the species layer is inconclusive -> LOW_N -> q_eff=1 (P_adj==P_motif), never silent
    suppression (option A, the safe low-N default)."""
    from intronIC.scoring.species_adjudicator import apply_pmotif_adjudication
    import pandas as pd
    p = tmp_path / "small.score_info.iic"
    _write_score_info(p, n_u2=80, n_u12=3, inject_nan=False)   # < min_u2=200 -> LOW_N
    res = apply_pmotif_adjudication(str(p), _tiny_ensemble(), params=P)
    out = pd.read_csv(p, sep="\t", keep_default_na=False)
    assert res.status is AdjStatus.LOW_N
    assert float(out["q"].iloc[0]) == 1.0
    pm = pd.to_numeric(out["P_motif"], errors="coerce").to_numpy()
    padj = pd.to_numeric(out["P_adj"], errors="coerce").to_numpy()
    fin = np.isfinite(pm) & np.isfinite(padj)
    assert np.allclose(padj[fin], pm[fin], atol=1e-9)
