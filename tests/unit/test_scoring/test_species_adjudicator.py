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
    margin, p_motif = _synth_genome(-2.0, 5.0, n_u2=100)
    r = adjudicate(margin, p_motif, P)
    assert r.status is AdjStatus.LOW_N
    assert not r.assessable
    assert r.q == 0.0


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
    assert "depth_tail" in ADJUDICATOR_PARAMS_VERSION


def test_bootstrap_is_reproducible():
    margin, p_motif = _synth_genome(-2.0, 5.0)
    r1 = adjudicate(margin, p_motif, P)
    r2 = adjudicate(margin, p_motif, P)
    assert r1.q_lo == r2.q_lo and r1.q_hi == r2.q_hi


def test_params_from_dict_ignores_unknown_keys():
    params = AdjudicatorParams.from_dict({"q_a": 3.0, "not_a_field": 99})
    assert params.q_a == 3.0
    assert not hasattr(params, "not_a_field")
