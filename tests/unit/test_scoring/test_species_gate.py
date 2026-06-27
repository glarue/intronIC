"""Unit tests for the raw_gated continuous species gate + support2_raw."""
import numpy as np
import pytest

from intronIC.core.intron import IntronScores
from intronIC.scoring.species_gate import (
    SpeciesGateParams, apply_species_gate, compute_gate_signals, _ramp,
)


# ---- support2_raw -------------------------------------------------------------------------
def test_support2_raw_all_ends_strong():
    s = IntronScores(five_raw_score=14.0, bp_raw_score=8.0, three_raw_score=2.0)
    # clipped-at-0 sorted = [2, 8, 14]; second-largest = 8
    assert s.support2_raw == pytest.approx(8.0)


def test_support2_raw_one_end_strong_is_zero():
    # strong 5' only (BP and 3' negative) -> second-largest of clipped = 0
    s = IntronScores(five_raw_score=15.0, bp_raw_score=-3.0, three_raw_score=-1.0)
    assert s.support2_raw == pytest.approx(0.0)


def test_support2_raw_missing_returns_none():
    assert IntronScores(five_raw_score=14.0, bp_raw_score=None, three_raw_score=2.0).support2_raw is None


def test_support2_raw_independent_of_zscores():
    # raw analog must not depend on z fields being present
    s = IntronScores(five_raw_score=10.0, bp_raw_score=6.0, three_raw_score=1.0,
                     five_z_score=None, bp_z_score=None, three_z_score=None)
    assert s.support2_raw == pytest.approx(6.0)


# ---- ramp ---------------------------------------------------------------------------------
def test_ramp_clamps():
    assert _ramp(0, 4, 50) == 0.0
    assert _ramp(100, 4, 50) == 1.0
    assert _ramp(27, 4, 50) == pytest.approx((27 - 4) / (50 - 4))


# ---- gate: loss vs bearer -----------------------------------------------------------------
def test_loss_species_gate_is_strict():
    # a loss species: motif-empty pool, no all-ends-strong core, a couple of stray one-end FPs
    rng = np.random.RandomState(0)
    n = 5000
    p = rng.uniform(0, 0.3, n)            # mostly low P
    s2 = np.zeros(n)                      # one-end-strong everywhere -> support2_raw ~ 0
    # 2 stray high-P introns but still one-end (support2 low) -> must NOT open the gate
    p[:2] = 0.999
    res = apply_species_gate(p, s2, SpeciesGateParams())
    assert res.g == pytest.approx(0.0)
    assert res.tau == pytest.approx(SpeciesGateParams().tau_high)


def test_loss_with_few_stray_cores_still_strict():
    # the hardening case: a handful of all-ends-strong high-P introns, but they don't form a mode
    rng = np.random.RandomState(1)
    n = 5000
    p = rng.uniform(0, 0.3, n); s2 = np.zeros(n)
    p[:3] = 0.99; s2[:3] = 7.0           # 3 stray strong cores
    res = apply_species_gate(p, s2, SpeciesGateParams())
    # core=3 < core_lo=4 -> g_core=0; and top-K gap dominated by the 0s -> g_gap=0 -> g=0
    assert res.g == pytest.approx(0.0)
    assert res.n_called <= 3


def test_bearer_species_gate_relaxes():
    # a bearer: a clear all-ends-strong U12 cluster on top of a U2 background
    rng = np.random.RandomState(2)
    n = 6000
    p = rng.uniform(0, 0.3, n); s2 = np.zeros(n)
    core_n = 120
    p[:core_n] = rng.uniform(0.9, 1.0, core_n)     # high-confidence
    s2[:core_n] = rng.uniform(6.0, 10.0, core_n)   # all-ends-strong
    res = apply_species_gate(p, s2, SpeciesGateParams())
    assert res.core >= 100
    assert res.gap >= 6.0
    assert res.g == pytest.approx(1.0, abs=1e-6)
    assert res.tau == pytest.approx(SpeciesGateParams().tau_low)
    # the cluster gets called at the relaxed threshold
    assert res.n_called >= core_n


def test_gate_signals_ignore_nan_probs():
    p = np.array([np.nan, np.nan, 0.99, 0.99]); s2 = np.array([0.0, 0.0, 7.0, 7.0])
    core, gap, g = compute_gate_signals(p, s2, SpeciesGateParams())
    assert core == 2  # NaN-P introns excluded


def test_params_from_dict_roundtrip():
    p = SpeciesGateParams.from_dict({"tau_high": 0.99, "tau_low": 0.8, "unknown_key": 1})
    assert p.tau_high == 0.99 and p.tau_low == 0.8
    assert SpeciesGateParams.from_dict(None) == SpeciesGateParams()


def test_shape_mismatch_raises():
    with pytest.raises(ValueError):
        apply_species_gate(np.zeros(3), np.zeros(4))
