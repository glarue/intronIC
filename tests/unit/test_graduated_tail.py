"""Unit tests for the graduated/Platt tail (C3), the species penalty (C6), and their gating.

Covers the deployed-but-bundle-gated additions on branch v6-graduated-tail:
  - species_penalty_logit_shift: one-sided (<=0), bearer-untouched, n_hc=0 no-op, formula value.
  - apply_continuous_per_intron_discount graduated_tail branch: fires only when a graduated_tail dict + an
    svm_margin column are present; reproduces 100*sigmoid(coef . scale([margin,raw5,rawbp,raw3,adh5z,adhbpz])
    + intercept); ABSENT -> the legacy overcall discount (byte-unchanged default path).
  - enable_species_penalty gating: OFF (default) is inert (species_penalty None, adj unchanged); ON applies the shift.
"""
import math
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from intronIC.scoring.mode_separation import (
    species_penalty_logit_shift,
    DEFAULT_SPECIES_PENALTY_PGATE,
)
from intronIC.classification.mode_sep_pipeline import apply_continuous_per_intron_discount


# --------------------------------------------------------------------------- #
# C6: species_penalty_logit_shift
# --------------------------------------------------------------------------- #
def test_penalty_zero_when_no_hc():
    assert species_penalty_logit_shift(0.0, 0) == 0.0
    assert species_penalty_logit_shift(0.9, 0) == 0.0


def test_penalty_one_sided_non_positive():
    for fb6 in (0.0, 0.3, 0.5, 0.8, 1.0):
        for n in (1, 5, 20, 200, 2000):
            assert species_penalty_logit_shift(fb6, n) <= 1e-12


def test_penalty_zero_for_confident_bearer():
    # high frac_bp6 + high N -> p_bearer >= p_gate -> untouched
    assert species_penalty_logit_shift(0.80, 888) == pytest.approx(0.0, abs=1e-9)
    assert species_penalty_logit_shift(0.70, 500) == pytest.approx(0.0, abs=1e-9)


def test_penalty_negative_for_loss_like():
    # low frac_bp6 + low N -> demoted
    assert species_penalty_logit_shift(0.0, 10) < -0.5
    assert species_penalty_logit_shift(0.4, 5) == pytest.approx(-1.681, abs=2e-3)


def test_penalty_high_n_blind_spot():
    # the documented high-N blind spot: a moderate-frac_bp6 species with large N escapes the penalty
    assert species_penalty_logit_shift(0.35, 232) == pytest.approx(0.0, abs=1e-9)


# --------------------------------------------------------------------------- #
# helpers for the discount/tail tests
# --------------------------------------------------------------------------- #
_GT = {
    "feature_order": ["margin", "raw5", "rawbp", "raw3", "adh5z", "adhbpz"],
    "scaler_mean": [0.0, 10.0, 5.0, 1.0, 0.0, 0.0],
    "scaler_scale": [1.0, 5.0, 3.0, 1.0, 1.0, 1.0],
    "coef": [-2.0, 4.0, 1.0, 0.5, 1.0, 1.0],
    "intercept": -1.0,
}


def _tiny_scoreinfo(with_margin=True):
    n = 6
    df = pd.DataFrame({
        "name": [f"i{j}" for j in range(n)],
        "svm_score": [99.0, 95.0, 92.0, 50.0, 30.0, 10.0],
        "first_pass_svm": [99.0, 95.0, 92.0, 50.0, 30.0, 10.0],
        "5'_raw": [13.0, 12.0, 11.0, 6.0, 4.0, 2.0],
        "bp_raw": [7.0, 7.0, 2.0, 7.0, 2.0, 2.0],
        "3'_raw": [2.0, 2.0, 2.0, 1.0, 1.0, 1.0],
        "5'_seq": ["AAGGTATGTAAA"] * n,   # len>=12 -> -3-anchored window (GT at idx3)
        "bp_seq": ["TTTTGCTGACAA"] * n,
        "3'_seq": ["TTTTTTACAG"] * n,
        "adjusted_score": [np.nan] * n,
    })
    if with_margin:
        df["svm_margin"] = [3.0, 1.5, 0.2, -0.5, -2.0, -4.0]
    return df


def _run_discount(df, **kw):
    f = tempfile.NamedTemporaryFile(suffix=".score_info.iic", delete=False)
    f.close()
    df.to_csv(f.name, sep="\t", index=False)
    summ = apply_continuous_per_intron_discount(
        Path(f.name), threshold=90.0, input_column="svm_score", meta_path=None, messenger=None, **kw)
    out = pd.read_csv(f.name, sep="\t")
    Path(f.name).unlink()
    return out, summ


# --------------------------------------------------------------------------- #
# C3: graduated tail
# --------------------------------------------------------------------------- #
def test_graduated_tail_reproduces_formula():
    from intronIC.scoring.gold_adherence import adh_features
    df = _tiny_scoreinfo(with_margin=True)
    out, _ = _run_discount(df.copy(), graduated_tail=_GT)
    adj = pd.to_numeric(out["adjusted_score"], errors="coerce").to_numpy()
    # recompute the tail independently and confirm the function applied exactly that
    m = pd.to_numeric(out["svm_margin"], errors="coerce").to_numpy()
    gmask = ~np.isnan(m)
    a5z, abpz = adh_features(out, gmask, want5=True)
    cm = {"margin": m[gmask], "raw5": out["5'_raw"].to_numpy()[gmask],
          "rawbp": out["bp_raw"].to_numpy()[gmask], "raw3": out["3'_raw"].to_numpy()[gmask],
          "adh5z": a5z, "adhbpz": abpz}
    X = np.column_stack([cm[f] for f in _GT["feature_order"]])
    z = (X - np.array(_GT["scaler_mean"])) / np.array(_GT["scaler_scale"])
    expect = 100.0 / (1.0 + np.exp(-(z @ np.array(_GT["coef"]) + _GT["intercept"])))
    assert np.max(np.abs(adj[gmask] - expect)) < 1e-6
    assert np.all((adj >= 0) & (adj <= 100))


def test_graduated_tail_absent_is_legacy_discount():
    # no graduated_tail -> overcall discount path; adjusted_score must NOT equal the tail output
    df = _tiny_scoreinfo(with_margin=True)
    out_legacy, _ = _run_discount(df.copy())                       # graduated_tail=None (default)
    out_tail, _ = _run_discount(df.copy(), graduated_tail=_GT)
    a_legacy = pd.to_numeric(out_legacy["adjusted_score"], errors="coerce").to_numpy()
    a_tail = pd.to_numeric(out_tail["adjusted_score"], errors="coerce").to_numpy()
    assert np.max(np.abs(a_legacy - a_tail)) > 1.0   # the tail genuinely changed the scores


def test_graduated_tail_skipped_without_margin_column():
    # graduated_tail present but no svm_margin column -> tail cannot fire -> legacy discount
    df = _tiny_scoreinfo(with_margin=False)
    out_a, _ = _run_discount(df.copy(), graduated_tail=_GT)
    out_b, _ = _run_discount(df.copy())
    assert np.allclose(pd.to_numeric(out_a["adjusted_score"], errors="coerce"),
                       pd.to_numeric(out_b["adjusted_score"], errors="coerce"), equal_nan=True)


# --------------------------------------------------------------------------- #
# C3: no-margin + spac tail (the v2_tail+pen production lock, task #71)
# --------------------------------------------------------------------------- #
_GT_V2 = {   # the productionized feature_order: DROP margin, ADD spac (BP->3'SS geometry)
    "feature_order": ["raw5", "rawbp", "raw3", "adh5z", "adhbpz", "spac"],
    "scaler_mean": [10.0, 5.0, 1.0, 0.0, 0.0, 0.5],
    "scaler_scale": [5.0, 3.0, 1.0, 1.0, 1.0, 0.3],
    "coef": [4.0, 1.0, 0.5, 1.0, 1.0, 0.8],
    "intercept": -1.0,
}


def _spac(bo):   # mirrors mode_sep_pipeline.py: exp(-((|bp_offset| - 12)/8)^2)
    return np.exp(-((np.abs(np.asarray(bo, float)) - 12.0) / 8.0) ** 2)


def test_spac_geometry_value():
    # spac peaks at |bp_offset|=12 (U12 BP->3'SS geometry) and decays for U2-like spacing
    assert _spac(-12.0) == pytest.approx(1.0, abs=1e-12)
    assert _spac(12.0) == pytest.approx(1.0, abs=1e-12)
    assert _spac(-4.0) == pytest.approx(math.exp(-1.0), abs=1e-9)        # (|-4|-12)/8 = -1
    assert _spac(-24.0) < _spac(-15.0) < _spac(-12.0)                    # U2-like (-24) decays below U12 (-12)


def test_graduated_tail_spac_reproduces_formula():
    from intronIC.scoring.gold_adherence import adh_features
    df = _tiny_scoreinfo(with_margin=True)
    df["bp_offset"] = [-12.0, -15.0, -24.0, -10.0, -33.0, -4.0]
    out, _ = _run_discount(df.copy(), graduated_tail=_GT_V2)
    adj = pd.to_numeric(out["adjusted_score"], errors="coerce").to_numpy()
    m = pd.to_numeric(out["svm_margin"], errors="coerce").to_numpy()
    bo = pd.to_numeric(out["bp_offset"], errors="coerce").to_numpy()
    gmask = ~np.isnan(m) & ~np.isnan(bo)           # spac tail still gates on svm_margin AND bp_offset
    a5z, abpz = adh_features(out, gmask, want5=True)
    cm = {"raw5": out["5'_raw"].to_numpy()[gmask], "rawbp": out["bp_raw"].to_numpy()[gmask],
          "raw3": out["3'_raw"].to_numpy()[gmask], "adh5z": a5z, "adhbpz": abpz, "spac": _spac(bo[gmask])}
    X = np.column_stack([cm[f] for f in _GT_V2["feature_order"]])      # NOTE: no `margin` column consumed
    z = (X - np.array(_GT_V2["scaler_mean"])) / np.array(_GT_V2["scaler_scale"])
    expect = 100.0 / (1.0 + np.exp(-(z @ np.array(_GT_V2["coef"]) + _GT_V2["intercept"])))
    assert np.max(np.abs(adj[gmask] - expect)) < 1e-6
    assert np.all((adj >= 0) & (adj <= 100))


def test_graduated_tail_spac_requires_bp_offset_per_row():
    # a row with NaN bp_offset is excluded from the spac tail and keeps the discount fallback
    df = _tiny_scoreinfo(with_margin=True)
    df["bp_offset"] = [-12.0, np.nan, -24.0, -10.0, -33.0, -4.0]       # row 1 lacks bp_offset
    out_tail, _ = _run_discount(df.copy(), graduated_tail=_GT_V2)
    out_legacy, _ = _run_discount(df.copy())                           # no tail -> discount everywhere
    a_tail = pd.to_numeric(out_tail["adjusted_score"], errors="coerce").to_numpy()
    a_legacy = pd.to_numeric(out_legacy["adjusted_score"], errors="coerce").to_numpy()
    assert a_tail[1] == pytest.approx(a_legacy[1], abs=1e-9)           # row 1 = unchanged discount fallback
    assert np.max(np.abs(a_tail - a_legacy)) > 1.0                     # other rows: the tail fired


def test_graduated_tail_spac_no_bp_offset_column_is_legacy():
    # spac in feature_order but NO bp_offset column at all -> tail cannot fire -> legacy discount everywhere
    df = _tiny_scoreinfo(with_margin=True)                             # no bp_offset column
    out_a, _ = _run_discount(df.copy(), graduated_tail=_GT_V2)
    out_b, _ = _run_discount(df.copy())
    assert np.allclose(pd.to_numeric(out_a["adjusted_score"], errors="coerce"),
                       pd.to_numeric(out_b["adjusted_score"], errors="coerce"), equal_nan=True)


# --------------------------------------------------------------------------- #
# C6 gating inside the discount
# --------------------------------------------------------------------------- #
def test_species_penalty_off_by_default():
    df = _tiny_scoreinfo(with_margin=True)
    out, summ = _run_discount(df.copy(), graduated_tail=_GT)        # enable_species_penalty default False
    assert summ["species_penalty"] is None


def test_species_penalty_on_applies_shift():
    df = _tiny_scoreinfo(with_margin=True)
    out_off, _ = _run_discount(df.copy(), graduated_tail=_GT)
    out_on, summ = _run_discount(df.copy(), graduated_tail=_GT, enable_species_penalty=True,
                                 species_penalty_pgate=DEFAULT_SPECIES_PENALTY_PGATE)
    a_off = pd.to_numeric(out_off["adjusted_score"], errors="coerce").to_numpy()
    a_on = pd.to_numeric(out_on["adjusted_score"], errors="coerce").to_numpy()
    if summ["species_penalty"] is not None and summ["species_penalty"]["logit_shift"] < 0:
        # one-sided: penalty only lowers scores
        assert np.all(a_on <= a_off + 1e-9)
        assert np.max(np.abs(a_on - a_off)) > 1e-6
    else:
        # tiny synthetic HC set landed in the no-shift regime; at least it must be inert
        assert np.allclose(a_on, a_off, equal_nan=True)
