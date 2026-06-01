"""Single-sequence / few-intron robustness (direct `-q` input is a first-class mode).

The whole scoring stack is population-based; this locks the graceful-degradation cascade at the
scoring-function layer (the full CLI can't run here due to the pickle/numpy-bundle block, but these
are the functions a 1–few-intron input actually hits). Asserts: no crash, the right fallback fires,
and the species-level terms (this workstream's gate + π_species) no-op cleanly at n→1.

Run: PYTHONPATH=<wt>/src <env>/bin/python3 tests/unit/test_low_n_robustness.py
"""
import math
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from intronIC.scoring.mode_separation import (
    fit_mode_separation, candidate_weight_from_svm, evaluate_gate,
)
from intronIC.scoring.cluster_validation import (
    compute_valley_depth, validate_u12_cluster, compute_adjusted_score,
)
from intronIC.classification.mode_sep_pipeline import apply_continuous_per_intron_discount


def test_fit_mode_separation_tiny_no_crash():
    """The 2-component fit runs upstream of the gate's n_floor guard, so it must not crash on
    n=1 (incl. an all-low-svm intron → all candidate weight ≈ 0)."""
    for svm_vals in ([99.0], [30.0], [99.0, 95.0], [99.0, 50.0, 20.0]):
        raw = np.array(svm_vals)
        w = candidate_weight_from_svm(np.array(svm_vals))
        s = fit_mode_separation(raw, w)  # must not raise
        assert math.isfinite(s.n_eff_candidates), svm_vals
        # mu_* may be NaN when weights collapse (all-low-svm); that's acceptable BECAUSE the gate
        # bails at below_n_floor first (NaN separation also can't pass the > checks). What matters
        # here is no exception.
        print(f"  fit n={len(svm_vals)} svm={svm_vals}: "
              f"n_eff={s.n_eff_candidates:.3f} sep={s.separation}")


def test_gate_below_n_floor_short_circuits():
    """Single candidate → n_eff < 5 → below_n_floor, returned BEFORE gap_fraction/csig are computed."""
    u12 = np.array([[3.0, 1.0, 1.0]])
    u2 = np.tile([0.0, 0.0, 0.0], (12, 1)).astype(float)
    g = evaluate_gate(u12, u2, mu_u12_5p_raw=15.7, n_eff_candidates=1.0,
                      separation_5p=2.0, valley_depth_fn=compute_valley_depth)
    assert not g.passes and g.reason == "below_n_floor", g
    assert g.gap_fraction is None and g.centroid_sigma is None  # never computed


def test_compute_valley_depth_tiny_insufficient():
    """< 3 U12 points → insufficient early-return: gap_fraction + UCL are NaN (π_species no-op)."""
    u2 = np.random.default_rng(0).normal(0, 1, (50, 3))
    for n in (1, 2):
        u12 = np.random.default_rng(1).normal(3.0, 0.5, (n, 3))
        r = compute_valley_depth(u12, u2)
        assert math.isnan(r["gap_fraction"]), (n, r["gap_fraction"])
        assert math.isnan(r["gap_fraction_ucl"]), (n, r["gap_fraction_ucl"])


def test_validate_u12_cluster_single_intron():
    """A single given intron → insufficient regime → gap_fraction_ucl NaN ⇒ fallback π_species no-ops."""
    r = validate_u12_cluster(
        five_z_scores=np.array([5.0]), bp_z_scores=np.array([3.0]),
        svm_scores=np.array([99.0]), type_ids=np.array(["u12"]),
        three_z_scores=np.array([1.0]),
    )
    assert r["regime"] == "insufficient", r["regime"]
    assert math.isnan(r["gap_fraction_ucl"]), r["gap_fraction_ucl"]


def test_compute_adjusted_score_noop_when_unassessable():
    """No population ⇒ UCL NaN ⇒ no species prior correction ⇒ exactly the σ-only score."""
    a = compute_adjusted_score(95.0, float("nan"), 5.0)
    sig_only = 100.0 / (1.0 + math.exp(-(math.log(0.95 / 0.05) - 3.0 * 0.05)))
    assert abs(a - sig_only) < 1e-6, (a, sig_only)


def test_continuous_discount_single_row():
    """Per-intron discount works on a 1-row score_info (the single-sequence FP guard)."""
    d = tempfile.mkdtemp()
    try:
        p = Path(d) / "one.score_info.iic"
        cols = ["name", "svm_score", "rel_score", "adjusted_score", "ensemble_sigma",
                "5'_raw", "bp_raw", "3'_raw", "first_pass_svm"]
        # high svm + weak motif (raw_sum ≈ 1) → a clear overcall the per-intron discount should bite
        row = ["solo_intron", 99.0, 9.0, "NA", 5.0, 1.0, 0.0, 0.0, 99.0]
        with open(p, "w") as f:
            f.write("\t".join(cols) + "\n")
            f.write("\t".join(str(x) for x in row) + "\n")
        summ = apply_continuous_per_intron_discount(p, threshold=90.0, input_column="svm_score")
        df = pd.read_csv(p, sep="\t")
        assert len(df) == 1
        adj = float(df["adjusted_score"].iloc[0])
        assert math.isfinite(adj) and adj < 90.0, adj  # overcall suppressed below the call line
        assert summ["n_called_pre_discount"] == 1 and summ["n_called_post_discount"] == 0
        print(f"  discount single-row: svm=99 weak-motif → adjusted={adj:.2f}")
    finally:
        shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"ALL {len(fns)} LOW-N ROBUSTNESS CHECKS PASSED")
