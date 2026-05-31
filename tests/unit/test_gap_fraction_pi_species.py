"""Step 3 (gate-gapfrac-graded): gap_fraction bootstrap-UCL → π_species.

The species prior replaces the unstable KDE valley_depth. It keys on the OPTIMISTIC
bootstrap upper-confidence-limit (UCL) of gap_fraction, so the low-N over-kill guard is
baked into the statistic itself:
  * compute_adjusted_score maps the UCL → π_species (graded; -inf suppresses; NaN no-op),
  * compute_valley_depth lifts the UCL above the point estimate when the bootstrap is wide
    (low-N) — protecting a real bearer whose point estimate is noisy — while a deep-overlap
    cluster keeps a low UCL and is suppressed.

Run directly (pytest may be blocked by a bundle numpy mismatch):
  PYTHONPATH=<wt>/src <env>/bin/python3 tests/unit/test_gap_fraction_pi_species.py
"""
import math
import numpy as np

from intronIC.scoring.cluster_validation import (
    compute_adjusted_score,
    compute_adjusted_scores_batch,
    compute_valley_depth,
    DEFAULT_GAP_FRACTION_MIDPOINT,
)

SVM = 95.0
SIGMA = 5.0


def _sigma_only(svm, sigma, k_sigma=3.0):
    logit = math.log((svm / 100.0) / (1 - svm / 100.0))
    return 100.0 / (1.0 + math.exp(-(logit - k_sigma * sigma / 100.0)))


def test_pi_species_regimes():
    hi = compute_adjusted_score(SVM, 0.8, SIGMA)               # clear separation (UCL high)
    lo = compute_adjusted_score(SVM, 0.0, SIGMA)               # poor separation (UCL at 0)
    lifted = compute_adjusted_score(SVM, 0.6, SIGMA)           # uncertainty lifted UCL > midpoint
    neg = compute_adjusted_score(SVM, float("-inf"), SIGMA)    # confident no-separation
    nan_adj = compute_adjusted_score(SVM, float("nan"), SIGMA)  # can't assess

    # Clear separation barely touched (only the σ penalty); poor separation strongly suppressed.
    assert abs(hi - SVM) < 5, hi
    assert lo < 15, lo
    assert hi > lo + 50, (hi, lo)

    # A UCL above the suppression midpoint (e.g. uncertainty lifted it there) ⇒ ~no suppression.
    assert lifted > 85, lifted
    assert DEFAULT_GAP_FRACTION_MIDPOINT == 0.35  # guard: the 0.6 case is above midpoint

    # Structural no-separation (-inf) → floor → near zero.
    assert neg < 5, neg

    # Can't assess (NaN) → no prior correction at all → exactly the σ-only score.
    assert abs(nan_adj - _sigma_only(SVM, SIGMA)) < 1e-6, (nan_adj, _sigma_only(SVM, SIGMA))


def test_monotonic_in_ucl():
    xs = [compute_adjusted_score(SVM, u, SIGMA)
          for u in (-5.0, -0.5, 0.0, 0.2, 0.35, 0.5, 0.8, 2.0)]
    assert all(b >= a - 1e-9 for a, b in zip(xs, xs[1:])), xs


def test_batch_matches_scalar():
    arr = compute_adjusted_scores_batch(np.array([SVM, SVM]), 0.0, np.array([SIGMA, SIGMA]))
    scalar = compute_adjusted_score(SVM, 0.0, SIGMA)
    assert abs(float(arr[0]) - scalar) < 1e-9, (arr[0], scalar)


def test_valley_depth_deterministic_and_order_invariant():
    rng = np.random.default_rng(0)
    u2 = rng.normal(0.0, 1.0, size=(2000, 3))
    u12 = rng.normal(3.0, 0.6, size=(40, 3))   # clearly separated

    r1 = compute_valley_depth(u12, u2)
    r2 = compute_valley_depth(u12, u2)
    assert r1["gap_fraction"] == r2["gap_fraction"]          # RNG-free point estimate
    assert r1["gap_fraction_ucl"] == r2["gap_fraction_ucl"]  # deterministic bootstrap

    perm = rng.permutation(len(u12))
    r3 = compute_valley_depth(u12[perm], u2)
    assert abs(r1["gap_fraction"] - r3["gap_fraction"]) < 1e-9
    assert abs(r1["gap_fraction_ucl"] - r3["gap_fraction_ucl"]) < 1e-9

    assert r1["gap_fraction"] > 0 and math.isfinite(r1["gap_fraction"])
    assert math.isfinite(r1["gap_fraction_ucl"])


def test_low_n_lifts_ucl_more_than_high_n():
    # Same separation, different N: the smaller sample has a wider bootstrap, so its UCL is
    # lifted further above the point estimate (the low-N protection).
    rng = np.random.default_rng(7)
    u2 = rng.normal(0.0, 1.0, size=(2000, 3))
    r_lo = compute_valley_depth(rng.normal(2.0, 1.0, size=(12, 3)), u2)
    r_hi = compute_valley_depth(rng.normal(2.0, 1.0, size=(500, 3)), u2)
    lift_lo = r_lo["gap_fraction_ucl"] - r_lo["gap_fraction"]
    lift_hi = r_hi["gap_fraction_ucl"] - r_hi["gap_fraction"]
    assert lift_lo > lift_hi, (lift_lo, lift_hi)
    assert lift_lo >= -1e-9, lift_lo   # UCL is optimistic (≥ point estimate)


def test_deep_overlap_ucl_suppresses():
    # U12 calls essentially drawn from U2 → even the optimistic UCL is below the midpoint,
    # so the species is suppressed (FP-suppression preserved despite ratio instability).
    rng = np.random.default_rng(9)
    u2 = rng.normal(0.0, 1.0, size=(2000, 3))
    u12 = rng.normal(0.1, 1.0, size=(60, 3))
    r = compute_valley_depth(u12, u2)
    ucl = r["gap_fraction_ucl"]
    assert ucl < DEFAULT_GAP_FRACTION_MIDPOINT, ucl
    assert compute_adjusted_score(95.0, ucl, 0.0) < 30, (ucl,)


def test_valley_depth_insufficient_is_nan():
    rng = np.random.default_rng(1)
    u2 = rng.normal(0.0, 1.0, size=(2000, 3))
    r = compute_valley_depth(rng.normal(3.0, 0.6, size=(2, 3)), u2)  # < 3 U12 points
    assert math.isnan(r["gap_fraction"])
    assert math.isnan(r["gap_fraction_ucl"])


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"ALL {len(fns)} CHECKS PASSED")
