"""Path C' species-level scaler gate for v4_aug deployment.

The gate decides per species whether to use adaptive, frozen, or strict
scoring based on the call-set asymmetry observed when scoring under BOTH
adaptive and frozen z-scaling. The rule is empirically derived from a
45-species reference panel; see docs in
``/mnt/data/u12/ipa/training_data/multispecies_v23/ARCHITECTURE_REASSESSMENT.md``
for the full analysis.

This module is pure logic — it operates on numpy arrays of probabilities
and returns routing decisions. The actual integration into the classify
pipeline (running both adaptive and frozen scoring) happens in
``cli/main.py``.

Currently not wired into the pipeline; Phase 1 plumbing only.
"""
from __future__ import annotations

from typing import Literal, TypedDict

import numpy as np


Route = Literal["adaptive", "frozen", "strict"]


class GateStats(TypedDict):
    """Per-species asymmetry counts used by the gate decision."""

    n_total: int
    n_both: int
    n_only_a: int
    n_only_f: int
    frac_only_a: float
    frac_only_f: float
    frac_both: float


def compute_asymmetry(
    scores_adaptive: np.ndarray,
    scores_frozen: np.ndarray,
    threshold: float = 90.0,
) -> GateStats:
    """Compute per-species call-set asymmetry at the given threshold.

    Defines four mutually exclusive sets of high-confidence calls:
      - both:    score_a >= thr AND score_f >= thr
      - only_a:  score_a >= thr AND score_f < thr
      - only_f:  score_a < thr  AND score_f >= thr
      - neither: both < thr (not counted; included implicitly by n_total below)

    n_total is the union of the three call sets (introns called by *some*
    mode); frac_* are normalized to n_total.

    Args:
        scores_adaptive: Per-intron probability * 100 under adaptive z.
        scores_frozen:   Per-intron probability * 100 under frozen z.
        threshold:       The call threshold (default 90 — production).

    Returns:
        GateStats dict.
    """
    if scores_adaptive.shape != scores_frozen.shape:
        raise ValueError(
            f"Adaptive/frozen score arrays must have same shape; "
            f"got {scores_adaptive.shape} vs {scores_frozen.shape}"
        )
    a = scores_adaptive >= threshold
    f = scores_frozen >= threshold
    n_both = int(np.sum(a & f))
    n_only_a = int(np.sum(a & ~f))
    n_only_f = int(np.sum(~a & f))
    n_total = n_both + n_only_a + n_only_f
    denom = max(n_total, 1)
    return {
        "n_total": n_total,
        "n_both": n_both,
        "n_only_a": n_only_a,
        "n_only_f": n_only_f,
        "frac_only_a": n_only_a / denom,
        "frac_only_f": n_only_f / denom,
        "frac_both": n_both / denom,
    }


def decide_route(
    stats: GateStats,
    only_f_thr: float = 0.10,
    only_a_thr: float = 0.50,
    n_both_max: int = 5,
    n_total_min: int = 5,
) -> Route:
    """Decide which scaler route applies for this species.

    Rule (revised May 14 after 45-species pressure test):
      1. frozen: frac_only_f >= only_f_thr AND n_total >= n_total_min
         (the species shows adaptive-demotion of borderline U12s; frozen
         recovers them).
      2. strict: frac_only_a >= only_a_thr AND n_both <= n_both_max
         (no n_total floor — adaptive over-calls into a U12-absent void
         even when total calls are few).
      3. adaptive: default for everything else, including quiet species
         (n_total < n_total_min) where neither pattern fires.

    The asymmetric n_total floor (no floor on strict, floor on frozen)
    handles low-call U12-absent species like ChlRei (n_total=2,
    only_a=100%, both=0) without over-firing the frozen route on noisy
    small samples.

    Args:
        stats:        Output of compute_asymmetry().
        only_f_thr:   frac_only_f threshold for frozen route. Default 0.10
                      (stable across 45-species panel; range [0.01, 0.19]
                      gives identical routing).
        only_a_thr:   frac_only_a threshold for strict route. Default 0.50
                      (range [≤0.30, ≥0.80] gives identical routing).
        n_both_max:   Max n_both allowed for strict route. Default 5.
        n_total_min:  Min n_total to apply frozen route. Default 5.

    Returns:
        "adaptive" | "frozen" | "strict".
    """
    if stats["frac_only_f"] >= only_f_thr and stats["n_total"] >= n_total_min:
        return "frozen"
    if stats["frac_only_a"] >= only_a_thr and stats["n_both"] <= n_both_max:
        return "strict"
    return "adaptive"


def apply_route(
    scores_adaptive: np.ndarray,
    scores_frozen: np.ndarray,
    route: Route,
    threshold: float = 90.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Produce final per-intron scores + call mask under the chosen route.

    Args:
        scores_adaptive: per-intron probability * 100 under adaptive z.
        scores_frozen:   per-intron probability * 100 under frozen z.
        route:           "adaptive" | "frozen" | "strict".
        threshold:       call threshold (default 90).

    Returns:
        (final_scores, called_mask):
          - final_scores: per-intron probability * 100 to populate svm_score.
            For adaptive/frozen, this is the chosen mode's score. For
            strict, this is min(adaptive, frozen) so the score reflects
            the agreement constraint.
          - called_mask: boolean array of which introns the rule calls.
    """
    if route == "adaptive":
        final = scores_adaptive
        called = scores_adaptive >= threshold
    elif route == "frozen":
        final = scores_frozen
        called = scores_frozen >= threshold
    elif route == "strict":
        final = np.minimum(scores_adaptive, scores_frozen)
        called = (scores_adaptive >= threshold) & (scores_frozen >= threshold)
    else:
        raise ValueError(f"Unknown route: {route!r}")
    return final, called
