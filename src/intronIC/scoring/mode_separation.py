"""Mode-separation z-scoring (intronIC v2.6).

Instead of standard `z = (raw - μ) / σ` (where σ captures only U2 spread),
mode-separation uses

    z = (raw - μ_U2_species) / (μ_U12_species - μ_U2_species)

so that z=0 ↔ U2 mode, z=1 ↔ U12 mode in every species. The classifier
sees a fixed target geometry; species-specific U12/U2 separation
compression (e.g., plants) is absorbed by the normalizer.

This module is pure transform logic. Two-pass classify orchestration,
gate decisions using compute_valley_depth(), and bundle loading live
elsewhere in the pipeline.

References
----------
- Empirical universal U12 anchors: 77-species sqrt-n weighted means
- Location-prior tolerance: 3·stdev + 0.5σ buffer = 3.6 raw PWM units
- Eligibility floor z_5p ≥ 0.3: legacy intuition, validated by Phase 0
  derivation (worst observed U12 lies at z_5p = 0.673; 0.37 margin)
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


# -----------------------------------------------------------------------------
# Defaults (frozen by Phase 0 derivation; see MODESEP_INTEGRATION_PLAN.md)
# -----------------------------------------------------------------------------

DEFAULT_THRESHOLD = 90.0
DEFAULT_N_FLOOR = 5
DEFAULT_VALLEY_MIN = 0.30
DEFAULT_Z5P_ELIGIBILITY = 0.30

DEFAULT_CANDIDATE_CENTER = 90.0
DEFAULT_CANDIDATE_STEEPNESS = 5.0

# Universal U12 anchors (raw PWM units), sqrt-n weighted across 77 species.
DEFAULT_UNIVERSAL_ANCHORS = {
    "five_raw": 15.671,
    "bp_raw": 6.704,
    "three_raw": 1.764,
}

# Location-prior tolerance on inferred μ_U12 5'_raw (deviation from anchor).
DEFAULT_MU_U12_TOLERANCE_5P = 3.6


# -----------------------------------------------------------------------------
# Data classes
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class ModeSeparationStats:
    """Per-feature U2/U12 mode estimates for one species."""
    mu_u2: float
    mu_u12: float
    spread_u2: float
    separation: float
    n_eff_candidates: float


@dataclass(frozen=True)
class GateDecision:
    """Outcome of the mode-separation gate."""
    passes: bool
    reason: str
    valley_depth: float | None = None
    mu_u12_offset: float | None = None
    n_eff_candidates: float | None = None


# -----------------------------------------------------------------------------
# Candidate weighting
# -----------------------------------------------------------------------------

def candidate_weight_from_svm(
    svm_scores: np.ndarray,
    center: float = DEFAULT_CANDIDATE_CENTER,
    steepness: float = DEFAULT_CANDIDATE_STEEPNESS,
) -> np.ndarray:
    """Logistic-gated soft candidate weight from first-pass SVM scores.

    With defaults: w ≈ 0 below svm=80, 0.5 at svm=90, ≈ 1.0 above svm=100.
    Smooth gating (not hard binarization) means weighted-median estimates
    of (μ_U2, μ_U12) degrade gracefully if first-pass calls are uncertain.
    """
    z = (np.asarray(svm_scores, dtype=float) - center) / steepness
    return 1.0 / (1.0 + np.exp(-z))


# -----------------------------------------------------------------------------
# Mode-separation z-transform
# -----------------------------------------------------------------------------

def fit_mode_separation(
    raw_scores: np.ndarray,
    candidate_weight: np.ndarray,
) -> ModeSeparationStats:
    """Estimate (μ_U2, μ_U12) for one feature column.

    Weighted-median is used for both modes — robust to first-pass
    misclassification at the tails. μ_U2 uses (1 - w) weights;
    μ_U12 uses w. Spread is the (1-w)-weighted IQR (kept for diagnostics).
    """
    raw = np.asarray(raw_scores, dtype=float)
    w = np.asarray(candidate_weight, dtype=float)
    if raw.shape != w.shape:
        raise ValueError(
            f"raw_scores and candidate_weight must have the same shape; "
            f"got {raw.shape} vs {w.shape}"
        )

    inv_w = 1.0 - w
    mu_u2 = _weighted_quantile(raw, inv_w, 0.5)
    q75 = _weighted_quantile(raw, inv_w, 0.75)
    q25 = _weighted_quantile(raw, inv_w, 0.25)
    mu_u12 = _weighted_quantile(raw, w, 0.5)

    return ModeSeparationStats(
        mu_u2=float(mu_u2),
        mu_u12=float(mu_u12),
        spread_u2=float(q75 - q25),
        separation=float(mu_u12 - mu_u2),
        n_eff_candidates=float(w.sum()),
    )


def apply_mode_separation_z(
    raw_scores: np.ndarray,
    stats: ModeSeparationStats,
) -> np.ndarray:
    """Apply z = (raw - μ_U2) / (μ_U12 - μ_U2). Caller must ensure separation > 0."""
    if stats.separation <= 0:
        raise ValueError(
            f"degenerate mode separation: μ_U12 ({stats.mu_u12}) ≤ μ_U2 "
            f"({stats.mu_u2}); the caller should have gated out this species"
        )
    return (np.asarray(raw_scores, dtype=float) - stats.mu_u2) / stats.separation


# -----------------------------------------------------------------------------
# Consensus feature
# -----------------------------------------------------------------------------

def compute_support2(z_5: np.ndarray, z_bp: np.ndarray, z_3: np.ndarray) -> np.ndarray:
    """Second-largest of clipped-positive z's per intron.

    Acts as a "two-of-three motifs strong" consensus feature — high
    when at least two of the three splice signals positively support U12.
    """
    z = np.stack([z_5, z_bp, z_3], axis=1)
    z_clip = np.clip(z, 0.0, None)
    z_sorted = np.sort(z_clip, axis=1)
    return z_sorted[:, 1]


# -----------------------------------------------------------------------------
# Gate decision
# -----------------------------------------------------------------------------

def evaluate_gate(
    u12_candidate_points: np.ndarray,
    u2_candidate_points: np.ndarray,
    mu_u12_5p_raw: float,
    n_eff_candidates: float,
    separation_5p: float,
    *,
    valley_depth_fn,
    n_floor: int = DEFAULT_N_FLOOR,
    valley_min: float = DEFAULT_VALLEY_MIN,
    mu_u12_prior: float = DEFAULT_UNIVERSAL_ANCHORS["five_raw"],
    mu_u12_tolerance: float = DEFAULT_MU_U12_TOLERANCE_5P,
) -> GateDecision:
    """Decide whether to apply mode-separation for the species.

    Three independent checks (cheap → expensive):

    1. **n_eff floor**: too few first-pass candidates → mode estimates are
       noisy. Falls back to first-pass scores.
    2. **μ_U12 location prior**: catches the failure mode where a noisy
       first-pass classifier confidently mis-locates μ_U12 in a U12-bearing-
       LIKE pattern. The inferred μ_U12_5'_raw must sit within
       ±`mu_u12_tolerance` of the cross-species universal anchor.
       Empirically (61 species, 12 phyla): max observed offset is 2.90
       raw PWM units; shipping tolerance 3.6 includes a 0.6-unit buffer
       for held-out species.
    3. **Density valley**: uses the multi-bandwidth Fisher-discriminant
       KDE valley detection from `cluster_validation.compute_valley_depth`
       (injected as `valley_depth_fn` so this module remains import-light
       and testable in isolation). Median valley depth ≥ `valley_min`
       implies a real bimodal U12/U2 separation in feature space.

    Note: the valley detection operates jointly in 3D over all motif features
    via Fisher's linear discriminant — NOT a cheap per-feature moment proxy.

    `valley_depth_fn` should be `cluster_validation.compute_valley_depth`
    (or a wrapper); we inject it instead of importing it directly to keep
    this module pure and avoid a scipy import in test contexts.
    """
    if n_eff_candidates < n_floor:
        return GateDecision(
            passes=False,
            reason="below_n_floor",
            n_eff_candidates=float(n_eff_candidates),
        )

    if separation_5p <= 0:
        return GateDecision(
            passes=False,
            reason="degenerate_separation",
            n_eff_candidates=float(n_eff_candidates),
        )

    mu_offset = float(mu_u12_5p_raw - mu_u12_prior)
    if abs(mu_offset) > mu_u12_tolerance:
        return GateDecision(
            passes=False,
            reason="u12_mode_outside_prior_range",
            mu_u12_offset=mu_offset,
            n_eff_candidates=float(n_eff_candidates),
        )

    valley_result = valley_depth_fn(u12_candidate_points, u2_candidate_points)
    valley_depth = valley_result.get("median_depth")
    if valley_depth is None or not np.isfinite(valley_depth) or valley_depth < valley_min:
        return GateDecision(
            passes=False,
            reason="no_kde_valley",
            valley_depth=(float(valley_depth)
                          if valley_depth is not None and np.isfinite(valley_depth)
                          else None),
            mu_u12_offset=mu_offset,
            n_eff_candidates=float(n_eff_candidates),
        )

    return GateDecision(
        passes=True,
        reason="ok",
        valley_depth=float(valley_depth),
        mu_u12_offset=mu_offset,
        n_eff_candidates=float(n_eff_candidates),
    )


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def _weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    """Weighted quantile via cumulative weight interpolation.

    Zero-weight input falls back to the plain quantile.
    """
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    total = weights.sum()
    if total <= 0:
        return float(np.quantile(values, q))
    order = np.argsort(values)
    cumw = np.cumsum(weights[order])
    cutoff = q * total
    idx = np.searchsorted(cumw, cutoff, side="left")
    idx = max(0, min(idx, len(values) - 1))
    return float(values[order][idx])


__all__ = [
    "DEFAULT_THRESHOLD",
    "DEFAULT_N_FLOOR",
    "DEFAULT_VALLEY_MIN",
    "DEFAULT_Z5P_ELIGIBILITY",
    "DEFAULT_CANDIDATE_CENTER",
    "DEFAULT_CANDIDATE_STEEPNESS",
    "DEFAULT_UNIVERSAL_ANCHORS",
    "DEFAULT_MU_U12_TOLERANCE_5P",
    "ModeSeparationStats",
    "GateDecision",
    "candidate_weight_from_svm",
    "fit_mode_separation",
    "apply_mode_separation_z",
    "compute_support2",
    "evaluate_gate",
]
