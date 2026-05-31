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

# Continuous per-intron discount defaults (v2.7+).
#
# Empirically tuned against the 14-species panel + Salpingoeca:
#   - k_overcall=2.0, tau_overcall=0.0 catches Salpingoeca-class
#     overcalls (svm_vs_naive ≈ +24 → penalty ≈ 48) without affecting
#     panel TPs (which have svm_vs_naive < 0; no overcall penalty fires)
#   - k_weakmot=0.0 DISABLES the weak-motif penalty by default. Empirically
#     it loses 4 IPA-validated borderline TPs (HomSap SUDS3 + ARPC5,
#     DanRer ap4e1, XenTro mios) — confirmed-by-conservation U12s where
#     the SVM correctly leverages non-motif features (bp_offset,
#     bp_scan_confidence, support2) to recover real U12s with decayed
#     primary motif evidence. The weak-motif penalty conflates these with
#     Salpingoeca-class noise even though the SVM is using non-motif
#     features differently in the two regimes.
#
# Empirical panel result with these defaults (vs raw v2.6 control 1950/3):
#   TP=1950  FN=3  FP_strong=0  FP_any=216  Salpingoeca = 29 → 6
DEFAULT_DISCOUNT_K_OVERCALL = 2.0
DEFAULT_DISCOUNT_TAU_OVERCALL = 0.0
DEFAULT_DISCOUNT_K_WEAKMOT = 0.0
DEFAULT_DISCOUNT_TAU_MOTIF = 10.0

# Boundary-mass diagnostic: fraction of eligible introns where second-pass
# mean P sits in [0.1, 0.9]. Subsample size for cost-efficient estimation.
# At BM=0.05, standard error ≈ 0.003 — sufficient precision for diagnostic.
DEFAULT_BM_SUBSAMPLE_SIZE = 5000
DEFAULT_BM_LOWER = 0.1
DEFAULT_BM_UPPER = 0.9


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
    gap_fraction: float | None = None          # gap_width/Δmean on the Fisher axis (gate decision stat)
    gap_fraction_ucl: float | None = None       # deterministic-bootstrap UCL (π_species input)


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
    valley_min: float = DEFAULT_VALLEY_MIN,        # legacy KDE-valley floor; no longer the decision (kept for diag)
    gap_fraction_min: float = 0.0,                 # lenient route floor: gap_fraction>0 ⇔ u12_p25 > u2_p99
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

    # Check 4 — separation gap. gap_fraction (gap_width/Δmean on the Fisher axis) REPLACES the
    # unstable KDE median_depth, which spuriously failed real low-N bearers (snRNA-confirmed: the
    # whole no_kde_valley bucket was real). This route check is deliberately LENIENT — route to
    # modesep whenever a clear gap exists (gap_fraction > gap_fraction_min). The GRADED suppression
    # of weak/uncertain separation happens downstream in π_species (keyed on the gap_fraction
    # bootstrap UCL). median_depth is still computed + surfaced as a diagnostic only. None
    # gap_fraction (compute_valley_depth early-returns: insufficient/coincident/overlap) ⇒ no gap ⇒ fail.
    valley_result = valley_depth_fn(u12_candidate_points, u2_candidate_points)
    gap_fraction = valley_result.get("gap_fraction")
    gap_fraction_ucl = valley_result.get("gap_fraction_ucl")
    valley_depth = valley_result.get("median_depth")
    _vd = float(valley_depth) if valley_depth is not None and np.isfinite(valley_depth) else None
    _gf = float(gap_fraction) if gap_fraction is not None and np.isfinite(gap_fraction) else None
    # Keep -inf (a valid "confident no-separation" UCL for π_species); drop only None/NaN.
    _gfucl = float(gap_fraction_ucl) if gap_fraction_ucl is not None and not np.isnan(gap_fraction_ucl) else None
    if _gf is None or _gf <= gap_fraction_min:
        return GateDecision(
            passes=False,
            reason="below_gap_fraction",
            valley_depth=_vd,
            mu_u12_offset=mu_offset,
            n_eff_candidates=float(n_eff_candidates),
            gap_fraction=_gf,
            gap_fraction_ucl=_gfucl,
        )

    return GateDecision(
        passes=True,
        reason="ok",
        valley_depth=_vd,
        mu_u12_offset=mu_offset,
        n_eff_candidates=float(n_eff_candidates),
        gap_fraction=_gf,
        gap_fraction_ucl=_gfucl,
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


# -----------------------------------------------------------------------------
# Continuous per-intron discount (v2.7+)
# -----------------------------------------------------------------------------

def continuous_discount_logit_shift(
    logit_svm: np.ndarray,
    raw_sum: np.ndarray,
    k_overcall: float = DEFAULT_DISCOUNT_K_OVERCALL,
    tau_overcall: float = DEFAULT_DISCOUNT_TAU_OVERCALL,
    k_weakmot: float = DEFAULT_DISCOUNT_K_WEAKMOT,
    tau_motif: float = DEFAULT_DISCOUNT_TAU_MOTIF,
) -> np.ndarray:
    """Compute the per-intron continuous discount in log-odds space.

    Two non-negative penalty terms:
      penalty_overcall = k_overcall * max(0, svm_vs_naive - tau_overcall)
      penalty_weakmot  = k_weakmot  * max(0, tau_motif - raw_sum)

    where svm_vs_naive = logit(p_svm) - raw_sum. Both terms are zero in the
    "healthy" regime (svm tracks naive log-LR + motif evidence is strong);
    they activate only when the SVM is overcalling relative to motif log-LR
    sum OR when motif evidence is weak.

    Returns the (non-positive) logit shift to apply: logit_final = logit_svm - shift.
    """
    logit_svm = np.asarray(logit_svm, dtype=float)
    raw_sum = np.asarray(raw_sum, dtype=float)
    svm_vs_naive = logit_svm - raw_sum
    penalty_oc = k_overcall * np.maximum(0.0, svm_vs_naive - tau_overcall)
    penalty_wm = k_weakmot * np.maximum(0.0, tau_motif - raw_sum)
    return penalty_oc + penalty_wm


def apply_continuous_discount(
    svm_scores: np.ndarray,
    raw_sums: np.ndarray,
    **kwargs,
) -> np.ndarray:
    """Apply continuous discount to svm_scores (0-100) given raw motif log-LR sums.

    Returns adjusted scores (0-100). Per-intron shift is non-positive, so
    adjusted_score <= svm_score always.

    kwargs are forwarded to continuous_discount_logit_shift().
    """
    svm = np.asarray(svm_scores, dtype=float)
    raw = np.asarray(raw_sums, dtype=float)
    eps = 1e-9
    p = np.clip(svm / 100.0, eps, 1 - eps)
    logit_p = np.log(p / (1 - p))
    shift = continuous_discount_logit_shift(logit_p, raw, **kwargs)
    logit_adj = logit_p - shift
    p_adj = 1.0 / (1.0 + np.exp(-logit_adj))
    return p_adj * 100.0


# -----------------------------------------------------------------------------
# Boundary mass (species-level OOD diagnostic, v2.7+)
# -----------------------------------------------------------------------------

def compute_boundary_mass(
    mean_p: np.ndarray,
    lower: float = DEFAULT_BM_LOWER,
    upper: float = DEFAULT_BM_UPPER,
) -> float:
    """Fraction of introns whose ensemble mean P sits in [lower, upper].

    High boundary mass signals an OOD species — the SVM hovers diffusely
    rather than committing cleanly to U2 vs U12 calls. Diagnostic only;
    not used in any gate decision in v2.7.
    """
    mp = np.asarray(mean_p, dtype=float)
    if len(mp) == 0:
        return float("nan")
    return float(((mp >= lower) & (mp <= upper)).mean())


def voting_fraction(per_model_probs: np.ndarray, threshold: float = 0.5) -> np.ndarray:
    """Per-intron fraction of ensemble sub-models voting U12 (P > threshold).

    Input shape: (n_models, n_introns). Returns 1-D array of length n_introns.
    Diagnostic only — surfaced alongside ensemble_sigma for users to inspect
    borderline cases where mean P is moderate but per-model agreement is high.
    """
    probs = np.asarray(per_model_probs, dtype=float)
    if probs.ndim != 2:
        raise ValueError(
            f"per_model_probs must be 2-D (n_models × n_introns); got "
            f"shape {probs.shape}"
        )
    return (probs > threshold).mean(axis=0)


__all__ = [
    "DEFAULT_THRESHOLD",
    "DEFAULT_N_FLOOR",
    "DEFAULT_VALLEY_MIN",
    "DEFAULT_Z5P_ELIGIBILITY",
    "DEFAULT_CANDIDATE_CENTER",
    "DEFAULT_CANDIDATE_STEEPNESS",
    "DEFAULT_UNIVERSAL_ANCHORS",
    "DEFAULT_MU_U12_TOLERANCE_5P",
    "DEFAULT_DISCOUNT_K_OVERCALL",
    "DEFAULT_DISCOUNT_TAU_OVERCALL",
    "DEFAULT_DISCOUNT_K_WEAKMOT",
    "DEFAULT_DISCOUNT_TAU_MOTIF",
    "DEFAULT_BM_SUBSAMPLE_SIZE",
    "DEFAULT_BM_LOWER",
    "DEFAULT_BM_UPPER",
    "ModeSeparationStats",
    "GateDecision",
    "candidate_weight_from_svm",
    "fit_mode_separation",
    "apply_mode_separation_z",
    "compute_support2",
    "evaluate_gate",
    "continuous_discount_logit_shift",
    "apply_continuous_discount",
    "compute_boundary_mass",
    "voting_fraction",
]
