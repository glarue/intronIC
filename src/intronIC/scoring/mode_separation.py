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
# DEPRECATED as a DECISION input (v3 gate-gapfrac Step 4): the absolute-anchor location check
# over-killed 24 snRNA-confirmed real U12 bearers (μ_U12 weighted-median is biased low at small
# n_eff, and divergent motifs score lower vs the reference PWM — both shift the absolute location
# but preserve species-internal separation). Superseded by the centroid_sigma floor below. The
# anchor (15.671) is retained only to compute the diagnostic mu_u12_offset.
DEFAULT_MU_U12_TOLERANCE_5P = 3.6

# Centroid-separation floor (U12↔U2 centroid distance in U2-σ units on the Fisher axis, i.e.
# `centroid_sigma` from compute_valley_depth) — the PHYLOGENY-NAIVE, species-internal replacement
# for the absolute μ_U12 anchor check. Guards against z-anchor corruption (a weak first-pass
# candidate cluster being mapped to z=1) using only the species' own intron pool. STEP-4 provisional
# value (gate-pass species observed ≥ 4.30; confirmed-absent cluster 3.3–4.6) — Step 6 calibrates it
# jointly with the π_species operating point on the snRNA + conservation anchors.
DEFAULT_CSIG_FLOOR = 4.0

# Core-presence bar (Step 7): the raw_sum (5'_raw + bp_raw + 3'_raw) above which a first-pass
# candidate counts as a genuine strong-motif U12 (a "core" intron). The species-level
# `core_fraction` = frac(candidate raw_sum > this bar) is the LABEL-FREE signal that distinguishes
# U12-LOSS species (core ~0-1% — their "U12 cluster" is chance weak-motif hits) from real bearers,
# however diverged (core 41-82%). Unlike the removed μ_U12 LOCATION anchor, a motif-strength bar is
# legitimately external: U12 identity is a UNIVERSAL biological constant (conserved U11/U12 snRNA
# binding), so a diverged real U12 still scores high motif log-odds. DERIVED, not arbitrary: ≈ the
# 5th percentile of the multispecies gold-U12 raw_sum distribution (conservation-labeled, IPA-
# independent of raw_sum). VALIDATED leave-clade-out: every held-out clade keeps 86-100% core under a
# bar derived from the OTHER clades (so it generalizes beyond its calibration lineages → not circular
# with PWM training); loss species stay at 1%. The wide 40-pt margin makes the exact value non-
# critical (anything in ~5-20 classifies identically). Re-derive per PWM retrain; travels w/ bundle.
#
# VALUE: set to 8.0 — deliberately the LOW end of the safe valley, NOT the gold-U12 p5 (≈13.9). The
# 12-diverged-bearer pressure test (Zoopago/Glomero/oomycete/amoebozoa/Mucorales) found the WIDEST
# loss-vs-real gap at bar≈8 (+35%): lowering 13.9→8 nearly doubles the worst diverged bearer's core
# margin (MucorLusit 22%→44%) at trivial loss-side cost (loss max 1%→9%). Choosing the low end is the
# explicit MITIGATION for the human-PWM anchoring of raw_sum — it maximizes diverged-real safety while
# every value in [5, 13.9] classifies loss-vs-real identically (loss core ≤9%, real ≥44%). See
# conservation_corpus/docs/calibration_plan.md §6i-6j.
DEFAULT_CORE_RAW_SUM_BAR = 8.0

# UNIFIED species-trust (experimental): if the first-pass candidate cluster has essentially no
# strong-motif core (core_fraction < this floor), route to FALLBACK instead of modesep — so a
# motif-empty species never enters the 2nd-pass z-anchor balloon (vs the current design, which
# balloons then claws back via the discount). 0.20 sits in the valley between loss species (≤0.09)
# and diverged real bearers (≥0.41) — a 30-pt margin. Set <0 to disable (= current behaviour).
DEFAULT_CORE_GATE_FLOOR = 0.20

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
# Species-level BPS penalty (frac_bp6 + logN), v6 (C6). Ported verbatim from the
# validated unified-species-trust branch. A gated, ONE-SIDED (<=0) additive logit
# shift applied to a species' calls AFTER the per-intron discount/graduated tail:
# a species whose post-scoring HC set looks non-bearer-like (low frac of U12-BPS,
# small N) gets demoted; confident bearers (p_bearer >= p_gate) are untouched.
# DISABLED by default (enable_species_penalty=False) — the shipped default bundle
# is byte-unchanged; the graduated bundle opts in only after the C6/V eval.
# Constants = full-precision 185-sp recal fit (doc-rounded values in comments).
# -----------------------------------------------------------------------------
DEFAULT_SPECIES_PENALTY_PGATE = 0.85
SPECIES_PENALTY_PI_TRAIN = 0.5
SPECIES_PENALTY_FLOOR = 0.001
SPECIES_PENALTY_INTERCEPT = -5.8433825296    # doc-rounded: -5.84338
SPECIES_PENALTY_COEF_FRACBP6 = 1.4928134895  # doc-rounded:  1.49281
SPECIES_PENALTY_COEF_LOGN = 1.9895442351     # doc-rounded:  1.98954


def species_penalty_logit_shift(
    frac_bp6: float,
    n_hc: int,
    *,
    p_gate: float = DEFAULT_SPECIES_PENALTY_PGATE,
    pi_train: float = SPECIES_PENALTY_PI_TRAIN,
    floor: float = SPECIES_PENALTY_FLOOR,
) -> float:
    """Gated species penalty → additive (<=0) logit shift to apply to adjusted_score.

    frac_bp6 = fraction of THIS species' HC (adjusted>=threshold) calls with bp_raw>=6 (U12-BPS present);
    n_hc = HC count. p_bearer = sigmoid(intercept + c_fb*frac_bp6 + c_logN*log1p(n_hc)). Confident
    bearers (p>=p_gate) → 0 shift (untouched); otherwise pi = floor + (pi_train-floor)*(p/p_gate) and the
    shift is log(pi/pi_train) < 0 (one-sided: never promotes). Returns 0.0 when n_hc==0 (nothing to penalize).
    Applied post-gate, so gate-catchable losses arrive with low HC, sidestepping the metric's high-N blind spot.
    """
    if n_hc <= 0:
        return 0.0
    z = (SPECIES_PENALTY_INTERCEPT
         + SPECIES_PENALTY_COEF_FRACBP6 * float(frac_bp6)
         + SPECIES_PENALTY_COEF_LOGN * float(np.log1p(n_hc)))
    p = 1.0 / (1.0 + float(np.exp(-z)))
    pi = pi_train if p >= p_gate else floor + (pi_train - floor) * (p / p_gate)
    return float(np.log(pi / pi_train))


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
    centroid_sigma: float | None = None         # U12↔U2 centroid distance in U2-σ (check #3 replacement)
    core_fraction: float | None = None          # frac(candidate raw_sum > bar): strong-motif U12 core presence


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
    mu_u12_tolerance: float = DEFAULT_MU_U12_TOLERANCE_5P,  # DEPRECATED: vestigial, no longer a check
    csig_floor: float = DEFAULT_CSIG_FLOOR,
    core_raw_sum_bar: float = DEFAULT_CORE_RAW_SUM_BAR,
    core_gate_floor: float = DEFAULT_CORE_GATE_FLOOR,
) -> GateDecision:
    """Decide whether to apply mode-separation for the species.

    Checks (cheap → expensive); all PHYLOGENY-NAIVE (species-internal — no taxonomy input):

    1. **n_eff floor**: too few first-pass candidates → mode estimates are
       noisy. Falls back to first-pass scores.
    2. **degenerate separation**: μ_U12 ≤ μ_U2 (inverted) → fall back.
    3/4. **Separation (Fisher discriminant, via `compute_valley_depth`)** — two species-internal
       conditions on the candidate clusters, both read from the single injected `valley_depth_fn`
       call (= `cluster_validation.compute_valley_depth`, injected to keep this module import-light):
         - **gap_fraction > `gap_fraction_min`** (lenient route floor): a positive gap between the
           U12 and U2 clusters on the Fisher axis. Graded suppression of weak/uncertain separation
           is deferred downstream to π_species (keyed on the gap_fraction bootstrap UCL).
         - **centroid_sigma ≥ `csig_floor`**: the U12↔U2 centroid distance in U2-σ units. This
           REPLACES the former absolute-anchor μ_U12 location check (v3 Step 4) — it guards the same
           failure mode (a weak first-pass candidate cluster being mapped to z=1 by the mode-sep
           z-transform → spurious 2nd-pass U12 calls) but species-internally, so it does not
           over-kill diverged real bearers whose absolute μ_U12 is shifted by estimator bias or
           motif divergence. `median_depth` (KDE valley) and `mu_u12_offset` (vs the universal
           anchor) are still computed + surfaced as DIAGNOSTICS only.

    Note: the separation metrics operate jointly in 3D over all motif features via Fisher's linear
    discriminant — NOT a cheap per-feature moment proxy.
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

    # mu_u12_offset vs the universal anchor — DIAGNOSTIC only now (v3 Step 4 removed the absolute
    # location FAIL; it over-killed diverged real bearers — see DEFAULT_CSIG_FLOOR). mu_u12_tolerance
    # is retained in the signature for back-compat but no longer gates anything.
    mu_offset = float(mu_u12_5p_raw - mu_u12_prior)

    # core_fraction (Step 7) — the LABEL-FREE species-level signal: fraction of first-pass candidates
    # whose summed motif log-odds (raw_sum = 5'_raw+bp_raw+3'_raw, = col-sum of the raw 3D candidate
    # points) clears the gold-derived motif bar. Loss species ≈ 0-1% (no genuine U12 core); real
    # bearers (however diverged) 41-82%. DIAGNOSTIC here; drives the graded π_species core-presence
    # term downstream (see cluster_validation._confidence_shrunk_pi_species). Computed once and carried
    # on every post-floor GateDecision.
    _core_fraction = None
    if u12_candidate_points is not None and len(u12_candidate_points) >= 1:
        _rs = u12_candidate_points.sum(axis=1)
        _core_fraction = float(np.mean(_rs > core_raw_sum_bar))

    # Checks 3/4 — separation on the Fisher discriminant, both from ONE compute_valley_depth call:
    #   (#4) gap_fraction > gap_fraction_min : lenient route floor (a positive U12/U2 gap); graded
    #        suppression of weak separation is deferred to π_species (gap_fraction bootstrap UCL).
    #   (#3) centroid_sigma ≥ csig_floor     : the species-internal replacement for the absolute
    #        μ_U12 anchor — guards z-anchor corruption without over-killing diverged real bearers.
    # median_depth is computed + surfaced as a diagnostic only. None gap_fraction (early-returns:
    # insufficient/coincident/overlap) ⇒ no gap ⇒ fail.
    valley_result = valley_depth_fn(u12_candidate_points, u2_candidate_points)
    gap_fraction = valley_result.get("gap_fraction")
    gap_fraction_ucl = valley_result.get("gap_fraction_ucl")
    centroid_sigma = valley_result.get("centroid_sigma")
    valley_depth = valley_result.get("median_depth")
    _vd = float(valley_depth) if valley_depth is not None and np.isfinite(valley_depth) else None
    _gf = float(gap_fraction) if gap_fraction is not None and np.isfinite(gap_fraction) else None
    # Keep -inf (a valid "confident no-separation" UCL for π_species); drop only None/NaN.
    _gfucl = float(gap_fraction_ucl) if gap_fraction_ucl is not None and not np.isnan(gap_fraction_ucl) else None
    _csig = float(centroid_sigma) if centroid_sigma is not None and np.isfinite(centroid_sigma) else None
    if _gf is None or _gf <= gap_fraction_min:
        return GateDecision(
            passes=False,
            reason="below_gap_fraction",
            valley_depth=_vd,
            mu_u12_offset=mu_offset,
            n_eff_candidates=float(n_eff_candidates),
            gap_fraction=_gf,
            gap_fraction_ucl=_gfucl,
            centroid_sigma=_csig,
            core_fraction=_core_fraction,
        )
    if _csig is None or _csig < csig_floor:
        return GateDecision(
            passes=False,
            reason="low_centroid_sigma",
            valley_depth=_vd,
            mu_u12_offset=mu_offset,
            n_eff_candidates=float(n_eff_candidates),
            gap_fraction=_gf,
            gap_fraction_ucl=_gfucl,
            centroid_sigma=_csig,
            core_fraction=_core_fraction,
        )

    # UNIFIED species-trust (experimental): a motif-empty cluster (no strong-motif U12 core) that
    # passes the geometric checks above is the z-anchor-corruption case (TetThe/Chlamydomonas) — it
    # would balloon in the 2nd pass then need clawing back. Route it to FALLBACK instead, so it never
    # balloons. Disabled when core_gate_floor < 0.
    if core_gate_floor >= 0 and _core_fraction is not None and _core_fraction < core_gate_floor:
        return GateDecision(
            passes=False,
            reason="low_core_fraction",
            valley_depth=_vd,
            mu_u12_offset=mu_offset,
            n_eff_candidates=float(n_eff_candidates),
            gap_fraction=_gf,
            gap_fraction_ucl=_gfucl,
            centroid_sigma=_csig,
            core_fraction=_core_fraction,
        )

    # RECALL TEST (env-gated, experimental): force this gate-passing species to fallback so we can
    # measure what the 2nd-pass mode-sep actually buys (modesep HC vs fallback HC = the 2nd-pass
    # recall contribution). All stats are computed above, so the fallback π_species still gets them.
    import os as _os
    if _os.environ.get("INTRONIC_FORCE_FALLBACK"):
        return GateDecision(
            passes=False,
            reason="forced_fallback",
            valley_depth=_vd,
            mu_u12_offset=mu_offset,
            n_eff_candidates=float(n_eff_candidates),
            gap_fraction=_gf,
            gap_fraction_ucl=_gfucl,
            centroid_sigma=_csig,
            core_fraction=_core_fraction,
        )

    return GateDecision(
        passes=True,
        reason="ok",
        valley_depth=_vd,
        mu_u12_offset=mu_offset,
        n_eff_candidates=float(n_eff_candidates),
        gap_fraction=_gf,
        gap_fraction_ucl=_gfucl,
        centroid_sigma=_csig,
        core_fraction=_core_fraction,
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
