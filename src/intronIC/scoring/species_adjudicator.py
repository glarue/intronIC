"""Species-level U12 adjudicator: decide whether a genome's motif-strong introns form a real, separated
U12 *population* (bearer) or are the upper tail of the unimodal U2 background (loss / recent-loss), and
turn that into a per-genome bearer probability ``q`` with a bootstrap CI.

Pipeline role (see ``docs/raw_gated_scoring.md`` §0b/§0c): the SVM gives a per-intron, species-agnostic
motif probability ``P_motif`` (calibrated ensemble margin). This module answers the per-*species* question
and produces the per-intron posterior ``P_adj = q * P_motif`` (functional-U12 posterior; clean because
``P(functional U12 | motif, not-a-bearer) ≈ 0``). The low-numbers uncertainty flows through ``q``'s CI.

Method (PRIMARY = z_excess population statistic + an EMPIRICAL gap-anchored gate). See
``docs/adjudicator_qdriver_postmortem.md`` for the full design + why the prior ``depth_tail -> q ->
adjusted_score`` chain was superseded (it baked in the 50% threshold and reported the high-confidence
count off an unidentified logistic slope; ``depth_tail`` is also count-blind and mis-ranks numerous-but-
shallow divergent bearers).
  - U2 reference = introns with ``P_motif < u2_threshold``; calls = ``P_motif >= call_threshold``.
  - **PRIMARY = ``z_excess``** = Poisson significance of the strong-call COUNT above what the genome's OWN
    U2 tail predicts: "is there a real U12 *population* beyond chance?" A count (so it sees numerous-but-
    shallow divergent bearers) referenced to the per-genome background (so it is robust at both size
    extremes, unlike raw count / density). Bearer/loss separate cleanly on the snRNA-corroborated panel.
  - **EMPIRICAL gap gate** (no fitted probability slope -> no threshold-baking): ``z_excess >=
    BEARER_FLOOR_Z`` -> DETECTED (a U12-motif population); ``<= LOSS_CEILING_Z`` -> NOT_DETECTED; in between
    (the empirical gap, the out-of-support region) -> INCONCLUSIVE (abstain). Anchors are frozen class extremes from
    the labelled panel, not tuning targets.
  - ``depth_tail``/``q`` and the EVT ``excess_z``/``p_gumbel`` are retained as labelled SECONDARY/back-compat
    diagnostics only.

CAVEAT (important): a high ``z_excess`` is *confirmatory* of a bearer only because it empirically correlates
with snRNA presence (0 of 96 panel bearers were high-z with searched-and-absent snRNAs). It is NOT a logical
guarantee. The runtime gate is MOTIF-ONLY, so DETECTED means "a U12-motif population BY MOTIF, corroborate
downstream" -- a high-z genome whose snRNAs are *searched and absent* must be flagged for investigation, not
auto-accepted. That snRNA cross-check lives in CALIBRATION (where it flags SUSPICIOUS) and in the database
layer downstream -- never in this module.

Calibration provenance: ``loss_ceiling_z=2.60`` = aspergillus coremiiformis (max snRNA-confirmed loss), frozen
from the cross-clade panel + 68 divergent bearers (Infernal E<=0.01, >=3-of-4, self-consistently defining-aware).
``bearer_floor_z=5.50`` is a TRUST threshold (2026-06-30c), widened from the original 4.00 (Mycotypha, lowest
corroborated bearer) after the tier1 v3 run surfaced a near-floor UNCERTAIN cluster the gate can't separate:
blyttiomyces (z=4.54, probable bearer) and monocercomonoides (z=4.24, probable loss) sit ADJACENT in z with
opposite truth and are unresolvable even by augmented clade-CMs + protein machinery. So z_excess is a trust
line, not a separator, in that band. 5.50 sits in the empty [4.54, 6.95] gap (below the lowest CONFIRMED bearer
batrachochytrium 6.95) -> near-floor calls -> INCONCLUSIVE/corroborate. Version-pinned by ``ADJUDICATOR_PARAMS_VERSION``;
see WtMTA ``snrna_cm/docs/zone_resolution_plan.md``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

import numpy as np

#: Version pin for the calibration constants (bundle-stamped in production). Bump when the anchors or the
#: z_excess/EVT recipe changes so a stale bundle can be detected. z_excess gap gate stays PRIMARY; 2026-07-03
#: the additive call-strength gate switched from the CI lower bound (cs_p95_lo >= 4.50) to the POINT cs_p95
#: (>= cs_point_threshold 5.0). 2026-07-03b: the strength gate's PRIMARY driver became the per-genome Gumbel
#: significance of the call upper tail (p_gumbel_p95 <= p_gumbel_threshold), with the absolute POINT cs_p95
#: retained as a null-free co-fallback; see the gate comment + docs/adjudicator_strength_evt.md.
ADJUDICATOR_PARAMS_VERSION = "zexcess_gap_pgumbel_cs_2026-07-03"

# Calibration constants frozen from the cross-clade panel + 68 divergent bearers + snRNA-vetted losses
# (eval_corpus/{corroborate_v2,refit_anchors_expanded,robust_sep_one_nu2}.py; see
# docs/adjudicator_qdriver_postmortem.md §4a). Bundle-stamped in production; defaults document the values.
DEFAULT_PLATT_A = 2.796       #: P_motif = sigma(PLATT_A * margin + PLATT_C)
DEFAULT_PLATT_C = -1.178
# z_excess (the PRIMARY population statistic): Poisson significance of the strong-call COUNT vs the genome's
# OWN U2-tail prediction. Empirical bearer/loss gap on the snRNA-corroborated panel; LCO 0 cross-errors.
DEFAULT_LOSS_CEILING_Z = 2.60   #: z_excess <= this => NOT_DETECTED (max snRNA-confirmed-loss; aspergillus)
DEFAULT_BEARER_FLOOR_Z = 5.50   #: z_excess >= this => DETECTED. TRUST threshold (not a separator): above the
                                #: near-floor UNCERTAIN cluster (~4.2-4.5; blyttiomyces/monocercomonoides — probable
                                #: bearer & probable loss ADJACENT in z, unresolvable even by augmented snRNA/protein),
                                #: below the lowest CONFIRMED bearer (batrachochytrium 6.95). Robust: any value in the
                                #: empty [4.54, 6.95] gap gives the identical split. Widened 2026-06-30c from 4.00; see
                                #: WtMTA snrna_cm/docs/zone_resolution_plan.md. Below this -> INCONCLUSIVE (corroborate).
# depth_tail/q are RETAINED as a labelled SECONDARY/back-compat diagnostic only (NOT the driver). The
# depth_tail->q logistic baked in the 50% threshold and reported HC off an unidentified slope; superseded.
DEFAULT_Q_A = 3.64            #: SECONDARY (back-compat): q = sigma(Q_A * depth_tail + Q_B)
DEFAULT_Q_B = -10.86

# Call-strength gate (PRIMARY, 2026-07; docs/call_strength_recovery.md, adjudicator_call_strength_plan.md,
# eval_corpus/STRENGTH_STRESS_FINDINGS.md). z_excess is a COUNT statistic: fooled by loss genomes that make
# MANY mediocre composition FPs (high count -> inflated z) and blind to divergent bearers with a FEW genuinely
# strong U12s (low count -> low z, in the INCONCLUSIVE gap or below). The un-clipped ensemble MARGIN preserves
# call STRENGTH that both z_excess and the Platt-saturated P_motif discard. cs_p95 = the 95th percentile of the
# call margins: a loss's calls top out ~bounded (composition can only mimic a U12 so far), a bearer's real U12s
# reach much higher; p95 (not max) discounts a lone relic-strength call.
#
# The gate is on the POINT cs_p95, NOT a bootstrap CI lower bound. A lower-bound gate self-calibrates to call
# count, but that is exactly wrong here: divergent bearers ARE few-call-but-strong, so the CI penalizes the
# signal we recover. Trichinella nativa has 5 strong calls (point 6.24) yet its lo@10 is 2.24 == the relic-loss
# Micromonas (2.08) -> a CI-lower-bound gate cannot separate "few but strong" (bearer) from "few and lucky"
# (relic loss) and suppresses the bearers (lo@10 recovers 8/14 vs point's 10/14 at the same 0 FP). Robustness
# against relic-loss inflation instead comes from three things that DON'T penalize the strong-but-sparse signal:
# p95 itself (ignores a lone above-95th relic), cs_min_calls>=3 (kills the degenerate k=2 case, e.g. Micromonas),
# and a fat threshold margin (5.0 sits ~1.36 above the k>=3 loss ceiling ~3.64 -> no real loss's p95 reaches it).
# Calibrated on the 58-bearer/38-loss/25-clade eval panel (index-labeled, snRNA-adjudicated): point cs_p95>=5.0
# recovers 10/14 z-missed divergent bearers at 0 loss FP, FLAT across [4.25,5.25]. Injection stress: a single
# spurious strong call must exceed ~5.5 (and only lifts low-k losses); two coordinated strong calls defeat any
# peak gate (-> downstream snRNA/IPA corroboration on DETECTED). cs_p95_lo/hi are still computed + logged as a
# per-genome confidence annotation (NOT the gate). See memory call-strength-divergent-recovery.
DEFAULT_CS_POINT_THRESHOLD = 5.0  #: cs_p95 (POINT p95 of call margins) >= this (AND >= cs_min_calls) => DETECTED
                                  #: — the null-free ABSOLUTE co-fallback (2026-07-03b: no longer the primary
                                  #: strength driver; p_gumbel_p95 is). Retained for robustness at marginal EVT fits.
DEFAULT_P_GUMBEL_THRESHOLD = 0.01 #: PRIMARY strength driver (2026-07-03b): DETECTED when the per-genome Gumbel
                                  #: tail-prob that the U2 background's OWN max reaches the call p95 is <= this
                                  #: (a ~1% one-sided outlier test on each genome's own U2 null; empirically ==
                                  #: the cs_p95>=5.0 line, but composition-ADAPTIVE — see docs/adjudicator_strength_evt.md)
DEFAULT_CS_MIN_CALLS = 3          #: below this many calls p95==max (outlier-sensitive) -> strength gate off
DEFAULT_CS_P95_PCT = 95.0         #: percentile of the (sorted) call margins defining cs_p95
DEFAULT_CS_P95_THRESHOLD = 4.50   #: LEGACY (pre-2026-07-03 CI-lower-bound gate on cs_p95_lo); no longer gates,
                                  #: retained only so bundles stamped with the old param load without error


def _sigmoid(z):
    return 1.0 / (1.0 + np.exp(-z))


class AdjStatus(str, Enum):
    """Operational status of one genome's adjudication (every result carries exactly one)."""
    ADJUDICATED = "ADJUDICATED"          #: full assessment succeeded (q + CI valid)
    UNDETERMINED = "UNDETERMINED"        #: q's 95% CI straddles the decision boundary -> bearer/loss unresolved
    LOW_N = "LOW_N"                      #: too few calls or species U2 to assess -> file side falls to q_eff=1
    DEGENERATE_TAIL = "DEGENERATE_TAIL"  #: U2 scale (MAD) ~0 / non-finite -> cannot standardize -> unassessable
    SCHEMA_FAIL = "SCHEMA_FAIL"          #: required inputs missing / non-finite / shape mismatch


class MotifCategory(str, Enum):
    """PRIMARY within-genome motif-population gate from the z_excess statistic. Names the MOTIF EVIDENCE
    for a U12 population, NOT the biological bearer/loss truth (that's the corroborated species ``u12_status``
    = motif x snRNA, a database-layer concept: U12_POSITIVE / U12_NEGATIVE / CONFLICT). ``DETECTED`` ==
    strong evidence of a U12 population, 'by motif — corroborate downstream' (see the module CAVEAT);
    ``NOT_DETECTED`` is NOT a loss call (motif-silent divergent bearers land here)."""
    DETECTED = "DETECTED"            #: z_excess >= bearer_floor_z -> a clear U12-motif population
    INCONCLUSIVE = "INCONCLUSIVE"    #: in the empirical gap -> motif can't decide bearer vs loss (undecidable, not "weak bearer")
    NOT_DETECTED = "NOT_DETECTED"    #: z_excess <= loss_ceiling_z -> calls consistent with the U2 background
    UNASSESSABLE = "UNASSESSABLE"    #: z_excess could not be computed (LOW_N / EVT-unfit / fail)


def classify_motif_category(z_excess: float, cs_p95: float, p_gumbel_p95: float, n_calls: int,
                            params: "AdjudicatorParams") -> "MotifCategory":
    """Empirical gap gate (COUNT) + call-strength gate (STRENGTH). ``z_excess`` places the genome against the
    frozen class extremes; the gap between them is the out-of-support / abstain (INCONCLUSIVE) zone. The
    call-strength gate ADDITIONALLY promotes a genome whose strongest calls are strong to DETECTED, resolving
    exactly the divergent bearers z_excess misses (a FEW strong real U12s -> low count -> low z, in the gap or
    below).

    STRENGTH driver = a PER-GENOME anomaly test on the call UPPER TAIL (2026-07-03b). ``p_gumbel_p95`` is the
    Gumbel tail-prob that the genome's OWN U2 background produces a maximum as strong as the call ``cs_p95`` — a
    ~1% one-sided outlier test referenced to each genome's own U2 null (``p_gumbel_p95 <= p_gumbel_threshold``).
    It is evaluated at the p95 (the upper tail), NOT the call median: divergent bearers are few-strong-among-many,
    so their median call is buried in the U2 background (the median-based ``excess_z``/``depth_tail`` miss them);
    the signal lives entirely in the tail (recovery 12->36/46 median->p95; see docs/adjudicator_strength_evt.md).
    The absolute POINT ``cs_p95 >= cs_point_threshold`` is retained as a null-free CO-FALLBACK (OR): the two agree
    on the panel but ``cs_p95`` stays stable at marginal EVT fits where the per-genome null is noisy. Both are
    behind the ``z_excess`` finiteness guard, so a genome that cannot reference its own U2 tail (EVT unfit / -q
    low-N) stays UNASSESSABLE either way. Relic-loss robustness: p95 (not max) + ``n_calls >= cs_min_calls``.
    Same DETECTED caveat: 'a U12-motif population BY MOTIF, corroborate downstream'. See
    docs/adjudicator_strength_evt.md + docs/call_strength_recovery.md + eval_corpus/STRENGTH_STRESS_FINDINGS.md."""
    if not np.isfinite(z_excess):
        return MotifCategory.UNASSESSABLE
    if z_excess >= params.bearer_floor_z:
        return MotifCategory.DETECTED
    if params.cs_gate_enabled and n_calls >= params.cs_min_calls and (
            (np.isfinite(p_gumbel_p95) and p_gumbel_p95 <= params.p_gumbel_threshold)   # PRIMARY: per-genome null
            or (np.isfinite(cs_p95) and cs_p95 >= params.cs_point_threshold)):           # CO-FALLBACK: null-free
        return MotifCategory.DETECTED
    if z_excess <= params.loss_ceiling_z:
        return MotifCategory.NOT_DETECTED
    return MotifCategory.INCONCLUSIVE


@dataclass(frozen=True)
class AdjudicatorParams:
    """Calibrated constants + method settings for the species adjudicator (bundle-stamped in production)."""
    platt_a: float = DEFAULT_PLATT_A    #: margin -> P_motif Platt slope
    platt_c: float = DEFAULT_PLATT_C    #: margin -> P_motif Platt intercept
    loss_ceiling_z: float = DEFAULT_LOSS_CEILING_Z  #: z_excess <= this -> NOT_DETECTED (PRIMARY gate)
    bearer_floor_z: float = DEFAULT_BEARER_FLOOR_Z  #: z_excess >= this -> DETECTED motif population (PRIMARY gate)
    q_a: float = DEFAULT_Q_A            #: SECONDARY/back-compat: depth_tail -> q logistic slope
    q_b: float = DEFAULT_Q_B            #: SECONDARY/back-compat: depth_tail -> q logistic intercept
    call_threshold: float = 0.90        #: P_motif >= this defines a "call"
    u2_threshold: float = 0.50          #: P_motif < this defines the U2 reference set
    tail_pct: float = 99.9              #: U2 percentile defining the tail edge for depth_tail (primary)
    evt_tail_pct: float = 90.0          #: U2 percentile above which the exponential tail is fit (secondary excess_z)
    evt_tail_frac: float = 0.10         #: 1 - evt_tail_pct/100 (fraction of U2 above the EVT tail threshold)
    min_evt_excesses: int = 20          #: below this many U2 exceedances -> secondary EVT unavailable (not a failure)
    bootstrap_n: int = 3000             #: bootstrap resamples for q's CI
    ci: tuple = (2.5, 97.5)             #: q CI percentiles
    ci_smooth_floor_frac: float = 0.3   #: smoothed-bootstrap kernel bandwidth floor, as a fraction of MAD_U2
                                        #: (prevents a degenerate width-0 CI on tied/few calls; count-aware)
    min_calls: int = 1                  #: below this -> LOW_N (no population to assess)
    min_u2: int = 200                   #: below this species U2 count -> LOW_N (cannot reference the tail).
                                        #: Lower it (config/CLI) to let smaller genomes self-adjudicate with
                                        #: their OWN (noisier) U2 tail — the opt-in low-N suppression path.
    random_seed: int = 0
    p_gumbel_threshold: float = DEFAULT_P_GUMBEL_THRESHOLD  #: PRIMARY strength driver: p_gumbel_p95 <= this => DETECTED
                                                        #: (per-genome Gumbel outlier test on the call p95; ~1% one-sided)
    cs_point_threshold: float = DEFAULT_CS_POINT_THRESHOLD  #: null-free CO-FALLBACK: POINT cs_p95 >= this (+ >= cs_min_calls) => DETECTED
    cs_min_calls: int = DEFAULT_CS_MIN_CALLS            #: min calls for the strength gate (below it p95==max, outlier-prone)
    cs_p95_pct: float = DEFAULT_CS_P95_PCT              #: percentile of the call margins defining cs_p95
    cs_ci_lo_pct: float = 10.0                          #: one-sided lower-bound percentile of the cs_p95 bootstrap
                                                        #: (LOGGED as cs_p95_lo, a per-genome confidence annotation;
                                                        #: NOT the gate since 2026-07-03 — the point cs_p95 gates)
    cs_p95_threshold: float = DEFAULT_CS_P95_THRESHOLD  #: LEGACY (old CI-lower-bound gate on cs_p95_lo); unused, kept
                                                        #: so old-stamped bundles load. The gate is cs_point_threshold.
    cs_gate_enabled: bool = True                        #: enable the call-strength DETECTED path (2026-07)
    params_version: str = ADJUDICATOR_PARAMS_VERSION

    @classmethod
    def from_dict(cls, d: Optional[dict]) -> "AdjudicatorParams":
        if not d:
            return cls()
        known = set(cls.__dataclass_fields__)
        return cls(**{k: v for k, v in d.items() if k in known})


@dataclass
class AdjudicatorResult:
    """Outcome of adjudicating one genome."""
    n_call: int                          #: number of P_motif-calls
    n_u2: int                            #: number of U2-reference introns
    depth_tail: float                    #: SECONDARY (back-compat): call-core depth beyond U2's 99.9th pct, MAD units
    q: float                             #: SECONDARY (back-compat): depth_tail->q P(bearer) point estimate
    q_lo: float                          #: SECONDARY: q lower CI (bootstrap)
    q_hi: float                          #: SECONDARY: q upper CI (bootstrap)
    status: AdjStatus = AdjStatus.ADJUDICATED
    z_excess: float = float("nan")       #: PRIMARY population statistic (Poisson count-excess over U2 tail)
    cs_p95: float = float("nan")         #: 95th-pct of the un-clipped call MARGINS — the null-free co-fallback gate value
    p_gumbel_p95: float = float("nan")   #: PRIMARY strength driver (2026-07-03b): per-genome Gumbel tail-prob at cs_p95
                                         #: (U2 background's own-max prob of reaching the call p95; <= p_gumbel_threshold => DETECTED)
    cs_p95_lo: float = float("nan")      #: cs_p95 bootstrap CI LOWER bound — logged confidence annotation (NOT the gate)
    cs_p95_hi: float = float("nan")      #: cs_p95 bootstrap CI upper bound — logged confidence annotation
    motif_category: MotifCategory = MotifCategory.UNASSESSABLE  #: PRIMARY motif gate (DETECTED/INCONCLUSIVE/NOT_DETECTED)
    excess_z: float = float("nan")       #: SECONDARY (size-aware significance): call-core beyond U2 expected-max
    p_gumbel: float = float("nan")       #: SECONDARY: Gumbel tail-prob that U2's own max reaches the call-core
    z_expmax: float = float("nan")       #: SECONDARY: U2 expected-max in robust-z units (the size-aware cutoff)
    secondary_available: bool = False    #: True when the EVT secondary (excess_z/p_gumbel) could be computed
    n_adj_called: Optional[int] = None   #: final P_adj>=0.5 U12 calls (set by the file-side application; None otherwise)
    n_motif_called: Optional[int] = None #: UNGATED motif calls (P_motif>=0.5, DISREGARDING motif_category) —
                                         #: == n_adj_called except for NOT_DETECTED genomes, where the gate
                                         #: suppresses these to u2 (file-side; None otherwise)
    params: AdjudicatorParams = field(repr=False, default_factory=AdjudicatorParams)

    @property
    def assessable(self) -> bool:
        """Back-compat: a usable q was produced (confident or undetermined, but not a hard failure)."""
        return self.status in (AdjStatus.ADJUDICATED, AdjStatus.UNDETERMINED)


# --------------------------------------------------------------------------------------------------
# Margin / P_motif helpers
# --------------------------------------------------------------------------------------------------
def ensemble_margin(ensemble_models, X: np.ndarray) -> np.ndarray:
    """Mean ``decision_function`` (the SVM margin) over the ensemble's calibrated sub-estimators.

    The bundle stores ``CalibratedClassifierCV`` models whose isotonic probability saturates at 0/1; the
    *margin* preserves the headroom the adjudicator needs. Averages over all (model x cv-fold) base SVCs.
    """
    subs = [sub.estimator for m in ensemble_models for sub in m.model.calibrated_classifiers_]
    acc = np.zeros(len(X))
    for est in subs:
        acc += est.decision_function(X)
    return acc / len(subs)


def p_motif_from_margin(margin: np.ndarray, params: AdjudicatorParams) -> np.ndarray:
    return _sigmoid(params.platt_a * np.asarray(margin, float) + params.platt_c)


# --------------------------------------------------------------------------------------------------
# U2 reference + the two separation statistics
# --------------------------------------------------------------------------------------------------
@dataclass
class _U2Reference:
    """Robust U2 location/scale + tail edges. ``z_expmax``/``evt_beta`` are NaN when the EVT tail can't be fit
    (too few exceedances) — the PRIMARY ``depth_tail`` does NOT depend on them, only the SECONDARY excess_z does."""
    med: float
    mad: float
    q_tail: float          #: the ``tail_pct`` (99.9th) U2 percentile — the depth reference (secondary depth_tail)
    z_expmax: float        #: U2 expected-max in robust-z units (secondary); NaN if EVT tail unfit
    evt_beta: float        #: Gumbel scale in robust-z units (secondary); NaN if EVT tail unfit
    finite: bool           #: True when (med, mad, q_tail) are usable (mad > 0 and finite)
    q_evt: float = float("nan")  #: the ``evt_tail_pct`` (90th) U2 percentile — the z_excess tail anchor
    lam: float = float("nan")    #: exponential U2-tail rate above q_evt (for z_excess); NaN if tail unfit
    exp_max_margin: float = float("nan")  #: Gumbel-expected U2 MAX on the raw-margin axis (q_evt + ln(frac*n_u2)/lam);
                                          #: NaN if tail unfit. Same quantity behind z_expmax; retained for the tail plot.


def _u2_reference(u2_margin: np.ndarray, n_u2: float, params: AdjudicatorParams) -> _U2Reference:
    """Robust (med, MAD), the 99.9th-pct tail edge (primary), and the EVT expected-max (secondary).

    PRIMARY ``depth_tail`` needs only (med, MAD, q_tail) — non-parametric in the tail SHAPE. The SECONDARY
    EVT (``z_expmax``/Gumbel) fits an exponential (xi=0) tail above ``evt_tail_pct`` by method-of-moments;
    that fit is fragile (needs >= ``min_evt_excesses`` exceedances) and is the labelled secondary only.
    """
    med = float(np.median(u2_margin))
    mad = 1.4826 * float(np.median(np.abs(u2_margin - med)))
    q_tail = float(np.percentile(u2_margin, params.tail_pct))

    finite = bool(np.isfinite(med) and np.isfinite(mad) and np.isfinite(q_tail) and mad > 1e-9)
    if not finite:
        return _U2Reference(med, mad, q_tail, np.nan, np.nan, finite=False)

    # ---- exponential U2-tail fit above the 90th pct: drives BOTH the PRIMARY z_excess (Poisson count-
    # excess over what this tail predicts) and the SECONDARY EVT expected-max. Fragile (needs
    # >= min_evt_excesses exceedances); when unfit, z_excess/excess_z are NaN -> the gate reads UNASSESSABLE.
    #
    # WHY EXPONENTIAL (xi=0) AND NOT A "BETTER-FITTING" TAIL (GPD/Weibull). The U2 margin tail IS mildly
    # lighter-than-exponential (cv_excess~0.87, GPD xi~-0.15 corpus-wide), so a GPD/Weibull fits the OBSERVED
    # tail better — yet both make z_excess WORSE (faithful gold-panel head-to-head, 1839 bearers/45 losses:
    # GPD/Weibull detonate divergent losses, e.g. Pristionchus 1.9->51; AUC 0.990->0.976; recover FEWER bearers
    # at the 0-loss-FP operating point; see eval_corpus/tail_model_proveout_2026-07/). This is not a paradox —
    # the exponential is not a tail DESCRIPTION, it is the NULL for a count-excess test, and it is the right
    # null from first principles (docs/adjudicator_design.md; memory tail-model-alternatives-closed):
    #   (1) EXTRAPOLATION, not fit: `call_core` sits far BEYOND the observed U2 data, so in-sample fit quality
    #       is irrelevant; given only the reliably-estimable mean excess, the exponential is the MAX-ENTROPY
    #       (least-committal, constant-hazard) extrapolation. A GPD xi<0 asserts an accelerating hazard / finite
    #       endpoint out where there is no data to justify it.
    #   (2) CONSERVATIVE by construction: the exponential OVER-predicts the null count (reality is lighter) so z
    #       UNDER-declares excess; a lighter/bounded tail sends pred_cnt->0, so z collapses to a raw COUNT — the
    #       loss-FP failure mode. Under the asymmetric loss (a false divergent-loss CALL is expensive; a
    #       conservative miss is cheap + rescuable by the strength gate) that conservative bias is the feature.
    #   (3) CROSS-GENOME COMPARABILITY: z is gated against a FIXED anchor (loss_ceiling/bearer_floor), which only
    #       works if z is stable across genomes. The 1-parameter mean-excess has low estimation variance; GPD's
    #       shape xi is high-variance at per-genome N, so its extrapolated count swings erratically by genome
    #       (why the losses explode NON-uniformly) and breaks the fixed-anchor comparison.
    # Boundary: this is NOT "simpler is always better" — for a per-intron p-value WITHIN the observed tail the
    # GPD would be correct. The exponential wins HERE from the conjunction {far-extrapolation, count-null,
    # asymmetric loss, fixed cross-genome anchor}; change any of those and re-derive.
    q_evt = float(np.percentile(u2_margin, params.evt_tail_pct))
    exc = u2_margin[u2_margin > q_evt] - q_evt
    z_expmax = np.nan
    evt_beta = np.nan
    lam = np.nan
    exp_max_margin = np.nan
    if len(exc) >= params.min_evt_excesses and exc.mean() > 0:
        lam = 1.0 / float(exc.mean())
        exp_max_margin = q_evt + np.log(max(params.evt_tail_frac * n_u2, 2.0)) / lam
        z_expmax = (exp_max_margin - med) / mad
        evt_beta = 1.0 / (lam * mad + 1e-9)
    return _U2Reference(med, mad, q_tail, z_expmax, evt_beta, finite=True, q_evt=q_evt, lam=lam,
                        exp_max_margin=float(exp_max_margin))


def compute_depth_tail(call_core_margin: float, ref: _U2Reference) -> float:
    """SECONDARY (back-compat): how far the call-core sits beyond U2's 99.9th pct, in MAD units. Count-blind
    (a central-tendency statistic) -> mis-ranks numerous-but-shallow divergent bearers. Not the gate driver."""
    return float((call_core_margin - ref.q_tail) / ref.mad)


def compute_z_excess(call_core_margin: float, calls: np.ndarray, n_u2: int, n_total: int,
                     ref: _U2Reference, params: AdjudicatorParams) -> float:
    """PRIMARY population statistic: Poisson significance of the strong-call COUNT above what the genome's
    OWN U2 tail predicts in the call region. Loss: the calls ARE the U2 tail (excess ~ 0 / negative).
    Bearer: a separate population the U2 tail can't explain (large positive excess). A COUNT (so it sees
    numerous-but-shallow divergent bearers depth_tail misses), referenced to the per-genome background (so a
    big loss reads ~0 and a small loss can't spike, unlike raw count / density). NaN if the U2 tail is unfit.
    Matches eval_corpus/robust_sep_one_nu2.py (the recipe the frozen anchors were calibrated with)."""
    if not (np.isfinite(ref.lam) and np.isfinite(ref.q_evt)) or n_total <= 0:
        return float("nan")
    frac_u2 = n_u2 / n_total
    s_u2 = frac_u2 * params.evt_tail_frac * np.exp(-ref.lam * (call_core_margin - ref.q_evt))
    pred_cnt = s_u2 * n_total                          # U2-predicted #calls above the call-core
    obs_cnt = float((calls > call_core_margin).sum())  # observed #calls above the call-core
    return float((obs_cnt - pred_cnt) / np.sqrt(pred_cnt + 1.0))


def q_from_depth_tail(depth_tail, params: AdjudicatorParams) -> np.ndarray:
    return _sigmoid(params.q_a * np.asarray(depth_tail, float) + params.q_b)


def compute_excess_z(call_core_margin: float, ref: _U2Reference) -> float:
    """SECONDARY size-aware significance: call-core beyond the U2 EXPECTED MAXIMUM, in robust-z units.

    NaN when the EVT tail could not be fit (``ref.z_expmax`` is NaN). Not the q-driver — diagnostic only.
    """
    robust_z = (call_core_margin - ref.med) / ref.mad
    return float(robust_z - ref.z_expmax)


def compute_p_gumbel(call_core_margin: float, ref: _U2Reference) -> float:
    """SECONDARY: Gumbel tail-prob that U2's own max reaches the call-core (NaN if EVT tail unfit)."""
    if not np.isfinite(ref.z_expmax) or not np.isfinite(ref.evt_beta):
        return float("nan")
    robust_z = (call_core_margin - ref.med) / ref.mad
    return float(1.0 - np.exp(-np.exp(-(robust_z - ref.z_expmax) / max(ref.evt_beta, 1e-6))))


def background_fdr(margin: np.ndarray, p_motif: np.ndarray, params: AdjudicatorParams) -> np.ndarray:
    """REPORT-ONLY per-intron background FDR: of the strong calls at least this strong, what fraction does
    the genome's OWN U2 tail already explain? A **third axis** — background-relative surprisal — deliberately
    NOT folded into ``P_motif`` or ``adjusted_score`` (see ``docs/adjudicator_qdriver_postmortem.md`` §3d:
    a species signal must not be multiplied into a calibrated per-intron probability, and the single-column
    ``rel_score > 0`` HC filter must keep working). It gates nothing.

    Definition (the per-intron decomposition of ``z_excess``, which evaluates the SAME arithmetic once at the
    call core). For each call with margin ``m_i``, using the identical deterministic ``_u2_reference`` fit:

        E_i    = evt_tail_frac * n_u2 * exp(-lam * (m_i - q_evt))   # U2 introns expected at or above m_i
        rank_i = #{j : m_j >= m_i}                                  # calls observed at or above m_i
        bg_fdr = min(E_i / rank_i, 1)                               # Benjamini-Hochberg q estimate

    with the standard BH step-up applied (a running minimum from the weakest call upward), so the column is
    MONOTONE in margin: a stronger call can never carry a worse ``bg_fdr`` than a weaker one. Without it the
    raw ratio is non-monotone (E and rank both fall with margin), which reads as noise in the output.

    Read it as: ``bg_fdr ~ 1`` -> indistinguishable from this genome's background; ``bg_fdr << 1`` -> the
    background cannot account for it. Computed only for the adjudicator's canonical call set
    (``P_motif >= call_threshold``) — the fitted domain of the exponential tail, and the only set for which
    the question is meaningful (non-calls ARE the background; their ranking is ``P_motif``). Populated even
    for a ``NOT_DETECTED`` genome, where the calling columns are zeroed: that is the point, since it
    separates a genome whose calls are uniformly background (Symbiodinium: every ``bg_fdr`` = 1) from one
    holding a genuine outlier (Phytophthora infestans: one call at 6.9e-3) — currently indistinguishable.

    ``rank_i`` uses a tie-max rule (``searchsorted``), so the result depends only on the multiset of margins,
    never on row order — required for streaming/in-memory bit-identical parity.

    CAVEAT (inherited, not new): the tail is fitted on the U2 side only, so the calls are missing from the
    exceedance set and ``lam`` is biased steep when they are a large share of it. Negligible at production
    ratios (Symbiodinium: 109 calls vs a 51,621-strong exceedance set = 0.2%), material only if calls
    approach ~10% of ``evt_tail_frac * n_u2``. ``z_excess`` uses the identical reference, so this adds no
    new exposure — but a genome that thin is already ``LOW_N``/``UNASSESSABLE`` territory.

    Returns a float array aligned to ``margin``; NaN for non-calls, unscorable rows, and whenever the U2 tail
    cannot be referenced (too few U2 / degenerate scale / EVT unfit), matching the adjudicator's own guards.
    """
    margin = np.asarray(margin, float)
    p_motif = np.asarray(p_motif, float)
    out = np.full(len(margin), np.nan)
    good = np.isfinite(margin) & np.isfinite(p_motif)
    u2 = margin[good & (p_motif < params.u2_threshold)]
    if len(u2) < params.min_u2:
        return out
    ref = _u2_reference(u2, len(u2), params)
    if not (ref.finite and np.isfinite(ref.lam) and np.isfinite(ref.q_evt)):
        return out
    is_call = good & (p_motif >= params.call_threshold)
    m = margin[is_call]
    if not len(m):
        return out
    e = params.evt_tail_frac * len(u2) * np.exp(-ref.lam * (m - ref.q_evt))
    # rank_i = #{j : m_j >= m_i}; ties share the max rank so the value is order-independent.
    rank = np.searchsorted(-np.sort(m)[::-1], -m, side="right")
    raw = e / rank
    # BH step-up: running minimum from the weakest call upward, so bg_fdr is monotone in margin.
    desc = np.argsort(-m, kind="stable")                 # strongest first
    adjusted = np.minimum.accumulate(raw[desc][::-1])[::-1]
    fdr = np.empty_like(raw)
    fdr[desc] = adjusted
    out[is_call] = np.minimum(fdr, 1.0)
    return out


# --------------------------------------------------------------------------------------------------
# The adjudicator
# --------------------------------------------------------------------------------------------------
def _fail(n_call, n_u2, status, params) -> AdjudicatorResult:
    """Construct a non-assessable result (q defaults to the loss side; CI collapsed). motif_category stays
    the default UNASSESSABLE (z_excess could not be computed): for SCHEMA_FAIL / DEGENERATE_TAIL, and for
    LOW_N when the species has too few U2 to reference its own tail (genuinely under-powered)."""
    return AdjudicatorResult(n_call=int(n_call), n_u2=int(n_u2), depth_tail=float("nan"),
                             q=0.0, q_lo=0.0, q_hi=0.0, status=status,
                             secondary_available=False, params=params)


def _no_population(n_u2, params) -> AdjudicatorResult:
    """Well-powered genome (>= ``min_u2`` U2 introns) with ZERO strong calls -> a definite ``NOT_DETECTED``,
    NOT ``UNASSESSABLE``. An ample U2 background with no U12-motif population is the CLEAREST loss: z_excess is
    0 by construction (no calls to exceed the U2 tail), so the gate verdict is unambiguous. The secondary
    q-pipeline is still degenerate (no calls to bootstrap) -> operational ``status=LOW_N``, but the PRIMARY
    ``motif_category`` is decisive. (Before zexcess_gap_2026-06-30b this fell through ``_fail`` -> UNASSESSABLE,
    wrongly pooling well-powered losses with truly under-powered low-U2 genomes.)"""
    return AdjudicatorResult(n_call=0, n_u2=int(n_u2), depth_tail=float("nan"),
                             q=0.0, q_lo=0.0, q_hi=0.0, status=AdjStatus.LOW_N,
                             z_excess=0.0, motif_category=MotifCategory.NOT_DETECTED,
                             secondary_available=False, params=params)


def _assess(calls: np.ndarray, n_u2: int, n_total: int, ref: _U2Reference,
            params: AdjudicatorParams) -> AdjudicatorResult:
    """Score a call set against the species U2 reference. PRIMARY: z_excess (Poisson count-excess) ->
    the empirical gap gate (motif_category). SECONDARY/back-compat: depth_tail -> q + bootstrap CI + status
    (ADJUDICATED if the bootstrap q-CI clears 0.5, else UNDETERMINED) and the EVT excess_z/p_gumbel."""
    call_core = float(np.median(calls))
    depth_tail = compute_depth_tail(call_core, ref)
    if not np.isfinite(depth_tail):
        return _fail(len(calls), n_u2, AdjStatus.DEGENERATE_TAIL, params)
    q_point = float(q_from_depth_tail(depth_tail, params))

    # PRIMARY point statistics: z_excess (count) + cs_p95 (call STRENGTH = the cs_p95_pct percentile of the
    # call MARGINS; calls is already sorted). The un-clipped margin retains headroom z_excess/P_motif discard;
    # the p95 (not max) is robust to a lone relic FP. The call-strength gate below is on this POINT cs_p95.
    z_excess = compute_z_excess(call_core, calls, n_u2, n_total, ref, params)
    cs_p95 = float(np.percentile(calls, params.cs_p95_pct))
    # PRIMARY strength driver: the same Gumbel EVT null as p_gumbel/excess_z, but evaluated at the call UPPER
    # TAIL (cs_p95) rather than the median (call_core). The upper tail is where the few-strong divergent-bearer
    # signal lives; the median is buried in the U2 background (see docs/adjudicator_strength_evt.md).
    p_gumbel_p95 = compute_p_gumbel(cs_p95, ref)

    excess_z = compute_excess_z(call_core, ref)
    p_gumbel = compute_p_gumbel(call_core, ref)
    secondary_available = bool(np.isfinite(excess_z))

    # COUNT-AWARE SMOOTHED bootstrap -> ONE resample serves BOTH CIs. A plain bootstrap of a median/percentile
    # is degenerate for tied/clustered calls (the resampled statistic takes few distinct values -> a width-0 CI
    # that silently defeats the confidence bands in exactly the low-call-count regime they exist for). Smoothing
    # each resampled call by a kernel (bandwidth = the robust call spread, floored to a fraction of the U2 scale)
    # (a) breaks the tie degeneracy and (b) makes BOTH CIs widen as the call count shrinks (sampling SD ~
    # bandwidth/sqrt(k)). It drives (i) the SECONDARY depth_tail->q CI (ADJUDICATED/UNDETERMINED) and (ii) the
    # cs_p95 CI (cs_p95_lo/hi), now LOGGED as a per-genome confidence annotation but NOT the gate (the point
    # cs_p95 gates; see the module comment on why a CI lower bound structurally suppresses divergent bearers).
    # Point estimates (q_point, cs_p95) use the UN-smoothed calls; smoothing only affects the CIs.
    rng = np.random.RandomState(params.random_seed)
    k = len(calls)
    mad_calls = 1.4826 * float(np.median(np.abs(calls - call_core)))
    bw = max(mad_calls, params.ci_smooth_floor_frac * ref.mad)
    qs = np.empty(params.bootstrap_n)
    cs_boot = np.empty(params.bootstrap_n)
    for b in range(params.bootstrap_n):
        samp = calls[rng.randint(0, k, k)] + rng.normal(0.0, bw, k)
        qs[b] = q_from_depth_tail(compute_depth_tail(float(np.median(samp)), ref), params)
        cs_boot[b] = float(np.percentile(samp, params.cs_p95_pct))
    q_lo, q_hi = (float(x) for x in np.percentile(qs, params.ci))
    # cs_p95_lo/hi = ONE-SIDED CI bounds (cs_ci_lo_pct, default 10% = 90% one-sided). Logged for provenance;
    # the gate is on the POINT cs_p95 (a CI lower bound penalizes call count = the divergent-bearer signal).
    cs_p95_lo = float(np.percentile(cs_boot, params.cs_ci_lo_pct))
    cs_p95_hi = float(np.percentile(cs_boot, 100.0 - params.cs_ci_lo_pct))

    # PRIMARY gate: z_excess (count) gap gate + call-strength gate — per-genome p_gumbel_p95 (primary) OR the
    # null-free POINT cs_p95 (co-fallback).
    motif_category = classify_motif_category(z_excess, cs_p95, p_gumbel_p95, len(calls), params)

    status = AdjStatus.ADJUDICATED if (q_lo > 0.5 or q_hi < 0.5) else AdjStatus.UNDETERMINED
    return AdjudicatorResult(
        n_call=int(len(calls)), n_u2=int(n_u2), depth_tail=float(depth_tail),
        q=q_point, q_lo=q_lo, q_hi=q_hi, status=status,
        z_excess=float(z_excess), cs_p95=cs_p95, p_gumbel_p95=float(p_gumbel_p95),
        cs_p95_lo=cs_p95_lo, cs_p95_hi=cs_p95_hi,
        motif_category=motif_category,
        excess_z=float(excess_z), p_gumbel=float(p_gumbel), z_expmax=float(ref.z_expmax),
        secondary_available=secondary_available, params=params)


def adjudicate(margin: np.ndarray, p_motif: np.ndarray,
               params: Optional[AdjudicatorParams] = None) -> AdjudicatorResult:
    """Run the species adjudicator on one genome's per-intron (margin, P_motif).

    Returns ``q`` = P(bearer) + bootstrap CI, driven by the size-invariant ``depth_tail``. ``q``'s CI widens
    for few/borderline calls (low-numbers uncertainty), so it stays the right confidence signal for
    arbitrarily-low-N bearers WITHOUT a count threshold (count enters only as CI width). The status field
    distinguishes a confident call (ADJUDICATED) from a CI that straddles the boundary (UNDETERMINED) and
    from the degenerate/low-data branches (LOW_N / DEGENERATE_TAIL / SCHEMA_FAIL) — never a silent NaN.

    When a species has fewer than ``min_u2`` U2 introns it cannot reference its own tail -> LOW_N (the
    file-side layer then defaults to ``q_eff = 1``, no suppression: option A, the safe default). Lowering
    ``min_u2`` (config/CLI) lets smaller genomes self-adjudicate against their OWN (noisier) U2 tail — the
    opt-in low-N suppression path (option B); a pooled GLOBAL U2 reference was evaluated and rejected (it
    over-suppresses clean-background genomes — see docs/raw_gated_scoring.md §0c).
    """
    params = params or AdjudicatorParams()

    # ---- schema / NaN guards ----
    margin = np.asarray(margin, float)
    p_motif = np.asarray(p_motif, float)
    if margin.shape != p_motif.shape or margin.ndim != 1 or len(margin) == 0:
        return _fail(0, 0, AdjStatus.SCHEMA_FAIL, params)
    good = np.isfinite(margin) & np.isfinite(p_motif)
    margin, p_motif = margin[good], p_motif[good]
    if len(margin) == 0:
        return _fail(0, 0, AdjStatus.SCHEMA_FAIL, params)

    # SORT the call set into a canonical order so the bootstrap CI is INDEPENDENT of the input row order.
    # The bootstrap resamples calls by index (rng.randint); without sorting, the streaming (per-contig) and
    # in-memory paths emit the same call *values* in different *order*, so a fixed-seed bootstrap draws
    # different elements -> different q CI (the point q is an order-independent median, so it matched, but
    # P_adj_lo/P_adj_hi diverged across the two paths). Sorting makes the whole adjudication order-invariant.
    calls = np.sort(margin[p_motif >= params.call_threshold])
    u2 = margin[p_motif < params.u2_threshold]
    n_total = len(p_motif)

    # ---- UNASSESSABLE: too few species U2 to reference the tail (genuinely under-powered) ----
    if len(u2) < params.min_u2:
        return _fail(len(calls), len(u2), AdjStatus.LOW_N, params)

    # ---- NOT_DETECTED: well-powered (>= min_u2 U2) but ZERO strong calls -> the clearest loss, not
    # "couldn't assess". A genome with an ample U2 background and no U12-motif population is a definite
    # NOT_DETECTED (z_excess 0 by construction); only the truly under-powered case above is UNASSESSABLE. ----
    if len(calls) < params.min_calls:
        return _no_population(len(u2), params)

    n_u2 = len(u2) / n_total * n_total   # full-population U2 count (== len(u2) unless subsampled upstream)
    ref = _u2_reference(u2, n_u2, params)

    # ---- DEGENERATE_TAIL: species U2 scale collapsed -> cannot standardize ----
    if not ref.finite:
        return _fail(len(calls), len(u2), AdjStatus.DEGENERATE_TAIL, params)

    return _assess(calls, len(u2), n_total, ref, params)


# --------------------------------------------------------------------------------------------------
# File-side post-process (the inference hook; mirrors species_gate.apply_raw_gated_postprocess)
# --------------------------------------------------------------------------------------------------
#: Raw-feature column order the ensemble was trained on (see eval_corpus/robust_sep_one.py). The 6th
#: feature, support2_raw, is the 2nd-largest of the clipped motif log-odds and is derived here.
_RAW_FEATURE_COLS = ("5'_raw", "bp_raw", "3'_raw", "bp_offset", "bp_scan_confidence")


def _effective_q(result: "AdjudicatorResult"):
    """PRIMARY (two-number design): a BINARY species gate from the ``z_excess`` tier — NOT a continuous
    multiplier. ``P_motif`` is the calibrated, post-hoc-thresholdable per-intron number; the species
    confidence is carried *separately* in ``motif_category`` and must NOT be folded into a per-intron
    probability (that conflation -- the superseded ``q*P_motif`` chain -- baked in a threshold and reported
    the confident count off an unidentified logistic slope; see ``docs/adjudicator_qdriver_postmortem.md``).

    - ``NOT_DETECTED`` -> ``q=0`` : suppress (the strong calls are consistent with the U2 tail / relic FPs).
    - everything else (DETECTED / INCONCLUSIVE / UNASSESSABLE) -> ``q=1`` : pass ``P_motif`` through unscaled.
      "P_motif is the call everywhere; suppress only the no-population case" (§0). INCONCLUSIVE/UNASSESSABLE
      calls are reported and flagged via ``motif_category`` for downstream filtering -- never silently boosted
      or deflated. (Note: DETECTED is 'by motif'; a high-z genome with searched-and-absent snRNAs must be
      investigated downstream, not auto-accepted -- see the module CAVEAT.)
    """
    return (0.0, 0.0, 0.0) if result.motif_category == MotifCategory.NOT_DETECTED else (1.0, 1.0, 1.0)


#: Sidecar file suffix for the per-species tail-model figure payload (consumed by
#: ``visualization.plots.tail_model_plot``). One JSON per genome, next to its ``score_info.iic``.
TAIL_MODEL_SIDECAR_SUFFIX = ".tail_model.iic.json"


def tail_model_sidecar_path(score_info_path) -> str:
    """Canonical sidecar path next to a ``score_info.iic`` (``base.score_info.iic`` ->
    ``base.tail_model.iic.json``). Reader (plot layer) and writer (below) derive it identically."""
    s = str(score_info_path)
    suf = ".score_info.iic"
    return (s[: -len(suf)] if s.endswith(suf) else s) + TAIL_MODEL_SIDECAR_SUFFIX


def _jf(v):
    """JSON-safe float: real NaN/inf -> None (JSON has no NaN; keep the file strictly valid)."""
    v = float(v)
    return v if np.isfinite(v) else None


def build_tail_model(margin: np.ndarray, p_motif: np.ndarray,
                     result: "AdjudicatorResult", params: AdjudicatorParams,
                     n_bins: int = 300) -> dict:
    """Assemble the per-species tail-model figure payload from the SAME ``(margin, P_motif)`` the
    adjudicator scored — production-faithful by construction: it re-derives the *deterministic* U2
    reference (``_u2_reference``) on the exact U2 split, so ``q90``/``lam``/``exp_max`` match the gate
    values with no re-fit or subsample drift (the trap that made an eval-cache genome read z=6.6 vs a
    production 4.2). Everything is stored in the native RAW ensemble-margin space; the plotter maps it to
    the ``logit(P_motif)`` axis via ``(platt_a, platt_c)`` so the exponential U2 tail is a straight line
    on log-y and the P=0.9 call line is a fixed reference.

    Returns a JSON-serializable dict. ``assessable`` is False (no fitted curve — bars only) when the U2
    tail can't be referenced (too few U2, degenerate scale, or too few EVT exceedances); the U2 histogram
    and the call margins are still emitted so the population stays visible.
    """
    margin = np.asarray(margin, float)
    p_motif = np.asarray(p_motif, float)
    good = np.isfinite(margin) & np.isfinite(p_motif)
    margin, p_motif = margin[good], p_motif[good]
    u2 = margin[p_motif < params.u2_threshold]
    calls = np.sort(margin[p_motif >= params.call_threshold])   # match adjudicate()'s canonical order

    payload = {
        "params_version": ADJUDICATOR_PARAMS_VERSION,
        "platt_a": float(params.platt_a), "platt_c": float(params.platt_c),
        "call_threshold": float(params.call_threshold), "u2_threshold": float(params.u2_threshold),
        "evt_tail_pct": float(params.evt_tail_pct), "evt_tail_frac": float(params.evt_tail_frac),
        "n_u2": int(len(u2)), "n_call": int(len(calls)),
        "z_excess": _jf(result.z_excess), "cs_p95": _jf(result.cs_p95),
        "p_gumbel_p95": _jf(result.p_gumbel_p95), "motif_category": result.motif_category.value,
        "call_margins": [float(m) for m in calls],   # small (few..few-hundred); exact call positions
    }
    if len(u2) < params.min_u2:
        payload.update(assessable=False, reason="LOW_N_U2")
        return payload

    ref = _u2_reference(u2, len(u2), params)   # deterministic -> bit-exact to the gate's fit
    # x-range: cover the U2 body + any calls + the modeled expected-max, so the extrapolated tail is visible.
    hi = float(u2.max())
    if len(calls):
        hi = max(hi, float(calls.max()))
    if np.isfinite(ref.exp_max_margin):
        hi = max(hi, float(ref.exp_max_margin))
    lo = float(np.percentile(u2, 0.1))
    pad = 0.5 * ref.mad if (ref.finite and np.isfinite(ref.mad)) else 0.0
    edges = np.linspace(lo, hi + pad, n_bins + 1)
    counts, _ = np.histogram(u2, bins=edges)
    payload.update({
        "assessable": bool(ref.finite and np.isfinite(ref.lam)),
        "reason": None if (ref.finite and np.isfinite(ref.lam)) else
                  ("DEGENERATE_TAIL" if not ref.finite else "EVT_UNFIT"),
        "u2_hist_edges": [float(e) for e in edges],
        "u2_hist_counts": [int(c) for c in counts],
        "med": _jf(ref.med), "mad": _jf(ref.mad), "q_tail": _jf(ref.q_tail),
        "q90": _jf(ref.q_evt), "lam": _jf(ref.lam), "exp_max": _jf(ref.exp_max_margin),
    })
    return payload


def apply_pmotif_adjudication(score_info_path, ensemble_models,
                              params: Optional[AdjudicatorParams] = None,
                              messenger=None) -> "AdjudicatorResult":
    """``pmotif_adjudicated`` post-classification: read score_info.iic, compute the per-intron motif
    probability ``P_motif`` from the ensemble MARGIN, run the per-species adjudicator, and write the two
    interpretable per-intron numbers + the call.

    TWO-NUMBER output (per-intron motif probability + species gate). Adds/overwrites in ``score_info.iic``:
      - ``P_motif``     : sigmoid(Platt(margin)) — the calibrated, species-agnostic, post-hoc-thresholdable
                          per-intron motif probability (threshold it at any certainty);
      - ``z_excess``    : the per-species population statistic (constant within a run);
      - ``motif_category``: DETECTED / INCONCLUSIVE / NOT_DETECTED / UNASSESSABLE — the motif-population gate
                          (DETECTED == 'a U12-motif population, by motif — corroborate downstream'; NOT a
                          confirmed bearer, and NOT_DETECTED is NOT a loss call);
      - ``type_id``     : ``u12`` iff ``P_motif >= 0.5`` AND ``motif_category != NOT_DETECTED`` (suppress only
                          the no-population case; INCONCLUSIVE/UNASSESSABLE calls are made and flagged);
      - ``adjusted_score`` = ``100 * q_eff * P_motif`` — the GATED call score (bed col5; ``>= 50`` iff u12);
      - ``rel_score`` = ``adjusted_score - 90`` — the same number recentred on the high-confidence threshold,
                          so ``rel_score > 0`` is a SELF-SUFFICIENT single-column HC filter that needs no
                          companion condition and no knowledge of the run's threshold. That identity is the
                          reason both are gated together: un-gating ``rel_score`` alone would break the
                          identity, and un-gating it at all would make ``rel_score > 0`` select the calls a
                          ``NOT_DETECTED`` genome explicitly does not make (Symbiodinium: 109 of them),
                          forcing every consumer into ``AND type_id == 'u12'``. The UNGATED per-intron
                          ranking is ``P_motif``, which is exactly what it is for.
      - ``bg_fdr``      : REPORT-ONLY background-relative surprisal for the call set (see :func:`background_fdr`).

    Removed 2026-07 (deterministic functions of ``P_motif`` + ``motif_category``, superseded ``q*P_motif``
    chain — see ``docs/adjudicator_qdriver_postmortem.md``): ``q``, ``P_adj``, ``P_adj_lo``, ``P_adj_hi``.
    ``P_adj_lo``/``P_adj_hi`` were bit-identical to ``P_adj`` by construction; ``q`` was a per-species
    constant. Use ``P_motif`` + ``motif_category``.

    Unscorable introns (NaN motif log-odds) keep NaN ``P_motif`` and are not called. Per-run = per-species
    (one genome). Returns the :class:`AdjudicatorResult` (z_excess + tier + secondary q/CI diagnostics).
    """
    import pandas as pd
    params = params or AdjudicatorParams()
    df = pd.read_csv(score_info_path, sep="\t", dtype=str, keep_default_na=False)

    def col(name):
        if name not in df.columns:
            return None
        # Coerce any non-numeric token to NaN. score_info encodes a missing value as "" for most
        # columns but as the literal "NA" for BP-scan fields (bp_offset / bp_scan_confidence) on introns
        # the branch-point scan produced no offset for — a plain `.replace("", "nan").astype(float)` chokes
        # on "NA" (ValueError). bp_offset/bp_scan_confidence NaN -> 0 downstream via np.nan_to_num, so "NA"
        # and "" are equivalent here; to_numeric(coerce) unifies them (and "nan"/"null"/etc.) safely.
        return pd.to_numeric(df[name], errors="coerce").to_numpy()

    feats = {c: col(c) for c in _RAW_FEATURE_COLS}
    if any(feats[c] is None for c in ("5'_raw", "bp_raw", "3'_raw")):
        raise ValueError("pmotif_adjudicated post-process: score_info missing 5'_raw/bp_raw/3'_raw")
    n = len(df)
    r5, rbp, r3 = feats["5'_raw"], feats["bp_raw"], feats["3'_raw"]
    boff = feats["bp_offset"] if feats["bp_offset"] is not None else np.zeros(n)
    bconf = feats["bp_scan_confidence"] if feats["bp_scan_confidence"] is not None else np.zeros(n)

    # support2_raw = 2nd-largest of the clipped [5', bp, 3'] motif log-odds (the trained 6th feature).
    s2 = np.sort(np.clip(np.column_stack([r5, rbp, r3]), 0.0, None), axis=1)[:, 1]
    X = np.column_stack([r5, rbp, r3, np.nan_to_num(boff), np.nan_to_num(bconf), s2])

    # margin + P_motif only on rows with finite motif log-odds (unscorable rows stay NaN, keep alignment).
    finite = np.isfinite(r5) & np.isfinite(rbp) & np.isfinite(r3)
    margin = np.full(n, np.nan)
    if finite.any():
        margin[finite] = ensemble_margin(ensemble_models, X[finite])
    p_motif = p_motif_from_margin(margin, params)   # NaN where margin NaN

    result = adjudicate(margin[finite], p_motif[finite], params)
    q_eff, _, _ = _effective_q(result)

    p_adj = q_eff * p_motif
    adjusted = 100.0 * p_adj
    # REPORT-ONLY third axis (gates nothing): per-intron background surprisal for the call set.
    bg_fdr = background_fdr(margin, p_motif, params)

    def fmt(arr):
        return [("nan" if not np.isfinite(v) else f"{v:.6g}") for v in arr]

    df["P_motif"] = fmt(p_motif)
    # PRIMARY species gate (constant within a run/species): the motif-only population statistic + category.
    # DETECTED == 'a U12-motif population, by motif — corroborate downstream' (CAVEAT); not a confirmed bearer.
    df["z_excess"] = f"{result.z_excess:.6g}"
    _g = lambda v: ("nan" if not np.isfinite(v) else f"{v:.6g}")
    df["cs_p95"] = _g(result.cs_p95)         # null-free co-fallback gate value (POINT p95 of the call margins)
    df["p_gumbel_p95"] = _g(result.p_gumbel_p95)  # PRIMARY strength driver: per-genome Gumbel tail-prob at cs_p95
    df["cs_p95_lo"] = _g(result.cs_p95_lo)   # bootstrap CI lower bound — logged confidence annotation (not the gate)
    df["cs_p95_hi"] = _g(result.cs_p95_hi)
    df["motif_category"] = result.motif_category.value
    # REPORT-ONLY: background-relative surprisal of each call vs THIS genome's own U2 tail. NaN for
    # non-calls. Never feeds a call — the gate stays in type_id/adjusted_score (see background_fdr).
    df["bg_fdr"] = fmt(bg_fdr)
    # The GATED call pair. rel_score is adjusted_score recentred on the HC threshold; keeping the identity
    # exact is what makes `rel_score > 0` a self-sufficient one-column HC filter (see the docstring).
    df["adjusted_score"] = fmt(adjusted)
    df["rel_score"] = fmt(adjusted - 90.0)
    called = np.isfinite(p_adj) & (p_adj >= 0.5)
    df["type_id"] = np.where(called, "u12", "u2")
    df.to_csv(score_info_path, sep="\t", index=False)

    # Per-species tail-model diagnostic sidecar: the U2 background margins, the fitted exponential tail, and
    # the call margins — production-faithful (same margins the adjudicator scored). Consumed by
    # visualization.plots.tail_model_plot. Best-effort: a plot payload must never break classification output.
    try:
        import json as _json
        payload = build_tail_model(margin[finite], p_motif[finite], result, params)
        with open(tail_model_sidecar_path(score_info_path), "w") as _fh:
            _json.dump(payload, _fh)
    except Exception as _tm_err:   # noqa: BLE001 — diagnostic sidecar is strictly optional
        if messenger is not None:
            messenger.warning(f"tail-model sidecar not written: {_tm_err}")

    result.n_adj_called = int(called.sum())
    # UNGATED motif calls: P_motif>=0.5 DISREGARDING the species motif_category gate. Equals n_adj_called for
    # every non-NOT_DETECTED genome; for a NOT_DETECTED genome (q_eff=0 -> p_adj=0, type_id=u2 for all) this
    # recovers the count the gate suppressed (the raw P_motif survives only in score_info, not meta/type_id).
    result.n_motif_called = int((np.isfinite(p_motif) & (p_motif >= 0.5)).sum())
    if messenger is not None:
        messenger.info(
            f"pmotif_adjudicated: status={result.status.value} motif_category={result.motif_category.value} "
            f"z_excess={result.z_excess:.2f} cs_p95={result.cs_p95:.2f} (lo={result.cs_p95_lo:.2f}) "
            f"n_call={result.n_call} "
            f"(secondary: depth_tail={result.depth_tail:.2f} q={q_eff:.3f}) "
            f"-> {int(called.sum())} U12 calls / {int(finite.sum())} scorable introns")
    return result
