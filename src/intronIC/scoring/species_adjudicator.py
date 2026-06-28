"""Species-level U12 adjudicator: decide whether a genome's motif-strong introns form a real, separated
U12 *population* (bearer) or are the upper tail of the unimodal U2 background (loss / recent-loss), and
turn that into a per-genome bearer probability ``q`` with a bootstrap CI.

Pipeline role (see ``docs/raw_gated_scoring.md`` §0b/§0c): the SVM gives a per-intron, species-agnostic
motif probability ``P_motif`` (calibrated ensemble margin). This module answers the per-*species* question
and produces the per-intron posterior ``P_adj = q * P_motif`` (functional-U12 posterior; clean because
``P(functional U12 | motif, not-a-bearer) ≈ 0``). The low-numbers uncertainty flows through ``q``'s CI.

Method (size-INVARIANT depth-beyond-U2-tail test on the margin):
  - U2 reference = introns with ``P_motif < u2_threshold``; calls = ``P_motif >= call_threshold``.
  - robust call-core = ``median(margin_calls)`` (median peels a fringe of edge-of-U2 mis-calls).
  - **PRIMARY q-driver = ``depth_tail`` = (call-core − U2's own 99.9th pct) / MAD_U2.** This is *size-
    invariant*: it asks "does a real U12 *population* sit in territory the U2 bulk does not reach?"
    (``Q_pop``), independent of genome size. ``q = sigma(q_a*depth_tail + q_b)``.
  - The earlier ``excess_z`` (call-core beyond the U2 EXPECTED MAXIMUM, which grows with the U2 count via a
    ``ln(N_u2)`` term) is a *size-aware significance* test (``Q_sig``); it is retained here as a **labelled
    SECONDARY** diagnostic, NOT the q-driver. See the 2026-06-27 committee/meta-review record in
    ``docs/raw_gated_scoring.md`` §0d for why the swap (the "#2 N-dependence" fix).

Calibration provenance: ``q``'s ``(q_a, q_b)`` are a **separation-safe Firth-penalized** logistic fit of
``depth_tail`` against bearer/loss truth on the 37-genome panel (``eval_corpus`` snRNA + IPA labels;
Symbiodinium recoded as a CONFLICT case and excluded from the loss class, not anchored as a hard loss). The
panel is cleanly separated so a plain MLE diverges; Firth (Jeffreys-prior penalty) is the standard fix.
**Leave-clade-out validated** (``eval_corpus/q_firth_leaveclade.py``): 26/27 genomes classify correctly
out-of-sample; the lone OOF miss is chlamydomonas (the deepest loss, 3 calls) tipping to q≈0.60 only when its
whole Chlorophyta clade is held out — and with 3 calls its bootstrap CI is wide, so it surfaces as
``UNDETERMINED`` (the undetermined band), not a confident bearer. The constants are version-pinned by
``ADJUDICATOR_PARAMS_VERSION``. (ship-blocker #2 of the 2026-06-27 committee verdict, §0d.)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

import numpy as np

#: Version pin for the calibration constants (bundle-stamped in production). Bump when (q_a, q_b) or the
#: depth_tail/EVT recipe changes so a stale bundle can be detected.
ADJUDICATOR_PARAMS_VERSION = "depth_tail_firth_2026-06-27"

# Calibration constants frozen from the 37-genome panel (eval_corpus/{robust_sep_one,q_firth_leaveclade}.py).
# These are bundle-stamped in production; the defaults document the prototype values.
DEFAULT_PLATT_A = 2.796       #: P_motif = sigma(PLATT_A * margin + PLATT_C)
DEFAULT_PLATT_C = -1.178
DEFAULT_Q_A = 3.64            #: q = sigma(Q_A * depth_tail + Q_B)   (Firth fit; q=0.5 at -Q_B/Q_A = +2.98, in the gap)
DEFAULT_Q_B = -10.86


def _sigmoid(z):
    return 1.0 / (1.0 + np.exp(-z))


class AdjStatus(str, Enum):
    """Operational status of one genome's adjudication (every result carries exactly one)."""
    ADJUDICATED = "ADJUDICATED"          #: full assessment succeeded (q + CI valid)
    UNDETERMINED = "UNDETERMINED"        #: q's 95% CI straddles the decision boundary -> bearer/loss unresolved
    LOW_N = "LOW_N"                      #: too few calls or species U2 to assess -> file side falls to q_eff=1
    DEGENERATE_TAIL = "DEGENERATE_TAIL"  #: U2 scale (MAD) ~0 / non-finite -> cannot standardize -> unassessable
    SCHEMA_FAIL = "SCHEMA_FAIL"          #: required inputs missing / non-finite / shape mismatch


@dataclass(frozen=True)
class AdjudicatorParams:
    """Calibrated constants + method settings for the species adjudicator (bundle-stamped in production)."""
    platt_a: float = DEFAULT_PLATT_A    #: margin -> P_motif Platt slope
    platt_c: float = DEFAULT_PLATT_C    #: margin -> P_motif Platt intercept
    q_a: float = DEFAULT_Q_A            #: depth_tail -> q logistic slope
    q_b: float = DEFAULT_Q_B            #: depth_tail -> q logistic intercept
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
    depth_tail: float                    #: PRIMARY: call-core depth beyond U2's 99.9th pct, in MAD units (>~3 => bearer)
    q: float                             #: P(genome is a U12 bearer), point estimate
    q_lo: float                          #: q lower CI (bootstrap)
    q_hi: float                          #: q upper CI (bootstrap)
    status: AdjStatus = AdjStatus.ADJUDICATED
    excess_z: float = float("nan")       #: SECONDARY (size-aware significance): call-core beyond U2 expected-max
    p_gumbel: float = float("nan")       #: SECONDARY: Gumbel tail-prob that U2's own max reaches the call-core
    z_expmax: float = float("nan")       #: SECONDARY: U2 expected-max in robust-z units (the size-aware cutoff)
    secondary_available: bool = False    #: True when the EVT secondary (excess_z/p_gumbel) could be computed
    n_adj_called: Optional[int] = None   #: final P_adj>=0.5 U12 calls (set by the file-side application; None otherwise)
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
    q_tail: float          #: the ``tail_pct`` (99.9th) U2 percentile — the primary depth reference
    z_expmax: float        #: U2 expected-max in robust-z units (secondary); NaN if EVT tail unfit
    evt_beta: float        #: Gumbel scale in robust-z units (secondary); NaN if EVT tail unfit
    finite: bool           #: True when (med, mad, q_tail) are usable (mad > 0 and finite)


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

    # ---- SECONDARY: exponential-tail EVT expected-max (size-aware; may be unavailable) ----
    q_evt = float(np.percentile(u2_margin, params.evt_tail_pct))
    exc = u2_margin[u2_margin > q_evt] - q_evt
    z_expmax = np.nan
    evt_beta = np.nan
    if len(exc) >= params.min_evt_excesses and exc.mean() > 0:
        lam = 1.0 / float(exc.mean())
        exp_max_margin = q_evt + np.log(max(params.evt_tail_frac * n_u2, 2.0)) / lam
        z_expmax = (exp_max_margin - med) / mad
        evt_beta = 1.0 / (lam * mad + 1e-9)
    return _U2Reference(med, mad, q_tail, z_expmax, evt_beta, finite=True)


def compute_depth_tail(call_core_margin: float, ref: _U2Reference) -> float:
    """PRIMARY size-invariant separation: how far the call-core sits beyond U2's 99.9th pct, in MAD units."""
    return float((call_core_margin - ref.q_tail) / ref.mad)


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


# --------------------------------------------------------------------------------------------------
# The adjudicator
# --------------------------------------------------------------------------------------------------
def _fail(n_call, n_u2, status, params) -> AdjudicatorResult:
    """Construct a non-assessable result (q defaults to the loss side; CI collapsed)."""
    return AdjudicatorResult(n_call=int(n_call), n_u2=int(n_u2), depth_tail=float("nan"),
                             q=0.0, q_lo=0.0, q_hi=0.0, status=status,
                             secondary_available=False, params=params)


def _assess(calls: np.ndarray, n_u2: int, ref: _U2Reference,
            params: AdjudicatorParams) -> AdjudicatorResult:
    """Score a call set against the species U2 reference: depth_tail -> q + bootstrap CI + status
    (ADJUDICATED if the bootstrap q-CI clears the 0.5 boundary, else UNDETERMINED). The secondary EVT
    diagnostic (excess_z/p_gumbel) is attached when the U2 tail could be fit."""
    call_core = float(np.median(calls))
    depth_tail = compute_depth_tail(call_core, ref)
    if not np.isfinite(depth_tail):
        return _fail(len(calls), n_u2, AdjStatus.DEGENERATE_TAIL, params)
    q_point = float(q_from_depth_tail(depth_tail, params))

    excess_z = compute_excess_z(call_core, ref)
    p_gumbel = compute_p_gumbel(call_core, ref)
    secondary_available = bool(np.isfinite(excess_z))

    # COUNT-AWARE SMOOTHED bootstrap of the call-core -> q CI. A plain bootstrap-of-the-median is degenerate
    # for tied/clustered calls (the resampled median takes few distinct values -> a width-0 CI that silently
    # defeats the UNDETERMINED band in exactly the low-call-count regime it exists for). Smoothing each
    # resampled call by a kernel (bandwidth = the robust call spread, floored to a fraction of the U2 scale)
    # (a) breaks the tie degeneracy and (b) makes the CI widen as the call count shrinks (the median's
    # sampling SD ~ bandwidth/sqrt(k)), so few/borderline calls correctly read UNDETERMINED. The point q
    # still uses the un-smoothed call-core, so this only affects the CI (and stays deterministic).
    rng = np.random.RandomState(params.random_seed)
    k = len(calls)
    mad_calls = 1.4826 * float(np.median(np.abs(calls - call_core)))
    bw = max(mad_calls, params.ci_smooth_floor_frac * ref.mad)
    qs = np.empty(params.bootstrap_n)
    for b in range(params.bootstrap_n):
        samp = calls[rng.randint(0, k, k)] + rng.normal(0.0, bw, k)
        core_b = float(np.median(samp))
        qs[b] = q_from_depth_tail(compute_depth_tail(core_b, ref), params)
    q_lo, q_hi = (float(x) for x in np.percentile(qs, params.ci))

    status = AdjStatus.ADJUDICATED if (q_lo > 0.5 or q_hi < 0.5) else AdjStatus.UNDETERMINED
    return AdjudicatorResult(
        n_call=int(len(calls)), n_u2=int(n_u2), depth_tail=float(depth_tail),
        q=q_point, q_lo=q_lo, q_hi=q_hi, status=status,
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

    # ---- LOW_N: too few calls or species U2 to assess a population ----
    if len(calls) < params.min_calls or len(u2) < params.min_u2:
        return _fail(len(calls), len(u2), AdjStatus.LOW_N, params)

    n_u2 = len(u2) / n_total * n_total   # full-population U2 count (== len(u2) unless subsampled upstream)
    ref = _u2_reference(u2, n_u2, params)

    # ---- DEGENERATE_TAIL: species U2 scale collapsed -> cannot standardize ----
    if not ref.finite:
        return _fail(len(calls), len(u2), AdjStatus.DEGENERATE_TAIL, params)

    return _assess(calls, len(u2), ref, params)


# --------------------------------------------------------------------------------------------------
# File-side post-process (the inference hook; mirrors species_gate.apply_raw_gated_postprocess)
# --------------------------------------------------------------------------------------------------
#: Raw-feature column order the ensemble was trained on (see eval_corpus/robust_sep_one.py). The 6th
#: feature, support2_raw, is the 2nd-largest of the clipped motif log-odds and is derived here.
_RAW_FEATURE_COLS = ("5'_raw", "bp_raw", "3'_raw", "bp_offset", "bp_scan_confidence")


def _effective_q(result: "AdjudicatorResult"):
    """Map an adjudication outcome to the per-intron multiplier ``q_eff`` (+ CI) applied to ``P_motif``.

    - ADJUDICATED / UNDETERMINED -> use the assessed ``q`` and its bootstrap CI.
    - LOW_N / DEGENERATE_TAIL / SCHEMA_FAIL -> the species layer is *inconclusive* (too few species U2, a
      malformed input, or no motif calls), so it must NOT suppress: default to ``q_eff = 1`` (P_adj =
      P_motif, the species-agnostic motif call — option A, the safe low-N default). A confirmed loss needs a
      successful adjudication to suppress; absent that ``P_motif`` is the call (matches the design's "P_motif
      everywhere; suppress only confirmed-depleted"). Users who want low-N suppression lower ``min_u2`` so a
      small genome reaches the ADJUDICATED path on its own U2 tail (option B).
    """
    if result.status in (AdjStatus.ADJUDICATED, AdjStatus.UNDETERMINED):
        return result.q, result.q_lo, result.q_hi
    return 1.0, 1.0, 1.0


def apply_pmotif_adjudication(score_info_path, ensemble_models,
                              params: Optional[AdjudicatorParams] = None,
                              messenger=None) -> "AdjudicatorResult":
    """``pmotif_adjudicated`` post-classification: read score_info.iic, compute the per-intron motif
    probability ``P_motif`` from the ensemble MARGIN, run the per-species adjudicator, and write the two
    interpretable per-intron numbers + the call.

    Adds/overwrites columns in ``score_info.iic``:
      - ``P_motif``  : sigmoid(Platt(margin)) — species-agnostic sequence-motif probability;
      - ``q``        : P(genome is a U12 bearer) (constant within a run / species);
      - ``P_adj``    : ``q_eff * P_motif`` — functional-U12 posterior, with ``P_adj_lo``/``P_adj_hi`` (CI);
      - ``adjusted_score`` = ``100*P_adj`` and ``rel_score`` = ``100*P_adj - 90`` (existing conventions:
        ``type_id == u12`` iff ``adjusted_score >= 50``; ``rel_score > 0`` iff high-confidence ``P_adj>=0.9``);
      - ``type_id``  : ``u12`` iff ``P_adj >= 0.5`` else ``u2``.

    Unscorable introns (NaN motif log-odds) keep NaN ``P_motif``/``P_adj`` and are not called. Per-run =
    per-species (one genome). Returns the :class:`AdjudicatorResult` (q + CI + status diagnostics).
    """
    import pandas as pd
    params = params or AdjudicatorParams()
    df = pd.read_csv(score_info_path, sep="\t", dtype=str, keep_default_na=False)

    def col(name):
        if name not in df.columns:
            return None
        return df[name].replace("", "nan").astype(float).to_numpy()

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
    q_eff, q_lo, q_hi = _effective_q(result)

    p_adj = q_eff * p_motif
    p_adj_lo = q_lo * p_motif
    p_adj_hi = q_hi * p_motif
    adjusted = 100.0 * p_adj

    def fmt(arr):
        return [("nan" if not np.isfinite(v) else f"{v:.6g}") for v in arr]

    df["P_motif"] = fmt(p_motif)
    df["q"] = f"{q_eff:.6g}"
    df["P_adj"] = fmt(p_adj)
    df["P_adj_lo"] = fmt(p_adj_lo)
    df["P_adj_hi"] = fmt(p_adj_hi)
    df["adjusted_score"] = fmt(adjusted)
    df["rel_score"] = fmt(adjusted - 90.0)
    called = np.isfinite(p_adj) & (p_adj >= 0.5)
    df["type_id"] = np.where(called, "u12", "u2")
    df.to_csv(score_info_path, sep="\t", index=False)

    result.n_adj_called = int(called.sum())
    if messenger is not None:
        messenger.info(
            f"pmotif_adjudicated: status={result.status.value} n_call={result.n_call} "
            f"depth_tail={result.depth_tail:.2f} q={q_eff:.3f} CI=[{q_lo:.2f},{q_hi:.2f}] "
            f"-> {int(called.sum())} U12 calls / {int(finite.sum())} scorable introns")
    return result
