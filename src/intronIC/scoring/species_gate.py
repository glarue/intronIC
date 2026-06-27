"""Continuous species gate for the ``raw_gated`` scoring mode.

The raw_gated architecture classifies on the background-corrected raw log-odds with a
single global operating point, then applies ONE per-species output threshold that slides
continuously between strict (suppress loss-species false positives) and relaxed (recover the
divergent weak-branch-point U12 tail). It replaces the z-normalizer + mode-separation +
continuous-discount stack. See ``eval_corpus/PROPOSED_ARCHITECTURE.md`` for the rationale and
the panel validation (the algorithm and its calibration were developed in
``eval_corpus/raw_gated_v2.py``).

The gate is driven by TWO per-species signals computed from the raw-feature SVM output, combined
as a PRODUCT (both must fire -> robust to a few stray strong introns in a loss species):

  core = #( support2_raw >= s_core  AND  P(U12) >= p_core )
         "Is there a U12 CORE?" — count of all-ends-strong, high-confidence introns.
  gap  = median support2_raw of the top-K introns by P(U12)
         "Are the species' BEST candidates U12-shaped (all-ends-strong)?" — real bearers ~6-10,
         loss-species FPs (one-end-strong) ~0-2. The bimodality / separated-real-mode check.

  g          = ramp(core, core_lo, core_hi) * ramp(gap, gap_lo, gap_hi)        # in [0, 1]
  tau        = tau_high * (1 - g) + tau_low * g                                # strict .. relaxed
  call       = P(U12) >= tau

`g=0` (no U12 core) -> tau_high (strict) -> loss species suppressed; `g=1` (strong core) ->
tau_low (relaxed) -> the divergent tail is recovered. Because a loss species cannot inflate
motif-empty introns into a U12 mode, the gate needs no discount to undo.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np


@dataclass(frozen=True)
class SpeciesGateParams:
    """Calibrated parameters for the continuous species gate.

    Defaults are the conservative first-pass values validated on a 48-species panel
    (FP-control prioritized over recall): clean on all loss species, ~99.75% bearer
    conservation-U12 recall. They are intended to be re-calibrated on a broader panel.
    """
    s_core: float = 5.0       #: support2_raw threshold defining an "all-ends-strong" intron
    p_core: float = 0.90      #: P(U12) threshold defining a "high-confidence" core intron
    topk: int = 100           #: number of top-P introns used for the `gap` signal
    core_lo: float = 4.0      #: core-count ramp low edge (g_core = 0 at/below)
    core_hi: float = 50.0     #: core-count ramp high edge (g_core = 1 at/above)
    gap_lo: float = 2.5       #: gap ramp low edge (g_gap = 0 at/below)
    gap_hi: float = 6.0       #: gap ramp high edge (g_gap = 1 at/above)
    tau_high: float = 0.985   #: strict per-species threshold (no U12 core)
    tau_low: float = 0.85     #: relaxed per-species threshold (strong U12 core)

    @classmethod
    def from_dict(cls, d: Optional[dict]) -> "SpeciesGateParams":
        if not d:
            return cls()
        known = {f for f in cls.__dataclass_fields__}
        return cls(**{k: v for k, v in d.items() if k in known})


@dataclass
class SpeciesGateResult:
    """Outcome of applying the gate to one species' introns."""
    core: int                 #: number of all-ends-strong high-confidence introns
    gap: float                #: median support2_raw of the top-K introns by P
    g: float                  #: combined gate signal in [0, 1]
    tau: float                #: per-species P(U12) threshold actually used
    calls: np.ndarray         #: boolean array, True where P(U12) >= tau (NaN P -> False)
    n_called: int             #: number of U12 calls
    params: SpeciesGateParams = field(repr=False, default_factory=SpeciesGateParams)


def _ramp(x: float, lo: float, hi: float) -> float:
    """Linear ramp clamped to [0, 1]: 0 at/below `lo`, 1 at/above `hi`."""
    if hi <= lo:
        return 1.0 if x >= hi else 0.0
    return float(min(1.0, max(0.0, (x - lo) / (hi - lo))))


def compute_gate_signals(probs: Sequence[float], support2_raw: Sequence[float],
                         params: SpeciesGateParams) -> tuple:
    """Compute (core, gap, g) for one species from its per-intron P(U12) and support2_raw.

    NaN probabilities (unscorable introns) are ignored. ``core`` is log-ramped (a U12-rich
    genome has hundreds of core introns; the count is genome-size dependent so the ramp is on
    log10), ``gap`` is linearly ramped.
    """
    p = np.asarray(probs, dtype=float)
    s2 = np.asarray(support2_raw, dtype=float)
    finite = ~np.isnan(p)
    core = int(np.sum(finite & (s2 >= params.s_core) & (p >= params.p_core)))
    if finite.any():
        k = min(params.topk, int(finite.sum()))
        # top-k indices by P among finite-P introns
        idx = np.argsort(np.where(finite, p, -np.inf))[-k:]
        gvals = s2[idx]
        gap = float(np.nanmedian(gvals)) if np.isfinite(gvals).any() else 0.0
    else:
        gap = 0.0
    g_core = _ramp(math.log10(core + 1), math.log10(params.core_lo), math.log10(params.core_hi))
    g_gap = _ramp(gap, params.gap_lo, params.gap_hi)
    g = g_core * g_gap
    return core, gap, g


def apply_species_gate(probs: Sequence[float], support2_raw: Sequence[float],
                       params: Optional[SpeciesGateParams] = None) -> SpeciesGateResult:
    """Apply the continuous species gate to one species' introns.

    Args:
        probs: per-intron P(U12) from the raw-feature classifier (NaN allowed for unscorable).
        support2_raw: per-intron :pyattr:`IntronScores.support2_raw` (NaN allowed).
        params: gate parameters (defaults to the calibrated :class:`SpeciesGateParams`).

    Returns:
        :class:`SpeciesGateResult` with the gate signals and the per-intron U12 calls.
    """
    params = params or SpeciesGateParams()
    p = np.asarray(probs, dtype=float)
    s2 = np.asarray(support2_raw, dtype=float)
    if p.shape != s2.shape:
        raise ValueError(f"probs and support2_raw shape mismatch: {p.shape} vs {s2.shape}")
    core, gap, g = compute_gate_signals(p, s2, params)
    tau = params.tau_high * (1.0 - g) + params.tau_low * g
    calls = np.where(np.isnan(p), False, p >= tau)
    return SpeciesGateResult(core=core, gap=gap, g=g, tau=tau, calls=calls,
                             n_called=int(calls.sum()), params=params)


def apply_raw_gated_postprocess(score_info_path, params: Optional[SpeciesGateParams] = None,
                                threshold: float = 90.0, messenger=None) -> SpeciesGateResult:
    """raw_gated post-classification: read score_info.iic, apply the species gate, write final calls.

    Replaces mode-separation + the continuous discount for ``scoring_mode == raw_gated`` bundles. Reads the
    raw motif log-odds (``5'_raw``/``bp_raw``/``3'_raw``) and the first-pass probability (``svm_score`` on
    a 0-100 scale), computes ``support2_raw`` and the per-species gate, and rewrites ``type_id`` (the final
    U12 call), ``adjusted_score`` (kept = the raw-model probability) and ``rel_score`` (``svm_score`` minus
    the per-species call threshold, so ``rel_score > 0`` iff called). Per-run = per-species (one genome).
    Returns the :class:`SpeciesGateResult` (gate diagnostics + calls).
    """
    import pandas as pd
    params = params or SpeciesGateParams()
    df = pd.read_csv(score_info_path, sep="\t", dtype=str, keep_default_na=False)

    def col(name):
        return df[name].replace("", "nan").astype(float).to_numpy() if name in df.columns else None
    r5, rbp, r3 = col("5'_raw"), col("bp_raw"), col("3'_raw")
    svm = col("svm_score")
    if r5 is None or rbp is None or r3 is None or svm is None:
        raise ValueError("raw_gated post-process: score_info missing 5'_raw/bp_raw/3'_raw/svm_score")
    s2 = np.sort(np.clip(np.column_stack([r5, rbp, r3]), 0.0, None), axis=1)[:, 1]
    P = svm / 100.0
    res = apply_species_gate(P, s2, params)

    tau_pct = 100.0 * res.tau
    df["type_id"] = np.where(res.calls, "u12", "u2")
    df["adjusted_score"] = [f"{v:.6g}" for v in svm]                 # keep the raw-model probability
    df["rel_score"] = [f"{v - tau_pct:.6g}" for v in svm]            # > 0 iff called
    df.to_csv(score_info_path, sep="\t", index=False)

    if messenger is not None:
        messenger.info(
            f"raw_gated gate: core={res.core} gap={res.gap:.1f} g={res.g:.2f} "
            f"tau={res.tau:.3f} -> {res.n_called} U12 calls / {len(P)} introns")
    return res
