"""Mode-separation second-pass orchestration (intronIC v2.6).

The first classify pass (cluster-aware ensemble) writes a score_info.iic
with raw PWM scores and first-pass svm_score per intron. This module:

1. Reads the score_info.iic
2. Derives soft candidate weights from first-pass svm_score
3. Estimates per-feature (μ_U2, μ_U12) by weighted median
4. Evaluates the gate (n_eff floor + μ_U12 location prior + 3D Fisher-KDE
   valley depth from cluster_validation.compute_valley_depth)
5. If gate passes: re-z-scores under mode-separation, runs the second-pass
   ensemble on eligible introns (z_5p ≥ floor), and rewrites svm_score
6. Always emits a diagnostic JSON with provenance + quality tier

Designed as a self-contained post-classify step. The existing classify
path scores with the bundle's first_pass_ensemble; this module owns the
recalibration phase.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from intronIC.scoring.mode_separation import (
    ModeSeparationStats,
    GateDecision,
    candidate_weight_from_svm,
    fit_mode_separation,
    apply_mode_separation_z,
    compute_support2,
    evaluate_gate,
    DEFAULT_THRESHOLD,
    DEFAULT_N_FLOOR,
    DEFAULT_VALLEY_MIN,
    DEFAULT_Z5P_ELIGIBILITY,
    DEFAULT_CANDIDATE_CENTER,
    DEFAULT_CANDIDATE_STEEPNESS,
)


# Feature column order matches v5_modesep training: 6D feature vector.
_FEATURE_COLS = ("five_z", "bp_z", "three_z",
                 "bp_offset", "bp_scan_confidence", "support2")


@dataclass
class ModeSepResult:
    """Outcome + diagnostics from a mode-sep run on one score_info file."""
    route: str                      # "modesep" | "first_pass_fallback"
    gate_reason: str
    n_introns: int
    n_eligible: int                 # introns that passed z_5p eligibility
    n_called_u12: int
    n_eff_candidates: float
    valley_depth: Optional[float]
    mu_u2_5p: float
    mu_u12_5p: float
    mu_u12_5p_offset: float
    median_ensemble_sigma_called: Optional[float]
    p90_ensemble_sigma_called: Optional[float]
    quality_tier: str               # "A" | "B" | "C" | "F"
    first_pass_model_id: str
    second_pass_model_id: str


def _f1_weighted_mean(probs: np.ndarray, ensemble) -> np.ndarray:
    """F1-weighted mean across ensemble sub-models (matches predictor.py)."""
    f1_scores = np.array(
        [float(getattr(m, "f1_score", 1.0)) for m in ensemble.models]
    )
    weights = f1_scores / f1_scores.sum()
    return np.dot(weights, probs)


def _score_second_pass_ensemble(ensemble, X: np.ndarray, n_jobs: int = -1):
    """Run all sub-models in parallel, return (mean × 100, std × 100)."""
    if X.shape[0] == 0:
        return np.array([]), np.array([])
    submodels = [m.model for m in ensemble.models]
    probas = Parallel(n_jobs=n_jobs, prefer="threads")(
        delayed(lambda m: m.predict_proba(X)[:, 1])(m) for m in submodels
    )
    probas = np.asarray(probas)  # (n_models, n_introns)
    mean_p = _f1_weighted_mean(probas, ensemble)
    std_p = probas.std(axis=0)
    return mean_p * 100.0, std_p * 100.0


def _assign_quality_tier(result: ModeSepResult) -> str:
    """Quality tier rubric from MODESEP_INTEGRATION_PLAN.md §3.2."""
    if result.route != "modesep":
        return "F"
    valley = result.valley_depth if result.valley_depth is not None else 0.0
    n_eff = result.n_eff_candidates
    med_sig = result.median_ensemble_sigma_called
    sig_ok_strict = med_sig is not None and med_sig <= 10.0
    sig_ok_loose = med_sig is not None and med_sig <= 15.0

    if valley >= 0.5 and n_eff >= 20 and sig_ok_strict:
        return "A"
    if (valley >= 0.3 or n_eff >= 10) and sig_ok_loose:
        return "B"
    return "C"


def apply_mode_separation_postprocess(
    score_info_path: Path,
    runtime: dict,
    *,
    valley_depth_fn,
    threshold: float = DEFAULT_THRESHOLD,
    diagnostics_path: Optional[Path] = None,
    messenger=None,
) -> ModeSepResult:
    """Apply the mode-sep second pass to a score_info.iic file in-place.

    Parameters
    ----------
    score_info_path : Path
        The score_info.iic written by the existing classify pipeline.
        Required columns: name, svm_score, 5'_raw, bp_raw, 3'_raw,
        bp_offset, bp_scan_confidence.
    runtime : dict
        Result of `normalize_model_bundle()` on a v5 modesep bundle.
        Expected keys: normalizer_mode='modesep', ensemble (first-pass),
        second_pass_ensemble, modesep_params.
    valley_depth_fn : callable
        Pass `intronIC.scoring.cluster_validation.compute_valley_depth`.
        Injected so this module doesn't drag scipy into test contexts.
    threshold : float
        SVM score threshold for the "called U12" count (default 90).
    diagnostics_path : Optional[Path]
        If set, writes the diagnostic JSON to this path.
    messenger : optional
        intronIC Messenger for status logging; falls back to stderr prints.

    Returns
    -------
    ModeSepResult — outcome + diagnostics. Caller can attach to summary.
    """
    if runtime.get("normalizer_mode") != "modesep":
        raise ValueError(
            "apply_mode_separation_postprocess called with non-modesep bundle"
        )

    def _log(msg):
        if messenger is not None:
            messenger.info(msg)
        else:
            print(msg)

    params = runtime.get("modesep_params", {})
    n_floor = int(params.get("n_floor_candidates", DEFAULT_N_FLOOR))
    valley_min = float(params.get("valley_depth_min", DEFAULT_VALLEY_MIN))
    z5p_floor = float(params.get("z_floor_eligibility", DEFAULT_Z5P_ELIGIBILITY))
    cw_center = float(params.get("candidate_weight_center", DEFAULT_CANDIDATE_CENTER))
    cw_steep = float(params.get("candidate_weight_steepness", DEFAULT_CANDIDATE_STEEPNESS))
    anchor_5p = float(params.get("universal_anchors", {}).get("five_raw", 15.671))
    tolerance = float(params.get("mu_u12_5p_tolerance", 3.6))
    first_pass_id = str(params.get("first_pass_model_id", "unknown"))
    second_pass_id = str(runtime.get("_v3_bundle", {}).get("model_id", "v5_modesep"))

    _log(f"[modesep] reading {score_info_path}")
    df = pd.read_csv(
        score_info_path, sep="\t", low_memory=False,
        usecols=["name", "svm_score",
                 "5'_raw", "bp_raw", "3'_raw",
                 "bp_offset", "bp_scan_confidence"],
    )
    for c in ("svm_score", "5'_raw", "bp_raw", "3'_raw",
              "bp_offset", "bp_scan_confidence"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df_valid = df.dropna(subset=["5'_raw", "bp_raw", "3'_raw"])
    n_total = len(df_valid)

    fp_svm = df_valid["svm_score"].fillna(0).to_numpy()
    cw = candidate_weight_from_svm(fp_svm, center=cw_center, steepness=cw_steep)

    s5 = fit_mode_separation(df_valid["5'_raw"].to_numpy(), cw)
    sbp = fit_mode_separation(df_valid["bp_raw"].to_numpy(), cw)
    s3 = fit_mode_separation(df_valid["3'_raw"].to_numpy(), cw)

    # Build feature points for the Fisher-KDE valley check.
    # u12_candidate_points: weight ≥ 0.5; u2_candidate_points: weight < 0.5.
    raw_3d = np.column_stack([
        df_valid["5'_raw"].to_numpy(),
        df_valid["bp_raw"].to_numpy(),
        df_valid["3'_raw"].to_numpy(),
    ])
    cand_mask = cw >= 0.5
    u12_points = raw_3d[cand_mask]
    u2_points = raw_3d[~cand_mask]

    gate = evaluate_gate(
        u12_candidate_points=u12_points,
        u2_candidate_points=u2_points,
        mu_u12_5p_raw=s5.mu_u12,
        n_eff_candidates=s5.n_eff_candidates,
        separation_5p=s5.separation,
        valley_depth_fn=valley_depth_fn,
        n_floor=n_floor,
        valley_min=valley_min,
        mu_u12_prior=anchor_5p,
        mu_u12_tolerance=tolerance,
    )

    if not gate.passes:
        _log(f"[modesep] GATE-FAIL ({gate.reason}); keeping first-pass scores")
        result = ModeSepResult(
            route="first_pass_fallback",
            gate_reason=gate.reason,
            n_introns=n_total,
            n_eligible=0,
            n_called_u12=int((fp_svm >= threshold).sum()),
            n_eff_candidates=s5.n_eff_candidates,
            valley_depth=gate.valley_depth,
            mu_u2_5p=s5.mu_u2,
            mu_u12_5p=s5.mu_u12,
            mu_u12_5p_offset=(gate.mu_u12_offset
                              if gate.mu_u12_offset is not None
                              else float(s5.mu_u12 - anchor_5p)),
            median_ensemble_sigma_called=None,
            p90_ensemble_sigma_called=None,
            quality_tier="F",
            first_pass_model_id=first_pass_id,
            second_pass_model_id=second_pass_id,
        )
        _maybe_write_diagnostics(result, diagnostics_path, score_info_path)
        return result

    z5 = apply_mode_separation_z(df_valid["5'_raw"].to_numpy(), s5)
    zbp = apply_mode_separation_z(df_valid["bp_raw"].to_numpy(), sbp)
    z3 = apply_mode_separation_z(df_valid["3'_raw"].to_numpy(), s3)
    support2 = compute_support2(z5, zbp, z3)
    X_all = np.column_stack([
        z5, zbp, z3,
        df_valid["bp_offset"].fillna(0).to_numpy(),
        df_valid["bp_scan_confidence"].fillna(0).to_numpy(),
        support2,
    ])

    elig = z5 >= z5p_floor
    n_elig = int(elig.sum())
    _log(f"[modesep] GATE-PASS μ_U2={s5.mu_u2:.2f} μ_U12={s5.mu_u12:.2f} "
         f"valley={gate.valley_depth:.3f}; scoring {n_elig:,}/{n_total:,} "
         f"introns (z_5p ≥ {z5p_floor})")

    second_pass_ensemble = runtime["second_pass_ensemble"]
    elig_idx = np.where(elig)[0]
    if n_elig > 0:
        mean_scores, std_scores = _score_second_pass_ensemble(
            second_pass_ensemble, X_all[elig]
        )
    else:
        mean_scores = np.array([])
        std_scores = np.array([])

    # Rebuild svm_score: first-pass for ineligible, second-pass for eligible.
    final_svm = fp_svm.copy().astype(float)
    final_svm[elig] = mean_scores

    # Per-intron ensemble σ (0 for ineligible since first-pass scored).
    ensemble_sigma_all = np.zeros_like(fp_svm, dtype=float)
    ensemble_sigma_all[elig] = std_scores

    called = final_svm >= threshold
    n_called = int(called.sum())

    # σ on called introns only — borderline calls are where σ matters.
    called_eligible_mask = called & elig
    if called_eligible_mask.any():
        sigmas_called = ensemble_sigma_all[called_eligible_mask]
        med_sigma = float(np.median(sigmas_called))
        p90_sigma = float(np.quantile(sigmas_called, 0.90))
    else:
        med_sigma = None
        p90_sigma = None

    # Rewrite score_info.iic with the updated svm_score column.
    # Preserve the original file order by realigning via the 'name' key.
    df_full = pd.read_csv(score_info_path, sep="\t", low_memory=False)
    update_keys = df_valid["name"].to_numpy()
    update_map = dict(zip(update_keys, final_svm))
    df_full["svm_score"] = df_full["name"].map(update_map).fillna(df_full["svm_score"])

    # Add diagnostic columns
    sigma_map = dict(zip(update_keys, ensemble_sigma_all))
    df_full["ensemble_sigma"] = df_full["name"].map(sigma_map).fillna(0.0)
    fp_map = dict(zip(update_keys, fp_svm))
    df_full["first_pass_svm"] = df_full["name"].map(fp_map)
    df_full["modesep_route"] = "modesep"
    df_full.loc[~df_full["name"].isin(update_keys), "modesep_route"] = "untouched"

    df_full.to_csv(score_info_path, sep="\t", index=False, float_format="%.6f")
    _log(f"[modesep] rewrote {score_info_path} (n_called={n_called}, "
         f"median σ on called = {('NA' if med_sigma is None else f'{med_sigma:.3f}')})")

    result = ModeSepResult(
        route="modesep",
        gate_reason=gate.reason,
        n_introns=n_total,
        n_eligible=n_elig,
        n_called_u12=n_called,
        n_eff_candidates=s5.n_eff_candidates,
        valley_depth=gate.valley_depth,
        mu_u2_5p=s5.mu_u2,
        mu_u12_5p=s5.mu_u12,
        mu_u12_5p_offset=(gate.mu_u12_offset
                          if gate.mu_u12_offset is not None
                          else float(s5.mu_u12 - anchor_5p)),
        median_ensemble_sigma_called=med_sigma,
        p90_ensemble_sigma_called=p90_sigma,
        quality_tier="",  # filled below
        first_pass_model_id=first_pass_id,
        second_pass_model_id=second_pass_id,
    )
    result = ModeSepResult(**{**asdict(result),
                              "quality_tier": _assign_quality_tier(result)})
    _maybe_write_diagnostics(result, diagnostics_path, score_info_path)
    return result


def _maybe_write_diagnostics(result: ModeSepResult,
                              diagnostics_path: Optional[Path],
                              score_info_path: Path) -> None:
    if diagnostics_path is None:
        return
    payload = asdict(result)
    payload["score_info_path"] = str(score_info_path)
    Path(diagnostics_path).write_text(json.dumps(payload, indent=2))


__all__ = ["apply_mode_separation_postprocess", "ModeSepResult"]
