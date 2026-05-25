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
    apply_continuous_discount,
    compute_boundary_mass,
    voting_fraction,
    DEFAULT_THRESHOLD,
    DEFAULT_N_FLOOR,
    DEFAULT_VALLEY_MIN,
    DEFAULT_Z5P_ELIGIBILITY,
    DEFAULT_CANDIDATE_CENTER,
    DEFAULT_CANDIDATE_STEEPNESS,
    DEFAULT_DISCOUNT_K_OVERCALL,
    DEFAULT_DISCOUNT_TAU_OVERCALL,
    DEFAULT_DISCOUNT_K_WEAKMOT,
    DEFAULT_DISCOUNT_TAU_MOTIF,
    DEFAULT_BM_SUBSAMPLE_SIZE,
    DEFAULT_BM_LOWER,
    DEFAULT_BM_UPPER,
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
    quality_tier: str               # "modesep_strong" | "modesep_standard" | "modesep_weak" | "first_pass_fallback"
    first_pass_model_id: str
    second_pass_model_id: str
    # v2.7+ additions (diagnostic-only species-level fields)
    boundary_mass: Optional[float] = None  # frac eligible w/ second-pass P in [0.1, 0.9]
    n_called_pre_discount: Optional[int] = None  # before continuous discount
    n_called_post_discount: Optional[int] = None  # after continuous discount
    continuous_discount_applied: bool = False


def _f1_weighted_mean(probs: np.ndarray, ensemble) -> np.ndarray:
    """F1-weighted mean across ensemble sub-models (matches predictor.py)."""
    f1_scores = np.array(
        [float(getattr(m, "f1_score", 1.0)) for m in ensemble.models]
    )
    weights = f1_scores / f1_scores.sum()
    return np.dot(weights, probs)


def _score_second_pass_ensemble(ensemble, X: np.ndarray, n_jobs: int = -1,
                                  return_probas: bool = False):
    """Run all sub-models in parallel.

    Default returns (mean × 100, std × 100). If `return_probas=True`,
    additionally returns the (n_models, n_introns) per-model probability
    matrix — used for downstream voting fraction and boundary-mass
    diagnostics.
    """
    if X.shape[0] == 0:
        if return_probas:
            return np.array([]), np.array([]), np.zeros((0, 0))
        return np.array([]), np.array([])
    submodels = [m.model for m in ensemble.models]
    probas = Parallel(n_jobs=n_jobs, prefer="threads")(
        delayed(lambda m: m.predict_proba(X)[:, 1])(m) for m in submodels
    )
    probas = np.asarray(probas)  # (n_models, n_introns)
    mean_p = _f1_weighted_mean(probas, ensemble)
    std_p = probas.std(axis=0)
    if return_probas:
        return mean_p * 100.0, std_p * 100.0, probas
    return mean_p * 100.0, std_p * 100.0


def _assign_quality_tier(result: ModeSepResult) -> str:
    """Quality tier rubric from MODESEP_INTEGRATION_PLAN.md §3.2.

    Tier strings (v2.7.1+): self-descriptive rename from the prior A/B/C/F
    school-grade codes. Mapping:
        modesep_strong       (was "A")
        modesep_standard     (was "B")
        modesep_weak         (was "C")
        first_pass_fallback  (was "F")
    """
    if result.route != "modesep":
        return "first_pass_fallback"
    valley = result.valley_depth if result.valley_depth is not None else 0.0
    n_eff = result.n_eff_candidates
    med_sig = result.median_ensemble_sigma_called
    sig_ok_strict = med_sig is not None and med_sig <= 10.0
    sig_ok_loose = med_sig is not None and med_sig <= 15.0

    if valley >= 0.5 and n_eff >= 20 and sig_ok_strict:
        return "modesep_strong"
    if (valley >= 0.3 or n_eff >= 10) and sig_ok_loose:
        return "modesep_standard"
    return "modesep_weak"


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
            quality_tier="first_pass_fallback",
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
        mean_scores, std_scores, probas_elig = _score_second_pass_ensemble(
            second_pass_ensemble, X_all[elig], return_probas=True
        )
    else:
        mean_scores = np.array([])
        std_scores = np.array([])
        probas_elig = np.zeros((0, 0))

    # Rebuild svm_score: first-pass for ineligible, second-pass for eligible.
    final_svm = fp_svm.copy().astype(float)
    final_svm[elig] = mean_scores

    # Per-intron ensemble σ (0 for ineligible since first-pass scored).
    ensemble_sigma_all = np.zeros_like(fp_svm, dtype=float)
    ensemble_sigma_all[elig] = std_scores

    # Per-intron voting fraction (v2.7 diagnostic): per-model votes among
    # second-pass models. Ineligible introns: NaN (no second-pass scoring).
    voting_all = np.full(len(fp_svm), np.nan, dtype=float)
    if probas_elig.size > 0:
        voting_all[elig] = voting_fraction(probas_elig)

    # Species-level boundary mass (v2.7 diagnostic): fraction of eligible
    # introns with second-pass mean P in [0.1, 0.9]. Subsample if large.
    bm_val: Optional[float] = None
    if probas_elig.size > 0:
        mean_p_elig = probas_elig.mean(axis=0)
        if len(mean_p_elig) > DEFAULT_BM_SUBSAMPLE_SIZE:
            rng = np.random.default_rng(0)
            sample = rng.choice(len(mean_p_elig),
                                size=DEFAULT_BM_SUBSAMPLE_SIZE, replace=False)
            bm_val = compute_boundary_mass(mean_p_elig[sample],
                                            lower=DEFAULT_BM_LOWER,
                                            upper=DEFAULT_BM_UPPER)
        else:
            bm_val = compute_boundary_mass(mean_p_elig,
                                            lower=DEFAULT_BM_LOWER,
                                            upper=DEFAULT_BM_UPPER)

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

    # v2.7 diagnostic columns: raw_sum, svm_vs_naive, voting_frac
    raw_sum_all = (df_valid["5'_raw"].to_numpy()
                   + df_valid["bp_raw"].to_numpy()
                   + df_valid["3'_raw"].to_numpy())
    eps = 1e-9
    p_clip = np.clip(final_svm / 100.0, eps, 1 - eps)
    logit_svm = np.log(p_clip / (1 - p_clip))
    svm_vs_naive_all = logit_svm - raw_sum_all

    raw_sum_map = dict(zip(update_keys, raw_sum_all))
    svn_map = dict(zip(update_keys, svm_vs_naive_all))
    voting_map = dict(zip(update_keys, voting_all))
    df_full["raw_sum"] = df_full["name"].map(raw_sum_map)
    df_full["svm_vs_naive"] = df_full["name"].map(svn_map)
    df_full["voting_frac"] = df_full["name"].map(voting_map)

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
        boundary_mass=bm_val,
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


def apply_continuous_per_intron_discount(
    score_info_path: Path,
    *,
    threshold: float = DEFAULT_THRESHOLD,
    k_overcall: float = DEFAULT_DISCOUNT_K_OVERCALL,
    tau_overcall: float = DEFAULT_DISCOUNT_TAU_OVERCALL,
    k_weakmot: float = DEFAULT_DISCOUNT_K_WEAKMOT,
    tau_motif: float = DEFAULT_DISCOUNT_TAU_MOTIF,
    input_column: str = "svm_score",
    messenger=None,
) -> dict:
    """Apply continuous per-intron discount to a score_info.iic file (v2.7).

    Reads `input_column` (default `svm_score`) and writes the result to
    `adjusted_score`:

        penalty_oc = k_overcall * max(0, svm_vs_naive - tau_overcall)
        penalty_wm = k_weakmot  * max(0, tau_motif - raw_sum)
        logit_adj  = logit(p_input) - penalty_oc - penalty_wm
        adjusted_score = sigmoid(logit_adj) * 100

    where svm_vs_naive = logit(p_input) - raw_sum and
          raw_sum     = 5'_raw + bp_raw + 3'_raw.

    `svm_score` is preserved. `adjusted_score` is overwritten with the
    discount applied.

    To chain discounts (e.g., legacy valley-depth discount on gate-fail
    species followed by this continuous discount), call this function with
    `input_column="adjusted_score"` after the legacy step. Otherwise the
    default `input_column="svm_score"` applies the discount directly to
    the raw / mode-sep-recalibrated SVM output.

    Returns a dict with summary stats (n_called pre/post, params used).
    """
    def _log(msg):
        if messenger is not None:
            messenger.info(msg)

    df = pd.read_csv(score_info_path, sep="\t", low_memory=False)
    for c in ("svm_score", "5'_raw", "bp_raw", "3'_raw", "adjusted_score"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    if input_column not in df.columns:
        raise ValueError(
            f"input_column '{input_column}' not in score_info; "
            f"available: {list(df.columns)}"
        )

    base = df[input_column].fillna(df.get("svm_score", 0)).to_numpy()

    valid = ~(df["5'_raw"].isna() | df["bp_raw"].isna() | df["3'_raw"].isna())
    raw_sum = np.where(
        valid,
        df["5'_raw"].fillna(0).to_numpy()
        + df["bp_raw"].fillna(0).to_numpy()
        + df["3'_raw"].fillna(0).to_numpy(),
        np.nan,
    )

    n_pre = int((base >= threshold).sum())
    # Apply discount only on valid rows; ineligible rows pass through
    adj = base.copy().astype(float)
    valid_idx = np.where(valid)[0]
    if len(valid_idx) > 0:
        adj_valid = apply_continuous_discount(
            base[valid_idx], raw_sum[valid_idx],
            k_overcall=k_overcall, tau_overcall=tau_overcall,
            k_weakmot=k_weakmot, tau_motif=tau_motif,
        )
        adj[valid_idx] = adj_valid
    n_post = int((adj >= threshold).sum())

    df["adjusted_score"] = adj

    # v2.7.1: unified labels (type_id / confidence / history)
    # See scoring/labeling.py for the rubric.
    from intronIC.scoring.labeling import (
        assign_labels, count_by_label, U2_STRONG_THRESHOLD,
        PROMOTED_DEMOTED_THRESHOLD,
    )
    # For Tier F species the first_pass_svm column is absent; fall back to
    # the input column (svm_score = first-pass output in that regime).
    if "first_pass_svm" in df.columns:
        adj_series = pd.Series(adj, index=df.index)
        fp_col = pd.to_numeric(df["first_pass_svm"], errors="coerce").fillna(adj_series)
    else:
        fp_col = pd.Series(base, index=df.index)
    fp_arr = fp_col.to_numpy()

    labels = [assign_labels(float(fp), float(a)) for fp, a in zip(fp_arr, adj)]
    df["type_id"] = [lab.type_id for lab in labels]
    df["confidence"] = [lab.confidence for lab in labels]
    df["history"] = [lab.history for lab in labels]

    df.to_csv(score_info_path, sep="\t", index=False, float_format="%.6f")
    _log(f"[continuous-discount] applied; n_called_pre={n_pre}, "
         f"n_called_post={n_post}, suppressed={n_pre - n_post}")

    label_counts = count_by_label(labels)
    _log(f"[labeling] u12={label_counts['u12_count']} "
         f"(strong={label_counts['u12_strong_count']}, "
         f"borderline={label_counts['u12_borderline_count']}, "
         f"promoted={label_counts['u12_promoted_count']}); "
         f"u2={label_counts['u2_count']} "
         f"(strong={label_counts['u2_strong_count']}, "
         f"borderline={label_counts['u2_borderline_count']}, "
         f"demoted={label_counts['u2_demoted_count']})")

    return {
        "n_called_pre_discount": n_pre,
        "n_called_post_discount": n_post,
        "n_suppressed": n_pre - n_post,
        "k_overcall": k_overcall,
        "tau_overcall": tau_overcall,
        "k_weakmot": k_weakmot,
        "tau_motif": tau_motif,
        **label_counts,
    }


__all__ = [
    "apply_mode_separation_postprocess",
    "apply_continuous_per_intron_discount",
    "ModeSepResult",
]
