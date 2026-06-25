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
from dataclasses import dataclass, asdict, replace
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
    species_penalty_logit_shift,
    DEFAULT_SPECIES_PENALTY_PGATE,
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
    # v3 gate-gapfrac (Step 3): separation statistic carried from the gate. gap_fraction_ucl
    # (the bootstrap upper-confidence-limit) drives the modesep-route π_species prior; both
    # are surfaced as diagnostics.
    gap_fraction: Optional[float] = None
    gap_fraction_ucl: Optional[float] = None
    centroid_sigma: Optional[float] = None
    core_fraction: Optional[float] = None       # frac(candidate raw_sum > bar): strong-motif U12 core
    # v6 C6 species-level BPS penalty: {frac_bp6, n_hc_pre_penalty, logit_shift, p_gate}. Populated
    # post-discount; None when the penalty is disabled or didn't fire (n_hc==0 / shift>=0). This is the
    # species-level FP-suppression lever — surfaced for the post-run watch-list (PRELAUNCH #7). Diagnostic only.
    species_penalty: Optional[dict] = None


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

    # C4: separable safe-fail + conservative_min routing.
    # (a) SAFE-FAIL: the gate checks only 5' separation, but apply_mode_separation_z runs on bp and 3' too and
    #     RAISES if their μ_U12 ≤ μ_U2 (separation ≤ 0) — a degenerate single feature on a sparse/diverged
    #     candidate cluster. Treat ANY degenerate feature as a gate failure so we fall back instead of crashing
    #     mid-scoring (task #50). Applies to ALL bundles (crash → graceful; strict improvement).
    # (b) CONSERVATIVE-MIN: under a graduated-tail bundle, a gate-FAIL but SEPARABLE species is NOT reverted to
    #     the permissive first-pass — it gets mode-sep + the graduated tail (route="conservative_min"). This
    #     matches the gate-FREE offline eval that validated the tail (FP@full-recall 11; every eligible intron
    #     scored by the tail regardless of gate), and protects recall on gate-failing divergent BEARERS
    #     (e.g. closterium_sp._nies-54, charophyte TP, fails low_centroid_sigma — fallback would drop its
    #     mode-sep recalibration). NO min(first-pass, mode-sep) is applied: the tail's negative margin coef +
    #     hard negatives supply the FP-suppression the isotonic-era min provided, and adding a min would deviate
    #     from the validated FP=11 (and risk the very divergent-bearer recall this recovers). Default
    #     (non-graduated) bundles keep the original gate-fail → first-pass fallback, byte-unchanged.
    separable = (s5.separation > 0 and sbp.separation > 0 and s3.separation > 0)
    if gate.passes and not separable:
        gate = replace(gate, passes=False, reason="degenerate_separation")
    conservative_min = (not gate.passes) and separable and bool(params.get("graduated_tail"))

    if not gate.passes and not conservative_min:
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
            gap_fraction=gate.gap_fraction,
            gap_fraction_ucl=gate.gap_fraction_ucl,
            centroid_sigma=gate.centroid_sigma,
            core_fraction=gate.core_fraction,
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
    # gate.valley_depth is Optional[float] (None on early gate returns — below_n_floor /
    # degenerate_separation — or when median_depth is non-finite). The conservative_min route
    # (gate-fail + separable + graduated) reaches this log line on a gate-fail, so None is
    # reachable here; format defensively (cf. the guard at module line ~142). Logging only —
    # no effect on scores.
    _vd_disp = gate.valley_depth if gate.valley_depth is not None else float("nan")
    _log(f"[modesep] {'CONSERVATIVE-MIN (gate-fail, separable, graduated)' if conservative_min else 'GATE-PASS'} "
         f"μ_U2={s5.mu_u2:.2f} μ_U12={s5.mu_u12:.2f} "
         f"valley={_vd_disp:.3f}; scoring {n_elig:,}/{n_total:,} "
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

    # C1: de-saturated 2nd-pass margin (mean base-SVC decision_function) for the graduated/Platt tail.
    # Computed ONLY when the bundle ships a graduated_tail (default bundles don't consume it -> skip the
    # extra decision_function pass). MUST match the SVC set the tail's scaler/coefs were fit on: the
    # training-feature builder proto_grad_features.py uses fold-[0] of each ensemble model
    # (BASE_SVCS = [m.model.calibrated_classifiers_[0].estimator for m in spp.models], n_models SVCs) — a
    # deliberate ~= mean-of-all-folds approximation. Averaging all CalibratedClassifierCV folds instead
    # (n_models × n_folds) shifts the margin by up to ~0.14 and silently miscalibrates the baked logistic
    # (V1b parity decomposition: all-folds -> adjusted_score max|Δ|≈6; fold-[0] -> max|Δ|≈5e-6). The OTHER
    # offline script proto_margin_extract.py uses all folds, but it did NOT feed phase4_build_graduated_bundle.
    # BACK-POCKET (needs a re-fit, not a hot-swap): averaging the FIRST K folds (K∈{2,3}) is a cost/quality dial
    # between fold-[0] and all-folds — variance falls ~1/K and is ~flat past 3 folds, so K=2-3 buys most of the
    # smoothing at <all-folds serve cost. Use first-K (folds are an exchangeable CV partition) for determinism;
    # NEVER a per-run random subset (would jitter calls + break -q reproducibility). Any K≠1 requires rebuilding
    # proto_grad_features.tsv + hard-neg margins at that K and re-baking/re-validating the bundle. See
    # GRADUATED_PRODUCTIONIZATION_PLAN.md "back-pocket".
    margin_all = np.full(len(fp_svm), np.nan, dtype=float)
    if params.get("graduated_tail") and n_elig > 0:
        base_svcs = [m.model.calibrated_classifiers_[0].estimator
                     for m in second_pass_ensemble.models]
        margin_all[elig] = np.mean([e.decision_function(X_all[elig]) for e in base_svcs], axis=0)

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
    if params.get("graduated_tail"):   # C1: persist the de-saturated margin for the graduated tail
        margin_map = dict(zip(update_keys, margin_all))
        df_full["svm_margin"] = df_full["name"].map(margin_map)
    df_full["modesep_route"] = "conservative_min" if conservative_min else "modesep"
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
        route="conservative_min" if conservative_min else "modesep",
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
        gap_fraction=gate.gap_fraction,
        gap_fraction_ucl=gate.gap_fraction_ucl,
        centroid_sigma=gate.centroid_sigma,
        core_fraction=gate.core_fraction,
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
    # Sanitize non-finite floats (gap_fraction = -inf for structural no-separation, or NaN)
    # to null so the diagnostics file stays valid, portable JSON. gate_reason carries the
    # categorical "why" for the gate-fail cases where these go null.
    payload = {k: (None if (isinstance(v, float) and not np.isfinite(v)) else v)
               for k, v in payload.items()}
    Path(diagnostics_path).write_text(json.dumps(payload, indent=2))


def _sync_meta_from_score_info(meta_path: Path, df: "pd.DataFrame",
                                messenger=None) -> None:
    """Propagate rel_score and type_id from score_info DataFrame to meta.iic.

    v2.7.1 introduces unified type_id labels and a discount-aware rel_score.
    meta.iic historically copies these from the in-memory Intron object,
    which is computed before the discount runs. This function syncs both
    columns from the on-disk score_info DataFrame after the discount has
    been applied.
    """
    import shutil
    import tempfile

    name_to_rel = dict(zip(df["name"].astype(str),
                           df["rel_score"].astype(float)))
    name_to_type = dict(zip(df["name"].astype(str),
                            df["type_id"].astype(str)))

    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".iic",
        dir=meta_path.parent, delete=False,
    )
    n_synced = 0
    try:
        with open(meta_path) as f_in:
            header_line = f_in.readline()
            tmp.write(header_line)
            hdr = header_line.rstrip("\n").split("\t")
            try:
                name_col = hdr.index("name")
                rel_col = hdr.index("rel_score")
                type_col = hdr.index("type_id")
            except ValueError:
                tmp.close()
                Path(tmp.name).unlink(missing_ok=True)
                return

            for line in f_in:
                parts = line.rstrip("\n").split("\t")
                if name_col < len(parts):
                    nm = parts[name_col]
                    if nm in name_to_rel:
                        parts[rel_col] = f"{name_to_rel[nm]:.4f}"
                        parts[type_col] = name_to_type[nm]
                        n_synced += 1
                tmp.write("\t".join(parts) + "\n")
        tmp.close()
        import os
        orig_mode = meta_path.stat().st_mode if meta_path.exists() else 0o644
        shutil.move(tmp.name, meta_path)
        os.chmod(meta_path, orig_mode & 0o777)
        if messenger is not None:
            messenger.info(
                f"[meta-sync] rewrote rel_score + type_id for "
                f"{n_synced} rows in {meta_path.name}"
            )
    except Exception:
        tmp.close()
        Path(tmp.name).unlink(missing_ok=True)
        raise


# Observability thresholds for the discount-clawback monitor (pure diagnostic; do not affect scores).
# WARN when the discount claws back a large svm-level reservoir (n_pre >> n_post) in BOTH absolute and
# relative terms — flags species whose final U12 set leans heavily on the discount (the tail-risk case
# under discount drift). Not an error: gate-passing bearers legitimately carry a reservoir.
RESERVOIR_WARN_ABS = 100     # absolute clawback floor (n_pre - n_post)
RESERVOIR_WARN_RATIO = 2.0   # and n_pre >= RATIO * n_post


def apply_continuous_per_intron_discount(
    score_info_path: Path,
    *,
    threshold: float = DEFAULT_THRESHOLD,
    k_overcall: float = DEFAULT_DISCOUNT_K_OVERCALL,
    tau_overcall: float = DEFAULT_DISCOUNT_TAU_OVERCALL,
    k_weakmot: float = DEFAULT_DISCOUNT_K_WEAKMOT,
    tau_motif: float = DEFAULT_DISCOUNT_TAU_MOTIF,
    input_column: str = "svm_score",
    meta_path: Optional[Path] = None,
    messenger=None,
    species_gap_fraction: Optional[float] = None,
    species_gap_fraction_ucl: Optional[float] = None,
    graduated_tail: Optional[dict] = None,
    enable_species_penalty: bool = False,
    species_penalty_pgate: float = DEFAULT_SPECIES_PENALTY_PGATE,
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

    def _warn(msg):
        if messenger is not None:
            messenger.warning(msg)

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

    # C3: graduated/Platt tail. When the bundle ships a `graduated_tail`, replace the discount on the
    # mode-sep-scored (margin-bearing) introns with the smooth logistic on [margin, raw5/bp/3, adh5z/adhbpz]:
    #     adjusted_score = 100 * sigmoid(coef . scale([margin, raw5, rawbp, raw3, adh5z, adhbpz]) + intercept)
    # Rows without `svm_margin` (ineligible or first_pass_fallback species) keep the discount computed above —
    # the unchanged fallback path. No graduated_tail -> default discount everywhere (byte-identical to before).
    if graduated_tail and "svm_margin" in df.columns:
        from intronIC.scoring.gold_adherence import adh_features
        gt = graduated_tail
        mcol = pd.to_numeric(df["svm_margin"], errors="coerce").to_numpy()
        gmask = ~np.isnan(mcol)
        if "spac" in gt["feature_order"]:          # spac (BP->3'SS geometry) feature needs bp_offset
            bo = (pd.to_numeric(df["bp_offset"], errors="coerce").to_numpy()
                  if "bp_offset" in df.columns else np.full(len(df), np.nan))
            gmask = gmask & ~np.isnan(bo)          # no bp_offset -> can't compute spac -> keep discount fallback
        if gmask.any():
            a5z, abpz = adh_features(df, gmask, want5=True)
            colmap = {"margin": mcol[gmask], "raw5": df["5'_raw"].to_numpy()[gmask],
                      "rawbp": df["bp_raw"].to_numpy()[gmask], "raw3": df["3'_raw"].to_numpy()[gmask],
                      "adh5z": a5z, "adhbpz": abpz}
            if "spac" in gt["feature_order"]:
                colmap["spac"] = np.exp(-((np.abs(bo[gmask]) - 12.0) / 8.0) ** 2)
            X = np.column_stack([colmap[f] for f in gt["feature_order"]])
            z = (X - np.asarray(gt["scaler_mean"], float)) / np.asarray(gt["scaler_scale"], float)
            adj[gmask] = 100.0 / (1.0 + np.exp(-(z @ np.asarray(gt["coef"], float) + float(gt["intercept"]))))

    # C6: species-level BPS penalty (frac_bp6 + logN), gated OFF by default. Applied AFTER the per-intron
    # discount/graduated tail, so the HC signature (frac_bp6, n_hc) is the species' POST-scoring call set —
    # gate/tail-suppressed losses arrive with low HC, sidestepping the metric's high-N blind spot. One-sided
    # (<=0): confident bearers (p_bearer>=p_gate) are untouched; loss-prone species get a negative logit shift
    # on adjusted_score. Adds species-level FP-suppression the per-intron tail structurally lacks (prove-out:
    # FP@full-recall 10→6, recall-positive, 0 TP lost on the v6.1 panel). The conservative_min min-of-discounts
    # is deliberately NOT ported — the graduated tail subsumes it (see C4). No graduated_tail / penalty disabled
    # -> this block is skipped (default bundle byte-unchanged). See GRADUATED_PRODUCTIONIZATION_PLAN.md.
    species_penalty = None
    if enable_species_penalty:
        hc_mask = adj >= threshold
        n_hc = int(hc_mask.sum())
        if n_hc > 0:
            bp_for_pen = pd.to_numeric(df.get("bp_raw"), errors="coerce").to_numpy()
            frac_bp6 = float(np.nanmean((bp_for_pen[hc_mask] >= 6).astype(float)))
            shift = species_penalty_logit_shift(frac_bp6, n_hc, p_gate=species_penalty_pgate)
            if shift < 0.0:
                p_adj = np.clip(adj / 100.0, 1e-9, 1 - 1e-9)
                adj = 100.0 / (1.0 + np.exp(-(np.log(p_adj / (1 - p_adj)) + shift)))
            species_penalty = {"frac_bp6": round(frac_bp6, 4), "n_hc_pre_penalty": n_hc,
                               "logit_shift": round(shift, 4), "p_gate": species_penalty_pgate}
            _log(f"[species-penalty] frac_bp6={frac_bp6:.3f} n_hc={n_hc} shift={shift:+.3f} "
                 f"-> HC {n_hc}→{int((adj >= threshold).sum())} (P_GATE={species_penalty_pgate})")

    n_post = int((adj >= threshold).sum())

    df["adjusted_score"] = adj
    # v2.7.1 fix: rel_score must track adjusted_score (the canonical
    # calling column post-discount). Without this, gate-pass species
    # ship a stale rel_score that reflects first_pass_svm or the
    # pre-discount mode-sep score, undercounting U12s in downstream
    # filters that use `rel_score > 0`.
    df["rel_score"] = np.round(adj - threshold, 4)

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

    # v3 gate-gapfrac (Step 5): surface the per-species separation statistic that drove the
    # π_species prior as constant columns (route-appropriate value passed by the caller — the
    # gate's UCL on the modesep route, validate_u12_cluster's on fallback). Diagnostic /
    # traceability only; the harness can also recompute from the 5'_z/bp_z/3'_z columns.
    if species_gap_fraction is not None:
        df["gap_fraction"] = species_gap_fraction
    if species_gap_fraction_ucl is not None:
        df["gap_fraction_ucl"] = species_gap_fraction_ucl

    df.to_csv(score_info_path, sep="\t", index=False, float_format="%.6f")
    _log(f"[continuous-discount] applied; n_called_pre={n_pre}, "
         f"n_called_post={n_post}, suppressed={n_pre - n_post}")

    # Discount-clawback monitor (pure observability — does not alter scores). Flags species whose
    # final U12 set leans heavily on the discount clawing back a large svm-level reservoir; the
    # tail-risk case if the discount ever drifts (see docs/threshold_shift_refactor_proposal.md).
    _reservoir = n_pre - n_post
    if _reservoir >= RESERVOIR_WARN_ABS and n_pre >= RESERVOIR_WARN_RATIO * max(n_post, 1):
        _warn(
            f"[discount-monitor] large clawback: svm-level {n_pre} → HC {n_post} "
            f"(reservoir {_reservoir}, {100 * _reservoir / max(n_pre, 1):.0f}% of svm-level calls); "
            f"this species' U12 calls rely heavily on the discount — verify if unexpected "
            f"(tail-risk under discount drift)."
        )

    # v2.7.1 fix: propagate updated rel_score and unified type_id to
    # meta.iic. Without this, meta.iic shipped stale rel_score (and a
    # legacy 50%-threshold type_id) that disagreed with score_info.iic
    # for every gate-pass species.
    if meta_path is not None and Path(meta_path).exists():
        _sync_meta_from_score_info(
            Path(meta_path), df, messenger=messenger
        )

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
        "species_penalty": species_penalty,
        **label_counts,
    }


__all__ = [
    "apply_mode_separation_postprocess",
    "apply_continuous_per_intron_discount",
    "ModeSepResult",
]
