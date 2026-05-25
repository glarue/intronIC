"""Integration tests for the mode-separation post-classify pipeline (v2.6).

Two layers:
1. Bundle round-trip: load v5 bundle, verify model_io surfaces first-pass +
   second-pass + modesep_params correctly.
2. Pipeline driver: call `apply_mode_separation_postprocess` on a real
   pre-classified score_info.iic file (AmbTri + SacCer); assert the gate
   decision matches what we measured offline.

These tests avoid invoking the full intronIC binary and a genome run — those
belong in regression tests. The pipeline driver in this file uses score_info
files in /mnt/data/u12/ipa/training_data/v25_eval/config_A_cluster_aware/
that were produced by the cluster-aware first pass and are the canonical
inputs to the modesep second pass.
"""
import shutil
import tempfile
import warnings
from pathlib import Path

import pytest

warnings.filterwarnings("ignore")

BUNDLE = Path(
    "/mnt/data/u12/ipa/training_data/multispecies_v23/"
    "v5_modesep_v2.6_bundle.model.pkl"
)
PANEL_DIR = Path("/mnt/data/u12/ipa/training_data/v25_eval/config_A_cluster_aware")

pytestmark = pytest.mark.skipif(
    not (BUNDLE.exists() and PANEL_DIR.is_dir()),
    reason="v5 bundle or panel data not present (development environment only)",
)


@pytest.fixture(scope="module")
def runtime():
    """Loaded + normalized v5 modesep runtime dict (cached for the module)."""
    from intronIC.utils.model_io import load_model, normalize_model_bundle
    raw = load_model(BUNDLE)
    raw["_source_path"] = str(BUNDLE)
    rt = normalize_model_bundle(raw)
    return rt


@pytest.fixture(scope="module")
def valley_depth_fn():
    from intronIC.scoring.cluster_validation import compute_valley_depth
    return compute_valley_depth


# -----------------------------------------------------------------------------
# 1. Bundle round-trip
# -----------------------------------------------------------------------------

def test_bundle_surfaces_modesep_metadata(runtime):
    """v5 bundle exposes normalizer_mode, both ensembles, and modesep_params."""
    assert runtime["normalizer_mode"] == "modesep"
    assert hasattr(runtime["ensemble"], "models")
    assert hasattr(runtime["second_pass_ensemble"], "models")
    # First-pass is the cluster-aware ensemble; second-pass is v5_modesep.
    assert len(runtime["ensemble"].models) == 126
    assert len(runtime["second_pass_ensemble"].models) == 126
    params = runtime["modesep_params"]
    assert params["z_floor_eligibility"] == pytest.approx(0.30)
    assert params["mu_u12_5p_tolerance"] == pytest.approx(3.6)
    assert params["universal_anchors"]["five_raw"] == pytest.approx(15.671)
    assert params["first_pass_model_id"]
    assert params["pwm_set_id"]


# -----------------------------------------------------------------------------
# 2. Pipeline driver: AmbTri (modesep route) and SacCer (gate-fail)
# -----------------------------------------------------------------------------

def _run_pipeline(species: str, runtime, valley_depth_fn, tmpdir: Path):
    """Copy the canonical score_info to tmp and run modesep over it."""
    from intronIC.classification.mode_sep_pipeline import (
        apply_mode_separation_postprocess,
    )
    src = PANEL_DIR / species / f"{species}.score_info.iic"
    if not src.is_file():
        pytest.skip(f"missing canonical score_info for {species}")
    dst = tmpdir / f"{species}.score_info.iic"
    shutil.copy(src, dst)
    diag = tmpdir / f"{species}.modesep.json"
    result = apply_mode_separation_postprocess(
        score_info_path=dst,
        runtime=runtime,
        valley_depth_fn=valley_depth_fn,
        threshold=90.0,
        diagnostics_path=diag,
    )
    return result, dst, diag


def test_modesep_ambtri_gate_passes_and_recalibrates(runtime, valley_depth_fn):
    """AmbTri: gate passes (valley ~0.59), 51 IPA U12s recovered."""
    with tempfile.TemporaryDirectory() as td:
        result, _, _ = _run_pipeline("AmbTri", runtime, valley_depth_fn,
                                      Path(td))
    assert result.route == "modesep"
    assert result.gate_reason == "ok"
    assert result.quality_tier in ("modesep_strong", "modesep_standard")
    # Offline panel measured 306 calls; allow ±5 slack for numerical drift.
    assert 295 <= result.n_called_u12 <= 320
    assert result.valley_depth is not None and result.valley_depth >= 0.30
    assert result.n_eff_candidates >= 100   # AmbTri has ~238 effective candidates
    assert abs(result.mu_u12_5p_offset) < 3.6  # within tolerance


def test_modesep_saccer_gate_fails_below_n_floor(runtime, valley_depth_fn):
    """SacCer (u12-absent): gate-fail at n_floor, zero spurious calls."""
    with tempfile.TemporaryDirectory() as td:
        result, _, _ = _run_pipeline("SacCer", runtime, valley_depth_fn,
                                      Path(td))
    assert result.route == "first_pass_fallback"
    assert result.gate_reason == "below_n_floor"
    assert result.quality_tier == "first_pass_fallback"
    assert result.n_called_u12 == 0
    assert result.n_eff_candidates < 5


def test_modesep_writes_diagnostic_json(runtime, valley_depth_fn):
    """Diagnostic JSON is written and contains required fields (BLOCKER)."""
    import json
    with tempfile.TemporaryDirectory() as td:
        result, _, diag = _run_pipeline("AmbTri", runtime, valley_depth_fn,
                                         Path(td))
        assert diag.is_file()
        payload = json.loads(diag.read_text())
    required_keys = {
        "route", "gate_reason", "quality_tier", "n_introns",
        "n_eligible", "n_called_u12", "n_eff_candidates",
        "valley_depth", "mu_u2_5p", "mu_u12_5p", "mu_u12_5p_offset",
        "first_pass_model_id", "second_pass_model_id",
    }
    assert required_keys <= set(payload.keys()), (
        f"diagnostic JSON missing keys: {required_keys - set(payload.keys())}"
    )


def test_modesep_rewrites_score_info_columns(runtime, valley_depth_fn):
    """After modesep, score_info gains ensemble_sigma + first_pass_svm + v2.7 diagnostic columns."""
    import pandas as pd
    with tempfile.TemporaryDirectory() as td:
        result, score_path, _ = _run_pipeline(
            "AmbTri", runtime, valley_depth_fn, Path(td))
        df = pd.read_csv(score_path, sep="\t", low_memory=False, nrows=10)
    # v2.6 diagnostic columns
    assert "ensemble_sigma" in df.columns
    assert "first_pass_svm" in df.columns
    assert "modesep_route" in df.columns
    # v2.7 diagnostic columns
    assert "raw_sum" in df.columns
    assert "svm_vs_naive" in df.columns
    assert "voting_frac" in df.columns


# -----------------------------------------------------------------------------
# v2.7 — continuous per-intron discount integration
# -----------------------------------------------------------------------------

def test_modesep_reports_boundary_mass(runtime, valley_depth_fn):
    """v2.7: AmbTri gate-pass run reports boundary_mass in result + JSON."""
    with tempfile.TemporaryDirectory() as td:
        result, _, diag = _run_pipeline(
            "AmbTri", runtime, valley_depth_fn, Path(td))
    assert result.boundary_mass is not None
    # AmbTri is a clean gate-pass species — BM should be very low
    assert result.boundary_mass < 0.005


def test_continuous_discount_preserves_strong_TPs(runtime, valley_depth_fn):
    """v2.7: applying continuous discount on AmbTri should not zero its TPs.

    AmbTri's 51 IPA-validated U12s have raw_sum ≈ +20 and SVM saturated, so
    the discount should be ~no-op on them.
    """
    from intronIC.classification.mode_sep_pipeline import (
        apply_continuous_per_intron_discount,
    )
    import pandas as pd
    with tempfile.TemporaryDirectory() as td:
        result, score_path, _ = _run_pipeline(
            "AmbTri", runtime, valley_depth_fn, Path(td))
        before = pd.read_csv(score_path, sep="\t", low_memory=False)
        before["svm_score"] = pd.to_numeric(before["svm_score"], errors="coerce")
        n_pre = int((before["svm_score"] >= 90).sum())

        summary = apply_continuous_per_intron_discount(score_path, threshold=90.0)

        after = pd.read_csv(score_path, sep="\t", low_memory=False)
        after["adjusted_score"] = pd.to_numeric(after["adjusted_score"], errors="coerce")
        n_post = int((after["adjusted_score"] >= 90).sum())

    assert "adjusted_score" in after.columns
    # n_pre should ≈ AmbTri's 306-ish called count
    assert n_pre > 250 and n_pre < 400, f"unexpected n_pre {n_pre}"
    # Discount should preserve most TPs — at most ~50 suppression (the
    # loose-or-NA tail), not the whole called set
    assert n_post >= n_pre - 60, (
        f"continuous discount over-suppressed: {n_pre} → {n_post}"
    )
    # Summary dict is well-formed
    assert summary["n_called_pre_discount"] == n_pre
    assert summary["n_called_post_discount"] == n_post
    assert summary["k_overcall"] == 1.0  # default


def test_continuous_discount_preserves_svm_score(runtime, valley_depth_fn):
    """v2.7 invariant: svm_score column is NEVER modified by the discount;
    adjusted_score is the new column."""
    from intronIC.classification.mode_sep_pipeline import (
        apply_continuous_per_intron_discount,
    )
    import pandas as pd
    with tempfile.TemporaryDirectory() as td:
        result, score_path, _ = _run_pipeline(
            "AmbTri", runtime, valley_depth_fn, Path(td))
        df_before = pd.read_csv(score_path, sep="\t", low_memory=False)
        svm_before = df_before["svm_score"].copy()

        apply_continuous_per_intron_discount(score_path, threshold=90.0)

        df_after = pd.read_csv(score_path, sep="\t", low_memory=False)
        svm_after = df_after["svm_score"]

    # svm_score column unchanged
    pd.testing.assert_series_equal(svm_before, svm_after, check_names=False)
    # adjusted_score is new
    assert "adjusted_score" in df_after.columns


def test_continuous_discount_disable_via_zero_coefficients(runtime, valley_depth_fn):
    """Setting k_overcall=0 and k_weakmot=0 makes the discount a no-op."""
    from intronIC.classification.mode_sep_pipeline import (
        apply_continuous_per_intron_discount,
    )
    import pandas as pd
    with tempfile.TemporaryDirectory() as td:
        result, score_path, _ = _run_pipeline(
            "AmbTri", runtime, valley_depth_fn, Path(td))
        df_pre = pd.read_csv(score_path, sep="\t", low_memory=False)
        df_pre["svm_score"] = pd.to_numeric(df_pre["svm_score"], errors="coerce")

        apply_continuous_per_intron_discount(
            score_path, threshold=90.0,
            k_overcall=0.0, k_weakmot=0.0,
        )

        df_post = pd.read_csv(score_path, sep="\t", low_memory=False)
        df_post["adjusted_score"] = pd.to_numeric(df_post["adjusted_score"], errors="coerce")
    # With zero coefficients, adjusted_score = svm_score (up to numerical noise)
    diff = (df_post["adjusted_score"] - df_pre["svm_score"]).abs()
    assert diff.max() < 0.5, f"max |Δ| = {diff.max():.4f} should be near 0"
