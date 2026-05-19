"""Phase 4.3 regression test for v2.6 mode-separation on the 14-species panel.

Runs the v5 modesep bundle's post-processing pipeline against each panel
species' canonical cluster-aware score_info.iic and asserts TP/FN/FP_strong/
FP_any counts match the validated baseline (±1 per species, ±5 totals).

Slow (~3-4 minutes for full panel). Skipped if panel data or gold labels
are missing. Run via:

    pytest tests/integration/test_modesep_panel_regression.py -v -s

The baseline numbers come from the offline panel evaluation in
`multispecies_v23/mode_separation/eval_panel_modesep.py` (and were used in
the consultant-facing brief MODESEP_BRIEF_FOR_CONSULT.md).
"""
import shutil
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

warnings.filterwarnings("ignore")

BUNDLE = Path("/mnt/data/u12/ipa/training_data/multispecies_v23/"
              "v5_modesep_v2.6_bundle.model.pkl")
PANEL_DIR = Path("/mnt/data/u12/ipa/training_data/v25_eval/config_A_cluster_aware")
GOLD_ALL = Path("/mnt/data/u12/ipa/training_data/gold_labels_v2.3.tsv")
PANEL_GOLD = Path("/mnt/data/u12/ipa/training_data/v25_eval/gold_panel/panel_gold.tsv")

pytestmark = pytest.mark.skipif(
    not (BUNDLE.exists() and PANEL_DIR.is_dir() and GOLD_ALL.exists()),
    reason="v5 bundle / panel data / gold labels not present",
)

# Baseline from validated offline panel run. Values: (TP, FN, FP_strong, FP_any,
# expected_route).
BASELINE = {
    "HomSap": (496, 0, 0,  34, "modesep"),
    "DanRer": (478, 1, 0,   8, "modesep"),
    "GalGal": (412, 1, 0,  72, "modesep"),
    "XenTro": (426, 0, 0,  86, "modesep"),
    "DroMel": (  8, 0, 0,   2, "modesep"),
    "AraTha": ( 47, 1, 0,   7, "modesep"),
    "AmbTri": ( 51, 0, 0,   7, "modesep"),
    "OrySat": ( 32, 0, 0,   0, "modesep"),
    "CaeEle": (  0, 0, 0,   0, "first_pass_fallback"),
    "AscSuu": (  0, 0, 0,   1, "first_pass_fallback"),
    "SchPom": (  0, 0, 0,   0, "first_pass_fallback"),
    "SacCer": (  0, 0, 0,   0, "first_pass_fallback"),
    "TetThe": (  0, 0, 0,   0, "first_pass_fallback"),
    "ChlRei": (  0, 0, 0,   0, "first_pass_fallback"),
}

TOLERANCE_TP_FN_PER_SPECIES = 1
TOLERANCE_FP_ANY_PER_SPECIES = 3


@pytest.fixture(scope="module")
def runtime_and_gold():
    from intronIC.utils.model_io import load_model, normalize_model_bundle
    raw = load_model(BUNDLE)
    raw["_source_path"] = str(BUNDLE)
    runtime = normalize_model_bundle(raw)

    gold = pd.read_csv(GOLD_ALL, sep="\t",
                       usecols=["intron_name", "species", "label", "alignment_id"])
    panel_gold = pd.read_csv(PANEL_GOLD, sep="\t")
    strong_set = set(panel_gold.loc[panel_gold["is_strong_u2"], "alignment_id"].dropna())
    return runtime, gold, strong_set


def _evaluate_species(sp: str, runtime, gold, strong_set,
                       apply_discount: bool = True):
    """Evaluate a species. apply_discount=True mirrors v2.7 production flow
    (mode-sep + continuous discount); False mirrors v2.6 (mode-sep only)."""
    from intronIC.classification.mode_sep_pipeline import (
        apply_mode_separation_postprocess,
        apply_continuous_per_intron_discount,
    )
    from intronIC.scoring.cluster_validation import compute_valley_depth

    src = PANEL_DIR / sp / f"{sp}.score_info.iic"
    if not src.is_file():
        pytest.skip(f"missing {sp} score_info")

    with tempfile.TemporaryDirectory() as td:
        dst = Path(td) / f"{sp}.score_info.iic"
        shutil.copy(src, dst)
        result = apply_mode_separation_postprocess(
            score_info_path=dst,
            runtime=runtime,
            valley_depth_fn=compute_valley_depth,
            threshold=90.0,
        )
        # v2.7: apply continuous discount after mode-sep.
        # On gate-fail, skip — the production cli/main.py applies the
        # legacy valley-depth discount first; this regression test exercises
        # only the modesep + continuous-discount portion.
        call_col = "svm_score"
        if apply_discount and result.route == "modesep":
            apply_continuous_per_intron_discount(
                score_info_path=dst, threshold=90.0,
                input_column="svm_score",
            )
            call_col = "adjusted_score"
        df = pd.read_csv(dst, sep="\t", low_memory=False)
    df["key"] = df["name"].astype(str).str.split(";", n=1).str[0].str.lower()
    df[call_col] = pd.to_numeric(df[call_col], errors="coerce")
    df["called"] = df[call_col] >= 90.0

    sp_gold = gold[gold["species"] == sp].copy()
    sp_gold["key"] = sp_gold["intron_name"].astype(str).str.lower()
    m = df.merge(sp_gold[["key", "label", "alignment_id"]], on="key", how="left")
    is_u12 = m["label"] == "u12"
    is_u2 = m["label"] == "u2"
    is_strong = is_u2 & m["alignment_id"].isin(strong_set)

    return {
        "TP": int((is_u12 & m["called"]).sum()),
        "FN": int((is_u12 & ~m["called"]).sum()),
        "FP_strong": int((is_strong & m["called"]).sum()),
        "FP_any": int((is_u2 & m["called"]).sum()),
        "route": result.route,
    }


@pytest.mark.parametrize("species", list(BASELINE.keys()))
def test_panel_species_regression(species, runtime_and_gold):
    runtime, gold, strong_set = runtime_and_gold
    metrics = _evaluate_species(species, runtime, gold, strong_set)
    exp_tp, exp_fn, exp_fps, exp_fpa, exp_route = BASELINE[species]

    assert metrics["route"] == exp_route, (
        f"{species}: route changed (was {exp_route!r}, now {metrics['route']!r})"
    )
    assert abs(metrics["TP"] - exp_tp) <= TOLERANCE_TP_FN_PER_SPECIES, (
        f"{species}: TP drift {exp_tp} → {metrics['TP']}"
    )
    assert abs(metrics["FN"] - exp_fn) <= TOLERANCE_TP_FN_PER_SPECIES, (
        f"{species}: FN drift {exp_fn} → {metrics['FN']}"
    )
    # FP_strong is a hard constraint — leakage onto strong-u2 clusters is
    # NEVER acceptable.
    assert metrics["FP_strong"] == exp_fps, (
        f"{species}: FP_strong drift {exp_fps} → {metrics['FP_strong']} "
        "(strong-u2 leakage is a hard regression)"
    )
    assert abs(metrics["FP_any"] - exp_fpa) <= TOLERANCE_FP_ANY_PER_SPECIES, (
        f"{species}: FP_any drift {exp_fpa} → {metrics['FP_any']}"
    )


def test_panel_aggregate_recall_holds(runtime_and_gold):
    """Aggregate recall across u12-bearing panel species ≥ 0.998 (Phase 0)."""
    runtime, gold, strong_set = runtime_and_gold
    u12_bearing = ["HomSap", "DanRer", "GalGal", "XenTro", "DroMel",
                   "AraTha", "AmbTri", "OrySat"]
    tp_sum = fn_sum = 0
    for sp in u12_bearing:
        m = _evaluate_species(sp, runtime, gold, strong_set)
        tp_sum += m["TP"]
        fn_sum += m["FN"]
    recall = tp_sum / (tp_sum + fn_sum)
    assert recall >= 0.998, (
        f"aggregate panel recall regressed: {recall:.4f} (baseline 0.9985)"
    )
