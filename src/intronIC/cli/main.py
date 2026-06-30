"""
Main CLI entry point for intronIC.

Orchestrates the complete intron classification pipeline.
"""

from __future__ import annotations

import gc
import json
import logging
import shutil
import sys
import tempfile
import time
from dataclasses import dataclass, replace
from multiprocessing import Pool
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

import joblib
import numpy as np
from smart_open import open as smart_open  # type: ignore[import-unresolved]


# Import pipeline components
from intronIC.core.intron import Intron, IntronMetadata, IntronSequences, OmissionReason
from intronIC.extraction.annotator import AnnotationHierarchyBuilder
from intronIC.extraction.filters import IntronFilter, prefilter_introns
from intronIC.extraction.intronator import IntronGenerator
from intronIC.extraction.sequences import SequenceExtractor
from intronIC.file_io.genome import GenomeReader
from intronIC.file_io.parsers import BEDParser, SequenceParser
from intronIC.file_io.writers import (
    BEDWriter,
    MetaWriter,
    ScoreWriter,
    SequenceWriter,
    count_calls_from_meta,
    summarize_boundaries_from_meta,
    summarize_introns,
)
from intronIC.scoring.scorer import IntronScorer
from intronIC.utils.coordinates import GenomicCoordinate
from intronIC.utils.metadata import (
    generate_pretrained_metadata,
    generate_training_metadata,
    write_metadata,
)
from intronIC.utils.model_io import (
    load_model,
    normalize_model_bundle,
    assert_scoreable_bundle,
    UnsupportedBundleError,
)
from intronIC.visualization.plots import plot_classification_results

from .args import IntronICArgumentParser
from .config import IntronICConfig, ScoringRegions
from .progress import IntronICProgressReporter

if TYPE_CHECKING:
    from rich.console import Console

    from intronIC.classification.trainer import SVMEnsemble
    from intronIC.cli.messenger import UnifiedMessenger
    from intronIC.scoring.pwm import PWMSet

# ============================================================================
# Helper Functions
# ============================================================================


def get_pwm_file(config: "IntronICConfig") -> Path:
    """
    Get PWM file path from config or use default.

    Args:
        config: IntronIC configuration

    Returns:
        Path to PWM file (supports .iic, .yaml, and .json formats)

    Raises:
        FileNotFoundError: If PWM file doesn't exist
    """
    if config.scoring.pwm_file is not None:
        pwm_file = config.scoring.pwm_file
    else:
        # Default to JSON format
        data_dir = Path(__file__).parent.parent / "data"
        pwm_file = data_dir / "intronIC_scoring_PWMs.json"

    if not pwm_file.exists():
        raise FileNotFoundError(f"PWM file not found: {pwm_file}")

    return pwm_file


def merge_pwm_sets(
    default_pwms: Dict[str, "PWMSet"], custom_pwms: Dict[str, "PWMSet"]
) -> Dict[str, "PWMSet"]:
    """
    Merge custom PWM matrices into defaults (custom overwrites defaults).

    Args:
        default_pwms: Default PWM sets
        custom_pwms: Custom PWM sets to merge in

    Returns:
        Merged PWM sets (modifies default_pwms in place and returns it)
    """
    for region in ["five", "bp", "three"]:
        if region in custom_pwms and region in default_pwms:
            custom_matrices = custom_pwms[region].matrices
            default_matrices = default_pwms[region].matrices

            # Update defaults with custom (custom overwrites)
            for key, pwm in custom_matrices.items():
                default_matrices[key] = pwm

    return default_pwms


def _create_background_accumulator(
    config: "IntronICConfig",
    pwm_sets: Dict[str, "PWMSet"],
) -> Tuple["SpeciesBackground", Dict[str, "PWMSet"], Dict[str, "PWMSet"], int, int]:
    """Create a SpeciesBackground accumulator and split PWMs into U12/U2.

    Shared setup used by all three BG correction paths (in-memory introns,
    SQLite readback, and streaming genome fetch).

    Args:
        config: Pipeline config with species_background settings
        pwm_sets: Original PWM sets (U12 + U2 combined)

    Returns:
        (bg, u12_sets, u2_sets, five_len, three_len)
    """
    from intronIC.scoring.background import BackgroundConfig, SpeciesBackground
    from intronIC.scoring.pwm import PWMSet

    bg_config = config.species_background

    # Split into U12 and U2 sets
    u12_sets: Dict[str, "PWMSet"] = {}
    u2_sets: Dict[str, "PWMSet"] = {}
    for region in ["five", "three", "bp"]:
        if region not in pwm_sets:
            continue
        u12_mats = {
            k: v for k, v in pwm_sets[region].matrices.items() if k[0] == "u12"
        }
        u2_mats = {
            k: v for k, v in pwm_sets[region].matrices.items() if k[0] == "u2"
        }
        u12_sets[region] = PWMSet(matrices=u12_mats)
        u2_sets[region] = PWMSet(matrices=u2_mats)

    five_len = (
        config.scoring.scoring_regions.five_end
        - config.scoring.scoring_regions.five_start
    )
    three_len = (
        config.scoring.scoring_regions.three_end
        - config.scoring.scoring_regions.three_start
    )

    core_config = BackgroundConfig(
        enabled=True,
        n0=bg_config.n0,
        trim_percentile=bg_config.trim_percentile,
        pseudocount_per_base=bg_config.pseudocount_per_base,
        n_iterations=bg_config.n_iterations,
        min_introns=bg_config.min_introns,
    )

    bg = SpeciesBackground(
        human_u2_pwm_sets=u2_sets,
        config=core_config,
        five_len=five_len,
        three_len=three_len,
        five_start=config.scoring.scoring_regions.five_start,
        three_start=config.scoring.scoring_regions.three_start,
    )

    return bg, u12_sets, u2_sets, five_len, three_len


def _finalize_background_correction(
    bg: "SpeciesBackground",
    u12_sets: Dict[str, "PWMSet"],
    n_accumulated: int,
    min_introns: int,
    n0: int,
    messenger: "UnifiedMessenger",
) -> Optional[Dict[str, "PWMSet"]]:
    """Build corrected PWMs if enough introns were accumulated.

    Args:
        bg: SpeciesBackground with accumulated frequencies
        u12_sets: U12 PWM sets (passed through unchanged)
        n_accumulated: Number of introns accumulated
        min_introns: Minimum required for correction
        n0: Bayesian shrinkage strength (for log message)
        messenger: For logging

    Returns:
        Corrected PWM sets, or None if insufficient introns
    """
    if n_accumulated < min_introns:
        messenger.info(
            f"Species background: {n_accumulated} introns < min_introns "
            f"({min_introns}), skipping correction"
        )
        return None

    messenger.info(
        f"Species background: computing empirical U2 from {n_accumulated} introns "
        f"(n0={n0})"
    )
    return bg.build_corrected_pwm_sets(u12_pwm_sets=u12_sets)


def apply_species_background(
    introns: List["Intron"],
    pwm_sets: Dict[str, "PWMSet"],
    config: "IntronICConfig",
    messenger: "UnifiedMessenger",
) -> Dict[str, "PWMSet"]:
    """Apply species-specific U2 background correction to PWM sets.

    Computes empirical per-position nucleotide frequencies from the species'
    own introns and blends them into the U2 PWMs. The U12 PWMs are unchanged.

    This must be called AFTER extraction (introns have sequences) and BEFORE
    scoring (scorer uses the returned corrected pwm_sets).

    Args:
        introns: Introns with extracted sequences
        pwm_sets: Original PWM sets (U12 + U2)
        config: Pipeline config with species_background settings
        messenger: For logging

    Returns:
        Corrected PWM sets (same structure, U2 entries replaced)
    """
    bg_config = config.species_background
    if not bg_config.enabled:
        return pwm_sets

    bg, u12_sets, _, five_len, three_len = _create_background_accumulator(
        config, pwm_sets
    )

    # Accumulate sequences from all introns
    five_start = config.scoring.scoring_regions.five_start
    five_end = config.scoring.scoring_regions.five_end
    three_start = config.scoring.scoring_regions.three_start
    three_end = config.scoring.scoring_regions.three_end
    bp_start = config.scoring.scoring_regions.bp_start
    bp_end = config.scoring.scoring_regions.bp_end

    include_duplicates = config.extraction.include_duplicates

    n_accumulated = 0
    for intron in introns:
        seqs = intron.sequences
        if seqs is None or seqs.seq is None:
            continue

        # Skip coordinate-duplicate isoforms (they share sequences with the
        # first occurrence via extract_sequences_with_deduplication; counting
        # them here would over-weight genes with many alternative
        # transcripts in the empirical U2 background and produce different
        # corrected PWMs than the streaming-classify path, which applies the
        # same dedup via IntronFilter after boundary correction).
        if not include_duplicates and intron.metadata and intron.metadata.duplicate:
            continue

        five_dnt = seqs.five_prime_dnt or (seqs.seq[:2] if seqs.seq else "")
        three_dnt = seqs.three_prime_dnt or (seqs.seq[-2:] if seqs.seq else "")
        if not five_dnt or not three_dnt:
            continue

        upstream = seqs.upstream_flank or ""
        seq = seqs.seq or ""
        downstream = seqs.downstream_flank or ""

        five_seq = upstream[five_start:] + seq[:five_end]
        three_seq = seq[three_start:] + downstream[:three_end]
        bp_region = seq[max(0, len(seq) + bp_start) : len(seq) + bp_end]

        bg.accumulate(
            intron.intron_id,
            five_dnt.upper(),
            three_dnt.upper(),
            five_seq.upper(),
            three_seq.upper(),
            bp_region.upper(),
        )
        n_accumulated += 1

    corrected = _finalize_background_correction(
        bg, u12_sets, n_accumulated, bg_config.min_introns, bg_config.n0, messenger
    )
    return corrected if corrected is not None else pwm_sets


def _sync_calls_to_meta_and_bed(score_path: Path, meta_path: Optional[Path],
                                bed_path: Optional[Path], messenger=None) -> None:
    """Propagate the final calls from a rewritten ``score_info.iic`` to ``meta.iic`` and ``bed.iic``.

    The file-side post-process scoring modes (``raw_gated``, ``pmotif_adjudicated``) compute their calls
    AFTER meta.iic/bed.iic are already on disk (both the streaming and in-memory paths write all output
    files from the in-memory introns' first-pass values, then run the shared post-classification pipeline).
    Without this sync those files keep first-pass ``type_id`` / ``rel_score`` (meta) and ``svm_score`` (bed
    col 4) — which silently corrupts the ``metrics.iic.json`` boundary tables (read from meta ``type_id``)
    whenever the post-process flips a call. Mirrors what the modesep/discount path already does.

    Updates: meta ``rel_score`` + ``type_id`` (via ``file_io.meta_sync.sync_meta_from_score_info``); bed
    col 4 (score) <- ``adjusted_score``, keyed by intron name. Omitted introns (absent from score_info)
    pass through unchanged.
    """
    import os
    import pandas as pd
    from intronIC.file_io.meta_sync import sync_meta_from_score_info

    df = pd.read_csv(score_path, sep="\t", dtype=str, keep_default_na=False)
    if "name" not in df.columns:
        return
    if meta_path is not None and meta_path.exists() and {"rel_score", "type_id"} <= set(df.columns):
        sync_meta_from_score_info(meta_path, df, messenger=messenger)
    if bed_path is not None and bed_path.exists() and "adjusted_score" in df.columns:
        name_to_score = dict(zip(df["name"].astype(str), df["adjusted_score"].astype(str)))
        tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".iic", dir=bed_path.parent, delete=False)
        n_synced = 0
        try:
            with open(bed_path) as f_in:
                for line in f_in:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) > 4 and parts[3] in name_to_score:   # col 3 = name, col 4 = score
                        parts[4] = name_to_score[parts[3]]
                        n_synced += 1
                    tmp.write("\t".join(parts) + "\n")
            tmp.close()
            orig_mode = bed_path.stat().st_mode if bed_path.exists() else 0o644
            shutil.move(tmp.name, bed_path)
            os.chmod(bed_path, orig_mode & 0o777)
            if messenger is not None:
                messenger.info(f"[bed-sync] rewrote score (col 4) for {n_synced} rows in {bed_path.name}")
        except Exception:
            tmp.close()
            Path(tmp.name).unlink(missing_ok=True)
            raise


def _train_and_save_raw_bundle(u12_scored, u2_scored, config, messenger):
    """Train a raw-feature SVM ensemble and save a ``scoring_mode``-stamped bundle (``raw_gated`` /
    ``pmotif_adjudicated``) — the reproducible ``main_train`` analog of ``scripts/build_raw_gated_bundle.py``.

    Drives off the already-scored reference introns (which carry ``five_raw_score``/``bp_raw_score``/
    ``three_raw_score``/``bp_offset``/``bp_scan_confidence``; ``support2_raw`` is a derived property) and
    deliberately bypasses ``ScoreNormalizer`` + ``IntronClassifier`` (both z-only). The bundle dict matches
    what the inference dispatch in ``_run_post_classification_pipeline`` expects (ensemble + scoring_mode +
    the mode's post-process params). Returns the saved model path.
    """
    import joblib
    from dataclasses import asdict
    from intronIC.classification.trainer import SVMTrainer, RAW_BASE_FEATURES
    from intronIC.classification.optimizer import SVMParameters

    mode = config.training.scoring_mode
    EXTRA = ("bp_offset", "bp_scan_confidence", "support2_raw")   # must match the inference feature order
    C = config.training.fixed_C or 200.0
    params = SVMParameters(
        C=C, calibration_method="isotonic", saturate_enabled=False, include_max=False,
        include_pairwise_mins=False, penalty="l2", class_weight_multiplier=1.0, loss="squared_hinge",
        kernel="rbf", gamma=0.001, extra_features=EXTRA, base_features=RAW_BASE_FEATURES,
    )
    trainer = SVMTrainer(
        n_models=config.training.n_models, random_state=config.training.seed, kernel="rbf",
        max_iter=config.training.max_iter, extra_feature_names=list(EXTRA),
        base_features=RAW_BASE_FEATURES,
    )
    messenger.info(
        f"Training {config.training.n_models}-model raw ensemble for scoring_mode={mode} "
        f"(C={C}, gamma=0.001, rbf, isotonic)"
    )
    ensemble = trainer.train_ensemble(
        list(u12_scored), list(u2_scored), params, subsample_u2=True, subsample_ratio=0.8,
    )

    bundle = {
        "version": f"{mode}_v1",
        "model_id": f"{mode}_{config.training.n_models}model_C{int(C)}_g0.001",
        "ensemble": ensemble,
        "scoring_mode": mode,
        "input_features": list(RAW_BASE_FEATURES) + list(EXTRA),
        "threshold": config.scoring.threshold,
        "provenance": {
            "n_u12": len(u12_scored), "n_u2": len(u2_scored),
            "n_models": config.training.n_models, "C": C,
            "note": "raw-feature bundle from main_train --scoring-mode; see docs/raw_gated_scoring.md §0c",
        },
    }
    if mode == "pmotif_adjudicated":
        from intronIC.scoring.species_adjudicator import AdjudicatorParams, ADJUDICATOR_PARAMS_VERSION
        bundle["adjudicator_params"] = asdict(AdjudicatorParams())
        bundle["provenance"]["adjudicator_params_version"] = ADJUDICATOR_PARAMS_VERSION
        bundle["provenance"]["calibration_provenance"] = (
            "DEFAULT (stamped from AdjudicatorParams; calibrated against the canonical eval_corpus raw "
            "ensemble). The Platt (margin->P_motif) and q (depth_tail->P(bearer)) constants are tied to a "
            "SPECIFIC trained ensemble's margin/depth_tail scale and do NOT transfer to a freshly-trained "
            "ensemble (verified: even same-corpus/same-n_models re-training shifts the scale). This bundle "
            "is NOT calibrated; it must be re-calibrated in eval_corpus (Platt OOF on the reference + q on "
            "the species panel) before production use."
        )
        # The frozen Platt/q constants are bound to a specific ensemble's scale (NOT reproduced by
        # re-training on the same corpus), so any main_train pmotif ensemble is uncalibrated by default.
        messenger.warning(
            "pmotif_adjudicated: this ensemble uses the DEFAULT Platt + q calibration constants, which are "
            "tied to the canonical eval_corpus ensemble and do NOT transfer to a freshly-trained ensemble. "
            "P_motif/q will be miscalibrated. Re-fit the calibration in eval_corpus (Platt OOF on the "
            "reference + q on the species panel) before using this bundle. See docs/raw_gated_scoring.md §0c."
        )
    else:  # raw_gated
        from intronIC.scoring.species_gate import SpeciesGateParams
        bundle["raw_gate_params"] = asdict(SpeciesGateParams())

    model_path = config.output.get_output_path(".model.pkl")
    joblib.dump(bundle, model_path, compress=3)
    messenger.success(
        f"Saved {mode} bundle to {model_path} "
        f"({len(ensemble.models)} models, {len(u12_scored)} U12 / {len(u2_scored)} U2)"
    )
    return model_path


@dataclass
class PostClassResult:
    """Artifacts produced by the shared post-classification pipeline.

    Returned by :func:`_run_post_classification_pipeline`. The raw-feature scoring
    modes (``raw_gated`` / ``pmotif_adjudicated``) compute their calls in-place on
    ``score_info.iic`` and only surface the called-U12 count here.
    """

    adjusted_hc_count: Optional[int] = None
    motif_category: Optional[str] = None  #: pmotif gate (DETECTED/BORDERLINE/NOT_DETECTED/UNASSESSABLE)
    z_excess: Optional[float] = None      #: pmotif population statistic (per-species)
    motif_called_u12: Optional[int] = None  #: ungated motif calls (P_motif>=0.5, disregarding motif_category)


def _run_post_classification_pipeline(
    *,
    five_z: "np.ndarray",
    bp_z: "np.ndarray",
    three_z: "np.ndarray",
    svm_s: "np.ndarray",
    type_ids: "np.ndarray",
    model_data: Any,
    config: "IntronICConfig",
    messenger: "UnifiedMessenger",
) -> PostClassResult:
    """Mode-sep (v2.6+) + continuous discount (v2.7) + unified labels (v2.7.1),
    with legacy cluster-validation as the gate-fail / non-modesep fallback.

    Shared by the streaming-classify and in-memory classify paths so the two
    routes produce equivalent output (the previous copy-paste drifted once —
    task #203). Callers pass already-filtered parallel score arrays; this rewrites
    ``adjusted_score`` / ``rel_score`` (and labels) into ``score_info.iic`` and
    ``meta.iic`` on disk, writes the ``modesep.json`` diagnostic sidecar when
    mode-sep runs, and returns the artifacts each caller needs for metrics
    assembly. An empty input (no classified introns cleared the caller's filter)
    is a no-op that returns an empty result.
    """
    result = PostClassResult()
    if len(five_z) == 0:
        return result

    score_path = config.output.get_output_path(".score_info.iic")
    meta_path = config.output.get_output_path(".meta.iic")

    # raw_gated scoring mode (transitional): replace mode-separation + the continuous discount with the
    # continuous species gate on the raw-feature classifier output. See scoring/species_gate.py and
    # eval_corpus/PROPOSED_ARCHITECTURE.md. Detected from the bundle's scoring_mode field.
    if isinstance(model_data, dict) and model_data.get("scoring_mode") == "raw_gated":
        from intronIC.scoring.species_gate import (
            apply_raw_gated_postprocess, SpeciesGateParams,
        )
        gate_params = SpeciesGateParams.from_dict(model_data.get("raw_gate_params"))
        gate_res = apply_raw_gated_postprocess(
            score_path, params=gate_params,
            threshold=config.scoring.threshold, messenger=messenger,
        )
        _sync_calls_to_meta_and_bed(
            score_path, meta_path, config.output.get_output_path(".bed.iic"), messenger=messenger)
        result.adjusted_hc_count = gate_res.n_called
        return result

    # pmotif_adjudicated scoring mode: per-intron P_motif (Platt-calibrated ensemble MARGIN) + the
    # per-species depth_tail adjudicator -> P_adj = q * P_motif. Replaces z-norm/mode-sep/discount. See
    # scoring/species_adjudicator.py and docs/raw_gated_scoring.md §0b/§0c. Runs in this shared helper so
    # streaming and in-memory are identical by construction.
    if isinstance(model_data, dict) and model_data.get("scoring_mode") == "pmotif_adjudicated":
        from intronIC.scoring.species_adjudicator import (
            apply_pmotif_adjudication, AdjudicatorParams,
        )
        adj_dict = dict(model_data.get("adjudicator_params") or {})
        if getattr(config.scoring, "adjudicator_min_u2", None) is not None:
            adj_dict["min_u2"] = config.scoring.adjudicator_min_u2   # opt-in low-N self-adjudication (option B)
        adj_params = AdjudicatorParams.from_dict(adj_dict)
        adj_res = apply_pmotif_adjudication(
            score_path, model_data["ensemble"].models, params=adj_params, messenger=messenger,
        )
        _sync_calls_to_meta_and_bed(
            score_path, meta_path, config.output.get_output_path(".bed.iic"), messenger=messenger)
        result.adjusted_hc_count = adj_res.n_adj_called
        result.motif_category = adj_res.motif_category.value
        result.z_excess = float(adj_res.z_excess) if adj_res.z_excess == adj_res.z_excess else None
        result.motif_called_u12 = adj_res.n_motif_called
        return result

    # Legacy z-normalization + mode-separation + continuous-discount scoring was removed
    # (supplant step 2). Only raw-feature bundles score (scoring_mode raw_gated /
    # pmotif_adjudicated), both handled above and returned. Reaching here means a legacy
    # zscore/modesep bundle slipped past the load-time guard (assert_scoreable_bundle).
    mode = model_data.get("scoring_mode") if isinstance(model_data, dict) else None
    raise UnsupportedBundleError(
        f"Unscoreable bundle reached post-classification (scoring_mode={mode!r}); the legacy "
        f"z-stack was removed (supplant step 2). Retrain with `--scoring-mode pmotif_adjudicated` "
        f"or use the `pre-zstack-removal` git tag for legacy bundles."
    )


def _finalize_classification_metrics(
    base_summary: dict,
    *,
    total_genes: int,
    total_introns_generated: int,
    total_scored: int,
    model_path: str,
    streaming_mode: str,
    normalizer_used: Any,
    meta_path: Optional[Path] = None,
    motif_category: Optional[str] = None,
    z_excess: Optional[float] = None,
    motif_called_u12: Optional[int] = None,
) -> dict:
    """Assemble the final per-species ``metrics.iic.json`` dict, shared by the
    streaming and in-memory classify paths so both emit equivalent metrics.

    SINGLE SOURCE OF TRUTH for all classification counts: every U12 / U2 /
    high-confidence count, percentage, and dinucleotide boundary is tallied here
    from the FINALIZED ``meta.iic`` (post-classification ``type_id`` / ``rel_score``)
    via ``count_calls_from_meta`` + ``summarize_boundaries_from_meta`` — so they
    cannot drift from the per-intron output and ``high_confidence_u12 <= u12_count``
    holds by construction. ``base_summary`` (from ``get_summary`` / ``summarize_introns``)
    supplies only the structural ``total_introns`` + ``threshold``; it carries no
    counts. The only intentional cross-mode difference is ``streaming_mode``. Caller
    must pass ``meta_path`` *after* post-classification has rewritten meta.iic.
    """
    summary = dict(base_summary)

    summary.update(
        {
            "total_genes": total_genes,
            "total_introns_generated": total_introns_generated,
            "total_scored": total_scored,
            "pretrained": True,
            "model_path": model_path,
            "streaming_mode": streaming_mode,
            "normalizer_used": normalizer_used,
        }
    )

    # ALL classification counts come from the finalized meta.iic (the authoritative
    # post-classification type_id/rel_score), via the two shared helpers — never the
    # write-time first-pass tally. u12_count + u2_count are the SCORED calls (type_id),
    # so omitted introns (type_id==NA: not_longest_isoform / short) are excluded from
    # both. high_confidence_u12 is the strong subset (rel_score > 0). Percentages are
    # over the scored set. Empty/zero if meta is unavailable.
    u12_boundaries, u2_boundaries, u12_hc_boundaries = {}, {}, {}
    n_u12 = n_hc = n_u2 = 0
    if meta_path is not None and Path(meta_path).exists():
        u12_boundaries, u2_boundaries, u12_hc_boundaries = summarize_boundaries_from_meta(
            meta_path, summary.get("threshold", 90.0)
        )
        n_u12, n_hc, n_u2 = count_calls_from_meta(meta_path)
    n_scored = n_u12 + n_u2
    summary["u12_count"] = n_u12
    summary["u2_count"] = n_u2
    summary["high_confidence_u12"] = n_hc
    summary["u12_percentage"] = (n_u12 / n_scored * 100) if n_scored else 0.0
    summary["high_confidence_percentage"] = (n_hc / n_scored * 100) if n_scored else 0.0
    summary["u12_boundaries"] = u12_boundaries
    summary["u2_boundaries"] = u2_boundaries
    # HIGH-CONFIDENCE dinucleotide breakdown (type_id==u12 AND rel_score>0): each count is a strict subset
    # of high_confidence_u12, so analyst tools can report e.g. AT-AC as a subset of the HC U12 headline
    # (the full u12_boundaries above is keyed on the wider type_id==u12 superset, P_motif>=0.5).
    summary["high_confidence_u12_boundaries"] = u12_hc_boundaries
    # PRIMARY two-number species call (pmotif_adjudicated): the motif-only gate + population statistic, so
    # consumers get the per-species determination from the metrics summary without parsing score_info.
    # DETECTED == 'a U12-motif population by motif, corroborate downstream' (snRNA cross-check is
    # database-layer); NOT_DETECTED is not a biological loss call (motif-silent divergent bearers land here).
    summary["motif_category"] = motif_category
    summary["z_excess"] = z_excess
    # UNGATED motif-call count (P_motif>=0.5, DISREGARDING motif_category). == u12_count for every
    # non-NOT_DETECTED genome; for NOT_DETECTED it's the count the species gate suppressed (u12_count -> 0),
    # so analysts can see what the motif alone called. Not derivable from the gated meta.iic tally.
    summary["motif_called_u12"] = motif_called_u12
    # The confident motif catalog: high-confidence U12 calls in a genome whose motif gate is DETECTED
    # (strong-motif calls AND motif_category==DETECTED). high_confidence_u12 itself is strong-motif calls
    # in any non-NOT_DETECTED genome (BORDERLINE calls included but flagged via motif_category).
    summary["confident_u12_motif"] = n_hc if motif_category == "DETECTED" else 0

    # Consistency guarantees (true by construction from the single finalized tally):
    # HC can never exceed the call count, and each dinucleotide breakdown is a subset of
    # its call count. These assert that no future change reintroduces a divergent path.
    assert n_hc <= n_u12, f"high_confidence_u12 ({n_hc}) > u12_count ({n_u12})"
    assert sum(u12_boundaries.values()) <= n_u12 and sum(u2_boundaries.values()) <= n_u2, (
        "dinucleotide boundaries exceed their call count"
    )
    assert sum(u12_hc_boundaries.values()) <= n_hc, (
        f"high_confidence_u12_boundaries ({sum(u12_hc_boundaries.values())}) > high_confidence_u12 ({n_hc})"
    )

    return summary


def load_pwms_with_fallback(
    config: "IntronICConfig", messenger: "UnifiedMessenger"
) -> Dict[str, "PWMSet"]:
    """
    Load PWM matrices with fallback to defaults for missing matrices.

    If a custom PWM file is provided, it is merged with the default matrices,
    allowing users to override only specific matrices without providing all of them.
    This matches v1.5.1 behavior (add_custom_matrices function).

    Args:
        config: IntronIC configuration
        messenger: Messenger for logging

    Returns:
        Dictionary mapping region to PWMSet (five, bp, three)

    Raises:
        FileNotFoundError: If PWM files don't exist
    """
    from intronIC.scoring.pwm import PWMLoader

    # Get default PWM file
    data_dir = Path(__file__).parent.parent / "data"
    default_pwm_file = data_dir / "intronIC_scoring_PWMs.json"

    # Load default matrices
    default_pwms = PWMLoader.load_from_file(
        default_pwm_file, pseudocount=config.scoring.pseudocount
    )

    # If no custom file, return defaults
    if config.scoring.pwm_file is None:
        return default_pwms

    # Load custom matrices
    custom_pwm_file = config.scoring.pwm_file
    if not custom_pwm_file.exists():
        raise FileNotFoundError(f"Custom PWM file not found: {custom_pwm_file}")

    custom_pwms = PWMLoader.load_from_file(
        custom_pwm_file, pseudocount=config.scoring.pseudocount
    )

    # Merge custom into defaults
    merge_pwm_sets(default_pwms, custom_pwms)

    # Log which matrices were customized
    custom_keys = []
    for region in ["five", "bp", "three"]:
        if region in custom_pwms:
            for key in custom_pwms[region].matrices.keys():
                custom_keys.append(f"{region}:{'-'.join(str(k) for k in key)}")

    if custom_keys:
        messenger.log_only(f"Custom PWM matrices: {', '.join(custom_keys)}")

    return default_pwms


# ============================================================================
# PHASE 3: Parallel Scoring Worker Function
# ============================================================================
# These functions must be at module level (not nested) to be picklable for multiprocessing

# Global variable to store PWMs in each worker process (set by initializer)
_worker_pwm_sets = None
_worker_scorer = None


def _init_worker(
    default_pwm_file: Path,
    custom_pwm_file: Optional[Path],
    five_coords: Tuple[int, int],
    bp_coords: Tuple[int, int],
    three_coords: Tuple[int, int],
    ignore_nc_dnts: bool,
    pseudocount: float,
):
    """
    Initialize worker process with PWMs and scorer.

    This runs once per worker process (not per task), so PWMs are loaded once
    and reused across all introns scored by this worker.

    Args:
        default_pwm_file: Path to default PWM matrices
        custom_pwm_file: Optional custom PWM file (overrides defaults)
        five_coords, bp_coords, three_coords: Scoring region coordinates
        ignore_nc_dnts: Whether to ignore non-canonical dinucleotides
        pseudocount: Pseudocount for PWM scoring
    """
    global _worker_pwm_sets, _worker_scorer

    from intronIC.scoring.pwm import PWMLoader
    from intronIC.scoring.scorer import IntronScorer

    # Load default PWMs
    pwm_sets = PWMLoader.load_from_file(default_pwm_file, pseudocount=pseudocount)

    # Merge custom PWMs if provided
    if custom_pwm_file is not None and custom_pwm_file.exists():
        custom_pwms = PWMLoader.load_from_file(custom_pwm_file, pseudocount=pseudocount)
        merge_pwm_sets(pwm_sets, custom_pwms)

    # Store in global for this worker process
    _worker_pwm_sets = pwm_sets

    # Create scorer once for this worker
    _worker_scorer = IntronScorer(
        pwm_sets=pwm_sets,
        five_coords=five_coords,
        bp_coords=bp_coords,
        three_coords=three_coords,
        ignore_nc_dnts=ignore_nc_dnts,
    )


def _init_worker_with_pwms(
    pwm_sets_dict: Dict,
    five_coords: Tuple[int, int],
    bp_coords: Tuple[int, int],
    three_coords: Tuple[int, int],
    ignore_nc_dnts: bool,
):
    """Initialize worker with pre-built PWM sets (e.g., after background correction).

    Unlike _init_worker which loads PWMs from disk, this receives PWM sets
    that have already been corrected with species-specific backgrounds.
    """
    global _worker_pwm_sets, _worker_scorer

    from intronIC.scoring.scorer import IntronScorer

    _worker_pwm_sets = pwm_sets_dict

    _worker_scorer = IntronScorer(
        pwm_sets=pwm_sets_dict,
        five_coords=five_coords,
        bp_coords=bp_coords,
        three_coords=three_coords,
        ignore_nc_dnts=ignore_nc_dnts,
    )


def _score_intron_worker_unpack(args):
    """Unpacking wrapper for imap_unordered compatibility with sequence tracking."""
    seq_idx, intron = args
    result, error = _score_intron_worker(intron)
    return seq_idx, result, error


def _score_intron_worker(intron: "Intron") -> Tuple[Optional["Intron"], Optional[str]]:
    """
    Worker function for parallel intron scoring.

    Uses PWMs and scorer initialized once per worker process (by _init_worker).
    This is much more efficient than loading PWMs for each intron.

    Args:
        intron: Intron to score

    Returns:
        (scored_intron, error_message) tuple
        If scoring succeeds: (intron, None)
        If scoring fails: (None, error_message)
    """
    global _worker_scorer

    try:
        # Use scorer initialized once for this worker process
        assert _worker_scorer is not None, "Worker scorer not initialized"
        scored = _worker_scorer.score_intron(intron)
        return (scored, None)

    except Exception as e:
        # Return error instead of crashing worker
        # Include intron ID in error message for debugging
        intron_id = intron.intron_id if hasattr(intron, "intron_id") else "unknown"
        return (None, f"{intron_id}: {str(e)}")


def log_data_block(
    logger: logging.Logger, header: str, lines: List[str], use_separator: bool = True
):
    """
    Log a block of data with a header as a single multi-line message.

    This ensures only the header gets a timestamp, while the data lines
    appear without individual timestamps for cleaner log formatting.

    Args:
        logger: Logger instance
        header: Header text (will get timestamp)
        lines: List of data lines to display
        use_separator: Whether to include separator lines (default: True)

    Example:
        log_data_block(
            logger,
            "Top 20 splice site boundaries (U12-type introns)",
            ["   1. GT-AG    11,684 (97.02%)", "   2. GC-AG       158 ( 1.31%)"]
        )

        Results in:
        [2025-11-13 00:18:48] INFO     Top 20 splice site boundaries (U12-type introns)
        ------------------------------------...
           1. GT-AG    11,684 (97.02%)
           2. GC-AG       158 ( 1.31%)
    """
    # Build the complete message starting with header
    parts = [header]

    if use_separator:
        parts.append("-" * 100)

    parts.extend(lines)

    # Log as single multi-line message (only first line gets timestamp)
    logger.info("\n".join(parts))


def clear_large_sequences_for_classification(introns: List[Intron]) -> List[Intron]:
    """
    Clear large sequence fields before classification to reduce memory.

    After scoring completes, the full intron sequence (seq) is no longer needed
    for classification. Only the small scored sequences (five_seq, three_seq,
    bp_seq, bp_seq_u2) are needed for classification.

    However, we KEEP upstream_flank, downstream_flank, and bp_region_seq because
    they are needed for meta.iic output (motif schematic and BP context).

    This function creates new intron objects with seq cleared, reducing memory
    usage by ~5-8 GB for 1M introns while preserving output quality.

    Clears:
    - seq (full intron sequence, ~500 bytes avg) - THE BIG ONE

    Keeps (for classification):
    - five_seq, three_seq, bp_seq, bp_seq_u2 (scored sequences, ~40 bytes total)
    - five_prime_dnt, three_prime_dnt (terminal dinucleotides, ~4 bytes total)

    Keeps (for meta.iic output):
    - upstream_flank, downstream_flank (needed for motif schematic, ~400 bytes total)
    - bp_region_seq (needed for BP context, ~50 bytes avg)
    - five_display_seq, three_display_seq (needed for motif schematic, ~50 bytes total)
    - bp_relative_coords (needed for BP context, ~16 bytes)

    Args:
        introns: List of scored introns with full sequences

    Returns:
        New list of introns with seq cleared (functional style)
    """
    from dataclasses import replace

    cleared = []
    for intron in introns:
        # Skip introns without sequences (pre-filtered)
        if intron.sequences is None:
            cleared.append(intron)
            continue

        # Create new intron with seq cleared (keep everything else for output)
        # Uses functional style - returns new object, original unchanged
        new_sequences = replace(
            intron.sequences,
            seq=None,  # Only clear the big one (~500 bytes avg)
        )
        cleared_intron = replace(intron, sequences=new_sequences)
        cleared.append(cleared_intron)

    return cleared


def format_count_with_percentage(count: int, total: int) -> str:
    """Format a count with percentage of total.

    Args:
        count: The count to format
        total: The total to calculate percentage from

    Returns:
        Formatted string like "1,234 (5.67%)" or "0 (0.00%)"

    Examples:
        >>> format_count_with_percentage(12074, 58933)
        '12,074 (20.49%)'
        >>> format_count_with_percentage(31, 12074)
        '31 (0.26%)'
    """
    if total == 0:
        percentage = 0.0
    else:
        percentage = (count / total) * 100
    return f"{count:,} ({percentage:.2f}%)"


def _extract_ensemble_coefficients(ensemble) -> dict:
    """Extract per-model coefficients and intercepts from trained ensemble.

    Returns a dict with feature names, per-model weights/intercepts,
    and summary statistics suitable for JSON metadata.

    For linear kernel SVC/LinearSVC: extracts coef_ and intercept_.
    For RBF kernel SVC: extracts n_support_, gamma, and support vector counts
    (no per-feature weights available for non-linear kernels).
    """
    import numpy as np

    base_features = ["s5", "sBP", "s3"]
    composite_map = {
        "min_5_bp": "min_5_bp", "min_5_3": "min_5_3", "min_bp_3": "min_bp_3",
        "min_all": "min_all", "sqrt_5_bp": "sqrt_5_bp",
        "corr_5_bp": "corr_5_bp", "corr_bp_3": "corr_bp_3", "corr_all": "corr_all",
        "gap_5_bp": "gap_5_bp", "gap_5_rest": "gap_5_rest",
        "absdiff_5_bp": "absdiff_5_bp", "absdiff_5_3": "absdiff_5_3",
        "absdiff_bp_3": "absdiff_bp_3",
        "neg_absdiff_5_bp": "neg_absdiff_5_bp", "neg_absdiff_5_3": "neg_absdiff_5_3",
        "neg_absdiff_bp_3": "neg_absdiff_bp_3", "max_5_bp": "max_5_bp", "max_5_3": "max_5_3",
    }

    all_coefs = []
    all_intercepts = []
    feature_names = None
    # For RBF kernel (no coef_ available)
    rbf_models = []

    for svm_model in ensemble.models:
        model = svm_model.model
        try:
            if hasattr(model, "calibrated_classifiers_"):
                fitted = model.calibrated_classifiers_[0].estimator
            else:
                fitted = model

            svc = None
            if hasattr(fitted, "named_steps"):
                for step in ["svc", "linearsvc", "classifier"]:
                    if step in fitted.named_steps:
                        svc = fitted.named_steps[step]
                        break

            if svc is None:
                continue

            # Check if this is a linear model with coef_ (LinearSVC or SVC(kernel='linear'))
            if hasattr(svc, "coef_"):
                coef = svc.coef_[0].tolist()
                intercept = float(svc.intercept_[0]) if hasattr(svc, "intercept_") else 0.0
                all_coefs.append(coef)
                all_intercepts.append(intercept)

                if feature_names is None:
                    transformer = fitted.named_steps.get("transform") if hasattr(fitted, "named_steps") else None
                    if transformer and hasattr(transformer, "features") and transformer.features:
                        feature_names = base_features + [
                            composite_map.get(f, f) for f in transformer.features
                        ]
                    else:
                        feature_names = base_features + [
                            f"feature_{j}" for j in range(3, len(coef))
                        ]
            else:
                # RBF or other non-linear kernel: store support vector metadata
                rbf_info = {
                    "n_support": svc.n_support_.tolist() if hasattr(svc, "n_support_") else [],
                    "total_support_vectors": int(sum(svc.n_support_)) if hasattr(svc, "n_support_") else 0,
                    "gamma": float(svc._gamma) if hasattr(svc, "_gamma") else getattr(svc, "gamma", "scale"),
                    "C": float(svc.C) if hasattr(svc, "C") else 0.0,
                    "kernel": getattr(svc, "kernel", "unknown"),
                }
                rbf_models.append(rbf_info)
        except Exception:
            continue

    # Return linear coefficient summary if available
    if all_coefs:
        coef_arr = np.array(all_coefs)
        int_arr = np.array(all_intercepts)

        return {
            "features": feature_names,
            "n_models": len(all_coefs),
            "per_model_weights": all_coefs,
            "per_model_intercepts": all_intercepts,
            "summary": {
                "mean_weights": coef_arr.mean(axis=0).tolist(),
                "std_weights": coef_arr.std(axis=0).tolist(),
                "mean_intercept": float(int_arr.mean()),
                "std_intercept": float(int_arr.std()),
                "relative_importance_pct": (
                    np.abs(coef_arr.mean(axis=0)) /
                    np.abs(coef_arr.mean(axis=0)).sum() * 100
                ).tolist(),
            },
        }

    # Return RBF model summary if available (no per-feature weights)
    if rbf_models:
        return {
            "kernel": rbf_models[0].get("kernel", "rbf"),
            "n_models": len(rbf_models),
            "per_model_support_vectors": rbf_models,
            "summary": {
                "mean_total_sv": float(np.mean([m["total_support_vectors"] for m in rbf_models])),
                "gamma": rbf_models[0].get("gamma", "scale"),
                "C": rbf_models[0].get("C", 0.0),
                "note": "Per-feature weights not available for non-linear kernels",
            },
        }

    return {}


def log_ensemble_coefficients(ensemble, messenger: "UnifiedMessenger"):
    """
    Extract and log learned coefficients from trained ensemble with detailed hyperparameters.

    Args:
        ensemble: SVMEnsemble with trained models
        messenger: UnifiedMessenger for logging
    """
    # Define feature names (matching BothEndsStrongTransformer canonical order)
    base_features = ["s5", "sBP", "s3"]
    composite_features_map = {
        "min_5_bp": "min_5_bp",
        "min_5_3": "min_5_3",
        "min_all": "min_all",
        "neg_absdiff_5_bp": "neg_absdiff_5_bp",
        "neg_absdiff_5_3": "neg_absdiff_5_3",
        "neg_absdiff_bp_3": "neg_absdiff_bp_3",
        "max_5_bp": "max_5_bp",
        "max_5_3": "max_5_3",
    }

    messenger.log_only("")
    messenger.log_only("=" * 80)
    messenger.log_only("MODEL DETAILS AND LEARNED COEFFICIENTS")
    messenger.log_only("=" * 80)

    all_coefficients = []
    all_intercepts = []
    all_feature_names = None
    model_hyperparams = []  # Store (C, penalty, loss, calibration, class_weight_mult)
    model_details = []  # Store full model details for grouped logging

    for i, svm_model in enumerate(ensemble.models):
        model = svm_model.model

        # Extract coefficients from CalibratedClassifierCV -> Pipeline -> SVC/LinearSVC
        try:
            # Navigate nested structure
            if hasattr(model, "calibrated_classifiers_"):
                calibrated_clf = model.calibrated_classifiers_[0]
                fitted_estimator = calibrated_clf.estimator

                # Get calibration method
                calib_method = getattr(model, "method", "unknown")

                # Get SVC/LinearSVC from pipeline
                linear_svc = None
                if hasattr(fitted_estimator, "named_steps"):
                    for step_name in ["svc", "linearsvc", "classifier"]:
                        if step_name in fitted_estimator.named_steps:
                            linear_svc = fitted_estimator.named_steps[step_name]
                            break

                # For RBF kernel, coef_ is not available - log support vector info instead
                if linear_svc is not None and not hasattr(linear_svc, "coef_"):
                    kernel = getattr(linear_svc, "kernel", "unknown")
                    n_sv = sum(linear_svc.n_support_) if hasattr(linear_svc, "n_support_") else "?"
                    gamma_val = getattr(linear_svc, "_gamma", getattr(linear_svc, "gamma", "?"))
                    C_val = getattr(linear_svc, "C", "?")

                    messenger.log_only(f"\nModel {i + 1}/{len(ensemble.models)}: SVC(kernel='{kernel}')")
                    messenger.log_only(f"  C={C_val}, gamma={gamma_val}, support_vectors={n_sv}")
                    messenger.log_only(f"  Calibration: {calib_method}")
                    if hasattr(linear_svc, "n_support_"):
                        messenger.log_only(f"  Per-class SVs: {linear_svc.n_support_.tolist()}")
                    continue  # Skip coef_ extraction for non-linear kernels

                if linear_svc is not None and hasattr(linear_svc, "coef_"):
                    coef = linear_svc.coef_[0]  # Shape (1, n_features)
                    intercept = (
                        linear_svc.intercept_[0]
                        if hasattr(linear_svc, "intercept_")
                        else 0.0
                    )

                    # Extract hyperparameters
                    C_param = getattr(linear_svc, "C", "unknown")
                    penalty = getattr(linear_svc, "penalty", "unknown")
                    loss = getattr(linear_svc, "loss", "unknown")

                    # Try to extract class_weight_multiplier if stored
                    class_weight_mult = "unknown"
                    if (
                        hasattr(linear_svc, "class_weight")
                        and linear_svc.class_weight is not None
                    ):
                        if (
                            isinstance(linear_svc.class_weight, dict)
                            and 1 in linear_svc.class_weight
                        ):
                            # Approximate multiplier from ratio (assuming balanced base)
                            class_weight_mult = f"~{linear_svc.class_weight[1]:.2f}"

                    model_hyperparams.append(
                        (C_param, penalty, loss, calib_method, class_weight_mult)
                    )

                    # Get feature list from transformer
                    transformer = fitted_estimator.named_steps.get("transform")
                    if (
                        transformer
                        and hasattr(transformer, "features")
                        and transformer.features is not None
                    ):
                        # Build feature names list - use transformer features directly
                        feature_names = base_features.copy()
                        for feat in transformer.features:
                            # Add features directly (transformer stores them without 'neg_' prefix)
                            # Map old naming convention to current if needed
                            if feat in composite_features_map:
                                feature_names.append(composite_features_map[feat])
                            else:
                                # Use feature name directly (handles both absdiff_X and neg_absdiff_X)
                                feature_names.append(feat)
                        all_feature_names = feature_names
                    else:
                        # Fallback: infer from coefficient count
                        n_features = len(coef)
                        if n_features == 3:
                            all_feature_names = base_features
                        else:
                            # Create generic names for composite features
                            all_feature_names = base_features + [
                                f"feature_{j}" for j in range(3, n_features)
                            ]

                    all_coefficients.append(coef)
                    all_intercepts.append(intercept)

                    # Store model details for grouped logging
                    model_details.append(
                        {
                            "index": i,
                            "hyperparams": (
                                C_param,
                                penalty,
                                loss,
                                calib_method,
                                class_weight_mult,
                            ),
                            "intercept": intercept,
                            "coefficients": coef,
                            "feature_names": all_feature_names,
                        }
                    )
        except Exception as e:
            messenger.log_only(
                f"Warning: Could not extract coefficients from model {i + 1}: {e}"
            )

    # Group models by identical hyperparameters and log details
    if model_details:
        from collections import defaultdict

        grouped = defaultdict(list)
        for detail in model_details:
            grouped[detail["hyperparams"]].append(detail)

        # Helper function to format model indices into ranges
        def format_model_indices(models):
            """Format list of model indices into compact range notation."""
            indices = sorted([m["index"] + 1 for m in models])  # 1-based for display
            if len(indices) == 1:
                return str(indices[0])

            ranges = []
            start = indices[0]
            end = indices[0]

            for i in range(1, len(indices)):
                if indices[i] == end + 1:
                    end = indices[i]
                else:
                    if start == end:
                        ranges.append(str(start))
                    else:
                        ranges.append(f"{start}-{end}")
                    start = indices[i]
                    end = indices[i]

            if start == end:
                ranges.append(str(start))
            else:
                ranges.append(f"{start}-{end}")

            return ", ".join(ranges)

        # Log each group
        for hyperparams, models in grouped.items():
            C_param, penalty, loss, calib_method, class_weight_mult = hyperparams

            # Header showing which models share these hyperparameters
            model_range = format_model_indices(models)
            messenger.log_only(
                f"\nModel(s) {model_range} (x{len(models)}) - Shared Hyperparameters:"
            )
            messenger.log_only(f"  C={C_param:.6f}, penalty={penalty}, loss={loss}")
            messenger.log_only(
                f"  Calibration={calib_method}, class_weight_mult={class_weight_mult}"
            )

            # Show individual model coefficients
            for model in models:
                messenger.log_only(
                    f"\n  Model {model['index'] + 1}/{len(ensemble.models)}:"
                )
                messenger.log_only(f"    Intercept: {model['intercept']:+.6f}")
                messenger.log_only("    Coefficients:")
                for feat_name, coef_val in zip(
                    model["feature_names"], model["coefficients"]
                ):
                    messenger.log_only(f"      {feat_name:20s}: {coef_val:+.6f}")

    # Compute and log ensemble statistics
    if all_coefficients:
        coef_array = np.array(all_coefficients)  # Shape: (n_models, n_features)
        intercept_array = np.array(all_intercepts)
        mean_coef = coef_array.mean(axis=0)
        std_coef = coef_array.std(axis=0)
        mean_intercept = intercept_array.mean()
        std_intercept = intercept_array.std()

        # Log hyperparameter summary
        messenger.log_only("\n" + "=" * 80)
        messenger.log_only("ENSEMBLE HYPERPARAMETERS SUMMARY")
        messenger.log_only("=" * 80)

        if model_hyperparams:
            # Extract unique values
            C_values = [hp[0] for hp in model_hyperparams if hp[0] != "unknown"]
            penalties = [hp[1] for hp in model_hyperparams if hp[1] != "unknown"]
            losses = [hp[2] for hp in model_hyperparams if hp[2] != "unknown"]
            calibs = [hp[3] for hp in model_hyperparams if hp[3] != "unknown"]

            if C_values:
                C_mean = np.mean(C_values)
                C_std = np.std(C_values)
                C_range = [min(C_values), max(C_values)]
                messenger.log_only(
                    f"C parameter: mean={C_mean:.6f}, std={C_std:.6f}, range={C_range}"
                )

            if penalties:
                unique_penalties = list(set(penalties))
                messenger.log_only(f"Penalty: {unique_penalties}")

            if losses:
                unique_losses = list(set(losses))
                messenger.log_only(f"Loss: {unique_losses}")

            if calibs:
                from collections import Counter

                calib_counts = Counter(calibs)
                messenger.log_only(f"Calibration methods: {dict(calib_counts)}")

        messenger.log_only("\n" + "=" * 80)
        messenger.log_only("ENSEMBLE COEFFICIENT STATISTICS")
        messenger.log_only("=" * 80)
        messenger.log_only(f"Ensemble size: {len(all_coefficients)} models")
        # all_feature_names is guaranteed to be set at this point if all_coefficients is truthy
        assert all_feature_names is not None, "Feature names not extracted from models"
        messenger.log_only(f"Feature dimensions: {len(all_feature_names)}")

        # Feature list
        if all_feature_names:
            messenger.log_only("\nFeatures used:")
            messenger.log_only(f"  Base features: {', '.join(base_features)}")
            composite_feats = [f for f in all_feature_names if f not in base_features]
            if composite_feats:
                messenger.log_only(
                    f"  Composite features: {', '.join(composite_feats)}"
                )
            else:
                messenger.log_only("  Composite features: none")

        messenger.log_only("\nIntercept Statistics:")
        messenger.log_only(f"  Mean: {mean_intercept:+.6f}")
        messenger.log_only(f"  Std:  {std_intercept:.6f}")
        messenger.log_only(
            f"  Range: [{intercept_array.min():+.6f}, {intercept_array.max():+.6f}]"
        )

        messenger.log_only("\nCoefficient Statistics:")
        messenger.log_only("-" * 80)
        messenger.log_only(
            f"{'Feature':<20s}  {'Mean':>12s}  {'Std':>12s}  {'Range':>24s}"
        )
        messenger.log_only("-" * 80)

        for feat_name, mean_val, std_val in zip(all_feature_names, mean_coef, std_coef):
            feat_min = coef_array[:, all_feature_names.index(feat_name)].min()
            feat_max = coef_array[:, all_feature_names.index(feat_name)].max()
            messenger.log_only(
                f"{feat_name:<20s}  {mean_val:+12.6f}  {std_val:12.6f}  "
                f"[{feat_min:+.6f}, {feat_max:+.6f}]"
            )

        # Feature importance ranking (by absolute mean coefficient)
        messenger.log_only("\nFeature Importance Ranking (by |mean coefficient|):")
        messenger.log_only("-" * 80)
        importance_list = [
            (name, abs(mean)) for name, mean in zip(all_feature_names, mean_coef)
        ]
        importance_list.sort(key=lambda x: x[1], reverse=True)
        for rank, (name, abs_coef) in enumerate(importance_list, 1):
            actual_coef = mean_coef[all_feature_names.index(name)]
            messenger.log_only(
                f"  {rank}. {name:20s}: {abs_coef:.6f} (actual: {actual_coef:+.6f})"
            )

        # Identify features with near-zero coefficients (L1-like sparsity)
        messenger.log_only("")
        threshold = 0.001
        near_zero = [
            (name, mean)
            for name, mean in zip(all_feature_names, mean_coef)
            if abs(mean) < threshold
        ]
        if near_zero:
            messenger.log_only(
                f"Features with near-zero coefficients (|mean| < {threshold}):"
            )
            for name, mean in near_zero:
                messenger.log_only(f"  {name}: {mean:+.6f}")
        else:
            messenger.log_only(
                f"All features have non-zero coefficients (|mean| >= {threshold})"
            )
    else:
        messenger.log_only("Warning: Could not extract coefficients from any models")

    messenger.log_only("")


def merge_scored_and_omitted_introns(
    scored_introns: List[Intron],
    all_introns: List[Intron],
    messenger: "UnifiedMessenger",
) -> List[Intron]:
    """
    Merge scored introns with unique omitted introns for complete meta output.

    Matches original intronIC behavior where .meta.iic includes both:
    - Scored introns (with classification results)
    - Omitted introns (with [o:X] tags but no scores)

    Duplicates are excluded from output (they're already filtered).

    Args:
        scored_introns: Introns that went through scoring/classification (have scores)
        all_introns: All extracted introns (including omitted)
        messenger: Unified messenger instance

    Returns:
        Combined list: scored introns + unique omitted introns (no duplicates)
    """
    # Create set of scored intron IDs for fast lookup
    scored_ids = {id(intron) for intron in scored_introns}

    # Find omitted introns that aren't duplicates
    # These should have metadata.omitted set and NOT be duplicates
    omitted_introns = [
        intron
        for intron in all_introns
        if id(intron) not in scored_ids
        and intron.metadata
        and intron.metadata.omitted != OmissionReason.NONE
        and not intron.metadata.duplicate
    ]

    total_output = len(scored_introns) + len(omitted_introns)
    total_all = len(all_introns)

    messenger.log_only(
        f"Merging output: {format_count_with_percentage(len(scored_introns), total_all)} scored + "
        f"{format_count_with_percentage(len(omitted_introns), total_all)} omitted = "
        f"{format_count_with_percentage(total_output, total_all)} total introns for output files"
    )

    # Return scored + omitted (duplicates already excluded)
    return scored_introns + omitted_introns


def calculate_minimum_intron_length(
    scoring_regions: ScoringRegions, bp_matrix_length: int
) -> int:
    """Calculate minimum intron length needed for scoring.

    Port from: intronIC.py:4600-4607

    The minimum length is determined by how much of the intron is consumed
    by the scoring regions. The key insight: we need the actual BP matrix
    length (e.g., 12bp), not the search window size (e.g., 50bp).

    Args:
        scoring_regions: Scoring region coordinates
        bp_matrix_length: Length of branch point PWM matrix

    Returns:
        Minimum intron length in bp
    """
    # Calculate intronic positions in 5' region (positions >= 0)
    # Port from: intronIC.py:4601
    five_range = range(scoring_regions.five_start, scoring_regions.five_end)
    intronic_five = len([e for e in five_range if e >= 0])

    # Calculate intronic positions in 3' region (positions < 0)
    # Port from: intronIC.py:4602
    three_range = range(scoring_regions.three_start, scoring_regions.three_end)
    intronic_three_from_score = len([e for e in three_range if e < 0])

    # BP region needs: BP_MATRIX_LENGTH + distance from 3' end
    # Port from: intronIC.py:4598, 4603-4604
    bp_margin = abs(scoring_regions.bp_end)
    intronic_three_from_bp = bp_matrix_length + bp_margin

    # Use whichever is larger
    intronic_three = max(intronic_three_from_bp, intronic_three_from_score)

    # Total minimum = positions needed at 5' end + positions needed at 3' end
    # Port from: intronIC.py:4607
    minimum_length = intronic_five + intronic_three

    return minimum_length


def setup_logging(config: IntronICConfig) -> tuple[logging.Logger, "Console"]:
    """Setup logging configuration with ANSI color support.

    Args:
        config: Pipeline configuration

    Returns:
        Tuple of (configured logger instance, Rich console for log file)
    """
    from rich.console import Console
    from rich.logging import RichHandler

    log_file = config.output.get_output_path(".iic.log")

    # Configure logging level based on flags
    # Both console and log file will use the same level for consistency
    if config.output.debug:
        log_level = logging.DEBUG
    elif config.output.quiet:
        log_level = logging.WARNING
    else:
        log_level = logging.INFO

    # Setup logger
    logger = logging.getLogger("intronIC")
    logger.setLevel(log_level)  # Logger level controls what gets through to handlers

    # Clear any existing handlers
    logger.handlers.clear()

    # Create Rich console for log file — plain text, no wrapping.
    # no_color: no ANSI escape codes, universally readable in any viewer
    # soft_wrap: no line wrapping, let the viewer handle it
    log_console = Console(
        file=open(log_file, "w", encoding="utf-8"),
        no_color=True,
        soft_wrap=True,
        highlight=False,
    )

    # File handler using Rich - preserves colors and formatting
    # NOTE: We only add a file handler. Console output is handled by UnifiedMessenger
    # to avoid duplicate messages. The logger is used ONLY for the log file.
    # File handler level matches the logger level (both controlled by --debug/--quiet)
    file_handler = RichHandler(
        console=log_console,
        show_time=True,
        show_level=True,
        show_path=False,
        markup=True,  # Enable Rich markup interpretation for proper ANSI formatting
        rich_tracebacks=True,
        tracebacks_show_locals=config.output.debug,
        level=log_level,  # Match logger level - file and console output are consistent
    )
    logger.addHandler(file_handler)

    # NO console handler - UnifiedMessenger handles console output to avoid duplication

    return logger, log_console


def load_reference_sequences(
    filepath: Path,
    max_count: Optional[int] = None,
    messenger: Optional["UnifiedMessenger"] = None,
) -> List[Intron]:
    """
    Load reference intron sequences from .iic.gz file.

    Format (tab-delimited):
    - intron_id
    - score
    - upstream_flank
    - intron_seq
    - downstream_flank
    - sources
    - source_count

    Args:
        filepath: Path to .iic.gz file
        max_count: Maximum number to load (None = all)
        messenger: Optional messenger instance

    Returns:
        List of Intron objects with sequences

    Note:
        No length filtering is applied here. Reference introns are filtered
        later using omit_check() after scoring regions are extracted, matching
        the original intronIC behavior.
    """
    introns = []

    with smart_open(filepath, "rt") as f:
        for line_num, line in enumerate(f, 1):
            # Skip comments
            if line.startswith("#"):
                continue

            # Parse line
            fields = line.strip().split("\t")
            if len(fields) < 5:
                continue

            intron_id = fields[0]
            upstream_flank = fields[2]
            intron_seq = fields[3]
            downstream_flank = fields[4]

            # Create dummy coordinates (not needed for scoring)
            coord = GenomicCoordinate(
                chromosome="ref",
                start=1,
                stop=len(intron_seq),
                strand="+",
                system="1-based",
            )

            # Extract terminal dinucleotides
            five_dnt = intron_seq[:2] if len(intron_seq) >= 2 else None
            three_dnt = intron_seq[-2:] if len(intron_seq) >= 2 else None

            # Create IntronSequences
            sequences = IntronSequences(
                seq=intron_seq,
                upstream_flank=upstream_flank,
                downstream_flank=downstream_flank,
                five_prime_dnt=five_dnt,
                three_prime_dnt=three_dnt,
            )

            # Create Intron (no length filtering - done later via omit_check)
            intron = Intron(intron_id=intron_id, coordinates=coord, sequences=sequences)
            introns.append(intron)

            # Check max count
            if max_count and len(introns) >= max_count:
                break

    if messenger:
        messenger.log_only(
            f"Loaded {len(introns)} reference sequences from {filepath.name}"
        )

    return introns


def load_genome(config: IntronICConfig, messenger: "UnifiedMessenger") -> GenomeReader:
    """Load genome file.

    Args:
        config: Pipeline configuration
        messenger: Unified messenger for output

    Returns:
        GenomeReader instance
    """
    assert config.input.genome is not None, "Genome path required"
    messenger.info(f"Loading genome: {config.input.genome}")
    # Use cached mode for faster repeated access
    reader = GenomeReader(config.input.genome, cached=True)
    if reader.cache:
        messenger.info(f"Loaded {len(reader.cache)} sequences into memory")
    return reader


def _process_contig_worker(
    contig_name: str, contig_introns: List[Intron], flank_len: int,
    u12_correction: bool, u12_correction_require_canonical: bool = True,
) -> List[Intron]:
    """
    Worker function for parallel contig processing.

    This function runs in a separate process and extracts sequences for a single
    contig using an indexed genome reader (initialized via Pool initializer).

    Args:
        contig_name: Name of the contig to process
        contig_introns: Pre-filtered introns on this contig (already filtered!)
        flank_len: Flanking sequence length
        u12_correction: Whether to apply U12 boundary corrections
        u12_correction_require_canonical: Only correct if result is canonical

    Returns:
        List of introns with sequences extracted

    Note:
        - Uses indexed FASTA for memory-efficient random access (~5-10 MB per worker)
        - No genome cache needed - each worker has lightweight file handle
        - Applies deduplication per-contig (duplicates can't span contigs)
    """
    from intronIC.extraction.sequences import SequenceExtractor
    from intronIC.file_io.indexed_genome import get_worker_genome

    try:
        # Get the worker's indexed genome reader
        indexed_genome = get_worker_genome()

        # Create extractor using indexed genome
        # Note: Using __new__ to inject custom genome reader; type ignores are needed
        extractor = SequenceExtractor.__new__(SequenceExtractor)
        extractor.genome_file = None  # type: ignore[assignment]
        extractor.genome_reader = indexed_genome  # type: ignore[assignment]
        extractor.use_cache = False  # No cache needed - using indexed access

        # Extract sequences with deduplication (per-contig deduplication)
        contig_with_seqs = list(
            extractor.extract_sequences_with_deduplication(
                contig_introns, flank_size=flank_len
            )
        )
    except Exception as e:
        import sys
        import traceback

        print(
            f"\n[ERROR] Worker failed on contig '{contig_name}': {e}", file=sys.stderr
        )
        print(f"Contig had {len(contig_introns)} introns", file=sys.stderr)
        traceback.print_exc()
        raise

    # Apply U12 corrections if enabled
    # Note: Call correct_intron_if_needed() on ALL introns - it handles the
    # noncanonical check internally via metadata.noncanonical flag.
    if u12_correction:
        from intronIC.extraction.boundary_correction import correct_intron_if_needed

        corrected_introns = []

        for intron in contig_with_seqs:
            corrected_intron, was_corrected = correct_intron_if_needed(
                intron, correction_enabled=True, use_strict_motif=True,
                require_canonical=u12_correction_require_canonical
            )
            if was_corrected:
                # Re-extract with new coordinates
                corrected_with_seq = list(
                    extractor.extract_sequences(
                        [corrected_intron], flank_size=flank_len
                    )
                )[0]
                corrected_introns.append(corrected_with_seq)
            else:
                corrected_introns.append(corrected_intron)

        contig_with_seqs = corrected_introns

    return contig_with_seqs


# Global for streaming worker initialization
_streaming_worker_db_path: str = ""
_streaming_worker_species_name: str = ""
_streaming_worker_simple_name: bool = False
_streaming_worker_no_abbreviate: bool = False
_streaming_worker_five_coords: tuple[int, int] = (0, 0)
_streaming_worker_bp_coords: tuple[int, int] = (0, 0)
_streaming_worker_three_coords: tuple[int, int] = (0, 0)
_streaming_worker_u12_correction_require_canonical: bool = True


def _init_streaming_worker(
    genome_path: str,
    db_path: str,
    species_name: str,
    simple_name: bool,
    no_abbreviate: bool,
    five_coords: tuple[int, int],
    bp_coords: tuple[int, int],
    three_coords: tuple[int, int],
    u12_correction_require_canonical: bool = True,
) -> None:
    """
    Initialize streaming worker process with genome and config.

    Called once per worker process by Pool initializer.
    Sets up both the genome reader and streaming-specific config.
    """
    from intronIC.file_io.indexed_genome import init_worker_genome

    # Initialize genome reader (reuse existing function)
    init_worker_genome(genome_path)

    # Store streaming-specific config in globals
    global _streaming_worker_db_path, _streaming_worker_species_name
    global _streaming_worker_simple_name, _streaming_worker_no_abbreviate
    global \
        _streaming_worker_five_coords, \
        _streaming_worker_bp_coords, \
        _streaming_worker_three_coords
    global _streaming_worker_u12_correction_require_canonical

    _streaming_worker_db_path = db_path
    _streaming_worker_species_name = species_name
    _streaming_worker_simple_name = simple_name
    _streaming_worker_no_abbreviate = no_abbreviate
    _streaming_worker_five_coords = five_coords
    _streaming_worker_bp_coords = bp_coords
    _streaming_worker_three_coords = three_coords
    _streaming_worker_u12_correction_require_canonical = u12_correction_require_canonical


def _process_contig_streaming_worker(
    worker_input: tuple,
) -> tuple[str, List[Intron], int, int]:
    """
    Worker function for parallel streaming contig processing.

    Like _process_contig_worker but also:
    1. Writes full sequences to SQLite (for later output)
    2. Extracts scoring motifs
    3. Returns lightweight introns (no full sequences)

    Args:
        worker_input: Tuple of (contig_name, contig_introns, flank_len,
            u12_correction, [u12_correction_require_canonical])

    Returns:
        Tuple of (contig_name, lightweight_introns, sequences_stored_count, corrections_count)
    """
    if len(worker_input) >= 5:
        contig_name, contig_introns, flank_len, u12_correction, _req_canonical = worker_input
    else:
        contig_name, contig_introns, flank_len, u12_correction = worker_input
        _req_canonical = True
    from intronIC.extraction.sequences import SequenceExtractor
    from intronIC.file_io.indexed_genome import get_worker_genome
    from intronIC.file_io.sequence_store import StreamingSequenceStore
    from intronIC.file_io.writers import generate_intron_name

    try:
        # Get the worker's indexed genome reader
        indexed_genome = get_worker_genome()

        # Create extractor using indexed genome
        extractor = SequenceExtractor.__new__(SequenceExtractor)
        extractor.genome_file = None  # type: ignore[assignment]
        extractor.genome_reader = indexed_genome  # type: ignore[assignment]
        extractor.use_cache = False

        # Extract sequences with deduplication
        contig_with_seqs = list(
            extractor.extract_sequences_with_deduplication(
                contig_introns, flank_size=flank_len
            )
        )
    except Exception as e:
        import sys
        import traceback

        print(
            f"\n[ERROR] Streaming worker failed on contig '{contig_name}': {e}",
            file=sys.stderr,
        )
        print(f"Contig had {len(contig_introns)} introns", file=sys.stderr)
        traceback.print_exc()
        raise

    # Apply U12 corrections if enabled
    corrections_count = 0
    if u12_correction:
        from intronIC.extraction.boundary_correction import correct_intron_if_needed

        corrected_introns = []

        for intron in contig_with_seqs:
            corrected_intron, was_corrected = correct_intron_if_needed(
                intron, correction_enabled=True, use_strict_motif=True,
                require_canonical=_streaming_worker_u12_correction_require_canonical
            )
            if was_corrected:
                corrected_with_seq = list(
                    extractor.extract_sequences(
                        [corrected_intron], flank_size=flank_len
                    )
                )[0]
                corrected_introns.append(corrected_with_seq)
                corrections_count += 1
            else:
                corrected_introns.append(corrected_intron)

        contig_with_seqs = corrected_introns

    # Filter introns with valid sequences for storage
    storable_introns = [i for i in contig_with_seqs if i.has_sequences]
    sequences_stored = len(storable_introns)

    # Write to SQLite (each worker opens its own connection)
    if storable_introns:
        store = StreamingSequenceStore(_streaming_worker_db_path)

        formatted_names = [
            generate_intron_name(
                intron,
                species_name=_streaming_worker_species_name,
                simple_name=_streaming_worker_simple_name,
                no_abbreviate=_streaming_worker_no_abbreviate,
            )
            for intron in storable_introns
        ]
        store.insert_batch(storable_introns, formatted_names=formatted_names)
        store.close()

    # Extract scoring motifs and create lightweight introns
    lightweight_introns = []
    for intron in contig_with_seqs:
        if intron.has_sequences:
            lightweight_intron = intron.extract_scoring_motifs(
                five_coords=_streaming_worker_five_coords,
                bp_coords=_streaming_worker_bp_coords,
                three_coords=_streaming_worker_three_coords,
            )
            lightweight_introns.append(lightweight_intron)
        else:
            lightweight_introns.append(intron)

    return contig_name, lightweight_introns, sequences_stored, corrections_count


def extract_introns_from_annotation(
    config: IntronICConfig,
    messenger: "UnifiedMessenger",
    reporter: IntronICProgressReporter,
) -> "tuple[List[Intron], int]":
    """
    Extract introns contig-by-contig with pre-filtering and deduplication.

    Returns ``(all_introns, n_genes)`` where ``n_genes`` is the total number of
    genes in the annotation hierarchy (used for metrics parity with streaming).

    This function implements the memory-optimized extraction pipeline:
    1. Parse annotations (coordinates only)
    2. Pre-filter before extraction (removes ~85-90% of introns)
    3. Process one contig at a time:
       - Extract sequences with deduplication
       - Apply U12 corrections
       - Free contig memory before next
    4. Return all introns WITH sequences for scoring

    Per-contig pipeline + prefilter keeps peak RSS bounded by the largest contig's
    intron set, not the whole genome's. Reference: ~5 GB peak on full human at -p 6.

    Note: A "contig" is a contiguous genomic sequence (chromosome, scaffold, or contig).
    This approach works for any assembly level and enables parallelization via -p flag.

    Note: Sequences will be written and cleared after scoring (done by caller).

    Args:
        config: Pipeline configuration
        genome_reader: Genome reader instance
        messenger: Unified messenger for output
        reporter: Progress reporter

    Returns:
        List of introns with sequences (ready for scoring)
    """
    messenger.info(f"Parsing annotation: {config.input.annotation}")

    # Build annotation hierarchy
    assert config.input.annotation is not None, "Annotation path required"
    builder = AnnotationHierarchyBuilder(
        child_features=["cds", "exon"],
        clean_names=config.output.clean_names,
        messenger=messenger,
    )
    genes = builder.build_from_file(str(config.input.annotation))

    # Report annotation statistics
    from intronIC.extraction.annotator import Exon, Transcript

    n_genes = len(genes)
    n_transcripts = sum(
        1 for f in builder.feature_index.values() if isinstance(f, Transcript)
    )
    n_exons = sum(1 for f in builder.feature_index.values() if isinstance(f, Exon))
    n_cds = sum(
        1
        for f in builder.feature_index.values()
        if isinstance(f, Exon)
        and f.attributes.get("_orig_feat_type", "").lower() == "cds"
    )

    messenger.info(
        f"Parsed annotation: {n_genes:,} genes, {n_transcripts:,} transcripts, "
        f"{n_cds:,} CDS, {n_exons:,} exons"
    )

    # Generate introns (coordinates only, NO sequences yet)
    messenger.log_only("Generating introns from exon pairs")
    generator = IntronGenerator(debug=config.output.debug, messenger=messenger)
    introns_iter = generator.generate_from_genes(genes, builder.feature_index)
    introns_all = list(introns_iter)
    messenger.log_only(f"Generated {len(introns_all):,} introns")

    # Report touching exons if any were found
    if generator.touching_exons_count > 0:
        messenger.log_only(
            f"Found {generator.touching_exons_count:,} touching/zero-length exon pairs "
            f"(annotation errors, counted in family size but omitted from output)"
        )

    # Filter by feature type (cds/exon/both)
    if config.extraction.feature_type == "cds":
        introns_list = [
            i
            for i in introns_all
            if i.metadata is not None and i.metadata.defined_by == "cds"
        ]
        messenger.log_only(f"Filtered to {len(introns_list):,} CDS-defined introns")
    elif config.extraction.feature_type == "exon":
        introns_list = [
            i
            for i in introns_all
            if i.metadata is not None and i.metadata.defined_by == "exon"
        ]
        messenger.log_only(f"Filtered to {len(introns_list):,} exon-defined introns")
    else:
        introns_list = introns_all

    # Free annotation hierarchy and intermediate lists (CRITICAL for memory!)
    # These can be several GB for human genome
    del introns_all
    del genes
    del builder
    gc.collect()
    messenger.log_only("Freed annotation hierarchy from memory")

    # Step 3: Pre-filter before extraction (MAJOR MEMORY SAVINGS!)
    # Pre-filtering happens silently - we'll report combined statistics after extraction
    messenger.log_only("Pre-filtering introns before sequence extraction")
    prefilter_result = prefilter_introns(
        introns=introns_list,
        min_length=config.extraction.min_intron_len,
        longest_only=not config.extraction.allow_multiple_isoforms,  # Inverse logic
        include_duplicates=config.extraction.include_duplicates,
    )

    extract_list = prefilter_result.extract_list
    skip_list = prefilter_result.skip_list

    # `or 1` guards 0/0 on intron-free genomes (e.g. trans-splicing kinetoplastids:
    # angomonas/strigomonas). The in-memory path reaches this prefilter with an empty
    # introns_list; numerators are 0 too, so the reported % is a correct 0.0. Logging only.
    _n_il = len(introns_list) or 1
    messenger.log_only(
        f"Pre-filter results: "
        f"extracting sequences for {len(extract_list):,} ({len(extract_list) / _n_il * 100:.1f}%), "
        f"skipping {len(skip_list):,} ({len(skip_list) / _n_il * 100:.1f}%)"
    )

    # Free original list
    del introns_list
    gc.collect()

    # Step 4: Extract sequences (delegate to helper for code reuse)
    all_introns = _extract_sequences_for_introns(
        extract_list, skip_list, config, messenger, reporter
    )

    # Return the gene count (n_genes, all genes in the hierarchy) alongside the
    # introns so the in-memory metrics match the streaming path's total_genes
    # (which sums len(contig_genes) per contig — i.e. all genes, incl. intron-less).
    return all_introns, n_genes


def _extract_sequences_for_introns(
    extract_list: List[Intron],
    skip_list: List[Intron],
    config: IntronICConfig,
    messenger: "UnifiedMessenger",
    reporter: IntronICProgressReporter,
) -> List[Intron]:
    """
    Helper function to extract sequences for introns.

    This is shared by both standard and streaming extraction modes.

    Args:
        extract_list: Introns that need sequences extracted
        skip_list: Introns to include without sequences
        config: Pipeline configuration
        messenger: Unified messenger
        reporter: Progress reporter

    Returns:
        All introns (extracted + skipped)
    """

    # Group introns by contig for contig-by-contig processing.
    #
    # Both the scored (extract_list) AND the omitted-but-written (skip_list)
    # introns are extracted, so omitted introns reach the output writers with
    # sequences populated (dnts, motif, bp_context, the noncanonical [n] tag).
    # This makes the standard in-memory path's meta.iic / introns.iic match the
    # streaming-classify path, which likewise feeds the full per-contig intron
    # set through extract_sequences_with_deduplication before filtering. The
    # dedup extractor shares one sequence object across coordinate-duplicate
    # introns, so this costs sequence extraction only once per unique locus.
    from collections import defaultdict

    introns_by_contig = defaultdict(list)
    for intron in (*extract_list, *skip_list):
        introns_by_contig[intron.coordinates.chromosome].append(intron)

    contigs = sorted(introns_by_contig.keys())

    # Load PWM matrices (for minimum length calculation)
    pwm_sets = load_pwms_with_fallback(config, messenger)
    bp_matrix_length = next(iter(pwm_sets["bp"].matrices.values())).length

    # Calculate minimum length (validated during pre-filtering)
    calculated_min = calculate_minimum_intron_length(
        config.scoring.scoring_regions, bp_matrix_length
    )
    # Note: actual_min_length is for logging; min_intron_len is enforced in pre-filter
    _ = max(config.extraction.min_intron_len, calculated_min)  # noqa: F841

    # Accumulator for all introns. skip_list introns are no longer appended
    # without sequences — they are grouped into introns_by_contig above and
    # extracted alongside extract_list, so every output intron has sequences.
    all_introns = []

    # Step 5: Determine parallel vs sequential mode
    n_processes = config.performance.processes
    use_parallel = n_processes > 1 and len(contigs) > 1

    # Initialize sequence extractor based on mode
    # Parallel mode: Don't load genome cache (workers use indexed FASTA)
    # Sequential mode: Load into cache for faster access
    if use_parallel:
        sequence_extractor = None  # Not used in parallel mode
    else:
        messenger.info(f"Loading genome: {config.input.genome}")
        sequence_extractor = SequenceExtractor(
            genome_file=str(config.input.genome), use_cache=True
        )
        if sequence_extractor.genome_reader.cache:
            messenger.info(
                f"Loaded {len(sequence_extractor.genome_reader.cache)} sequences into memory"
            )

    if use_parallel:
        # Parallel mode: Process multiple contigs concurrently using indexed FASTA
        messenger.info(
            f"Extracting sequences for {len(contigs)} contigs in parallel using {n_processes} processes"
        )

        # Check if index exists before creating it
        import os

        index_path = str(config.input.genome) + ".fxi"
        index_existed = os.path.exists(index_path)

        # IMPORTANT: pyfastx requires the index to be created in the main process
        # before workers can access it. Create a temporary reader to ensure index exists.
        from intronIC.file_io.indexed_genome import IndexedGenomeReader

        _temp_reader = IndexedGenomeReader(str(config.input.genome), use_cache=False)
        # Access .fasta to trigger index creation if needed
        _ = _temp_reader.fasta
        del _temp_reader  # Close and clean up

        # Now log what happened
        if index_existed:
            index_size = os.path.getsize(index_path)
            if index_size < 1024**2:
                size_str = f"{index_size / 1024:.1f} KB"
            elif index_size < 1024**3:
                size_str = f"{index_size / (1024**2):.1f} MB"
            else:
                size_str = f"{index_size / (1024**3):.1f} GB"
            messenger.info(f"Using existing genome index ({size_str})")
        else:
            # Index was just created
            index_size = os.path.getsize(index_path)
            if index_size < 1024**2:
                size_str = f"{index_size / 1024:.1f} KB"
            elif index_size < 1024**3:
                size_str = f"{index_size / (1024**2):.1f} MB"
            else:
                size_str = f"{index_size / (1024**3):.1f} GB"
            messenger.info(f"Created genome index ({size_str})")

        # Get contig lengths for length-weighted progress reporting
        from intronIC.file_io.indexed_genome import get_contig_lengths

        assert config.input.genome is not None, "Genome path required"
        contig_lengths = get_contig_lengths(config.input.genome)

        # Prepare length-weighted progress tracking
        contig_length_list = [contig_lengths[c] for c in contigs]
        cumulative_lengths = np.cumsum(contig_length_list)
        total_length = cumulative_lengths[-1]

        # Report total genome size
        if total_length < 1e6:
            size_str = f"{total_length / 1e3:.1f} Kb"
        elif total_length < 1e9:
            size_str = f"{total_length / 1e6:.1f} Mb"
        else:
            size_str = f"{total_length / 1e9:.2f} Gb"
        messenger.info(f"Total genome size: {total_length:,} bp ({size_str})")

        # Prepare inputs for worker processes (no genome cache - workers use indexed FASTA!)
        worker_inputs = [
            (
                contig,
                introns_by_contig[contig],
                config.extraction.flank_len,
                config.extraction.u12_boundary_correction,
                config.extraction.u12_correction_require_canonical,
            )
            for contig in contigs
        ]

        # Import worker initializer
        from intronIC.file_io.indexed_genome import init_worker_genome

        # Process contigs in parallel with progress tracking
        completed = 0
        total_introns_extracted = 0
        last_reported_percent = 0

        with Pool(
            processes=n_processes,
            initializer=init_worker_genome,
            initargs=(str(config.input.genome),),
        ) as pool:
            # Use starmap to process all contigs
            try:
                for contig_introns_with_seqs in pool.starmap(
                    _process_contig_worker, worker_inputs
                ):
                    completed += 1
                    introns_count = len(contig_introns_with_seqs)
                    total_introns_extracted += introns_count

                    # Debug: Check if result is None
                    if contig_introns_with_seqs is None:
                        raise RuntimeError(
                            f"Worker returned None for contig {completed}/{len(contigs)}. "
                            f"This indicates the worker function failed to return a result."
                        )

                    all_introns.extend(contig_introns_with_seqs)

                    # Log progress every 10% or on completion (length-weighted)
                    completed_length = cumulative_lengths[
                        completed - 1
                    ]  # -1 because completed is 1-indexed
                    current_percent = int((completed_length / total_length) * 100)
                    # Report when we cross a 10% boundary or complete
                    if (
                        current_percent // 10 > last_reported_percent // 10
                    ) or completed == len(contigs):
                        messenger.info(
                            f"Progress: {completed}/{len(contigs)} contigs ({current_percent}% of genome) - "
                            f"{total_introns_extracted:,} sequences extracted"
                        )
                        last_reported_percent = current_percent
            except Exception as e:
                messenger.error(
                    f"Parallel processing failed after {completed}/{len(contigs)} contigs"
                )
                messenger.error(f"Error: {e}")

                # Check if it's a pyfastx concurrency issue
                if "Fasta" in str(e) or "get seq count" in str(e):
                    messenger.error(
                        "This appears to be a pyfastx concurrency issue with this genome file"
                    )
                    messenger.error(
                        "pyfastx may not support concurrent access for some files"
                    )
                    messenger.error(
                        "Try running in sequential mode (remove -p flag) as a workaround"
                    )

                import traceback

                traceback.print_exc()
                raise

        messenger.info(
            f"Parallel processing complete: {len(contigs)} contigs processed"
        )

    else:
        # Sequential mode: Process one contig at a time (no parallelization)
        mode = (
            "sequentially"
            if not use_parallel
            else f"using 1 process (n_processes={n_processes})"
        )
        messenger.info(f"Extracting sequences for {len(contigs)} contigs {mode}")

        for contig_idx, contig in enumerate(contigs, 1):
            contig_introns = introns_by_contig[contig]
            messenger.info(
                f"[{contig_idx}/{len(contigs)}] Processing {contig}: "
                f"{len(contig_introns):,} introns"
            )

            # 5a: Extract sequences with deduplication
            assert sequence_extractor is not None, "Sequence extractor not initialized"
            contig_with_seqs = list(
                sequence_extractor.extract_sequences_with_deduplication(
                    contig_introns, flank_size=config.extraction.flank_len
                )
            )

            # 5b: Apply U12 corrections if enabled
            # Note: Call correct_intron_if_needed() on ALL introns - it handles the
            # noncanonical check internally via metadata.noncanonical flag.
            if config.extraction.u12_boundary_correction:
                from intronIC.extraction.boundary_correction import (
                    correct_intron_if_needed,
                )

                corrected_count = 0
                corrected_contig_introns = []

                for intron in contig_with_seqs:
                    corrected_intron, was_corrected = correct_intron_if_needed(
                        intron, correction_enabled=True, use_strict_motif=True,
                        require_canonical=_streaming_worker_u12_correction_require_canonical
                    )
                    if was_corrected:
                        # Re-extract with new coordinates
                        corrected_with_seq = list(
                            sequence_extractor.extract_sequences(
                                [corrected_intron],
                                flank_size=config.extraction.flank_len,
                            )
                        )[0]
                        corrected_contig_introns.append(corrected_with_seq)
                        corrected_count += 1
                    else:
                        corrected_contig_introns.append(corrected_intron)

                contig_with_seqs = corrected_contig_introns
                if corrected_count > 0:
                    messenger.log_only(
                        f"  Applied U12 corrections to {corrected_count} introns"
                    )

            # 5c: Add to accumulator (WITH sequences - needed for scoring)
            all_introns.extend(contig_with_seqs)

            # Free memory from this contig
            del contig_introns
            del contig_with_seqs
            gc.collect()

    # Free contig grouping
    del introns_by_contig
    del extract_list
    gc.collect()

    return all_introns


def extract_introns_streaming(
    config: IntronICConfig,
    messenger: "UnifiedMessenger",
    reporter: IntronICProgressReporter,
) -> Tuple[List[Intron], Path]:
    """
    Extract introns in streaming mode (roughly halves peak memory vs in-memory).

    This mode processes introns contig-by-contig and:
    1. Writes full sequences to SQLite for later output (preserving insertion order)
    2. Extracts minimal scoring motifs and clears full sequences from memory
    3. Returns lightweight introns with only motifs (for scoring)

    The SQLite database path is returned so main_classify() can read sequences
    back after scoring to write the final .introns.iic file with scores.

    Args:
        config: Pipeline configuration
        messenger: Unified messenger for output
        reporter: Progress reporter

    Returns:
        Tuple of (introns_with_motifs, db_path) where:
        - introns_with_motifs: List of Intron objects with motifs only (no full sequences)
        - db_path: Path to SQLite database containing full sequences
    """
    from collections import defaultdict

    from intronIC.extraction.filters import prefilter_introns
    from intronIC.extraction.intronator import IntronGenerator
    from intronIC.file_io.parsers import BioGLAnnotationParser
    from intronIC.file_io.sequence_store import StreamingSequenceStore

    messenger.info(f"Parsing annotation: {config.input.annotation}")
    messenger.info("Using streaming mode (roughly half the peak memory of --in-memory)")

    # Parse annotation and group by contig
    messenger.log_only("Parsing annotation and grouping by contig")
    assert config.input.annotation is not None, "Annotation path required"
    parser = BioGLAnnotationParser(clean_names=config.output.clean_names)
    annotations_by_contig = defaultdict(list)

    for ann_line in parser.parse_file(str(config.input.annotation)):
        annotations_by_contig[ann_line.region].append(ann_line)

    contigs = sorted(annotations_by_contig.keys())
    messenger.log_only(f"Found {len(contigs)} contigs")

    # Generate introns contig-by-contig for memory efficiency
    all_introns_no_seq = []
    total_genes = 0
    total_introns_generated = 0

    for contig_idx, contig in enumerate(contigs, 1):
        contig_annotations = annotations_by_contig[contig]

        if not contig_annotations or all(
            a.feat_type == "region" for a in contig_annotations
        ):
            continue

        if len(contig_annotations) > 10:
            messenger.log_only(
                f"Processing contig {contig_idx}/{len(contigs)}: {contig} ({len(contig_annotations)} features)"
            )

        builder = AnnotationHierarchyBuilder(
            child_features=["cds", "exon"],
            clean_names=config.output.clean_names,
            messenger=messenger,
        )

        try:
            contig_genes = builder.build_from_annotations(contig_annotations)
        except ValueError as e:
            if "Could not establish parent-child relationships" in str(e):
                continue
            raise

        if not contig_genes:
            continue

        generator = IntronGenerator(debug=config.output.debug, messenger=messenger)
        contig_introns = list(
            generator.generate_from_genes(contig_genes, builder.feature_index)
        )

        # Filter by feature type
        if config.extraction.feature_type == "cds":
            contig_introns = [
                i
                for i in contig_introns
                if i.metadata is not None and i.metadata.defined_by == "cds"
            ]
        elif config.extraction.feature_type == "exon":
            contig_introns = [
                i
                for i in contig_introns
                if i.metadata is not None and i.metadata.defined_by == "exon"
            ]

        total_genes += len(contig_genes)
        total_introns_generated += len(contig_introns)
        all_introns_no_seq.extend(contig_introns)

        del contig_genes
        del builder
        del contig_introns
        gc.collect()

    del annotations_by_contig
    gc.collect()

    messenger.info(
        f"Processed {total_genes:,} genes from {len(contigs)} contigs, "
        f"generated {total_introns_generated:,} introns"
    )

    # Pre-filter before extraction
    messenger.log_only("Pre-filtering introns before sequence extraction")
    prefilter_result = prefilter_introns(
        introns=all_introns_no_seq,
        min_length=config.extraction.min_intron_len,
        longest_only=not config.extraction.allow_multiple_isoforms,
        include_duplicates=config.extraction.include_duplicates,
    )

    extract_list = prefilter_result.extract_list
    skip_list = prefilter_result.skip_list

    messenger.log_only(
        f"Pre-filter results: "
        f"extracting sequences for {len(extract_list):,}, skipping {len(skip_list):,}"
    )

    del all_introns_no_seq
    gc.collect()

    # Create temporary SQLite store for sequences
    temp_dir = tempfile.mkdtemp(prefix="intronIC_sequences_")
    db_path = Path(temp_dir) / "sequences.db"
    store = StreamingSequenceStore.create(db_path)

    # Group introns by contig for extraction
    introns_by_contig = defaultdict(list)
    for intron in extract_list:
        introns_by_contig[intron.coordinates.chromosome].append(intron)

    contigs = sorted(introns_by_contig.keys())

    # Get scoring region coordinates
    five_coords = (
        config.scoring.scoring_regions.five_start,
        config.scoring.scoring_regions.five_end,
    )
    bp_coords = (
        config.scoring.scoring_regions.bp_start,
        config.scoring.scoring_regions.bp_end,
    )
    three_coords = (
        config.scoring.scoring_regions.three_start,
        config.scoring.scoring_regions.three_end,
    )

    # Get contig lengths for length-weighted progress reporting
    import numpy as np

    from intronIC.file_io.indexed_genome import get_contig_lengths

    assert config.input.genome is not None, "Genome path required"
    contig_lengths = get_contig_lengths(config.input.genome)

    # Filter out contigs missing from the genome FASTA
    missing = [c for c in contigs if c not in contig_lengths]
    if missing:
        messenger.warning(
            f"Skipping {len(missing)} contig(s) not found in genome FASTA: "
            + ", ".join(missing[:5])
            + (f" ... and {len(missing) - 5} more" if len(missing) > 5 else "")
        )
        contigs = [c for c in contigs if c in contig_lengths]

    # Prepare length-weighted progress tracking
    contig_length_list = [contig_lengths[c] for c in contigs]
    cumulative_lengths = np.cumsum(contig_length_list)
    total_length = cumulative_lengths[-1]

    # Determine parallelization
    n_processes = config.performance.processes
    use_parallel = n_processes > 1 and len(contigs) > 1

    if use_parallel:
        messenger.info(
            f"Extracting sequences in parallel ({n_processes} processes, {len(contigs)} contigs)"
        )
    else:
        messenger.info(f"Extracting sequences sequentially ({len(contigs)} contigs)")

    # Close the store we created - workers will open their own connections
    store.close()

    # Accumulators
    all_introns_with_motifs = []
    all_introns_with_motifs.extend(skip_list)  # Skipped introns (no sequences)
    total_sequences_stored = 0
    total_corrected = 0

    if use_parallel:
        # Parallel processing using Pool
        from multiprocessing import Pool

        # Prepare worker inputs and create contig index mapping
        worker_inputs = [
            (
                contig,
                introns_by_contig[contig],
                config.extraction.flank_len,
                config.extraction.u12_boundary_correction,
                config.extraction.u12_correction_require_canonical,
            )
            for contig in contigs
        ]
        contig_to_index = {contig: idx for idx, contig in enumerate(contigs)}

        # Process contigs in parallel with progress bar
        completed = 0
        completed_length = 0

        # Create progress bar
        progress = reporter.create_progress()

        with Pool(
            processes=n_processes,
            initializer=_init_streaming_worker,
            initargs=(
                str(config.input.genome),
                str(db_path),
                config.output.species_name,
                config.output.uninformative_naming,
                config.output.no_abbreviate,
                five_coords,
                bp_coords,
                three_coords,
                config.extraction.u12_correction_require_canonical,
            ),
        ) as pool, progress:
            # Add task with total = total genome length for smooth progress
            task = progress.add_task(
                "[cyan]Extracting sequences (parallel streaming)...", total=total_length
            )

            try:
                # Use imap_unordered to get results as soon as any worker completes
                # This provides better visual feedback - small contigs update progress immediately
                for (
                    contig_name,
                    lightweight_introns,
                    seqs_stored,
                    corrections,
                ) in pool.imap_unordered(
                    _process_contig_streaming_worker, worker_inputs
                ):
                    completed += 1
                    all_introns_with_motifs.extend(lightweight_introns)
                    total_sequences_stored += seqs_stored
                    total_corrected += corrections

                    # Update progress bar based on this contig's length
                    contig_idx = contig_to_index[contig_name]
                    contig_length = contig_length_list[contig_idx]
                    completed_length += contig_length
                    progress.update(task, completed=completed_length)
            except Exception as e:
                messenger.error(
                    f"Parallel streaming extraction failed after {completed}/{len(contigs)} contigs"
                )
                messenger.error(f"Error: {e}")
                raise
    else:
        # Sequential processing (original behavior)
        messenger.info(f"Loading genome: {config.input.genome}")
        sequence_extractor = SequenceExtractor(
            genome_file=str(config.input.genome), use_cache=True
        )
        if sequence_extractor.genome_reader.cache:
            messenger.info(
                f"Loaded {len(sequence_extractor.genome_reader.cache)} sequences into memory"
            )

        # Reopen store for sequential writes
        store = StreamingSequenceStore(db_path)

        for contig_idx, contig in enumerate(contigs, 1):
            contig_introns = introns_by_contig[contig]
            messenger.log_only(
                f"[{contig_idx}/{len(contigs)}] Processing {contig}: {len(contig_introns):,} introns"
            )

            # Extract sequences
            contig_with_seqs = list(
                sequence_extractor.extract_sequences_with_deduplication(
                    contig_introns, flank_size=config.extraction.flank_len
                )
            )

            # Apply U12 corrections if enabled
            if config.extraction.u12_boundary_correction:
                from intronIC.extraction.boundary_correction import (
                    correct_intron_if_needed,
                )

                corrected_count = 0
                corrected_contig_introns = []

                for intron in contig_with_seqs:
                    corrected_intron, was_corrected = correct_intron_if_needed(
                        intron, correction_enabled=True, use_strict_motif=True,
                        require_canonical=_streaming_worker_u12_correction_require_canonical
                    )
                    if was_corrected:
                        corrected_with_seq = list(
                            sequence_extractor.extract_sequences(
                                [corrected_intron],
                                flank_size=config.extraction.flank_len,
                            )
                        )[0]
                        corrected_contig_introns.append(corrected_with_seq)
                        corrected_count += 1
                    else:
                        corrected_contig_introns.append(corrected_intron)

                contig_with_seqs = corrected_contig_introns
                total_corrected += corrected_count

            # Filter and store sequences
            storable_introns = [i for i in contig_with_seqs if i.has_sequences]

            if storable_introns:
                from intronIC.file_io.writers import generate_intron_name

                formatted_names = [
                    generate_intron_name(
                        intron,
                        species_name=config.output.species_name,
                        simple_name=config.output.uninformative_naming,
                        no_abbreviate=config.output.no_abbreviate,
                    )
                    for intron in storable_introns
                ]
                store.insert_batch(storable_introns, formatted_names=formatted_names)
                total_sequences_stored += len(storable_introns)

            # Extract scoring motifs
            for intron in contig_with_seqs:
                if intron.has_sequences:
                    lightweight_intron = intron.extract_scoring_motifs(
                        five_coords=five_coords,
                        bp_coords=bp_coords,
                        three_coords=three_coords,
                    )
                    all_introns_with_motifs.append(lightweight_intron)
                else:
                    all_introns_with_motifs.append(intron)

            del contig_introns
            del contig_with_seqs
            gc.collect()

        store.close()
        del sequence_extractor
        gc.collect()

    # Free remaining resources
    del introns_by_contig
    del extract_list
    gc.collect()

    messenger.info(
        f"Streaming extraction complete: {total_sequences_stored:,} sequences stored, "
        f"{total_corrected:,} U12 corrections applied"
    )

    return all_introns_with_motifs, db_path


def extract_introns_from_bed(
    config: IntronICConfig,
    genome_reader: GenomeReader,
    messenger: "UnifiedMessenger",
    reporter: IntronICProgressReporter,
) -> List[Intron]:
    """Extract introns from BED file.

    Args:
        config: Pipeline configuration
        genome_reader: Genome reader instance
        messenger: Unified messenger for output
        reporter: Progress reporter (for progress bars only)

    Returns:
        List of extracted introns
    """
    messenger.info(f"Reading BED file: {config.input.bed}")

    # Parse BED file
    assert config.input.bed is not None, "BED file path required"
    parser = BEDParser()
    bed_lines = list(parser.parse_file(str(config.input.bed)))
    messenger.log_only(f"Parsed {len(bed_lines)} introns from BED")

    # Convert BEDLine objects to Intron objects with GenomicCoordinate
    # BED format uses 0-based half-open coordinates, convert to 1-based
    from intronIC.utils.coordinates import bed_to_internal

    introns_no_seq: List[Intron] = []
    for i, bed_line in enumerate(bed_lines):
        coord = bed_to_internal(
            chrom=bed_line.chrom,
            start=bed_line.start,
            stop=bed_line.stop,
            strand=bed_line.strand if bed_line.strand in ("+", "-") else "+",
        )
        # Use BED name if provided, otherwise just the index number
        # (output formatting adds "i" prefix, so "intron_" is redundant)
        intron_id = bed_line.name if bed_line.name != "." else str(i + 1)
        intron = Intron(intron_id=intron_id, coordinates=coord, metadata=IntronMetadata())
        introns_no_seq.append(intron)

    # Extract sequences
    messenger.info("Extracting sequences from genome")
    sequence_extractor = SequenceExtractor(
        genome_file=str(config.input.genome), use_cache=True
    )
    introns_with_seq = sequence_extractor.extract_sequences(
        introns_no_seq, flank_size=config.extraction.flank_len
    )
    # Materialize generator to list
    introns_all = list(introns_with_seq)
    messenger.log_only(f"Extracted sequences for {len(introns_all)} introns")

    # Free memory from coordinate-only list (no longer needed)
    del introns_no_seq
    gc.collect()

    # Apply U12 boundary correction (if enabled)
    if config.extraction.u12_boundary_correction:
        from intronIC.extraction.boundary_correction import correct_intron_if_needed

        messenger.info(
            "Checking non-canonical introns for U12-type boundary corrections"
        )
        corrected_count = 0
        corrected_introns = []

        for intron in introns_all:
            corrected_intron, was_corrected = correct_intron_if_needed(
                intron, correction_enabled=True, use_strict_motif=True,
                require_canonical=config.extraction.u12_correction_require_canonical
            )

            if was_corrected:
                corrected_with_seq = sequence_extractor.extract_sequences(
                    [corrected_intron], flank_size=config.extraction.flank_len
                )
                corrected_intron = list(corrected_with_seq)[0]
                corrected_count += 1

            corrected_introns.append(corrected_intron)

        introns_all = corrected_introns
        if corrected_count > 0:
            total_introns = len(introns_all)
            messenger.log_only(
                f"Applied U12-type boundary corrections to "
                f"{format_count_with_percentage(corrected_count, total_introns)} non-canonical introns"
            )

    # Load PWM matrices to get BP matrix length for minimum calculation
    # Port from: intronIC.py:4591-4592
    pwm_sets = load_pwms_with_fallback(config, messenger)
    # Get any U12 BP matrix to determine length (all should be same length)
    bp_matrix_length = next(iter(pwm_sets["bp"].matrices.values())).length
    messenger.log_only(f"BP matrix length: {bp_matrix_length}bp", level="debug")

    # Calculate actual minimum length needed for scoring regions
    # Port from: intronIC.py:4600-4617
    calculated_min = calculate_minimum_intron_length(
        config.scoring.scoring_regions, bp_matrix_length
    )
    actual_min_length = max(config.extraction.min_intron_len, calculated_min)

    messenger.log_only(
        f"Minimum intron length: {actual_min_length}bp "
        f"(user: {config.extraction.min_intron_len}bp, "
        f"scoring regions: {calculated_min}bp)"
    )

    # Keep ALL introns at extraction time (matching original behavior)
    # Let IntronFilter decide what to omit during the filtering phase
    # This ensures short introns are marked as omitted and written to output files
    # NOTE: The original keeps all introns through extraction, then marks short ones
    # as omitted during filtering (see intronIC.py lines 4772-4786)
    introns = introns_all

    messenger.log_only(
        f"Extracted {len(introns):,} introns (length filtering during scoring filter phase)"
    )

    return introns


def scored_window_slices(seq, upstream, downstream, scoring_regions):
    """Build the (five_seq, three_seq) SCORED slices from the intron + flanks,
    mirroring extraction/sequences.py. These are the EXACT regions that feed the
    5'/3' PWMs, so .score_info.iic reflects what was scored (not a display window).
    Invariant (flanks long enough): len(five_seq) == five_end - five_start;
                                    len(three_seq) == three_end - three_start.
    """
    up = upstream or ""
    down = downstream or ""
    s = seq or ""
    scoring_seq = up + s + down
    ul, dl = len(up), len(down)
    fs, fe = scoring_regions.five_start + ul, scoring_regions.five_end + ul
    ts, te = scoring_regions.three_start - dl, scoring_regions.three_end - dl
    fs = fs if fs != 0 else None
    fe = fe if fe != 0 else None
    ts = ts if ts != 0 else None
    te = te if te != 0 else None
    return scoring_seq[fs:fe], scoring_seq[ts:te]


def load_introns_from_sequences(
    config: IntronICConfig, messenger: "UnifiedMessenger"
) -> List[Intron]:
    """Load introns from pre-extracted sequences.

    Args:
        config: Pipeline configuration
        messenger: Unified messenger for output

    Returns:
        List of introns with sequences
    """
    import re as _re
    from intronIC.core.intron import (
        Intron,
        IntronFlags,
        IntronMetadata,
        IntronScores,
        IntronSequences,
        OmissionReason,
    )
    from intronIC.extraction.sequences import CANONICAL_PAIRS
    from intronIC.utils.coordinates import GenomicCoordinate

    messenger.info(f"Loading sequences: {config.input.sequence_file}")

    # Parse sequence file
    assert config.input.sequence_file is not None, "Sequence file path required"
    parser = SequenceParser()
    sequence_lines = list(parser.parse_file(str(config.input.sequence_file)))
    messenger.log_only(f"Loaded {len(sequence_lines)} sequences from file")

    # Convert SequenceLine objects to Intron objects
    introns = []
    for seq_line in sequence_lines:
        # Create minimal genomic coordinate (placeholder - we don't have real coords from sequence file)
        # These are not used for anything meaningful, just required by the Intron dataclass
        # The SEQUENCE_ONLY flag in metadata signals to skip coordinate-based operations
        coord = GenomicCoordinate(
            chromosome="",
            start=1,
            stop=len(seq_line.sequence) + 1 if seq_line.sequence else 1,
            strand="+",  # Placeholder (validation requires + or -)
            system="1-based",
        )

        # Extract terminal dinucleotides for sequences
        seq = seq_line.sequence
        five_prime_dnt = seq[:2] if seq and len(seq) >= 4 else None
        three_prime_dnt = seq[-2:] if seq and len(seq) >= 4 else None

        # Build five_seq/three_seq as the ACTUAL scored slices so .score_info.iic
        # reflects exactly what was scored — not a display 10mer. (Was: five_seq=seq[:10],
        # which dropped the exonic 5' positions -3..-1; scoring was always correct via
        # on-the-fly re-extraction, so this is an OUTPUT-only fix.)
        _five_seq, _three_seq = scored_window_slices(
            seq, seq_line.upstream_flank, seq_line.downstream_flank,
            config.scoring.scoring_regions,
        )
        sequences = IntronSequences(
            seq=seq_line.sequence,
            upstream_flank=seq_line.upstream_flank,
            downstream_flank=seq_line.downstream_flank,
            five_seq=_five_seq or (seq[:10] if seq and len(seq) >= 10 else seq),
            three_seq=_three_seq or (seq[-10:] if seq and len(seq) >= 10 else seq),
            bp_seq=None,  # Will be found during scoring
            bp_region_seq=None,  # Will be extracted during scoring
            five_prime_dnt=five_prime_dnt,
            three_prime_dnt=three_prime_dnt,
        )

        # Create IntronScores if score is present
        scores = None
        if seq_line.score is not None:
            scores = IntronScores(svm_score=seq_line.score)

        # Equivalence fix 1: recover NON-CANONICAL status so re-ingest masks the same
        # terminal dinucleotides the original run did (otherwise non-canonical termini are
        # scored against the canonical PWM → large 5'/3' raw errors). Prefer the original
        # determination preserved in the name tag (;[n]); fall back to the dinucleotides for
        # files written without it.
        _flags = IntronFlags.SEQUENCE_ONLY | IntronFlags.LONGEST_ISOFORM
        _name = seq_line.name or ""
        if _re.search(r"\[n\]", _name) or (
            "[n]" not in _name  # no tag present → derive from termini
            and five_prime_dnt and three_prime_dnt
            and (five_prime_dnt.upper(), three_prime_dnt.upper()) not in CANONICAL_PAIRS
        ):
            _flags |= IntronFlags.NONCANONICAL
        # Equivalence fix 2: honor the original omission recorded in the name tag (;[o:CODE]),
        # so re-ingest reproduces the exact scored pool that fed background correction +
        # z-normalization (the SHORT re-filter on the stored body length diverges otherwise).
        _om = _re.search(r"\[o:([a-z])\]", _name)
        _omitted = OmissionReason.from_code(_om.group(1)) if _om else OmissionReason.NONE

        # Create minimal metadata with defaults to avoid None errors in filters
        # But leave parent/grandparent as None so generate_intron_name() falls back
        # to using intron_id directly (line 256-258 in writers.py), preserving
        # the original name from the input file
        metadata = IntronMetadata(
            parent=None,
            grandparent=None,
            index=None,
            family_size=None,
            parent_length=None,
            line_number=None,
            defined_by=None,
            # Mark as sequence-only input (no real genomic coordinates)
            # This signals to skip duplicate detection and BED output
            flags=_flags,
            omitted=_omitted,
        )

        # Create Intron object
        intron = Intron(
            intron_id=seq_line.name,
            coordinates=coord,
            sequences=sequences,
            scores=scores,
            metadata=metadata,
        )
        introns.append(intron)

    messenger.log_only(f"Converted {len(introns)} sequences to Intron objects")
    return introns


def score_introns(
    introns: List[Intron],
    config: IntronICConfig,
    messenger: "UnifiedMessenger",
    reporter: IntronICProgressReporter,
    precomputed_pwm_sets: Optional[Dict] = None,
) -> List[Intron]:
    """Score introns with PWM matrices (parallel or sequential).

    PHASE 3: Restored parallel scoring from original intronIC v1.5.1
    Uses multiprocessing.Pool when config.performance.processes > 1

    Args:
        introns: List of introns to score
        config: Pipeline configuration
        messenger: Unified messenger for output
        reporter: Progress reporter (for progress bars only)
        precomputed_pwm_sets: If provided, use these PWMs instead of loading
            from disk. Used when species background correction has already
            been computed (e.g., from streaming SQLite readback).

    Returns:
        List of introns with raw scores
    """
    if precomputed_pwm_sets is not None:
        messenger.info("Using precomputed PWM matrices (species background corrected)")
        pwm_sets = precomputed_pwm_sets
    else:
        messenger.info("Loading PWM matrices")

        # Load PWM matrices with fallback support for custom matrices
        pwm_sets = load_pwms_with_fallback(config, messenger)

        # Apply species-specific U2 background correction if enabled
        # (only needed when PWMs weren't precomputed from streaming path)
        pwm_sets = apply_species_background(introns, pwm_sets, config, messenger)

    # Extract scoring configuration
    five_coords = (
        config.scoring.scoring_regions.five_start,
        config.scoring.scoring_regions.five_end,
    )
    bp_coords = (
        config.scoring.scoring_regions.bp_start,
        config.scoring.scoring_regions.bp_end,
    )
    three_coords = (
        config.scoring.scoring_regions.three_start,
        config.scoring.scoring_regions.three_end,
    )
    ignore_nc_dnts = config.scoring.ignore_nc_dnts
    n_workers = config.performance.processes

    # ========================================================================
    # PARALLEL SCORING (n_workers > 1)
    # ========================================================================
    if n_workers > 1:
        messenger.info(f"Calculating PWM scores (parallel, {n_workers} workers)")

        progress = reporter.create_progress()
        scored_introns = []
        failed_count = 0

        # Choose worker initializer: if background correction was applied,
        # pass pre-built pwm_sets directly to workers. Otherwise, have
        # workers load PWMs from disk (original behavior).
        if config.species_background.enabled:
            worker_init = _init_worker_with_pwms
            worker_initargs = (
                pwm_sets,
                five_coords,
                bp_coords,
                three_coords,
                ignore_nc_dnts,
            )
        else:
            data_dir = Path(__file__).parent.parent / "data"
            default_pwm_file = data_dir / "intronIC_scoring_PWMs.json"
            custom_pwm_file = config.scoring.pwm_file
            worker_init = _init_worker
            worker_initargs = (
                default_pwm_file,
                custom_pwm_file,
                five_coords,
                bp_coords,
                three_coords,
                ignore_nc_dnts,
                config.scoring.pseudocount,
            )

        with progress:
            task = progress.add_task("[cyan]Scoring introns...", total=len(introns))

            with Pool(
                processes=n_workers,
                initializer=worker_init,
                initargs=worker_initargs,
            ) as pool:
                try:
                    # Use imap_unordered with chunking for smooth progress updates
                    # This streams results back as they complete, not in order
                    # Chunksize: balance between overhead and memory usage
                    chunksize = max(1, min(1000, len(introns) // (n_workers * 4)))

                    # Use imap_unordered for smooth progress monitoring across all workers
                    # Add sequence indices to track original order, then sort at the end
                    # This gives us both smooth progress AND correct ordering!
                    # Workers now only need (seq_idx, intron) - PWMs loaded once per worker
                    results_iter = pool.imap_unordered(
                        _score_intron_worker_unpack,
                        zip(
                            range(len(introns)),  # Add sequence index
                            introns,
                        ),
                        chunksize=chunksize,
                    )

                    # Collect results with their sequence indices
                    # Store in dict to handle out-of-order completion
                    results_dict = {}
                    for seq_idx, scored_intron, error in results_iter:
                        if error is not None:
                            messenger.warning(f"Failed to score intron: {error}")
                            failed_count += 1
                        else:
                            results_dict[seq_idx] = scored_intron

                        progress.update(task, advance=1)

                    # Restore original order by sorting by sequence index
                    scored_introns = [
                        results_dict[i] for i in sorted(results_dict.keys())
                    ]

                except KeyboardInterrupt:
                    messenger.warning("User interrupt - terminating workers")
                    pool.terminate()
                    raise
                finally:
                    pool.close()
                    pool.join()

    # ========================================================================
    # SEQUENTIAL SCORING (n_workers == 1)
    # ========================================================================
    else:
        messenger.info("Calculating PWM scores (sequential)")

        # PWMs already loaded with custom matrix fallback at function start
        # Create scorer
        scorer = IntronScorer(
            pwm_sets=pwm_sets,
            five_coords=five_coords,
            bp_coords=bp_coords,
            three_coords=three_coords,
            ignore_nc_dnts=ignore_nc_dnts,
        )

        # Score introns sequentially with error handling
        progress = reporter.create_progress()
        scored_introns = []
        failed_count = 0

        with progress:
            task = progress.add_task("[cyan]Scoring introns...", total=len(introns))

            for intron in introns:
                try:
                    scored = scorer.score_intron(intron)
                    scored_introns.append(scored)
                except Exception as e:
                    # Log but continue - don't let one bad intron crash the pipeline
                    messenger.warning(
                        f"Failed to score intron {intron.intron_id}: {str(e)}. Skipping."
                    )
                    failed_count += 1

                progress.update(task, advance=1)

    # ========================================================================
    # COMPLETION (both modes)
    # ========================================================================
    total_attempted = len(introns)
    messenger.log_only(
        f"Scored {format_count_with_percentage(len(scored_introns), total_attempted)} introns successfully"
    )
    if failed_count > 0:
        messenger.warning(
            f"Failed to score {format_count_with_percentage(failed_count, total_attempted)} introns "
            f"(see warnings above)"
        )

    return scored_introns



def _apply_prior_adjustment(
    classified_introns: List[Intron],
    training_prior: float,
    target_prior: float,
    threshold: float,
    messenger: "UnifiedMessenger",
) -> List[Intron]:
    """
    Apply Bayesian prior adjustment to classification probabilities.

    This function adjusts probabilities via Bayes' rule to account for different
    U12 base rates between training and target species. This is particularly
    important for U12-absent species where the human prior significantly
    overestimates the true prior.

    Args:
        classified_introns: Introns with initial probabilities
        training_prior: U12 prior in training data (e.g., 0.005 for human)
        target_prior: Expected U12 prior in target species (e.g., 1e-6 for C. elegans)
        threshold: Classification threshold (for diagnostics)
        messenger: For logging diagnostics

    Returns:
        Introns with prior-adjusted probabilities

    Note:
        This is independent of margin alignment and can be applied separately or
        in combination.
    """
    from dataclasses import replace

    from intronIC.scoring.prior_adjustment import (
        adjust_probabilities_for_prior,
        compute_prior_adjustment_diagnostics,
    )

    messenger.log_only(
        f"Adjusting probabilities: training π={training_prior:.2e} → target π={target_prior:.2e}"
    )

    # Step 1: Extract probabilities (convert from 0-100 to 0-1)
    # Filter introns that have valid scores
    probs = np.array(
        [
            i.scores.svm_score / 100.0
            for i in classified_introns
            if i.scores is not None and i.scores.svm_score is not None
        ]
    )

    # Step 2: Apply prior adjustment
    probs_adj = adjust_probabilities_for_prior(
        probabilities=probs, training_prior=training_prior, target_prior=target_prior
    )

    # Step 3: Convert back to 0-100 scale
    # Ensure probs_adj is an array (adjust_probabilities_for_prior should return ndarray)
    probs_adj_array = np.asarray(probs_adj)
    probs_adj_scaled = probs_adj_array * 100.0

    # Step 4: Compute diagnostics
    diag = compute_prior_adjustment_diagnostics(
        probabilities=probs,
        adjusted_probabilities=probs_adj_array,
        training_prior=training_prior,
        target_prior=target_prior,
        threshold=threshold / 100.0,  # Convert threshold to 0-1 scale
    )

    # Step 5: Update introns with adjusted probabilities
    updated_introns = []
    for intron, new_prob in zip(classified_introns, probs_adj_scaled):
        if intron.scores is None or intron.metadata is None:
            # Skip introns without scores/metadata
            updated_introns.append(intron)
            continue
        # Update scores object
        new_scores = replace(intron.scores, svm_score=float(new_prob))
        # Reclassify based on new probability
        new_type_id = "u12" if new_prob >= threshold else "u2"
        new_metadata = replace(intron.metadata, type_id=new_type_id)
        # Update intron
        new_intron = replace(intron, scores=new_scores, metadata=new_metadata)
        updated_introns.append(new_intron)

    # Log effect
    messenger.log_only(
        f"Prior adjustment effect: {diag['n_u12_before']} → {diag['n_u12_after']} U12 predictions "
        f"({100 * diag['frac_u12_after']:.3f}%)"
    )
    messenger.log_only(
        f"Mean probability: {diag['mean_prob_before']:.4f} → {diag['mean_prob_after']:.4f}"
    )
    messenger.log_only(f"Odds ratio: {diag['odds_ratio']:.4e}")

    return updated_introns


# =============================================================================
# Parallel streaming classification worker
# =============================================================================

# Streaming-classification worker state. multiprocessing.Pool requires the
# initializer to populate a module-level name in each child interpreter; one
# global holding a dataclass satisfies that as well as N scalars, while making
# the worker's contract explicit and easy to trace.
@dataclass
class StreamingClassifyWorkerContext:
    genome_path: str
    annotation_db_path: str
    ensemble: Any
    scaler: Any
    pwm_sets: Any
    config: dict


_streaming_classify_ctx: Optional[StreamingClassifyWorkerContext] = None

# ── Parallel BG accumulation worker ──────────────────────────────────

_bg_worker_genome_path: str = ""
_bg_worker_annotation_db_path: str = ""
_bg_worker_config: dict = {}


def _init_bg_worker(genome_path: str, annotation_db_path: str, config_dict: dict) -> None:
    """Initialize BG accumulation worker process."""
    from intronIC.file_io.indexed_genome import init_worker_genome

    init_worker_genome(genome_path)

    global _bg_worker_genome_path, _bg_worker_annotation_db_path, _bg_worker_config
    _bg_worker_genome_path = genome_path
    _bg_worker_annotation_db_path = annotation_db_path
    _bg_worker_config = config_dict


def _process_contigs_bg_worker(
    contig_batch: list,
) -> tuple:
    """Accumulate BG nucleotide frequencies for a batch of contigs.

    Calls the shared :func:`_streaming_extract_and_filter_contig` helper
    so the BG accumulator sees the same intron set that the adaptive-fit
    pre-pass and the main classify pass will see (and that the in-memory
    pipeline accumulates over). This is required for streaming and
    in-memory to produce identical classifications.

    Args:
        contig_batch: List of (contig_name, count) tuples

    Returns:
        (five_counts, three_counts, bp_counts, n_accumulated)
        where *_counts are dicts suitable for merge_counts().
    """
    from intronIC.file_io.indexed_genome import get_worker_genome
    from intronIC.scoring.background import _RegionAccumulator, _BPSAccumulator

    config_dict = _bg_worker_config
    five_len = config_dict['five_len']
    three_len = config_dict['three_len']
    five_start = config_dict['five_start']
    five_end = config_dict['five_end']
    three_start = config_dict['three_start']
    three_end = config_dict['three_end']
    bp_start = config_dict['bp_start']
    bp_end = config_dict['bp_end']
    include_duplicates = config_dict.get('include_duplicates', False)

    five_acc = _RegionAccumulator(five_len)
    three_acc = _RegionAccumulator(three_len)
    bp_acc = _BPSAccumulator()
    n_accumulated = 0

    indexed_genome = get_worker_genome()

    for contig_name, _count in contig_batch:
        _, scorable, _ = _streaming_extract_and_filter_contig(
            contig=contig_name,
            annotation_db_path=_bg_worker_annotation_db_path,
            config=config_dict,
            indexed_genome=indexed_genome,
        )

        for intron in scorable:
            seqs = intron.sequences
            if seqs is None or seqs.seq is None:
                continue
            if (
                not include_duplicates
                and intron.metadata
                and intron.metadata.duplicate
            ):
                # Belt-and-suspenders: helper already drops duplicates,
                # but mirror in-memory's apply_species_background guard.
                continue

            five_dnt = seqs.five_prime_dnt or (seqs.seq[:2] if seqs.seq else "")
            three_dnt = seqs.three_prime_dnt or (seqs.seq[-2:] if seqs.seq else "")
            if not five_dnt or not three_dnt:
                continue

            upstream = seqs.upstream_flank or ""
            seq = seqs.seq or ""
            downstream = seqs.downstream_flank or ""

            five_seq = (upstream[five_start:] + seq[:five_end]).upper()
            three_seq = (seq[three_start:] + downstream[:three_end]).upper()
            bp_region = seq[
                max(0, len(seq) + bp_start) : len(seq) + bp_end
            ].upper()

            five_acc.add(five_dnt.upper(), five_seq)
            three_acc.add(three_dnt.upper(), three_seq)
            bp_acc.add(five_dnt.upper(), bp_region)
            n_accumulated += 1

    return (five_acc.export_counts(), three_acc.export_counts(),
            bp_acc.export_counts(), n_accumulated)


def _streaming_extract_and_filter_contig(
    contig: str,
    annotation_db_path: str,
    config: dict,
    indexed_genome: Any,
) -> tuple[List[Intron], List[Intron], dict]:
    """Extract, sequence, U12-correct, and filter introns for one contig.

    Shared by every per-contig pass in the streaming-classify pipeline:
    BG correction, adaptive normalizer fit, and the main classify pass.
    Centralising it guarantees all three see the same intron set, which
    is required for streaming results to match the in-memory pipeline.

    Returns:
        (filtered_introns, scorable, stats)
        - filtered_introns: ALL introns including omitted (for output).
        - scorable: subset that passed filtering, has sequences, and is
          eligible for scoring (``omitted == NONE`` and not a
          coordinate-duplicate when ``include_duplicates`` is False).
        - stats: dict with the standard streaming stats keys.
    """
    from intronIC.extraction.annotator import AnnotationHierarchyBuilder
    from intronIC.extraction.filters import FilterStats, IntronFilter
    from intronIC.extraction.intronator import IntronGenerator
    from intronIC.extraction.sequences import SequenceExtractor
    from intronIC.file_io.annotation_store import StreamingAnnotationStore

    stats: dict = {
        "genes": 0,
        "introns_generated": 0,
        "scored": 0,
        "filter_stats": FilterStats(),
        "duplicate_map": {},
        "overlap_map": {},
    }

    annotation_store = StreamingAnnotationStore(annotation_db_path)
    contig_annotations = annotation_store.get_annotations_for_contig(contig)
    annotation_store.close()

    if not contig_annotations or all(
        a.feat_type == "region" for a in contig_annotations
    ):
        return [], [], stats

    builder = AnnotationHierarchyBuilder(
        child_features=["cds", "exon"],
        clean_names=config["clean_names"],
        messenger=None,
    )
    try:
        contig_genes = builder.build_from_annotations(contig_annotations)
    except ValueError as e:
        if "Could not establish parent-child relationships" in str(e):
            return [], [], stats
        raise
    del contig_annotations

    if not contig_genes:
        return [], [], stats

    stats["genes"] = len(contig_genes)

    generator = IntronGenerator(debug=config["debug"], messenger=None)
    contig_introns = list(
        generator.generate_from_genes(contig_genes, builder.feature_index)
    )

    if config["feature_type"] == "cds":
        contig_introns = [
            i for i in contig_introns
            if i.metadata is not None and i.metadata.defined_by == "cds"
        ]
    elif config["feature_type"] == "exon":
        contig_introns = [
            i for i in contig_introns
            if i.metadata is not None and i.metadata.defined_by == "exon"
        ]

    stats["introns_generated"] = len(contig_introns)

    if not contig_introns:
        return [], [], stats

    extractor = SequenceExtractor.__new__(SequenceExtractor)
    extractor.genome_file = None  # type: ignore[assignment]
    extractor.genome_reader = indexed_genome  # type: ignore[assignment]
    extractor.use_cache = False

    contig_with_seqs = list(
        extractor.extract_sequences_with_deduplication(
            contig_introns, flank_size=config["flank_len"]
        )
    )

    if config["u12_boundary_correction"]:
        from intronIC.extraction.boundary_correction import correct_intron_if_needed

        corrected_contig_introns = []
        for intron in contig_with_seqs:
            corrected_intron, was_corrected = correct_intron_if_needed(
                intron, correction_enabled=True, use_strict_motif=True,
                require_canonical=config.get("u12_correction_require_canonical", True),
            )
            if was_corrected:
                corrected_with_seq = list(
                    extractor.extract_sequences(
                        [corrected_intron], flank_size=config["flank_len"]
                    )
                )[0]
                corrected_contig_introns.append(corrected_with_seq)
            else:
                corrected_contig_introns.append(corrected_intron)
        contig_with_seqs = corrected_contig_introns

    intron_filter = IntronFilter(
        min_length=config["min_intron_len"],
        bp_matrix_length=7,
        scoring_regions=["five", "three"],
        allow_noncanonical=not config["exclude_noncanonical"],
        allow_overlap=not config["no_intron_overlap"],
        longest_only=True,
        include_duplicates=config["include_duplicates"],
    )
    filtered_introns = intron_filter.filter_introns(contig_with_seqs)

    # Match the in-memory standard pipeline: prefilter_introns drops
    # coordinate-duplicates before scoring (and the pre-write filter drops
    # them again before output). IntronFilter only TAGS duplicates here,
    # so we drop them explicitly. Without this, streaming would score and
    # emit each isoform of a shared intron, while in-memory writes each
    # unique locus once — producing different counts on identical input.
    if not config["include_duplicates"]:
        filtered_introns = [
            i for i in filtered_introns
            if i.metadata is None or not i.metadata.duplicate
        ]

    stats["filter_stats"] = intron_filter.stats
    stats["duplicate_map"] = intron_filter.get_duplicate_map()
    stats["overlap_map"] = intron_filter.get_overlap_map()

    scorable = [
        i for i in filtered_introns
        if i.has_sequences
        and (i.metadata is None or i.metadata.omitted == OmissionReason.NONE)
    ]

    return filtered_introns, scorable, stats


def _streaming_extract_and_score_contig(
    contig: str,
    annotation_db_path: str,
    pwm_sets: Any,
    config: dict,
    indexed_genome: Any,
) -> tuple[List[Intron], List[Intron], dict]:
    """Extract, filter, and PWM-score introns for one contig.

    Thin wrapper around :func:`_streaming_extract_and_filter_contig` that
    additionally scores the scorable introns with the given (possibly
    BG-corrected) PWMs. Returns introns with raw scores populated; the
    caller handles normalization and classification.
    """
    from intronIC.scoring.scorer import IntronScorer

    filtered_introns, scorable, stats = _streaming_extract_and_filter_contig(
        contig=contig,
        annotation_db_path=annotation_db_path,
        config=config,
        indexed_genome=indexed_genome,
    )

    if not scorable:
        return filtered_introns, [], stats

    scorer = IntronScorer(
        pwm_sets=pwm_sets,
        five_coords=(config["five_start"], config["five_end"]),
        bp_coords=(config["bp_start"], config["bp_end"]),
        three_coords=(config["three_start"], config["three_end"]),
        ignore_nc_dnts=config["ignore_nc_dnts"],
    )

    scored_scorable = [scorer.score_intron(i) for i in scorable]
    stats["scored"] = len(scored_scorable)

    return filtered_introns, scored_scorable, stats

def _init_streaming_classify_worker(
    genome_path: str,
    annotation_db_path: str,
    ensemble: Any,
    scaler: Any,
    pwm_sets: Any,
    config_dict: dict,
) -> None:
    """
    Initialize streaming classification worker process.

    Called once per worker process by Pool initializer.
    Sets up genome reader and stores classification components.
    """
    from intronIC.file_io.indexed_genome import init_worker_genome

    # Initialize genome reader
    init_worker_genome(genome_path)

    # Store classification components in the single worker-context global.
    global _streaming_classify_ctx
    _streaming_classify_ctx = StreamingClassifyWorkerContext(
        genome_path=genome_path,
        annotation_db_path=annotation_db_path,
        ensemble=ensemble,
        scaler=scaler,
        pwm_sets=pwm_sets,
        config=config_dict,
    )


def _process_contig_streaming_classify_worker(
    contig_input: tuple[str, int],
) -> tuple[str, List[Intron], List[Intron], dict]:
    """
    Worker function for parallel streaming classification.

    Processes a single contig: extracts sequences, scores, normalizes, and classifies.
    Returns contig name, classified introns, all filtered introns, and statistics.

    Args:
        contig_input: Tuple of (contig_name, contig_annotation_count)

    Returns:
        Tuple of (contig_name, classified_introns, filtered_introns, stats_dict)
        - classified_introns: Introns that were scored and classified
        - filtered_introns: ALL introns including omitted ones (for output files)
    """
    contig, _count = contig_input
    from intronIC.classification.predictor import classify_introns_batch
    from intronIC.file_io.indexed_genome import get_worker_genome

    # Access worker context
    ctx = _streaming_classify_ctx
    config = ctx.config
    ensemble = ctx.ensemble
    pwm_sets = ctx.pwm_sets

    indexed_genome = get_worker_genome()
    filtered_introns, scored_scorable, base_stats = _streaming_extract_and_score_contig(
        contig=contig,
        annotation_db_path=ctx.annotation_db_path,
        pwm_sets=pwm_sets,
        config=config,
        indexed_genome=indexed_genome,
    )

    stats = {
        **base_stats,
        "classified": 0,
    }

    if not scored_scorable:
        return contig, [], filtered_introns, stats

    # Raw-feature bundles classify directly from the background-corrected raw motif scores
    # (z-normalization removed — supplant 2b).
    classified_introns = classify_introns_batch(
        scored_scorable,
        ensemble,
        threshold=config["threshold"],
    )
    stats["classified"] = len(classified_introns)

    # Dinucleotide boundary tallies are accumulated downstream by the
    # StreamingOutputWriter (deduplicated, type_id-keyed) — no per-contig
    # pre-dedup counters needed here.
    return contig, classified_introns, filtered_introns, stats


def classify_streaming_per_contig(
    config: IntronICConfig,
    messenger: "UnifiedMessenger",
    reporter: IntronICProgressReporter,
) -> Tuple[int, dict]:
    """True streaming classification - processes one contig at a time.

    This function provides ~90% memory savings compared to standard mode by:
    1. Processing introns per-contig instead of accumulating all in memory
    2. Using pre-trained model's frozen scaler (no fitting needed)
    3. Writing outputs immediately after classification
    4. Freeing memory after each contig

    Peak memory: O(largest_chromosome) instead of O(all_introns)

    REQUIREMENTS:
    - Only works with pre-trained models (--pretrained-model required)
    - Only works in annotation mode (not BED or sequences input)

    Args:
        config: Pipeline configuration
        messenger: Unified messenger for console and log output
        reporter: Progress reporter

    Returns:
        Tuple of (total_classified, classification_summary_dict)

    Raises:
        ValueError: If pretrained model not specified
        ValueError: If not in annotation mode
    """
    from intronIC.file_io.annotation_store import StreamingAnnotationStore
    from intronIC.file_io.indexed_genome import IndexedGenomeReader
    from intronIC.file_io.parsers import BioGLAnnotationParser
    from intronIC.file_io.writers import MappingWriter, StreamingOutputWriter

    # Validate requirements
    if config.training.pretrained_model_path is None:
        raise ValueError(
            "True streaming mode requires a pre-trained model. "
            "Use --pretrained-model to specify the model file."
        )

    if config.input.mode != "annotation":
        raise ValueError(
            "True streaming mode currently only supports annotation input mode. "
            "Use standard mode for BED or sequence inputs."
        )

    messenger.info("Streaming mode: processing per-contig")

    # Load pre-trained model and extract components
    messenger.info(
        f"Loading pretrained model from {config.training.pretrained_model_path}"
    )
    model_data = normalize_model_bundle(
        load_model(config.training.pretrained_model_path)
    )
    assert_scoreable_bundle(model_data)  # fail fast on legacy z-stack bundles (supplant step 2)

    # Handle both old and new model format
    if isinstance(model_data, dict):
        ensemble = model_data["ensemble"]
        saved_normalizer = model_data.get("normalizer", None)
    else:
        ensemble = model_data
        saved_normalizer = None

    messenger.log_only(f"Loaded ensemble with {len(ensemble.models)} models")

    # Raw-feature bundles (the only scoreable kind after the z-stack removal) classify
    # directly from the background-corrected raw motif scores — no per-species
    # z-normalization, so the per-contig workers need no scaler. The z-only
    # --load-normalizer / --save-normalizer options no longer apply.
    if config.scoring.load_normalizer is not None or config.scoring.save_normalizer:
        messenger.warning(
            "--load-normalizer / --save-normalizer are ignored: z-normalization was removed "
            "(supplant 2b); raw-feature bundles score directly from raw motif scores."
        )
    scaler = None
    scaler_source = "none (raw features)"

    # Load PWM matrices
    messenger.info("Loading PWM matrices")
    pwm_sets = load_pwms_with_fallback(config, messenger)

    # =========================================================================
    # MEMORY OPTIMIZATION 1: Parse annotations into SQLite (single pass)
    # This avoids loading all annotations into memory at once
    # =========================================================================
    messenger.info(f"Indexing annotation: {config.input.annotation}")
    assert config.input.annotation is not None, "Annotation path required"

    # Create temporary SQLite store for annotations
    temp_dir = tempfile.mkdtemp(prefix="intronIC_annotations_")
    annotation_db_path = Path(temp_dir) / "annotations.db"

    parser = BioGLAnnotationParser(clean_names=config.output.clean_names)
    annotation_store = StreamingAnnotationStore.create_from_file(
        annotation_path=config.input.annotation,
        db_path=annotation_db_path,
        parser=parser,
    )

    contigs_with_counts = annotation_store.get_contigs_with_counts()
    total_annotations = annotation_store.get_total_annotations()
    messenger.log_only(
        f"Indexed {total_annotations:,} annotations across {len(contigs_with_counts)} contigs"
    )

    # =========================================================================
    # MEMORY OPTIMIZATION 2: Indexed genome access
    # Use pyfastx for efficient random access without loading entire genome
    # =========================================================================
    messenger.info(f"Using indexed genome access: {config.input.genome}")
    genome_reader = IndexedGenomeReader(str(config.input.genome), use_cache=False)

    # Force the pyfastx .fxi index to be built single-threaded in the main
    # process before any worker pool starts. Without this, workers in the
    # BG-correction phase (or the adaptive-fit phase, or the classify phase)
    # all race to build the same .fxi file, leaving it 0-byte and
    # corrupting later access. get_contig_lengths is the lightest call
    # that materialises the index — we need contig lengths later for
    # progress tracking anyway, so capture them once here and reuse.
    from intronIC.file_io.indexed_genome import get_contig_lengths
    assert config.input.genome is not None, "Genome path required"
    contig_lengths = get_contig_lengths(config.input.genome)
    messenger.log_only(
        f"Built pyfastx index covering {len(contig_lengths):,} contigs"
    )

    # Species-specific U2 background correction in streaming classify mode.
    # Lightweight first-pass: parse annotations → generate introns → fetch
    # only the short scored motif windows (~83 bytes/intron) → accumulate
    # nucleotide frequencies. Parallelized across contigs when -p > 1.
    if config.species_background.enabled:
        bg, u12_sets_bg, _, five_len, three_len = _create_background_accumulator(
            config, pwm_sets
        )

        # Build a single config dict with everything the shared
        # extract+filter helper needs, plus the BG-specific lengths.
        bg_config_dict = {
            'five_start': config.scoring.scoring_regions.five_start,
            'five_end': config.scoring.scoring_regions.five_end,
            'three_start': config.scoring.scoring_regions.three_start,
            'three_end': config.scoring.scoring_regions.three_end,
            'bp_start': config.scoring.scoring_regions.bp_start,
            'bp_end': config.scoring.scoring_regions.bp_end,
            'five_len': five_len,
            'three_len': three_len,
            'clean_names': config.output.clean_names,
            'debug': config.output.debug,
            'include_duplicates': config.extraction.include_duplicates,
            'min_intron_len': config.extraction.min_intron_len,
            # Required by _streaming_extract_and_filter_contig:
            'feature_type': config.extraction.feature_type,
            'flank_len': config.extraction.flank_len,
            'u12_boundary_correction': (
                config.extraction.u12_boundary_correction
            ),
            'u12_correction_require_canonical': (
                config.extraction.u12_correction_require_canonical
            ),
            'exclude_noncanonical': config.scoring.exclude_noncanonical,
            'no_intron_overlap': config.extraction.no_intron_overlap,
        }

        n_bg_introns = 0
        n_bg_workers = config.performance.processes

        # Both paths use the same worker function (_process_contigs_bg_worker)
        # to guarantee identical BG correction regardless of -p value.
        # The worker creates its own genome reader, avoiding pyfastx state
        # corruption from many sequential fetch() calls on a single handle.
        bg_init_args = (
            str(config.input.genome),
            str(annotation_store.db_path),
            bg_config_dict,
        )

        if n_bg_workers > 1:
            messenger.info(
                f"Species background: accumulating frequencies "
                f"(parallel, {n_bg_workers} workers)..."
            )

            batch_size = max(1, len(contigs_with_counts) // (n_bg_workers * 2))
            contig_batches = [
                contigs_with_counts[i:i + batch_size]
                for i in range(0, len(contigs_with_counts), batch_size)
            ]

            with Pool(
                processes=n_bg_workers,
                initializer=_init_bg_worker,
                initargs=bg_init_args,
            ) as pool:
                for five_counts, three_counts, bp_counts, n_acc in pool.imap_unordered(
                    _process_contigs_bg_worker, contig_batches
                ):
                    bg.merge_worker_counts(five_counts, three_counts, bp_counts)
                    n_bg_introns += n_acc

        else:
            messenger.info("Species background: accumulating frequencies across contigs...")
            _init_bg_worker(*bg_init_args)
            five_counts, three_counts, bp_counts, n_acc = (
                _process_contigs_bg_worker(contigs_with_counts)
            )
            bg.merge_worker_counts(five_counts, three_counts, bp_counts)
            n_bg_introns = n_acc

        bg_config = config.species_background
        corrected = _finalize_background_correction(
            bg, u12_sets_bg, n_bg_introns, bg_config.min_introns,
            bg_config.n0, messenger,
        )
        if corrected is not None:
            pwm_sets = corrected

    # Initialize output writer
    output_writer = StreamingOutputWriter(
        output_dir=config.output.output_dir,
        base_name=config.output.base_filename,
        species_name=config.output.species_name,
        simple_name=config.output.uninformative_naming,
        no_abbreviate=config.output.no_abbreviate,
        write_bed=True,
        write_sequences=True,
        write_scores=True,
        no_headers=config.output.no_headers,
    )
    output_writer.threshold = config.scoring.threshold

    # Statistics tracking
    total_genes = 0
    total_introns_generated = 0
    total_scored = 0
    total_classified = 0

    # Lightweight score accumulation for post-hoc cluster validation
    # Only stores 3 floats per classified intron (~24 bytes each)
    accumulated_five_z: list = []
    accumulated_bp_z: list = []
    accumulated_three_z: list = []
    accumulated_svm_scores: list = []
    accumulated_type_ids: list = []

    # Accumulated filter statistics
    from intronIC.extraction.filters import FilterStats

    accumulated_filter_stats = FilterStats()

    # Accumulated duplicate and overlap maps
    accumulated_duplicate_map: Dict[str, Set[str]] = {}
    accumulated_overlap_map: Dict[str, Set[str]] = {}

    # Use contig lengths captured up-front (used here for length-weighted
    # progress reporting and for filtering out missing contigs).
    import numpy as np

    # Filter out annotation contigs missing from the genome FASTA
    # (e.g. organellar genes referencing contigs not in nuclear assemblies)
    missing_contigs = [
        (contig, count) for contig, count in contigs_with_counts
        if contig not in contig_lengths
    ]
    if missing_contigs:
        total_missing_annotations = sum(c for _, c in missing_contigs)
        missing_names = [c for c, _ in missing_contigs]
        messenger.warning(
            f"Skipping {len(missing_contigs)} contig(s) not found in genome FASTA "
            f"({total_missing_annotations:,} annotations): "
            + ", ".join(missing_names[:5])
            + (f" ... and {len(missing_names) - 5} more" if len(missing_names) > 5 else "")
        )
        contigs_with_counts = [
            (contig, count) for contig, count in contigs_with_counts
            if contig in contig_lengths
        ]

    # Prepare length-weighted progress tracking
    contig_names = [contig for contig, _ in contigs_with_counts]
    contig_length_list = [contig_lengths[c] for c in contig_names]
    cumulative_lengths = np.cumsum(contig_length_list)
    total_length = cumulative_lengths[-1]

    # Determine parallelization strategy.
    # Both paths use the same per-contig worker function
    # (_process_contig_streaming_classify_worker) to guarantee identical
    # classification logic regardless of -p value.
    n_processes = config.performance.processes
    use_parallel = n_processes > 1 and len(contigs_with_counts) > 1

    # Build config dict consumed by the worker function
    worker_config = {
        "clean_names": config.output.clean_names,
        "debug": config.output.debug,
        "feature_type": config.extraction.feature_type,
        "flank_len": config.extraction.flank_len,
        "u12_boundary_correction": config.extraction.u12_boundary_correction,
        "u12_correction_require_canonical": config.extraction.u12_correction_require_canonical,
        "min_intron_len": config.extraction.min_intron_len,
        "exclude_noncanonical": config.scoring.exclude_noncanonical,
        "no_intron_overlap": config.extraction.no_intron_overlap,
        "include_duplicates": config.extraction.include_duplicates,
        "threshold": config.scoring.threshold,
        "ignore_nc_dnts": config.scoring.ignore_nc_dnts,
        "five_start": config.scoring.scoring_regions.five_start,
        "five_end": config.scoring.scoring_regions.five_end,
        "bp_start": config.scoring.scoring_regions.bp_start,
        "bp_end": config.scoring.scoring_regions.bp_end,
        "three_start": config.scoring.scoring_regions.three_start,
        "three_end": config.scoring.scoring_regions.three_end,
    }

    # Close annotation store in main process — workers (or direct calls in
    # sequential mode) open their own read-only connections.
    annotation_store.close()

    # Free main-process genome reader; workers create their own
    del genome_reader
    gc.collect()

    # Prepare worker inputs and contig index mapping
    worker_inputs = [(contig, count) for contig, count in contigs_with_counts]
    contig_to_index = {
        contig: idx for idx, (contig, _) in enumerate(contigs_with_counts)
    }

    worker_init_args = (
        str(config.input.genome),
        str(annotation_db_path),
        ensemble,
        scaler,
        pwm_sets,
        worker_config,
    )

    # Set up dispatch: Pool for parallel, direct calls for sequential
    from contextlib import nullcontext

    if use_parallel:
        messenger.info(
            f"Processing {len(contigs_with_counts)} contigs in parallel "
            f"({n_processes} processes)"
        )
        pool_mgr = Pool(
            processes=n_processes,
            initializer=_init_streaming_classify_worker,
            initargs=worker_init_args,
        )
    else:
        messenger.info(f"Processing {len(contigs_with_counts)} contigs sequentially")
        _init_streaming_classify_worker(*worker_init_args)
        pool_mgr = nullcontext()

    # Process contigs and accumulate results
    all_classified_introns = []
    all_filtered_introns = []
    completed_length = 0

    progress = reporter.create_progress()

    with pool_mgr as pool, progress:
        task = progress.add_task(
            "[cyan]Classifying introns...", total=total_length
        )

        if use_parallel:
            contig_iter = pool.imap_unordered(
                _process_contig_streaming_classify_worker, worker_inputs
            )
        else:
            contig_iter = (
                _process_contig_streaming_classify_worker(ci)
                for ci in worker_inputs
            )

        try:
            for contig_name, classified_introns, filtered_introns, stats in contig_iter:
                all_classified_introns.extend(classified_introns)
                all_filtered_introns.extend(filtered_introns)

                # Accumulate statistics
                total_genes += stats["genes"]
                total_introns_generated += stats["introns_generated"]
                total_scored += stats["scored"]
                total_classified += stats["classified"]

                # Accumulate filter statistics
                filter_stats = stats["filter_stats"]
                accumulated_filter_stats.duplicates += filter_stats.duplicates
                accumulated_filter_stats.short += filter_stats.short
                accumulated_filter_stats.ambiguous += filter_stats.ambiguous
                accumulated_filter_stats.noncanonical += filter_stats.noncanonical
                accumulated_filter_stats.overlap += filter_stats.overlap
                accumulated_filter_stats.isoform += filter_stats.isoform
                accumulated_filter_stats.total_introns += filter_stats.total_introns
                accumulated_filter_stats.kept_introns += filter_stats.kept_introns

                # Accumulate duplicate and overlap maps
                accumulated_duplicate_map.update(stats["duplicate_map"])
                accumulated_overlap_map.update(stats["overlap_map"])

                # Accumulate the scored intron set for the post-classification pipeline. The
                # raw-feature post-process keys off this set's COUNT (it reads scores from
                # score_info.iic on disk); the parallel arrays carry raw motif scores so the
                # in-memory and streaming paths feed the pipeline an identical intron set.
                for intron in classified_introns:
                    if (
                        intron.scores
                        and intron.scores.svm_score is not None
                        and intron.metadata
                        and intron.metadata.type_id in ('u12', 'u2')
                    ):
                        accumulated_five_z.append(intron.scores.five_raw_score)
                        accumulated_bp_z.append(intron.scores.bp_raw_score)
                        accumulated_three_z.append(intron.scores.three_raw_score)
                        accumulated_svm_scores.append(intron.scores.svm_score)
                        accumulated_type_ids.append(intron.metadata.type_id)

                # Update progress bar based on this contig's length
                contig_idx = contig_to_index[contig_name]
                contig_length = contig_length_list[contig_idx]
                completed_length += contig_length
                progress.update(task, completed=completed_length)
        except Exception as e:
            messenger.error(f"Streaming classification failed: {e}")
            raise

    # Merge classified introns with omitted introns for complete output
    all_introns_for_output = merge_scored_and_omitted_introns(
        all_classified_introns, all_filtered_introns, messenger
    )

    # Sort by coordinates for consistent, deterministic output order
    all_introns_for_output.sort(
        key=lambda i: (
            i.coordinates.chromosome,
            i.coordinates.start,
            i.coordinates.stop,
        )
    )

    with output_writer:
        for intron in all_introns_for_output:
            output_writer.write_intron(intron)

    del all_classified_introns, all_filtered_introns, all_introns_for_output
    gc.collect()

    # Free global resources and cleanup temp files
    annotation_store.cleanup()  # Delete SQLite annotation database
    gc.collect()

    # =========================================================================
    # POST-HOC CLUSTER VALIDATION AND SCORE ADJUSTMENT
    # 1. Compute valley depth from 2D (5'z, BPz) KDE bimodality
    # 2. Apply log-odds score adjustment: valley prior + σ penalty
    # 3. Rewrite output files with adjusted_score and updated rel_score
    #
    # Mode-separation (v2.6+) replaces this branch with a two-pass
    # recalibration when the loaded bundle is in modesep mode.
    # =========================================================================
    _post = _run_post_classification_pipeline(
        five_z=np.array(accumulated_five_z),
        bp_z=np.array(accumulated_bp_z),
        three_z=np.array(accumulated_three_z),
        svm_s=np.array(accumulated_svm_scores),
        type_ids=np.array(accumulated_type_ids),
        model_data=model_data,
        config=config,
        messenger=messenger,
    )

    # Free accumulated score data
    del accumulated_five_z, accumulated_bp_z, accumulated_three_z, accumulated_svm_scores, accumulated_type_ids

    # Build summary via the shared finalizer (same dict in-memory mode produces).
    summary = _finalize_classification_metrics(
        output_writer.get_summary(),
        total_genes=total_genes,
        total_introns_generated=total_introns_generated,
        total_scored=total_scored,
        model_path=str(config.training.pretrained_model_path),
        streaming_mode="per_contig",
        normalizer_used=scaler_source,
        meta_path=config.output.get_output_path(".meta.iic"),
        motif_category=_post.motif_category,
        z_excess=_post.z_excess,
        motif_called_u12=_post.motif_called_u12,
    )

    messenger.info(
        f"Streaming classification complete: {total_classified:,} introns classified"
    )
    # Surface the normalizer used so users see at a glance whether adaptive
    # fitting actually ran or whether we fell back to frozen (e.g. sparse
    # input).
    _NORMALIZER_LABEL = {
        "adaptive_fit": "adaptive (fitted per-species)",
        "user_supplied": "user-supplied (--load-normalizer)",
        "bundled_saved": "bundled saved scaler (training-species)",
        "frozen_fallback_min_introns": "frozen fallback (input below adaptive-fit threshold)",
        "frozen_fallback_no_introns": "frozen fallback (no scoreable introns produced)",
    }
    _norm_label = _NORMALIZER_LABEL.get(scaler_source, str(scaler_source))
    if scaler_source and scaler_source.startswith("frozen_fallback"):
        messenger.warning(f"Normalizer used: {_norm_label}")
    else:
        messenger.info(f"Normalizer used: {_norm_label}")
    messenger.log_only(
        f"Total genes: {total_genes:,}, introns generated: {total_introns_generated:,}"
    )

    # Display filtering summary
    messenger.print_filtering_summary(
        total=accumulated_filter_stats.total_introns,
        short=accumulated_filter_stats.short,
        ambiguous=accumulated_filter_stats.ambiguous,
        noncanonical=accumulated_filter_stats.noncanonical,
        isoform=accumulated_filter_stats.isoform,
        overlap=accumulated_filter_stats.overlap,
        duplicates=accumulated_filter_stats.duplicates,
        kept=accumulated_filter_stats.kept_introns,
        include_duplicates=config.extraction.include_duplicates,
        include_isoforms=False,  # Streaming always uses longest_only=True
        exclude_noncanonical=config.scoring.exclude_noncanonical,
        exclude_overlap=config.extraction.no_intron_overlap,
    )

    # Display classification summary from the FINAL summary (post-discount), so the
    # console counts + dinucleotide breakdowns match metrics.iic.json and the
    # in-memory path (not the writer's stale first-pass high_confidence_u12).
    _u12_boundaries = summary.get("u12_boundaries", {})
    _u2_boundaries = summary.get("u2_boundaries", {})
    _u12_total = summary.get("u12_count", 0)
    messenger.print_classification_results(
        total=total_classified,
        u12_count=_u12_total,
        u2_count=total_classified - _u12_total,
        atac_count=_u12_boundaries.get("AT-AC", 0),
        threshold=config.scoring.threshold,
    )

    # Display boundary statistics
    if _u12_boundaries:
        sorted_boundaries = sorted(_u12_boundaries.items(), key=lambda x: (-x[1], x[0]))
        messenger.print_dinucleotide_boundaries(
            intron_type="U12-type", boundaries=sorted_boundaries, top_n=20
        )

    if _u2_boundaries:
        sorted_boundaries = sorted(_u2_boundaries.items(), key=lambda x: (-x[1], x[0]))
        messenger.print_dinucleotide_boundaries(
            intron_type="U2-type", boundaries=sorted_boundaries, top_n=20
        )

    # Write duplicate and overlap mapping files if there are any mappings
    if accumulated_duplicate_map:
        dupe_map_path = (
            config.output.output_dir / f"{config.output.base_filename}.dupe_map.iic"
        )
        with MappingWriter(dupe_map_path) as dupe_writer:
            dupe_writer.write_mappings(accumulated_duplicate_map)
        messenger.log_only(
            f"Wrote {dupe_writer.mappings_written:,} duplicate mappings to {dupe_map_path}"
        )

    if accumulated_overlap_map:
        overlap_map_path = (
            config.output.output_dir / f"{config.output.base_filename}.overlap_map.iic"
        )
        with MappingWriter(overlap_map_path) as overlap_writer:
            overlap_writer.write_mappings(accumulated_overlap_map)
        messenger.log_only(
            f"Wrote {overlap_writer.mappings_written:,} overlap mappings to {overlap_map_path}"
        )

    # Clean up temporary annotation database
    temp_dir = annotation_db_path.parent
    if temp_dir.name.startswith("intronIC_annotations_"):
        shutil.rmtree(temp_dir)
        messenger.log_only(f"Cleaned up temporary directory: {temp_dir}")

    return total_classified, summary


def classify_with_pretrained_model(
    introns: List[Intron],
    model_path: Path,
    config: IntronICConfig,
    messenger: "UnifiedMessenger",
    reporter: IntronICProgressReporter,
) -> Tuple[List[Intron], dict]:
    """Classify introns using a pretrained model with cross-species domain adaptation.

    This function enables applying a trained SVM to new species without curated references.
    It implements unsupervised domain adaptation by fitting the normalizer on unlabeled
    experimental data, which is statistically valid because:
    - 99.5% of introns are U2-type → robust estimators learn species-specific U2 baseline
    - No label leakage (normalization uses only marginal feature distribution)
    - Corrects covariate shift from species-specific sequence characteristics

    Args:
        introns: Experimental introns (scored, not normalized)
        model_path: Path to pretrained model file (.model.pkl)
        config: Pipeline configuration
        messenger: Unified messenger for console and log output
        reporter: Progress reporter

    Returns:
        Tuple of (classified introns, classification metrics)

    Note:
        The saved normalizer from training is NOT used. Instead, we fit a new normalizer
        on the experimental data to correct for species-specific score distributions.
    """
    import joblib

    # Handle empty introns list
    if not introns:
        messenger.warning("No introns to classify")
        return [], {
            "pretrained": True,
            "model_path": str(model_path),
            "n_classified": 0,
        }, None

    # Drop coordinate-duplicate introns up-front so the adaptive normalizer
    # is fit on the same intron set the streaming-classify path uses (which
    # filters duplicates in its per-contig pre-pass) and so the same set is
    # transformed and classified. Without this, in-memory fits on a larger
    # set with duplicate observations skewing the median/IQR, then writes
    # different z-scores than streaming for the same locus.
    if not config.extraction.include_duplicates:
        n_before = len(introns)
        introns = [i for i in introns if not (i.metadata and i.metadata.duplicate)]
        n_dropped = n_before - len(introns)
        if n_dropped > 0:
            messenger.log_only(
                f"Dropping {n_dropped:,} coordinate-duplicate introns before "
                "classification (use -d to include)"
            )

    messenger.info(f"Loading pretrained model from {model_path}")

    # Load model bundle
    if not model_path.exists():
        raise FileNotFoundError(f"Pretrained model not found: {model_path}")

    model_data = normalize_model_bundle(load_model(model_path))
    assert_scoreable_bundle(model_data)  # fail fast on legacy z-stack bundles (supplant step 2)

    # Handle both old format (SVMEnsemble directly) and new format (dict bundle)
    if isinstance(model_data, dict):
        # New format: {'ensemble': ..., 'normalizer': ..., 'threshold': ...}
        # v3 multispecies bundles are normalized into this shape upstream.
        ensemble = model_data["ensemble"]
        saved_normalizer = model_data.get("normalizer", None)
        _saved_threshold = model_data.get("threshold", config.scoring.threshold)  # noqa: F841
        if "_v3_bundle" in model_data:
            v3_id = model_data["_v3_bundle"].get("model_id", "<unknown>")
            messenger.log_only(f"Loaded v3 multispecies bundle (model_id={v3_id})")
        else:
            messenger.log_only("Loaded model bundle (dict format)")
    else:
        # Old format: SVMEnsemble directly (backward compatibility)
        ensemble = model_data
        saved_normalizer = None
        _saved_threshold = config.scoring.threshold  # Reserved for future use
        messenger.log_only(
            "Loaded model ensemble (legacy format - backward compatibility)"
        )

    messenger.log_only(f"Loaded ensemble with {len(ensemble.models)} models")
    messenger.log_only(f"Using threshold: {config.scoring.threshold}")

    # PROBE (2b-2): skip z-normalization — raw-feature bundles predict directly from raw scores.
    normalizer_source = "none (raw features)"
    normalized_introns = introns

    # Classify using loaded ensemble
    messenger.info("Classifying with pretrained model")

    from intronIC.classification.predictor import SVMPredictor

    predictor = SVMPredictor(
        threshold=config.scoring.threshold,  # Use config threshold (can override saved)
        n_jobs=config.performance.processes,
    )

    classified_introns = list(predictor.predict(ensemble, normalized_introns))
    messenger.log_only(f"Classified {len(classified_introns)} introns")

    # Apply prior adjustment if species-specific prior was provided
    if config.scoring.species_prior is not None:
        training_prior = model_data.get("training_prior")
        if training_prior is not None:
            messenger.info(
                f"Applying prior adjustment for target species (π={config.scoring.species_prior:.2e})"
            )
            classified_introns = _apply_prior_adjustment(
                classified_introns=classified_introns,
                training_prior=training_prior,
                target_prior=config.scoring.species_prior,
                threshold=config.scoring.threshold,
                messenger=messenger,
            )
        else:
            messenger.warning(
                "Cannot apply prior adjustment: model lacks training prior information"
            )
            messenger.warning(
                "Retrain model with current version to enable this feature"
            )

    # Create metrics (limited since we skipped training)
    metrics = {
        "optimized_C": "N/A (pretrained)",
        "cv_score": "N/A (pretrained)",
        "n_models": len(ensemble.models),
        "pretrained": True,
        "model_path": str(model_path),
        "normalizer_used": normalizer_source,
    }
    # Surface the normalizer choice so users see it next to other run-level
    # info; warn when we fell back rather than fitted adaptively.
    _NORMALIZER_LABEL = {
        "adaptive_fit": "adaptive (fitted per-species)",
        "user_supplied": "user-supplied (--load-normalizer)",
        "bundled_human": "bundled human/training-species scaler",
        "frozen_fallback_min_introns": "frozen fallback (input below adaptive-fit threshold)",
    }
    _norm_label = _NORMALIZER_LABEL.get(normalizer_source, str(normalizer_source))
    if normalizer_source and normalizer_source.startswith("frozen_fallback"):
        messenger.warning(f"Normalizer used: {_norm_label}")
    else:
        messenger.info(f"Normalizer used: {_norm_label}")

    # Generate metadata for pretrained model usage
    run_metadata_path = config.output.get_output_path(".run_metadata.json")
    messenger.log_only("Recording pretrained model usage")
    run_metadata = generate_pretrained_metadata(
        model_path=model_path, threshold=config.scoring.threshold
    )
    write_metadata(run_metadata, run_metadata_path)
    messenger.log_only(f"Run metadata saved to {run_metadata_path}")

    # Return model_data so the in-memory caller can dispatch mode-sep + v2.7
    # discount + unified labels in its post-classification block (same pipeline
    # the streaming path runs). Without this the in-memory path silently falls
    # back to the v2.3-era cluster_validation only — see task #203.
    return classified_introns, metrics, model_data



def write_outputs(
    introns: List[Intron],
    config: IntronICConfig,
    messenger: "UnifiedMessenger",
    reporter: IntronICProgressReporter,
    scored_only: Optional[List[Intron]] = None,
    skip_sequences: bool = False,
):
    """Write output files.

    Port from: intronIC.py:4820-4912 (filter_introns_write_files) and 5232-5267 (main)

    Args:
        introns: All introns for .bed.iic, .meta.iic, .introns.iic (scored + omitted)
        config: Pipeline configuration
        messenger: Unified messenger for console and log output
        reporter: Progress reporter
        scored_only: Introns for .score_info.iic (scored only, no omitted). If None, uses introns.
        skip_sequences: If True, skip writing .introns.iic (already written earlier)

    Notes:
        Original intronIC writes different intron sets to different files:
        - .introns.iic/.seqs.iic: ALL introns except duplicates (unless -d)
        - .bed.iic, .meta.iic: Scored + omitted non-duplicates
        - .score_info.iic: ONLY scored introns (no omitted)

        This is achieved by:
        1. Two-phase writing for .bed/.meta (omitted in filter function, scored in main)
        2. Single write for .seqs (in filter function)
        3. Single write for .score_info (in main, only finalized_introns)
    """
    messenger.info("Writing output files")

    # Filter duplicates from BOTH the all-introns list AND the score-info
    # source list. Without filtering ``scored_only``, coordinate-duplicate
    # introns get an SVM score and end up in score_info.iic even though
    # they're absent from .bed/.meta/.introns. That breaks streaming-vs-
    # in-memory equivalence; --include_duplicates is the explicit opt-in.
    # Port from: intronIC.py:4806-4807
    if not config.extraction.include_duplicates:
        original_count = len(introns)
        introns = [i for i in introns if not (i.metadata and i.metadata.duplicate)]
        filtered_count = original_count - len(introns)
        if filtered_count > 0:
            messenger.log_only(
                f"Filtered out {filtered_count} duplicate introns (use -d to include)"
            )
        if scored_only is not None:
            scored_only = [
                i for i in scored_only
                if not (i.metadata and i.metadata.duplicate)
            ]

    output_dir = config.output.output_dir
    base_name = config.output.base_filename
    species_name = config.output.species_name
    simple_name = config.output.uninformative_naming
    no_abbreviate = config.output.no_abbreviate

    # Write BED file (scored + omitted non-duplicates)
    # Skip BED output for sequence-only introns (no real genomic coordinates)
    # Port from: intronIC.py:4823-4835 (omitted) + 5232-5237 (scored)
    from intronIC.core.intron import IntronFlags

    has_real_coordinates = any(
        intron.metadata is not None
        and IntronFlags.SEQUENCE_ONLY not in intron.metadata.flags
        for intron in introns
    )

    if has_real_coordinates:
        bed_path = output_dir / f"{base_name}.bed.iic"
        messenger.log_only(f"Writing BED file: {bed_path}")
        bed_writer = BEDWriter(bed_path)
        with bed_writer:
            for intron in introns:
                bed_writer.write_intron(
                    intron,
                    species_name=species_name,
                    simple_name=simple_name,
                    no_abbreviate=no_abbreviate,
                )
        messenger.log_only(f"Wrote {len(introns)} introns to BED file")
    else:
        messenger.log_only(
            "Skipping BED file (sequence-only input has no genomic coordinates)"
        )
        bed_path = None

    # Write metadata file (scored + omitted non-duplicates)
    # Port from: intronIC.py:4823-4835 (omitted) + 5240, 5262-5264 (scored)
    meta_path = output_dir / f"{base_name}.meta.iic"
    messenger.log_only(f"Writing metadata file: {meta_path}")
    meta_writer = MetaWriter(meta_path)
    with meta_writer:
        if not config.output.no_headers:
            meta_writer.write_header()
        for intron in introns:
            meta_writer.write_intron(
                intron,
                species_name=species_name,
                simple_name=simple_name,
                no_abbreviate=no_abbreviate,
            )
    messenger.log_only(f"Wrote metadata for {len(introns)} introns")

    # Write sequences file (all introns except duplicates, unless -d)
    # Port from: intronIC.py:4842-4845
    seq_path = output_dir / f"{base_name}.introns.iic"
    if not skip_sequences:
        messenger.log_only(f"Writing sequences file: {seq_path}")
        seq_writer = SequenceWriter(seq_path)
        with seq_writer:
            for intron in introns:
                seq_writer.write_intron(
                    intron,
                    species_name=species_name,
                    simple_name=simple_name,
                    no_abbreviate=no_abbreviate,
                )
        messenger.log_only(f"Wrote sequences for {len(introns)} introns")
    else:
        messenger.log_only(
            f"Sequences already written during scoring phase: {seq_path}"
        )

    # Write score info file (ONLY scored introns, no omitted)
    # Port from: intronIC.py:5239-5261
    # CRITICAL: Only write scored introns, not omitted ones
    score_introns = scored_only if scored_only is not None else introns
    score_path = output_dir / f"{base_name}.score_info.iic"
    messenger.log_only(f"Writing score info file: {score_path}")
    score_writer = ScoreWriter(score_path)
    with score_writer:
        if not config.output.no_headers:
            score_writer.write_header()
        for intron in score_introns:
            score_writer.write_intron(
                intron,
                species_name=species_name,
                simple_name=simple_name,
                no_abbreviate=no_abbreviate,
            )
    messenger.log_only(f"Wrote score info for {len(score_introns)} introns")

    output_files: dict[str, str] = {
        "Metadata": str(meta_path),
        "Sequences": str(seq_path),
        "Scores": str(score_path),
        "Log": str(config.output.get_output_path(".iic.log")),
    }

    # Only include BED file if we wrote one (requires real genomic coordinates)
    if bed_path:
        output_files["BED"] = str(bed_path)

    messenger.print_file_tree(output_files)
    messenger.success("All output files written successfully")


def main_train(config: IntronICConfig):
    """Train a model on reference data only (no genome/annotation needed).

    This is a pure training workflow:
    1. Load reference U12/U2 sequences
    2. Score reference sequences with PWM matrices
    3. Normalize scores
    4. Train ensemble of models with cross-validation
    5. Save model and metadata to disk

    Args:
        config: Training configuration

    No genome/annotation required!
    """
    # Track start time
    start_time = time.time()

    # Setup logging
    logger, log_console = setup_logging(config)
    reporter = IntronICProgressReporter(quiet=config.output.quiet)

    from .messenger import UnifiedMessenger

    messenger = UnifiedMessenger(
        console=reporter.console,
        log_console=log_console,
        logger=logger,
        quiet=config.output.quiet,
    )

    # Print header
    reporter.console.print(
        "\n[bold cyan]═══════════════════════════════════════════════════════[/bold cyan]"
    )
    reporter.console.print(
        "[bold cyan]                  intronIC TRAINING MODE                  [/bold cyan]"
    )
    reporter.console.print(
        "[bold cyan]═══════════════════════════════════════════════════════[/bold cyan]\n"
    )
    messenger.info(f"Training model: {config.output.species_name}")
    messenger.info("No genome/annotation needed - using reference sequences only")

    # Log command and config for reproducibility
    import sys

    messenger.log_only("=" * 80)
    messenger.log_only("COMMAND AND CONFIGURATION")
    messenger.log_only("=" * 80)

    # Format command for easy copy/paste (one line, let terminal wrap)
    messenger.log_only(f"Command: {' '.join(sys.argv)}")
    messenger.log_only(f"Working directory: {Path.cwd()}")
    messenger.log_only("")

    # Log reference input files with full paths
    messenger.log_only("Input Files:")
    if config.scoring.reference_u12s:
        messenger.log_only(
            f"  U12 reference: {config.scoring.reference_u12s.absolute()}"
        )
    if config.scoring.reference_u2s:
        messenger.log_only(f"  U2 reference: {config.scoring.reference_u2s.absolute()}")
    messenger.log_only("")

    # Log config file path if loaded
    if config.config_path:
        messenger.log_only(f"Config file: {config.config_path.absolute()}")
        messenger.log_only("")

    # Log key configuration parameters in condensed format
    messenger.log_only("Configuration Parameters:")
    messenger.log_only(f"  Run name: {config.output.species_name or 'N/A'}")
    messenger.log_only(f"  Random seed: {config.training.seed}")
    messenger.log_only(f"  Classification threshold: {config.scoring.threshold}%")
    messenger.log_only("")

    messenger.log_only("Training:")
    messenger.log_only(f"  n_models: {config.ensemble.n_models}")
    messenger.log_only(f"  max_iter: {config.ensemble.max_iter}")
    messenger.log_only(f"  eval_mode: {config.training.eval_mode}")
    if config.training.fixed_C:
        messenger.log_only(f"  C: {config.training.fixed_C:.6e} (fixed)")
    else:
        messenger.log_only("  C: optimized via grid search")
        messenger.log_only(f"  n_optimization_rounds: {config.optimizer.n_rounds}")
        messenger.log_only(f"  n_cv_folds: {config.optimizer.cv_folds}")
    messenger.log_only("")

    # Log optimizer-specific config
    opt = config.optimizer
    messenger.log_only("Optimizer:")
    messenger.log_only(f"  penalty_options: {list(opt.penalty_options)}")
    if opt.features:
        features_str = ", ".join(opt.features)
        messenger.log_only(f"  features: [{features_str}]")
    else:
        messenger.log_only("  features: default (4D)")
    if opt.gamma_imbalance_options:
        messenger.log_only(
            f"  gamma_imbalance_options: {list(opt.gamma_imbalance_options)}"
        )
    messenger.log_only(
        f"  class_weight_multipliers: {list(opt.class_weight_multipliers)}"
    )
    messenger.log_only("")

    messenger.log_only("=" * 80)
    messenger.log_only("")

    pipeline_steps = [
        "Load reference data",
        "Score reference sequences",
        "Normalize scores",
        "Train classifier",
    ]
    messenger.console_only("")
    reporter.print_pipeline_steps(pipeline_steps)

    try:
        # Steps 1-3: Load, score, and normalize reference data
        messenger.step(1, "Load Reference Data", pipeline_steps)
        messenger.step(2, "Score Reference Sequences", pipeline_steps)
        messenger.step(3, "Normalize Scores", pipeline_steps)

        # Load and process reference sequences
        # This duplicates some logic from normalize_scores() but avoids issues with empty experimental introns
        data_dir = Path(__file__).parent.parent / "data"

        u12_file = config.scoring.reference_u12s or (
            data_dir / "u12_reference_multispecies.introns.iic.gz"
        )
        u2_file = config.scoring.reference_u2s or (
            data_dir / "u2_reference_multispecies.introns.iic.gz"
        )

        if not u12_file.exists() or not u2_file.exists():
            raise FileNotFoundError(
                f"Reference data not found. U12: {u12_file}, U2: {u2_file}"
            )

        messenger.log_only("Loading reference sequences")
        u12_reference = load_reference_sequences(u12_file, messenger=messenger)
        u2_reference = load_reference_sequences(u2_file, messenger=messenger)
        messenger.log_only(
            f"Loaded {len(u12_reference)} U12 and {len(u2_reference)} U2 reference introns"
        )

        # Score reference introns
        messenger.log_only("Scoring reference sequences")
        all_reference = u12_reference + u2_reference
        scored_reference = score_introns(all_reference, config, messenger, reporter)

        # Split back into U12 and U2
        u12_scored = scored_reference[: len(u12_reference)]
        u2_scored = scored_reference[len(u12_reference) :]

        # Raw-feature scoring modes (raw_gated / pmotif_adjudicated): skip z-normalization and the
        # z-only IntronClassifier path; train the raw ensemble directly and stamp the bundle so inference
        # runs the species-gate / pmotif adjudicator. See docs/raw_gated_scoring.md §0c.
        if config.training.scoring_mode in ("raw_gated", "pmotif_adjudicated"):
            messenger.success(
                f"Loaded and scored {len(u12_scored)} U12 and {len(u2_scored)} U2 reference introns "
                f"(raw features; z-normalization skipped for scoring_mode={config.training.scoring_mode})"
            )
            messenger.step(4, "Train Classifier", pipeline_steps)
            _train_and_save_raw_bundle(u12_scored, u2_scored, config, messenger)
            messenger.success("Model trained and saved")
        else:
            raise ValueError(
                f"scoring_mode={config.training.scoring_mode!r} is not trainable: the z-normalization "
                f"training path (ScoreNormalizer + IntronClassifier) was removed (supplant 2b). Train "
                f"with --scoring-mode pmotif_adjudicated (default) or raw_gated; check out the "
                f"pre-zstack-removal git tag to train a legacy zscore bundle."
            )

        # Print final summary
        elapsed = time.time() - start_time
        hours, remainder = divmod(int(elapsed), 3600)
        minutes, seconds = divmod(remainder, 60)

        messenger.console_only("")
        if hours > 0:
            time_str = f"{hours}h {minutes}m {seconds}s"
        elif minutes > 0:
            time_str = f"{minutes}m {seconds}s"
        else:
            time_str = f"{seconds}s"

        messenger.success(f"Training complete! (Runtime: {time_str})")

    except Exception as e:
        messenger.error(f"Training failed: {str(e)}")
        raise


def main_extract(config: IntronICConfig):
    """Extract intron sequences without classification.

    This workflow extracts introns and writes sequence files but does not
    perform scoring or classification. Useful for:
    1. Creating custom reference sets for training
    2. Extracting introns for external analysis
    3. Preparing data for manual curation

    Workflow:
    1. Load genome and annotation/BED
    2. Extract introns
    3. Filter introns
    4. Write sequence outputs (no scores or classification)

    Args:
        config: Pipeline configuration
    """
    # Track start time
    start_time = time.time()

    # Setup logging
    logger, log_console = setup_logging(config)
    reporter = IntronICProgressReporter(quiet=config.output.quiet)

    from .messenger import UnifiedMessenger

    messenger = UnifiedMessenger(
        console=reporter.console,
        log_console=log_console,
        logger=logger,
        quiet=config.output.quiet,
    )

    # Print startup banner
    import sys

    messenger.print_startup_banner(
        species_name=config.output.species_name,
        input_mode=config.input.mode,
        output_dir=str(config.output.output_dir.absolute()),
        threshold=None,  # No threshold for extract mode
        command=" ".join(sys.argv),
        working_dir=str(Path.cwd()),
        model_path=None,  # No model for extract mode
        genome_path=str(config.input.genome.absolute())
        if config.input.mode in ["annotation", "bed"] and config.input.genome
        else None,
        annotation_path=str(config.input.annotation.absolute())
        if config.input.mode == "annotation" and config.input.annotation
        else None,
        bed_path=str(config.input.bed.absolute())
        if config.input.mode == "bed" and config.input.bed
        else None,
        sequences_path=None,
    )

    # Note: Extract mode only supports annotation and bed modes (not pre-extracted sequences)
    if config.input.mode not in ["annotation", "bed"]:
        messenger.error(
            f"Extract mode only supports annotation or BED input, not '{config.input.mode}'"
        )
        raise ValueError(
            f"Extract mode requires genome + annotation or genome + BED, got mode: {config.input.mode}"
        )

    # Extract introns (this includes filtering)
    messenger.info("=" * 80)
    messenger.info("EXTRACTION MODE - Extracting intron sequences")
    messenger.info("=" * 80)

    # Load genome
    messenger.info(f"Loading genome: {config.input.genome}")
    genome_reader = GenomeReader(str(config.input.genome))

    # Count sequences (simple count, no caching)
    seq_count = sum(1 for _ in genome_reader.stream())
    messenger.success(f"Genome loaded: {seq_count:,} sequences")

    # Extract based on input mode
    if config.input.mode == "annotation":
        if config.performance.streaming:
            # Streaming extraction
            introns, db_path = extract_introns_streaming(config, messenger, reporter)
            messenger.success(f"Extracted {len(introns):,} introns (streaming mode)")
        else:
            # Standard extraction (gene count unused in extract-only mode)
            introns, _ = extract_introns_from_annotation(config, messenger, reporter)
            messenger.success(f"Extracted {len(introns):,} introns")
    elif config.input.mode == "bed":
        # BED extraction
        introns = extract_introns_from_bed(config, genome_reader, messenger, reporter)
        messenger.success(f"Extracted {len(introns):,} introns from BED")

    # Write outputs (sequences only, no scores/classification)
    messenger.info("=" * 80)
    messenger.info("Writing output files")
    messenger.info("=" * 80)

    # Write intron sequences
    seq_path = config.output.get_output_path(".introns.iic")
    meta_path = config.output.get_output_path(".meta.iic")
    bed_path = config.output.get_output_path(".bed.iic")

    from intronIC.file_io.writers import BEDWriter, MetaWriter, SequenceWriter

    # If streaming mode, read sequences from SQLite
    if config.performance.streaming and config.input.mode == "annotation":
        messenger.info(f"Writing sequences from SQLite: {seq_path}")
        from intronIC.file_io.sequence_store import StreamingSequenceStore

        store = StreamingSequenceStore(db_path)
        seq_writer = SequenceWriter(seq_path)
        with seq_writer:
            for row in store.iter_all():
                # Write using raw sequence data from SQLite
                seq_writer.write_from_row(
                    intron_id=row.intron_id,
                    formatted_name=row.formatted_name,
                    upstream_flank=row.upstream_flank,
                    seq=row.seq,
                    downstream_flank=row.downstream_flank,
                    terminal_dnts=row.terminal_dnts,
                    svm_score=None,  # No scores in extract mode
                )
        store.close()

        # Clean up temporary directory
        temp_dir = db_path.parent
        if temp_dir.name.startswith("intronIC_sequences_"):
            shutil.rmtree(temp_dir)
            messenger.log_only(f"Cleaned up temporary directory: {temp_dir}")
    else:
        # Standard mode - write from memory
        # Filter introns with sequences
        introns_with_seqs = [i for i in introns if i.has_sequences]
        messenger.info(
            f"Writing sequences: {seq_path} ({len(introns_with_seqs):,} introns with sequences)"
        )
        seq_writer = SequenceWriter(seq_path)
        with seq_writer:
            seq_writer.write_introns(
                introns_with_seqs,
                species_name=config.output.species_name,
                simple_name=not config.output.no_abbreviate,
                no_abbreviate=config.output.no_abbreviate,
                include_score=False,  # No scores in extract mode
            )

    # Write metadata
    messenger.info(f"Writing metadata: {meta_path}")
    meta_writer = MetaWriter(meta_path)
    with meta_writer:
        meta_writer.write_introns(
            introns,
            species_name=config.output.species_name,
            simple_name=not config.output.no_abbreviate,
        )

    # Write BED
    messenger.info(f"Writing BED: {bed_path}")
    bed_writer = BEDWriter(bed_path)
    with bed_writer:
        bed_writer.write_introns(
            introns,
            species_name=config.output.species_name,
            simple_name=not config.output.no_abbreviate,
            no_abbreviate=config.output.no_abbreviate,
        )

    # Calculate runtime
    elapsed_seconds = time.time() - start_time
    hours, remainder = divmod(int(elapsed_seconds), 3600)
    minutes, seconds = divmod(remainder, 60)

    if hours > 0:
        runtime_str = f"{hours}h {minutes}m {seconds}s"
    elif minutes > 0:
        runtime_str = f"{minutes}m {seconds}s"
    else:
        runtime_str = f"{seconds}s"

    messenger.success(f"Extraction complete! (Runtime: {runtime_str})")
    messenger.info(f"Output files written to: {config.output.output_dir}")


def main_classify(config: IntronICConfig):
    """Run the complete intronIC classification pipeline.

    This is the standard workflow:
    1. Load input data (genome/annotation/bed/sequences)
    2. Extract introns
    3. Filter introns
    4. Score introns
    5. Load or train model
    6. Classify introns
    7. Write outputs

    Args:
        config: Pipeline configuration
    """
    # Track start time for runtime reporting
    start_time = time.time()

    # Setup logging and reporting with ANSI color support in log files
    logger, log_console = setup_logging(config)
    reporter = IntronICProgressReporter(quiet=config.output.quiet)

    # Create unified messenger for synchronized console + log output
    # Both destinations get Rich formatting with ANSI colors
    from .messenger import UnifiedMessenger

    messenger = UnifiedMessenger(
        console=reporter.console,
        log_console=log_console,
        logger=logger,
        quiet=config.output.quiet,
    )

    # Print unified startup banner (both console and log)
    # Convert paths to absolute for logging
    import sys

    messenger.print_startup_banner(
        species_name=config.output.species_name,
        input_mode=config.input.mode,
        output_dir=str(config.output.output_dir.absolute()),
        threshold=config.scoring.threshold,
        command=" ".join(sys.argv),
        working_dir=str(Path.cwd()),
        model_path=str(config.training.pretrained_model_path.absolute())
        if config.training.pretrained_model_path
        else None,
        genome_path=str(config.input.genome.absolute())
        if config.input.mode in ["annotation", "bed"] and config.input.genome
        else None,
        annotation_path=str(config.input.annotation.absolute())
        if config.input.mode == "annotation" and config.input.annotation
        else None,
        bed_path=str(config.input.bed.absolute())
        if config.input.mode == "bed" and config.input.bed
        else None,
        sequences_path=str(config.input.sequence_file.absolute())
        if config.input.mode == "sequences" and config.input.sequence_file
        else None,
    )

    # Define pipeline steps
    pipeline_steps = [
        "Load and extract introns",
        "Score introns with PWMs",
        "Normalize scores",
        "Train and apply classifier",
        "Write output files",
    ]

    try:
        # TRUE STREAMING MODE: Per-contig processing with immediate output
        # This mode provides ~90% memory savings but requires:
        # 1. Pre-trained model (has frozen scaler)
        # 2. Annotation input mode
        # Species background correction is handled via a lightweight first-pass
        # in classify_streaming_per_contig() (genome fetch, not full extraction).
        if (
            config.performance.streaming
            and config.training.pretrained_model_path
            and config.input.mode == "annotation"
        ):
            total_classified, summary = classify_streaming_per_contig(
                config, messenger, reporter
            )

            # Save classification metrics
            metrics_path = config.output.get_output_path(".metrics.iic.json")
            messenger.log_only(f"Saving classification metrics to {metrics_path}")
            with open(metrics_path, "w") as f:
                json.dump(summary, f, indent=2)

            # Generate visualization plots from output file
            messenger.log_only("Generating visualization plots")
            try:
                from intronIC.visualization.plots import (
                    plot_classification_results_from_file,
                )

                score_file = config.output.get_output_path(".score_info.iic")
                if score_file.exists():
                    plot_classification_results_from_file(
                        score_file=score_file,
                        output_dir=config.output.output_dir,
                        species_name=config.output.base_filename,
                        threshold=config.scoring.threshold,
                        fig_dpi=300,
                    )
                    messenger.log_only("Successfully generated classification plots")
                else:
                    messenger.warning(f"Score file not found: {score_file}")
            except Exception as plot_error:
                import traceback

                messenger.warning(f"Failed to generate plots: {plot_error}")
                messenger.warning(f"Traceback: {traceback.format_exc()}")

            # Calculate and log total runtime
            elapsed_seconds = time.time() - start_time
            hours, remainder = divmod(int(elapsed_seconds), 3600)
            minutes, seconds = divmod(remainder, 60)

            if hours > 0:
                runtime_str = f"{hours}h {minutes}m {seconds}s"
            elif minutes > 0:
                runtime_str = f"{minutes}m {seconds}s"
            else:
                runtime_str = f"{seconds}s"

            messenger.success(f"Pipeline complete! (Runtime: {runtime_str})")
            return

        # Step 1: Load and extract introns
        messenger.step(1, "Load and Extract Introns", pipeline_steps)

        # Track streaming mode state for later sequence output
        streaming_db_path = None

        _inmem_total_genes = 0  # set by annotation/standard extraction below
        if config.input.mode == "annotation":
            # Don't load genome here - extraction handles it internally
            # (parallel mode uses indexed access, sequential uses cache)
            if config.performance.streaming:
                # Streaming mode: roughly half the peak memory of in-memory; stores sequences in SQLite
                introns, streaming_db_path = extract_introns_streaming(
                    config, messenger, reporter
                )

                # If species background correction is enabled, accumulate
                # frequencies from the SQLite sequence store before scoring.
                # This avoids loading all sequences into memory — we just
                # iterate the store and count bases at scored positions.
                if config.species_background.enabled and streaming_db_path:
                    from intronIC.file_io.sequence_store import StreamingSequenceStore

                    pwm_sets_for_bg = load_pwms_with_fallback(config, messenger)
                    bg, u12_sets_bg, _, five_len, three_len = (
                        _create_background_accumulator(config, pwm_sets_for_bg)
                    )

                    messenger.info("Species background: reading sequences from SQLite...")
                    five_start = config.scoring.scoring_regions.five_start
                    five_end = config.scoring.scoring_regions.five_end
                    three_start = config.scoring.scoring_regions.three_start
                    three_end = config.scoring.scoring_regions.three_end
                    bp_start = config.scoring.scoring_regions.bp_start
                    bp_end = config.scoring.scoring_regions.bp_end

                    store = StreamingSequenceStore(streaming_db_path)
                    n_acc = 0
                    for row in store.iter_all():
                        seq = row.seq or ""
                        up = row.upstream_flank or ""
                        dn = row.downstream_flank or ""
                        dnts = row.terminal_dnts or ""
                        if "-" not in dnts or len(seq) < 30:
                            continue
                        five_dnt, three_dnt = dnts.split("-")

                        five_seq = (up[five_start:] + seq[:five_end]).upper()
                        three_seq = (seq[three_start:] + dn[:three_end]).upper()
                        bp_seq = seq[
                            max(0, len(seq) + bp_start) : len(seq) + bp_end
                        ].upper()

                        bg.accumulate(
                            row.intron_id,
                            five_dnt.upper(),
                            three_dnt.upper(),
                            five_seq[:five_len],
                            three_seq[-three_len:],
                            bp_seq,
                        )
                        n_acc += 1
                    store.close()

                    bg_cfg = config.species_background
                    _streaming_corrected_pwm_sets = _finalize_background_correction(
                        bg, u12_sets_bg, n_acc, bg_cfg.min_introns, bg_cfg.n0,
                        messenger,
                    )
                else:
                    _streaming_corrected_pwm_sets = None
            else:
                _streaming_corrected_pwm_sets = None
                # Standard processing: faster but uses more memory
                introns, _inmem_total_genes = extract_introns_from_annotation(
                    config, messenger, reporter
                )
        elif config.input.mode == "bed":
            _streaming_corrected_pwm_sets = None
            genome_reader = load_genome(config, messenger)
            introns = extract_introns_from_bed(
                config, genome_reader, messenger, reporter
            )
        elif config.input.mode == "sequences":
            _streaming_corrected_pwm_sets = None
            introns = load_introns_from_sequences(config, messenger)
        else:
            raise ValueError(f"Unknown input mode: {config.input.mode}")

        # Capture generation-time counts for metrics parity with the streaming
        # path (which accumulates the same totals per contig). `introns` here is
        # the full generated set (extract + skip), pre-filter; `_inmem_total_genes`
        # was set by the annotation extraction above (all genes in the hierarchy).
        _inmem_total_generated = len(introns)

        # Filter introns before scoring (duplicates, short introns, longest isoform)
        # This matches original intronIC behavior where filtering happens BEFORE scoring
        # to avoid scoring 5x more introns than necessary (which causes O(n²) slowdown)
        messenger.log_only("Filtering introns for scoring")

        # Create filter with scoring-appropriate settings:
        # - longest_only=True: Only score longest isoform per gene (filters ~8k introns)
        # - include_duplicates=False: Don't score duplicates (filters ~38k introns)
        # - min_length: Filter short introns
        # - allow_noncanonical: Based on exclude_noncanonical flag
        # - allow_overlap: Based on no_intron_overlap flag
        # Note: For sequence file inputs, all introns have unique placeholder coordinates
        # and no grandparent/parent metadata, so longest_only and duplicate filtering
        # effectively become no-ops (all pass through as "longest isoform")
        intron_filter = IntronFilter(
            min_length=config.extraction.min_intron_len,
            bp_matrix_length=7,  # Default from original
            scoring_regions=["five", "three"],  # Check these for ambiguous bases
            allow_noncanonical=not config.scoring.exclude_noncanonical,
            allow_overlap=not config.extraction.no_intron_overlap,
            longest_only=True,  # For sequences: no-op (no grandparent info)
            include_duplicates=False,  # For sequences: no-op (unique coords)
        )

        filtered_introns = intron_filter.filter_introns(introns)

        # Report comprehensive filtering statistics (combines pre-filter + post-extraction)
        stats = intron_filter.stats

        # Report filtering results using unified method (console + log)
        # Pass counts and user options so messenger can show in correct column
        messenger.print_filtering_summary(
            total=stats.total_introns,
            short=stats.short,
            ambiguous=stats.ambiguous,
            noncanonical=stats.noncanonical,
            isoform=stats.isoform,
            overlap=stats.overlap,
            duplicates=stats.duplicates,
            kept=stats.kept_introns,
            # User options that affect which column counts appear in
            include_duplicates=config.extraction.include_duplicates,
            include_isoforms=config.extraction.allow_multiple_isoforms,
            exclude_noncanonical=config.scoring.exclude_noncanonical,
            exclude_overlap=config.extraction.no_intron_overlap,
        )

        messenger.success(f"Processed {len(introns):,} introns from annotation")

        # Important: Use filtered_introns for scoring, but keep original introns list
        # for potential output (user may want duplicates via -d flag).
        #
        # Drop omitted introns (AMBIGUOUS, SHORT, NONCANONICAL, etc.) from the
        # scoring list — IntronFilter has already tagged them and they will be
        # re-merged with NA scores by merge_scored_and_omitted_introns. Without
        # this filter, all-N motifs get scored with raw=0 and produce garbage
        # z-scores / SVM calls in score_info, while the streaming-classify path
        # correctly omits them (see _streaming_extract_and_filter_contig).
        introns_for_scoring = [
            i for i in filtered_introns
            if i.metadata is None or i.metadata.omitted == OmissionReason.NONE
        ]

        # Step 2: Score introns
        messenger.step(2, "Score Introns with PWMs", pipeline_steps)

        scored_introns = score_introns(
            introns_for_scoring, config, messenger, reporter,
            precomputed_pwm_sets=_streaming_corrected_pwm_sets,
        )

        # Defer sequence writing until AFTER classification for both modes
        # so we can include SVM scores. Sequences remain in memory for standard mode
        # until writing is complete (uses more memory but produces correct output).
        # Streaming mode: Sequences are in SQLite, nothing to clear yet.
        if streaming_db_path is None:
            messenger.log_only(
                "Standard mode: sequences will be written after classification"
            )
        else:
            messenger.log_only(
                "Streaming mode: sequences will be written from SQLite after classification"
            )

        # Check if using pretrained model
        if config.training.pretrained_model_path:
            # Skip normalization/training - use pretrained model
            messenger.step(3, "Classify with Pretrained Model", pipeline_steps)
            classified_introns, metrics, model_data = classify_with_pretrained_model(
                scored_introns,
                config.training.pretrained_model_path,
                config,
                messenger,
                reporter,
            )
        else:
            # The no-pretrained-model "normalize + train + classify in one pass" flow used
            # ScoreNormalizer + IntronClassifier (z-only), removed in supplant 2b. Classify now
            # requires a raw-feature bundle (the bundled pmotif default is used when --model is
            # omitted, so this branch is only reached if --model is explicitly cleared).
            raise ValueError(
                "intronIC classify requires a pretrained raw-feature bundle. The no-model "
                "train-and-classify flow used z-normalization + IntronClassifier, which were "
                "removed (supplant 2b). Omit --model to use the bundled pmotif default, or pass "
                "--model with a pmotif_adjudicated / raw_gated bundle."
            )

        # The classification summary + dinucleotide boundary tables are printed
        # later, after post-classification finalizes type_id in meta.iic (see
        # below). Reporting them here would show the first-pass (pre-discount)
        # call and disagree with metrics.iic.json; the streaming path reports
        # after its finalizer for the same reason.

        # Generate visualization plots
        messenger.log_only("Generating visualization plots")
        try:
            plot_classification_results(
                introns=classified_introns,
                output_dir=config.output.output_dir,
                species_name=config.output.base_filename,
                threshold=config.scoring.threshold,
                fig_dpi=300,
            )
            messenger.log_only("Successfully generated classification plots")
        except Exception as plot_error:
            import traceback

            messenger.warning(f"Failed to generate plots: {plot_error}")
            messenger.warning(f"Traceback: {traceback.format_exc()}")
            # Continue even if plotting fails

        # Step 5: Write outputs
        messenger.step(5, "Write Output Files", pipeline_steps)

        # STREAMING MODE: Write sequences from SQLite with SVM scores now that classification is done
        if streaming_db_path is not None:
            messenger.info("Writing intron sequences to file (from SQLite)")
            from intronIC.file_io.sequence_store import StreamingSequenceStore

            seq_output_path = (
                config.output.output_dir / f"{config.output.base_filename}.introns.iic"
            )
            seq_writer = SequenceWriter(seq_output_path)

            # Build lookup of intron_id -> svm_score from classified introns
            score_lookup = {}
            for intron in classified_introns:
                if intron.scores and intron.scores.svm_score is not None:
                    score_lookup[intron.intron_id] = intron.scores.svm_score

            # Build set of duplicate intron_ids to filter (matches standard mode behavior)
            # The metadata is preserved in the lightweight intron objects
            duplicate_ids = set()
            if not config.extraction.include_duplicates:
                for intron in introns:
                    if intron.metadata and intron.metadata.duplicate:
                        duplicate_ids.add(intron.intron_id)

            store = StreamingSequenceStore(streaming_db_path)
            introns_written = 0
            duplicates_skipped = 0

            with seq_writer:
                for row in store.iter_all():
                    # Skip duplicates if not -d flag (matches standard mode)
                    if row.intron_id in duplicate_ids:
                        duplicates_skipped += 1
                        continue

                    # Get score if available (None for omitted introns)
                    svm_score = score_lookup.get(row.intron_id)

                    # Write using raw sequence data from SQLite
                    # Use the pre-computed formatted_name from SQLite
                    seq_writer.write_from_row(
                        intron_id=row.intron_id,
                        formatted_name=row.formatted_name,
                        upstream_flank=row.upstream_flank,
                        seq=row.seq,
                        downstream_flank=row.downstream_flank,
                        terminal_dnts=row.terminal_dnts,
                        svm_score=svm_score,
                    )
                    introns_written += 1

            store.cleanup()  # Delete SQLite database

            # Clean up temporary directory
            temp_dir = streaming_db_path.parent
            if temp_dir.name.startswith("intronIC_sequences_"):
                shutil.rmtree(temp_dir)
                messenger.log_only(f"Cleaned up temporary directory: {temp_dir}")

            if duplicates_skipped > 0:
                messenger.log_only(f"Skipped {duplicates_skipped} duplicate introns")
            messenger.log_only(
                f"Wrote sequences for {introns_written} introns to {seq_output_path.name} (from SQLite)"
            )

        # Merge classified introns with omitted introns for complete meta output
        # This matches original intronIC behavior where .meta.iic includes all introns
        # (scored + omitted), not just the ones that went through classification
        # Note: Use filtered_introns (not introns) because omission reasons are set during filtering
        all_introns_for_output = merge_scored_and_omitted_introns(
            classified_introns, filtered_introns, messenger
        )

        # STANDARD MODE: Write sequences with SVM scores now that classification is done
        if streaming_db_path is None:
            messenger.info("Writing intron sequences to file")
            seq_output_path = (
                config.output.output_dir / f"{config.output.base_filename}.introns.iic"
            )
            seq_writer = SequenceWriter(seq_output_path)

            # Write from all_introns_for_output to include omitted introns
            # Omitted introns have sequences but no scores
            introns_to_write = all_introns_for_output
            introns_written = 0
            with seq_writer:
                for intron in introns_to_write:
                    # Filter duplicates if not -d flag
                    if not config.extraction.include_duplicates:
                        if intron.metadata and intron.metadata.duplicate:
                            continue

                    # Skip introns without sequences (extraction failed or not attempted)
                    if not intron.has_sequences:
                        continue

                    seq_writer.write_intron(
                        intron,
                        species_name=config.output.species_name,
                        simple_name=config.output.uninformative_naming,
                        no_abbreviate=config.output.no_abbreviate,
                        include_score=True,  # Include SVM score (now available after classification)
                    )
                    introns_written += 1
            messenger.log_only(
                f"Wrote sequences for {introns_written} introns to {seq_output_path.name}"
            )

            # Now clear sequences to free memory
            all_introns_for_output = clear_large_sequences_for_classification(
                all_introns_for_output
            )
            introns = clear_large_sequences_for_classification(introns)
            gc.collect()

        # Write outputs with different intron sets for different files:
        # - all_introns_for_output: for .bed.iic, .meta.iic, .introns.iic (scored + omitted)
        # - classified_introns: for .score_info.iic (scored only, no omitted)
        # Port from: intronIC.py writes finalized_introns to .score_info (line 5239-5261)
        write_outputs(
            all_introns_for_output,
            config,
            messenger,
            reporter,
            scored_only=classified_introns,
            skip_sequences=True,  # Sequences already written (standard mode) or just above (streaming mode)
        )

        # Post-classification pipeline: mode-sep (v2.6+) + continuous discount
        # (v2.7) + unified labels (v2.7.1), with legacy cluster-validation as
        # the gate-fail fallback. Mirrors the streaming-classify path so the
        # two routes produce equivalent output (task #203 — May 27 2026,
        # in-memory previously fell through to the v2.3 legacy-only path).
        if classified_introns:
            # Filter to introns with all three z-scores + svm_score present AND
            # a resolved u12/u2 type_id, so the parallel arrays line up. This guard
            # is identical to the streaming accumulator filter
            # (classify_streaming_per_contig) so both paths feed the
            # post-classification pipeline the same intron set — previously the
            # in-memory path admitted metadata-less introns as a fabricated 'u2',
            # which could diverge from streaming on some species.
            cv_introns = [
                i for i in classified_introns
                if i.scores
                and i.scores.svm_score is not None
                and i.metadata
                and i.metadata.type_id in ('u12', 'u2')
            ]
            five_z = np.array([i.scores.five_raw_score for i in cv_introns])
            bp_z = np.array([i.scores.bp_raw_score for i in cv_introns])
            three_z = np.array([i.scores.three_raw_score for i in cv_introns])
            svm_s = np.array([i.scores.svm_score for i in cv_introns])
            type_ids = np.array([i.metadata.type_id for i in cv_introns])

            if len(five_z) > 0:
                # Raw_gated / pmotif_adjudicated post-classification, shared with the
                # streaming-classify path (single source of truth).
                _post = _run_post_classification_pipeline(
                    five_z=five_z, bp_z=bp_z, three_z=three_z,
                    svm_s=svm_s, type_ids=type_ids,
                    model_data=model_data, config=config, messenger=messenger,
                )
                # Rebuild metrics through the shared finalizer so this path's
                # metrics.iic.json matches the streaming path (the only intentional
                # difference is streaming_mode). base_summary supplies only
                # total_introns; all U12/U2/HC counts + boundaries are tallied from the
                # finalized meta.iic inside the finalizer (the single source of truth).
                _base_summary = summarize_introns(
                    all_introns_for_output, config.scoring.threshold
                )
                metrics = _finalize_classification_metrics(
                    _base_summary,
                    total_genes=_inmem_total_genes,
                    total_introns_generated=_inmem_total_generated,
                    total_scored=len(classified_introns),
                    model_path=str(config.training.pretrained_model_path),
                    streaming_mode="in_memory",
                    normalizer_used=(metrics.get("normalizer_used") if metrics else None),
                    meta_path=config.output.get_output_path(".meta.iic"),
                    motif_category=_post.motif_category,
                    z_excess=_post.z_excess,
                    motif_called_u12=_post.motif_called_u12,
                )

        # Save classification metrics to JSON file (after score adjustment)
        if metrics:
            metrics_path = config.output.get_output_path(".metrics.iic.json")
            messenger.log_only(f"Saving classification metrics to {metrics_path}")
            with open(metrics_path, "w") as f:
                json.dump(metrics, f, indent=2)

            # Display the FINAL (post-discount) classification summary + boundary
            # tables from the finalized metrics, so the console matches
            # metrics.iic.json and the streaming path (vs the first-pass call that
            # the objects still carry at this point).
            _u12_total = metrics.get("u12_count", 0)
            _total_scored = metrics.get("total_scored", _u12_total + metrics.get("u2_count", 0))
            messenger.print_classification_results(
                total=_total_scored,
                u12_count=_u12_total,
                u2_count=_total_scored - _u12_total,
                atac_count=metrics.get("u12_boundaries", {}).get("AT-AC", 0),
                threshold=config.scoring.threshold,
            )
            for label, key in (("U12-type", "u12_boundaries"),
                               ("U2-type", "u2_boundaries")):
                bnd = metrics.get(key) or {}
                if bnd:
                    messenger.print_dinucleotide_boundaries(
                        intron_type=label,
                        boundaries=sorted(bnd.items(), key=lambda x: (-x[1], x[0])),
                        top_n=20,
                    )

        # Calculate and log total runtime
        elapsed_seconds = time.time() - start_time
        hours, remainder = divmod(int(elapsed_seconds), 3600)
        minutes, seconds = divmod(remainder, 60)

        if hours > 0:
            runtime_str = f"{hours}h {minutes}m {seconds}s"
        elif minutes > 0:
            runtime_str = f"{minutes}m {seconds}s"
        else:
            runtime_str = f"{seconds}s"

        messenger.success(f"Pipeline complete! (Runtime: {runtime_str})")

    except Exception as e:
        logger.exception("Pipeline failed with error")
        messenger.error(f"Pipeline failed: {str(e)}")
        raise


def main_test(args):
    """Run installation test with bundled test data.

    Args:
        args: Parsed arguments from argparse
    """
    from rich.console import Console

    console = Console()

    # Find bundled test data.
    #
    # Human chr19 (Ensembl 91, ~58 Mb, ~12k introns) is the smoke fixture.
    # Running adaptive RobustScaler on chr19's introns alone compresses
    # the projection axis (chr19 under-samples the genome-wide U2 IQR),
    # so the multi-bandwidth valley metric under-reports even when the
    # U12 cluster is real. We work around that by passing
    # --load-normalizer pointed at the bundled multispecies scaler
    # (v3_fallback_normalizer.pkl) — the broad-based scaler keeps the
    # U12 lower tail separated from the U2 right tail, valley fires
    # cleanly at depth ~0.79.
    data_dir = Path(__file__).parent.parent / "data" / "test_data"
    genome_file = data_dir / "Homo_sapiens.Chr19.Ensembl_91.fa.gz"
    annotation_file = data_dir / "Homo_sapiens.Chr19.Ensembl_91.gff3.gz"
    bundled_normalizer = (
        Path(__file__).parent.parent / "data" / "v3_fallback_normalizer.pkl"
    )

    # Check if test data exists
    if not genome_file.exists() or not annotation_file.exists():
        console.print(
            "[red]Error: Bundled test data not found![/red]", style="bold"
        )
        console.print(f"Expected location: {data_dir}")
        console.print("\nTest data should include:")
        console.print(f"  - {genome_file.name}")
        console.print(f"  - {annotation_file.name}")
        return 1

    # Show test data location
    console.print("\n[bold cyan]intronIC Installation Test[/bold cyan]")
    console.print(f"Test data location: [green]{data_dir}[/green]")
    console.print(f"  Genome:     {genome_file.name}")
    console.print(f"  Annotation: {annotation_file.name}")

    # If --show-only, exit here
    if getattr(args, "show_only", False):
        console.print("\n[bold]To run test manually:[/bold]")
        console.print(f"  intronIC classify -g {genome_file} \\")
        console.print(f"                    -a {annotation_file} \\")
        console.print(f"                    -n homo_sapiens_chr19 \\")
        console.print(f"                    --load-normalizer {bundled_normalizer} -p 4")
        return 0

    # Run quick classification test
    console.print("\n[bold]Running classification test...[/bold]")

    # Use temporary directory or user-specified output directory
    if args.output_dir:
        output_dir = args.output_dir
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        import tempfile
        output_dir = Path(tempfile.mkdtemp(prefix="intronic_test_"))

    console.print(f"Output directory: [green]{output_dir}[/green]")

    # Build command-line args for classify mode.
    # --load-normalizer wires the bundled multispecies scaler in so the
    # adaptive pre-pass is skipped — chr19 alone is too narrow to fit
    # adaptive without compressing the U2 IQR.
    test_args = [
        "classify",
        "-g", str(genome_file),
        "-a", str(annotation_file),
        "-n", "homo_sapiens_chr19",
        "-o", str(output_dir),
        "-p", str(args.processes),
        "--load-normalizer", str(bundled_normalizer),
    ]

    # Run classification
    start_time = time.time()
    try:
        main(test_args)
        elapsed = time.time() - start_time

        # Check results
        metrics_file = output_dir / "homo_sapiens_chr19.metrics.iic.json"
        if metrics_file.exists():
            # Parse metrics to get counts
            import json
            total_introns = "?"
            u12_introns = "?"
            try:
                with open(metrics_file) as f:
                    metrics = json.load(f)
                    # Use total_scored (introns after filtering) not total_introns_generated
                    total_introns = metrics.get("total_scored", "?")
                    # Use high_confidence_u12 (matches what's shown in results table)
                    u12_introns = metrics.get("high_confidence_u12", "?")
            except Exception:
                pass

            console.print(f"\n[bold green]✓ Test completed successfully![/bold green]")
            console.print(f"  Runtime: {elapsed:.1f}s")
            console.print(f"  Total introns: {total_introns:,}" if isinstance(total_introns, int) else f"  Total introns: {total_introns}")
            console.print(f"  U12 introns: {u12_introns:,}" if isinstance(u12_introns, int) else f"  U12 introns: {u12_introns}")

            if not args.output_dir:
                console.print(f"\n[dim]Output saved to: {output_dir}[/dim]")
        else:
            console.print("[yellow]Warning: Test completed but results file not found[/yellow]")

    except Exception as e:
        console.print(f"[red]Test failed: {str(e)}[/red]")
        raise


def main(args=None):
    """Main entry point for intronIC CLI.

    Args:
        args: Optional list of arguments (defaults to sys.argv)

    Returns:
        Exit code (0 for success, 1 for failure)
    """
    try:
        # Parse arguments
        parser = IntronICArgumentParser()
        parsed_args = parser.parse_args(args)

        # Handle --generate-config
        if getattr(parsed_args, "generate_config", False):
            # Copy the built-in config to current directory
            install_dir = Path(__file__).parent.parent.parent.parent
            source_config = install_dir / "config" / "config.yaml"
            dest_config = Path(".intronIC.yaml")

            if dest_config.exists():
                print(f"Config file already exists: {dest_config}")
                print("Remove it first if you want to regenerate.")
                return 1

            shutil.copy(source_config, dest_config)
            print(f"Generated configuration file: {dest_config}")
            print("Edit this file to customize your intronIC settings.")
            return 0

        # Route to appropriate entry point based on command
        command = getattr(parsed_args, "command", "classify")

        # Test mode doesn't need config - handle it directly
        if command == "test":
            return main_test(parsed_args)

        # Create unified configuration from YAML config + CLI args
        # CLI args take precedence over YAML values
        # This creates a single source of truth for the entire run
        config = IntronICConfig.from_yaml_and_args(parsed_args)

        if command == "train":
            # Train mode: Train model on reference sequences
            main_train(config)
        elif command == "extract":
            # Extract mode: Extract intron sequences without classification
            main_extract(config)
        elif command == "classify":
            # Classify mode: Standard classification pipeline
            main_classify(config)
        else:
            raise ValueError(f"Unknown command: {command}")

        return 0

    except KeyboardInterrupt:
        print("\n\nPipeline interrupted by user", file=sys.stderr)
        return 130  # Standard exit code for SIGINT

    except Exception as e:
        print(f"\n\nError: {str(e)}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())


if __name__ == "__main__":
    sys.exit(main())
if __name__ == "__main__":
    sys.exit(main())
if __name__ == "__main__":
    sys.exit(main())
