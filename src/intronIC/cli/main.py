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
from dataclasses import replace
from multiprocessing import Pool
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

import joblib
import numpy as np
from smart_open import open as smart_open  # type: ignore[import-unresolved]

from intronIC.classification.classifier import IntronClassifier

# Import pipeline components
from intronIC.core.intron import Intron, IntronSequences, OmissionReason
from intronIC.extraction.annotator import AnnotationHierarchyBuilder
from intronIC.extraction.filters import IntronFilter, prefilter_introns
from intronIC.extraction.intronator import IntronGenerator
from intronIC.extraction.sequences import SequenceExtractor
from intronIC.file_io.genome import GenomeReader
from intronIC.file_io.parsers import BEDParser, SequenceParser
from intronIC.file_io.writers import BEDWriter, MetaWriter, ScoreWriter, SequenceWriter
from intronIC.scoring.normalizer import ScoreNormalizer
from intronIC.scoring.scorer import IntronScorer
from intronIC.utils.coordinates import GenomicCoordinate
from intronIC.utils.metadata import (
    generate_pretrained_metadata,
    generate_training_metadata,
    write_metadata,
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

        # Extract coefficients from CalibratedClassifierCV -> Pipeline -> LinearSVC
        try:
            # Navigate nested structure
            if hasattr(model, "calibrated_classifiers_"):
                calibrated_clf = model.calibrated_classifiers_[0]
                fitted_estimator = calibrated_clf.estimator

                # Get calibration method
                calib_method = getattr(model, "method", "unknown")

                # Get LinearSVC from pipeline
                linear_svc = None
                if hasattr(fitted_estimator, "named_steps"):
                    for step_name in ["svc", "linearsvc", "classifier"]:
                        if step_name in fitted_estimator.named_steps:
                            linear_svc = fitted_estimator.named_steps[step_name]
                            break

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

    # Create Rich console for log file with ANSI color support
    # force_terminal=True ensures ANSI codes are written even to files
    # highlight=False disables automatic syntax highlighting to avoid false matches (e.g. "PUT" in "OUTPUT")
    log_console = Console(
        file=open(log_file, "w", encoding="utf-8"),
        force_terminal=True,
        width=120,
        legacy_windows=False,
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
    contig_name: str, contig_introns: List[Intron], flank_len: int, u12_correction: bool
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
                intron, correction_enabled=True, use_strict_motif=True
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


def _init_streaming_worker(
    genome_path: str,
    db_path: str,
    species_name: str,
    simple_name: bool,
    no_abbreviate: bool,
    five_coords: tuple[int, int],
    bp_coords: tuple[int, int],
    three_coords: tuple[int, int],
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

    _streaming_worker_db_path = db_path
    _streaming_worker_species_name = species_name
    _streaming_worker_simple_name = simple_name
    _streaming_worker_no_abbreviate = no_abbreviate
    _streaming_worker_five_coords = five_coords
    _streaming_worker_bp_coords = bp_coords
    _streaming_worker_three_coords = three_coords


def _process_contig_streaming_worker(
    worker_input: tuple[str, List[Intron], int, bool],
) -> tuple[str, List[Intron], int, int]:
    """
    Worker function for parallel streaming contig processing.

    Like _process_contig_worker but also:
    1. Writes full sequences to SQLite (for later output)
    2. Extracts scoring motifs
    3. Returns lightweight introns (no full sequences)

    Args:
        worker_input: Tuple of (contig_name, contig_introns, flank_len, u12_correction)

    Returns:
        Tuple of (contig_name, lightweight_introns, sequences_stored_count, corrections_count)
    """
    contig_name, contig_introns, flank_len, u12_correction = worker_input
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
                intron, correction_enabled=True, use_strict_motif=True
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
) -> List[Intron]:
    """
    Extract introns contig-by-contig with pre-filtering and deduplication.

    This function implements the memory-optimized extraction pipeline:
    1. Parse annotations (coordinates only)
    2. Pre-filter before extraction (removes ~85-90% of introns)
    3. Process one contig at a time:
       - Extract sequences with deduplication
       - Apply U12 corrections
       - Free contig memory before next
    4. Return all introns WITH sequences for scoring

    Memory savings from pre-filtering: 28 GB → ~4-5 GB peak (82-85% reduction)

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

    messenger.log_only(
        f"Pre-filter results: "
        f"extracting sequences for {len(extract_list):,} ({len(extract_list) / len(introns_list) * 100:.1f}%), "
        f"skipping {len(skip_list):,} ({len(skip_list) / len(introns_list) * 100:.1f}%)"
    )

    # Free original list
    del introns_list
    gc.collect()

    # Step 4: Extract sequences (delegate to helper for code reuse)
    all_introns = _extract_sequences_for_introns(
        extract_list, skip_list, config, messenger, reporter
    )

    return all_introns


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

    # Group introns by contig for contig-by-contig processing
    from collections import defaultdict

    introns_by_contig = defaultdict(list)
    for intron in extract_list:
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

    # Accumulator for all introns
    all_introns = []
    all_introns.extend(skip_list)  # Add skipped introns (no sequences)

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
                        intron, correction_enabled=True, use_strict_motif=True
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
    Extract introns in streaming mode (~85% memory savings).

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
    messenger.info("Using streaming mode (~85% memory savings)")

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
                        intron, correction_enabled=True, use_strict_motif=True
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
        intron = Intron(intron_id=intron_id, coordinates=coord)
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
                intron, correction_enabled=True, use_strict_motif=True
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
    from intronIC.core.intron import (
        Intron,
        IntronFlags,
        IntronMetadata,
        IntronScores,
        IntronSequences,
    )
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

        # Create IntronSequences with the loaded data
        # Populate five_seq and three_seq to satisfy has_sequences() check
        # These are approximations since we don't know exact scoring regions,
        # but the scorer will extract the proper regions during scoring
        sequences = IntronSequences(
            seq=seq_line.sequence,
            upstream_flank=seq_line.upstream_flank,
            downstream_flank=seq_line.downstream_flank,
            # Populate with dummy values to pass has_sequences() check
            # The actual values used for scoring are extracted on-the-fly
            five_seq=seq[:10] if seq and len(seq) >= 10 else seq,
            three_seq=seq[-10:] if seq and len(seq) >= 10 else seq,
            bp_seq=None,  # Will be found during scoring
            bp_region_seq=None,  # Will be extracted during scoring
            five_prime_dnt=five_prime_dnt,
            three_prime_dnt=three_prime_dnt,
        )

        # Create IntronScores if score is present
        scores = None
        if seq_line.score is not None:
            scores = IntronScores(svm_score=seq_line.score)

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
            flags=IntronFlags.SEQUENCE_ONLY | IntronFlags.LONGEST_ISOFORM,
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
) -> List[Intron]:
    """Score introns with PWM matrices (parallel or sequential).

    PHASE 3: Restored parallel scoring from original intronIC v1.5.1
    Uses multiprocessing.Pool when config.performance.processes > 1

    Args:
        introns: List of introns to score
        config: Pipeline configuration
        messenger: Unified messenger for output
        reporter: Progress reporter (for progress bars only)

    Returns:
        List of introns with raw scores
    """
    messenger.info("Loading PWM matrices")

    # Load PWM matrices with fallback support for custom matrices
    pwm_sets = load_pwms_with_fallback(config, messenger)

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

        # Get PWM file paths for workers (they load PWMs in separate processes)
        data_dir = Path(__file__).parent.parent / "data"
        default_pwm_file = data_dir / "intronIC_scoring_PWMs.json"
        custom_pwm_file = config.scoring.pwm_file

        progress = reporter.create_progress()
        scored_introns = []
        failed_count = 0

        with progress:
            task = progress.add_task("[cyan]Scoring introns...", total=len(introns))

            # Use initializer to load PWMs once per worker process (not per intron)
            with Pool(
                processes=n_workers,
                initializer=_init_worker,
                initargs=(
                    default_pwm_file,
                    custom_pwm_file,
                    five_coords,
                    bp_coords,
                    three_coords,
                    ignore_nc_dnts,
                    config.scoring.pseudocount,
                ),
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


def normalize_scores(
    introns: List[Intron],
    config: IntronICConfig,
    messenger: "UnifiedMessenger",
    reporter: IntronICProgressReporter,
) -> Tuple[List[Intron], List[Intron], List[Intron], "ScoreNormalizer"]:
    """Normalize intron scores using z-score transformation.

    Args:
        introns: List of introns with raw scores
        config: Pipeline configuration
        messenger: Unified messenger for console and log output
        reporter: Progress reporter

    Returns:
        Tuple of (normalized experimental introns, u12 reference, u2 reference, normalizer)
    """
    messenger.info("Loading reference sequences")

    # Load reference data - use custom paths if provided, otherwise use defaults
    data_dir = Path(__file__).parent.parent / "data"

    if config.scoring.reference_u12s:
        u12_file = config.scoring.reference_u12s
    else:
        u12_file = data_dir / "u12_reference.introns.iic.gz"

    if config.scoring.reference_u2s:
        u2_file = config.scoring.reference_u2s
    else:
        u2_file = data_dir / "u2_reference.introns.iic.gz"

    if not u12_file.exists() or not u2_file.exists():
        raise FileNotFoundError(
            f"Reference data not found. U12: {u12_file}, U2: {u2_file}"
        )

    u12_reference = load_reference_sequences(u12_file, messenger=messenger)
    u2_reference = load_reference_sequences(u2_file, messenger=messenger)

    messenger.log_only(
        f"Loaded {len(u12_reference)} U12-type and {len(u2_reference)} U2-type reference introns"
    )

    # Score reference introns
    messenger.info("Scoring reference sequences")

    # Load PWM matrices (U2 BP matrix now included in main file)
    pwm_sets = load_pwms_with_fallback(config, messenger)

    scorer = IntronScorer(
        pwm_sets=pwm_sets,
        five_coords=(
            config.scoring.scoring_regions.five_start,
            config.scoring.scoring_regions.five_end,
        ),
        bp_coords=(
            config.scoring.scoring_regions.bp_start,
            config.scoring.scoring_regions.bp_end,
        ),
        three_coords=(
            config.scoring.scoring_regions.three_start,
            config.scoring.scoring_regions.three_end,
        ),
        ignore_nc_dnts=config.scoring.ignore_nc_dnts,
    )

    # Score reference introns
    u12_scored = [scorer.score_intron(i) for i in u12_reference]
    u2_scored = [scorer.score_intron(i) for i in u2_reference]
    reference_introns = u12_scored + u2_scored

    # Normalize scores
    messenger.info("Normalizing scores with z-score transformation")
    messenger.log_only("Fitting normalizer on reference data")

    normalizer = ScoreNormalizer()
    normalizer.fit(reference_introns, dataset_type="reference")

    # Transform reference introns (needed for classification)
    # Note: transform() takes an iterable and returns an iterator, so we materialize with list()
    u12_normalized = list(normalizer.transform(u12_scored, dataset_type="reference"))
    u2_normalized = list(normalizer.transform(u2_scored, dataset_type="reference"))

    # Transform experimental introns
    normalized_introns = list(
        normalizer.transform(introns, dataset_type="experimental")
    )

    messenger.log_only(f"Normalized {len(normalized_introns)} experimental introns")
    messenger.log_only(
        f"Normalized {len(u12_normalized)} U12-type and {len(u2_normalized)} U2-type reference introns"
    )
    return normalized_introns, u12_normalized, u2_normalized, normalizer


def _apply_margin_alignment(
    classified_introns: List[Intron],
    ensemble: "SVMEnsemble",
    human_margin_stats: dict,
    threshold: float,
    messenger: "UnifiedMessenger",
) -> List[Intron]:
    """
    Apply margin-space alignment to classified introns (post-processing approach).

    This function re-calibrates probabilities by aligning target species margins
    to the human U2 margin distribution, then re-applying Platt calibration with
    adjusted parameters.

    Args:
        classified_introns: Introns with initial probabilities from ensemble
        ensemble: Trained SVM ensemble
        human_margin_stats: Dict with 'margin_median', 'margin_iqr', 'n_samples'
        threshold: Classification threshold (for diagnostics)
        messenger: For logging diagnostics

    Returns:
        Introns with margin-aligned probabilities

    Note:
        This is mathematically equivalent to modifying the predictor, but simpler
        to implement as a post-processing step.
    """
    from dataclasses import replace

    from scipy.special import expit as sigmoid

    from intronIC.scoring.margin_alignment import (
        compute_margin_alignment,
        fold_alignment_into_platt,
    )

    messenger.log_only("Computing margin alignment parameters...")

    # Step 1: Extract margins from classified introns
    # The ensemble stores calibrated probabilities, we need to reverse-engineer margins
    # Alternative: Re-compute margins from features
    # We'll use the second approach for accuracy

    # Extract features from introns (require scores to be present)
    features = np.array(
        [
            [i.scores.five_z_score, i.scores.bp_z_score, i.scores.three_z_score]
            for i in classified_introns
            if i.scores is not None
        ]
    )

    # Get first model's base pipeline (for margin computation)
    first_model = ensemble.models[0].model
    if hasattr(first_model, "calibrated_classifiers_"):
        base_pipeline = first_model.calibrated_classifiers_[0].estimator
    else:
        base_pipeline = first_model

    # Compute margins (before calibration)
    target_margins = base_pipeline.decision_function(features)

    # Step 2: Compute alignment parameters
    scale, shift, align_stats = compute_margin_alignment(
        target_margins=target_margins,
        human_median=human_margin_stats["margin_median"],
        human_iqr=human_margin_stats["margin_iqr"],
        human_n=human_margin_stats["n_samples"],
    )

    # Log diagnostics
    messenger.log_only(align_stats.format_summary())

    # Step 3: Check calibration method and apply alignment
    if hasattr(first_model, "calibrated_classifiers_"):
        # For CalibratedClassifierCV with cv='prefit', there's one calibrated classifier
        calibrator = first_model.calibrated_classifiers_[0].calibrators[0]

        # Check if it's Platt (sigmoid) or isotonic calibration
        if hasattr(calibrator, "a_"):
            # Platt/sigmoid calibration - can fold alignment into parameters
            A = calibrator.a_
            B = calibrator.b_

            # Step 4: Fold alignment into Platt parameters
            A_new, B_new = fold_alignment_into_platt(A, B, scale, shift)

            messenger.log_only(f"Platt parameters: A={A:.4f}, B={B:.4f}")
            messenger.log_only(f"Aligned Platt:    A'={A_new:.4f}, B'={B_new:.4f}")

            # Step 5: Re-calibrate probabilities with aligned parameters
            # p_new = sigmoid(A_new * margin + B_new)
            probs_new = (
                sigmoid(A_new * target_margins + B_new) * 100.0
            )  # Convert to 0-100
        else:
            # Isotonic regression - cannot fold alignment, apply transform then re-calibrate
            messenger.log_only(
                "Model uses isotonic calibration - applying alignment then re-calibrating"
            )

            # Transform margins
            aligned_margins = scale * target_margins + shift

            # Re-calibrate with isotonic regressor
            probs_new = calibrator.predict(aligned_margins) * 100.0  # Convert to 0-100
    else:
        # Shouldn't happen, but fallback
        messenger.warning(
            "Could not extract calibration parameters - skipping alignment"
        )
        return classified_introns

    # Step 6: Update introns with new probabilities
    updated_introns = []
    for intron, new_prob in zip(classified_introns, probs_new):
        if intron.scores is None or intron.metadata is None:
            # Skip introns without scores/metadata (shouldn't happen but be defensive)
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
    n_u12_before = sum(
        1
        for i in classified_introns
        if i.metadata is not None and i.metadata.type_id == "u12"
    )
    n_u12_after = sum(
        1
        for i in updated_introns
        if i.metadata is not None and i.metadata.type_id == "u12"
    )
    messenger.log_only(
        f"Margin alignment effect: {n_u12_before} → {n_u12_after} U12 predictions "
        f"({100 * n_u12_after / len(updated_introns):.3f}%)"
    )

    return updated_introns


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

# Globals for streaming classification worker initialization
_streaming_classify_worker_genome_path: str = ""
_streaming_classify_worker_annotation_db_path: str = ""
_streaming_classify_worker_ensemble: Any = None
_streaming_classify_worker_scaler: Any = None
_streaming_classify_worker_pwm_sets: Any = None
_streaming_classify_worker_config: dict = {}


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

    # Store classification components in globals
    global _streaming_classify_worker_genome_path
    global _streaming_classify_worker_annotation_db_path
    global _streaming_classify_worker_ensemble
    global _streaming_classify_worker_scaler
    global _streaming_classify_worker_pwm_sets
    global _streaming_classify_worker_config

    _streaming_classify_worker_genome_path = genome_path
    _streaming_classify_worker_annotation_db_path = annotation_db_path
    _streaming_classify_worker_ensemble = ensemble
    _streaming_classify_worker_scaler = scaler
    _streaming_classify_worker_pwm_sets = pwm_sets
    _streaming_classify_worker_config = config_dict


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
    contig, contig_annotation_count = contig_input
    from collections import Counter

    from intronIC.classification.predictor import classify_introns_batch
    from intronIC.extraction.annotator import AnnotationHierarchyBuilder
    from intronIC.extraction.filters import IntronFilter
    from intronIC.extraction.intronator import IntronGenerator
    from intronIC.extraction.sequences import SequenceExtractor
    from intronIC.file_io.annotation_store import StreamingAnnotationStore
    from intronIC.file_io.indexed_genome import get_worker_genome
    from intronIC.scoring.scorer import IntronScorer, score_and_normalize_batch

    # Access worker globals
    config = _streaming_classify_worker_config
    ensemble = _streaming_classify_worker_ensemble
    scaler = _streaming_classify_worker_scaler
    pwm_sets = _streaming_classify_worker_pwm_sets

    # Statistics for this contig
    from intronIC.extraction.filters import FilterStats

    stats = {
        "genes": 0,
        "introns_generated": 0,
        "scored": 0,
        "classified": 0,
        "boundaries_u12": Counter(),
        "boundaries_u2": Counter(),
        "filter_stats": FilterStats(),  # Initialize empty, will be updated if filtering occurs
        "duplicate_map": {},  # Initialize empty, will be updated if filtering occurs
        "overlap_map": {},  # Initialize empty, will be updated if filtering occurs
    }

    # Open annotation store (read-only, each worker gets own connection)
    annotation_store = StreamingAnnotationStore(
        _streaming_classify_worker_annotation_db_path
    )

    # Get annotations for this contig
    contig_annotations = annotation_store.get_annotations_for_contig(contig)
    annotation_store.close()

    # Skip empty contigs
    if not contig_annotations or all(
        a.feat_type == "region" for a in contig_annotations
    ):
        return contig, [], [], stats

    # Build gene hierarchy
    builder = AnnotationHierarchyBuilder(
        child_features=["cds", "exon"],
        clean_names=config["clean_names"],
        messenger=None,  # No logging from workers
    )

    try:
        contig_genes = builder.build_from_annotations(contig_annotations)
    except ValueError as e:
        if "Could not establish parent-child relationships" in str(e):
            return contig, [], [], stats
        raise

    del contig_annotations

    if not contig_genes:
        return contig, [], [], stats

    stats["genes"] = len(contig_genes)

    # Generate introns
    generator = IntronGenerator(debug=config["debug"], messenger=None)
    contig_introns = list(
        generator.generate_from_genes(contig_genes, builder.feature_index)
    )

    # Filter by feature type
    if config["feature_type"] == "cds":
        contig_introns = [
            i
            for i in contig_introns
            if i.metadata is not None and i.metadata.defined_by == "cds"
        ]
    elif config["feature_type"] == "exon":
        contig_introns = [
            i
            for i in contig_introns
            if i.metadata is not None and i.metadata.defined_by == "exon"
        ]

    stats["introns_generated"] = len(contig_introns)

    if not contig_introns:
        del contig_genes, builder
        return contig, [], [], stats

    # Get worker's genome reader and create extractor
    indexed_genome = get_worker_genome()
    extractor = SequenceExtractor.__new__(SequenceExtractor)
    extractor.genome_file = None  # type: ignore[assignment]
    extractor.genome_reader = indexed_genome  # type: ignore[assignment]
    extractor.use_cache = False

    # Extract sequences
    contig_with_seqs = list(
        extractor.extract_sequences_with_deduplication(
            contig_introns, flank_size=config["flank_len"]
        )
    )

    # Apply U12 corrections if enabled
    if config["u12_boundary_correction"]:
        from intronIC.extraction.boundary_correction import correct_intron_if_needed

        corrected_contig_introns = []
        for intron in contig_with_seqs:
            corrected_intron, was_corrected = correct_intron_if_needed(
                intron, correction_enabled=True, use_strict_motif=True
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

    # Create filter
    intron_filter = IntronFilter(
        min_length=config["min_intron_len"],
        bp_matrix_length=7,
        scoring_regions=["five", "three"],
        allow_noncanonical=not config["exclude_noncanonical"],
        allow_overlap=not config["no_intron_overlap"],
        longest_only=True,
        include_duplicates=config["include_duplicates"],
    )

    # Filter introns
    filtered_introns = intron_filter.filter_introns(contig_with_seqs)

    # Capture filter statistics and maps
    stats["filter_stats"] = intron_filter.stats
    stats["duplicate_map"] = intron_filter.get_duplicate_map()
    stats["overlap_map"] = intron_filter.get_overlap_map()

    # Filter to only scorable introns
    scorable = [
        i
        for i in filtered_introns
        if i.has_sequences
        and (i.metadata is None or i.metadata.omitted == OmissionReason.NONE)
    ]

    if not scorable:
        del contig_genes, builder, contig_introns, contig_with_seqs
        # Return filtered_introns even when no scorable - these are omitted introns
        return contig, [], filtered_introns, stats

    # Create scorer
    scorer = IntronScorer(
        pwm_sets=pwm_sets,
        five_coords=(config["five_start"], config["five_end"]),
        bp_coords=(config["bp_start"], config["bp_end"]),
        three_coords=(config["three_start"], config["three_end"]),
        ignore_nc_dnts=config["ignore_nc_dnts"],
    )

    # Score and normalize
    scored_introns = score_and_normalize_batch(scorable, scorer, scaler)
    stats["scored"] = len(scored_introns)

    # Classify
    classified_introns = classify_introns_batch(
        scored_introns,
        ensemble,
        threshold=config["threshold"],
    )
    stats["classified"] = len(classified_introns)

    # Track boundary statistics
    for intron in classified_introns:
        if (
            intron.metadata
            and intron.sequences
            and intron.sequences.terminal_dinucleotides
        ):
            dnts = intron.sequences.terminal_dinucleotides
            if (
                intron.scores
                and intron.scores.svm_score is not None
                and intron.scores.svm_score >= config["threshold"]
            ):
                stats["boundaries_u12"][dnts] += 1
            else:
                stats["boundaries_u2"][dnts] += 1

    # Clean up
    del contig_genes, builder, contig_introns, contig_with_seqs
    del scorable, scored_introns

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
    from collections import Counter

    from intronIC.classification.predictor import classify_introns_batch
    from intronIC.extraction.intronator import IntronGenerator
    from intronIC.file_io.annotation_store import StreamingAnnotationStore
    from intronIC.file_io.indexed_genome import IndexedGenomeReader
    from intronIC.file_io.parsers import BioGLAnnotationParser
    from intronIC.file_io.writers import MappingWriter, StreamingOutputWriter
    from intronIC.scoring.scorer import score_and_normalize_batch

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
    model_data = joblib.load(config.training.pretrained_model_path)

    # Handle both old and new model format
    if isinstance(model_data, dict):
        ensemble = model_data["ensemble"]
        saved_normalizer = model_data.get("normalizer", None)
    else:
        ensemble = model_data
        saved_normalizer = None

    messenger.log_only(f"Loaded ensemble with {len(ensemble.models)} models")

    # Extract frozen scaler from normalizer
    if saved_normalizer is None:
        raise ValueError(
            "Model bundle does not contain a normalizer. "
            "Cannot use true streaming mode without frozen scaler. "
            "Use standard --streaming mode instead."
        )

    scaler = saved_normalizer.get_frozen_scaler()
    messenger.log_only("Extracted frozen scaler from model normalizer")

    # Load PWM matrices
    messenger.info("Loading PWM matrices")
    pwm_sets = load_pwms_with_fallback(config, messenger)

    # Configure scorer
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

    scorer = IntronScorer(
        pwm_sets=pwm_sets,
        five_coords=five_coords,
        bp_coords=bp_coords,
        three_coords=three_coords,
        ignore_nc_dnts=config.scoring.ignore_nc_dnts,
    )

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
    boundaries_u12: Counter = Counter()
    boundaries_u2: Counter = Counter()

    # Accumulated filter statistics
    from intronIC.extraction.filters import FilterStats

    accumulated_filter_stats = FilterStats()

    # Accumulated duplicate and overlap maps
    accumulated_duplicate_map: Dict[str, Set[str]] = {}
    accumulated_overlap_map: Dict[str, Set[str]] = {}

    # Get contig lengths for length-weighted progress reporting
    import numpy as np

    from intronIC.file_io.indexed_genome import get_contig_lengths

    assert config.input.genome is not None, "Genome path required"
    contig_lengths = get_contig_lengths(config.input.genome)

    # Prepare length-weighted progress tracking
    contig_names = [contig for contig, _ in contigs_with_counts]
    contig_length_list = [contig_lengths[c] for c in contig_names]
    cumulative_lengths = np.cumsum(contig_length_list)
    total_length = cumulative_lengths[-1]

    # Determine parallelization
    n_processes = config.performance.processes
    use_parallel = n_processes > 1 and len(contigs_with_counts) > 1

    if use_parallel:
        messenger.info(
            f"Processing {len(contigs_with_counts)} contigs in parallel "
            f"({n_processes} processes)"
        )

        # Build config dict for workers
        worker_config = {
            "clean_names": config.output.clean_names,
            "debug": config.output.debug,
            "feature_type": config.extraction.feature_type,
            "flank_len": config.extraction.flank_len,
            "u12_boundary_correction": config.extraction.u12_boundary_correction,
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

        # Close annotation store in main process - workers will open their own
        annotation_store.close()

        # Prepare worker inputs and create contig index mapping
        worker_inputs = [(contig, count) for contig, count in contigs_with_counts]
        contig_to_index = {
            contig: idx for idx, (contig, _) in enumerate(contigs_with_counts)
        }

        # Process contigs in parallel with progress bar
        all_classified_introns = []
        all_filtered_introns = []  # Includes omitted introns for output files
        completed = 0
        completed_length = 0

        # Create progress bar
        progress = reporter.create_progress()

        with Pool(
            processes=n_processes,
            initializer=_init_streaming_classify_worker,
            initargs=(
                str(config.input.genome),
                str(annotation_db_path),
                ensemble,
                scaler,
                pwm_sets,
                worker_config,
            ),
        ) as pool, progress:
            # Add task with total = total genome length for smooth progress
            task = progress.add_task(
                "[cyan]Classifying introns (parallel streaming)...", total=total_length
            )

            try:
                # Use imap_unordered to get results as soon as any worker completes
                # This provides better visual feedback - small contigs update progress immediately
                for contig_name, classified_introns, filtered_introns, stats in pool.imap_unordered(
                    _process_contig_streaming_classify_worker, worker_inputs
                ):
                    completed += 1
                    all_classified_introns.extend(classified_introns)
                    all_filtered_introns.extend(filtered_introns)

                    # Accumulate statistics
                    total_genes += stats["genes"]
                    total_introns_generated += stats["introns_generated"]
                    total_scored += stats["scored"]
                    total_classified += stats["classified"]
                    boundaries_u12.update(stats["boundaries_u12"])
                    boundaries_u2.update(stats["boundaries_u2"])

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

                    # Update progress bar based on this contig's length
                    contig_idx = contig_to_index[contig_name]
                    contig_length = contig_length_list[contig_idx]
                    completed_length += contig_length
                    progress.update(task, completed=completed_length)
            except Exception as e:
                messenger.error(
                    f"Parallel streaming classification failed after "
                    f"{completed}/{len(contigs_with_counts)} contigs"
                )
                messenger.error(f"Error: {e}")
                raise

        # Merge classified introns with omitted introns for complete output
        # classified_introns have scores, filtered_introns includes omitted ones without scores
        all_introns_for_output = merge_scored_and_omitted_introns(
            all_classified_introns, all_filtered_introns, messenger
        )

        # Write outputs sequentially (maintains deterministic order by sorting)
        # Sort by coordinates for consistent output order
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

    else:
        # Sequential processing (original behavior)
        messenger.info(f"Processing {len(contigs_with_counts)} contigs sequentially")

        # Create IntronFilter for per-contig filtering
        intron_filter = IntronFilter(
            min_length=config.extraction.min_intron_len,
            bp_matrix_length=7,
            scoring_regions=["five", "three"],
            allow_noncanonical=not config.scoring.exclude_noncanonical,
            allow_overlap=not config.extraction.no_intron_overlap,
            longest_only=True,
            include_duplicates=config.extraction.include_duplicates,
        )

        # Process each contig with progress tracking
        progress = reporter.create_progress()

        with output_writer, progress:
            task = progress.add_task(
                "[cyan]Processing annotations...", total=total_annotations
            )

            for contig_idx, (contig, contig_annotation_count) in enumerate(
                contigs_with_counts, 1
            ):
                # Get annotations for this contig from SQLite
                contig_annotations = annotation_store.get_annotations_for_contig(contig)

                # Skip empty contigs or region-only contigs
                if not contig_annotations or all(
                    a.feat_type == "region" for a in contig_annotations
                ):
                    progress.update(task, advance=contig_annotation_count)
                    continue

                # Build gene hierarchy for this contig
                builder = AnnotationHierarchyBuilder(
                    child_features=["cds", "exon"],
                    clean_names=config.output.clean_names,
                    messenger=messenger,
                )

                try:
                    contig_genes = builder.build_from_annotations(contig_annotations)
                except ValueError as e:
                    if "Could not establish parent-child relationships" in str(e):
                        progress.update(task, advance=contig_annotation_count)
                        continue
                    raise

                del contig_annotations

                if not contig_genes:
                    progress.update(task, advance=contig_annotation_count)
                    continue

                total_genes += len(contig_genes)

                # Generate introns for this contig
                generator = IntronGenerator(
                    debug=config.output.debug, messenger=messenger
                )
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

                total_introns_generated += len(contig_introns)

                if not contig_introns:
                    progress.update(task, advance=contig_annotation_count)
                    del contig_genes, builder
                    gc.collect()
                    continue

                # Create sequence extractor using indexed genome reader
                sequence_extractor = SequenceExtractor.from_indexed_reader(
                    str(config.input.genome), genome_reader
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

                    corrected_contig_introns = []
                    for intron in contig_with_seqs:
                        corrected_intron, was_corrected = correct_intron_if_needed(
                            intron, correction_enabled=True, use_strict_motif=True
                        )
                        if was_corrected:
                            corrected_with_seq = list(
                                sequence_extractor.extract_sequences(
                                    [corrected_intron],
                                    flank_size=config.extraction.flank_len,
                                )
                            )[0]
                            corrected_contig_introns.append(corrected_with_seq)
                        else:
                            corrected_contig_introns.append(corrected_intron)
                    contig_with_seqs = corrected_contig_introns

                # Filter introns (duplicates, short, noncanonical, etc.)
                filtered_introns = intron_filter.filter_introns(contig_with_seqs)

                # Accumulate filter statistics from this contig
                accumulated_filter_stats.duplicates += intron_filter.stats.duplicates
                accumulated_filter_stats.short += intron_filter.stats.short
                accumulated_filter_stats.ambiguous += intron_filter.stats.ambiguous
                accumulated_filter_stats.noncanonical += (
                    intron_filter.stats.noncanonical
                )
                accumulated_filter_stats.overlap += intron_filter.stats.overlap
                accumulated_filter_stats.isoform += intron_filter.stats.isoform
                accumulated_filter_stats.total_introns += (
                    intron_filter.stats.total_introns
                )
                accumulated_filter_stats.kept_introns += (
                    intron_filter.stats.kept_introns
                )

                # Accumulate duplicate and overlap maps
                accumulated_duplicate_map.update(intron_filter.get_duplicate_map())
                accumulated_overlap_map.update(intron_filter.get_overlap_map())

                # Filter to only scorable introns (have sequences and aren't omitted)
                scorable = [
                    i
                    for i in filtered_introns
                    if i.has_sequences
                    and (
                        i.metadata is None or i.metadata.omitted == OmissionReason.NONE
                    )
                ]

                if not scorable:
                    # Still write omitted introns even when no scorable introns
                    for intron in filtered_introns:
                        if intron.metadata and intron.metadata.omitted != OmissionReason.NONE:
                            output_writer.write_intron(intron)
                    progress.update(task, advance=contig_annotation_count)
                    del (
                        contig_genes,
                        builder,
                        contig_introns,
                        contig_with_seqs,
                        filtered_introns,
                    )
                    gc.collect()
                    continue

                # Score and normalize using frozen scaler
                scored_introns = score_and_normalize_batch(scorable, scorer, scaler)
                total_scored += len(scored_introns)

                # Classify
                classified_introns = classify_introns_batch(
                    scored_introns,
                    ensemble,
                    threshold=config.scoring.threshold,
                )
                total_classified += len(classified_introns)

                # Get omitted introns from filtered_introns (those not in scorable)
                scored_ids = {id(intron) for intron in classified_introns}
                omitted_introns = [
                    intron
                    for intron in filtered_introns
                    if id(intron) not in scored_ids
                    and intron.metadata
                    and intron.metadata.omitted != OmissionReason.NONE
                ]

                # Write classified introns and track statistics
                for intron in classified_introns:
                    output_writer.write_intron(intron)

                    # Track boundary statistics
                    if (
                        intron.metadata
                        and intron.sequences
                        and intron.sequences.terminal_dinucleotides
                    ):
                        dnts = intron.sequences.terminal_dinucleotides
                        if (
                            intron.scores
                            and intron.scores.svm_score is not None
                            and intron.scores.svm_score >= config.scoring.threshold
                        ):
                            boundaries_u12[dnts] += 1
                        else:
                            boundaries_u2[dnts] += 1

                # Write omitted introns (those with omission reasons)
                for intron in omitted_introns:
                    output_writer.write_intron(intron)

                # Free contig memory
                del contig_genes, builder, contig_introns, contig_with_seqs
                del filtered_introns, scorable, scored_introns, classified_introns, omitted_introns
                gc.collect()

                progress.update(task, advance=contig_annotation_count)

        del genome_reader
        gc.collect()

    # Free global resources and cleanup temp files
    annotation_store.cleanup()  # Delete SQLite annotation database
    gc.collect()

    # Build summary
    summary = output_writer.get_summary()
    summary.update(
        {
            "total_genes": total_genes,
            "total_introns_generated": total_introns_generated,
            "total_scored": total_scored,
            "pretrained": True,
            "model_path": str(config.training.pretrained_model_path),
            "streaming_mode": "per_contig",
            "u12_boundaries": dict(boundaries_u12.most_common(20)),
            "u2_boundaries": dict(boundaries_u2.most_common(20)),
        }
    )

    messenger.info(
        f"Streaming classification complete: {total_classified:,} introns classified"
    )
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

    # Display classification summary
    messenger.print_classification_results(
        total=total_classified,
        u12_count=output_writer.high_confidence_u12,
        u2_count=total_classified - output_writer.high_confidence_u12,
        atac_count=boundaries_u12.get("AT-AC", 0),
        threshold=config.scoring.threshold,
    )

    # Display boundary statistics
    if boundaries_u12:
        sorted_boundaries = sorted(boundaries_u12.items(), key=lambda x: (-x[1], x[0]))
        messenger.print_dinucleotide_boundaries(
            intron_type="U12-type", boundaries=sorted_boundaries, top_n=20
        )

    if boundaries_u2:
        sorted_boundaries = sorted(boundaries_u2.items(), key=lambda x: (-x[1], x[0]))
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
        }

    messenger.info(f"Loading pretrained model from {model_path}")

    # Load model bundle
    if not model_path.exists():
        raise FileNotFoundError(f"Pretrained model not found: {model_path}")

    model_data = joblib.load(model_path)

    # Handle both old format (SVMEnsemble directly) and new format (dict bundle)
    if isinstance(model_data, dict):
        # New format: {'ensemble': ..., 'normalizer': ..., 'threshold': ...}
        ensemble = model_data["ensemble"]
        saved_normalizer = model_data.get("normalizer", None)
        _saved_threshold = model_data.get("threshold", config.scoring.threshold)  # noqa: F841
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

    # Determine normalizer mode
    normalizer_mode = config.scoring.normalizer_mode
    if normalizer_mode == "auto":
        # Use human scaler if available, otherwise adaptive
        if saved_normalizer is not None:
            normalizer_mode = "human"
            messenger.log_only("Auto mode: Using saved human scaler (recommended)")
        else:
            normalizer_mode = "adaptive"
            messenger.log_only("Auto mode: Falling back to adaptive (no saved scaler)")

    # Apply normalizer strategy
    if normalizer_mode == "human":
        if saved_normalizer is None:
            raise ValueError(
                "Normalizer mode 'human' requested but model has no saved scaler. "
                "Retrain model or use '--normalizer_mode adaptive'."
            )
        messenger.info("Using human-trained normalizer (scaler from training species)")
        messenger.log_only("This preserves composition bias correction across species")
        normalizer = saved_normalizer
    else:  # adaptive
        # Check if user wants to load a saved normalizer
        if config.scoring.load_normalizer is not None:
            messenger.info(
                f"Loading saved normalizer from {config.scoring.load_normalizer}"
            )
            normalizer = joblib.load(config.scoring.load_normalizer)
            messenger.log_only("Using saved normalizer for reproducible normalization")
        else:
            # Fit normalizer on experimental data (feature re-scaling)
            messenger.info(
                "Fitting normalizer on experimental data (domain adaptation)"
            )
            messenger.log_only("Re-normalizing features to target species distribution")
            messenger.log_only(
                f"Fitting normalizer on {len(introns)} experimental introns"
            )
            normalizer = ScoreNormalizer()
            normalizer.fit(introns, dataset_type="unlabeled")

            # Save normalizer if requested
            if config.scoring.save_normalizer:
                normalizer_path = config.output.get_output_path(".normalizer.pkl")
                messenger.info(f"Saving fitted normalizer to {normalizer_path}")
                joblib.dump(normalizer, normalizer_path, compress=3)
                messenger.log_only(
                    "Future runs can reuse this normalizer with --load-normalizer"
                )

    normalized_introns = list(
        normalizer.transform(introns, dataset_type="experimental")
    )
    messenger.log_only(f"Normalized {len(normalized_introns)} experimental introns")

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
    }

    # Generate metadata for pretrained model usage
    run_metadata_path = config.output.get_output_path(".run_metadata.json")
    messenger.log_only("Recording pretrained model usage")
    run_metadata = generate_pretrained_metadata(
        model_path=model_path, threshold=config.scoring.threshold
    )
    write_metadata(run_metadata, run_metadata_path)
    messenger.log_only(f"Run metadata saved to {run_metadata_path}")

    return classified_introns, metrics


def classify_introns(
    introns: List[Intron],
    u12_reference: List[Intron],
    u2_reference: List[Intron],
    normalizer: "ScoreNormalizer",
    config: IntronICConfig,
    messenger: "UnifiedMessenger",
    reporter: IntronICProgressReporter,
) -> Tuple[List[Intron], dict]:
    """Classify introns as U2 or U12 type.

    Args:
        introns: List of introns with z-scores
        u12_reference: U12 reference introns (scored and normalized)
        u2_reference: U2 reference introns (scored and normalized)
        normalizer: Fitted normalizer used for score transformation
        config: Pipeline configuration
        messenger: Unified messenger for console and log output
        reporter: Progress reporter

    Returns:
        Tuple of (classified introns, classification metrics)
    """
    messenger.info("Training SVM classifier")

    # All optimizer/ensemble settings are pre-loaded in config.optimizer and config.ensemble
    # These were populated from YAML during IntronICConfig.from_yaml_and_args()
    # No need to re-load the YAML file here!

    opt = config.optimizer
    ens = config.ensemble

    # Log config source
    if config.config_path:
        messenger.log_only(f"Using configuration from: {config.config_path}")
    else:
        messenger.log_only("Using default configuration")

    # Log optimizer settings
    messenger.log_only("Optimizer configuration:")
    messenger.log_only(f"  Optimization rounds: {opt.n_rounds}")
    messenger.log_only(f"  Initial grid points: {opt.n_points_initial}")
    messenger.log_only(f"  CV folds: {opt.cv_folds}")
    messenger.log_only(f"  Parallel jobs: {opt.n_jobs}")
    messenger.log_only(f"  Max iterations: {opt.max_iter}")
    messenger.log_only(f"  Scoring metric: {opt.scoring_metric}")
    messenger.log_only(f"  Penalty options: {list(opt.penalty_options)}")
    messenger.log_only(f"  Loss options: {list(opt.loss_options)}")
    messenger.log_only(
        f"  Class weight multipliers: {list(opt.class_weight_multipliers)}"
    )
    if opt.features:
        messenger.log_only(f"  Features: {list(opt.features)}")
    else:
        messenger.log_only("  Features: default (4D)")
    if opt.gamma_imbalance_options:
        messenger.log_only(
            f"  Gamma imbalance options: {list(opt.gamma_imbalance_options)}"
        )
    if opt.param_grid_override:
        messenger.log_only(
            f"  Parameter grid: {len(opt.param_grid_override)} hyperparameter sets"
        )
    messenger.log_only(
        f"  C bounds: eff_C_pos_range={opt.eff_C_pos_range}, eff_C_neg_max={opt.eff_C_neg_max}"
    )

    # Log ensemble settings
    messenger.log_only("Ensemble configuration:")
    messenger.log_only(f"  n_models: {ens.n_models}")
    messenger.log_only(f"  subsample_u2: {ens.subsample_u2}")
    messenger.log_only(f"  subsample_ratio: {ens.subsample_ratio}")
    messenger.log_only(f"  max_iter: {ens.max_iter}")
    messenger.log_only(f"  random_state: {ens.random_state}")

    # Create classifier with correct parameter names
    # IntronClassifier API uses:
    # - classification_threshold (not threshold)
    # - n_ensemble_models (not n_models)
    # - fixed_c (not fixed_C)
    # - random_state (not seed)
    # - cv_processes (for cross-validation parallelization)
    # - classification_processes (for prediction parallelization)
    # - param_grid_override (optional custom parameter grid)
    # - n_points_initial (initial grid points for optimization)
    # - subsample_u2, subsample_ratio (ensemble diversity settings)
    # - max_iter (for LinearSVC convergence during ensemble training)

    # Get use_fold_averaged_params: CLI arg takes precedence over YAML config
    # Check yaml_config for training.use_fold_averaged_params
    yaml_use_fold_averaged = config.yaml_config.get("training", {}).get(
        "use_fold_averaged_params"
    )
    if config.training.use_fold_averaged_params is not None:
        # CLI arg was explicitly provided
        use_fold_averaged_params = config.training.use_fold_averaged_params
    elif yaml_use_fold_averaged is not None:
        # Use YAML value as fallback
        use_fold_averaged_params = yaml_use_fold_averaged
    else:
        # Default
        use_fold_averaged_params = False

    classifier = IntronClassifier(
        n_optimization_rounds=opt.n_rounds,
        classification_threshold=config.scoring.threshold,
        n_ensemble_models=ens.n_models,
        subsample_u2=ens.subsample_u2,
        subsample_ratio=ens.subsample_ratio,
        fixed_c=config.training.fixed_C,
        optimize_c=(config.training.fixed_C is None),
        random_state=ens.random_state,
        cv_processes=opt.n_jobs,
        classification_processes=config.performance.processes,
        max_iter=ens.max_iter,
        eval_mode=config.training.eval_mode,
        n_cv_folds=opt.cv_folds,
        test_fraction=config.training.test_fraction,
        scoring_metric=opt.scoring_metric,
        penalty_options=list(opt.penalty_options),
        loss_options=list(opt.loss_options),
        class_weight_multipliers=list(opt.class_weight_multipliers),
        use_multiplier_tiebreaker=opt.use_multiplier_tiebreaker,
        features_list=list(opt.features) if opt.features else None,
        gamma_imbalance_options=list(opt.gamma_imbalance_options)
        if opt.gamma_imbalance_options
        else None,
        param_grid_override=opt.param_grid_override,
        n_points_initial=opt.n_points_initial,
        eff_C_pos_range=opt.eff_C_pos_range,
        eff_C_neg_max=opt.eff_C_neg_max,
        use_fold_averaged_params=use_fold_averaged_params,
    )

    # Run complete classification pipeline (optimize + train + classify)
    messenger.log_only(f"Running classification on {len(introns)} experimental introns")
    messenger.log_only(
        f"Reference data: {len(u12_reference)} U12-type, {len(u2_reference)} U2-type"
    )

    result = classifier.classify(
        u12_reference=u12_reference, u2_reference=u2_reference, experimental=introns
    )

    # Extract metrics from result
    metrics = {
        "optimized_C": result.parameters.C,
        "cv_score": result.parameters.cv_score,
        "n_models": len(result.ensemble.models),
    }

    # Add evaluation metrics if nested CV or split eval was performed
    if result.eval_result is not None:
        # Check if it's NestedCVResult (has mean_f1) or SplitEvalResult (has f1_score)
        if hasattr(result.eval_result, "mean_f1"):
            # Nested CV result
            metrics["mean_f1"] = result.eval_result.mean_f1
            metrics["std_f1"] = result.eval_result.std_f1
            metrics["mean_pr_auc"] = result.eval_result.mean_pr_auc
            metrics["std_pr_auc"] = result.eval_result.std_pr_auc
            metrics["n_cv_folds"] = result.eval_result.n_folds
        elif hasattr(result.eval_result, "f1_score"):
            # Split evaluation result
            metrics["f1"] = result.eval_result.f1_score
            metrics["pr_auc"] = result.eval_result.pr_auc

    messenger.log_only("Classification complete")
    messenger.log_only(f"  Optimized C: {metrics['optimized_C']:.6e}")
    messenger.log_only(f"  Models trained: {metrics['n_models']}")

    # Log evaluation metrics if available
    if "mean_f1" in metrics:
        messenger.log_only(f"  Nested CV results ({metrics['n_cv_folds']} folds):")
        messenger.log_only(
            f"    Mean F1: {metrics['mean_f1']:.4f} ± {metrics['std_f1']:.4f}"
        )
        messenger.log_only(
            f"    Mean PR-AUC: {metrics['mean_pr_auc']:.4f} ± {metrics['std_pr_auc']:.4f}"
        )

        # Add per-fold results table to log
        eval_result = result.eval_result
        if eval_result is not None and hasattr(eval_result, "fold_results"):
            messenger.print_nested_cv_results(
                n_folds=eval_result.n_folds,
                mean_f1=eval_result.mean_f1,
                std_f1=eval_result.std_f1,
                mean_pr_auc=eval_result.mean_pr_auc,
                std_pr_auc=eval_result.std_pr_auc,
                fold_results=eval_result.fold_results,
            )
    elif "f1" in metrics:
        messenger.log_only("  Test set evaluation:")
        messenger.log_only(f"    F1: {metrics['f1']:.4f}")
        messenger.log_only(f"    PR-AUC: {metrics['pr_auc']:.4f}")

    # Generate training reference plots if evaluation was performed
    if result.eval_result is not None:
        messenger.log_only("Generating training reference plots")
        try:
            from intronIC.visualization.plots import plot_training_results

            # Extract normalized z-scores from reference introns
            u2_scores = np.array(
                [
                    [i.scores.five_z_score, i.scores.bp_z_score]
                    for i in u2_reference
                    if i.scores
                    and i.scores.five_z_score is not None
                    and i.scores.bp_z_score is not None
                ]
            )
            u12_scores = np.array(
                [
                    [i.scores.five_z_score, i.scores.bp_z_score]
                    for i in u12_reference
                    if i.scores
                    and i.scores.five_z_score is not None
                    and i.scores.bp_z_score is not None
                ]
            )

            # Get PR curves and AUC based on evaluation type
            if hasattr(result.eval_result, "mean_f1"):
                # Nested CV - multiple curves
                pr_curves = result.eval_result.pr_curves
                pr_auc = result.eval_result.mean_pr_auc
            else:
                # Split eval - single curve
                pr_curves = [(result.eval_result.precision, result.eval_result.recall)]
                pr_auc = result.eval_result.test_pr_auc

            plot_training_results(
                u2_scores=u2_scores,
                u12_scores=u12_scores,
                pr_curves=pr_curves,
                pr_auc=pr_auc,
                output_dir=config.output.output_dir,
                species_name=config.output.base_filename,
                fig_dpi=300,
            )
            messenger.log_only("Successfully generated training reference plots")
        except Exception as plot_error:
            import traceback

            messenger.warning(f"Failed to generate training plots: {plot_error}")
            messenger.warning(f"Traceback: {traceback.format_exc()}")
            # Continue even if plotting fails

    # Log learned coefficients from ensemble
    log_ensemble_coefficients(result.ensemble, messenger)

    # Save trained model with human-trained normalizer for cross-species classification
    model_path = config.output.get_output_path(".model.pkl")
    messenger.log_only(f"Saving trained model to {model_path}")

    # Compute human U2 margin statistics for margin alignment in adaptive mode
    # Use first model in ensemble as representative
    messenger.log_only("Computing U2 margin statistics for cross-species adaptation")
    first_model_obj = result.ensemble.models[0]
    first_model = first_model_obj.model  # This is the sklearn Pipeline

    # Extract features from U2 reference introns
    u2_features = np.array(
        [
            [i.scores.five_z_score, i.scores.bp_z_score, i.scores.three_z_score]
            for i in u2_reference
            if i.scores is not None
        ]
    )

    # Get the transformer and base estimator from the Pipeline
    # Pipeline structure: scale -> transform -> svc -> calibration
    # We need to transform features then get decision_function from base SVM
    try:
        # Transform features through the pipeline up to (but not including) calibration
        # The pipeline has: scale, transform, svc, and CalibratedClassifierCV wraps the whole thing
        if hasattr(first_model, "calibrated_classifiers_"):
            # Model is CalibratedClassifierCV - get the base pipeline
            base_pipeline = first_model.calibrated_classifiers_[0].estimator
        else:
            # Direct pipeline
            base_pipeline = first_model

        # Transform features and get decision function
        # Note: decision_function is before Platt calibration
        u2_margins = base_pipeline.decision_function(u2_features)

        # Compute robust statistics (median and IQR)
        mu_u2 = float(np.median(u2_margins))
        q25, q75 = np.percentile(u2_margins, [25, 75])
        sigma_u2 = float(q75 - q25)  # IQR

        messenger.log_only(
            f"  U2 margin stats: median={mu_u2:.3f}, IQR={sigma_u2:.3f}, N={len(u2_margins):,}"
        )

        human_negative_stats = {
            "margin_median": mu_u2,
            "margin_iqr": sigma_u2,
            "n_samples": len(u2_margins),
        }
    except Exception as e:
        messenger.warning(f"Failed to compute margin statistics: {e}")
        messenger.warning("Model will not support margin-aligned adaptive mode")
        human_negative_stats = None

    # Compute training prior for prior adjustment
    n_u12 = len(u12_reference)
    n_u2 = len(u2_reference)
    training_prior = n_u12 / (n_u12 + n_u2)
    messenger.log_only(
        f"  Training U12 prior: {training_prior:.4f} ({n_u12}/{n_u12 + n_u2})"
    )

    # Build model bundle
    model_bundle = {
        "ensemble": result.ensemble,
        "normalizer": normalizer,  # Save human-trained scaler for cross-species use
        "threshold": config.scoring.threshold,
        "human_negative_stats": human_negative_stats,  # For margin alignment (NEW)
        "training_prior": training_prior,  # For prior adjustment (NEW)
    }
    joblib.dump(model_bundle, model_path, compress=3)
    messenger.log_only(
        "Model saved successfully with cross-species adaptation statistics"
    )

    # Generate and save training metadata
    metadata_path = model_path.with_suffix(".metadata.json")
    messenger.log_only("Generating training metadata")
    metadata = generate_training_metadata(
        model_name=model_path.stem,
        u12_reference_path=config.scoring.reference_u12s,
        u2_reference_path=config.scoring.reference_u2s,
        u12_introns=u12_reference,
        u2_introns=u2_reference,
        optimized_C=result.parameters.C,
        calibration_method=result.parameters.calibration_method,
        cv_score=result.parameters.cv_score,
        n_models=len(result.ensemble.models),
        threshold=config.scoring.threshold,
        eval_result=result.eval_result,
        max_iter=config.training.max_iter,
        kernel="linear",
        seed=config.training.seed,
    )
    write_metadata(metadata, metadata_path)
    messenger.log_only(f"Training metadata saved to {metadata_path}")

    return list(result.classified_introns), metrics


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

    # Filter duplicates if not including them
    # Port from: intronIC.py:4806-4807
    if not config.extraction.include_duplicates:
        original_count = len(introns)
        introns = [i for i in introns if not (i.metadata and i.metadata.duplicate)]
        filtered_count = original_count - len(introns)
        if filtered_count > 0:
            messenger.log_only(
                f"Filtered out {filtered_count} duplicate introns (use -d to include)"
            )

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
            data_dir / "u12_reference.introns.iic.gz"
        )
        u2_file = config.scoring.reference_u2s or (
            data_dir / "u2_reference.introns.iic.gz"
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

        # Normalize scores
        messenger.log_only("Normalizing scores with z-score transformation")
        from intronIC.scoring.normalizer import ScoreNormalizer

        normalizer = ScoreNormalizer()
        normalizer.fit(scored_reference, dataset_type="reference")

        u12_ref_norm = list(normalizer.transform(u12_scored, dataset_type="reference"))
        u2_ref_norm = list(normalizer.transform(u2_scored, dataset_type="reference"))

        messenger.success(
            f"Loaded, scored, and normalized {len(u12_ref_norm)} U12 and {len(u2_ref_norm)} U2 reference introns"
        )

        # Step 4: Train classifier (model is saved internally by classify_introns)
        messenger.step(4, "Train Classifier", pipeline_steps)
        # Pass empty list for experimental introns - we're only training on references
        # classify_introns() will train the model and save it to disk
        classified_introns, metrics = classify_introns(
            introns=[],  # No experimental introns in train mode
            u12_reference=u12_ref_norm,
            u2_reference=u2_ref_norm,
            normalizer=normalizer,
            config=config,
            messenger=messenger,
            reporter=reporter,
        )
        messenger.success("Model trained and saved")

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
            # Standard extraction
            introns = extract_introns_from_annotation(config, messenger, reporter)
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
        # The function handles extraction, scoring, classification, and output internally.
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

        if config.input.mode == "annotation":
            # Don't load genome here - extraction handles it internally
            # (parallel mode uses indexed access, sequential uses cache)
            if config.performance.streaming:
                # Streaming mode: ~85% memory savings, stores sequences in SQLite
                introns, streaming_db_path = extract_introns_streaming(
                    config, messenger, reporter
                )
            else:
                # Standard processing: faster but uses more memory
                introns = extract_introns_from_annotation(config, messenger, reporter)
        elif config.input.mode == "bed":
            genome_reader = load_genome(config, messenger)
            introns = extract_introns_from_bed(
                config, genome_reader, messenger, reporter
            )
        elif config.input.mode == "sequences":
            introns = load_introns_from_sequences(config, messenger)
        else:
            raise ValueError(f"Unknown input mode: {config.input.mode}")

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
        # for potential output (user may want duplicates via -d flag)
        introns_for_scoring = filtered_introns

        # Step 2: Score introns
        messenger.step(2, "Score Introns with PWMs", pipeline_steps)

        scored_introns = score_introns(introns_for_scoring, config, messenger, reporter)

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
            classified_introns, metrics = classify_with_pretrained_model(
                scored_introns,
                config.training.pretrained_model_path,
                config,
                messenger,
                reporter,
            )
        else:
            # Normal flow: normalize + train + classify
            # Step 3: Normalize scores
            messenger.step(3, "Normalize Scores", pipeline_steps)
            normalized_introns, u12_reference, u2_reference, normalizer = (
                normalize_scores(scored_introns, config, messenger, reporter)
            )
            messenger.success("Scores normalized")

            # Step 4: Train and apply classifier
            messenger.step(4, "Train and Apply Classifier", pipeline_steps)
            classified_introns, metrics = classify_introns(
                normalized_introns,
                u12_reference,
                u2_reference,
                normalizer,
                config,
                messenger,
                reporter,
            )

        # Count classifications based on threshold (for reporting "high confidence" U12s)
        # Note: type_id is based on raw classifier (>50%), but reporting uses threshold
        threshold = config.scoring.threshold
        u12_count = sum(
            1
            for i in classified_introns
            if i.scores is not None
            and i.scores.svm_score is not None
            and i.scores.svm_score >= threshold
        )
        u2_count = len(classified_introns) - u12_count

        # Count AT-AC introns among high-confidence U12s (score >= threshold)
        atac_count = sum(
            1
            for i in classified_introns
            if (
                i.scores is not None
                and i.scores.svm_score is not None
                and i.scores.svm_score >= config.scoring.threshold
                and i.sequences
                and i.sequences.terminal_dinucleotides == "AT-AC"
            )
        )

        # Display classification summary with unified formatting
        total_classified = len(classified_introns)
        messenger.print_classification_results(
            total=total_classified,
            u12_count=u12_count,
            u2_count=u2_count,
            atac_count=atac_count,
            threshold=config.scoring.threshold,
        )

        # Collect and log splice site boundary statistics (separate by U12/U2)
        from collections import Counter

        boundaries_u12 = Counter()
        boundaries_u2 = Counter()

        for intron in classified_introns:
            # Count all boundaries (canonical and non-canonical) by type
            # Use threshold to match the main counts, not raw classifier type_id
            if (
                intron.metadata
                and intron.sequences
                and intron.sequences.terminal_dinucleotides
            ):
                dnts = intron.sequences.terminal_dinucleotides
                if (
                    intron.scores is not None
                    and intron.scores.svm_score is not None
                    and intron.scores.svm_score >= config.scoring.threshold
                ):
                    boundaries_u12[dnts] += 1
                else:
                    boundaries_u2[dnts] += 1

        # Display U12 boundary statistics with unified formatting
        if boundaries_u12:
            # Sort by count (descending), then alphabetically by dinucleotide
            sorted_boundaries = sorted(
                boundaries_u12.items(), key=lambda x: (-x[1], x[0])
            )
            messenger.print_dinucleotide_boundaries(
                intron_type="U12-type", boundaries=sorted_boundaries, top_n=20
            )

        # Display U2 boundary statistics with unified formatting
        if boundaries_u2:
            # Sort by count (descending), then alphabetically by dinucleotide
            sorted_boundaries = sorted(
                boundaries_u2.items(), key=lambda x: (-x[1], x[0])
            )
            messenger.print_dinucleotide_boundaries(
                intron_type="U2-type", boundaries=sorted_boundaries, top_n=20
            )

        # Save classification metrics to JSON file
        if metrics:
            metrics_path = config.output.get_output_path(".metrics.iic.json")
            messenger.log_only(f"Saving classification metrics to {metrics_path}")
            with open(metrics_path, "w") as f:
                json.dump(metrics, f, indent=2)

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

    # Find bundled test data
    data_dir = Path(__file__).parent.parent / "data" / "test_data"
    genome_file = data_dir / "Homo_sapiens.Chr19.Ensembl_91.fa.gz"
    annotation_file = data_dir / "Homo_sapiens.Chr19.Ensembl_91.gff3.gz"

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
        console.print(f"                    -n homo_sapiens_chr19 -p 4")
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

    # Build command-line args for classify mode
    test_args = [
        "classify",
        "-g", str(genome_file),
        "-a", str(annotation_file),
        "-n", "test_chr19",
        "-o", str(output_dir),
        "-p", str(args.processes),
    ]

    # Run classification
    start_time = time.time()
    try:
        main(test_args)
        elapsed = time.time() - start_time

        # Check results
        meta_file = output_dir / "test_chr19.meta.iic"
        if meta_file.exists():
            # Parse metadata to get counts
            import json
            total_introns = "?"
            u12_introns = "?"
            with open(meta_file) as f:
                for line in f:
                    if line.startswith("#"):
                        try:
                            meta = json.loads(line[1:].strip())
                            total_introns = meta.get("total_introns", "?")
                            u12_introns = meta.get("u12_introns", "?")
                            break
                        except:
                            pass

            console.print(f"\n[bold green]✓ Test completed successfully![/bold green]")
            console.print(f"  Runtime: {elapsed:.1f}s")
            console.print(f"  Total introns: {total_introns}")
            console.print(f"  U12 introns: {u12_introns}")

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
