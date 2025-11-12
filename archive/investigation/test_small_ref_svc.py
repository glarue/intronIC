#!/usr/bin/env python3
"""
Test SVC with CalibratedClassifierCV using a small reference subset.

This uses only 50 U12s and 5000 U2s to speed up testing while still
maintaining the same pipeline structure.
"""

import sys
import time
import random
from pathlib import Path
from argparse import Namespace

# Add intronIC to path
sys.path.insert(0, str(Path(__file__).parent / "intronIC"))

from intronIC import (
    get_reference_introns,
    get_raw_scores,
    scale_scores,
    optimize_svm,
    train_svm,
    parallel_svm_score,
    average_svm_score_info
)

def log(msg):
    """Log with timestamp."""
    timestamp = time.strftime("%H:%M:%S")
    print(f"[{timestamp}] {msg}", flush=True)

def main():
    log("="*80)
    log("Testing SVC + CalibratedClassifierCV with SMALL reference subset")
    log("="*80)

    # Setup args
    args = Namespace()
    args.u12_type = 'gt-ag'
    args.reference_u12s = "intronIC/data/u12_reference.introns.iic.gz"
    args.reference_u2s = "intronIC/data/u2_reference.introns.iic.gz"
    args.pwms = "intronIC/data/scoring_matrices.fasta.iic"
    args.five_score_coords = (-3, 9)
    args.three_score_coords = (-6, 4)
    args.bp_search_coords = (-55, -5)
    args.no_nc = False
    args.feature_set = 'include_both'
    args.threshold = 90.0
    args.C = None
    args.processes = 4
    args.cv_processes = 4
    args.verbose = True

    # Extract experimental introns from test data
    log("Extracting experimental introns from test data...")
    args.annotation = "intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz"
    args.genome = "intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz"
    args.feature = 'both'
    args.min_intron_len = 30
    args.no_intron_overlap = False
    args.allow_multiple_isoforms = True

    from intronIC import introns_from_annotation
    experimental_introns = list(introns_from_annotation(args))
    log(f"Extracted {len(experimental_introns)} experimental introns")

    # Load full reference set
    log("\nLoading FULL reference sequences...")
    ref_introns_full = get_reference_introns(args)
    u12_full = [i for i in ref_introns_full if i.type_id == 'u12']
    u2_full = [i for i in ref_introns_full if i.type_id == 'u2']
    log(f"Full reference: {len(u12_full)} U12, {len(u2_full)} U2")

    # Sample smaller subset
    log("\nSampling SMALL reference subset...")
    random.seed(42)
    u12_small = random.sample(u12_full, min(50, len(u12_full)))
    u2_small = random.sample(u2_full, min(5000, len(u2_full)))
    ref_introns_small = u12_small + u2_small
    log(f"Small reference: {len(u12_small)} U12, {len(u2_small)} U2")

    # Score all introns
    log("\nScoring experimental introns with PWMs...")
    start = time.time()
    get_raw_scores(
        experimental_introns,
        args.pwms,
        args.five_score_coords,
        args.three_score_coords,
        args.bp_search_coords,
        args.u12_type,
        build_u2_bp=False,
        bp_region=None
    )
    log(f"Scored experimental introns in {time.time()-start:.1f}s")

    log("\nScoring reference introns with PWMs...")
    start = time.time()
    get_raw_scores(
        ref_introns_small,
        args.pwms,
        args.five_score_coords,
        args.three_score_coords,
        args.bp_search_coords,
        args.u12_type,
        build_u2_bp=False,
        bp_region=None
    )
    log(f"Scored reference introns in {time.time()-start:.1f}s")

    # Normalize scores
    log("\nNormalizing scores...")
    start = time.time()
    scale_scores(experimental_introns, ref_introns_small)
    log(f"Normalized scores in {time.time()-start:.1f}s")

    # Classification
    log("\n" + "="*80)
    log("CLASSIFICATION PIPELINE")
    log("="*80)

    # Stage 1: Optimize hyperparameters
    log("\n[Stage 1/4] Hyperparameter optimization...")
    log("This will take several minutes with the small reference set")
    start = time.time()

    try:
        best_params = optimize_svm(
            ref_introns_small,
            args,
            use_probability=False  # Test with CalibratedClassifierCV
        )
        opt_time = time.time() - start
        log(f"✓ Optimization complete in {opt_time:.1f}s")
        log(f"  Best C: {best_params['C']:.2e}")
        if 'method' in best_params:
            log(f"  Best calibration method: {best_params['method']}")
    except Exception as e:
        log(f"✗ Optimization FAILED: {e}")
        import traceback
        traceback.print_exc()
        return 1

    # Stage 2: Train models
    log("\n[Stage 2/4] Training SVM models...")
    start = time.time()

    try:
        models = train_svm(
            ref_introns_small,
            args,
            best_params=best_params,
            num_models=3,  # Use 3 models for ensemble
            subsample_u2=False,  # Don't subsample the already-small set
            use_probability=False  # Use CalibratedClassifierCV
        )
        train_time = time.time() - start
        log(f"✓ Trained {len(models)} models in {train_time:.1f}s")
        for i, (model, metrics) in enumerate(models):
            log(f"  Model {i+1}: F1={metrics['f1_score']:.4f}")
    except Exception as e:
        log(f"✗ Training FAILED: {e}")
        import traceback
        traceback.print_exc()
        return 1

    # Stage 3: Score experimental introns
    log("\n[Stage 3/4] Scoring experimental introns...")
    start = time.time()

    try:
        parallel_svm_score(
            experimental_introns,
            models,
            processes=args.processes,
            method='predict_proba'
        )
        score_time = time.time() - start
        log(f"✓ Scored {len(experimental_introns)} introns in {score_time:.1f}s")
    except Exception as e:
        log(f"✗ Scoring FAILED: {e}")
        import traceback
        traceback.print_exc()
        return 1

    # Stage 4: Aggregate scores
    log("\n[Stage 4/4] Aggregating ensemble scores...")
    start = time.time()

    try:
        average_svm_score_info(
            experimental_introns,
            models,
            threshold=args.threshold
        )
        agg_time = time.time() - start
        log(f"✓ Aggregated scores in {agg_time:.1f}s")
    except Exception as e:
        log(f"✗ Aggregation FAILED: {e}")
        import traceback
        traceback.print_exc()
        return 1

    # Results
    log("\n" + "="*80)
    log("RESULTS")
    log("="*80)

    u12_predicted = sum(1 for i in experimental_introns if i.type_id == 'u12')
    u2_predicted = sum(1 for i in experimental_introns if i.type_id == 'u2')

    log(f"Total experimental introns: {len(experimental_introns)}")
    log(f"  Predicted U12: {u12_predicted} ({100*u12_predicted/len(experimental_introns):.2f}%)")
    log(f"  Predicted U2:  {u2_predicted} ({100*u2_predicted/len(experimental_introns):.2f}%)")

    # Timing summary
    log("\n" + "="*80)
    log("TIMING SUMMARY")
    log("="*80)
    log(f"Optimization: {opt_time:6.1f}s")
    log(f"Training:     {train_time:6.1f}s")
    log(f"Scoring:      {score_time:6.1f}s")
    log(f"Aggregation:  {agg_time:6.1f}s")
    log(f"Total:        {opt_time+train_time+score_time+agg_time:6.1f}s")

    # Get top U12 predictions
    log("\n" + "="*80)
    log("TOP 10 U12 PREDICTIONS")
    log("="*80)

    scored = [i for i in experimental_introns if hasattr(i, 'svm_score')]
    scored.sort(key=lambda x: x.svm_score, reverse=True)

    for i, intron in enumerate(scored[:10], 1):
        log(f"{i:2d}. {intron.label:40s} {intron.svm_score:5.1f}%  {intron.type_id}")

    log("\n" + "="*80)
    log("TEST COMPLETE!")
    log("="*80)

    return 0

if __name__ == "__main__":
    sys.exit(main())
