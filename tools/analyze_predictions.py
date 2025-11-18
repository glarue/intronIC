#!/usr/bin/env python3
"""
Analyze prediction patterns from score_info.iic files.

This tool performs post-hoc analysis of model predictions by:
1. Computing all BothEndsStrong features manually from z-scores
2. Identifying patterns (one-end-strong, negative scores, imbalance)
3. Helping diagnose false positives and model behavior

Usage:
    # Analyze top predicted U12s
    python tools/analyze_predictions.py run_tests/species.score_info.iic

    # Show only top 5
    python tools/analyze_predictions.py run_tests/species.score_info.iic -n 5

    # Show both high-scoring and low-scoring predictions
    python tools/analyze_predictions.py run_tests/species.score_info.iic --show-extremes

    # Filter by score range
    python tools/analyze_predictions.py run_tests/species.score_info.iic --min-score 5.0

See: docs/POST_HOC_ANALYSIS.md for methodology and interpretation
"""
import sys
import argparse
import numpy as np
from pathlib import Path


def compute_features(s5, sBP, s3):
    """
    Compute all BothEndsStrong features from base z-scores.

    This replicates the feature computation done by the transformer,
    allowing post-hoc analysis of what features contributed to a prediction.

    Args:
        s5: 5' splice site z-score
        sBP: Branch point z-score
        s3: 3' splice site z-score

    Returns:
        Dict of feature_name: value for all possible features
    """
    # Base features
    features = {
        's5': s5,
        'sBP': sBP,
        's3': s3
    }

    # Pairwise min features (may or may not be in model)
    features['min_5_bp'] = min(s5, sBP)
    features['min_5_3'] = min(s5, s3)
    features['min_bp_3'] = min(sBP, s3)

    # 3-way min (recommended feature)
    features['min_all'] = min(s5, sBP, s3)

    # Pairwise max features (may or may not be in model)
    features['max_5_bp'] = max(s5, sBP)
    features['max_5_3'] = max(s5, s3)
    features['max_bp_3'] = max(sBP, s3)
    features['max_all'] = max(s5, sBP, s3)

    # Imbalance penalties (should be in model)
    features['neg_absdiff_5_bp'] = -abs(s5 - sBP)
    features['neg_absdiff_5_3'] = -abs(s5 - s3)
    features['neg_absdiff_bp_3'] = -abs(sBP - s3)

    return features


def identify_pattern(s5, sBP, s3, features):
    """
    Identify problematic patterns in the z-scores.

    Returns:
        List of warning strings describing issues
    """
    warnings = []

    # Check for negative scores
    negative_scores = []
    if s5 < 0:
        negative_scores.append('s5')
    if sBP < 0:
        negative_scores.append('sBP')
    if s3 < 0:
        negative_scores.append('s3')

    if negative_scores:
        warnings.append(f"NEGATIVE scores: {', '.join(negative_scores)}")

    # Check for one-end-strong pattern
    # Defined as: exactly one score > 1.0, others < 1.0
    high_scores = sum(1 for score in [s5, sBP, s3] if score > 1.0)
    if high_scores == 1:
        warnings.append("ONE-END-STRONG pattern detected")

    # Check for high imbalance
    if abs(features['neg_absdiff_5_bp']) > 1.0:
        warnings.append(f"High 5'/BP imbalance (|diff|={abs(features['neg_absdiff_5_bp']):.2f})")
    if abs(features['neg_absdiff_5_3']) > 1.0:
        warnings.append(f"High 5'/3' imbalance (|diff|={abs(features['neg_absdiff_5_3']):.2f})")
    if abs(features['neg_absdiff_bp_3']) > 1.0:
        warnings.append(f"High BP/3' imbalance (|diff|={abs(features['neg_absdiff_bp_3']):.2f})")

    # Check if min_all is negative (should strongly reject)
    if features['min_all'] < 0:
        warnings.append("min_all is NEGATIVE (should strongly penalize)")

    # Check if all scores are weak but balanced
    if all(abs(score) < 0.5 for score in [s5, sBP, s3]):
        warnings.append("All scores weak (<0.5)")

    return warnings


def parse_score_file(score_file):
    """
    Parse score_info.iic file.

    Format (tab-delimited):
        name  rel_score  svm_score  5'_seq  5'_raw  5'_z
        bp_seq  bp_seq_u2  bp_raw  bp_z  3'_seq  3'_raw  3'_z
        min(5,bp)  min(5,3)  max(5,bp)  max(5,3)  decision_dist

    Returns:
        List of dicts with parsed data
    """
    introns = []
    with open(score_file) as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 13:
                continue

            try:
                # Parse fields
                name = fields[0]
                rel_score = float(fields[1])
                svm_score = float(fields[2])
                five_z = float(fields[5])
                bp_z = float(fields[9])
                three_z = float(fields[12])

                introns.append({
                    'name': name,
                    'rel_score': rel_score,
                    'svm_score': svm_score,
                    's5': five_z,
                    'sBP': bp_z,
                    's3': three_z
                })
            except (ValueError, IndexError):
                continue

    return introns


def print_intron_analysis(intron_data, index, total):
    """Print detailed analysis for one intron."""
    print(f"\n{'='*80}")
    print(f"#{index}/{total}: {intron_data['name'][:70]}")
    print(f"  Relative score: {intron_data['rel_score']:.2f}")
    print(f"  SVM probability: {intron_data['svm_score']:.1f}%")

    # Compute all features
    features = compute_features(
        intron_data['s5'],
        intron_data['sBP'],
        intron_data['s3']
    )

    # Print base features
    print(f"\n  Base z-scores:")
    print(f"    s5:  {features['s5']:>7.4f}")
    print(f"    sBP: {features['sBP']:>7.4f}")
    print(f"    s3:  {features['s3']:>7.4f}")

    # Print min features (model may use these)
    print(f"\n  Min features (require all/both strong):")
    print(f"    min_all:   {features['min_all']:>7.4f}  ← ALL THREE must be strong")
    print(f"    min_5_bp:  {features['min_5_bp']:>7.4f}  ← Both 5' and BP")
    print(f"    min_5_3:   {features['min_5_3']:>7.4f}  ← Both 5' and 3'")
    print(f"    min_bp_3:  {features['min_bp_3']:>7.4f}  ← Both BP and 3'")

    # Print max features (model may use these - usually bad sign)
    print(f"\n  Max features (at least one strong - can reward FPs):")
    print(f"    max_all:   {features['max_all']:>7.4f}")
    print(f"    max_5_bp:  {features['max_5_bp']:>7.4f}")
    print(f"    max_5_3:   {features['max_5_3']:>7.4f}")
    print(f"    max_bp_3:  {features['max_bp_3']:>7.4f}")

    # Print penalties (should be in model)
    print(f"\n  Imbalance penalties (should be small for real U12s):")
    print(f"    neg_absdiff_5_bp: {features['neg_absdiff_5_bp']:>7.4f}  ← 5'/BP balance")
    print(f"    neg_absdiff_5_3:  {features['neg_absdiff_5_3']:>7.4f}  ← 5'/3' balance")
    print(f"    neg_absdiff_bp_3: {features['neg_absdiff_bp_3']:>7.4f}  ← BP/3' balance")

    # Identify patterns
    warnings = identify_pattern(
        intron_data['s5'],
        intron_data['sBP'],
        intron_data['s3'],
        features
    )

    if warnings:
        print(f"\n  ⚠ Pattern warnings:")
        for warning in warnings:
            print(f"    • {warning}")
    else:
        print(f"\n  ✓ No obvious problems detected")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze predictions from score_info.iic file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        'score_file',
        type=Path,
        help='Path to .score_info.iic file'
    )
    parser.add_argument(
        '-n', '--top-n',
        type=int,
        default=10,
        help='Number of top predictions to analyze (default: 10)'
    )
    parser.add_argument(
        '--min-score',
        type=float,
        default=0.0,
        help='Minimum relative score to include (default: 0.0, i.e., predicted U12s)'
    )
    parser.add_argument(
        '--max-score',
        type=float,
        default=None,
        help='Maximum relative score to include (default: none)'
    )
    parser.add_argument(
        '--show-extremes',
        action='store_true',
        help='Show both highest and lowest scoring predictions'
    )
    parser.add_argument(
        '--summary',
        action='store_true',
        help='Show summary table only, without detailed analysis'
    )

    args = parser.parse_args()

    if not args.score_file.exists():
        print(f"Error: {args.score_file} not found", file=sys.stderr)
        sys.exit(1)

    # Parse file
    print(f"Reading {args.score_file}...")
    introns = parse_score_file(args.score_file)

    # Filter by score range
    filtered = [
        i for i in introns
        if i['rel_score'] >= args.min_score and
           (args.max_score is None or i['rel_score'] <= args.max_score)
    ]

    # Sort by relative score
    filtered.sort(key=lambda x: x['rel_score'], reverse=True)

    print(f"\nFound {len(introns)} total introns")
    print(f"After filtering: {len(filtered)} introns")

    if not filtered:
        print("No introns match criteria")
        return

    # Show summary table
    print(f"\n{'='*80}")
    print("SUMMARY TABLE")
    print(f"{'='*80}")
    print(f"{'#':<4} {'RelScore':>9} {'SVM%':>6} {'s5':>7} {'sBP':>7} {'s3':>7} {'min_all':>8} {'Pattern':<20}")
    print("-" * 80)

    n_show = min(args.top_n, len(filtered))
    for i, intron in enumerate(filtered[:n_show], 1):
        features = compute_features(intron['s5'], intron['sBP'], intron['s3'])
        warnings = identify_pattern(intron['s5'], intron['sBP'], intron['s3'], features)
        pattern = warnings[0][:20] if warnings else "OK"

        print(f"{i:<4} {intron['rel_score']:>9.2f} {intron['svm_score']:>6.1f} "
              f"{intron['s5']:>7.2f} {intron['sBP']:>7.2f} {intron['s3']:>7.2f} "
              f"{features['min_all']:>8.2f} {pattern:<20}")

    # Show detailed analysis unless --summary
    if not args.summary:
        print(f"\n{'='*80}")
        print(f"DETAILED ANALYSIS (Top {n_show})")
        print(f"{'='*80}")

        for i, intron in enumerate(filtered[:n_show], 1):
            print_intron_analysis(intron, i, n_show)

    # Show extremes if requested
    if args.show_extremes and len(filtered) > args.top_n:
        print(f"\n{'='*80}")
        print(f"LOWEST SCORING PREDICTIONS")
        print(f"{'='*80}")

        for i, intron in enumerate(filtered[-args.top_n:], 1):
            print_intron_analysis(intron, i, args.top_n)

    print(f"\n{'='*80}")
    print("Analysis complete")
    print(f"{'='*80}")


if __name__ == '__main__':
    main()
