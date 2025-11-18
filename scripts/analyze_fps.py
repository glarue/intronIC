#!/usr/bin/env python3
"""
Analyze top-scoring false positives to understand feature contributions.
"""
import sys
import numpy as np
from pathlib import Path

def compute_features(s5, sBP, s3):
    """
    Compute all BothEndsStrong features from base z-scores.

    Returns dict of feature_name: value
    """
    # Base features
    features = {
        's5': s5,
        'sBP': sBP,
        's3': s3
    }

    # Pairwise min features
    features['min_5_bp'] = min(s5, sBP)
    features['min_5_3'] = min(s5, s3)
    features['min_bp_3'] = min(sBP, s3)

    # 3-way min
    features['min_all'] = min(s5, sBP, s3)

    # Pairwise max features
    features['max_5_bp'] = max(s5, sBP)
    features['max_5_3'] = max(s5, s3)
    features['max_bp_3'] = max(sBP, s3)
    features['max_all'] = max(s5, sBP, s3)

    # Imbalance penalties
    features['neg_absdiff_5_bp'] = -abs(s5 - sBP)
    features['neg_absdiff_5_3'] = -abs(s5 - s3)
    features['neg_absdiff_bp_3'] = -abs(sBP - s3)

    return features

def analyze_score_file(score_file, n_top=10):
    """
    Analyze top-scoring introns from score_info file.
    """
    print(f"Reading {score_file}...")

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

                # Only look at predicted U12s (positive scores)
                if rel_score > 0:
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

    # Sort by relative score descending
    introns.sort(key=lambda x: x['rel_score'], reverse=True)

    print(f"\nFound {len(introns)} predicted U12s")
    print(f"Analyzing top {n_top} highest-scoring introns:\n")

    # Analyze top N
    for i, intron in enumerate(introns[:n_top], 1):
        print(f"\n{'='*80}")
        print(f"#{i}: {intron['name'][:60]}")
        print(f"  Relative score: {intron['rel_score']:.2f}")
        print(f"  SVM probability: {intron['svm_score']:.1f}%")

        # Compute all features
        features = compute_features(intron['s5'], intron['sBP'], intron['s3'])

        # Print base features
        print(f"\n  Base z-scores:")
        print(f"    s5:  {features['s5']:>7.4f}")
        print(f"    sBP: {features['sBP']:>7.4f}")
        print(f"    s3:  {features['s3']:>7.4f}")

        # Print min features
        print(f"\n  Min features (require all strong):")
        print(f"    min_all:   {features['min_all']:>7.4f}  ← ALL THREE must be strong")
        print(f"    min_5_bp:  {features['min_5_bp']:>7.4f}")
        print(f"    min_5_3:   {features['min_5_3']:>7.4f}")
        print(f"    min_bp_3:  {features['min_bp_3']:>7.4f}")

        # Print max features
        print(f"\n  Max features (at least one strong):")
        print(f"    max_all:   {features['max_all']:>7.4f}")
        print(f"    max_5_bp:  {features['max_5_bp']:>7.4f}")
        print(f"    max_5_3:   {features['max_5_3']:>7.4f}")
        print(f"    max_bp_3:  {features['max_bp_3']:>7.4f}")

        # Print penalties
        print(f"\n  Imbalance penalties (should be small for real U12s):")
        print(f"    neg_absdiff_5_bp: {features['neg_absdiff_5_bp']:>7.4f}  ← 5'/BP balance")
        print(f"    neg_absdiff_5_3:  {features['neg_absdiff_5_3']:>7.4f}  ← 5'/3' balance")
        print(f"    neg_absdiff_bp_3: {features['neg_absdiff_bp_3']:>7.4f}  ← BP/3' balance")

        # Identify pattern
        print(f"\n  Pattern analysis:")
        negative_scores = [k for k, v in [('s5', features['s5']), ('sBP', features['sBP']), ('s3', features['s3'])] if v < 0]
        weak_scores = [k for k, v in [('s5', features['s5']), ('sBP', features['sBP']), ('s3', features['s3'])] if -0.5 < v < 0.5]

        if negative_scores:
            print(f"    ⚠ NEGATIVE scores: {', '.join(negative_scores)}")
        if len([f for f in [features['s5'], features['sBP'], features['s3']] if f > 1.0]) == 1:
            print(f"    ⚠ ONE-END-STRONG pattern detected")
        if abs(features['neg_absdiff_5_bp']) > 1.0 or abs(features['neg_absdiff_5_3']) > 1.0 or abs(features['neg_absdiff_bp_3']) > 1.0:
            print(f"    ⚠ High imbalance penalties (>1.0)")
        if features['min_all'] < 0:
            print(f"    ⚠ min_all is NEGATIVE (should strongly penalize)")

if __name__ == '__main__':
    score_file = Path('run_tests/caenorhabditis_elegans.pretrained.score_info.iic')

    if not score_file.exists():
        print(f"Error: {score_file} not found")
        sys.exit(1)

    analyze_score_file(score_file, n_top=10)
