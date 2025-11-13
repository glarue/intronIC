#!/usr/bin/env python3
"""
Analyze false positive patterns in C. elegans compared to reference U12 introns.
"""
import gzip
import sys
from pathlib import Path
import statistics

def parse_score_file(filepath):
    """Parse .score_info.iic file and extract z-scores."""
    scores = []

    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 14:
                continue

            try:
                intron_name = parts[0]
                rel_score = float(parts[1])
                svm_score = float(parts[2])

                # Extract z-scores (columns 6, 10, 13)
                z_5prime = float(parts[5])
                z_bps = float(parts[9])
                z_3prime = float(parts[12])

                scores.append({
                    'name': intron_name,
                    'rel_score': rel_score,
                    'svm_score': svm_score,
                    'z_5prime': z_5prime,
                    'z_bps': z_bps,
                    'z_3prime': z_3prime,
                    'min_z': min(z_5prime, z_bps, z_3prime),
                    'max_z': max(z_5prime, z_bps, z_3prime),
                    'asymmetry': abs(z_5prime - z_3prime),
                    'sum_ends': z_5prime + z_3prime,
                })
            except (ValueError, IndexError):
                continue

    return scores

def analyze_celegans_fps():
    """Analyze C. elegans false positives."""
    celegans_file = Path('archive/test_outputs/celegans/caenorhabditis_elegans.score_info.iic')

    if not celegans_file.exists():
        print(f"Error: {celegans_file} not found")
        return

    all_scores = parse_score_file(celegans_file)

    # Separate U12-type calls (false positives) from U2-type
    fps = [s for s in all_scores if s['rel_score'] > 0]
    u2s = [s for s in all_scores if s['rel_score'] <= 0]

    print("=" * 80)
    print("C. elegans False Positive Analysis")
    print("=" * 80)
    print(f"\nTotal introns: {len(all_scores):,}")
    print(f"Called as U12-type (FALSE POSITIVES): {len(fps)}")
    print(f"Called as U2-type: {len(u2s):,}")
    print(f"False positive rate: {100 * len(fps) / len(all_scores):.4f}%")

    print("\n" + "-" * 80)
    print("FALSE POSITIVE SCORE PATTERNS:")
    print("-" * 80)
    print(f"{'Intron':<45} {'SVM%':>6} {'5z':>7} {'BPz':>7} {'3z':>7} {'MinZ':>7} {'|5-3|':>7}")
    print("-" * 80)

    for fp in sorted(fps, key=lambda x: x['svm_score'], reverse=True):
        print(f"{fp['name'][:45]:<45} {fp['svm_score']:>6.1f} "
              f"{fp['z_5prime']:>7.2f} {fp['z_bps']:>7.2f} {fp['z_3prime']:>7.2f} "
              f"{fp['min_z']:>7.2f} {fp['asymmetry']:>7.2f}")

    print("\n" + "=" * 80)
    print("KEY OBSERVATIONS:")
    print("=" * 80)

    # Count introns with negative minimum z-score
    neg_min = [fp for fp in fps if fp['min_z'] < 0]
    print(f"\n1. Introns with at least ONE NEGATIVE z-score: {len(neg_min)}/{len(fps)}")
    for fp in neg_min:
        which_neg = []
        if fp['z_5prime'] < 0:
            which_neg.append(f"5'={fp['z_5prime']:.2f}")
        if fp['z_bps'] < 0:
            which_neg.append(f"BPS={fp['z_bps']:.2f}")
        if fp['z_3prime'] < 0:
            which_neg.append(f"3'={fp['z_3prime']:.2f}")
        print(f"   - {fp['name'][:50]}: {', '.join(which_neg)}")

    # Count introns with weak minimum z-score (< 0.5)
    weak_min = [fp for fp in fps if fp['min_z'] < 0.5]
    print(f"\n2. Introns with WEAK minimum z-score (< 0.5): {len(weak_min)}/{len(fps)}")

    # Count introns with high asymmetry
    high_asym = [fp for fp in fps if fp['asymmetry'] > 2.0]
    print(f"\n3. Introns with HIGH ASYMMETRY (|5'-3'| > 2.0): {len(high_asym)}/{len(fps)}")
    for fp in high_asym:
        print(f"   - {fp['name'][:50]}: 5'={fp['z_5prime']:.2f}, 3'={fp['z_3prime']:.2f}, diff={fp['asymmetry']:.2f}")

    # Extreme single-feature cases
    extreme = [fp for fp in fps if fp['max_z'] > 5 and fp['min_z'] < 1]
    print(f"\n4. EXTREME single-feature cases (max>5, min<1): {len(extreme)}/{len(fps)}")
    for fp in extreme:
        print(f"   - {fp['name'][:50]}: 5'={fp['z_5prime']:.2f}, BPS={fp['z_bps']:.2f}, 3'={fp['z_3prime']:.2f}")

    print("\n" + "=" * 80)
    print("STATISTICS:")
    print("=" * 80)

    if fps:
        print(f"\nMean z-scores (false positives):")
        z5_vals = [fp['z_5prime'] for fp in fps]
        zbps_vals = [fp['z_bps'] for fp in fps]
        z3_vals = [fp['z_3prime'] for fp in fps]
        min_vals = [fp['min_z'] for fp in fps]
        asym_vals = [fp['asymmetry'] for fp in fps]

        print(f"  5' SS:  {statistics.mean(z5_vals):.3f} ± {statistics.stdev(z5_vals) if len(z5_vals) > 1 else 0:.3f}")
        print(f"  BPS:    {statistics.mean(zbps_vals):.3f} ± {statistics.stdev(zbps_vals) if len(zbps_vals) > 1 else 0:.3f}")
        print(f"  3' SS:  {statistics.mean(z3_vals):.3f} ± {statistics.stdev(z3_vals) if len(z3_vals) > 1 else 0:.3f}")
        print(f"\nMin z-score: {statistics.mean(min_vals):.3f} ± {statistics.stdev(min_vals) if len(min_vals) > 1 else 0:.3f}")
        print(f"Asymmetry:   {statistics.mean(asym_vals):.3f} ± {statistics.stdev(asym_vals) if len(asym_vals) > 1 else 0:.3f}")

    return fps, u2s

if __name__ == '__main__':
    analyze_celegans_fps()
