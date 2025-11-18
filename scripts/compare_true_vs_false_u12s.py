#!/usr/bin/env python3
"""
Compare true U12-type introns (H. sapiens) vs false positives (C. elegans).
"""
import statistics
from pathlib import Path

def parse_u12_scores(filepath):
    """Parse U12-type introns from score_info.iic file."""
    scores = []

    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 14:
                continue

            try:
                rel_score = float(parts[1])
                if rel_score <= 0:
                    continue  # Skip U2-type

                intron_name = parts[0]
                svm_score = float(parts[2])
                z_5prime = float(parts[5])
                z_bps = float(parts[9])
                z_3prime = float(parts[12])

                scores.append({
                    'name': intron_name,
                    'svm_score': svm_score,
                    'z_5prime': z_5prime,
                    'z_bps': z_bps,
                    'z_3prime': z_3prime,
                    'min_z': min(z_5prime, z_bps, z_3prime),
                    'max_z': max(z_5prime, z_bps, z_3prime),
                    'asymmetry': abs(z_5prime - z_3prime),
                    'sum_ends': z_5prime + z_3prime,
                    'all_positive': (z_5prime > 0 and z_bps > 0 and z_3prime > 0),
                    'weak_feature': (z_5prime < 0.5 or z_bps < 0.5 or z_3prime < 0.5)
                })
            except (ValueError, IndexError):
                continue

    return scores

def print_statistics(name, data):
    """Print statistics for a dataset."""
    if not data:
        print(f"No data for {name}")
        return

    print(f"\n{'='*80}")
    print(f"{name}")
    print(f"{'='*80}")
    print(f"Sample size: {len(data)}")

    # Z-score statistics
    z5_vals = [d['z_5prime'] for d in data]
    zbps_vals = [d['z_bps'] for d in data]
    z3_vals = [d['z_3prime'] for d in data]
    min_vals = [d['min_z'] for d in data]
    max_vals = [d['max_z'] for d in data]
    asym_vals = [d['asymmetry'] for d in data]
    sum_vals = [d['sum_ends'] for d in data]

    print(f"\nZ-score ranges (mean ± std):")
    print(f"  5' SS:  {statistics.mean(z5_vals):>7.3f} ± {statistics.stdev(z5_vals) if len(z5_vals) > 1 else 0:>6.3f}  "
          f"[{min(z5_vals):>7.3f}, {max(z5_vals):>7.3f}]")
    print(f"  BPS:    {statistics.mean(zbps_vals):>7.3f} ± {statistics.stdev(zbps_vals) if len(zbps_vals) > 1 else 0:>6.3f}  "
          f"[{min(zbps_vals):>7.3f}, {max(zbps_vals):>7.3f}]")
    print(f"  3' SS:  {statistics.mean(z3_vals):>7.3f} ± {statistics.stdev(z3_vals) if len(z3_vals) > 1 else 0:>6.3f}  "
          f"[{min(z3_vals):>7.3f}, {max(z3_vals):>7.3f}]")

    print(f"\nCombined metrics:")
    print(f"  Min z:      {statistics.mean(min_vals):>7.3f} ± {statistics.stdev(min_vals) if len(min_vals) > 1 else 0:>6.3f}  "
          f"[{min(min_vals):>7.3f}, {max(min_vals):>7.3f}]")
    print(f"  Max z:      {statistics.mean(max_vals):>7.3f} ± {statistics.stdev(max_vals) if len(max_vals) > 1 else 0:>6.3f}  "
          f"[{min(max_vals):>7.3f}, {max(max_vals):>7.3f}]")
    print(f"  Asymmetry:  {statistics.mean(asym_vals):>7.3f} ± {statistics.stdev(asym_vals) if len(asym_vals) > 1 else 0:>6.3f}  "
          f"[{min(asym_vals):>7.3f}, {max(asym_vals):>7.3f}]")
    print(f"  Sum (5'+3'): {statistics.mean(sum_vals):>7.3f} ± {statistics.stdev(sum_vals) if len(sum_vals) > 1 else 0:>6.3f}  "
          f"[{min(sum_vals):>7.3f}, {max(sum_vals):>7.3f}]")

    # Key diagnostic counts
    neg_count = sum(1 for d in data if d['min_z'] < 0)
    weak_count = sum(1 for d in data if d['weak_feature'])
    all_pos_count = sum(1 for d in data if d['all_positive'])
    high_asym = sum(1 for d in data if d['asymmetry'] > 2.0)
    extreme = sum(1 for d in data if d['max_z'] > 5 and d['min_z'] < 1)

    print(f"\nDiagnostic counts:")
    print(f"  At least ONE negative z-score:      {neg_count:>3} / {len(data)} ({100*neg_count/len(data):>5.1f}%)")
    print(f"  At least ONE weak feature (z<0.5):  {weak_count:>3} / {len(data)} ({100*weak_count/len(data):>5.1f}%)")
    print(f"  ALL three features positive (z>0):  {all_pos_count:>3} / {len(data)} ({100*all_pos_count/len(data):>5.1f}%)")
    print(f"  High asymmetry (|5'-3'| > 2.0):     {high_asym:>3} / {len(data)} ({100*high_asym/len(data):>5.1f}%)")
    print(f"  Extreme single-feature (max>5, min<1): {extreme:>3} / {len(data)} ({100*extreme/len(data):>5.1f}%)")

def main():
    # Parse data
    celegans_file = Path('archive/test_outputs/celegans/caenorhabditis_elegans.score_info.iic')
    human_file = Path('homo_sapiens.score_info.iic')

    if not celegans_file.exists():
        print(f"Error: {celegans_file} not found")
        return
    if not human_file.exists():
        print(f"Error: {human_file} not found")
        return

    fps = parse_u12_scores(celegans_file)
    true_u12s = parse_u12_scores(human_file)

    print("\n" + "="*80)
    print("COMPARISON: TRUE U12s (H. sapiens) vs FALSE POSITIVES (C. elegans)")
    print("="*80)

    # Print statistics for both
    print_statistics("TRUE U12-TYPE INTRONS (Homo sapiens)", true_u12s)
    print_statistics("FALSE POSITIVES (C. elegans - no real U12s)", fps)

    # Direct comparison
    print("\n" + "="*80)
    print("KEY DIFFERENCES")
    print("="*80)

    if true_u12s and fps:
        print(f"\n1. MINIMUM Z-SCORE (weakest feature):")
        true_min_mean = statistics.mean([d['min_z'] for d in true_u12s])
        fp_min_mean = statistics.mean([d['min_z'] for d in fps])
        print(f"   True U12s:       {true_min_mean:>7.3f} (average)")
        print(f"   False positives: {fp_min_mean:>7.3f} (average)")
        print(f"   → Difference:    {true_min_mean - fp_min_mean:>7.3f}")

        print(f"\n2. ALL FEATURES POSITIVE:")
        true_all_pos = sum(1 for d in true_u12s if d['all_positive']) / len(true_u12s) * 100
        fp_all_pos = sum(1 for d in fps if d['all_positive']) / len(fps) * 100
        print(f"   True U12s:       {true_all_pos:>5.1f}%")
        print(f"   False positives: {fp_all_pos:>5.1f}%")

        print(f"\n3. ASYMMETRY (|5' - 3'|):")
        true_asym_mean = statistics.mean([d['asymmetry'] for d in true_u12s])
        fp_asym_mean = statistics.mean([d['asymmetry'] for d in fps])
        print(f"   True U12s:       {true_asym_mean:>7.3f} (average)")
        print(f"   False positives: {fp_asym_mean:>7.3f} (average)")

        print(f"\n4. CONSISTENCY (sum of 5' + 3'):")
        true_sum_mean = statistics.mean([d['sum_ends'] for d in true_u12s])
        fp_sum_mean = statistics.mean([d['sum_ends'] for d in fps])
        print(f"   True U12s:       {true_sum_mean:>7.3f} (average)")
        print(f"   False positives: {fp_sum_mean:>7.3f} (average)")

        print(f"\n5. AT LEAST ONE WEAK FEATURE (z < 0.5):")
        true_weak = sum(1 for d in true_u12s if d['weak_feature']) / len(true_u12s) * 100
        fp_weak = sum(1 for d in fps if d['weak_feature']) / len(fps) * 100
        print(f"   True U12s:       {true_weak:>5.1f}%")
        print(f"   False positives: {fp_weak:>5.1f}% ← ALL OF THEM!")

    # Recommendations
    print("\n" + "="*80)
    print("RECOMMENDATIONS FOR REDUCING FALSE POSITIVES")
    print("="*80)

    print("\n1. **Add consistency features (as suggested by expert):**")
    print("   - Feature: |z_5' - z_3'| (asymmetry penalty)")
    print("   - Feature: z_5' + z_3' (consistency reward)")
    print("   → Linear model can learn to penalize one-end-strong introns")

    print("\n2. **Require minimum z-score threshold:**")
    if true_u12s and fps:
        true_min_5th = sorted([d['min_z'] for d in true_u12s])[int(0.05 * len(true_u12s))]
        print(f"   → 5th percentile for true U12s: {true_min_5th:.3f}")
        print(f"   → Could add post-filter: min_z > {max(0.0, true_min_5th):.1f}")

    print("\n3. **Zero-anchored robust spread (ZAR) normalization:**")
    print("   - Compute per-species robust spread: σ̂ = median(min(|s|, Q_0.995(|s|)))")
    print("   - Normalize: s̃ = s / σ̂ (no centering, zero stays zero)")
    print("   → Prevents species-specific score inflation")

    print("\n4. **Ensemble with feature importance:**")
    print("   - Current: Equal weight to all 3 features")
    print("   - Better: Weight BPS more heavily (strongest U12 signal)")
    print("   → Or use learned feature importance from training")

if __name__ == '__main__':
    main()
