#!/usr/bin/env python3
"""
Investigate why refactored generates 58,962 introns vs original's 58,933 (+29 difference).

Strategy:
1. Run both versions with duplicate inclusion (-d flag for original)
2. Count duplicates per coordinate in both versions
3. Identify coordinates with different duplicate counts
4. Analyze patterns: which transcripts, CDS vs exon, etc.
"""

import subprocess
import sys
from collections import defaultdict, Counter
from pathlib import Path

print("="*80)
print("DUPLICATE GENERATION INVESTIGATION")
print("="*80)

# Step 1: Generate outputs with ALL duplicates included
print("\n[Step 1] Running original with duplicate inclusion (-d flag)...")
result = subprocess.run(
    ["python", "intronIC/intronIC_patched.py",
     "-g", "intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz",
     "-a", "intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz",
     "-n", "orig_with_dupes",
     "-s", "-d"],  # sequences-only + include duplicates
    capture_output=True,
    text=True
)
if result.returncode != 0:
    print(f"ERROR: {result.stderr}")
    sys.exit(1)

# Extract intron count from log
orig_count = None
for line in result.stdout.split('\n'):
    if 'introns found' in line:
        # Format: "[#] [58933] introns found in..."
        import re
        match = re.search(r'\[(\d+)\].*introns found', line)
        if match:
            orig_count = int(match.group(1))
            break

print(f"✓ Original generated {orig_count:,} introns (including duplicates)")

print("\n[Step 2] Running refactored (duplicates included by default)...")
result = subprocess.run(
    ["python", "-m", "intronIC_refactored.cli.main",
     "-g", "intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz",
     "-a", "intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz",
     "-n", "ref_with_dupes",
     "-o", "comparison_outputs/refactored",
     "-s"],  # sequences-only (duplicates included by default)
    capture_output=True,
    text=True
)
if result.returncode != 0:
    print(f"ERROR: {result.stderr}")
    sys.exit(1)

# Extract intron count from log (check both stdout and stderr)
import re
ref_count = None
output_text = result.stdout + result.stderr
for line in output_text.split('\n'):
    if 'Generated' in line and 'introns' in line:
        # Format: "INFO: Generated 58962 introns"
        match = re.search(r'Generated (\d+) introns', line)
        if match:
            ref_count = int(match.group(1))
            break

if ref_count is None:
    print("ERROR: Could not extract intron count from refactored output")
    print("Output sample:")
    print(output_text[:500])
    sys.exit(1)

print(f"✓ Refactored generated {ref_count:,} introns (including duplicates)")

print(f"\n[Step 3] Analyzing duplicate patterns...")
print(f"  Difference: {ref_count - orig_count:+,} introns")

# Step 2: Parse BED files and count duplicates per coordinate
print("\n[Step 4] Counting duplicates per coordinate...")

def parse_bed_file(filename):
    """Parse BED file and group introns by coordinates."""
    coord_to_introns = defaultdict(list)

    with open(filename) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            chrom, start, stop, name, score, strand = parts[:6]
            coord = f"{chrom}:{start}-{stop}:{strand}"

            # Extract transcript and intron index from name
            # Original: "OriFin-gene:ENSG@transcript:ENST-intron_N(M)"
            # Refactored: "transcript:ENST_N[tags]"
            coord_to_introns[coord].append({
                'name': name,
                'line': line.strip()
            })

    return coord_to_introns

orig_coords = parse_bed_file('orig_with_dupes.bed.iic')
ref_coords = parse_bed_file('comparison_outputs/refactored/ref_with_dupes.bed.iic')

# Count duplicates per coordinate
print(f"\nOriginal coordinate groups: {len(orig_coords):,}")
print(f"Refactored coordinate groups: {len(ref_coords):,}")

# Find coordinates with different duplicate counts
print(f"\n[Step 5] Finding coordinates with different duplicate counts...")

diff_counts = []
for coord in set(orig_coords.keys()) | set(ref_coords.keys()):
    orig_count = len(orig_coords.get(coord, []))
    ref_count = len(ref_coords.get(coord, []))

    if orig_count != ref_count:
        diff_counts.append({
            'coord': coord,
            'orig_count': orig_count,
            'ref_count': ref_count,
            'diff': ref_count - orig_count
        })

diff_counts.sort(key=lambda x: abs(x['diff']), reverse=True)

print(f"\nFound {len(diff_counts)} coordinates with different duplicate counts")
print(f"\nTop 20 coordinates with most difference:")
print(f"{'Coordinate':<40} {'Original':<10} {'Refactored':<12} {'Diff':<8}")
print("-"*80)

for item in diff_counts[:20]:
    print(f"{item['coord']:<40} {item['orig_count']:<10} {item['ref_count']:<12} {item['diff']:+8}")

# Analyze which have more duplicates in refactored
more_in_ref = [d for d in diff_counts if d['diff'] > 0]
more_in_orig = [d for d in diff_counts if d['diff'] < 0]

print(f"\nSummary:")
print(f"  Coordinates with MORE duplicates in refactored: {len(more_in_ref)}")
print(f"  Coordinates with FEWER duplicates in refactored: {len(more_in_orig)}")
print(f"  Total extra introns from 'more': {sum(d['diff'] for d in more_in_ref):+,}")
print(f"  Total fewer introns from 'less': {sum(d['diff'] for d in more_in_orig):+,}")
print(f"  Net difference: {sum(d['diff'] for d in diff_counts):+,}")

# Step 3: Analyze specific examples
print(f"\n[Step 6] Analyzing specific examples...")

if more_in_ref:
    print(f"\nExample coordinate with MORE duplicates in refactored:")
    example = more_in_ref[0]
    coord = example['coord']

    print(f"\nCoordinate: {coord}")
    print(f"  Original count: {example['orig_count']}")
    print(f"  Refactored count: {example['ref_count']}")

    print(f"\n  Original introns:")
    for intron in orig_coords.get(coord, [])[:5]:
        print(f"    {intron['name'][:80]}")

    print(f"\n  Refactored introns:")
    for intron in ref_coords.get(coord, [])[:5]:
        print(f"    {intron['name'][:80]}")

if more_in_orig:
    print(f"\nExample coordinate with FEWER duplicates in refactored:")
    example = more_in_orig[0]
    coord = example['coord']

    print(f"\nCoordinate: {coord}")
    print(f"  Original count: {example['orig_count']}")
    print(f"  Refactored count: {example['ref_count']}")

    print(f"\n  Original introns:")
    for intron in orig_coords.get(coord, [])[:5]:
        print(f"    {intron['name'][:80]}")

    print(f"\n  Refactored introns:")
    for intron in ref_coords.get(coord, [])[:5]:
        print(f"    {intron['name'][:80]}")

# Step 4: Check if related to 1bp features
print(f"\n[Step 7] Checking if differences involve very short introns...")

very_short_diffs = [d for d in diff_counts if 'diff' in d]
short_coords = []

for item in diff_counts[:50]:  # Check top 50
    coord = item['coord']
    try:
        start, stop = coord.split(':')[1].split('-')
        length = int(stop) - int(start)
        if length <= 10:
            short_coords.append((coord, length, item['diff']))
    except:
        pass

if short_coords:
    print(f"\nVery short introns (≤10bp) with duplicate count differences:")
    for coord, length, diff in short_coords[:10]:
        print(f"  {coord:<40} {length:3}bp  diff={diff:+3}")

print("\n" + "="*80)
print("Investigation complete. Check output files:")
print("  - orig_with_dupes.bed.iic")
print("  - comparison_outputs/refactored/ref_with_dupes.bed.iic")
print("="*80)
