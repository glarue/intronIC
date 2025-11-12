#!/usr/bin/env python3
"""
Investigate why refactored scores 7,691 introns vs original's 12,074 (~4,400 difference).

Strategy:
1. Run both versions WITH scoring enabled
2. Parse meta files to see omission codes
3. Compare which introns are scored vs omitted
4. Identify patterns in the ~4,400 missing scored introns
"""

import subprocess
import sys
import re
from collections import Counter, defaultdict
from pathlib import Path

print("="*80)
print("SCORING DISCREPANCY INVESTIGATION")
print("="*80)

# Step 1: Run original WITH scoring (no -s flag)
print("\n[Step 1] Running original WITH scoring...")
result = subprocess.run(
    ["python", "intronIC/intronIC_patched.py",
     "-g", "intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz",
     "-a", "intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz",
     "-n", "orig_scored"],
    capture_output=True,
    text=True,
    timeout=300
)

if result.returncode != 0:
    print(f"ERROR: {result.stderr}")
    sys.exit(1)

# Extract counts from original output
orig_total = None
orig_scored = None
for line in result.stdout.split('\n'):
    if 'introns found' in line:
        match = re.search(r'\[(\d+)\].*introns found', line)
        if match:
            orig_total = int(match.group(1))
    if 'introns kept for scoring' in line or 'unique introns for scoring' in line:
        match = re.search(r'\[(\d+)\]', line)
        if match:
            orig_scored = int(match.group(1))

print(f"✓ Original: {orig_total:,} total → {orig_scored:,} scored")

# Step 2: Run refactored WITH scoring (no -s flag)
print("\n[Step 2] Running refactored WITH scoring...")
result = subprocess.run(
    ["python", "-m", "intronIC_refactored.cli.main",
     "-g", "intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz",
     "-a", "intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz",
     "-n", "ref_scored",
     "-o", "comparison_outputs/refactored"],
    capture_output=True,
    text=True,
    timeout=300
)

if result.returncode != 0:
    print(f"ERROR: {result.stderr}")
    sys.exit(1)

# Extract counts from refactored output
ref_total = None
ref_scored = None
output_text = result.stdout + result.stderr
for line in output_text.split('\n'):
    if 'Generated' in line and 'introns' in line:
        match = re.search(r'Generated (\d+) introns', line)
        if match:
            ref_total = int(match.group(1))
    if 'Loaded' in line and 'introns' in line:
        match = re.search(r'Loaded (\d+) introns', line)
        if match:
            ref_scored = int(match.group(1))

print(f"✓ Refactored: {ref_total:,} total → {ref_scored:,} scored")

# Calculate difference
print(f"\n[Step 3] Analyzing difference...")
print(f"  Difference in scored introns: {orig_scored - ref_scored:,}")

# Step 3: Parse meta files to analyze omission patterns
print(f"\n[Step 4] Parsing metadata files to analyze omission patterns...")

def parse_meta_file(filename):
    """Parse meta.iic file and categorize introns."""
    scored_coords = set()
    omitted_by_code = defaultdict(set)

    with open(filename) as f:
        header = f.readline()  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue

            name = parts[0]
            # Extract coordinates from name
            # Format: "species-gene:ID@transcript:ID-intron_N(M);[tags];type"
            # or "transcript:ID_N[tags]"

            # Try to extract coordinate info from other columns if available
            # For now, use the name field

            # Check if intron was omitted by looking for omission tags
            # Original: [i]=isoform, [n]=noncanonical, [s]=short, etc.
            # Refactored: Similar tags in name

            omitted = False
            omission_code = None

            # Look for omission indicators in the name
            if '[i]' in name:
                omission_code = 'i'  # not longest isoform
                omitted = True
            elif 'not_longest_isoform' in line:
                omission_code = 'i'
                omitted = True

            if '[n]' in name:
                omission_code = 'n'  # noncanonical
                omitted = True
            elif 'noncanonical' in line.lower():
                if omission_code != 'i':  # Don't override if already marked
                    omission_code = 'n'
                omitted = True

            if '[s]' in name:
                omission_code = 's'  # short
                omitted = True

            if not omitted:
                scored_coords.add(name)
            else:
                omitted_by_code[omission_code or 'unknown'].add(name)

    return scored_coords, omitted_by_code

# Parse original metadata
orig_scored_coords, orig_omitted = parse_meta_file('orig_scored.meta.iic')
print(f"\nOriginal omission breakdown:")
for code, introns in sorted(orig_omitted.items()):
    code_name = {
        'i': 'not longest isoform',
        'n': 'noncanonical',
        's': 'short',
        'a': 'ambiguous sequence',
        'v': 'coordinate overlap',
        'unknown': 'unknown'
    }.get(code, code)
    print(f"  [{code}] {code_name}: {len(introns):,}")
print(f"  Scored (no omission): {len(orig_scored_coords):,}")

# Parse refactored metadata
ref_scored_coords, ref_omitted = parse_meta_file('comparison_outputs/refactored/ref_scored.meta.iic')
print(f"\nRefactored omission breakdown:")
for code, introns in sorted(ref_omitted.items()):
    code_name = {
        'i': 'not longest isoform',
        'n': 'noncanonical',
        's': 'short',
        'a': 'ambiguous sequence',
        'v': 'coordinate overlap',
        'unknown': 'unknown'
    }.get(code, code)
    print(f"  [{code}] {code_name}: {len(introns):,}")
print(f"  Scored (no omission): {len(ref_scored_coords):,}")

# Compare omission patterns
print(f"\n[Step 5] Comparing omission patterns...")
print(f"\n{'Category':<30} {'Original':<12} {'Refactored':<12} {'Diff':<10}")
print("-"*70)

all_codes = set(orig_omitted.keys()) | set(ref_omitted.keys()) | {'scored'}
for code in sorted(all_codes):
    if code == 'scored':
        orig_count = len(orig_scored_coords)
        ref_count = len(ref_scored_coords)
        label = "Scored (no omission)"
    else:
        orig_count = len(orig_omitted.get(code, set()))
        ref_count = len(ref_omitted.get(code, set()))
        code_name = {
            'i': '[i] not longest isoform',
            'n': '[n] noncanonical',
            's': '[s] short',
            'a': '[a] ambiguous sequence',
            'v': '[v] coordinate overlap',
        }.get(code, f'[{code}]')
        label = code_name

    diff = ref_count - orig_count
    print(f"{label:<30} {orig_count:<12,} {ref_count:<12,} {diff:+10,}")

print(f"\n[Step 6] Summary")
print(f"  Original total: {orig_total:,}")
print(f"  Refactored total: {ref_total:,}")
print(f"  Original scored: {orig_scored:,}")
print(f"  Refactored scored: {ref_scored:,}")
print(f"  Scoring gap: {orig_scored - ref_scored:,} fewer introns scored in refactored")

print("\n" + "="*80)
print("Investigation complete. Check output files:")
print("  - orig_scored.meta.iic")
print("  - comparison_outputs/refactored/ref_scored.meta.iic")
print("="*80)
