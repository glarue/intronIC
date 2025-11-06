#!/usr/bin/env python3
"""
Investigate the remaining 236 intron discrepancy in scoring.

Strategy:
1. Compare which introns are scored vs omitted in each version
2. Identify specific genes/transcripts causing the difference
3. Analyze sorting/tagging patterns
"""

import re
from collections import defaultdict, Counter
from pathlib import Path

print("="*80)
print("INVESTIGATING 236 INTRON DISCREPANCY")
print("="*80)

# Parse original scored introns
print("\n[Step 1] Loading original scored introns...")
original_scored = set()
original_omitted_reasons = {}

with open('orig_scored.bed.iic') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 6:
            # Extract intron identifier from name field
            # Format: OriSco-gene:ENSG00000196218@transcript:ENST00000355481-intron_69(104);[o:i];[i]
            name = parts[3]
            # Parse to get gene, transcript, intron index
            match = re.search(r'gene:([^@]+)@transcript:([^-]+)-intron_(\d+)\((\d+)\)', name)
            if match:
                gene, transcript, intron_idx, family_size = match.groups()
                intron_id = f"{gene}:{transcript}:{intron_idx}"
                original_scored.add(intron_id)

print(f"  Found {len(original_scored):,} scored introns in original")

# Parse original omitted introns from meta file
with open('orig_scored.meta.iic') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 1:
            name = parts[0]
            # Check for omission code
            if '[o:' in name:
                match_omit = re.search(r'\[o:([a-z])\]', name)
                if match_omit:
                    omit_code = match_omit.group(1)
                    # Parse intron ID
                    match = re.search(r'gene:([^@]+)@transcript:([^-]+)-intron_(\d+)\((\d+)\)', name)
                    if match:
                        gene, transcript, intron_idx, family_size = match.groups()
                        intron_id = f"{gene}:{transcript}:{intron_idx}"
                        original_omitted_reasons[intron_id] = omit_code

print(f"  Found {len(original_omitted_reasons):,} omitted introns in original")

# Parse refactored scored introns
print("\n[Step 2] Loading refactored scored introns...")
refactored_scored = set()
refactored_omitted_reasons = {}

# Try the most recent test output
test_files = [
    'comparison_outputs/refactored/ref_test_isfals.bed.iic',
    'comparison_outputs/refactored/ref_scored_final.bed.iic',
    'comparison_outputs/refactored/ref_scored_v3.bed.iic'
]

ref_bed_file = None
for f in test_files:
    if Path(f).exists():
        ref_bed_file = f
        break

if not ref_bed_file:
    print("  ERROR: No refactored BED file found!")
    print(f"  Tried: {test_files}")
    exit(1)

print(f"  Using: {ref_bed_file}")

with open(ref_bed_file) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 6:
            name = parts[3]
            # Parse refactored format: transcript:ENST00000476247_1[i]
            # Need to extract gene from metadata file
            match = re.search(r'transcript:([^_]+)_(\d+)', name)
            if match:
                transcript, intron_idx = match.groups()
                # We'll need to map transcript->gene
                refactored_scored.add(f"{transcript}:{intron_idx}")

print(f"  Found {len(refactored_scored):,} scored introns in refactored")

# Parse refactored meta for gene mapping and omission reasons
print("\n[Step 3] Loading refactored metadata...")
ref_meta_file = ref_bed_file.replace('.bed.', '.meta.')
transcript_to_gene = {}
refactored_omitted_full = {}

if Path(ref_meta_file).exists():
    with open(ref_meta_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 14:
                name = parts[0]
                transcript_col = parts[6]  # transcript:ENST...
                gene_col = parts[7]  # gene:ENSG...
                omission_col = parts[14]  # last column with omission info

                # Extract transcript and gene IDs
                match_tx = re.search(r'transcript:([^_]+)_(\d+)', name)
                match_gene = re.search(r'gene:(.+)', gene_col)

                if match_tx and match_gene:
                    transcript, intron_idx = match_tx.groups()
                    gene = match_gene.group(1)
                    transcript_to_gene[transcript] = gene

                    intron_id = f"{gene}:{transcript}:{intron_idx}"

                    # Check omission reason
                    if 'not_longest_isoform' in omission_col:
                        refactored_omitted_reasons[intron_id] = 'i'
                    elif 'short' in omission_col:
                        refactored_omitted_reasons[intron_id] = 's'

print(f"  Mapped {len(transcript_to_gene):,} transcripts to genes")
print(f"  Found {len(refactored_omitted_reasons):,} omitted introns in refactored")

# Convert refactored scored set to use full IDs with genes
print("\n[Step 4] Reconciling intron IDs...")
refactored_scored_full = set()
for tx_intron in refactored_scored:
    transcript, intron_idx = tx_intron.split(':')
    if transcript in transcript_to_gene:
        gene = transcript_to_gene[transcript]
        intron_id = f"{gene}:{transcript}:{intron_idx}"
        refactored_scored_full.add(intron_id)

print(f"  Reconciled {len(refactored_scored_full):,} refactored scored introns")

# Compare
print("\n[Step 5] Comparing scored introns...")
only_in_original = original_scored - refactored_scored_full
only_in_refactored = refactored_scored_full - original_scored
in_both = original_scored & refactored_scored_full

print(f"  In both: {len(in_both):,}")
print(f"  Only in original (should be scored but isn't): {len(only_in_original):,}")
print(f"  Only in refactored (scored but shouldn't be): {len(only_in_refactored):,}")

# Analyze the differences
print("\n[Step 6] Analyzing introns scored in original but not refactored...")
if only_in_original:
    # Group by gene
    by_gene = defaultdict(list)
    for intron_id in only_in_original:
        gene = intron_id.split(':')[0]
        by_gene[gene].append(intron_id)

    print(f"  Affected genes: {len(by_gene)}")
    print(f"\n  Top 10 genes by affected introns:")
    for gene, introns in sorted(by_gene.items(), key=lambda x: len(x[1]), reverse=True)[:10]:
        print(f"    {gene}: {len(introns)} introns")
        for intron in sorted(introns)[:3]:
            _, transcript, idx = intron.split(':')
            reason = refactored_omitted_reasons.get(intron, 'unknown')
            print(f"      {transcript} intron {idx} - refactored marked as [{reason}]")

    # Check omission reasons
    print(f"\n  Omission reasons in refactored for these introns:")
    reason_counts = Counter()
    for intron_id in only_in_original:
        reason = refactored_omitted_reasons.get(intron_id, 'not_found')
        reason_counts[reason] += 1

    for reason, count in reason_counts.most_common():
        print(f"    [{reason}]: {count}")

print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print(f"\nOriginal scores: {len(original_scored):,}")
print(f"Refactored scores: {len(refactored_scored_full):,}")
print(f"Difference: {len(original_scored) - len(refactored_scored_full):,}")
print(f"\nIntrons scored in original but not refactored: {len(only_in_original):,}")
print("These introns are being incorrectly marked as omitted in the refactored version.")
print("\nNext step: Investigate why these specific introns are marked as 'not longest isoform'")
print("when they should be kept.")
