#!/usr/bin/env python3
"""
Compare longest transcript selections using correct file understanding:
- .introns.iic has ALL unique introns (scored + omitted)
- Introns WITHOUT [i] tag are from longest isoform
- Introns WITH [i] tag are NOT from longest isoform
"""

import re
from collections import defaultdict, Counter

print("="*80)
print("COMPARING LONGEST TRANSCRIPT SELECTIONS")
print("="*80)

# Parse original .introns.iic to identify longest transcripts
print("\n[Step 1] Parsing original introns...")
original_longest = {}  # gene -> transcript
original_intron_counts = Counter()  # transcript -> intron count
original_all_transcripts = defaultdict(set)  # gene -> all transcripts

with open('orig_scored.introns.iic') as f:
    for line in f:
        # Format: OriSco-gene:ENSG00000181143@transcript:ENST00000397910-intron_1(83)
        match = re.search(r'gene:([A-Z0-9]+)@transcript:([A-Z0-9]+)-intron_\d+\((\d+)\)', line)
        if match:
            gene, transcript, family_size = match.groups()
            original_all_transcripts[gene].add(transcript)
            original_intron_counts[transcript] = int(family_size)

            # Check if [i] tag is present
            if ';[i]' not in line:
                # This intron is from the longest isoform
                if gene not in original_longest:
                    original_longest[gene] = transcript
                elif original_longest[gene] != transcript:
                    print(f"  WARNING: Gene {gene} has multiple 'longest' transcripts: {original_longest[gene]} and {transcript}")

print(f"  Identified {len(original_longest):,} genes with longest transcript")
print(f"  Total unique transcripts: {len(original_intron_counts):,}")
print(f"  Total genes: {len(original_all_transcripts):,}")

# Parse refactored - wait for current run or use latest available
print("\n[Step 2] Parsing refactored introns...")
import os
import time

# Wait for current run to finish writing
ref_introns_file = 'comparison_outputs/refactored/ref_final_comparison.introns.iic'
if not os.path.exists(ref_introns_file):
    print(f"  Waiting for {ref_introns_file} to be created...")
    for _ in range(30):  # Wait up to 30 seconds
        time.sleep(1)
        if os.path.exists(ref_introns_file):
            break
    else:
        print(f"  File not ready, using older output")
        ref_introns_file = 'comparison_outputs/refactored/ref_with_dupes.introns.iic'

print(f"  Using: {ref_introns_file}")

refactored_longest = {}  # gene -> transcript
refactored_intron_counts = Counter()
refactored_all_transcripts = defaultdict(set)

with open(ref_introns_file) as f:
    for line in f:
        # Format: transcript:ENST00000476247_1[i]
        parts = line.strip().split('\t')
        if len(parts) >= 8:
            name = parts[0]
            transcript_col = parts[6]
            gene_col = parts[7]

            # Extract transcript and gene
            match_tx = re.search(r'transcript:([A-Z0-9]+)', transcript_col)
            match_gene = re.search(r'gene:([A-Z0-9]+)', gene_col)
            match_name = re.search(r'transcript:([A-Z0-9]+)_\d+', name)

            if match_tx and match_gene and match_name:
                transcript = match_tx.group(1)
                gene = match_gene.group(1)
                refactored_all_transcripts[gene].add(transcript)

                # Count introns per transcript
                refactored_intron_counts[transcript] += 1

                # Check if marked as NOT longest isoform
                if '[i]' not in name:
                    # This is from longest isoform
                    if gene not in refactored_longest:
                        refactored_longest[gene] = transcript

print(f"  Identified {len(refactored_longest):,} genes with longest transcript")
print(f"  Total unique transcripts: {len(refactored_intron_counts):,}")
print(f"  Total genes: {len(refactored_all_transcripts):,}")

# Compare
print("\n[Step 3] Comparing longest transcript selections...")
genes_in_both = set(original_longest.keys()) & set(refactored_longest.keys())
print(f"  Genes present in both: {len(genes_in_both):,}")

different_longest = []
for gene in genes_in_both:
    orig_tx = original_longest[gene]
    ref_tx = refactored_longest[gene]
    if orig_tx != ref_tx:
        different_longest.append((gene, orig_tx, ref_tx))

print(f"  Genes with DIFFERENT longest transcript: {len(different_longest):,}")

if different_longest:
    print(f"\n  Sample (first 30):")
    print(f"  {'Gene':<18} {'Original TX':<18} {'Refactored TX':<18} {'Orig Introns':<15} {'Ref Introns':<15}")
    print(f"  {'-'*18} {'-'*18} {'-'*18} {'-'*15} {'-'*15}")
    for gene, orig_tx, ref_tx in different_longest[:30]:
        orig_count = original_intron_counts.get(orig_tx, 0)
        ref_count = refactored_intron_counts.get(ref_tx, 0)
        print(f"  {gene:<18} {orig_tx:<18} {ref_tx:<18} {orig_count:<15} {ref_count:<15}")

# Estimate intron impact
print("\n[Step 4] Estimating intron count impact...")
affected_introns = 0
for gene, orig_tx, ref_tx in different_longest:
    # Introns from orig_tx are now marked as [i] in refactored
    # Count how many introns orig_tx had
    affected_introns += original_intron_counts.get(orig_tx, 0)

print(f"  Introns from transcripts that are no longer 'longest': ~{affected_introns:,}")
print(f"  Target discrepancy: 232 introns")
print(f"  Ratio: {affected_introns / 232:.1f}x (some introns may be duplicates)")

print("\n" + "="*80)
print("CONCLUSION")
print("="*80)
print(f"\n{len(different_longest):,} genes have different 'longest' transcript selections.")
print(f"This affects approximately {affected_introns:,} intron assignments.")
print(f"\nThe discrepancy of 232 introns is likely caused by these {len(different_longest):,} genes")
print("where a different transcript was selected as 'longest' due to missing")
print("'line_number' tiebreaker in the hierarchical sort.")
