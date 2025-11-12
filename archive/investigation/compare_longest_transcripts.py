#!/usr/bin/env python3
"""
Compare which transcripts are selected as "longest" between original and refactored.

This should reveal why 232 more introns are being marked as "not longest isoform".
"""

import re
from collections import defaultdict, Counter

print("="*80)
print("COMPARING LONGEST TRANSCRIPT SELECTIONS")
print("="*80)

# Parse original to identify which transcripts are "longest" per gene
print("\n[Step 1] Identifying longest transcripts in original...")
original_longest = {}  # gene -> transcript
original_gene_transcripts = defaultdict(set)  # gene -> set of all transcripts

# From orig_scored.bed.iic, introns NOT marked with [i] are from longest isoform
with open('orig_scored.bed.iic') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 4:
            name = parts[3]
            # Parse: OriSco-gene:ENSG00000196218@transcript:ENST00000355481-intron_69(104);[o:i];[i]
            match = re.search(r'gene:([^@]+)@transcript:([^-]+)', name)
            if match:
                gene, transcript = match.groups()
                original_gene_transcripts[gene].add(transcript)

                # Check if [i] tag is present (NOT longest isoform)
                if ';[i]' not in name:
                    # This is from longest isoform
                    if gene not in original_longest:
                        original_longest[gene] = transcript
                    elif original_longest[gene] != transcript:
                        # Multiple transcripts marked as longest? Shouldn't happen
                        pass

print(f"  Found {len(original_longest):,} genes with longest transcript identified")
print(f"  Total genes: {len(original_gene_transcripts):,}")

# Parse refactored to identify longest transcripts
print("\n[Step 2] Identifying longest transcripts in refactored...")
refactored_longest = {}  # gene -> transcript
refactored_gene_transcripts = defaultdict(set)  # gene -> set of all transcripts

# Wait for the current run to finish, or use ref_with_dupes which has the hierarchical sorting
# Let me check if output is ready
import os
if os.path.exists('comparison_outputs/refactored/ref_final_comparison.meta.iic'):
    meta_file = 'comparison_outputs/refactored/ref_final_comparison.meta.iic'
elif os.path.exists('comparison_outputs/refactored/ref_with_dupes.meta.iic'):
    meta_file = 'comparison_outputs/refactored/ref_with_dupes.meta.iic'
    print("  Using ref_with_dupes (may not have latest fixes)")
else:
    print("  ERROR: No refactored meta file found")
    exit(1)

print(f"  Using: {meta_file}")

with open(meta_file) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 15:
            name = parts[0]
            transcript_col = parts[6]  # transcript:ENST...
            gene_col = parts[7]  # gene:ENSG...
            omission_col = parts[14]  # omission reason

            # Extract IDs
            match_tx = re.search(r'transcript:([A-Z0-9]+)', transcript_col)
            match_gene = re.search(r'gene:([A-Z0-9]+)', gene_col)

            if match_tx and match_gene:
                transcript = match_tx.group(1)
                gene = match_gene.group(1)
                refactored_gene_transcripts[gene].add(transcript)

                # Check if marked as longest isoform
                if 'not_longest_isoform' not in omission_col:
                    # This is from longest isoform (or at least not marked as NOT longest)
                    if gene not in refactored_longest:
                        refactored_longest[gene] = transcript

print(f"  Found {len(refactored_longest):,} genes with longest transcript identified")
print(f"  Total genes: {len(refactored_gene_transcripts):,}")

# Compare
print("\n[Step 3] Comparing longest transcript selections...")
genes_in_both = set(original_longest.keys()) & set(refactored_longest.keys())
print(f"  Genes in both: {len(genes_in_both):,}")

different_selections = []
for gene in genes_in_both:
    orig_tx = original_longest[gene]
    ref_tx = refactored_longest[gene]
    if orig_tx != ref_tx:
        different_selections.append((gene, orig_tx, ref_tx))

print(f"  Genes with DIFFERENT longest transcript: {len(different_selections):,}")

if different_selections:
    print(f"\n  Sample of genes with different selections (first 20):")
    print(f"  {'Gene':<20} {'Original':<20} {'Refactored':<20}")
    print(f"  {'-'*20} {'-'*20} {'-'*20}")
    for gene, orig_tx, ref_tx in different_selections[:20]:
        print(f"  {gene:<20} {orig_tx:<20} {ref_tx:<20}")

# Count how many introns this affects
print("\n[Step 4] Estimating impact on intron counts...")
affected_intron_count = 0

# For each gene with different selection, count introns in the non-selected transcripts
for gene, orig_tx, ref_tx in different_selections:
    # In refactored, introns from orig_tx are now marked as "not longest"
    # We need to count how many introns are in orig_tx
    # This is approximate without full intron data
    pass

print(f"\n  Genes with different longest transcript: {len(different_selections):,}")
print(f"  This could affect the ~232 intron discrepancy")

# Analyze patterns
print("\n[Step 5] Analyzing patterns...")
# Are there genes where multiple transcripts have same length?
# We'd need parent_length data to check this

print("\n" + "="*80)
print("HYPOTHESIS")
print("="*80)
print(f"\n{len(different_selections):,} genes have different 'longest' transcript selections.")
print("\nThis suggests the hierarchical sort is producing a different order")
print("for transcripts with identical sort key attributes.")
print("\nWithout 'line_number' as a final tiebreaker, Python's sort may")
print("order these transcripts differently than the original, leading to")
print("different 'first seen' assignments for ~232 introns across these genes.")
