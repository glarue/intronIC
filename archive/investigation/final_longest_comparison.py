#!/usr/bin/env python3
"""
Final comparison of longest transcript selections.
Uses the correct file formats:
- Original: .introns.iic (all introns), check for [i] tag
- Refactored: .meta.iic (all introns), check column 15 for "not_longest_isoform"
"""

import re
from collections import defaultdict, Counter

print("="*80)
print("FINAL LONGEST TRANSCRIPT COMPARISON")
print("="*80)

# Parse original
print("\n[Step 1] Parsing original introns (orig_scored.introns.iic)...")
original_longest = {}  # gene -> transcript
original_all_introns = []  # (gene, transcript, index)
original_gene_transcripts = defaultdict(set)

with open('orig_scored.introns.iic') as f:
    for line in f:
        # Format: OriSco-gene:ENSG00000181143@transcript:ENST00000397910-intron_1(83)
        match = re.search(r'gene:([A-Z0-9]+)@transcript:([A-Z0-9]+)-intron_(\d+)', line)
        if match:
            gene, transcript, index = match.groups()
            original_all_introns.append((gene, transcript, index))
            original_gene_transcripts[gene].add(transcript)

            # Check if this is from longest isoform
            if ';[i]' not in line:  # No [i] tag = longest isoform
                if gene not in original_longest:
                    original_longest[gene] = transcript
                elif original_longest[gene] != transcript:
                    # Multiple transcripts from same gene without [i] tag - shouldn't happen
                    pass

print(f"  Total introns: {len(original_all_introns):,}")
print(f"  Genes with longest identified: {len(original_longest):,}")
print(f"  Total genes: {len(original_gene_transcripts):,}")

# Parse refactored
print("\n[Step 2] Parsing refactored introns (ref_seqs_only.meta.iic)...")
refactored_longest = {}  # gene -> transcript
refactored_all_introns = []  # (gene, transcript, index)
refactored_gene_transcripts = defaultdict(set)

with open('comparison_outputs/refactored/ref_seqs_only.meta.iic') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 15:
            name = parts[0]  # transcript:ENST00000477565_1
            transcript_col = parts[6]  # transcript:ENST00000477565
            gene_col = parts[7]  # gene:ENSG00000254901
            omission_col = parts[14]  # "not_longest_isoform" or empty

            # Extract IDs
            match_name = re.search(r'transcript:([A-Z0-9]+)_(\d+)', name)
            match_gene = re.search(r'gene:([A-Z0-9]+)', gene_col)

            if match_name and match_gene:
                transcript, index = match_name.groups()
                gene = match_gene.group(1)

                refactored_all_introns.append((gene, transcript, index))
                refactored_gene_transcripts[gene].add(transcript)

                # Check if from longest isoform
                if omission_col != 'not_longest_isoform':
                    if gene not in refactored_longest:
                        refactored_longest[gene] = transcript

print(f"  Total introns: {len(refactored_all_introns):,}")
print(f"  Genes with longest identified: {len(refactored_longest):,}")
print(f"  Total genes: {len(refactored_gene_transcripts):,}")

# Compare
print("\n[Step 3] Comparing longest transcript selections...")
genes_in_both = set(original_longest.keys()) & set(refactored_longest.keys())
print(f"  Genes present in both: {len(genes_in_both):,}")

different_longest = []
for gene in sorted(genes_in_both):
    orig_tx = original_longest[gene]
    ref_tx = refactored_longest[gene]
    if orig_tx != ref_tx:
        different_longest.append((gene, orig_tx, ref_tx))

print(f"  Genes with DIFFERENT longest transcript: {len(different_longest):,}")

if different_longest:
    print(f"\n  Sample (first 40):")
    print(f"  {'Gene':<20} {'Original TX':<20} {'Refactored TX':<20}")
    print(f"  {'-'*20} {'-'*20} {'-'*20}")
    for gene, orig_tx, ref_tx in different_longest[:40]:
        print(f"  {gene:<20} {orig_tx:<20} {ref_tx:<20}")

# Count affected introns
print("\n[Step 4] Counting affected introns...")
# For each gene with different selection, count introns from the transcript
# that is marked as longest in original but NOT in refactored
affected_intron_count = 0
for gene, orig_tx, ref_tx in different_longest:
    # Count introns from orig_tx (these are now marked as [i] in refactored)
    count = sum(1 for g, t, idx in original_all_introns if g == gene and t == orig_tx)
    affected_intron_count += count

print(f"  Introns affected by different longest selection: {affected_intron_count:,}")
print(f"  Target discrepancy: 232 introns (more marked as not_longest in refactored)")
print(f"  Match: {affected_intron_count == 232}")

print("\n" + "="*80)
print("ANALYSIS")
print("="*80)
print(f"\n{len(different_longest):,} genes have different 'longest' transcript selections.")
print(f"This causes {affected_intron_count:,} introns to be incorrectly marked as 'not longest isoform'.")
print(f"\nRoot cause: Missing 'line_number' tiebreaker in hierarchical sort.")
print(f"When transcripts have identical (defined_by, parent_length, parent, family_size, index),")
print(f"Python's sort order may differ from original, causing different 'first seen' assignments.")
