#!/usr/bin/env python3
"""
Analyze the 30 missing introns between original and refactored intronIC.
Generate a detailed comparison report.
"""

# First, generate fresh outputs with duplicates included
import subprocess
import sys

print("=" * 80)
print("GENERATING FRESH OUTPUTS")
print("=" * 80)

# Original with duplicates
print("\n[1/2] Running original intronIC with duplicates...")
result = subprocess.run(
    ["python", "/home/user/intronIC/intronIC/intronIC_patched.py",
     "-g", "/home/user/intronIC/intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz",
     "-a", "/home/user/intronIC/intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz",
     "-n", "original_analysis",
     "-s", "-d"],
    capture_output=True,
    text=True,
    cwd="/home/user/intronIC/comparison_outputs/original"
)
if result.returncode != 0:
    print(f"ERROR: {result.stderr}")
    sys.exit(1)
print(f"✓ Original: extracted introns")

# Refactored
print("\n[2/2] Running refactored intronIC...")
result = subprocess.run(
    ["python", "-m", "intronIC_refactored.cli.main",
     "-g", "/home/user/intronIC/intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz",
     "-a", "/home/user/intronIC/intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz",
     "-n", "refactored_analysis",
     "-o", "/home/user/intronIC/comparison_outputs/refactored",
     "-s"],
    capture_output=True,
    text=True
)
if result.returncode != 0:
    print(f"ERROR: {result.stderr}")
    sys.exit(1)
print(f"✓ Refactored: extracted introns")

# Now analyze the outputs
print("\n" + "=" * 80)
print("ANALYZING INTRON DIFFERENCES")
print("=" * 80)

# Load original introns
orig_introns = {}
with open("/home/user/intronIC/comparison_outputs/original/original_analysis.bed.iic") as f:
    for line in f:
        fields = line.strip().split('\t')
        chrom, start, stop, name, score, strand = fields[:6]
        key = f"{chrom}:{start}-{stop}:{strand}"
        orig_introns[key] = name

# Load refactored introns
refact_introns = {}
with open("/home/user/intronIC/comparison_outputs/refactored/refactored_analysis.bed.iic") as f:
    for line in f:
        fields = line.strip().split('\t')
        chrom, start, stop, name, score, strand = fields[:6]
        key = f"{chrom}:{start}-{stop}:{strand}"
        refact_introns[key] = name

print(f"\nOriginal introns:    {len(orig_introns):,}")
print(f"Refactored introns:  {len(refact_introns):,}")
print(f"Difference:          {len(orig_introns) - len(refact_introns):,}")

# Find missing and extra
missing = set(orig_introns.keys()) - set(refact_introns.keys())
extra = set(refact_introns.keys()) - set(orig_introns.keys())

print(f"\nMissing in refactored: {len(missing)}")
print(f"Extra in refactored:   {len(extra)}")

# Analyze missing introns
print("\n" + "=" * 80)
print("MISSING INTRONS ANALYSIS")
print("=" * 80)

print(f"\nFirst 20 missing intron coordinates:")
for i, coord in enumerate(sorted(missing)[:20]):
    name = orig_introns[coord]
    print(f"  {i+1:2}. {coord:30} {name}")

# Check if there's a pattern by transcript
from collections import Counter
missing_transcripts = Counter()
for coord in missing:
    name = orig_introns[coord]
    # Extract transcript ID from name (format: OriAna-gene:XXX@transcript:YYY-intron_N)
    if '@transcript:' in name:
        trans = name.split('@transcript:')[1].split('-')[0]
        missing_transcripts[trans] += 1

print(f"\nTranscripts with most missing introns:")
for trans, count in missing_transcripts.most_common(10):
    print(f"  {trans}: {count} introns")

# Analyze extra introns
print("\n" + "=" * 80)
print("EXTRA INTRONS ANALYSIS")
print("=" * 80)

print(f"\nFirst 20 extra intron coordinates:")
for i, coord in enumerate(sorted(extra)[:20]):
    name = refact_introns[coord]
    print(f"  {i+1:2}. {coord:30} {name}")

# Check for coordinate shifts
print("\n" + "=" * 80)
print("COORDINATE SHIFT ANALYSIS")
print("=" * 80)
print("\nChecking if missing/extra introns are just coordinate shifts...")

# For each missing intron, check if there's a nearby extra intron (same stop, different start)
shifts = []
for missing_coord in sorted(missing)[:10]:
    chrom, coords, strand = missing_coord.split(':')
    m_start, m_stop = coords.split('-')

    # Look for extra introns with same chromosome, strand, and stop
    matches = []
    for extra_coord in extra:
        e_chrom, e_coords, e_strand = extra_coord.split(':')
        e_start, e_stop = e_coords.split('-')

        if e_chrom == chrom and e_strand == strand and e_stop == m_stop:
            shift = int(e_start) - int(m_start)
            matches.append((extra_coord, shift))

    if matches:
        shifts.append((missing_coord, matches[0]))
        print(f"\n  Missing: {missing_coord}")
        print(f"  Extra:   {matches[0][0]} (shift: {matches[0][1]:+d}bp)")

print(f"\n\nFound {len(shifts)} coordinate shifts in first 10 missing introns")

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print(f"""
Original extracted:     {len(orig_introns):,} introns
Refactored extracted:   {len(refact_introns):,} introns
Difference:             {len(orig_introns) - len(refact_introns):,} introns

Missing in refactored:  {len(missing)} introns
Extra in refactored:    {len(extra)} introns

This suggests the introns are not just "missing" but have different coordinates.
The two-pass CDS/exon merging may be creating introns with slightly different
boundaries compared to the original.

Next steps:
1. Examine specific transcripts with discrepancies
2. Check if missing introns are all exon-only (no CDS coverage)
3. Verify overlap detection logic in two-pass algorithm
""")
