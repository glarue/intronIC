#!/usr/bin/env python3
"""
Investigate why refactored filters to 7,636 introns for scoring vs original's 12,074.

Strategy:
1. Extract introns from both versions (we already have the outputs)
2. Compare filtering/omission patterns
3. Identify which filtering criteria differ
"""

import re
from collections import Counter, defaultdict
from pathlib import Path

print("="*80)
print("FILTERING DISCREPANCY INVESTIGATION")
print("="*80)

# We know:
# - Both extract 58,933 introns → 20,252 unique
# - Original: 12,074 introns included in scoring (8,178 filtered out)
# - Refactored: 7,636 introns included in scoring (12,616 filtered out)
# - Difference: 4,438 fewer introns in refactored

print("\nKnown counts:")
print(f"  Extraction: 58,933 → 20,252 unique (both versions)")
print(f"  Original scoring: 12,074 introns (8,178 filtered)")
print(f"  Refactored scoring: 7,636 introns (12,616 filtered)")
print(f"  Discrepancy: 4,438 more introns filtered in refactored")

# Strategy: Look at the refactored filtering output to see breakdown
print("\n[Step 1] Analyzing refactored filtering output from log...")

# Parse the refactored log for filtering stats
log_path = "/tmp/ref_scored.log"
if Path(log_path).exists():
    with open(log_path) as f:
        log_text = f.read()

    # Look for filtering statistics
    # Expected format: "Omitted: X short, Y ambiguous, Z non-canonical, W non-longest isoform, V overlapping"
    match = re.search(
        r'Omitted: (\d+) short, (\d+) ambiguous, (\d+) non-canonical, (\d+) non-longest isoform, (\d+) overlapping',
        log_text
    )

    if match:
        ref_short = int(match.group(1))
        ref_ambiguous = int(match.group(2))
        ref_noncanonical = int(match.group(3))
        ref_isoform = int(match.group(4))
        ref_overlap = int(match.group(5))

        print(f"\nRefactored filtering breakdown:")
        print(f"  Short: {ref_short:,}")
        print(f"  Ambiguous: {ref_ambiguous:,}")
        print(f"  Non-canonical: {ref_noncanonical:,}")
        print(f"  Non-longest isoform: {ref_isoform:,}")
        print(f"  Overlapping: {ref_overlap:,}")
        print(f"  Total filtered: {ref_short + ref_ambiguous + ref_noncanonical + ref_isoform + ref_overlap:,}")

        # Look for duplicates line
        dup_match = re.search(r'Duplicates marked: (\d+)', log_text)
        if dup_match:
            ref_duplicates = int(dup_match.group(1))
            print(f"  Duplicates: {ref_duplicates:,}")
    else:
        print("  Could not find filtering stats in refactored log")
else:
    print(f"  Log file not found: {log_path}")

# Now let's parse the original log for comparison
print("\n[Step 2] Analyzing original filtering from log...")

orig_log_path = "orig_scored.log.iic"
if Path(orig_log_path).exists():
    with open(orig_log_path) as f:
        orig_log = f.read()

    # Extract counts from original log
    # Format varies, look for key lines
    counts = {}
    for line in orig_log.split('\n'):
        if 'introns found' in line:
            match = re.search(r'\[(\d+)\].*introns found', line)
            if match:
                counts['total'] = int(match.group(1))

        if 'redundant coordinates' in line:
            match = re.search(r'\[(\d+)\].*redundant coordinates', line)
            if match:
                counts['duplicates'] = int(match.group(1))

        if 'introns included in scoring' in line:
            match = re.search(r'\[(\d+)\].*introns included', line)
            if match:
                counts['scoring'] = int(match.group(1))

        if 'introns omitted' in line:
            match = re.search(r'\[(\d+)\].*introns omitted', line)
            if match:
                counts['omitted'] = int(match.group(1))

    print(f"\nOriginal log counts:")
    for key, value in counts.items():
        print(f"  {key}: {value:,}")

    # Calculate derived values
    if 'total' in counts and 'duplicates' in counts:
        unique = counts['total'] - counts['duplicates']
        print(f"  unique (calculated): {unique:,}")

    if 'unique' in locals() and 'scoring' in counts:
        filtered_orig = unique - counts['scoring']
        print(f"  filtered (calculated): {filtered_orig:,}")

else:
    print(f"  Log file not found: {orig_log_path}")

# Key hypothesis: The bp_region check may be filtering differently
print("\n[Step 3] Hypothesis: Branch point region validation differs")
print("\nThe refactored code checks bp_region_seq for:")
print("  1. Must be at least 7bp long")
print("  2. Must have at least 7bp contiguous stretch of valid bases (no N's)")
print("\nThe original may not have this strict check, or checks it differently.")
print("\nThis could explain the ~4,400 difference if many introns have:")
print("  - Short bp regions (< 7bp due to intron being near 5' end)")
print("  - Ambiguous bases in bp region")

print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print("\nFiltering discrepancy: Refactored filters 4,438 MORE introns than original")
print("\nMost likely causes:")
print("  1. Branch point region validation (bp_region_seq length/quality)")
print("  2. Ambiguous base checking in bp region")
print("  3. Definition of 'longest isoform' (first-seen vs actual length)")
print("\nNext steps:")
print("  1. Check original's bp region validation logic")
print("  2. Compare bp_region extraction between versions")
print("  3. Verify 'longest isoform' identification matches")
print("="*80)
