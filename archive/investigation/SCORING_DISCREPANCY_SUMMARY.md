# Scoring Discrepancy Investigation Summary

**Date**: 2025-11-06
**Session**: Continued from context limit
**Status**: Partially resolved - root cause identified, fix in progress

---

## Problem Statement

The refactored version scores significantly fewer introns than the original:
- **Original**: 12,074 introns kept for scoring
- **Refactored**: 7,636 introns kept for scoring (before fix)
- **After fix attempt**: 5,704 introns kept for scoring
- **Discrepancy**: ~6,400 fewer introns (getting worse!)

---

## Root Cause Identified

### Primary Issue: Longest Isoform Identification

The original intronIC sorts introns by **hierarchical_sort_attrs** which includes:
1. `defined_by` (CDS before exon)
2. **`parent_length * -1`** (DESCENDING by transcript length!) ← KEY!
3. `parent` (transcript ID)
4. `family_size * -1` (descending)
5. `index` (intron index)
6. `line_number` (annotation order)

This means the "first seen wins" logic actually becomes **"longest transcript wins"** because introns from longer transcripts are processed first!

### Original Code Reference
```python
# intronIC/intronIC_patched.py, line 2633
def hierarchical_sort_attrs(intron):
    sort_features = [
        intron.defined_by,
        intron.parent_length * -1,  # NEGATIVE = descending sort!
        intron.parent,
        intron.family_size * -1,
        intron.index,
        intron.line_number]
    return tuple(sort_features)

# Line 2674 - introns sorted before processing
for intron in sorted(
    introns_by_region[region_name],
    key=hierarchical_sort_attrs):
```

### Refactored Issue

The refactored version used simple "first seen in annotation file" order, which produces completely different longest isoform assignments!

---

## Fix Applied

Updated `/home/user/intronIC/intronIC_refactored/extraction/filters.py` line 130-165:

```python
def _identify_longest_isoforms(self, introns: List[Intron]) -> None:
    """
    First pass: identify "longest" transcript per gene.

    Matches original's hierarchical_sort_attrs sorting.
    """
    # Sort introns by parent_length descending
    sorted_introns = sorted(
        introns,
        key=lambda i: (
            i.metadata.defined_by or '',         # CDS before exon
            -(i.metadata.parent_length or 0),   # Descending by parent length (KEY!)
            i.metadata.parent or '',             # Transcript ID
            -(i.metadata.family_size or 0),     # Descending by family size
            i.metadata.index or 0                # Intron index
        )
    )

    # First transcript seen (after sorting) is the longest
    for intron in sorted_introns:
        grandparent = intron.metadata.grandparent
        parent = intron.metadata.parent

        if grandparent and grandparent not in self.longest_isoforms:
            self.longest_isoforms[grandparent] = parent
```

---

## Testing Results

### Test 1: Logic Verification
Created `debug_longest_isoform.py` to test sorting logic:
- ✅ Sorting works correctly
- ✅ Longest transcripts identified correctly
- ✅ Logic matches original behavior

### Test 2: Full Pipeline
Ran refactored with fix:
- ❌ **Result**: 5,704 introns kept (WORSE than before!)
- **Before fix**: 7,636 introns
- **After fix**: 5,704 introns
- **Expected**: 12,074 introns

---

## Current Status

The fix is **partially correct** but produces unexpected results:
1. ✅ Sorting logic is correct (verified in isolation)
2. ✅ Longest isoform identification logic is correct
3. ❌ Full pipeline produces wrong count (5,704 vs 12,074)
4. ⚠️ Something else is filtering out additional introns

---

## Hypothesis for Remaining Issue

Possible causes for the 5,704 vs 12,074 discrepancy:

1. **Order of operations**: The refactored may be applying filters in a different order than the original
2. **Omission checking**: The `_check_omission` method may be too aggressive
3. **BP region validation**: The branch point region check may be filtering too many introns
4. **Duplicate handling**: Duplicates may be interacting with longest_isoform filtering incorrectly

---

## Detailed Counts

### Original (from logs)
```
58,933 total introns extracted
38,681 duplicates
20,252 unique introns
8,178 omitted (8,112 for [i] not longest isoform, 66 for [s] short)
12,074 kept for scoring
```

### Refactored (before fix)
```
58,863 total introns extracted
38,677 duplicates
20,186 unique introns
50,893 marked as "non-longest isoform" (includes duplicates!)
7,636 kept for scoring
```

### Refactored (after fix)
```
58,863 total introns extracted
38,677 duplicates
20,186 unique introns
46,365 marked as "non-longest isoform" (includes duplicates!)
5,704 kept for scoring  ← WORSE!
```

---

## Next Steps

1. **Debug the filtering process**:
   - Add logging to show which introns are being filtered and why
   - Compare the `longest_isoforms` dictionary between original and refactored
   - Check if the sorted order is being maintained throughout filtering

2. **Verify all attributes**:
   - Ensure `parent_length` is populated correctly for all introns
   - Check if any introns have `None` values that affect sorting

3. **Compare specific transcripts**:
   - Pick a specific gene with multiple isoforms
   - Trace through both versions to see which transcript is marked as "longest"
   - Verify that the same transcript is chosen in both versions

4. **Check omission logic**:
   - The bp_region validation might be interacting poorly with the new sorting
   - Verify that omission checking happens in the same order as the original

---

## Files Modified

- `/home/user/intronIC/intronIC_refactored/extraction/filters.py` (lines 130-165)
  - Updated `_identify_longest_isoforms` to sort by parent_length descending

## Investigation Scripts Created

- `/home/user/intronIC/investigate_filtering_discrepancy.py` - Analyzes filtering patterns
- `/home/user/intronIC/debug_longest_isoform.py` - Tests sorting logic in isolation

---

## Key Learnings

1. **"Longest isoform" != "first seen"**: The original uses length-based sorting, not annotation order
2. **Sort order matters**: Different processing orders produce different "longest" assignments
3. **Context in code comments can be misleading**: Code said "first seen wins" but meant "first in SORTED order wins"

---

## Recommendation

The issue is subtle and requires deeper investigation. Suggested approach:

1. Add detailed logging to both versions showing:
   - Which transcript is selected as "longest" for each gene
   - Why specific introns are being omitted
   - The order introns are processed in

2. Create a minimal test case with just one gene that has multiple isoforms

3. Step through both versions side-by-side to find where they diverge

This will likely require another session to fully resolve.
