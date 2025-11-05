# Gold Standard Validation - Critical Issues Found

## Issue 1: Missing 30 Introns (58,933 → 58,903)

**Problem:** Refactored version extracts 30 fewer introns than original.

**Original:** 58,933 introns found
**Refactored:** 58,903 introns found
**Difference:** -30 introns (0.05%)

**Hypothesis:** Need to investigate:
- Whether overlapping exon handling differs
- Whether exon pair generation logic differs
- Whether we're creating introns from both CDS and exon features correctly

**Status:** Requires investigation

---

## Issue 2: Incorrect Reference Intron Filtering (4 introns skipped)

**Problem:** Refactored version skips 4 U2 reference introns that original keeps.

**Original:** 387 U12 + 20,690 U2 = 21,077 total (none reported skipped)
**Refactored:** 387 U12 + 20,686 U2 = 21,073 total (4 skipped as "shorter than 55bp")

**Root Cause:**
In `cli/main.py:load_reference_sequences()` line 182:
```python
if len(intron_seq) < min_length:
    skipped_short += 1
    continue
```

This checks raw intron sequence length BEFORE extracting scoring regions.

**Original's approach:**
The original loads all introns, extracts scoring regions, THEN checks:
1. Total intron length >= min_length (default 30bp)
2. Scoring regions have sufficient valid (non-ambiguous) sequence
3. BP region has at least `bp_matrix_length` valid characters

So introns with length 50-54bp might pass if they have valid scoring regions, but we're rejecting them upfront.

**Fix:** Remove the length check from `load_reference_sequences()`. Let the omit_check logic (which runs after scoring region extraction) handle filtering.

---

## Issue 3: Missing Filtering Before Scoring (CRITICAL - 5x slowdown)

**Problem:** Refactored version scores ALL introns instead of filtering first.

**Original filtering pipeline:**
```
58,933 extracted
- 38,681 duplicates
- 8,178 omitted (66 short + 8,112 non-longest isoform)
= 12,074 sent to scoring
```

**Refactored (BROKEN):**
```
58,903 extracted
- 40 omitted (short introns only)
= 58,863 sent to scoring ❌
```

**Impact:**
- 5x more introns scored (58,863 vs 12,074)
- SVM optimization is O(n²) → **25x+ longer runtime**
- Original: 8.7 minutes
- Refactored: 30+ minutes, crashed

**Root Cause:**
Filtering logic in `cli/main.py:write_outputs()` happens AFTER scoring, but should happen BEFORE.

**Fix Required:**
1. Implement duplicate detection before scoring
2. Implement longest-isoform tagging before scoring
3. Filter out omitted introns before passing to scoring pipeline
4. Only apply `include_duplicates` flag at OUTPUT time, not scoring time

---

## Priority

1. **CRITICAL:** Fix Issue #3 (filtering before scoring) - blocks validation
2. **HIGH:** Fix Issue #2 (reference intron filtering) - affects classification accuracy
3. **MEDIUM:** Investigate Issue #1 (missing 30 introns) - affects completeness

## Next Steps

1. Fix reference intron loading (remove premature length check)
2. Implement duplicate detection and longest-isoform filtering BEFORE scoring
3. Re-run validation with fixes
4. Investigate the 30 missing introns
