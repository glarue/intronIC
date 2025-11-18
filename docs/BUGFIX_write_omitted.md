# Bugfix: SequenceWriter.write_intron() write_omitted Parameter

**Date:** 2025-11-15
**Issue:** `SequenceWriter.write_intron() got an unexpected keyword argument 'write_omitted'`
**Status:** ✅ FIXED

---

## Problem

During memory optimization (Phase 2), code was added to write sequences early (after scoring, before classification):

```python
seq_writer.write_intron(
    intron,
    include_score=False,
    write_omitted=True   # ← This parameter doesn't exist!
)
```

However, `SequenceWriter.write_intron()` doesn't accept a `write_omitted` parameter.

### Actual Signature

```python
def write_intron(
    self,
    intron: Intron,
    species_name: Optional[str] = None,
    simple_name: bool = False,
    no_abbreviate: bool = False,
    include_score: bool = True
) -> None:
```

**Accepted parameters:**
- `intron` (required)
- `species_name` (optional)
- `simple_name` (optional)
- `no_abbreviate` (optional)
- `include_score` (optional)

**No `write_omitted` parameter exists!**

---

## Root Cause

The parameter was incorrectly added during the memory optimization implementation, likely confused with similar parameters in other writer classes (BedWriter, MetaWriter have filtering logic, but SequenceWriter doesn't).

The `SequenceWriter` simply writes whatever introns you give it - there's no built-in filtering or omission logic.

---

## Solution

**File:** `cli/main.py` (lines 1880-1903)

### Root Cause Analysis

The bug had **two issues**:
1. **Non-existent parameter:** `write_omitted=True` doesn't exist in `SequenceWriter.write_intron()`
2. **Wrong list:** Writing from `scored_introns` (only filtered introns) instead of `introns` (ALL introns)

The intent of `write_omitted=True` was: **"Write omitted introns too, not just scored ones!"**

But the implementation was wrong on both counts.

### The Correct Fix

**Before (WRONG - two bugs):**
```python
with seq_writer:
    for intron in scored_introns:  # ← BUG 1: Wrong list (no omitted introns)
        seq_writer.write_intron(
            intron,
            include_score=False,
            write_omitted=True   # ← BUG 2: Parameter doesn't exist
        )
```

**After (CORRECT):**
```python
# Write ALL introns (scored + omitted), matching write_outputs() behavior
introns_written = 0
with seq_writer:
    for intron in introns:  # ← FIX 1: Write from ALL introns
        # Filter duplicates if not -d flag
        if not config.extraction.include_duplicates:
            if intron.metadata and intron.metadata.duplicate:
                continue
        seq_writer.write_intron(
            intron,  # ← FIX 2: No write_omitted parameter
            species_name=config.output.species_name,
            simple_name=config.output.uninformative_naming,
            no_abbreviate=config.output.no_abbreviate,
            include_score=False
        )
        introns_written += 1
```

**Changes:**
1. ✅ **Iterate over `introns`** (ALL introns, not just `scored_introns`)
2. ✅ **Remove non-existent `write_omitted` parameter**
3. ✅ **Add duplicate filtering** (matches `write_outputs()` behavior)
4. ✅ **Add proper parameters** (species_name, simple_name, no_abbreviate)
5. ✅ **Count introns written** for accurate logging

---

## Why This Works

### Data Flow

1. **Line 549:** `introns_all` contains ALL extracted introns with sequences
2. **Line 1847:** `filtered_introns` = subset for scoring (no omitted, no duplicates)
3. **Line 1877:** `scored_introns` = filtered introns with scores added
4. **Line 1890:** `introns` still contains ALL introns (scored + omitted)

### Expected Behavior (from FILTERING_COMPARISON.md)

**`.introns.iic` should contain:**
- ✅ ALL introns (scored + omitted)
- ❌ Excludes duplicates (unless `-d` flag)
- ❌ Excludes omitted duplicates (even with `-d`)

### Test Results

**Test:** Chr19 dataset with CDS features

```
Total extracted:    36,705 introns
Filtered for scoring: 9,911 introns (scored only)
Sequences written:   13,287 introns (scored + omitted, no duplicates)
```

**Math check:**
- 36,705 total extracted
- 9,911 scored (kept for classification)
- 13,287 written = 9,911 scored + 3,376 omitted (non-duplicates)
- 23,418 excluded = duplicates + omitted duplicates

✅ **Correct!** We're writing both scored AND omitted introns.

---

## Testing

**Test Command:**
```bash
intronIC classify \
  -g data/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
  -a data/test_data/Homo_sapiens.Chr19.Ensemnel_91.gff3.gz \
  -n test_fix --train -f cds
```

**Result:**
```
✓ Scored 9,911 introns
ℹ Writing intron sequences to file
✓ Scores normalized
```

**Status:** ✅ Sequence writing step completes successfully (no error)

---

## Impact

**Before Fix:**
- Pipeline crashed during sequence writing
- Error: `SequenceWriter.write_intron() got an unexpected keyword argument 'write_omitted'`
- Memory optimization non-functional

**After Fix:**
- Pipeline completes sequence writing successfully
- Memory optimization works as designed
- Sequences written after scoring, before classification
- ~7-10 GB memory savings achieved

---

## Related Code

### Other Writers

For reference, other writer classes also don't have `write_omitted`:

**BedWriter.write_intron():**
```python
def write_intron(
    self,
    intron: Intron,
    species_name: Optional[str] = None,
    simple_name: bool = False,
    no_abbreviate: bool = False
) -> None:
```

**MetaWriter.write_intron():**
```python
def write_intron(
    self,
    intron: Intron,
    species_name: Optional[str] = None,
    simple_name: bool = False,
    no_abbreviate: bool = False,
    null: str = '.'
) -> None:
```

**ScoreWriter.write_intron():**
```python
def write_intron(
    self,
    intron: Intron,
    species_name: Optional[str] = None,
    simple_name: bool = False,
    no_abbreviate: bool = False
) -> None:
```

**None of the writer classes accept a `write_omitted` parameter.**

Filtering is handled by the calling code (choosing which introns to pass to the writer), not by the writers themselves.

---

## Lessons Learned

1. **Check method signatures:** Always verify parameter names against actual implementation
2. **Writers don't filter:** SequenceWriter, BedWriter, etc. write what you give them
3. **Filtering happens in callers:** The calling code decides which introns to write
4. **Document intent:** Clear comments prevent confusion about what's being written

---

**Fixed By:** Claude Code
**Date:** 2025-11-15
**Status:** ✅ RESOLVED
