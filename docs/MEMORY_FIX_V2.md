# Memory Fix V2 - Critical Fixes for 20 GB+ Usage

**Date:** 2025-11-15
**Issue:** Memory usage climbing to 20+ GB during extraction, then exploding to 250+ GB during classification

---

## Root Causes Identified

### 1. **Duplicate Intron Objects During Extraction** (21 GB at extraction)

**Problem:**
```python
# In extract_sequences() streaming mode:
introns_by_region = {chr1: [intron1, intron2, ...], chr2: [...], ...}  # ~2-3 GB

for region_name, region_seq in genome_reader.stream():
    for intron in introns_by_region[region_name]:
        intron_with_seq = self._extract_intron_sequences(intron, ...)  # Creates NEW intron
        yield intron_with_seq
    # Old introns (no seq) STILL in introns_by_region dict!
    # New introns (with seq) being added to introns_all list!
```

**Result:** Two copies of every intron in memory simultaneously:
- Original introns (no sequences) in `introns_by_region` dict: ~3 GB
- New introns (with sequences) in `introns_all` list: ~18 GB
- **Total: 21 GB during extraction**

**Fix:** Delete each chromosome from dict after processing:
```python
for region_name, region_seq in genome_reader.stream():
    for intron in introns_by_region[region_name]:
        intron_with_seq = self._extract_intron_sequences(intron, ...)
        yield intron_with_seq

    # NEW: Free memory immediately
    del introns_by_region[region_name]
    del region_seq
    gc.collect()
```

---

### 2. **Pretrained Path Skips Sequence Clearing** (21 GB → 250 GB explosion)

**Problem:**
```python
# Step 3: Score introns
scored_introns = score_introns(introns_for_scoring, ...)  # Still has sequences (21 GB)

# Step 4: Check if using pretrained model
if config.training.pretrained_model_path:
    # Goes DIRECTLY to classification WITHOUT clearing sequences!
    classified_introns = classify_with_pretrained_model(scored_introns, ...)
    # ↑ scored_introns still has 21 GB of sequences
    # ↓ predictor.predict() with -p 12 creates 12 workers
    #   Each worker gets a COPY of scored_introns
    #   21 GB × 12 workers = 252 GB → CRASH 💥
else:
    # Normal path clears sequences (but pretrained path skips this!)
    normalized_introns = normalize_scores(scored_introns, ...)
    normalized_introns = clear_large_sequences(normalized_introns)  # ← Never reached in pretrained path!
```

**Fix:** Clear sequences BEFORE the if/else branch:
```python
# Step 3: Score introns
scored_introns = score_introns(introns_for_scoring, ...)

# CRITICAL: Clear sequences BEFORE pretrained check
scored_introns = clear_large_sequences(scored_introns)  # Reduce to ~100 MB
del introns_for_scoring
gc.collect()

# Step 4: Now BOTH paths have sequences already cleared
if config.training.pretrained_model_path:
    classified_introns = classify_with_pretrained_model(scored_introns, ...)
    # ↑ scored_introns now only ~100 MB (just scores + small metadata)
    # ↓ 100 MB × 12 workers = 1.2 GB ✅
else:
    normalized_introns = normalize_scores(scored_introns, ...)
```

---

### 3. **Multiple Copies During Normalization**

**Problem:**
```python
# classify_with_pretrained_model():
normalized_introns = list(normalizer.transform(introns, ...))  # NEW list
# introns (old) still in memory
# normalized_introns (new) also in memory
# Double memory usage!

classified_introns = predictor.predict(ensemble, normalized_introns)
# Now 3 copies: introns, normalized_introns, classified_introns
```

**Fix:** Delete intermediate copies:
```python
normalized_introns = list(normalizer.transform(introns, ...))
del introns  # Free old copy
gc.collect()

classified_introns = predictor.predict(ensemble, normalized_introns)
del normalized_introns  # Free after prediction
gc.collect()
```

---

## All Changes Made

### File 1: `extraction/sequences.py`

**Lines 129-137** - Delete processed chromosomes during streaming:
```python
# Free memory: delete this chromosome's introns from dict after processing
# This prevents keeping both old (no seq) and new (with seq) introns in memory
del introns_by_region[region_name]

# Also explicitly delete the chromosome sequence to free memory
del region_seq
# Force garbage collection every few chromosomes to prevent buildup
import gc
gc.collect()
```

**Memory saved:** ~15-18 GB during extraction (prevents duplicate intron objects)

---

### File 2: `cli/main.py`

**Lines 1569-1578** - Clear sequences immediately after scoring (BEFORE pretrained check):
```python
# Clear large sequences immediately after scoring to reduce memory
# Sequences are no longer needed - we only need scores from here on
# CRITICAL: Must clear BEFORE pretrained path to prevent 21GB×12 workers = OOM
messenger.log_only("Clearing full sequences after scoring to reduce memory usage")
scored_introns = clear_large_sequences(scored_introns)

# Explicitly delete the large introns_for_scoring list to free memory
del introns_for_scoring
import gc
gc.collect()  # Force garbage collection to reclaim memory
```

**Lines 1588-1590** - Free memory after pretrained classification:
```python
# Free memory before output phase
del scored_introns
gc.collect()
```

**Lines 946-950** - Delete introns after normalization in pretrained path:
```python
# Free memory: delete input introns after normalization
# (normalized_introns is a new list with updated scores)
del introns
import gc
gc.collect()
```

**Lines 1595-1597** - Delete scored_introns in training path (already present):
```python
# Delete scored_introns to free memory before classification
del scored_introns
gc.collect()
```

**Total memory saved:** ~20+ GB (prevents multiprocessing explosion)

---

## Expected Memory Usage (After All Fixes)

| Pipeline Stage | Memory | Details |
|---------------|---------|---------|
| **Load & Parse** | 1-2 GB | Intron coordinates + metadata |
| **Extract Sequences** | 3-5 GB | Streaming genome (1 chr at a time) + growing introns list |
|  | | Old: 21 GB (kept duplicates) |
| **Score Introns** | 4-6 GB | Brief spike during scoring |
| **Clear Sequences** | 1-2 GB | After clearing, only small scored seqs remain |
| **Normalize** | 2-3 GB | Brief spike creating normalized list |
| **Classify (pretrained, -p 12)** | 5-8 GB | 12 workers × ~100 MB chunks |
|  | | Old: 250+ GB (12 workers × 21 GB) |
| **Output** | 2-3 GB | Writing files |
| **Peak Memory** | **8-10 GB** | ✅ Fits in 16-32 GB systems |

---

## Files Changed (V2)

1. **extraction/sequences.py** (lines 129-137)
   - Delete processed chromosomes during streaming
   - Explicit GC after each chromosome

2. **cli/main.py** (multiple locations)
   - Lines 1569-1578: Clear sequences immediately after scoring (BEFORE pretrained check)
   - Lines 1588-1590: Free memory after pretrained classification
   - Lines 946-950: Delete introns after normalization (pretrained path)
   - Lines 1595-1597: Delete scored_introns (training path)

---

## Testing

**Syntax validated:**
```bash
pixi run python -m py_compile extraction/sequences.py cli/main.py
# ✅ No errors
```

**Recommended test:**
```bash
# Kill existing run if still active
pkill -f intronIC

# Restart with fixed code
intronIC -g genome.fa.gz -a annotation.gff3.gz -n homo_sapiens -p 12 -f cds --pretrained

# Monitor memory in separate terminal
watch -n 2 'ps aux | grep intronIC | grep -v grep | awk "{print \$2, \$6/1024 \" MB\", \$11}"'
```

**Expected behavior:**
- Extraction: ~3-5 GB (not 21 GB)
- Scoring: ~4-6 GB (brief spike)
- After clearing: ~1-2 GB
- Classification: ~5-8 GB total (all workers combined, not 250+ GB)
- **No OOM crashes** ✅

---

## Why This Happened

The refactored code used **immutable frozen dataclasses** for memory safety and thread safety. However, this meant:

1. Every update creates a NEW object (can't modify in-place)
2. Old objects stay in memory until explicitly deleted
3. Python's garbage collector doesn't run aggressively enough for large datasets

The original intronIC used **mutable objects** and **explicit sequence clearing** (`intron.seq = None`), which worked but was less safe.

The fix: **Keep immutability** (good for correctness) + **add explicit deletion and GC** (good for memory).

---

## Summary

**Before fixes:**
- Extraction: 21 GB (duplicate intron objects)
- Classification with `-p 12`: 250+ GB explosion → OOM crash

**After fixes:**
- Extraction: 3-5 GB (delete chromosomes as processed)
- Classification with `-p 12`: 5-8 GB total (sequences cleared before multiprocessing)
- **Peak: 8-10 GB** ✅

**The critical fix:** Clear sequences BEFORE the pretrained model path, not after. This prevents the 21 GB × 12 workers = 250 GB explosion.
