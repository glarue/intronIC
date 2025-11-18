# Memory Analysis: Original intronIC vs Refactored

**Date:** 2025-11-15
**Goal:** Understand how original achieved ~5 GB peak vs our ~14 GB

---

## Original intronIC Memory Strategy

### Key Insight: **NEVER materialize full list with sequences**

**Code flow (lines 4668-4847):**

```python
def filter_introns_write_files(...):
    final_introns = []  # Will contain introns WITHOUT full sequences

    # Open output files
    seq_file = open('seqs.iic', 'w')
    meta_file = open('meta.iic', 'w')
    bed_file = open('bed.iic', 'w')

    # GENERATOR: yields one intron at a time (NOT a list!)
    for intron in get_sub_seqs(all_introns, genome, ...):  # Line 4737

        # 1. Check if should be omitted
        intron.omit_check(...)  # Line 4748

        # 2. Tag duplicates/isoforms
        intron = add_tags(intron, ...)  # Line 4772

        # 3. Only add to final_introns if NOT omitted AND NOT duplicate
        if not intron.omitted and intron.duplicate is False:
            final_introns.append(intron)  # Line 4793

        # 4. Write to files IMMEDIATELY (while seq still in memory)
        if not NO_SEQS:
            seq_string = output_format(intron, 'SEQ', ...)
            seq_file.write(seq_string + '\n')  # Line 4844-4845

        # 5. CLEAR SEQUENCE IMMEDIATELY (!!!)
        intron.seq = None  # Line 4847 ← KEY LINE!

        # Loop continues - intron goes out of scope, garbage collected
        # Next intron loaded with sequences, processed, written, cleared...

    # Close files
    seq_file.close()

    # Return final_introns - these have seq=None!
    return final_introns
```

---

## Memory Comparison

### Original intronIC

| Stage | In Memory | Size | Details |
|-------|-----------|------|---------|
| **Iteration Start** | Current intron (with seq) | ~10 KB | One at a time from generator |
| | `final_introns` (no seq) | ~100 MB | Only non-omitted, non-dupes, seq=None |
| **After Writing** | Current intron cleared | ~100 bytes | seq=None |
| | `final_introns` growing | ~100 MB → ~1 GB | Accumulating cleared introns |
| **Peak During Loop** | ~1-2 GB | | One chromosome in genome + accumulating list |
| **After Loop Complete** | `final_introns` (no seq) | ~1-2 GB | ~200k introns × 5 KB each (no full seq) |
| **During Scoring** | Same introns + scores | ~2-3 GB | Only motif seqs needed |
| **Peak Total** | **~5 GB** ✅ | | Matches reported usage |

### Refactored intronIC (Current)

| Stage | In Memory | Size | Details |
|-------|-----------|------|---------|
| **Before Extraction** | `introns_list` (no seq) | ~3 GB | All introns, coordinates only |
| **During Extraction** | `introns_list` (no seq) | ~3 GB | Old list still in memory |
| | `introns_all` (growing) | 0 → ~13 GB | NEW list with sequences |
| **After Extraction** | `introns_all` (with seq) | ~13 GB | ALL introns, 1M × 13 KB |
| | (OLD: introns_list deleted) | ~0 GB | After our fix |
| **During Scoring** | `scored_introns` (with seq) | ~13 GB | Still have full sequences |
| **Peak Total** | **~13-14 GB** ❌ | | 3× higher than original |

---

## Key Differences

### 1. **Generator vs Materialized List**

**Original:**
```python
for intron in get_sub_seqs(...):  # Generator - one at a time
    process(intron)
    write(intron)
    intron.seq = None  # Clear immediately
    # Intron goes out of scope, GC reclaims memory
```

**Memory during loop:** ~100 MB (current intron) + ~1 GB (accumulating final_introns without seq)

**Refactored:**
```python
introns_with_seq = sequence_extractor.extract_sequences(...)  # Generator
introns_all = list(introns_with_seq)  # ← MATERIALIZES ENTIRE LIST!
# Now ALL introns with sequences are in memory
```

**Memory after materialization:** ~13 GB (1M introns × 13 KB with sequences)

---

### 2. **Immediate vs Deferred Sequence Clearing**

**Original:**
```python
for intron in get_sub_seqs(...):
    write_sequence_to_file(intron)  # Write while seq in memory
    intron.seq = None  # Clear IMMEDIATELY
    final_introns.append(intron)  # Append intron WITHOUT seq
```

**Result:** `final_introns` contains introns with `seq = None` → small memory footprint

**Refactored:**
```python
introns_all = list(extract_sequences(...))  # ALL introns with sequences
# ... process ...
scored_introns = score_introns(introns_all)  # Still have sequences
# ... later ...
write_outputs(scored_introns, ...)  # THEN write sequences
# Only cleared sequences in dual-track approach for classification
```

**Result:** Keep full sequences until output writing → large memory footprint

---

### 3. **Filtering Before vs After Sequence Extraction**

**Original:** Partial filtering during extraction

```python
for intron in get_sub_seqs(...):
    intron.omit_check(...)  # Mark as omitted

    if not intron.omitted and not intron.duplicate:
        final_introns.append(intron)  # Only keep non-omitted

    # But still writes ALL introns to seq file (including omitted)
    write_sequence(intron)
    intron.seq = None
```

**Result:** `final_introns` contains ~200k introns (1M - 800k omitted/dupes)

**Refactored:**
```python
introns_all = extract_sequences(all_introns)  # Extract for ALL 1M introns
filtered_introns = filter(introns_all)  # Filter AFTER extraction
```

**Result:** Extract and keep sequences for all 1M introns, then filter

---

## Memory Breakdown: Where Our 14 GB Goes

For **1 million introns** with full sequences:

| Component | Per Intron | Total (1M introns) |
|-----------|------------|-------------------|
| **Coordinates** | ~200 bytes | 200 MB |
| **Metadata** | ~500 bytes | 500 MB |
| **Full intron sequence** | ~500 bytes | **500 MB** |
| **Upstream flank** | ~200 bytes | **200 MB** |
| **Downstream flank** | ~200 bytes | **200 MB** |
| **BP region sequence** | ~50 bytes | **50 MB** |
| **Motif sequences** | ~30 bytes | 30 MB |
| **Scores** | ~100 bytes | 100 MB |
| **Python object overhead** | ~12 KB | **12 GB** ← Biggest! |
| **Total per intron** | ~13-14 KB | **13-14 GB** |

**Python object overhead is massive:**
- Dict for `__dict__` (even with slots)
- Reference counting
- Type pointers
- Garbage collection metadata

**Original intronIC objects WITHOUT sequences:**

| Component | Per Intron | Total (200k kept) |
|-----------|------------|------------------|
| Coordinates | ~200 bytes | 40 MB |
| Metadata | ~500 bytes | 100 MB |
| **seq = None** | 8 bytes | 1.6 MB |
| Motif sequences | ~30 bytes | 6 MB |
| Scores | ~100 bytes | 20 MB |
| Python overhead | ~4 KB | **800 MB** |
| **Total per intron** | ~5 KB | **~1 GB** |

Plus ~1 GB for omitted introns (written to file but not kept in final_introns)

**Total:** ~2-3 GB for final_introns, ~5 GB peak during operations

---

## Why Our "Dual-Track" Approach Still Uses More Memory

**Our current approach:**
1. Extract sequences for all 1M introns: **13 GB**
2. Score all: **13 GB** (same introns)
3. Create cleared copy for classification: **+100 MB**
4. Keep originals for output: **Still 13 GB**
5. Write outputs: **13 GB**

**We're keeping full sequences until output writing!**

---

## Solution Options

### Option A: **Immediate Write & Clear** (Matches Original)

**Approach:** Write sequences to file during extraction, clear immediately

```python
# In extract_introns_from_annotation():
seq_file = open(f'{species}.introns.iic', 'w')

for intron in sequence_extractor.extract_sequences(introns_list, ...):
    # Write sequence to file IMMEDIATELY
    write_sequence(seq_file, intron)

    # Clear large sequences, keep only motif sequences
    intron = clear_large_sequences_keep_motifs(intron)

    # Now intron is small (~5 KB instead of 13 KB)
    introns_all.append(intron)

seq_file.close()

# introns_all now contains 1M introns WITHOUT full sequences
# Memory: 1M × 5 KB = 5 GB (instead of 13 GB)
```

**Memory savings:** 13 GB → 5 GB ✅

**Pros:**
- Matches original memory profile
- Works for any genome size
- Sequences written once, never kept in memory

**Cons:**
- Must re-read sequences from file if needed later (not needed in practice)
- Slightly more complex flow (write during extraction)

---

### Option B: **Clear After Filtering** (Partial Fix)

**Approach:** Filter first (identify omitted/dupes), only extract sequences for kept introns

```python
# 1. Create introns without sequences
introns_no_seq = extract_intron_coordinates(...)  # 3 GB

# 2. Filter BEFORE sequence extraction
filtered_introns = filter(introns_no_seq)  # Remove omitted/dupes
# Result: ~200k introns to extract (instead of 1M)

# 3. Extract sequences only for filtered introns
introns_with_seq = extract_sequences(filtered_introns, ...)
# Memory: 200k × 13 KB = 2.6 GB (instead of 13 GB)
```

**Memory savings:** 13 GB → 3 GB ✅

**Pros:**
- Don't waste time/memory extracting sequences for introns we'll discard
- Simpler than Option A

**Cons:**
- Can't write sequences for omitted introns (output file difference)
- May break if user wants sequences for all introns

---

### Option C: **Generator-Based Pipeline** (Most Like Original)

**Approach:** Never materialize full list, process one chromosome at a time

```python
def process_genome_streaming(genome, annotation, ...):
    seq_file = open('seqs.iic', 'w')
    final_introns = []

    # Process one chromosome at a time
    for chrom_name, chrom_seq in stream_genome(genome):
        chrom_introns = get_introns_for_chromosome(annotation, chrom_name)

        for intron in extract_sequences_for_chromosome(chrom_introns, chrom_seq):
            # Filter
            if should_omit(intron):
                intron.mark_omitted()

            # Write to file
            write_sequence(seq_file, intron)

            # Clear sequence
            intron = clear_large_sequences(intron)

            # Keep only non-omitted
            if not intron.omitted:
                final_introns.append(intron)

        # Chromosome done, memory freed automatically
        del chrom_introns, chrom_seq

    seq_file.close()
    return final_introns  # ~5 GB
```

**Memory savings:** 13 GB → 5 GB ✅

**Pros:**
- Matches original exactly
- Constant memory regardless of genome size
- Most memory-efficient

**Cons:**
- Requires significant refactoring
- Harder to test individual components
- Generator-based flow is complex

---

## Recommendation

**For immediate fix:** Implement **Option A** (Immediate Write & Clear)

**Why:**
1. Minimal code changes (modify extraction functions only)
2. Gets us to original's ~5 GB memory usage
3. Maintains all functionality
4. Write sequences during extraction (when they're in memory anyway)
5. Clear immediately after writing
6. Keep small introns (~5 KB) for scoring

**Implementation:**
1. Modify `extract_introns_from_annotation()` to open sequence file
2. Write each intron to file as it's extracted
3. Clear large sequences immediately after writing
4. Keep motif sequences for scoring
5. Result: `introns` list contains small introns (seq=None)

**Expected memory:**
- During extraction: 3-5 GB (streaming genome + accumulating small introns)
- During scoring: 4-6 GB (small introns + scores)
- During classification: 5-7 GB (dual-track approach)
- **Peak: ~7 GB** (vs current 14 GB, vs original 5 GB)

Close to original but with better code organization!

---

## Next Steps

1. Implement Option A (write sequences during extraction)
2. Verify memory usage drops to ~7 GB
3. If still too high, investigate Python object overhead (consider numpy arrays for scores)
4. If needed, implement Option C (full streaming pipeline)

The key insight: **Don't materialize the full list with sequences!**
