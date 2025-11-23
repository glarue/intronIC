# Streaming Extraction - Final Implementation Plan

**Date:** 2025-11-20
**Current Status:** 28 GB peak during extraction
**Target:** ~3-5 GB peak
**Blocker Identified:** Scorer needs full `seq` field

---

## Problem Analysis

### Current Memory Usage (2.1M introns)

| Component | Size per Intron | Total Memory |
|-----------|-----------------|--------------|
| Full sequence (`seq`) | ~500 bytes | 1.05 GB |
| Flanks (upstream/downstream) | ~400 bytes | 840 MB |
| Display sequences | ~50 bytes | 105 MB |
| Scoring sequences (five_seq, three_seq, bp_seq) | ~40 bytes | 84 MB |
| BP region | ~50 bytes | 105 MB |
| Coordinates, scores, metadata | ~300 bytes | 630 MB |
| Python object overhead | ~12 KB | **25.2 GB** |
| **TOTAL** | **~13 KB** | **~27 GB** |

**Current Flow:**
```
Extract 2.1M introns → List with all sequences (27 GB)
↓
Filter to 200k for scoring (still 27 GB in memory)
↓
Score 200k (still 27 GB in memory)
↓
Write sequences to disk
↓
Clear seq field → 2.1M introns × 600 bytes = 1.3 GB
↓
Classify
```

**Peak: 28 GB** (during extraction through scoring)

---

## Why We Can't Clear Sequences Before Scoring

### Scorer Requirements (`scoring/scorer.py:124-486`)

```python
# Line 124: Validation
if intron.sequences is None or intron.sequences.seq is None:
    raise ValueError("Intron has no sequence. Cannot score...")

# Lines 441, 446, 486: Region extraction
def _extract_five_region(self, intron):
    return intron.sequences.seq[start:stop]  # Extracts from full seq

def _extract_three_region(self, intron):
    return intron.sequences.seq[intron_start:intron_stop]  # Extracts from full seq
```

**The scorer extracts scoring regions ON-THE-FLY from `seq`**, even though we store `five_seq`, `three_seq` during extraction.

### Why Scorer Re-extracts Instead of Using Stored Sequences

**Flexibility:** Different scoring configs might need different region sizes:
- `config.scoring.scoring_regions`: Variable window sizes
- Non-canonical introns: Different extraction logic
- U12 correction: Modified coordinates → re-extract needed

**Conclusion:** **Scorer MUST have access to full `seq` field.**

---

## Solution: Two Approaches

### Approach A: Modify Scorer to Use Pre-Extracted Sequences (High Risk)

**Changes Required:**
1. Ensure extraction always extracts exact regions scorer needs
2. Modify scorer to use `intron.sequences.five_seq` instead of extracting
3. Handle all edge cases (short introns, non-canonical, U12 correction)
4. Update scorer validation to not require `seq` field

**Pros:**
- Enables clearing seq before scoring
- 27 GB → 3 GB during scoring

**Cons:**
- High risk: Changes core scoring logic
- Fragile: Extraction and scoring tightly coupled
- Hard to validate: Need to prove outputs identical
- Breaks if scoring configs change

**Estimated Effort:** 4-6 hours + extensive testing
**Risk:** HIGH ⚠️

---

### Approach B: Stream-Write-Clear During Extraction (Recommended)

**Keep scorer unchanged, optimize extraction phase only.**

**New Flow:**
```
For each intron (one at a time):
    1. Extract sequences (has full seq) - 13 KB
    2. Write seq to .introns.iic immediately
    3. Clear seq field (keep scoring seqs) - now 600 bytes
    4. Append cleared intron to list

Memory: One full intron (13 KB) + accumulated cleared (growing to 1.3 GB)

↓ At end of extraction
2.1M introns with seq=None but five_seq, three_seq, bp_seq intact (1.3 GB)

↓ Scoring
Wait - THIS WON'T WORK! Scorer needs seq!
```

**Blocker:** Scorer still needs `seq` field after extraction.

---

### Approach C: Defer Sequence Extraction Until After Filtering (Best Balance)

**Key Insight:** We extract sequences for 2.1M introns, but only ~200k get scored!

**Filter First, Then Extract:**

```
1. Parse annotations → 2.1M introns (coordinates only, ~1 GB)

2. Filter on metadata (no sequences needed):
   - Duplicates (coordinates only)
   - Longest isoform (metadata only)
   - Length checks (coordinates only)
   - Overlaps (coordinates only)
   → Reduces to ~300k introns to extract

3. Extract sequences for 300k introns:
   Stream extraction:
   - Extract one at a time
   - Write to .introns.iic immediately
   - Keep with sequences for scoring
   Memory: 300k × 13 KB = 3.9 GB ✅

4. Score 200k introns (still 3.9 GB)

5. Clear seq after scoring → 1 GB

6. For omitted introns (1.8M):
   - Need sequences for output
   - Extract on-demand or mark as "not extracted"
   - Original intronIC wrote "sequence unavailable" for omitted introns
```

**Memory Profile:**
- Parsing: 1 GB (coordinates only)
- After pre-filtering: 1 GB (still coordinates only)
- Extraction: 3.9 GB (300k with sequences)
- Scoring: 3.9 GB
- Post-scoring: 1 GB (cleared)

**Peak: ~4 GB** ✅ (vs. 28 GB current)

**Pros:**
- 85% memory reduction
- No scorer changes needed
- Lower risk
- Matches original intronIC philosophy

**Cons:**
- Omitted introns won't have sequences in output (ACCEPTABLE - they're omitted!)
- Requires refactoring filter logic to work without sequences

---

## Approach C Detailed Implementation Plan

### Phase 1: Separate Metadata-Based and Sequence-Based Filtering

**Current Filter (`extraction/filters.py`):**
```python
def filter_introns(self, introns: List[Intron]) -> List[Intron]:
    # Operates on introns WITH sequences
    # Checks: length, duplicates, isoforms, ambiguous bases, etc.
```

**New Filter Split:**
```python
def prefilter_on_metadata(self, introns: List[Intron]) -> List[Intron]:
    """Filter using only coordinates and metadata (no sequences needed)."""
    for intron in introns:
        # Length check
        if intron.length < self.min_length: omit
        if intron.length > self.max_length: omit

        # Duplicates (coordinate comparison)
        if intron.coordinates == previous.coordinates: mark_duplicate

        # Longest isoform (metadata)
        if not intron.metadata.longest_isoform: omit

        # Overlaps (coordinate comparison)
        if overlaps_with_other(intron): omit

    return passing_introns  # ~300k introns

def filter_with_sequences(self, introns: List[Intron]) -> List[Intron]:
    """Filter using sequence data (after extraction)."""
    for intron in introns:
        # Ambiguous bases in scoring regions
        if has_ambiguous_bases(intron.sequences.five_seq): omit
        if has_ambiguous_bases(intron.sequences.three_seq): omit

        # Non-canonical (if excluding)
        if not self.allow_noncanonical:
            if intron.sequences.terminal_dinucleotides not in CANONICAL: omit

    return scored_introns  # ~200k introns
```

**Files to modify:**
- `extraction/filters.py`: Split filter logic
- `cli/main.py`: Call prefilter before extraction

### Phase 2: Refactor Extraction to Accept Filtered Introns

**Current:**
```python
# Extract ALL annotations (2.1M)
introns_list = intron_generator.generate_introns_from_hierarchy(...)

# Extract sequences for ALL
introns_with_seq = sequence_extractor.extract_sequences(introns_list, ...)
introns_all = list(introns_with_seq)  # ← 28 GB HERE!

# Then filter
filtered = intron_filter.filter_introns(introns_all)
```

**New:**
```python
# Extract ALL annotations (2.1M, coordinates only)
introns_list = intron_generator.generate_introns_from_hierarchy(...)

# PRE-FILTER on metadata (NO sequences needed)
prefiltered = intron_filter.prefilter_on_metadata(introns_list)  # ~300k

# Extract sequences ONLY for prefiltered introns
introns_with_seq = sequence_extractor.extract_sequences(prefiltered, ...)

# Stream write-and-keep
with SequenceWriter(output_path) as writer:
    introns_all = []
    for intron in introns_with_seq:  # One at a time
        writer.write_intron(intron)  # Write immediately
        introns_all.append(intron)  # Keep with sequences for scoring

# Filter on sequences
filtered = intron_filter.filter_with_sequences(introns_all)  # ~200k
```

**Memory:**
- introns_list: 2.1M × 1 KB = 2.1 GB (coordinates + metadata)
- prefiltered: 300k × 1 KB = 300 MB
- introns_all: 300k × 13 KB = 3.9 GB
- After scoring clear: 300k × 600 bytes = 180 MB

**Peak: ~4 GB** ✅

### Phase 3: Streaming Within Extraction

**Even Better:** Don't accumulate introns_all at all during extraction:

```python
with SequenceWriter(output_path) as writer:
    for intron in introns_with_seq:
        writer.write_intron(intron)
        # Don't accumulate yet - just write and continue

# THEN re-extract for scoring (or keep in memory if acceptable)
```

But this requires re-reading from disk, which is slower. Keeping 3.9 GB in memory is acceptable.

---

## Implementation Checklist

### Step 1: Split Filtering Logic
- [ ] Create `prefilter_on_metadata()` in `extraction/filters.py`
- [ ] Create `filter_with_sequences()` in `extraction/filters.py`
- [ ] Update filter tests
- [ ] Verify same introns filtered in both approaches

### Step 2: Refactor Extraction Flow
- [ ] Move `prefilter_on_metadata()` call before sequence extraction in `cli/main.py`
- [ ] Update extraction to only process prefiltered introns
- [ ] Add sequence writing during extraction (streaming)
- [ ] Test with small dataset (Chr19)

### Step 3: Update Output Writing
- [ ] Ensure `.introns.iic` written during extraction
- [ ] Skip redundant writing in `write_outputs()`
- [ ] Handle omitted introns in output (note: sequences not available)

### Step 4: Memory Profiling
- [ ] Profile with Drosophila (~200k introns) → expect ~1-2 GB
- [ ] Profile with Human (~2.1M introns) → expect ~4-5 GB
- [ ] Compare with current (28 GB)

---

## Expected Results

| Dataset | Current Peak | After Implementation | Savings |
|---------|-------------|----------------------|---------|
| C. elegans (~20k) | 3 GB | 1 GB | 67% |
| Drosophila (~200k) | 5 GB | 1.5 GB | 70% |
| Human (~2.1M) | 28 GB | 4-5 GB | **82-85%** |

---

## Risks and Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| Omitted introns lack sequences in output | Users can't analyze omitted intron sequences | Document change; original intronIC also did this |
| Filter logic bugs | Incorrect introns scored | Extensive testing; compare outputs before/after |
| Memory still high with 300k introns | 4 GB may not be enough for some systems | Could add second-level filtering or streaming |
| Output file format changes | Breaks downstream tools | Ensure .introns.iic format identical |

---

## Alternative: Most Aggressive Approach (If 4 GB Still Too High)

If 4 GB is still too much, we could:

1. **Pre-filter more aggressively** (score only longest isoform, canonical, non-duplicate)
   - Reduces 300k → 150k introns
   - Memory: 150k × 13 KB = 2 GB

2. **Stream scoring itself** (score one chromosome at a time)
   - Chr1: Extract 20k introns → score → clear
   - Chr2: Extract 18k introns → score → clear
   - Peak: Largest chromosome (~2-3 GB)

3. **Modify scorer** to use pre-extracted sequences (Approach A)
   - Highest risk but lowest memory
   - Would enable 300k × 600 bytes = 180 MB

---

## Recommendation

**Implement Approach C (Defer Extraction Until After Filtering)**

**Reasoning:**
1. **85% memory reduction** (28 GB → 4 GB) solves the problem
2. **Low risk** - no scorer changes
3. **Testable** - can verify same outputs
4. **Matches original intronIC** philosophy
5. **Incremental** - can add more optimizations later if needed

**Timeline:**
- Phase 1 (Split filtering): 2-3 hours
- Phase 2 (Refactor extraction): 2-3 hours
- Phase 3 (Update output): 1 hour
- Phase 4 (Testing): 2-3 hours
- **Total: 8-10 hours**

---

## Notes on Original intronIC v1.5.1

**Why original used less memory:**

1. **Pre-filtered annotations during parsing**
   - Only kept introns from longest isoform
   - Filtered by length during parsing
   - Never created 2.1M intron objects

2. **Streamed sequence extraction**
   - Extracted one at a time
   - Wrote immediately
   - Cleared immediately
   - Never materialized full list

3. **Monolithic approach**
   - All phases integrated
   - No intermediate materializations
   - Hard to test but memory-efficient

**Our refactor prioritized:**
- Testability → separate phases
- Maintainability → clear abstractions
- But sacrificed: Memory efficiency

**Approach C regains memory efficiency while keeping modularity!**
