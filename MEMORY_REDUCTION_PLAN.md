# Memory Reduction Plan: 30 GB → 5-10 GB

## Problem Statement

**Current state:** Processing human genome uses ~30 GB peak memory
**Target state:** Reduce to ~5-10 GB (similar to v1.5.1)
**Impact:** Enable processing on machines with 16-32 GB RAM

## Root Cause Analysis

### Critical Bottleneck: `cli/main.py:673`
```python
introns_with_seq = sequence_extractor.extract_sequences(introns)
introns_all = list(introns_with_seq)  # ← Materializes 1M introns × 13 KB = 13 GB
```

### Memory Breakdown (30 GB total)
- **13 GB**: All intron sequences materialized into list
- **7-9 GB**: Sequences for introns that get filtered (80% omitted/duplicated)
- **3-5 GB**: Sequence duplication during parallel classification (12 workers)
- **3 GB**: Genome chromosome in memory (already optimized with `use_cache=False`)

### Per-Intron Memory Footprint
- `sequences.seq`: ~500 bytes (full intron sequence)
- `sequences.upstream_flank`: ~200 bytes
- `sequences.downstream_flank`: ~200 bytes
- `sequences.bp_region_seq`: ~50 bytes
- Python object overhead: ~12 KB
- **Total: ~13 KB per intron**
- **For 1M introns: 13 GB**

### Why v1.5.1 Used Only 5-10 GB

**Key difference: Generator pattern with immediate clearing**

`archive/v1.5.1/intronIC/intronIC.py:4668-4847`:
```python
for intron in get_sub_seqs(all_introns, genome, ...):  # Generator
    intron.omit_check(...)
    intron = add_tags(intron, ...)

    # Write sequence immediately
    seq_file.write(output_format(intron, 'SEQ', ...) + '\n')

    # Clear sequence immediately
    intron.seq = None  # ← KEY: Only keep metadata

    final_introns.append(intron)  # No sequences
```

**Memory during loop:**
- Current intron with sequences: ~10 KB
- Accumulated introns without sequences: ~1 GB
- **Peak: ~5-10 GB** ✅

## Solution Strategy

**Replace:** `Extract ALL → Filter ALL → Score ALL → Write ALL → Clear ALL`
**With:** `Extract ONE → Write → Clear → Repeat` (streaming pattern)

## Implementation Phases

### Phase 4: Implement Sequence Clearing Method [Foundation]
**File:** `core/intron.py`
**Complexity:** Low
**Impact:** Enables Phases 1 & 3

Add method to clear large sequence fields while preserving small scoring fields:

```python
def clear_sequences(self) -> 'Intron':
    """Return a new Intron with large sequences cleared.

    Clears: seq, upstream_flank, downstream_flank, bp_region_seq (~950 bytes)
    Keeps: five_seq, three_seq, bp_seq, bp_seq_u12 (~40 bytes - needed for scoring)

    Returns:
        New Intron instance with sequences cleared
    """
```

**Tasks:**
- [ ] Add `clear_sequences()` method to Intron class
- [ ] Add `clear_all_sequences()` method (aggressive version)
- [ ] Add unit tests for both methods
- [ ] Verify frozen dataclass works with replace()

**Testing:**
```python
# Test that clearing works
intron_with_seq = create_test_intron(seq="ATCG"*100)
intron_cleared = intron_with_seq.clear_sequences()
assert intron_cleared.sequences.seq is None
assert intron_cleared.sequences.five_seq == intron_with_seq.sequences.five_seq  # Kept
```

---

### Phase 5: Add Streaming Writer Support [Foundation]
**File:** `file_io/writers.py`
**Complexity:** Low
**Impact:** Enables Phase 1

Add context manager and single-intron writing to writers:

```python
class SequenceWriter:
    def __enter__(self):
        self.file_handle = open(self.file_path, 'w')
        return self

    def write_single(self, intron: Intron):
        """Write a single intron's sequence immediately."""
        seq_line = self._format_sequence(intron)
        self.file_handle.write(seq_line + '\n')

    def __exit__(self, *args):
        if self.file_handle:
            self.file_handle.close()
```

**Tasks:**
- [ ] Add `__enter__` and `__exit__` to SequenceWriter
- [ ] Add `write_single()` method to SequenceWriter
- [ ] Update MetaWriter similarly if needed
- [ ] Add unit tests for streaming writes
- [ ] Verify output format matches batch write format

**Testing:**
```python
# Test that streaming produces identical output
with SequenceWriter(path) as writer:
    for intron in introns:
        writer.write_single(intron)

# Compare with original batch write
assert files_identical(streaming_output, batch_output)
```

---

### Phase 1: Restructure Extraction Pipeline [18 GB savings]
**File:** `cli/main.py:668-680`
**Complexity:** Medium
**Impact:** PRIMARY - Eliminates 13 GB materialization

Replace list materialization with streaming generator pattern:

**Current:**
```python
introns_with_seq = sequence_extractor.extract_sequences(introns)
introns_all = list(introns_with_seq)  # ← REMOVE
```

**New:**
```python
introns_with_seq = sequence_extractor.extract_sequences(filtered_introns)

final_introns = []
with SequenceWriter(seq_output_path) as seq_writer:
    for intron in introns_with_seq:
        # Write sequence immediately
        if not config.no_seqs:
            seq_writer.write_single(intron)

        # Clear large sequences immediately
        intron = intron.clear_sequences()

        # Keep only metadata
        final_introns.append(intron)

# Continue with scoring using final_introns (no large sequences)
```

**Tasks:**
- [ ] Identify all `list(introns_with_seq)` calls (lines 673, 802, 899)
- [ ] Replace with streaming loop at line 673
- [ ] Open SequenceWriter in context manager
- [ ] Write sequences immediately in loop
- [ ] Clear sequences after writing
- [ ] Update downstream code to use `final_introns`
- [ ] Verify scoring still works without large sequences
- [ ] Integration test with small genome (C. elegans)

**Testing:**
```bash
# Run with small dataset
timeout 60 pixi run intronIC train \
    -g data/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
    -a data/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz \
    -n test_streaming \
    -o test_streaming_output

# Compare outputs
diff test_streaming_output/*.iic archive/expected_output/*.iic
```

**Memory profiling:**
```bash
# Before and after comparison
mprof run pixi run intronIC train <args>
mprof plot  # Should show ~13 GB reduction
```

---

### Phase 3: Clear Sequences Before Classification [3-5 GB savings]
**File:** `cli/main.py:1262-1274`
**Complexity:** Low
**Impact:** Eliminates duplication across parallel workers

Clear sequences before parallel classification to reduce per-worker memory:

**Current:**
```python
classified_introns = list(predictor.predict(ensemble, normalized_introns))
# Sequences still in memory during parallel classification
```

**New:**
```python
# Clear any remaining sequences before parallelization
introns_no_seq = [i.clear_all_sequences() for i in normalized_introns]
classified_introns = list(predictor.predict(ensemble, introns_no_seq))
```

**Tasks:**
- [ ] Add sequence clearing before `predictor.predict()` call (line 1274)
- [ ] Verify classification doesn't need sequence data
- [ ] Test with full pipeline
- [ ] Profile memory reduction during parallel phase

**Testing:**
```bash
# Memory profiling should show reduction during classification phase
mprof run pixi run intronIC classify <args>
# Peak during classification should be ~400 MB instead of 4-6 GB
```

---

### Phase 2: Pre-Filter Before Sequence Extraction [7-9 GB savings]
**Files:** `extraction/filters.py`, `cli/main.py:650-670`
**Complexity:** Medium
**Impact:** OPTIONAL - Further optimization

Split filtering into metadata-based (pre-extraction) and sequence-based (post-extraction):

**Metadata filters (no sequences needed):**
- Length checks (min/max intron size)
- Coordinate validation
- Basic canonical checks (if splice sites in coordinates)

**Sequence filters (need sequences):**
- Duplicate detection (sequence comparison)
- Overlap detection
- N-content checks

**Tasks:**
- [ ] Analyze filter logic in `extraction/filters.py`
- [ ] Split into `prefilter_on_metadata()` and `filter_with_sequences()`
- [ ] Add pre-filtering step in main.py before extraction
- [ ] Ensure duplicate/overlap detection still works
- [ ] Integration test to verify same filtering results

**Implementation:**
```python
# In cli/main.py
introns_prefiltered = filter_engine.prefilter_on_metadata(introns)
# Now extract ~300k instead of 1M
introns_with_seq = sequence_extractor.extract_sequences(introns_prefiltered)

# Continue with sequence-based filtering
for intron in introns_with_seq:
    filter_engine.apply_sequence_filters(intron)
    ...
```

**Testing:**
```python
# Verify same introns filtered before vs after
original_filtered = old_filter(all_introns)
new_prefiltered = prefilter_on_metadata(all_introns)
new_final = filter_with_sequences(extract_sequences(new_prefiltered))

assert set(original_filtered) == set(new_final)
```

---

## Expected Memory Improvements

| Component | Current | After Phases 4+5+1 | After All | Savings |
|-----------|---------|-------------------|-----------|---------|
| Extraction | 13 GB (materialized) | 13 KB (streaming) | 13 KB | **13 GB** |
| Filtered waste | 9 GB (80% unused) | 9 GB | 3 GB | **6 GB** |
| Classification | 5 GB (worker duplication) | 5 GB | 400 MB | **5 GB** |
| Genome | 3 GB | 3 GB | 3 GB | 0 |
| **TOTAL** | **~30 GB** | **~17 GB** | **~6 GB** | **~24 GB** |

**Priority implementation:** Phases 4+5+1 → 18 GB savings (43% reduction)
**Full implementation:** All phases → 24 GB savings (80% reduction)

## Testing Strategy

### Unit Tests
- `test_intron_clear_sequences()` - Verify clearing methods work
- `test_streaming_writers()` - Verify single-write produces correct output
- `test_generator_pipeline()` - Verify streaming doesn't lose data

### Integration Tests
- Small genome (Chr19): Verify identical output files
- C. elegans: Verify identical metrics/classifications
- Drosophila: Memory profiling to verify <10 GB

### Memory Profiling
```bash
# Install memory profiler
pip install memory-profiler

# Profile memory usage
mprof run pixi run intronIC train \
    -g data/genomes/human.fa.gz \
    -a data/annotations/human.gff3.gz \
    -n human_test \
    -o human_test_output

# Generate plot
mprof plot --output memory_profile.png

# Check peak memory
mprof peak  # Should be <10 GB
```

### Validation Criteria
✅ Output files identical (sequences, scores, classifications)
✅ No sequences lost during streaming
✅ Peak memory <10 GB for human genome
✅ All existing tests pass
✅ Runtime increase <20%

## Implementation Order

1. **Phase 4** - Add `clear_sequences()` method (30 min)
2. **Phase 5** - Add streaming writers (45 min)
3. **Phase 1** - Restructure extraction pipeline (2-3 hours)
4. **Integration testing** - Verify correctness (1 hour)
5. **Memory profiling** - Verify improvements (30 min)
6. **Phase 3** - Clear before classification (30 min)
7. **Phase 2** (optional) - Pre-filtering (2-3 hours)

**Estimated time:** 5-6 hours for critical path (Phases 4+5+1+3)
**Estimated time:** 8-10 hours for full implementation (all phases)

## Risks and Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| Sequences needed after clearing | Breaks scoring | Keep small scoring fields in `clear_sequences()` |
| File I/O slowdown | Runtime increase | Use buffered writes, benchmark I/O |
| Generator state bugs | Lost/duplicate introns | Comprehensive integration tests |
| Breaking existing tests | Development friction | Run test suite after each phase |
| Memory profiling inaccurate | Can't verify improvements | Use multiple profiling tools (mprof, tracemalloc) |

## Rollback Plan

Each phase is independent and can be reverted:
- Phase 4: Remove `clear_sequences()` method
- Phase 5: Remove streaming writer methods
- Phase 1: Restore `list(introns_with_seq)` materialization
- Phase 3: Remove clearing before classification
- Phase 2: Restore original filter order

Git branch strategy:
```bash
# Create feature branch
git checkout -b feature/memory-reduction

# Commit after each phase
git commit -m "Phase 4: Add sequence clearing methods"
git commit -m "Phase 5: Add streaming writer support"
git commit -m "Phase 1: Implement streaming extraction pipeline"

# If issues arise, revert specific commits
git revert <commit-hash>
```

## Success Criteria

### Must Have
- ✅ Peak memory <10 GB for human genome (from 30 GB)
- ✅ All output files identical to current implementation
- ✅ All existing tests pass
- ✅ No data loss (same number of introns processed)

### Nice to Have
- ✅ Runtime increase <20%
- ✅ Memory <6 GB with Phase 2 implementation
- ✅ Clear memory profiling graphs showing improvements
- ✅ Documentation updated with new architecture

## Future Optimizations

If still >10 GB after all phases, consider:

1. **Chromosome-by-chromosome processing** (v1.5.1 approach)
   - Process chr1 completely, free memory
   - Process chr2 completely, free memory
   - Requires annotation reorganization
   - Could reach 2-3 GB peak

2. **Disk-backed storage** for large intermediate data
   - Store filtered introns on disk between phases
   - Memory-map score matrices
   - Use SQLite for intron metadata

3. **Streaming classification** without batching
   - Process introns one-at-a-time through classifier
   - Tradeoff: Slower but lower memory

4. **Compressed in-memory storage**
   - Compress sequence strings with zlib
   - Decompress only when needed
   - Tradeoff: CPU for memory

## References

- Original investigation: (from Task agent analysis)
- v1.5.1 implementation: `archive/v1.5.1/intronIC/intronIC.py:4668-4847`
- Current memory optimization: `docs/MEMORY_OPTIMIZATION.md`
- Enum-based storage: `ENUM_BASED_STORAGE_IMPLEMENTATION.md`

---

**Document created:** 2025-11-20
**Target completion:** TBD
**Status:** Planning complete, ready for implementation
