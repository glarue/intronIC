# Memory Optimization & Performance Enhancement Plan

**Date:** 2025-11-15
**Branch:** `claude/fix-scaler-centering-01C2BuWBX24F7n3cCBX1QUWu`
**Goal:** Reduce memory from ~24 GB to ~7-10 GB and restore parallel scoring

---

## Current State

### Memory Usage (Full Human Genome)
- **Current:** ~24 GB during extraction
- **Target:** ~7-10 GB peak
- **Original intronIC:** ~5 GB peak

### Performance
- **Scoring:** Sequential (single-threaded)
- **Original:** Parallelized with `-p` flag

### Breakdown of Current 24 GB
| Component | Memory | Issue |
|-----------|--------|-------|
| Genome (cached) | ~4 GB | `use_cache=True` |
| `introns_list` (coords only) | ~3 GB | **Never deleted after extraction** |
| `introns_all` (with sequences) | ~14 GB | Materialized list with full sequences |
| Python overhead | ~3 GB | Object overhead for 1M introns |
| **Total** | **~24 GB** | ❌ |

---

## Three-Phase Implementation Plan

---

## PHASE 1: Quick Memory Fix (~3 GB savings)

**Goal:** Delete `introns_list` immediately after extraction
**Risk:** Low
**Memory saved:** ~3 GB
**Files changed:** 1 (cli/main.py)

### Pre-Implementation Checks

1. **Verify no references to `introns_list` after extraction**
   - Search for all uses of `introns_list` in `extract_introns_from_annotation()`
   - Confirm last use is at line 418 (passed to `extract_sequences()`)
   - Confirm `introns_all` (line 422) replaces it entirely

2. **Verify same pattern in `extract_introns_from_bed()`**
   - Check if similar variable exists
   - Apply same fix if applicable

### Implementation Steps

1. **Modify `cli/main.py::extract_introns_from_annotation()`**
   ```python
   # After line 422:
   introns_all = list(introns_with_seq)
   messenger.log_only(f"Extracted sequences for {len(introns_all)} introns")

   # NEW: Delete old list to free memory
   del introns_list
   gc.collect()
   messenger.log_only("Freed memory from coordinate-only intron list")
   ```

2. **Add import if needed**
   ```python
   import gc  # At top of file if not present
   ```

3. **Apply same fix to `extract_introns_from_bed()` if applicable**
   - Check for similar pattern
   - Apply same deletion logic

### Testing

1. **Syntax validation**
   ```bash
   pixi run python -m py_compile cli/main.py
   ```

2. **Functional test (Chr19)**
   ```bash
   intronIC -g data/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
            -a data/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz \
            -n test_phase1 -p 4 -f cds
   ```
   - Verify all outputs generated
   - Verify no errors in log
   - Compare metrics to baseline

3. **Memory monitoring (if possible)**
   ```bash
   # In separate terminal during run:
   watch -n 2 'ps aux | grep intronIC | grep -v grep'
   ```

### Success Criteria

- ✅ Code compiles without errors
- ✅ Chr19 test completes successfully
- ✅ All output files generated
- ✅ No functional regressions
- ✅ Memory reduced by ~3 GB (if measurable)

### Rollback Plan

If issues arise:
```bash
git diff cli/main.py  # Review changes
git checkout cli/main.py  # Revert if needed
```

---

## PHASE 2: Medium Memory Fix (~7-10 GB savings)

**Goal:** Write sequences during extraction, clear immediately, keep only scoring regions
**Risk:** Medium
**Memory saved:** ~7-10 GB
**Files changed:** 2-3 (cli/main.py, possibly new helper function)

### Pre-Implementation Analysis

1. **Understand sequence lifecycle**
   - Extraction: Full sequences populated
   - Boundary correction: Needs full sequences (upstream_flank, seq, downstream_flank)
   - Scoring: Needs full sequences (extracts regions on-demand)
   - Classification: Only needs scoring regions (five_seq, three_seq, bp_seq)
   - Output: Needs full sequences for `.introns.iic`, scoring regions for `.score_info.iic`

2. **Identify "write-during-extraction" pattern from original**
   - Original wrote sequences in `filter_introns_write_files()` generator loop
   - Set `intron.seq = None` immediately after writing (line 4847 in original)
   - Result: `final_introns` list had seq=None, small memory footprint

3. **Determine clearing strategy**
   - **Option A:** Write during initial extraction, clear immediately
     - Problem: Boundary correction needs sequences
   - **Option B:** Write after boundary correction, clear before scoring
     - Problem: Scorer needs full sequences to extract regions on-demand
   - **Option C:** Write after scoring, clear before classification
     - Best option: Scoring complete, classification only needs scored regions

### Implementation Strategy

**Approach: Early Write + Late Clear**

1. Write sequences to `.introns.iic` after scoring completes
2. Clear large sequences before classification
3. Keep only: `five_seq`, `three_seq`, `bp_seq`, `bp_seq_u2`, `five_prime_dnt`, `three_prime_dnt`, `bp_relative_coords`
4. Use dual-track approach for classification (cleared copy) and output (keep originals OR re-read from file)

### Implementation Steps

#### Step 2.1: Create sequence clearing function

**File:** `cli/main.py`

```python
def clear_large_sequences_for_classification(introns: List[Intron]) -> List[Intron]:
    """
    Clear large sequence fields before classification to reduce memory.

    Clears:
    - seq (full intron sequence)
    - upstream_flank (exonic context)
    - downstream_flank (exonic context)
    - bp_region_seq (branch point search region)

    Keeps:
    - five_seq, three_seq, bp_seq, bp_seq_u2 (scored sequences)
    - five_prime_dnt, three_prime_dnt (terminal dinucleotides)
    - bp_relative_coords (branch point position)
    - All scores and metadata

    Args:
        introns: List of scored introns with full sequences

    Returns:
        New list of introns with large sequences cleared
    """
    from dataclasses import replace

    cleared = []
    for intron in introns:
        # Create new intron with cleared sequences
        new_sequences = replace(
            intron.sequences,
            seq=None,
            upstream_flank=None,
            downstream_flank=None,
            bp_region_seq=None
        )
        cleared_intron = replace(intron, sequences=new_sequences)
        cleared.append(cleared_intron)

    return cleared
```

#### Step 2.2: Write sequences after scoring

**File:** `cli/main.py` - in `run_training_mode()` function

Find where scoring completes (around line 1519), add sequence writing:

```python
# After scoring completes (line 1519+)
scored_introns = score_introns(introns_for_scoring, config, messenger, reporter)

# NEW: Write sequences to file NOW (while still in memory)
messenger.info("Writing intron sequences to file")
seq_output_path = output_dir / f"{base_name}.introns.iic"
from output.sequence_writer import SequenceWriter

seq_writer = SequenceWriter(seq_output_path)
with seq_writer:
    for intron in scored_introns:
        seq_writer.write_intron(
            intron,
            include_score=False,  # Don't write scores yet (not classified)
            write_omitted=True    # Write all for now
        )
messenger.log_only(f"Wrote sequences for {len(scored_introns)} introns")
```

#### Step 2.3: Clear sequences before classification

**File:** `cli/main.py` - before classification step

```python
# NEW: Clear large sequences before classification
messenger.log_only("Clearing large sequences to reduce memory for classification")
scored_introns_for_classification = clear_large_sequences_for_classification(scored_introns)

# Use cleared introns for classification
if pretrained_model:
    classified = classify_with_pretrained_model(
        scored_introns_for_classification,  # Use cleared introns
        ...
    )
```

#### Step 2.4: Update output writing to skip sequences

**File:** `cli/main.py` - in output section

```python
# In write_outputs() call:
write_outputs(
    classified_introns,
    output_dir,
    base_name,
    config,
    messenger,
    skip_sequences=True,  # NEW: Skip because already written
    ...
)
```

Modify `write_outputs()` to handle `skip_sequences`:

```python
def write_outputs(..., skip_sequences: bool = False, ...):
    ...
    if not skip_sequences:
        # Write sequences (sequences-only mode or not already written)
        seq_writer.write_all(...)
    else:
        messenger.log_only(f"Sequences already written, skipping: {seq_path}")
```

### Testing

1. **Syntax validation**
   ```bash
   pixi run python -m py_compile cli/main.py
   ```

2. **Functional test (Chr19)**
   ```bash
   intronIC -g data/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
            -a data/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz \
            -n test_phase2 -p 4 -f cds
   ```
   - Verify `.introns.iic` written after scoring
   - Verify classification completes successfully
   - Verify all output files generated
   - Compare output to Phase 1 output (should be identical)

3. **Memory monitoring**
   - Peak should be ~10-12 GB (after Phase 1 + Phase 2)

### Success Criteria

- ✅ Code compiles without errors
- ✅ Chr19 test completes successfully
- ✅ Sequences written after scoring
- ✅ Classification uses cleared introns
- ✅ All output files generated and match Phase 1
- ✅ Memory reduced by additional ~7-10 GB
- ✅ No functional regressions

### Rollback Plan

```bash
git diff cli/main.py  # Review changes
git stash  # Save changes if needed
# Or commit before Phase 3 to create checkpoint
```

---

## PHASE 3: Parallel Scoring Restoration

**Goal:** Restore multiprocessing for scoring when `-p > 1`
**Risk:** Low-Medium
**Performance gain:** 4-8x speedup with `-p 12`
**Memory impact:** ~120-150 MB additional (negligible)
**Files changed:** 2-3 (cli/main.py, possibly scoring/scorer.py)

### Pre-Implementation Analysis

1. **Review original parallel scoring**
   - Used `multiprocessing.Pool.starmap()`
   - Passed one intron per worker invocation
   - Shared PWM matrices via `itertools.repeat()`
   - Memory-efficient: no large chunks, matrices shared

2. **Check current scorer interface**
   - `IntronScorer.score_intron(intron)` - scores single intron
   - Returns new intron with scores populated
   - Functional style (no side effects)

3. **Identify challenges**
   - Need to pass scorer configuration to workers
   - PWM matrices might not pickle (check)
   - Progress bar updates need thread-safe mechanism
   - Error handling per intron

### Implementation Strategy

**Approach: Pool.starmap with wrapper function**

1. Create a module-level worker function (picklable)
2. Pass scorer configuration (PWM paths, coords) instead of scorer instance
3. Each worker creates its own scorer instance
4. Use `pool.starmap()` with one intron per call
5. Collect results, handle errors gracefully

### Implementation Steps

#### Step 3.1: Create picklable worker function

**File:** `cli/main.py` (at module level, before any classes/main functions)

```python
def _score_intron_worker(
    intron: Intron,
    pwm_file: Path,
    u2_bp_file: Path,
    five_coords: Tuple[int, int],
    bp_coords: Tuple[int, int],
    three_coords: Tuple[int, int],
    ignore_nc_dnts: bool,
    pseudocount: float
) -> Tuple[Optional[Intron], Optional[str]]:
    """
    Worker function for parallel scoring.

    Returns:
        (scored_intron, error_message) tuple
        If scoring fails, scored_intron is None and error_message is set
    """
    try:
        # Import here to avoid issues with multiprocessing pickling
        from scoring.scorer import IntronScorer
        from scoring.pwm import PWMLoader, PWMSet
        from pathlib import Path

        # Load PWMs (cached by OS file system for efficiency)
        pwm_sets = PWMLoader.load_from_file(pwm_file, pseudocount=pseudocount)

        # Load U2 BP matrix if available
        if u2_bp_file.exists():
            u2_bp_matrices = PWMLoader.load_from_file(u2_bp_file, pseudocount=pseudocount)
            if 'bp' in u2_bp_matrices and u2_bp_matrices['bp'].u2_gtag:
                updated_matrices = dict(pwm_sets['bp'].matrices)
                updated_matrices[('u2', 'gtag')] = u2_bp_matrices['bp'].u2_gtag
                pwm_sets['bp'] = PWMSet(matrices=updated_matrices)

        # Create scorer
        scorer = IntronScorer(
            pwm_sets=pwm_sets,
            five_coords=five_coords,
            bp_coords=bp_coords,
            three_coords=three_coords,
            ignore_nc_dnts=ignore_nc_dnts
        )

        # Score intron
        scored = scorer.score_intron(intron)
        return (scored, None)

    except Exception as e:
        # Return error instead of crashing worker
        return (None, f"{intron.intron_id}: {str(e)}")
```

#### Step 3.2: Modify score_introns() to use parallel scoring

**File:** `cli/main.py` - function `score_introns()`

```python
def score_introns(
    introns: List[Intron],
    config: IntronICConfig,
    messenger: 'UnifiedMessenger',
    reporter: IntronICProgressReporter
) -> List[Intron]:
    """Score introns with PWM matrices (parallel or sequential)."""

    messenger.info("Loading PWM matrices")

    # Determine PWM file paths
    pwm_file = Path(__file__).parent.parent / "data" / "scoring_matrices.fasta.iic"
    if not pwm_file.exists():
        raise FileNotFoundError(f"PWM file not found: {pwm_file}")

    u2_bp_file = Path(__file__).parent.parent / "data" / "u2.conserved_empirical_bp_pwm.iic"

    # Extract scoring configuration
    five_coords = (
        config.scoring.scoring_regions.five_start,
        config.scoring.scoring_regions.five_end
    )
    bp_coords = (
        config.scoring.scoring_regions.bp_start,
        config.scoring.scoring_regions.bp_end
    )
    three_coords = (
        config.scoring.scoring_regions.three_start,
        config.scoring.scoring_regions.three_end
    )
    ignore_nc_dnts = config.scoring.ignore_nc_dnts
    pseudocount = config.scoring.pseudocount

    # Parallel or sequential scoring
    n_workers = config.runtime.n_workers

    if n_workers > 1:
        # PARALLEL SCORING
        from multiprocessing import Pool
        from itertools import repeat

        messenger.info(f"Calculating PWM scores (parallel, {n_workers} workers)")

        progress = reporter.create_progress()
        scored_introns = []
        failed_count = 0

        with progress:
            task = progress.add_task(
                "[cyan]Scoring introns...",
                total=len(introns)
            )

            with Pool(processes=n_workers) as pool:
                try:
                    # starmap returns results in order
                    results = pool.starmap(
                        _score_intron_worker,
                        zip(
                            introns,
                            repeat(pwm_file),
                            repeat(u2_bp_file),
                            repeat(five_coords),
                            repeat(bp_coords),
                            repeat(three_coords),
                            repeat(ignore_nc_dnts),
                            repeat(pseudocount)
                        )
                    )

                    # Process results
                    for scored_intron, error in results:
                        if error is not None:
                            messenger.warning(f"Failed to score intron: {error}")
                            failed_count += 1
                        else:
                            scored_introns.append(scored_intron)

                        progress.update(task, advance=1)

                except KeyboardInterrupt:
                    messenger.warning("User interrupt - terminating workers")
                    pool.terminate()
                    raise
                finally:
                    pool.close()
                    pool.join()

    else:
        # SEQUENTIAL SCORING (original logic)
        from scoring.scorer import IntronScorer
        from scoring.pwm import PWMLoader, PWMSet

        messenger.info("Calculating PWM scores (sequential)")

        # Load PWMs once
        pwm_sets = PWMLoader.load_from_file(pwm_file, pseudocount=pseudocount)

        if u2_bp_file.exists():
            u2_bp_matrices = PWMLoader.load_from_file(u2_bp_file, pseudocount=pseudocount)
            if 'bp' in u2_bp_matrices and u2_bp_matrices['bp'].u2_gtag:
                updated_matrices = dict(pwm_sets['bp'].matrices)
                updated_matrices[('u2', 'gtag')] = u2_bp_matrices['bp'].u2_gtag
                pwm_sets['bp'] = PWMSet(matrices=updated_matrices)
                messenger.log_only("Loaded conserved U2-type BP matrix")

        # Create scorer
        scorer = IntronScorer(
            pwm_sets=pwm_sets,
            five_coords=five_coords,
            bp_coords=bp_coords,
            three_coords=three_coords,
            ignore_nc_dnts=ignore_nc_dnts
        )

        # Score sequentially
        progress = reporter.create_progress()
        scored_introns = []
        failed_count = 0

        with progress:
            task = progress.add_task(
                "[cyan]Scoring introns...",
                total=len(introns)
            )

            for intron in introns:
                try:
                    scored = scorer.score_intron(intron)
                    scored_introns.append(scored)
                except Exception as e:
                    messenger.warning(
                        f"Failed to score intron {intron.intron_id}: {str(e)}. Skipping."
                    )
                    failed_count += 1

                progress.update(task, advance=1)

    # Common completion logic
    total_attempted = len(introns)
    total_scored = len(scored_introns)

    if failed_count > 0:
        messenger.warning(
            f"Scoring completed with {failed_count} failures "
            f"({total_scored}/{total_attempted} succeeded)"
        )
    else:
        messenger.log_only(f"Successfully scored {total_scored} introns")

    return scored_introns
```

#### Step 3.3: Add necessary imports

**File:** `cli/main.py` (top of file)

```python
from multiprocessing import Pool  # If not already present
from itertools import repeat  # If not already present
```

### Testing

#### Test 3.1: Sequential scoring (baseline)

```bash
intronIC -g data/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
         -a data/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz \
         -n test_phase3_seq -p 1 -f cds
```

- Verify uses sequential path
- Time the run
- Save outputs as baseline

#### Test 3.2: Parallel scoring (4 workers)

```bash
intronIC -g data/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
         -a data/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz \
         -n test_phase3_par4 -p 4 -f cds
```

- Verify uses parallel path
- Time the run (should be ~2-3x faster)
- Compare outputs to sequential (should be identical)

#### Test 3.3: Parallel scoring (12 workers)

```bash
intronIC -g data/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
         -a data/test_data/Homo_sapiens.Chr19.Ensembling_91.gff3.gz \
         -n test_phase3_par12 -p 12 -f cds
```

- Verify uses parallel path
- Time the run
- Compare outputs to sequential

#### Test 3.4: Error handling

Create a test with intentionally malformed introns to verify error handling.

### Success Criteria

- ✅ Code compiles without errors
- ✅ Sequential mode works (p=1)
- ✅ Parallel mode works (p=4, p=12)
- ✅ Parallel outputs identical to sequential
- ✅ Parallel mode shows speedup (2-4x with p=4)
- ✅ Error handling works (doesn't crash on bad introns)
- ✅ Memory overhead negligible (~100-200 MB max)

### Rollback Plan

```bash
git diff cli/main.py  # Review changes
git stash  # Save if needed
```

---

## Integration & Final Testing

After all three phases complete:

### Test 1: Full pipeline with all optimizations (Chr19)

```bash
intronIC -g data/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
         -a data/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz \
         -n test_all_phases -p 12 -f cds \
         --pretrained_model run_tests/homo_sapiens.default.model.pkl
```

**Expected results:**
- Memory: ~10-12 GB peak (down from 24 GB)
- Runtime: Faster than before (parallel scoring)
- All outputs generated correctly
- Metrics match baseline

### Test 2: Memory monitoring during full run

```bash
# Terminal 1: Run intronIC
intronIC -g <full_genome> -a <full_annotation> -n test_full -p 12 -f cds

# Terminal 2: Monitor memory
watch -n 2 'ps aux | grep intronIC | grep -v grep | awk "{print \$6/1024 \" MB\"}"'
```

**Monitor at each stage:**
- Initial extraction: Should stay under 15 GB
- Boundary correction: Should stay under 15 GB
- Scoring: Should stay under 12 GB (parallel workers)
- Classification: Should stay under 10 GB (cleared sequences)

### Test 3: Performance benchmarking

Run same dataset with different worker counts:

```bash
for p in 1 4 8 12; do
    echo "Testing with -p $p"
    time intronIC -g data.fa.gz -a data.gff3.gz -n test_p${p} -p $p -f cds
done
```

Expected scaling:
- p=1: Baseline time
- p=4: ~3-3.5x faster
- p=8: ~5-6x faster
- p=12: ~6-8x faster

### Success Criteria (All Phases)

- ✅ Memory reduced from ~24 GB to ~10 GB (60% reduction)
- ✅ Scoring parallelized successfully
- ✅ 4-8x speedup with 12 workers
- ✅ All outputs identical to baseline
- ✅ No functional regressions
- ✅ Code well-documented and maintainable

---

## Documentation Updates

After all phases complete, update:

1. **docs/MEMORY_OPTIMIZATION.md**
   - Update "Current Memory Usage" section with new numbers
   - Mark optimizations as implemented
   - Update recommendations

2. **CLAUDE.md**
   - Update performance stats
   - Note parallel scoring capability

3. **This file (MEMORY_AND_PERFORMANCE_FIX_PLAN.md)**
   - Add final results section
   - Document actual memory/performance gains
   - Add lessons learned

---

## Risk Mitigation

### General Strategies

1. **Commit after each phase**
   ```bash
   git add -A
   git commit -m "Phase X: [description]"
   ```

2. **Keep baseline outputs for comparison**
   ```bash
   # Before Phase 1
   intronIC -g chr19.fa.gz -a chr19.gff3.gz -n baseline -p 4 -f cds
   mv baseline.* baseline_outputs/
   ```

3. **Test incrementally**
   - Don't proceed to next phase if current phase fails
   - Fix issues before moving forward

4. **Monitor logs carefully**
   - Check for warnings/errors after each test
   - Compare log structure to baseline

### Phase-Specific Risks

**Phase 1:**
- Risk: Accidentally delete wrong variable
- Mitigation: Careful code review, search for all uses

**Phase 2:**
- Risk: Clearing sequences too early breaks scorer
- Mitigation: Clear only AFTER scoring completes, extensive testing

**Phase 3:**
- Risk: Multiprocessing errors, pickle issues
- Mitigation: Worker function at module level, careful imports, error handling

---

## Timeline Estimate

- **Phase 1:** 15-30 minutes (simple deletion)
- **Phase 2:** 45-90 minutes (write logic + clearing)
- **Phase 3:** 60-90 minutes (parallel scoring + testing)
- **Integration testing:** 30-60 minutes
- **Total:** 2.5-4.5 hours

---

## Appendix: Memory Calculations

### Current State (24 GB)
```
Genome cached:        4 GB
introns_list:         3 GB  ← Phase 1 removes
introns_all (seqs):  14 GB  ← Phase 2 reduces to ~4 GB
Python overhead:      3 GB
---------------------------------
Total:               24 GB
```

### After Phase 1 (21 GB)
```
Genome cached:        4 GB
introns_all (seqs):  14 GB
Python overhead:      3 GB
---------------------------------
Total:               21 GB  (-3 GB)
```

### After Phase 2 (10 GB)
```
Genome cached:        4 GB
introns (cleared):    4 GB  ← Full seqs cleared
Scoring regions:      1 GB  ← Small sequences only
Python overhead:      1 GB  ← Reduced
---------------------------------
Total:               10 GB  (-11 GB from Phase 1)
```

### After Phase 3 (10-11 GB)
```
Genome cached:        4 GB
introns (cleared):    4 GB
Scoring regions:      1 GB
Python overhead:      1 GB
Worker overhead:    0.5 GB  ← New, minimal
---------------------------------
Total:               10.5 GB  (+0.5 GB, acceptable)
```

**Final savings: 24 GB → 10.5 GB = 56% reduction**

---

## Success Definition

This plan succeeds if:

1. ✅ Memory usage reduced to ~10 GB peak (from 24 GB)
2. ✅ Scoring parallelized and shows 4-8x speedup
3. ✅ All outputs match baseline (no functional changes)
4. ✅ Code remains maintainable and well-documented
5. ✅ No regressions in error handling or edge cases
6. ✅ Changes are backward compatible with existing configs

Ready to proceed with Phase 1!
