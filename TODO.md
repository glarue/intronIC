# intronIC Development TODO

This document tracks planned features, enhancements, and known issues for future development.

---

## Priority Levels

- 🔴 **HIGH** - Critical for next release or blocks other work
- 🟡 **MEDIUM** - Important but not blocking
- 🟢 **LOW** - Nice to have, implement when convenient
- 💭 **RESEARCH** - Requires investigation or proof of concept

---

## Feature Requests

### 🟢 Flexible Feature Selection for Training
**Status:** Deferred (post-v2.0)
**Complexity:** Moderate-High
**Issue:** Currently requires all three motifs (5'SS, BP, 3'SS). Cannot train with subset.

**Use Cases:**
- Train using only 5'SS scores (1D model)
- Train using 5'SS + BP (2D model, skip 3'SS)
- Train using 5'SS + 3'SS (2D model, skip BP)
- Feature ablation studies to understand motif importance

**Implementation Plan:** See `FEATURE_SELECTION_IMPLEMENTATION_PLAN.md`

**Scope:**
- Configuration changes
- Scoring/normalization updates
- Transformer refactor
- Training/prediction updates
- Comprehensive testing

**Files Affected:** 16 files, ~600-800 LOC changes

**Blocking Issues:** None

**Decision Criteria to Proceed:**
- User explicitly requests single-feature training
- Research question requires feature ablation
- Performance issues with BP scoring in certain species
- After v2.0 release and architecture stabilization

**Related Files:**
- Implementation plan: `FEATURE_SELECTION_IMPLEMENTATION_PLAN.md`
- Config: `config/config.yaml`, `src/intronIC/cli/config.py`
- Core logic: `src/intronIC/classification/transformers.py`

---

## Known Issues

*(No known issues currently)*

---

## Enhancements

### 🟡 [Add medium priority enhancements here]

---

## Documentation

### 🟢 [Add documentation TODOs here]

---

## Performance Optimizations

### 🟢 [Add optimization ideas here]

---

## Testing

### 🟡 Add Unit Tests for Streaming Mode
**Status:** Planned
**Plan:** [STREAMING_UNIT_TEST_PLAN.md](docs/STREAMING_UNIT_TEST_PLAN.md)

Need tests for:
- `StreamingSequenceStore` - insert, retrieve, duplicate handling, WAL mode, cleanup
- `ScoringMotifs` - creation from Intron objects, memory footprint
- Dual-mode scoring in `IntronScorer` and `BranchPointScorer`
- End-to-end streaming vs standard output comparison
- Edge cases: short introns, AT-AC boundaries, BP region exclusion

**Test Data:** Use `src/intronIC/data/test_data/Homo_sapiens.Chr19.*` (smaller, faster than full genomes)

**Priority:**
1. Core unit tests (StreamingSequenceStore, ScoringMotifs, dual-mode scoring)
2. Integration tests (streaming vs standard parity)
3. Performance tests (memory profiling, optional)

### 🟡 [Add testing improvements here]

---

## Technical Debt

### 🟢 [Add refactoring needs here]

---

## Research Ideas

### 💭 [Add research questions here]

---

## Completed Items

Items are moved here after implementation.

### ✅ Test Suite Audit and PWM Data Migration (v2.0)
**Completed:** 2025-12-01
**Branch:** `refactor/src-layout`

Audited test suite and migrated PWM data files to JSON format.

**Test Suite Changes:**
- Archived obsolete `test_pwm_format_equivalence.py` (was testing non-existent YAML format)
- Removed 3 redundant tests from `test_error_handling.py` (already covered in `test_coordinates.py`)
- Created new `test_pwm_format_equivalence.py` with 12 tests comparing legacy .iic with JSON format
- Fixed 2 edge case tests that were skipped due to API changes (`test_z_score_normalization_with_zero_variance`, `test_z_score_with_single_reference`)
- Fixed `conftest.py` data_dir path (was pointing to non-existent root `data/` instead of `src/intronIC/data/`)

**PWM Data Migration:**
- Moved `scoring_matrices.fasta.iic` to `data/archive/` (superseded by JSON format)
- `intronIC_scoring_PWMs.json` is now the default PWM file
- Added `data/archive/README.md` documenting all archived PWM files
- Updated `wiki/Training-data-and-PWMs.md` to reflect new JSON format

**Test Results:** 509 passed, 2 skipped (down from 492 passed, 22 skipped)

**Remaining Skips (Legitimate):**
- `test_normalizer.py:367` - RobustScaler math differs from StandardScaler
- `test_scorer.py:624` - Test PWM files don't have complete U2/U12 pairs

### ✅ U12 Boundary Correction Consistency Fix (v2.0)
**Completed:** 2025-11-30
**Branch:** `refactor/src-layout`
**Documentation:** [STREAMING_MODE_CONSOLIDATION.md](docs/STREAMING_MODE_CONSOLIDATION.md)

Fixed inconsistency where streaming mode reported different U12 counts than standard mode for C. elegans.

**Root Cause:** The `extract_scoring_motifs()` method in streaming mode was not properly excluding the 5' splice site region when extracting the BP search region for short introns, causing different BP scores.

**Fix Applied:**
- Added 5' splice site exclusion logic to `extract_scoring_motifs()` in `src/intronIC/core/intron.py` (matches `BranchPointScorer._extract_search_region()`)
- Both methods now use: `bp_start_idx = max(bp_start_idx, five_end)` to exclude 5'SS region

**Verified:** Both modes now produce identical results for C. elegans:
- 1 U12 intron at >90% probability (`WBGene00014227@ZK1128.1.1_6` at 99.47%)
- 7 total U12 calls (below threshold)
- 109,823 U2 introns

### ✅ Per-Contig Mode Removal (v2.0)
**Completed:** 2025-11-30
**Branch:** `refactor/src-layout`
**Documentation:** [STREAMING_MODE_CONSOLIDATION.md](docs/STREAMING_MODE_CONSOLIDATION.md)

Removed the deprecated `--per-contig` mode, consolidating into streaming mode.

**Rationale:** Per-contig mode (~18% memory savings) was superseded by streaming mode (~85% memory savings). Both use contig-by-contig processing internally, but streaming also stores sequences in SQLite.

**Changes:**
- Removed `--per-contig` CLI argument from `args.py`
- Removed `per_contig` field from `PerformanceConfig` in `config.py`
- Removed `extract_introns_from_annotation_per_contig()` function (~160 lines) from `main.py`
- Simplified routing logic in `main_classify()` to only have streaming vs standard modes
- Updated comments/docstrings referencing per-contig mode

**Result:** Two clean extraction paths:
1. **Standard mode** (default) - Fast, loads all introns into memory
2. **Streaming mode** (`--streaming`) - Memory-efficient (~85% savings), uses SQLite

### ✅ Streaming Mode Implementation (v2.0)
**Completed:** 2025-11-30
**Branch:** `refactor/src-layout`

Implemented true streaming extraction with ~85% memory savings:

**Completed Components:**
1. ✅ `ScoringMotifs` dataclass - lightweight structure (~100-150 bytes) for PWM scoring
2. ✅ `StreamingSequenceStore` - SQLite-backed storage with WAL mode for parallel workers
3. ✅ `IntronScorer` dual-mode - checks `intron.motifs` first (streaming), falls back to `intron.sequences`
4. ✅ `BranchPointScorer` dual-mode - same pattern as IntronScorer
5. ✅ `extract_introns_streaming()` in main.py - full streaming extraction function
6. ✅ CLI `--streaming` flag - enables streaming mode
7. ✅ Deferred sequence writing - both modes write introns.iic AFTER classification so SVM scores are available
8. ✅ `formatted_name` column in SQLite - stores pre-computed names for proper output format
9. ✅ Duplicate filtering - streaming mode filters duplicates when writing from SQLite

**Key Files Modified:**
- `src/intronIC/cli/main.py` - `extract_introns_streaming()`, routing in `main_classify()`
- `src/intronIC/file_io/sequence_store.py` - `StreamingSequenceStore` with `formatted_name` support
- `src/intronIC/file_io/writers.py` - `write_from_row()` updated signature
- `src/intronIC/scoring/intron_scorer.py` - dual-mode scoring
- `src/intronIC/scoring/branch_point.py` - dual-mode BP scoring
- `src/intronIC/core/intron.py` - `ScoringMotifs` dataclass

**Verified:**
- Both standard and streaming modes produce identical line counts (109,830 for C. elegans)
- SVM scores properly written to introns.iic in both modes
- Streaming mode stores 109,955 sequences in SQLite, outputs 109,830 after duplicate filtering

---

## Notes

- Keep this file updated as new issues/features are identified
- Link to detailed planning documents for complex features
- Update priority levels as project needs change
- Mark items as completed with date and commit hash
