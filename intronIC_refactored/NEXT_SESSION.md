# Next Session: M1.3 - Extraction Pipeline

**Date:** TBD
**Branch:** `refactor/phase-0-foundation`
**Current Status:** M1.2 COMPLETE ✅

---

## What Was Just Completed (M1.2)

✅ **I/O Layer** - Full implementation with 262 passing tests
- GenomeReader (streaming + cached modes)
- Parsers (GFF3/GTF, BED, sequences)
- Writers (BED, metadata, sequences, scores, mappings)
- Tag/Attribute mapping system (hybrid compact/verbose)
- Format matches original intronIC exactly

**Key commits:**
- `a86fb28` - M1.2 COMPLETE: I/O Layer implementation
- `eec3f90` - Updated MILESTONES.md

---

## Next Milestone: M1.3 - Extraction Pipeline

**Goal:** Implement the intron extraction logic from annotations

### Components to Build

#### 1. `extraction/hierarchy.py` (~300 lines)
**Purpose:** Parse annotation hierarchy (genes → transcripts → exons)

**Key functions:**
- `build_annotation_hierarchy()` - Create parent-child relationships
- `parse_annotation_file()` - Use BioGLAnnotationParser
- Return structured hierarchy (genes contain transcripts contain exons)

**Reference:** Original `annotation_hierarchy()` at line 1567

#### 2. `extraction/intronator.py` (~200 lines)
**Purpose:** Generate introns from exon pairs

**Key functions:**
- `extract_introns_from_transcript()` - Process one transcript
- `introns_from_exon_pairs()` - Create Intron objects from adjacent exons
- Use `Intron.from_exon_pair()` factory method

**Reference:** Original `intronator()` at line 1955

#### 3. `extraction/sequences.py` (~250 lines)
**Purpose:** Extract sequences from genome

**Key functions:**
- `add_sequences_to_intron()` - Fetch from GenomeReader
- `extract_splice_sites()` - 5' and 3' sequences
- `extract_branch_point_region()` - BP search region
- `extract_flanks()` - Upstream/downstream exonic flanks
- Handle strand correctly (reverse complement for minus strand)

**Reference:** Original `get_sub_seqs()` at line 2645

#### 4. `extraction/filters.py` (~300 lines)
**Purpose:** Quality control, deduplication, filtering

**Key functions:**
- `check_omission()` - Apply omission criteria (length, ambiguous bases, etc.)
- `find_duplicates()` - Coordinate-based duplicate detection
- `mark_longest_isoforms()` - Identify longest transcript per gene
- `detect_overlaps()` - Find coordinate overlaps
- Set `intron.metadata.omitted` field appropriately

**Reference:** Original `omit_check()` at line 671, `add_tags()` at line 1897

#### 5. Integration Test
**File:** `tests/integration/test_extraction_pipeline.py`

**Tests:**
- Parse chr19 annotation → extract introns
- Verify intron count matches gold standard (20,252)
- Check coordinate ranges are valid
- Verify sequences are correct
- Test filtering logic

---

## Test-Driven Development Approach

### Step 1: Write Tests First
Create `tests/unit/test_extraction.py` with tests for:
- Hierarchy parsing
- Intron generation from exons
- Sequence extraction
- Filtering logic

### Step 2: Implement Components
Build each module to pass its tests

### Step 3: Integration Test
Full pipeline: annotation file → intron objects with sequences

---

## Key Design Principles

1. **Use existing abstractions:**
   - `GenomicCoordinate` for all coordinates
   - `Intron.from_exon_pair()` for intron creation
   - `GenomeReader` for sequence fetching
   - Composition pattern for adding data

2. **Generator-friendly:**
   - Process transcripts one at a time
   - Yield introns rather than collecting all in memory
   - Clear sequences after use if needed

3. **Type hints throughout:**
   - Full type annotations
   - Use Protocol for interfaces if needed

4. **Match original logic:**
   - Preserve extraction algorithms from original
   - Use same coordinate calculations
   - Apply same filtering criteria

---

## Success Criteria

✅ Extract 20,252 introns from chr19 annotation (matches gold standard)
✅ All coordinates within valid ranges
✅ Sequences match genome (spot checks)
✅ Filtering correctly identifies omissions
✅ Duplicate detection works
✅ All tests pass (target: 300+ total tests)

---

## Estimated Effort

- **Time:** 3-4 hours (test writing + implementation)
- **Complexity:** Medium-High (complex hierarchy logic, state management)
- **Risk:** Medium (need to preserve original extraction logic exactly)

---

## References

**Original intronIC functions:**
- `annotation_hierarchy()` - line 1567
- `intronator()` - line 1955  
- `get_sub_seqs()` - line 2645
- `omit_check()` - line 671
- `add_tags()` - line 1897
- `filter_introns_write_files()` - line 4668

**Documentation:**
- CLAUDE.md - Complete original code documentation
- REFACTORING_PLAN.md - Architecture decisions
- SESSION_NOTES.md - Implementation history

---

## Quick Start Commands

```bash
# Navigate to project
cd /home/glarue/code/intronIC/intronIC_refactored

# Run tests
PYTHONPATH=. pixi run pytest

# Check test coverage
PYTHONPATH=. pixi run pytest --cov=extraction --cov-report=term-missing

# Run specific test file
PYTHONPATH=. pixi run pytest tests/unit/test_extraction.py -xvs
```

---

## Current State

**Branch:** `refactor/phase-0-foundation`
**Tests:** 262/262 passing
**Modules Complete:**
- ✅ utils/ (coordinates, sequences)
- ✅ core/ (models, intron)
- ✅ file_io/ (genome, parsers, writers)
- ⏳ extraction/ (NEXT)

**Git clean:** Yes, all changes committed
