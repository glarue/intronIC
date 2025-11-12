# Test Suite Improvement Plan
**Date:** 2025-11-11
**Project:** intronIC_refactored
**Current Status:** 7.5/10 - Good foundation, needs targeted improvements

---

## Table of Contents

1. [Overview](#overview)
2. [Current State Summary](#current-state-summary)
3. [Implementation Phases](#implementation-phases)
4. [Detailed Task Breakdown](#detailed-task-breakdown)
5. [Success Criteria](#success-criteria)
6. [Timeline](#timeline)

---

## Overview

### Goals
1. **Complete incomplete tests** - No skipped or TODO tests in critical paths
2. **Add robust error handling tests** - Handle corrupted files, edge cases gracefully
3. **Improve test organization** - Move misplaced files, add markers, improve structure
4. **Enhance maintainability** - Reduce redundancy, improve clarity, add documentation

### Non-Goals
- ❌ Gold standard regression testing (codebase has evolved beyond direct comparison)
- ❌ Performance optimization testing (focus on correctness first)
- ❌ Extensive mutation testing (overkill for current stage)

---

## Current State Summary

### Test Statistics
- **Total Lines:** ~9,163 lines across 24 test files
- **Organization:** Unit + Integration tests properly separated
- **Coverage:** Good core coverage, gaps in integration/error handling

### Strengths ✅
1. **ML Integrity Tests** - Exemplary prevention of data leakage (test_normalizer.py)
2. **Coordinate System Tests** - Outstanding coverage preventing off-by-one errors
3. **Model Tests** - Comprehensive data structure validation
4. **Real Data Integration** - Uses actual chr19 annotation and reference sequences

### Critical Issues 🔴
1. Skipped tests in `test_scorer.py` (real PWM testing)
2. Incomplete extraction pipeline tests (many TODOs)
3. Root-level test files (6 files misplaced)
4. No error handling tests (corrupted files, disk errors, etc.)
5. Limited edge case coverage (extreme sizes, unusual inputs)

---

## Implementation Phases

### Phase 1: Critical Fixes (Priority: HIGH) 🔴
**Goal:** Complete all incomplete tests and fix organizational issues
**Time Estimate:** 3-4 hours

**Tasks:**
1. Implement skipped PWM scoring tests
2. Complete extraction pipeline TODOs
3. Move root-level test files to proper location
4. Add basic error handling tests

**Impact:** HIGH - Closes critical gaps in test coverage

---

### Phase 2: Robustness & Edge Cases (Priority: MEDIUM) 🟡
**Goal:** Make tests more comprehensive and robust
**Time Estimate:** 4-5 hours

**Tasks:**
1. Add file corruption/malformation tests
2. Add edge case tests (very short/long introns, unusual genomes)
3. Add writer error tests (disk full, permissions)
4. Add scorer coordinate mismatch tests

**Impact:** MEDIUM - Improves robustness for production use

---

### Phase 3: Quality & Organization (Priority: LOW) 🟢
**Goal:** Improve test maintainability and clarity
**Time Estimate:** 2-3 hours

**Tasks:**
1. Add pytest markers (slow, integration, requires_data)
2. Reduce test redundancy with parametrize
3. Add test documentation/comments
4. Set up code coverage reporting

**Impact:** LOW - Quality of life improvements

---

## Detailed Task Breakdown

### Phase 1 Tasks

#### 1.1: Implement Skipped PWM Scoring Tests
**File:** `tests/unit/test_scoring/test_scorer.py`
**Current Issue:** Tests marked with `@pytest.mark.skip` and TODO comments

**Action:**
```python
# Remove @pytest.mark.skip decorators
# Implement tests with actual PWM matrices from data/

def test_score_u12_intron_with_real_pwm():
    """Test scoring a known U12 intron with real U12 PWM."""
    # Load actual U12 PWM from data/scoring_matrices.fasta.iic
    # Score a known U12 intron
    # Assert U12 score > U2 score
    pass

def test_score_u2_intron_with_real_pwm():
    """Test scoring a known U2 intron with real U2 PWM."""
    # Load actual U2 PWM
    # Score a known U2 intron
    # Assert U2 score > U12 score
    pass

def test_branch_point_detection_accuracy():
    """Test that BP detection finds correct position."""
    # Use known U12 intron with TCCTTAAC motif at specific position
    # Verify detected BP matches expected position
    pass
```

**Validation:** All scorer tests pass without skips

---

#### 1.2: Complete Extraction Pipeline Tests
**File:** `tests/integration/test_extraction_pipeline.py`
**Current Issue:** Many incomplete tests with TODO/commented code

**Action:**
```python
# Complete these TODOs:

def test_filtering_removes_short_introns():
    """Test that introns below minimum length are filtered."""
    # Create test annotation with introns of various lengths
    # Apply filter with min_length=30
    # Assert short introns marked as omitted with 'o:s'
    pass

def test_filtering_removes_ambiguous_bases():
    """Test that introns with N's are filtered."""
    # Create intron with ambiguous bases in scoring regions
    # Assert marked as omitted with 'o:a'
    pass

def test_duplicate_detection_and_priority():
    """Test duplicate resolution follows correct priority."""
    # Create duplicates: longest isoform vs shorter isoform
    # Assert longest isoform kept, others marked duplicate
    pass

def test_overlap_detection():
    """Test overlapping intron coordinates are detected."""
    # Create introns with overlapping coordinates
    # Assert marked as omitted with 'o:v'
    pass

def test_noncanonical_handling():
    """Test non-canonical introns handled correctly."""
    # Create GT-TG, AT-AG introns
    # Test both with and without --no_nc flag
    pass
```

**Validation:** All extraction tests pass, no TODOs remain

---

#### 1.3: Reorganize Root-Level Test Files
**Current Issue:** 6 test files in project root instead of tests/ directory

**Files to Move:**
```
intronIC_refactored/
├── test_pwm_fix.py          → tests/unit/test_scoring/test_pwm_fix.py
├── test_matrix_conversion.py → tests/unit/test_scoring/test_matrix_conversion.py
├── test_bp_scoring.py        → tests/unit/test_scoring/test_bp_scoring.py
├── test_coordinate_fix.py    → tests/unit/test_io/test_coordinate_fix.py
├── test_omission_tagging.py  → tests/unit/test_extraction/test_omission_tagging.py
└── test_type_id_fix.py       → tests/unit/test_classification/test_type_id_fix.py
```

**Action:**
1. Create subdirectories if needed
2. Move files with `git mv`
3. Update imports if necessary
4. Run tests to verify they still work
5. Remove from root .gitignore if present

**Validation:** `pytest tests/` runs all tests, none in root

---

#### 1.4: Add Basic Error Handling Tests
**New File:** `tests/unit/test_io/test_error_handling.py`

**Tests to Add:**
```python
def test_corrupted_fasta_raises_error():
    """Test that corrupted FASTA files raise appropriate error."""
    # Create truncated FASTA
    # Attempt to read
    # Assert raises ParseError with informative message
    pass

def test_malformed_gff3_raises_error():
    """Test that malformed GFF3 raises appropriate error."""
    # Create GFF3 with missing columns
    # Assert raises ParseError
    pass

def test_missing_file_raises_error():
    """Test that missing input file raises FileNotFoundError."""
    # Attempt to read non-existent file
    # Assert raises FileNotFoundError
    pass

def test_empty_annotation_raises_error():
    """Test that empty annotation file raises appropriate error."""
    # Create empty GFF3
    # Assert raises ValueError or EmptyAnnotationError
    pass

def test_writer_handles_disk_full():
    """Test writers handle disk full gracefully."""
    # Mock disk full condition
    # Attempt write
    # Assert raises IOError with clear message
    pass
```

**Validation:** Error handling tests pass

---

### Phase 2 Tasks

#### 2.1: File Corruption Tests
**File:** `tests/unit/test_io/test_parsers.py` (extend)

**Tests to Add:**
```python
def test_truncated_gzipped_fasta():
    """Test handling of truncated gzipped FASTA."""
    pass

def test_gff3_with_wrong_encoding():
    """Test GFF3 with non-UTF-8 encoding."""
    pass

def test_bed_with_invalid_coordinates():
    """Test BED file with start > stop."""
    pass

def test_mixed_line_endings():
    """Test files with mixed \\r\\n and \\n."""
    pass
```

---

#### 2.2: Edge Case Tests
**File:** `tests/unit/test_extraction/test_edge_cases.py` (new)

**Tests to Add:**
```python
def test_single_exon_gene():
    """Test gene with only one exon (no introns)."""
    pass

def test_very_short_intron():
    """Test intron of 1bp (should be filtered)."""
    pass

def test_very_long_intron():
    """Test intron of 1Mbp (memory handling)."""
    pass

def test_intron_at_chromosome_boundary():
    """Test intron at start/end of chromosome."""
    pass

def test_overlapping_genes_different_strands():
    """Test genes on opposite strands with overlapping coordinates."""
    pass

def test_nested_genes():
    """Test gene entirely within intron of another gene."""
    pass
```

---

#### 2.3: Writer Error Tests
**File:** `tests/unit/test_io/test_writers.py` (extend)

**Tests to Add:**
```python
def test_writer_permission_denied():
    """Test writer handles permission errors."""
    pass

def test_writer_path_does_not_exist():
    """Test writer handles non-existent directory."""
    pass

def test_writer_disk_quota_exceeded():
    """Test writer handles disk quota errors."""
    pass

def test_writer_invalid_characters_in_filename():
    """Test writer handles special characters in output names."""
    pass
```

---

#### 2.4: Scorer Coordinate Mismatch Tests
**File:** `tests/unit/test_scoring/test_scorer.py` (extend)

**Tests to Add:**
```python
def test_scorer_intron_too_short_for_pwm():
    """Test scorer handles intron shorter than PWM length."""
    pass

def test_scorer_pwm_extends_past_intron_boundary():
    """Test scorer handles PWM overlapping intron boundaries."""
    pass

def test_scorer_branch_point_region_too_short():
    """Test scorer handles insufficient BP search space."""
    pass
```

---

### Phase 3 Tasks

#### 3.1: Add Pytest Markers
**File:** `pyproject.toml` (update)

**Add Markers:**
```toml
[tool.pytest.ini_options]
markers = [
    "slow: marks tests as slow (>1 second) (deselect with '-m \"not slow\"')",
    "integration: marks tests as integration tests",
    "requires_data: marks tests that require chr19 data files",
    "unit: marks tests as unit tests",
]
```

**Tag Tests:**
- Add `@pytest.mark.slow` to tests taking >1 second
- Add `@pytest.mark.integration` to integration tests
- Add `@pytest.mark.requires_data` to tests needing chr19 data

---

#### 3.2: Reduce Redundancy with Parametrize
**Target Files:** Various (especially test_parsers.py, test_writers.py)

**Example Refactor:**
```python
# Before (redundant):
def test_parse_canonical_gtag():
    assert parse_dinucleotides("GTXXXXAG") == "GT-AG"

def test_parse_canonical_gcag():
    assert parse_dinucleotides("GCXXXXAG") == "GC-AG"

def test_parse_canonical_atac():
    assert parse_dinucleotides("ATXXXXAC") == "AT-AC"

# After (parametrized):
@pytest.mark.parametrize("sequence,expected", [
    ("GTXXXXAG", "GT-AG"),
    ("GCXXXXAG", "GC-AG"),
    ("ATXXXXAC", "AT-AC"),
])
def test_parse_canonical_dinucleotides(sequence, expected):
    assert parse_dinucleotides(sequence) == expected
```

**Target:** Reduce test code by ~15-20% through parametrization

---

#### 3.3: Add Test Documentation
**Action:** Add module-level docstrings to all test files

**Template:**
```python
"""
Tests for [module name].

Coverage:
- [Aspect 1]: 95%
- [Aspect 2]: 90%
- [Aspect 3]: 85%

Critical Tests:
- test_critical_behavior: Ensures [critical property]
- test_edge_case_X: Handles [edge case]

Known Gaps:
- [Gap 1]: Not tested because [reason]
- [Gap 2]: Future work

Author: intronIC test suite
Date: 2025-11-11
"""
```

---

#### 3.4: Set Up Code Coverage Reporting
**File:** `pyproject.toml` (extend)

**Add Coverage Config:**
```toml
[tool.pytest.ini_options]
addopts = [
    "--cov=.",
    "--cov-report=html",
    "--cov-report=term-missing",
]

[tool.coverage.run]
source = ["."]
omit = [
    "*/tests/*",
    "*/test_*.py",
    "*/__pycache__/*",
]

[tool.coverage.report]
precision = 2
show_missing = true
skip_covered = false
```

**Usage:**
```bash
pixi run pytest --cov
# View HTML report: open htmlcov/index.html
```

**Target Coverage Levels:**
- Core models: >95%
- Scoring/Classification: >90%
- I/O modules: >85%
- CLI: >75%

---

## Success Criteria

### Phase 1 (Critical)
- ✅ Zero skipped tests in core modules
- ✅ Zero TODO comments in test files
- ✅ All tests in tests/ directory (none in root)
- ✅ Basic error handling for all I/O operations
- ✅ All tests pass: `pytest tests/ -v`

### Phase 2 (Robustness)
- ✅ Corruption tests for all input formats
- ✅ Edge case tests for unusual genomes/introns
- ✅ Writer error handling for common failures
- ✅ Scorer handles coordinate mismatches gracefully

### Phase 3 (Quality)
- ✅ All tests properly marked (slow, integration, etc.)
- ✅ Code coverage >85% for core modules
- ✅ Test code reduced by 15-20% through parametrization
- ✅ All test files have comprehensive docstrings

---

## Timeline

### Week 1: Phase 1 (Critical Fixes)
- **Day 1-2:** Implement skipped PWM tests, complete extraction tests
- **Day 3:** Reorganize root-level test files
- **Day 4:** Add basic error handling tests
- **Day 5:** Review, fix any failures, ensure all Phase 1 criteria met

### Week 2: Phase 2 (Robustness)
- **Day 1-2:** File corruption and malformation tests
- **Day 3:** Edge case tests (short/long introns, boundaries)
- **Day 4:** Writer error tests
- **Day 5:** Scorer coordinate mismatch tests, review

### Week 3: Phase 3 (Quality)
- **Day 1:** Add pytest markers to all tests
- **Day 2:** Refactor redundant tests with parametrize
- **Day 3:** Add test documentation
- **Day 4:** Set up coverage reporting, fill coverage gaps
- **Day 5:** Final review, documentation, celebration 🎉

**Total Estimated Time:** ~9-12 hours of focused work

---

## Notes & Considerations

### What We're NOT Doing (and Why)
1. **Gold Standard Regression Tests:** Codebase has evolved with improvements (better ML integrity, enhanced output). Direct comparison to old output is no longer meaningful or desirable.
2. **Performance Benchmarks:** Current focus is correctness. Performance optimization is future work.
3. **Exhaustive Mutation Testing:** Overkill for current stage; focus on high-value test cases.

### Dependencies
- Real chr19 data in `test_data/` directory
- Reference U12/U2 sequences in `data/` directory
- PWM matrices in `data/scoring_matrices.fasta.iic`

### Risk Mitigation
- Run full test suite after each phase
- Keep test changes small and incremental
- Use git branches for each phase
- Review coverage reports regularly

---

## Appendix: Test Organization Structure

### Final Structure (After Phase 1)
```
tests/
├── __init__.py
├── conftest.py                    # Shared fixtures
├── unit/
│   ├── __init__.py
│   ├── test_models.py             # Core data structures
│   ├── test_sequences.py          # Sequence utilities
│   ├── test_coordinates.py        # Coordinate systems ⭐
│   ├── test_io/
│   │   ├── __init__.py
│   │   ├── test_parsers.py        # Input file parsing
│   │   ├── test_genome.py         # Genome reading
│   │   ├── test_writers.py        # Output writing
│   │   ├── test_error_handling.py # NEW: Error handling
│   │   └── test_coordinate_fix.py # MOVED from root
│   ├── test_scoring/
│   │   ├── __init__.py
│   │   ├── test_pwm.py            # PWM operations
│   │   ├── test_branch_point.py   # BP detection
│   │   ├── test_scorer.py         # Intron scoring (FIX skips)
│   │   ├── test_pwm_fix.py        # MOVED from root
│   │   ├── test_bp_scoring.py     # MOVED from root
│   │   └── test_matrix_conversion.py # MOVED from root
│   ├── test_normalization/
│   │   ├── __init__.py
│   │   └── test_normalizer.py     # Z-score normalization ⭐⭐
│   ├── test_classification/
│   │   ├── __init__.py
│   │   ├── test_trainer.py        # SVM training
│   │   ├── test_classifier.py     # Full classifier ⭐
│   │   ├── test_optimizer.py      # Hyperparameter optimization
│   │   ├── test_predictor.py      # Ensemble prediction
│   │   └── test_type_id_fix.py    # MOVED from root
│   └── test_extraction/
│       ├── __init__.py
│       ├── test_annotator.py      # Annotation hierarchy
│       ├── test_intronator.py     # Intron generation
│       ├── test_filters.py        # Quality filters
│       ├── test_edge_cases.py     # NEW: Edge cases
│       └── test_omission_tagging.py # MOVED from root
└── integration/
    ├── __init__.py
    ├── test_extraction_pipeline.py # COMPLETE TODOs
    ├── test_classification_pipeline.py
    └── test_cli.py                 # Command-line interface

⭐ = Exemplary test file (use as template)
```

---

**Document Version:** 1.0
**Last Updated:** 2025-11-11
**Author:** intronIC development team
