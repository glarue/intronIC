# TEST COVERAGE GAPS - PRIORITIZED

## Summary
- **Current**: 457 tests (excellent unit coverage)
- **Gap**: Limited end-to-end integration and error handling tests
- **Recommendation**: Add 8-10 critical integration tests before production

---

## CRITICAL GAPS (Must Address Before Production)

### 1. Full Pipeline Test ⭐⭐⭐ (HIGHEST PRIORITY)
**Gap**: No test runs the complete workflow end-to-end with real data.

**What's needed**:
```python
def test_full_annotation_pipeline_chr19():
    """Complete pipeline: genome + annotation → classified introns."""
    # Load Chr19 test data (we already have this!)
    # Extract → Score → Normalize → Classify → Write outputs
    # Verify:
    #   - All output files created (.meta, .bed, .seqs, .scores, .log)
    #   - Some U12s found (Chr19 has ~10 known U12s)
    #   - Output files have correct format
    #   - Intron counts match across files
```

**Why critical**: This is what users actually do. If this fails, nothing works.

---

### 2. Error Recovery Tests ⭐⭐ (HIGH PRIORITY)

**Gap**: 96 error handlers in code, only 70 tested. Missing tests for:

```python
# Common user errors that should fail gracefully:

def test_annotation_references_missing_chromosome():
    """Annotation has chr1, genome only has chr2."""
    # Should: Log warning, skip those features, continue
    # Current: Untested - might crash

def test_empty_annotation_file():
    """GFF3 file has headers but no features."""
    # Should: Error with helpful message
    # Current: Untested

def test_all_introns_filtered_out():
    """All introns too short or ambiguous."""
    # Should: Warn user, continue with empty set
    # Current: Might crash downstream

def test_insufficient_reference_sequences():
    """Only 2 U12 sequences in reference (need ~50+)."""
    # Should: Error before training
    # Current: SVM might fail cryptically
```

**Why critical**: Production tools must handle bad inputs gracefully.

---

### 3. Output Validation ⭐⭐ (HIGH PRIORITY)

**Gap**: We write output files but don't validate format correctness.

```python
def test_output_files_valid_format():
    """All output files can be parsed back."""
    # Write outputs, then:
    #   - Parse .bed with BEDParser
    #   - Parse .meta as TSV
    #   - Verify intron IDs match across files

def test_output_coordinates_consistent():
    """Same intron has same coords in .bed and .meta."""

def test_bed_file_is_valid_bed():
    """Can load into IGV/UCSC Genome Browser."""
```

---

## IMPORTANT GAPS (Should Address Soon)

### 4. Cross-Component Integration ⭐

**Gap**: Components tested individually, but interactions not fully covered.

```python
def test_filtering_affects_scoring():
    """Verify filtered introns don't get scored."""

def test_duplicate_introns_output_once():
    """Duplicates tagged but only representative scored/output."""

def test_normalizer_uses_only_reference_data():
    """Z-scores fitted on reference, not experimental."""
    # This is Issue #1 - has some coverage but could be more explicit
```

---

### 5. Real-World Workflows ⭐

**Gap**: Advanced features not tested in realistic scenarios.

```python
def test_recursive_training_improves_results():
    """Two-pass training on distant species."""
    # Run once, extract confident U12s, rebuild matrices, re-run

def test_custom_pwm_matrices():
    """User provides their own PWMs."""

def test_multiple_isoforms_flag():
    """With/without -i produces different intron sets."""
```

---

## NICE-TO-HAVE GAPS (Lower Priority)

### 6. Performance Tests
- Parallel processing produces identical results
- Memory usage stays reasonable for large genomes
- Progress bars update correctly

### 7. Edge Cases
- Non-canonical boundary correction
- Zero U12s found (valid for some species)
- Very short chromosomes
- Circular chromosomes (bacteria)

---

## RECOMMENDATION: Implementation Plan

### Phase 1 (Before "Production Ready") - 4-6 hours
1. `test_full_pipeline_chr19()` - End-to-end with test data
2. `test_annotation_genome_chromosome_mismatch()` - Common error
3. `test_output_files_valid_format()` - Output validation
4. `test_empty_annotation_graceful_failure()` - Error handling

### Phase 2 (Polish) - 3-4 hours
5. `test_all_introns_filtered_out()`
6. `test_output_coordinate_consistency()`
7. `test_duplicate_handling()`
8. `test_custom_pwm_workflow()`

### Phase 3 (Optional Enhancement)
9. Performance/scalability tests
10. Advanced workflow tests

---

## Current Test Distribution (457 tests)

| Category | Count | Coverage |
|----------|-------|----------|
| M1.1 Core Models | 96 | Excellent ✓ |
| M1.2 I/O Layer | 103 | Excellent ✓ |
| M1.3 Extraction | 54 | Good ✓ |
| M1.4 Scoring | 114 | Excellent ✓ |
| M1.5 Classification | 81 | Excellent ✓ |
| M2.1 CLI | 23 | Good ✓ |
| **Integration** | **21** | **Needs work** |
| Other | 6 | Minimal |

**Analysis**: Strong unit coverage, weak integration coverage.

---

## What We Can Skip (Already Well-Covered)

✓ Unit tests for individual components (excellent)
✓ Data structure validation (Intron, GenomicCoordinate)
✓ Parser/Writer functionality (103 tests)
✓ PWM scoring logic (comprehensive)
✓ SVM training mechanics (thorough)
✓ CLI argument parsing (well-covered)

**The main gap is integration and error handling, not unit logic.**

---

## Decision Point

Should we:
1. **Add Phase 1 tests now** (4-6 hours) - Recommended before moving forward
2. **Skip for now** and continue with next milestones
3. **Prioritize differently** based on project needs

Your call!
