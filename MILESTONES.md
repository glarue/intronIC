# intronIC Refactoring: Milestone Tracker

**Branch:** `refactor/phase-0-foundation`
**Started:** 2025-11-01
**Target Completion:** 2-3 weeks
**Current Phase:** Phase 0 - Foundation & Safety Net

---

## Overview

This document tracks progress through the intronIC refactoring project. Each phase builds incrementally on the previous, maintaining a working state throughout.

**Key Documents:**
- **REFACTORING_PLAN.md** - Comprehensive architecture and module design
- **REFACTORING_EXECUTION.md** - Detailed chunk-by-chunk execution plan
- **MAJOR_CONCERNS.md** - Critical issues and risk assessment
- **FIELD_APPROACHES.md** - Field-proven improvements to incorporate
- **CLAUDE.md** - Complete codebase documentation

---

## Phase 0: Foundation & Safety Net ⏳

**Status:** IN PROGRESS
**Priority:** CRITICAL - Must complete before any refactoring
**Goal:** Establish regression testing and core abstractions

### Milestones

#### M0.1: Environment Setup ✅ COMPLETE
- [x] Create `intronIC_refactored/` directory structure
- [x] Set up Pixi environment with dependencies
- [x] Fix biogl API compatibility issues
- [x] Create pyproject.toml and pixi.toml
- [x] Verify environment works with test run

**Completed:** 2025-11-01
**Deliverables:**
- Working Pixi environment (Python 3.12.12)
- All dependencies installed and compatible
- Test data symlinked
- Clean project structure

---

#### M0.2: Gold Standard Baseline ✅ COMPLETE
- [x] Run original intronIC on chr19 test data
- [x] Generate complete output set (20,252 introns)
- [x] Validate outputs (BED, meta, introns, scores, plots)
- [x] Document classification results (35 U12 introns detected)
- [x] Archive baseline for comparison

**Completed:** 2025-11-01
**Runtime:** 5.222 minutes
**Deliverables:**
- `gold_standard_outputs/` with all output files
- `VALIDATION_SUMMARY.md` with complete metrics
- Classification baseline: 35 U12s (0.17%), F1=1.0, PR-AUC=1.0

---

#### M0.3: Regression Test Suite ✅ COMPLETE
- [x] Create `tests/test_gold_standard.py`
- [x] Implement 6 core regression tests
- [x] Verify all tests pass (6/6)
- [x] Add pytest to dependencies
- [x] Document test coverage

**Completed:** 2025-11-01
**Test Results:**
```
tests/test_gold_standard.py::test_gold_standard_exists        PASSED
tests/test_gold_standard.py::test_intron_count                PASSED
tests/test_gold_standard.py::test_bed_format                  PASSED
tests/test_gold_standard.py::test_u12_introns_detected        PASSED
tests/test_gold_standard.py::test_classification_scores       PASSED
tests/test_gold_standard.py::test_log_file_completeness       PASSED
```

---

#### M0.4: Git Branch & Documentation ⏳ IN PROGRESS
- [x] Create `refactor/phase-0-foundation` branch
- [ ] Create MILESTONES.md (this file)
- [ ] Commit planning documents
  - [ ] REFACTORING_PLAN.md
  - [ ] REFACTORING_EXECUTION.md
  - [ ] MAJOR_CONCERNS.md
  - [ ] FIELD_APPROACHES.md
  - [ ] CLAUDE.md
- [ ] Commit refactored environment
  - [ ] intronIC_refactored/ structure
  - [ ] pixi.toml and pyproject.toml
  - [ ] tests/test_gold_standard.py
  - [ ] SETUP_COMPLETE.md
  - [ ] VALIDATION_SUMMARY.md
- [ ] Create initial commit: "Phase 0: Environment and baseline established"

**Target:** 2025-11-01
**Status:** Setting up now

---

#### M0.5: Coordinate System Tests 🔜 NEXT
- [ ] Create `tests/unit/test_coordinates.py`
- [ ] Test 0-based ↔ 1-based conversion
- [ ] Test BED format parsing (0-based input)
- [ ] Test internal coordinate handling (1-based)
- [ ] Test BED output conversion (back to 0-based)
- [ ] Test strand handling in coordinates
- [ ] Test validation (start < stop)
- [ ] All coordinate tests pass

**Priority:** 🔴 CRITICAL
**Reason:** Coordinate bugs are #1 risk area (see MAJOR_CONCERNS.md)
**Target:** Next session

---

#### M0.6: Core Coordinate Abstraction 🔜 UPCOMING
- [ ] Create `utils/coordinates.py`
- [ ] Implement `GenomicCoordinate` dataclass
- [ ] Add coordinate system enum (0-based, 1-based)
- [ ] Add validation (start < stop, valid strand)
- [ ] Implement conversion functions
- [ ] Full type hints
- [ ] Docstrings with examples
- [ ] Unit tests pass

**Priority:** 🔴 CRITICAL
**Reason:** Foundation for all coordinate handling
**Target:** After M0.5

---

## Phase 1: Core Infrastructure 📋 PLANNED

**Status:** Not Started
**Prerequisites:** Phase 0 complete
**Duration:** ~1 week (2-3 hours per chunk)

### M1.1: Core Data Models
- [ ] Create `core/models.py` (GenomeFeature hierarchy)
- [ ] Create `core/intron.py` with composition pattern
  - [ ] IntronScores dataclass
  - [ ] IntronSequences dataclass
  - [ ] IntronMetadata dataclass
  - [ ] Intron class using composition
- [ ] Create `utils/sequences.py` (rev_comp, helpers)
- [ ] Unit tests for all models
- [ ] Can create Intron objects successfully

**Priority:** Medium
**Risk:** Low (new code, well-defined)

---

### M1.2: I/O Layer
- [ ] Create `io/genome.py` (GenomeReader with caching)
- [ ] Create `io/parsers.py`
  - [ ] AnnotationParser (GFF3/GTF)
  - [ ] BEDParser (with coordinate conversion)
  - [ ] SequenceParser (.iic files)
- [ ] Create `io/writers.py` (all output formats)
- [ ] Create `io/formats.py` (format specifications)
- [ ] Test with chr19 data
- [ ] Compare outputs to gold standard

**Priority:** High
**Risk:** Medium (coordinate conversions at boundaries)

---

### M1.3: Extraction Pipeline
- [ ] Create `extraction/annotator.py` (hierarchy parsing)
- [ ] Create `extraction/intronator.py` (exon pairs → introns)
- [ ] Create `extraction/sequences.py` (sequence extraction)
- [ ] Create `extraction/filters.py` (deduplication, filtering)
- [ ] Integration test: annotation → introns
- [ ] Verify intron count matches original (20,252)
- [ ] Compare coordinates to gold standard

**Priority:** High
**Risk:** High (complex logic, state management)

---

## Phase 2: Scoring & Classification 📋 PLANNED

**Status:** Not Started
**Prerequisites:** Phase 1 complete
**Duration:** ~1 week

### M2.1: PWM Scoring
- [ ] Create `scoring/pwm.py` (matrix loading)
- [ ] Create `scoring/scorer.py` (seq_score, bp_score)
- [ ] Create `scoring/features.py` (5'ss, BP, 3'ss extraction)
- [ ] Create `scoring/normalization.py` (z-scores)
- [ ] **NEW:** Create `scoring/ppt.py` (PPT detection)
- [ ] Unit tests for all scoring functions
- [ ] Compare scores to gold standard (within tolerance)

**Priority:** High
**Risk:** Low (well-defined math, testable)

---

### M2.2: SVM Classification
- [ ] Create `classification/svm.py` (training, prediction)
- [ ] Create `classification/optimization.py` (hyperparameter tuning)
- [ ] Create `classification/ensemble.py` (multi-model)
- [ ] Integration test: full scoring → classification
- [ ] Verify F1=1.0, PR-AUC=1.0 on training data
- [ ] Compare predictions to gold standard

**Priority:** High
**Risk:** Low (sklearn handles complexity)

---

### M2.3: New Validators & Features ⭐
- [ ] **NEW:** Create `classification/validators.py`
  - [ ] `validate_bp_distance()` (8-30 nt constraint)
  - [ ] `calculate_confidence_level()` (multi-tier)
  - [ ] `calculate_score_differentials()` (U12 vs U2)
- [ ] Add BP distance validation to pipeline
- [ ] Add confidence levels to output
- [ ] Add PPT strength to features
- [ ] Test with known U12 introns
- [ ] Document new features

**Priority:** High (field-proven improvements)
**Risk:** Low (additions, not modifications)

---

### M2.4: Advanced Classification
- [ ] Create `classification/recursive.py` (recursive training)
- [ ] Refactor `u12_correction()` to return new object
- [ ] Add validation before accepting corrections
- [ ] Test recursive training
- [ ] Compare to gold standard

**Priority:** Medium
**Risk:** High (complex state management)

---

## Phase 3: CLI & Configuration 📋 PLANNED

**Status:** Not Started
**Prerequisites:** Phase 2 complete
**Duration:** ~2-3 days

### M3.1: Configuration System
- [ ] Create `config/defaults.py` (dataclasses)
- [ ] Create `config/schemas.py` (validation)
- [ ] Add YAML config file support
- [ ] Test config loading and validation

**Priority:** Medium
**Risk:** Low

---

### M3.2: Command-Line Interface
- [ ] Create `cli/parser.py` (argparse wrapper)
- [ ] Create `cli/main.py` (entry point)
- [ ] Create `__main__.py` for `python -m intronIC_refactored`
- [ ] Preserve original CLI compatibility
- [ ] End-to-end test from command line
- [ ] Compare outputs to gold standard

**Priority:** High
**Risk:** Low

---

## Phase 4: Analysis & Visualization 📋 PLANNED

**Status:** Not Started
**Prerequisites:** Phase 3 complete
**Duration:** ~2-3 days

### M4.1: Analysis Tools
- [ ] Create `analysis/plotting.py` (scatter, PR-AUC, etc.)
- [ ] Create `analysis/statistics.py` (summary stats)
- [ ] Generate all plots
- [ ] Compare plots to gold standard visually

**Priority:** Medium
**Risk:** Low

---

## Phase 5: Testing & Documentation 📋 PLANNED

**Status:** Not Started
**Prerequisites:** Phases 1-4 complete
**Duration:** ~3-4 days

### M5.1: Comprehensive Test Suite
- [ ] Unit tests for all modules (>85% coverage)
- [ ] Integration tests (full pipeline)
- [ ] Property-based tests (hypothesis)
- [ ] Performance benchmarks
- [ ] All tests pass

**Priority:** Critical
**Risk:** N/A

---

### M5.2: Documentation
- [ ] API documentation (all modules)
- [ ] Tutorial (getting started)
- [ ] Migration guide (v1 → v2)
- [ ] Algorithm documentation
- [ ] Contributing guide
- [ ] Type check with mypy (100% coverage)
- [ ] Format with black

**Priority:** High
**Risk:** Low

---

### M5.3: CI/CD Setup
- [ ] Create `.github/workflows/tests.yml`
- [ ] Create `.github/workflows/lint.yml`
- [ ] Test across Python versions (3.8-3.12)
- [ ] Automated testing on push/PR
- [ ] Coverage reporting

**Priority:** Medium
**Risk:** Low

---

## Phase 6: Release Preparation 📋 PLANNED

**Status:** Not Started
**Prerequisites:** Phase 5 complete
**Duration:** ~2-3 days

### M6.1: Final Validation
- [ ] Run full pipeline on chr19 (final comparison)
- [ ] Run on diverse test datasets
- [ ] Performance benchmarking (vs original)
- [ ] Memory profiling
- [ ] Fix any remaining issues

**Priority:** Critical
**Risk:** N/A

---

### M6.2: Release
- [ ] Version bump to 2.0.0
- [ ] Create CHANGELOG.md
- [ ] Tag release
- [ ] Create GitHub release
- [ ] Update README with v2 info
- [ ] Announce to users

**Priority:** High
**Risk:** N/A

---

## Success Criteria

### Phase 0 (Foundation) ✅
- [x] Gold standard outputs generated
- [x] Regression tests pass (6/6)
- [x] Environment working
- [ ] Coordinate tests written
- [ ] GenomicCoordinate abstraction created
- [ ] Git branch established with planning docs

### Phase 1-2 (Core + Scoring)
- [ ] Full pipeline functional
- [ ] All outputs match gold standard
- [ ] New features (BP distance, PPT, confidence) working
- [ ] >85% test coverage for new code
- [ ] No performance regression

### Phase 3-4 (CLI + Viz)
- [ ] Command-line interface working
- [ ] All plots generated correctly
- [ ] Original CLI arguments supported

### Phase 5-6 (Test + Release)
- [ ] >90% overall test coverage
- [ ] All documentation complete
- [ ] CI/CD functional
- [ ] Ready for public use

---

## Risk Register

### Critical Risks (Must Mitigate)

#### R1: Coordinate Conversion Errors 🔴
**Likelihood:** High (historical bugs)
**Impact:** Critical (wrong results)
**Mitigation:**
- Extensive tests for coordinate handling (M0.5)
- Explicit GenomicCoordinate abstraction (M0.6)
- Validation at every conversion
- Review all BED I/O carefully

**Status:** Mitigating with M0.5 and M0.6

---

#### R2: Regression in Classification Accuracy 🔴
**Likelihood:** Medium
**Impact:** Critical (defeats purpose)
**Mitigation:**
- Comprehensive regression tests
- Compare every output to gold standard
- Validate on multiple datasets
- Track F1, PR-AUC metrics

**Status:** Protected by M0.2 and M0.3

---

#### R3: Performance Regression 🟡
**Likelihood:** Medium
**Impact:** High (user experience)
**Mitigation:**
- Benchmark each phase
- Profile bottlenecks
- Optimize before release
- Acceptance: <10% slowdown

**Status:** Will benchmark in Phase 6

---

### Medium Risks

#### R4: Scope Creep 🟡
**Likelihood:** High
**Impact:** Medium (timeline)
**Mitigation:**
- Strict adherence to roadmap
- Defer nice-to-haves to v2.1
- Focus on MVP first

**Status:** Controlled with this milestone document

---

#### R5: State Management Bugs 🟡
**Likelihood:** Medium
**Impact:** Medium
**Mitigation:**
- Use composition (IntronScores, etc.)
- Reduce mutability
- Clear state transitions
- Extensive testing

**Status:** Addressed in architecture (see REFACTORING_PLAN.md)

---

## Timeline

### Optimistic (Full-Time Focus)
- **Phase 0:** 2-3 days ← WE ARE HERE
- **Phase 1:** 4-5 days
- **Phase 2:** 4-5 days
- **Phase 3:** 2-3 days
- **Phase 4:** 2-3 days
- **Phase 5:** 3-4 days
- **Phase 6:** 2-3 days
- **Total:** 2-3 weeks

### Realistic (Part-Time)
- **Phase 0:** 1 week ← WE ARE HERE
- **Phase 1-2:** 2 weeks
- **Phase 3-4:** 1 week
- **Phase 5-6:** 1 week
- **Total:** 5-6 weeks

### Current Status
- **Started:** 2025-11-01
- **Current Phase:** Phase 0 (Foundation)
- **Progress:** M0.1-M0.3 complete (60% of Phase 0)
- **Next:** M0.4-M0.6 (coordinate handling)

---

## Notes

### Key Decisions Made
1. **Pixi for environment management** - Better than conda/venv for reproducibility
2. **Composition over inheritance** - IntronScores, IntronSequences, IntronMetadata
3. **Explicit coordinate systems** - GenomicCoordinate with validation
4. **Test-first approach** - Regression tests before refactoring
5. **Incremental delivery** - Each phase produces working code

### Lessons Learned
1. Coordinate bugs are the #1 risk - extensive testing required
2. Gold standard baseline is essential - can't refactor without it
3. API compatibility matters - biogl version caused issues
4. Test environment early - catch dependency problems

### Questions for Future Sessions
1. Should we use dataclasses with `frozen=True` for immutability?
2. Type hints: strict (mypy) or lenient?
3. Performance target: match original or 10% speedup?
4. Support Python 3.8+ or only 3.10+?

---

## Quick Reference

### Commands

```bash
# Run tests
pixi run python -m pytest tests/ -v

# Run gold standard test
pixi run test-gold

# Check current branch
git branch

# View milestone progress
cat MILESTONES.md | grep -E "(\[x\]|\[ \])" | head -20

# Compare outputs
diff homo_sapiens_test.bed.iic gold_standard_baseline/homo_sapiens_gold.bed.iic
```

### Files to Track
- This file: `MILESTONES.md`
- Planning: `REFACTORING_PLAN.md`, `REFACTORING_EXECUTION.md`
- Issues: `MAJOR_CONCERNS.md`
- Validation: `VALIDATION_SUMMARY.md`
- Setup: `SETUP_COMPLETE.md`

---

**Last Updated:** 2025-11-01
**Next Review:** After M0.6 complete
**Maintained By:** Project Team
