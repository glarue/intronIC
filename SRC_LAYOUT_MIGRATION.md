# src/ Layout Migration Plan

**Goal**: Migrate intronIC from flat layout to src layout for better packaging isolation.

**Status**: Planning Phase

---

## Migration Overview

### Current Structure
```
intronIC/  (repo root)
  cli/
  extraction/
  file_io/
  scoring/
  classification/
  visualization/
  utils/
  core/
  data/
  tests/
  scripts/
  pyproject.toml
```

### Target Structure
```
intronIC/  (repo root)
  src/
    intronIC/
      cli/
      extraction/
      file_io/
      scoring/
      classification/
      visualization/
      utils/
      core/
      data/
  tests/
  scripts/
  pyproject.toml
```

### Import Changes
All imports change from:
```python
from cli.main import run
from extraction.annotator import Annotator
```
to:
```python
from intronIC.cli.main import run
from intronIC.extraction.annotator import Annotator
```

---

## Phase 1: Configuration Files

### pyproject.toml
- [ ] Update `[project.scripts]` entry point: `cli.main:main` → `intronIC.cli.main:main`
- [ ] Update `[tool.hatch.build.targets.wheel]` packages list
- [ ] Add `where = ["src"]` to package discovery
- [ ] Update `[tool.pytest.ini_options]` pythonpath

---

## Phase 2: Source Code Files

### Core Module (3 files)
- [ ] `core/__init__.py`
- [ ] `core/intron.py` - Has imports
- [ ] `core/models.py` - Has imports

### CLI Module (8 files)
- [ ] `cli/__init__.py`
- [ ] `cli/args.py`
- [ ] `cli/args_old.py`
- [ ] `cli/config.py` - Uses `__file__`
- [ ] `cli/config_loader.py`
- [ ] `cli/main.py` - Uses `__file__`, has many imports
- [ ] `cli/messenger.py`
- [ ] `cli/progress.py`

### Extraction Module (6 files)
- [ ] `extraction/__init__.py`
- [ ] `extraction/annotator.py` - Has imports
- [ ] `extraction/boundary_correction.py`
- [ ] `extraction/filters.py` - Has imports
- [ ] `extraction/intronator.py`
- [ ] `extraction/sequences.py` - Has imports

### File I/O Module (6 files)
- [ ] `file_io/__init__.py`
- [ ] `file_io/annotation_index.py`
- [ ] `file_io/genome.py` - Has imports
- [ ] `file_io/indexed_genome.py`
- [ ] `file_io/parsers.py`
- [ ] `file_io/writers.py`

### Scoring Module (5 files)
- [ ] `scoring/__init__.py`
- [ ] `scoring/branch_point.py` - Has imports
- [ ] `scoring/normalizer.py`
- [ ] `scoring/pwm.py`
- [ ] `scoring/scorer.py` - Has imports

### Classification Module (15 files)
- [ ] `classification/__init__.py`
- [ ] `classification/calibration_validator.py`
- [ ] `classification/classifier.py` - Has imports
- [ ] `classification/clipping.py`
- [ ] `classification/config_loader.py` - Uses `__file__`, has imports
- [ ] `classification/fast_test_config.py`
- [ ] `classification/model_inspector.py` - Has imports
- [ ] `classification/nested_cv.py` - Has imports
- [ ] `classification/optimizer.py` - Has imports
- [ ] `classification/predictor.py` - Has imports
- [ ] `classification/progress_tracker.py`
- [ ] `classification/saturating.py`
- [ ] `classification/split_eval.py` - Has imports
- [ ] `classification/trainer.py` - Has imports
- [ ] `classification/transformers.py`

### Visualization Module (2 files)
- [ ] `visualization/__init__.py`
- [ ] `visualization/plots.py`

### Utils Module (7 files)
- [ ] `utils/__init__.py`
- [ ] `utils/convert_pwm_format.py`
- [ ] `utils/coordinates.py`
- [ ] `utils/logging_utils.py`
- [ ] `utils/metadata.py` - Uses `__file__`
- [ ] `utils/model_io.py`
- [ ] `utils/sequences.py`

**Source Subtotal: 58 files**

---

## Phase 3: Test Files

### Unit Tests - Core (3 files)
- [ ] `tests/unit/test_intron.py` - Uses `__file__`, has imports
- [ ] `tests/unit/test_models.py` - Has imports
- [ ] `tests/unit/test_sequences.py` - Has imports

### Unit Tests - File I/O (2 files)
- [ ] `tests/unit/test_genome.py` - Uses `__file__`, has imports
- [ ] `tests/unit/test_parsers.py` - Has imports
- [ ] `tests/unit/test_writers.py` - Has imports

### Unit Tests - Extraction (2 files)
- [ ] `tests/unit/test_coordinates.py` - Has imports
- [ ] `tests/unit/test_filters.py` - Has imports

### Unit Tests - Scoring (16 files)
- [ ] `tests/unit/test_scoring/__init__.py`
- [ ] `tests/unit/test_scoring/test_branch_point.py` - Has imports
- [ ] `tests/unit/test_scoring/test_correct_matrices.py` - Uses `__file__`, has imports
- [ ] `tests/unit/test_scoring/test_log_ratio.py` - Uses `__file__`, has imports
- [ ] `tests/unit/test_scoring/test_log_ratio_deterministic.py` - Has imports
- [ ] `tests/unit/test_scoring/test_matrix_conversion.py` - Uses `__file__`, has imports
- [ ] `tests/unit/test_scoring/test_matrix_selection.py` - Uses `__file__`, has imports
- [ ] `tests/unit/test_scoring/test_new_scaling_architecture.py` - Has imports
- [ ] `tests/unit/test_scoring/test_normalizer.py` - Has imports
- [ ] `tests/unit/test_scoring/test_pwm.py` - Has imports
- [ ] `tests/unit/test_scoring/test_pwm_debug.py` - Uses `__file__`, has imports
- [ ] `tests/unit/test_scoring/test_pwm_fix.py` - Uses `__file__`, has imports
- [ ] `tests/unit/test_scoring/test_pwm_format_equivalence.py` - Has imports
- [ ] `tests/unit/test_scoring/test_pwm_indexing.py` - Has imports
- [ ] `tests/unit/test_scoring/test_scorer.py` - Has imports

### Unit Tests - Classification (5 files)
- [ ] `tests/unit/test_classification/__init__.py`
- [ ] `tests/unit/test_classification/test_classifier.py` - Has imports
- [ ] `tests/unit/test_classification/test_optimizer.py` - Has imports
- [ ] `tests/unit/test_classification/test_predictor.py` - Has imports
- [ ] `tests/unit/test_classification/test_trainer.py` - Has imports

### Unit Tests - File I/O (1 file)
- [ ] `tests/unit/test_file_io/test_writer_errors.py` - Has imports

### Unit Tests - Misc (3 files)
- [ ] `tests/unit/test_edge_cases.py` - Has imports
- [ ] `tests/unit/test_error_handling.py` - Has imports

### Integration Tests (4 files)
- [ ] `tests/integration/__init__.py`
- [ ] `tests/integration/test_classification_pipeline.py` - Uses `__file__`, has imports
- [ ] `tests/integration/test_cli.py` - Has imports
- [ ] `tests/integration/test_extraction_pipeline.py` - Has imports
- [ ] `tests/integration/test_parser_writer_pipeline.py` - Has imports

### Test Support (2 files)
- [ ] `tests/conftest.py` - Uses `__file__`
- [ ] `tests/test_gold_standard.py` - Has imports

**Test Subtotal: 39 files**

---

## Phase 4: Scripts and Utilities

### Scripts (10 files)
- [ ] `scripts/__init__.py`
- [ ] `scripts/analyze_false_positives.py`
- [ ] `scripts/analyze_fps.py`
- [ ] `scripts/check_version.py` - Uses `__file__`
- [ ] `scripts/compare_true_vs_false_u12s.py`
- [ ] `scripts/debug_scores.py` - Uses `__file__`, has imports
- [ ] `scripts/test_bothends_smoke.py` - Uses `__file__`, has imports
- [ ] `scripts/test_chr19_fast.py` - Uses `__file__`
- [ ] `scripts/test_corrected_architecture.py` - Has imports
- [ ] `scripts/test_transformer.py` - Has imports

### Root Scripts (1 file)
- [ ] `debug_cycles.py` - Uses `__file__`, has imports

**Scripts Subtotal: 11 files**

---

## Phase 5: Special Cases Requiring Extra Attention

### Files Using `__file__` (31 files)
These need path adjustments after move to src/intronIC/:

**Source files:**
1. `cli/config.py` - Likely loading config files
2. `cli/main.py` - Loading default PWM/model files from data/
3. `classification/config_loader.py` - Loading config files
4. `utils/metadata.py` - Package metadata loading

**Test files:**
5. `tests/conftest.py` - Test fixtures and paths
6. `tests/unit/test_genome.py` - Test data paths
7. `tests/unit/test_intron.py` - Test data paths
8. `tests/unit/test_scoring/test_correct_matrices.py` - Test data paths
9. `tests/unit/test_scoring/test_log_ratio.py` - Test data paths
10. `tests/unit/test_scoring/test_matrix_conversion.py` - Test data paths
11. `tests/unit/test_scoring/test_matrix_selection.py` - Test data paths
12. `tests/unit/test_scoring/test_pwm_debug.py` - Test data paths
13. `tests/unit/test_scoring/test_pwm_fix.py` - Test data paths
14. `tests/integration/test_classification_pipeline.py` - Test data paths

**Scripts:**
15. `debug_cycles.py` - Path resolution
16. `scripts/check_version.py` - Version file path
17. `scripts/debug_scores.py` - Data file paths
18. `scripts/test_bothends_smoke.py` - Test data paths
19. `scripts/test_chr19_fast.py` - Test data paths

### Common `__file__` Patterns to Update:

**Pattern 1: Data files relative to source**
```python
# Before:
data_dir = Path(__file__).parent.parent / "data"

# After:
data_dir = Path(__file__).parent.parent.parent / "data"  # Add extra .parent
# OR use package resources (more robust)
```

**Pattern 2: Test fixtures**
```python
# Before:
fixtures_dir = Path(__file__).parent / "fixtures"

# After (no change - tests/ stays at root):
fixtures_dir = Path(__file__).parent / "fixtures"
```

---

## Validation Checklist

After migration:

### Build & Install
- [ ] `pixi clean cache --pypi` - Clear package cache
- [ ] `pixi install` - Reinstall with new structure
- [ ] `pixi run intronIC --help` - Verify CLI works
- [ ] `pixi run intronIC --version` - Verify version detection

### Import Tests
- [ ] `pixi run python -c "from intronIC.cli.main import main"` - Test new imports
- [ ] `pixi run python -c "from intronIC.scoring.pwm import PWMLoader"` - Test module imports
- [ ] `pixi run python -c "from intronIC.core.intron import Intron"` - Test core imports

### Unit Tests
- [ ] `pixi run pytest tests/unit/test_scoring/ -v` - Scoring tests
- [ ] `pixi run pytest tests/unit/test_classification/ -v` - Classification tests
- [ ] `pixi run pytest tests/unit/ -v` - All unit tests

### Integration Tests
- [ ] `pixi run pytest tests/integration/ -v` - Integration tests
- [ ] `pixi run pytest tests/test_gold_standard.py -v` - Gold standard test

### Full Pipeline Tests
- [ ] Run small dataset (C. elegans) - Verify end-to-end classification
- [ ] Run training pipeline - Verify model training works
- [ ] Verify output files generated correctly

---

## Execution Strategy

1. **Create feature branch**:
   ```bash
   git checkout -b refactor/src-layout
   ```

2. **Create directory structure**:
   ```bash
   mkdir -p src/intronIC
   ```

3. **Move source packages** (keeps git history):
   ```bash
   git mv cli src/intronIC/
   git mv core src/intronIC/
   git mv extraction src/intronIC/
   git mv file_io src/intronIC/
   git mv scoring src/intronIC/
   git mv classification src/intronIC/
   git mv visualization src/intronIC/
   git mv utils src/intronIC/
   git mv data src/intronIC/
   ```

4. **Update pyproject.toml** (Phase 1)

5. **Automated import updates**:
   ```bash
   # Use sed/awk or Python script to update all imports
   # Example for manual review first:
   grep -r "^from \(cli\|extraction\|file_io\|scoring\|classification\|visualization\|utils\|core\)" \
     src/ tests/ scripts/ *.py
   ```

6. **Fix `__file__` paths** (Phase 5 - manual review)

7. **Reinstall and test**:
   ```bash
   pixi clean cache --pypi
   pixi install
   pixi run pytest tests/
   ```

8. **Commit**:
   ```bash
   git add -A
   git commit -m "Refactor: Migrate to src/ layout for better packaging isolation

   - Move all source packages to src/intronIC/
   - Update all imports to use intronIC.* namespace
   - Update pyproject.toml for src layout
   - Fix __file__ paths for new structure
   - All tests passing
   "
   ```

---

## Rollback Plan

If issues arise:
```bash
git checkout main  # Return to previous branch
git branch -D refactor/src-layout  # Delete migration branch
```

All changes isolated to feature branch until validated.

---

## File Count Summary

- **Source files**: 58
- **Test files**: 39
- **Scripts**: 11
- **Config files**: 1 (pyproject.toml)
- **Special attention** (`__file__`): 31 files

**Total files to update**: ~109 files

---

## Notes

- Tests stay at root level (`tests/`) - no need to move
- Scripts stay at root level (`scripts/`) - import updates only
- Data moves to `src/intronIC/data/` (with source)
- Archive directory unaffected (not imported)
