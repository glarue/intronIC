# intronIC v2.0.0 Migration Plan

**Date:** 2025-11-11
**Goal:** Replace v1.5.1 implementation with refactored v2.0.0 codebase
**Reference:** See `MIGRATION_MANIFEST.md` for detailed file listing

---

## Phase 1: Pre-Migration Preparation

### 1.1 Create Git Tags & Branches for v1.5.1

```bash
# Tag the current state as v1.5.1
git tag -a v1.5.1 -m "Final release of monolithic implementation"

# Create archive branch for v1.5.1
git branch v1-archive
```

**Purpose:** Preserve the original implementation for reference and potential rollback.

---

### 1.2 Create Archive Directory Structure

```bash
mkdir -p archive/v1.5.1
mkdir -p archive/investigation
mkdir -p archive/refactoring
mkdir -p archive/test_outputs
mkdir -p archive/comparison_test
mkdir -p archive/classification_backup
mkdir -p docs
```

**Purpose:** Organize files we want to keep but move out of active development.

---

## Phase 2: Archive Old Files

### 2.1 Archive v1.5.1 Implementation

```bash
# Original implementation
git mv intronIC/intronIC.py archive/v1.5.1/intronIC/
git mv intronIC/pwms_from_seqs.py archive/v1.5.1/intronIC/
git mv intronIC/matrices_from_seqs.py archive/v1.5.1/intronIC/
git mv intronIC/_version.py archive/v1.5.1/intronIC/
git mv intronIC/__init__.py archive/v1.5.1/intronIC/

# Old packaging
git mv setup.py archive/v1.5.1/
git mv setup.cfg archive/v1.5.1/
git mv versioneer.py archive/v1.5.1/
git mv MANIFEST.in archive/v1.5.1/  # Will be replaced with new version

# Old test_data location (will be moved to data/test_data/ later)
# Keep for now, will remove after copying to new location in Phase 3
```

---

### 2.2 Archive Investigation Files

```bash
# Analysis scripts
git mv analyze_*.py archive/investigation/
git mv compare_*.py archive/investigation/
git mv debug_*.py archive/investigation/
git mv diagnose_*.py archive/investigation/
git mv investigate_*.py archive/investigation/
git mv profile_*.py archive/investigation/
git mv test_*.py archive/investigation/  # Root-level test scripts only

# Investigation documentation
git mv *_INVESTIGATION.md archive/investigation/
git mv *_ANALYSIS.md archive/investigation/
git mv *_AUDIT*.md archive/investigation/
git mv *_DISCREPANCY*.md archive/investigation/
git mv *_RESULTS.md archive/investigation/
git mv *_SUMMARY.md archive/investigation/  # Session summaries
git mv NEXT_STEPS*.md archive/investigation/
git mv OPTIMIZATION_CALL_FLOW.md archive/investigation/
git mv FILESYSTEM_REFERENCE.md archive/investigation/

# Log and profile files
git mv *.log archive/investigation/
git mv *.pstats archive/investigation/
git mv *.txt archive/investigation/ # Output text files

# Reference data copies (root level)
git mv u2_reference_small.introns.iic.gz archive/investigation/
git mv u12_reference_small.introns.iic.gz archive/investigation/
```

---

### 2.3 Archive Refactoring Session Files

```bash
cd intronIC_refactored/

# Session documentation
git mv SESSION_*.md ../archive/refactoring/
git mv SETUP*.md ../archive/refactoring/
git mv VALIDATION*.md ../archive/refactoring/
git mv PROGRESS*.md ../archive/refactoring/
git mv SCORING_ANALYSIS.md ../archive/refactoring/
git mv U12_CORRECTION_*.md ../archive/refactoring/
git mv PORTING_CHECKLIST.md ../archive/refactoring/
git mv REVIEW_FINDINGS.md ../archive/refactoring/
git mv TODO_FROM_REVIEW.md ../archive/refactoring/
git mv DYNAMIC_TAGGING_IMPLEMENTATION.md ../archive/refactoring/
git mv OUTPUT_FORMAT_REFACTOR_PLAN.md ../archive/refactoring/
git mv NEXT_SESSION.md ../archive/refactoring/

# Test outputs
git mv *.bed.iic ../archive/test_outputs/ 2>/dev/null || true
git mv *.introns.iic ../archive/test_outputs/ 2>/dev/null || true
git mv *.meta.iic ../archive/test_outputs/ 2>/dev/null || true
git mv *.score_info.iic ../archive/test_outputs/ 2>/dev/null || true
git mv *.plot.*.png ../archive/test_outputs/ 2>/dev/null || true
git mv *.log ../archive/test_outputs/ 2>/dev/null || true
git mv *.metrics.json ../archive/test_outputs/ 2>/dev/null || true
git mv *.tar.gz ../archive/test_outputs/ 2>/dev/null || true
git mv drosophila/ ../archive/test_outputs/ 2>/dev/null || true

# Keep homo_sapiens.model.pkl - we need this for default pretrained model
# (will be moved to data/ in Phase 3)

# Backup directories
git mv classification_backup/ ../archive/classification_backup/

cd ..
```

---

### 2.4 Archive Comparison Test Directory

```bash
git mv comparison_test/ archive/comparison_test/
```

---

### 2.5 Clean Up Build Artifacts (Delete, Don't Archive)

```bash
cd intronIC_refactored/

# Remove build artifacts
rm -rf __pycache__/
rm -rf .pytest_cache/
rm -rf .pixi/  # Will be regenerated
rm -f .gitignore  # Redundant with root
rm -f test_run.log  # Empty file

cd ..
```

---

## Phase 3: Prepare Data Directory

### 3.1 Create Clean Data Directory Structure

```bash
# The intronIC/data/ directory should already exist from v1
# We'll update it with files from both old and new versions
# Nest test_data under data/ for better organization

mkdir -p intronIC/data
mkdir -p intronIC/data/test_data
```

---

### 3.2 Copy Production Data Files

```bash
# PWM matrices (from original - these are the authoritative versions)
cp intronIC/data/scoring_matrices.fasta.iic intronIC/data/
cp intronIC/data/scoring_matrices.2+_pubs.fasta.iic intronIC/data/
cp intronIC/data/scoring_matrices.2+.pub_required.fasta.iic intronIC/data/
cp intronIC/data/u2.conserved_empirical_bp_pwm.iic intronIC/data/
cp intronIC/data/empirical.u12s.2+.pub_required.matrices intronIC/data/

# Reference sequences (from original)
cp intronIC/data/u12_reference.introns.iic.gz intronIC/data/
cp intronIC/data/u2_reference.introns.iic.gz intronIC/data/

# Test data (from original, now nested under data/)
cp intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz intronIC/data/test_data/
cp intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz intronIC/data/test_data/

# Default pretrained model (from refactored outputs)
cp intronIC_refactored/homo_sapiens.model.pkl intronIC/data/default_pretrained.model.pkl
```

**Note:** These are `cp` commands (not `git mv`) because the data directory structure stays as-is - we're just ensuring all needed files are present.

---

### 3.3 Remove Test-Only Data Files

```bash
cd intronIC/data/

# Remove test/debug reference files
rm -f u12_reference_small.introns.iic.gz
rm -f u2_reference_small.introns.iic.gz
rm -f u12_reference_micro.introns.iic.gz
rm -f u2_reference_micro.introns.iic.gz
rm -f u12_reference_tiny.introns.iic.gz
rm -f u2_reference_tiny.introns.iic.gz

cd ../..

# Remove old test_data directory (now nested under data/)
rm -rf intronIC/test_data/
```

---

## Phase 4: Move Refactored Code to Root

### 4.1 Move Module Directories

```bash
cd intronIC_refactored/

# Move all module directories to root
git mv cli/ ../
git mv core/ ../
git mv extraction/ ../
git mv file_io/ ../
git mv scoring/ ../
git mv classification/ ../
git mv utils/ ../
git mv visualization/ ../
git mv tests/ ../
git mv scripts/ ../
```

---

### 4.2 Move Configuration Files

```bash
# Package configuration
git mv pyproject.toml ../
git mv pixi.toml ../
git mv pixi.lock ../
git mv MANIFEST.in ../
git mv __init__.py ../
git mv __main__.py ../

# Installation helpers
git mv install.sh ../
git mv install.bat ../

# Testing configuration
git mv pytest.ini ../
```

---

### 4.3 Move Documentation

```bash
# User-facing documentation (to root)
git mv README.md ../
git mv INSTALL.md ../
git mv QUICKSTART.md ../
# LICENSE already exists in root - skip

# Developer documentation (to docs/)
git mv PACKAGING_STRATEGY.md ../docs/
git mv VERSIONING_STRATEGY.md ../docs/
git mv VERSIONING_COMPARISON.md ../docs/
git mv RELEASE_PROCESS.md ../docs/
git mv REFACTOR_SUMMARY.md ../docs/
git mv TEST_IMPROVEMENT_PLAN.md ../docs/
git mv IMPLEMENTATION_PLAN.md ../docs/
git mv PERFORMANCE_IMPROVEMENTS.md ../docs/
git mv SETUP.md ../docs/

cd ..
```

---

### 4.4 Handle Special Files

```bash
# The test_data symlink in intronIC_refactored/ can be removed
# (it pointed to ../intronIC/test_data which we're keeping)
rm intronIC_refactored/test_data

# Remove the now-empty intronIC subdirectory in intronIC_refactored
# (the data files were already in the root-level intronIC/data/)
# This subdirectory just held the symlink to test_data
rmdir intronIC_refactored/intronIC 2>/dev/null || true

# Remove examples directory if it exists
rm -rf intronIC_refactored/examples
```

---

### 4.5 Remove Empty Refactored Directory

```bash
# At this point, intronIC_refactored/ should be empty (or nearly so)
# Verify it's empty:
ls -la intronIC_refactored/

# If empty, remove it:
rmdir intronIC_refactored/
```

---

## Phase 5: Code Modifications

### 5.1 Add Default Pretrained Model Support

**File:** `cli/args.py`

Add default for `--pretrained_model`:

```python
training_group.add_argument(
    '--pretrained_model',
    type=Path,
    default=None,  # Will be set in main.py if None
    help='Path to pretrained model file (.model.pkl). Skips training and uses saved model for classification'
)
```

**File:** `cli/main.py`

Add logic to use default pretrained model:

```python
def get_default_pretrained_model_path():
    """Get path to default pretrained model in data directory."""
    default_path = Path(__file__).parent.parent / "intronIC" / "data" / "default_pretrained.model.pkl"
    if default_path.exists():
        return default_path
    return None

# In main() function, before model training section:
# Check if user wants to use pretrained model
pretrained_model_path = config.training.pretrained_model
if pretrained_model_path is None and config.training.use_pretrained:
    # Use default if available
    pretrained_model_path = get_default_pretrained_model_path()
    if pretrained_model_path:
        logger.info(f"Using default pretrained model: {pretrained_model_path}")
```

Add CLI argument for automatic pretrained model usage:

```python
# In args.py
training_group.add_argument(
    '--use_pretrained',
    action='store_true',
    help='Use default pretrained model from data directory (skip training)'
)
```

---

### 5.2 Update Version References

Verify version is set to `2.0.0` in:
- `pyproject.toml`: `version = "2.0.0"`
- `pixi.toml`: `version = "2.0.0"`
- `__init__.py`: Uses `importlib.metadata.version()` (dynamic)

---

### 5.3 Update README for v2.0.0

**File:** `README.md`

Update version badges, installation instructions, and changelog to reflect v2.0.0.

---

## Phase 6: Verification

### 6.1 Verify Directory Structure

```bash
# Check that root-level structure looks correct
ls -la

# Should see:
# - Module directories: cli/, core/, extraction/, file_io/, scoring/, classification/, utils/, visualization/, tests/
# - Data: intronIC/data/ (with test_data/ nested inside)
# - Config: pyproject.toml, pixi.toml, pytest.ini
# - Docs: README.md, INSTALL.md, QUICKSTART.md, docs/
# - Archive: archive/
```

---

### 6.2 Verify Data Files

```bash
# Check data directory has all required files
ls -lh intronIC/data/

# Should include:
# - 5 PWM matrix files
# - 2 reference sequence files (.iic.gz)
# - 1 pretrained model (default_pretrained.model.pkl)
# - test_data/ subdirectory

ls -lh intronIC/data/test_data/

# Should include:
# - Homo_sapiens.Chr19.Ensembl_91.fa.gz
# - Homo_sapiens.Chr19.Ensembl_91.gff3.gz
```

---

### 6.3 Regenerate Pixi Environment

```bash
# Remove old environment
rm -rf .pixi/

# Regenerate from pixi.toml
pixi install

# Verify installation
pixi run intronIC --version
# Should show: intronIC 2.0.0
```

---

### 6.4 Run Tests

```bash
# Run full test suite
pixi run pytest tests/ -v

# Run quick integration test with Chr19 data
pixi run intronIC \
  -g intronIC/data/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
  -a intronIC/data/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz \
  -n test_migration \
  --eval_mode split

# Verify output files were created
ls -lh test_migration.*
```

---

### 6.5 Test Pretrained Model

```bash
# Test using default pretrained model
pixi run intronIC \
  -g intronIC/data/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
  -a intronIC/data/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz \
  -n test_pretrained \
  --use_pretrained

# Should skip training and use default_pretrained.model.pkl
```

---

## Phase 7: Git Commit

### 7.1 Stage All Changes

```bash
# Add all changes
git add -A

# Review what's being committed
git status
git diff --cached --stat
```

---

### 7.2 Create Migration Commit

```bash
git commit -m "Migrate to v2.0.0 refactored architecture

BREAKING CHANGE: Complete architectural refactor from monolithic to modular design

Major Changes:
- Replace monolithic intronIC.py (6,093 lines) with modular architecture
- 9 new modules: cli, core, extraction, file_io, scoring, classification, utils, visualization, tests
- Modern packaging: pyproject.toml + pixi.toml (replace setup.py + versioneer)
- Comprehensive test suite: 76+ tests with pytest
- Add default pretrained model support
- Improved documentation: INSTALL.md, QUICKSTART.md, docs/

Architecture:
- Separation of concerns: each module has single responsibility
- Immutable data structures using dataclasses
- Type hints throughout
- Zero-Anchored Robust (ZAR) normalization for LLR features
- Explicit ML integrity safeguards (prevent data leakage)

Data:
- Consolidated data directory: intronIC/data/
- Include default pretrained model: default_pretrained.model.pkl
- Preserve all production PWM matrices and reference sequences
- Remove test-only data files (small/micro/tiny references)

Archived:
- v1.5.1 implementation → archive/v1.5.1/
- Investigation files → archive/investigation/
- Refactoring docs → archive/refactoring/
- Test outputs → archive/test_outputs/

Version: 2.0.0
Previous: v1.5.1 (tagged and archived)"
```

---

### 7.3 Tag New Version

```bash
git tag -a v2.0.0 -m "Release v2.0.0: Modular refactored architecture"
```

---

## Phase 8: Cleanup and Documentation

### 8.1 Update GitHub Branch

```bash
# Push to current branch
git push origin claude/investigate-ref-performance-011CUt9YD4yNb4AzXF2TKAU9

# Create PR to merge into main
gh pr create --title "Release v2.0.0: Modular Architecture Refactor" \
  --body "Complete migration to refactored v2.0.0 architecture. See MIGRATION_PLAN.md for details."
```

---

### 8.2 Create Release Notes

Create `RELEASE_NOTES_v2.0.0.md` summarizing:
- Major architectural changes
- Breaking changes from v1.5.1
- Migration guide for users
- New features and improvements

---

### 8.3 Update Documentation

Ensure all documentation reflects new structure:
- Update links in README.md
- Update QUICKSTART.md examples
- Update INSTALL.md for new package structure

---

## Phase 9: Optional PyPI Publication

See `docs/RELEASE_PROCESS.md` for PyPI publication steps:

```bash
# Build distribution
python -m build

# Upload to PyPI
python -m twine upload dist/*
```

---

## Rollback Plan

If migration fails or issues are discovered:

```bash
# Revert to v1.5.1
git reset --hard v1.5.1

# Or use archive branch
git checkout v1-archive
```

---

## Post-Migration Checklist

- [ ] Phase 1: Git tags and branches created
- [ ] Phase 2: Old files archived
- [ ] Phase 3: Data directory cleaned and organized
- [ ] Phase 4: Refactored code moved to root
- [ ] Phase 5: Code modifications completed
- [ ] Phase 6: All verification tests passed
- [ ] Phase 7: Changes committed and tagged
- [ ] Phase 8: Documentation updated
- [ ] Phase 9: (Optional) Published to PyPI

---

## Notes

- All `git mv` commands preserve file history
- Archive directories allow easy reference to old code
- Data files are consolidated in `intronIC/data/`
- Default pretrained model enables quick classification without training
- v1.5.1 remains accessible via git tag and branch
