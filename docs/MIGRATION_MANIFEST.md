# intronIC v2.0.0 Migration Manifest

## Overview
This document defines which files to keep and their final destinations when migrating from the refactored subdirectory to the root-level structure.

---

## 1. Core Application Code

**Source:** `intronIC_refactored/`
**Destination:** `/` (root)

### Module Directories (keep all)
- `cli/` → `/cli/`
- `core/` → `/core/`
- `extraction/` → `/extraction/`
- `file_io/` → `/file_io/`
- `scoring/` → `/scoring/`
- `classification/` → `/classification/`
- `utils/` → `/utils/`
- `visualization/` → `/visualization/`
- `tests/` → `/tests/`
- `scripts/` → `/scripts/`

---

## 2. Data Files

**Source:** `intronIC/data/` and `intronIC_refactored/intronIC/data/`
**Destination:** `/intronIC/data/`

### Required Files (production data only):

#### PWM Matrices
- `scoring_matrices.fasta.iic` → `/intronIC/data/scoring_matrices.fasta.iic`
- `scoring_matrices.2+_pubs.fasta.iic` → `/intronIC/data/scoring_matrices.2+_pubs.fasta.iic`
- `scoring_matrices.2+.pub_required.fasta.iic` → `/intronIC/data/scoring_matrices.2+.pub_required.fasta.iic`
- `u2.conserved_empirical_bp_pwm.iic` → `/intronIC/data/u2.conserved_empirical_bp_pwm.iic`
- `empirical.u12s.2+.pub_required.matrices` → `/intronIC/data/empirical.u12s.2+.pub_required.matrices`

#### Reference Sequences (full datasets only)
- `u12_reference.introns.iic.gz` → `/intronIC/data/u12_reference.introns.iic.gz`
- `u2_reference.introns.iic.gz` → `/intronIC/data/u2_reference.introns.iic.gz`

#### Pretrained Models
- `homo_sapiens.model.pkl` → `/intronIC/data/default_pretrained.model.pkl`

### Test Data

**Source:** `intronIC/test_data/`
**Destination:** `/intronIC/data/test_data/`

- `Homo_sapiens.Chr19.Ensembl_91.fa.gz` → `/intronIC/data/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz`
- `Homo_sapiens.Chr19.Ensembl_91.gff3.gz` → `/intronIC/data/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz`

### Excluded Data Files (test/debug only, not needed):
- `u12_reference_small.introns.iic.gz` ❌ (testing only)
- `u2_reference_small.introns.iic.gz` ❌ (testing only)
- `u12_reference_micro.introns.iic.gz` ❌ (testing only)
- `u2_reference_micro.introns.iic.gz` ❌ (testing only)
- `u12_reference_tiny.introns.iic.gz` ❌ (testing only)
- `u2_reference_tiny.introns.iic.gz` ❌ (testing only)

---

## 3. Configuration & Packaging Files

**Source:** `intronIC_refactored/`
**Destination:** `/` (root)

### Package Configuration
- `pyproject.toml` → `/pyproject.toml` ✅
- `pixi.toml` → `/pixi.toml` ✅
- `pixi.lock` → `/pixi.lock` ✅
- `MANIFEST.in` → `/MANIFEST.in` ✅
- `__init__.py` → `/__init__.py` ✅
- `__main__.py` → `/__main__.py` ✅

### Installation Helpers
- `install.sh` → `/install.sh` ✅
- `install.bat` → `/install.bat` ✅

### Configuration Files from Root (keep existing)
- `.gitignore` → `/.gitignore` ✅ (merge if needed)
- `.gitattributes` → `/.gitattributes` ✅ (keep existing)

---

## 4. Documentation Files

**Source:** `intronIC_refactored/`
**Destination:** `/` (root)

### User Documentation
- `README.md` → `/README.md` ✅
- `INSTALL.md` → `/INSTALL.md` ✅
- `QUICKSTART.md` → `/QUICKSTART.md` ✅
- `LICENSE` → `/LICENSE` ✅ (keep existing from root)

### Developer Documentation
- `PACKAGING_STRATEGY.md` → `/docs/PACKAGING_STRATEGY.md`
- `VERSIONING_STRATEGY.md` → `/docs/VERSIONING_STRATEGY.md`
- `VERSIONING_COMPARISON.md` → `/docs/VERSIONING_COMPARISON.md`
- `RELEASE_PROCESS.md` → `/docs/RELEASE_PROCESS.md`
- `REFACTOR_SUMMARY.md` → `/docs/REFACTOR_SUMMARY.md`
- `TEST_IMPROVEMENT_PLAN.md` → `/docs/TEST_IMPROVEMENT_PLAN.md`
- `IMPLEMENTATION_PLAN.md` → `/docs/IMPLEMENTATION_PLAN.md`
- `PERFORMANCE_IMPROVEMENTS.md` → `/docs/PERFORMANCE_IMPROVEMENTS.md`

### Testing Configuration
- `pytest.ini` → `/pytest.ini` ✅

---

## 5. Files to Archive

**Destination:** `/archive/` (new directory)

### v1.5.1 Original Implementation
- `intronIC/intronIC.py` → `/archive/v1.5.1/intronIC/intronIC.py`
- `intronIC/pwms_from_seqs.py` → `/archive/v1.5.1/intronIC/pwms_from_seqs.py`
- `intronIC/matrices_from_seqs.py` → `/archive/v1.5.1/intronIC/matrices_from_seqs.py`
- `intronIC/_version.py` → `/archive/v1.5.1/intronIC/_version.py`
- `intronIC/__init__.py` → `/archive/v1.5.1/intronIC/__init__.py`

### Old Packaging Files (replaced by pyproject.toml)
- `setup.py` → `/archive/v1.5.1/setup.py`
- `setup.cfg` → `/archive/v1.5.1/setup.cfg`
- `versioneer.py` → `/archive/v1.5.1/versioneer.py`

### Investigation/Debug Files
All root-level `.py` and `.md` files related to debugging:
- `analyze_*.py` → `/archive/investigation/`
- `compare_*.py` → `/archive/investigation/`
- `debug_*.py` → `/archive/investigation/`
- `diagnose_*.py` → `/archive/investigation/`
- `investigate_*.py` → `/archive/investigation/`
- `profile_*.py` → `/archive/investigation/`
- `test_*.py` (root level only, not tests/) → `/archive/investigation/`
- `*_INVESTIGATION.md` → `/archive/investigation/`
- `*_ANALYSIS.md` → `/archive/investigation/`
- `*_AUDIT*.md` → `/archive/investigation/`
- `*_DISCREPANCY*.md` → `/archive/investigation/`
- `*_SUMMARY.md` (session summaries) → `/archive/investigation/`
- `*_RESULTS.md` → `/archive/investigation/`
- `*.log` files → `/archive/investigation/`
- `*.pstats` files → `/archive/investigation/`

### Refactored Session Files
- `intronIC_refactored/SESSION_*.md` → `/archive/refactoring/`
- `intronIC_refactored/SETUP*.md` → `/archive/refactoring/`
- `intronIC_refactored/VALIDATION*.md` → `/archive/refactoring/`
- `intronIC_refactored/PROGRESS*.md` → `/archive/refactoring/`
- `intronIC_refactored/SCORING_ANALYSIS.md` → `/archive/refactoring/`
- `intronIC_refactored/U12_CORRECTION_*.md` → `/archive/refactoring/`
- `intronIC_refactored/PORTING_CHECKLIST.md` → `/archive/refactoring/`
- `intronIC_refactored/REVIEW_FINDINGS.md` → `/archive/refactoring/`
- `intronIC_refactored/TODO_FROM_REVIEW.md` → `/archive/refactoring/`
- `intronIC_refactored/DYNAMIC_TAGGING_IMPLEMENTATION.md` → `/archive/refactoring/`
- `intronIC_refactored/OUTPUT_FORMAT_REFACTOR_PLAN.md` → `/archive/refactoring/`

### Test Output Files
- `intronIC_refactored/*.bed.iic` → `/archive/test_outputs/`
- `intronIC_refactored/*.introns.iic` → `/archive/test_outputs/`
- `intronIC_refactored/*.meta.iic` → `/archive/test_outputs/`
- `intronIC_refactored/*.score_info.iic` → `/archive/test_outputs/`
- `intronIC_refactored/*.plot.*.png` → `/archive/test_outputs/`
- `intronIC_refactored/*.log` → `/archive/test_outputs/`
- `intronIC_refactored/*.metrics.json` → `/archive/test_outputs/`
- `intronIC_refactored/*.model.pkl` → `/archive/test_outputs/`
- `intronIC_refactored/*.tar.gz` → `/archive/test_outputs/`
- `intronIC_refactored/drosophila/` → `/archive/test_outputs/drosophila/`

### Comparison Test Data
- `comparison_test/` → `/archive/comparison_test/`

### Backup Directories
- `intronIC_refactored/classification_backup/` → `/archive/classification_backup/`

---

## 6. Files to Delete (not needed)

### Build Artifacts
- `intronIC_refactored/__pycache__/` ❌
- `intronIC_refactored/.pytest_cache/` ❌
- `intronIC_refactored/.pixi/` ❌ (will be regenerated)
- Root `__pycache__/` ❌
- `*.pyc` files ❌

### Temporary/Generated Files
- `intronIC_refactored/.gitignore` ❌ (redundant with root)
- `intronIC_refactored/test_run.log` ❌ (empty)
- Root level `*.log` files ❌ (already archived)
- Root level `*.pstats` files ❌ (already archived)

### Old Reference Files (in root, superseded)
- `u2_reference_small.introns.iic.gz` (root level) ❌
- `u12_reference_small.introns.iic.gz` (root level) ❌

---

## 7. Files to Keep in Root (untouched)

### Project Metadata
- `CLAUDE.md` ✅ (keep existing)
- `README.md` ✅ (replace with refactored version)
- `LICENSE` ✅ (keep existing)
- `_config.yml` ✅ (GitHub Pages config)
- `.gitattributes` ✅
- `.gitignore` ✅
- `.git/` ✅
- `.claude/` ✅

### Reference Directories
- `images/` ✅ (project images)

---

## Summary Counts

### Keep & Move
- **Code modules:** 9 directories
- **Data files:** 7 files (PWMs + references)
- **Test data:** 2 files (Chr19)
- **Config files:** 8 files
- **Documentation:** 13 files

### Archive
- **v1.5.1 code:** ~10 files
- **Investigation files:** ~50+ files
- **Refactoring docs:** ~15 files
- **Test outputs:** ~30+ files

### Delete
- **Build artifacts:** 3 directories
- **Temporary files:** ~10+ files
