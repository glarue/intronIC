# intronIC Filesystem Reference & Testing Guide

**Date:** 2025-11-07
**Purpose:** Eliminate confusion about paths, file locations, and how to run tests

---

## Directory Structure

```
/home/user/
├── intronIC/                              # GIT REPO ROOT (THIS IS THE REPO)
│   │
│   ├── .git/                              # Git metadata
│   │
│   ├── intronIC/                          # ORIGINAL CODE (Python package)
│   │   ├── intronIC.py                    # Monolithic original (6093 lines)
│   │   ├── __init__.py
│   │   ├── _version.py
│   │   ├── data/                          # Training matrices & reference data
│   │   │   ├── scoring_matrices.fasta.iic
│   │   │   ├── u12_reference.introns.iic.gz  # 387 U12 introns
│   │   │   ├── u2_reference.introns.iic.gz   # 20,690 U2 introns
│   │   │   ├── u12_reference_small.introns.iic.gz  # Smaller test set
│   │   │   └── u2_reference_small.introns.iic.gz
│   │   └── test_data/                     # Test datasets
│   │       ├── Homo_sapiens.Chr19.Ensembl_91.fa.gz     # 16MB
│   │       └── Homo_sapiens.Chr19.Ensembl_91.gff3.gz   # 2.2MB
│   │
│   ├── intronIC_refactored/               # REFACTORED CODE (modular Python package)
│   │   ├── cli/                           # Command-line interface
│   │   │   └── main.py                    # Entry point for refactored version
│   │   ├── classification/                # SVM training & prediction
│   │   │   ├── optimizer.py               # Hyperparameter optimization
│   │   │   ├── trainer.py                 # Model training
│   │   │   └── predictor.py               # Classification
│   │   ├── core/                          # Core data structures
│   │   │   └── intron.py                  # Intron class definition
│   │   ├── extraction/                    # Intron extraction logic
│   │   ├── file_io/                       # File reading/writing
│   │   ├── scoring/                       # PWM scoring
│   │   ├── utils/                         # Utilities
│   │   └── intronIC/                      # DATA DIRECTORY (shared with original)
│   │       └── data/                      # Same reference data as original
│   │           ├── scoring_matrices.fasta.iic
│   │           ├── u12_reference.introns.iic.gz
│   │           ├── u2_reference.introns.iic.gz
│   │           └── ... (same as intronIC/data/)
│   │
│   ├── comparison_test/                   # Test outputs (INSIDE repo)
│   │   ├── original_output/               # Original version results
│   │   ├── refactored_output/             # Refactored version results
│   │   ├── original_run.log               # Original logs
│   │   └── refactored_run.log             # Refactored logs
│   │
│   ├── test_*.py                          # Test scripts (various)
│   ├── LINEAR_SVC_INVESTIGATION.md        # LinearSVC investigation
│   ├── PERFORMANCE_COMPARISON_RESULTS.md  # Performance test results
│   └── ... (various documentation files)
│
└── comparison_test/                       # DUPLICATE test outputs (OUTSIDE repo)
    ├── original_output/
    ├── refactored_output/
    ├── refactored_bestpractices/
    ├── refactored_fixed/
    └── ... (scattered test runs)
```

---

## Key File Locations

### Original Code
| Item | Absolute Path |
|------|---------------|
| **Main script** | `/home/user/intronIC/intronIC/intronIC.py` |
| **Test data** | `/home/user/intronIC/intronIC/test_data/` |
| **Training data** | `/home/user/intronIC/intronIC/data/` |
| **Package root** | `/home/user/intronIC/intronIC/` |

### Refactored Code
| Item | Absolute Path |
|------|---------------|
| **Main entry** | `/home/user/intronIC/intronIC_refactored/cli/main.py` |
| **Test data** | `/home/user/intronIC/intronIC/test_data/` (SHARED with original!) |
| **Training data** | `/home/user/intronIC/intronIC_refactored/intronIC/data/` |
| **Package root** | `/home/user/intronIC/intronIC_refactored/` |
| **Modules** | `/home/user/intronIC/intronIC_refactored/{cli,classification,core,etc}/` |

### Test Outputs (RECOMMENDED LOCATION)
| Item | Absolute Path |
|------|---------------|
| **Original output** | `/home/user/intronIC/comparison_test/original_output/` |
| **Refactored output** | `/home/user/intronIC/comparison_test/refactored_output/` |
| **Test logs** | `/home/user/intronIC/comparison_test/*.log` |

---

## Running the Original Version

### Working Directory: `/home/user/intronIC/intronIC` (package root)

```bash
cd /home/user/intronIC/intronIC

# Run with test data
python intronIC.py \
    -n homo_sapiens_original \
    -g test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
    -a test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz \
    --reference_u12s data/u12_reference.introns.iic.gz \
    --reference_u2s data/u2_reference.introns.iic.gz \
    -o ../comparison_test/original_output \
    -p 4
```

**Key Points:**
- Must run from `/home/user/intronIC/intronIC/` directory
- Uses relative paths: `test_data/`, `data/`, `../comparison_test/`
- Output goes to `/home/user/intronIC/comparison_test/original_output/`
- Command: `python intronIC.py` (direct script invocation)

---

## Running the Refactored Version

### Working Directory: `/home/user/intronIC/intronIC_refactored` (package root)

```bash
cd /home/user/intronIC/intronIC_refactored

# Run with test data
python -m cli.main \
    -n homo_sapiens_refactored \
    -g ../intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
    -a ../intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz \
    --reference_u12s intronIC/data/u12_reference.introns.iic.gz \
    --reference_u2s intronIC/data/u2_reference.introns.iic.gz \
    -o ../comparison_test/refactored_output \
    -p 4
```

**Key Points:**
- Must run from `/home/user/intronIC/intronIC_refactored/` directory
- Test data: `../intronIC/test_data/` (goes UP to intronIC/, then into intronIC/)
- Training data: `intronIC/data/` (local intronIC/ directory)
- Output: `../comparison_test/refactored_output/`
- Command: `python -m cli.main` (module invocation)

---

## Data Sharing Architecture

### Test Data (Genome & Annotation)
**Location:** `/home/user/intronIC/intronIC/test_data/`
**Shared:** YES - Both versions use the SAME test data
**Why:** Test datasets are large (18MB), don't duplicate them

### Training Data (Reference Introns)
**Original:** `/home/user/intronIC/intronIC/data/`
**Refactored:** `/home/user/intronIC/intronIC_refactored/intronIC/data/`
**Shared:** NO - Separate copies (but identical content)
**Why:** Each version has its own data directory for independence

---

## Running Complete Comparison Tests

### Option 1: Using Single Test Script (TODO - create this)

```bash
cd /home/user/intronIC

# Run both versions and compare
./run_comparison_test.sh
```

### Option 2: Manual Sequential Runs

```bash
cd /home/user/intronIC

# Clean output directories
rm -rf comparison_test/original_output/*
rm -rf comparison_test/refactored_output/*

# Run original
cd intronIC
time python intronIC.py \
    -n homo_sapiens_original \
    -g test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
    -a test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz \
    --reference_u12s data/u12_reference.introns.iic.gz \
    --reference_u2s data/u2_reference.introns.iic.gz \
    -o ../comparison_test/original_output \
    -p 4 \
    2>&1 | tee ../comparison_test/original_run.log

# Run refactored
cd ../intronIC_refactored
time python -m cli.main \
    -n homo_sapiens_refactored \
    -g ../intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
    -a ../intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz \
    --reference_u12s intronIC/data/u12_reference.introns.iic.gz \
    --reference_u2s intronIC/data/u2_reference.introns.iic.gz \
    -o ../comparison_test/refactored_output \
    -p 4 \
    2>&1 | tee ../comparison_test/refactored_run.log

# Compare results
cd ..
echo "=== Original Results ===" && tail -50 comparison_test/original_run.log
echo -e "\n=== Refactored Results ===" && tail -50 comparison_test/refactored_run.log
```

---

## Common Pitfalls & Solutions

### ❌ Problem: "No module named cli"
**Cause:** Running from wrong directory
**Solution:** Must be in `/home/user/intronIC/intronIC_refactored/`

### ❌ Problem: "File not found: intronIC/data/..."
**Cause:** Wrong path specification from current directory
**Solution:** Check working directory, adjust paths accordingly

### ❌ Problem: Output goes to unexpected location
**Cause:** Using absolute path that points somewhere else
**Solution:** Use relative paths from documented working directories

### ❌ Problem: Test data not found
**Cause:** Test data is in `intronIC/intronIC/test_data/`, not `intronIC/test_data/`
**Solution:** Use correct path: `../intronIC/test_data/` from refactored directory

### ❌ Problem: Two comparison_test directories
**Status:** KNOWN ISSUE - There are TWO:
1. `/home/user/intronIC/comparison_test/` (PREFERRED - inside repo)
2. `/home/user/comparison_test/` (outside repo, should be cleaned up)

**Recommendation:**
- Use ONLY `/home/user/intronIC/comparison_test/` (inside repo)
- Clean up `/home/user/comparison_test/` when convenient
- All test commands above use the in-repo location

---

## Quick Reference Commands

### Check Current Working Directory
```bash
pwd
# Expected: /home/user/intronIC/intronIC_refactored or /home/user/intronIC/intronIC
```

### Verify File Existence
```bash
# From refactored directory
ls -lh ../intronIC/test_data/
ls -lh intronIC/data/

# From original directory
ls -lh test_data/
ls -lh data/
```

### Check Git Branch
```bash
cd /home/user/intronIC
git branch
git status
```

### View Recent Test Logs
```bash
cd /home/user/intronIC
tail -100 comparison_test/original_run.log
tail -100 comparison_test/refactored_run.log
```

---

## Test Output Files

When tests complete successfully, you should see these files:

### Original Output
```
/home/user/intronIC/comparison_test/original_output/
├── homo_sapiens_original.bed.iic
├── homo_sapiens_original.meta.iic
├── homo_sapiens_original.introns.iic
├── homo_sapiens_original.score_info.iic
├── homo_sapiens_original.log
└── homo_sapiens_original.*.png (plots)
```

### Refactored Output
```
/home/user/intronIC/comparison_test/refactored_output/
├── homo_sapiens_refactored.bed.iic
├── homo_sapiens_refactored.meta.iic
├── homo_sapiens_refactored.introns.iic
├── homo_sapiens_refactored.score_info.iic
├── homo_sapiens_refactored.log
└── homo_sapiens_refactored.*.png (plots)
```

---

## Performance Testing Results Location

All performance investigation documents are in the repo root:

```
/home/user/intronIC/
├── PERFORMANCE_COMPARISON_RESULTS.md     # Latest comparison (4.5x speedup)
├── PERFORMANCE_TEST_RESULTS.md           # Component tests
├── LINEAR_SVC_INVESTIGATION.md           # Current investigation
├── INVESTIGATION_COMPLETE.md             # Previous investigation summary
└── PERFORMANCE_TESTING_STRATEGY.md       # Testing methodology
```

---

## Environment Variables (None Required)

intronIC does not require any special environment variables. All paths are specified via command-line arguments.

---

## Python Environment

```bash
python --version
# Expected: Python 3.11+

pip list | grep -E "scikit-learn|numpy|pandas"
# Expected:
#   scikit-learn  1.3.x+
#   numpy         1.24.x+ (but < 2.0)
#   pandas        2.x+
```

---

## Git Workflow

**Current Branch:** `claude/investigate-ref-performance-011CUt9YD4yNb4AzXF2TKAU9`

```bash
cd /home/user/intronIC

# View current branch
git branch

# View recent commits
git log --oneline -10

# View changed files
git status

# Push changes
git push -u origin claude/investigate-ref-performance-011CUt9YD4yNb4AzXF2TKAU9
```

---

## Cleanup Recommendations

### Immediate
- Remove `/home/user/comparison_test/` (duplicate, outside repo)
- Move all test outputs to `/home/user/intronIC/comparison_test/`

### Future
- Create `run_comparison_test.sh` script to automate testing
- Add `.gitignore` entry for `comparison_test/` outputs
- Create separate `test_outputs/` directory for experimental runs

---

## Summary Cheat Sheet

| Task | Working Dir | Command |
|------|-------------|---------|
| **Run Original** | `/home/user/intronIC/intronIC` | `python intronIC.py -n NAME -g test_data/... -a test_data/... --reference_u12s data/... --reference_u2s data/... -o ../comparison_test/original_output` |
| **Run Refactored** | `/home/user/intronIC/intronIC_refactored` | `python -m cli.main -n NAME -g ../intronIC/test_data/... -a ../intronIC/test_data/... --reference_u12s intronIC/data/... --reference_u2s intronIC/data/... -o ../comparison_test/refactored_output` |
| **View Logs** | `/home/user/intronIC` | `tail -100 comparison_test/{original,refactored}_run.log` |
| **Check Status** | `/home/user/intronIC` | `git status && git log --oneline -5` |

---

**Last Updated:** 2025-11-07
**Maintainer:** Development team working on performance investigation
