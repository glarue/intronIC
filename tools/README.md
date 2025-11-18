# intronIC Analysis Tools

Post-hoc analysis and debugging tools for model predictions and behavior.

---

## Tools

### analyze_predictions.py

**Purpose:** Analyze prediction patterns from score_info.iic files

**What it does:**
- Manually computes all BothEndsStrong features from z-scores
- Identifies problematic patterns (one-end-strong, negative scores, imbalance)
- Helps diagnose false positives and understand model behavior

**Usage:**
```bash
# Basic usage - analyze top 10 predicted U12s
python tools/analyze_predictions.py run_tests/species.score_info.iic

# Show only top 5
python tools/analyze_predictions.py run_tests/species.score_info.iic -n 5

# Summary table only (faster)
python tools/analyze_predictions.py run_tests/species.score_info.iic --summary

# Filter by score range
python tools/analyze_predictions.py run_tests/species.score_info.iic --min-score 5.0

# Show high and low scoring predictions
python tools/analyze_predictions.py run_tests/species.score_info.iic --show-extremes
```

**When to use:**
- Investigating false positives
- Comparing model versions (before/after fix)
- Understanding species-specific issues
- Diagnosing why certain predictions scored high

**See:** `docs/POST_HOC_ANALYSIS.md` for detailed methodology

---

### inspect_model.py

**Purpose:** Inspect trained model coefficients and configuration

**What it does:**
- Loads pickled model files
- Shows feature names and learned coefficients
- Displays transformer configuration (include_max, include_pairwise_mins)
- Shows all models in ensemble

**Usage:**
```bash
# Default: inspect run_tests/hsapiens/homo_sapiens.model.pkl
python tools/inspect_model.py

# Inspect specific model
python tools/inspect_model.py path/to/model.pkl
```

**Output:**
```
Features (7): ['s5', 'sBP', 's3', 'min_all', ...]

Coefficients (7 features):
  s5                  : +1.249035
  sBP                 : +0.543210
  s3                  : +0.321098
  min_all             : +0.654321  ← Should be non-zero!
  neg_absdiff_5_bp    : +0.234567
  neg_absdiff_5_3     : +0.345678
  neg_absdiff_bp_3    : +0.456789
```

**When to use:**
- After training to verify feature selection
- Diagnosing L1 penalty issues (zeroed features)
- Checking if critical features (min_all, neg_absdiff) are active
- Understanding which features the model relies on

**Red flags:**
- min_all coefficient = 0 (L1 zeroed it)
- All neg_absdiff coefficients = 0 (no imbalance penalty)
- max features present with positive coefficients (rewarding one-end-strong)
- sBP or s3 coefficients = 0 (only using 5'SS)

---

## Typical Workflow: Diagnosing False Positives

### 1. Check predictions
```bash
python tools/analyze_predictions.py species.score_info.iic -n 20 --summary
```

Look for patterns in warnings column:
- Many "NEGATIVE scores" → model not enforcing minimum thresholds
- Many "ONE-END-STRONG" → model rewarding imbalance
- Many "min_all is NEGATIVE" → min_all not enforced

### 2. Inspect model
```bash
python tools/inspect_model.py path/to/model.pkl
```

Check if critical features are active:
- Is min_all coefficient non-zero?
- Are neg_absdiff coefficients non-zero?
- Are max features absent (or zero)?

### 3. Diagnose root cause

**If min_all = 0:**
- L1 penalty zeroed it due to redundancy with min_5_bp, min_5_3
- Fix: Remove redundant features, use L2 penalty

**If neg_absdiff = 0:**
- L1 penalty didn't see them as useful (no imbalance in training data)
- Fix: Force L2 to ensure they contribute

**If max features present:**
- Rewarding "at least one strong" instead of "all strong"
- Fix: Remove max features from transformer

### 4. Verify fix
```bash
# Retrain with fix
intronIC train -n species -p 12

# Re-analyze predictions
python tools/analyze_predictions.py species_v2.score_info.iic -n 20 --summary

# Compare
diff <(python tools/analyze_predictions.py old.score_info.iic --summary) \
     <(python tools/analyze_predictions.py new.score_info.iic --summary)
```

Expected improvements:
- Fewer total predictions
- No predictions with negative min_all
- No one-end-strong patterns
- Smaller imbalance values

---

## Example: L1 Penalty Investigation

**Scenario:** Model predicting false positives in C. elegans

**Step 1: Analyze predictions**
```bash
$ python tools/analyze_predictions.py caenorhabditis_elegans.score_info.iic

#1: RelScore=10.00, min_all=-1.90, Pattern="NEGATIVE scores: s3"
#5: RelScore=7.96, min_all=-0.15, Pattern="ONE-END-STRONG"
```

**Finding:** 6/8 predictions have negative min_all (should be rejected!)

**Step 2: Inspect model**
```bash
$ python tools/inspect_model.py homo_sapiens.model.pkl

min_all             : +0.000000  ← ZEROED!
max_5_3             : +0.396311  ← REWARDING IMBALANCE!
```

**Finding:** L1 penalty zeroed min_all and kept max_5_3

**Step 3: Root cause**
- L1 saw min_all as redundant with min_5_bp, min_5_3
- L1 kept pairwise mins, zeroed the most important one
- L1 kept max_5_3 which rewards one-end-strong patterns

**Step 4: Fix**
- Remove redundant features (keep only min_all)
- Force L2 penalty
- Remove max features

See: `FP_ROOT_CAUSE_ANALYSIS.md` for complete analysis

---

## Files Created

- `tools/analyze_predictions.py` - Main analysis tool
- `tools/inspect_model.py` - Model coefficient inspector
- `docs/POST_HOC_ANALYSIS.md` - Detailed methodology guide
- `tools/README.md` - This file

---

## Dependencies

Both tools use only standard libraries and intronIC components:
- numpy
- pathlib
- argparse
- joblib (for model loading)

No additional installation required if intronIC is already set up.

---

**Created:** 2025-11-16
**Purpose:** Post-hoc diagnostics for model predictions and false positive investigation
