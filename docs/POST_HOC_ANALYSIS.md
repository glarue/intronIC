# Post-Hoc Prediction Analysis Methodology

**Purpose:** Diagnose false positives and understand model behavior by manually computing features from z-scores

**Tool:** `tools/analyze_predictions.py`

---

## When to Use This Analysis

### 1. Investigating False Positives

**Scenario:** Model predicts U12s in species known to have none/few

**Example:**
```bash
# C. elegans has ~0 known U12s, but model predicts 8
python tools/analyze_predictions.py \
    run_tests/caenorhabditis_elegans.score_info.iic
```

**What to look for:**
- ⚠ Negative base scores (s5, sBP, or s3 < 0)
- ⚠ One-end-strong patterns (one score >>1, others weak/negative)
- ⚠ High imbalance (|difference| > 2-3 standard deviations)
- ⚠ Negative min_all (should strongly reject)

### 2. Comparing Model Versions

**Scenario:** Testing if new model fixes FP issues

**Example:**
```bash
# Old model
python tools/analyze_predictions.py \
    run_tests/species.old_model.score_info.iic \
    --summary > old_summary.txt

# New model
python tools/analyze_predictions.py \
    run_tests/species.new_model.score_info.iic \
    --summary > new_summary.txt

# Compare
diff old_summary.txt new_summary.txt
```

**What to compare:**
- Number of predictions with warnings
- Pattern of min_all values
- Magnitude of imbalance penalties

### 3. Understanding Training Issues

**Scenario:** Model coefficients seem wrong (e.g., features zeroed by L1)

**Example:**
```bash
# Check reference data first
python tools/analyze_predictions.py \
    data/u12_reference.score_info.iic \
    --min-score -10 --max-score 10

# Look for:
# - Are reference U12s balanced? (small neg_absdiff)
# - Are they all strong? (min_all > 1.0)
# - Any edge cases in training data?
```

### 4. Species-Specific Issues

**Scenario:** Model works on vertebrates but fails on invertebrates

**Example:**
```bash
# Compare top predictions across species
python tools/analyze_predictions.py human.score_info.iic -n 5 --summary
python tools/analyze_predictions.py celegans.score_info.iic -n 5 --summary
python tools/analyze_predictions.py drosophila.score_info.iic -n 5 --summary
```

**What to look for:**
- Different z-score distributions
- Different patterns of imbalance
- Consistent vs. inconsistent warnings

---

## How to Interpret Output

### Summary Table

```
#    RelScore   SVM%      s5     sBP      s3  min_all Pattern
1        10.00  100.0    0.45    1.37   -1.90    -1.90 NEGATIVE scores: s3
2        10.00  100.0    0.28    1.12    0.41     0.28 OK
3         7.96   98.0   -0.15    0.30   11.45    -0.15 ONE-END-STRONG
```

**Columns:**
- `RelScore`: Distance from decision boundary (>0 = predicted U12)
- `SVM%`: Calibrated probability (90% = threshold by default)
- `s5, sBP, s3`: Base z-scores (LLR, zero-anchored)
- `min_all`: Strictest "all strong" constraint
- `Pattern`: First warning detected (or "OK")

**Interpretation:**
- **#1 (FP):** s3 is NEGATIVE (-1.90) → should reject, but min_all not enforced
- **#2 (OK):** All positive, balanced (small differences)
- **#3 (FP):** One-end-strong (s3=11.45 very high, s5 negative) → clear FP

### Detailed Analysis

```
  Base z-scores:
    s5:   0.4450
    sBP:  1.3708
    s3:  -1.9010

  Min features (require all/both strong):
    min_all:   -1.9010  ← ALL THREE must be strong
    ...

  Imbalance penalties (should be small for real U12s):
    neg_absdiff_5_bp: -0.9258  ← 5'/BP balance
    neg_absdiff_5_3:  -2.3460  ← 5'/3' balance
    neg_absdiff_bp_3: -3.2718  ← BP/3' balance

  ⚠ Pattern warnings:
    • NEGATIVE scores: s3
    • ONE-END-STRONG pattern detected
    • High BP/3' imbalance (|diff|=3.27)
    • min_all is NEGATIVE (should strongly penalize)
```

**Red flags:**
1. **min_all negative** → At least one signal is weak/negative (should reject)
2. **Large imbalance** (>2-3) → Inconsistent signals (likely FP)
3. **One-end-strong** → Classic FP pattern (one motif matches by chance)
4. **Negative base scores** → Below background level (should reject)

**Green flags (real U12s):**
1. **All z-scores > 1.0** → All signals strong
2. **min_all > 1.0** → Consistent strength
3. **Small imbalances** (<0.5) → Balanced signals
4. **No warnings** → Looks like real U12

### Feature Computation Details

The tool computes ALL possible BothEndsStrong features:

**Base (3):** s5, sBP, s3

**Min features (4):**
- `min_5_bp = min(s5, sBP)` - Both 5' and BP strong
- `min_5_3 = min(s5, s3)` - Both 5' and 3' strong
- `min_bp_3 = min(sBP, s3)` - Both BP and 3' strong
- `min_all = min(s5, sBP, s3)` - ALL THREE strong (strictest)

**Max features (4):**
- `max_5_bp = max(s5, sBP)` - At least one of 5'/BP strong
- `max_5_3 = max(s5, s3)` - At least one of 5'/3' strong
- `max_bp_3 = max(sBP, s3)` - At least one of BP/3' strong
- `max_all = max(s5, sBP, s3)` - At least one strong (weakest)

**Imbalance penalties (3):**
- `neg_absdiff_5_bp = -|s5 - sBP|` - Penalty for 5'/BP imbalance
- `neg_absdiff_5_3 = -|s5 - s3|` - Penalty for 5'/3' imbalance
- `neg_absdiff_bp_3 = -|sBP - s3|` - Penalty for BP/3' imbalance

**Note:** Model may only use subset of these (e.g., 7 features: base + min_all + penalties)

---

## Common Diagnosis Workflows

### Workflow 1: Why are there false positives?

```bash
# Step 1: Analyze top predictions
python tools/analyze_predictions.py species.score_info.iic -n 20

# Step 2: Count warning types
python tools/analyze_predictions.py species.score_info.iic --summary | \
    grep -oE "(NEGATIVE|ONE-END|imbalance|min_all)" | sort | uniq -c

# Step 3: Check model coefficients
python tools/inspect_model.py

# Step 4: Compare to expected behavior
# - Are critical features (min_all, neg_absdiff) in the model?
# - Do they have non-zero coefficients?
# - Are max features present (bad - reward imbalance)?
```

**Common findings:**
- L1 penalty zeroed min_all → no "all strong" enforcement
- L1 penalty zeroed neg_absdiff → no imbalance penalty
- max features present with positive coefficients → rewarding one-end-strong
- Missing features → model lacks capacity to reject FPs

### Workflow 2: Did my fix work?

```bash
# Before fix
python tools/analyze_predictions.py old_model.score_info.iic --summary > before.txt

# After fix
python tools/analyze_predictions.py new_model.score_info.iic --summary > after.txt

# Compare
echo "=== BEFORE ==="
head -20 before.txt
echo "=== AFTER ==="
head -20 after.txt

# Check specific metrics
echo "FPs with negative min_all (BEFORE):"
grep "min_all" before.txt | awk '$7 < 0 {count++} END {print count}'

echo "FPs with negative min_all (AFTER):"
grep "min_all" after.txt | awk '$7 < 0 {count++} END {print count}'
```

**Expected after fix:**
- Fewer total predictions
- No predictions with negative min_all
- No predictions with extreme imbalance (>5.0)
- No one-end-strong patterns

### Workflow 3: Understanding species differences

```bash
# Check if patterns are consistent
for species in human mouse celegans drosophila; do
    echo "=== $species ==="
    python tools/analyze_predictions.py ${species}.score_info.iic \
        --summary -n 3 | tail -5
done
```

**What to look for:**
- Are FP patterns similar across species? (suggests model issue)
- Are FP patterns different? (suggests species-specific issue)
- Do vertebrates have different patterns than invertebrates?

---

## Case Study: L1 Penalty Zeroing Critical Features

**Problem:** C. elegans model predicted 8 U12s, all false positives

**Step 1: Analyze predictions**
```bash
python tools/analyze_predictions.py \
    caenorhabditis_elegans.pretrained.score_info.iic -n 10
```

**Findings:**
```
#1: min_all = -1.90 (NEGATIVE!), max_5_3 = 0.45
#4: min_all = -2.76 (NEGATIVE!), max_5_3 = 0.36
#5: min_all = -0.15 (NEGATIVE!), max_5_3 = 11.45 (HUGE!)
```

6 out of 8 FPs had **negative min_all** → should be strongly rejected

**Step 2: Inspect model**
```bash
python tools/inspect_model.py
```

**Findings:**
```
min_all:         +0.000  ← ZEROED! (L1 penalty)
neg_absdiff_5_3: +0.000  ← ZEROED!
max_5_3:         +0.396  ← ACTIVE! (rewards imbalance)
```

**Diagnosis:** L1 zeroed the features designed to reject FPs!

**Step 3: Understand why**
- `min_all = min(min_5_bp, min_5_3)` → redundant
- L1 saw redundancy, kept pairwise mins, zeroed min_all
- But min_all is the STRICTEST constraint (most important!)
- L1 made arbitrary choice based on training data, not domain knowledge

**Solution:**
1. Remove redundant features (keep only min_all)
2. Force L2 penalty (ensures all features contribute)
3. Remove max features (don't reward imbalance)

**Step 4: Verify fix**
```bash
# Retrain model
intronIC train -n homo_sapiens -p 12

# Re-analyze
python tools/analyze_predictions.py \
    caenorhabditis_elegans.new_model.score_info.iic -n 10
```

**Expected:** 0-2 predictions (min_all and neg_absdiff should reject FPs)

---

## Advanced: Manual Feature Contribution Calculation

If you want to understand EXACTLY why a prediction scored high, calculate the linear combination:

```python
# From model inspection
coefficients = {
    's5': 1.249,
    'sBP': 0.000,
    's3': 0.000,
    'min_5_bp': 0.528,
    'min_5_3': 0.172,
    'min_all': 0.000,
    'neg_absdiff_5_bp': 0.000,
    'neg_absdiff_5_3': 0.000,
    'neg_absdiff_bp_3': 0.000,
    'max_5_3': 0.396
}
intercept = -0.199

# From analyze_predictions.py output
features = {
    's5': -0.154,
    'sBP': 0.301,
    's3': 11.451,
    'min_5_bp': -0.154,
    'min_5_3': -0.154,
    'min_all': -0.154,
    'neg_absdiff_5_bp': -0.454,
    'neg_absdiff_5_3': -11.605,
    'neg_absdiff_bp_3': -11.150,
    'max_5_3': 11.451
}

# Calculate score
score = intercept
for feat_name, feat_val in features.items():
    coef = coefficients.get(feat_name, 0)
    contribution = coef * feat_val
    print(f"{feat_name:20s}: {coef:+7.3f} × {feat_val:+7.3f} = {contribution:+7.3f}")
    score += contribution

print(f"\nTotal score: {score:.3f}")
```

This shows EXACTLY which features contributed positive/negative scores.

---

## Tips

1. **Always check model first:** Use `inspect_model.py` to see what features are active
2. **Start with summary:** Use `--summary` to get overview before detailed analysis
3. **Compare to reference:** Analyze reference U12s to see expected patterns
4. **Look for consistent patterns:** If all FPs have same warning, that's the issue
5. **Combine with plots:** Scatter plots show separation, this tool shows WHY

---

## Files

- **Tool:** `tools/analyze_predictions.py`
- **Model inspector:** `tools/inspect_model.py`
- **Example analysis:** `FP_ROOT_CAUSE_ANALYSIS.md`
- **This guide:** `docs/POST_HOC_ANALYSIS.md`

---

**Created:** 2025-11-16
**Use case:** Post-hoc diagnosis of false positives and model behavior
