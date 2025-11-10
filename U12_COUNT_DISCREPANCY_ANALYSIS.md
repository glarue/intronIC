# U12 Count Discrepancy Analysis

## Problem Statement

When classifying the same experimental data, different total numbers of U12 introns are found depending on the execution mode:
1. **Pretrained model** (`--pretrained_model`)
2. **Fixed C parameter** (`-C` with value from previous optimization)
3. **Full optimization** (normal training from scratch)

## Root Cause: Different Normalization Strategies

### The Two Code Paths

#### Path 1: Pretrained Model (`cli/main.py:1061-1067`)
```python
if config.training.pretrained_model_path:
    classified_introns, metrics = classify_with_pretrained_model(
        scored_introns, config.training.pretrained_model_path, config, reporter, logger
    )
```

**Normalization in `classify_with_pretrained_model()` (lines 723-735):**
```python
# Fit normalizer on experimental data (cross-species domain adaptation)
normalizer = ScoreNormalizer()
normalizer.fit(introns, dataset_type='unlabeled')  # ← EXPERIMENTAL DATA
normalized_introns = list(normalizer.transform(introns, dataset_type='experimental'))
```

#### Path 2: Normal Training (`cli/main.py:1068-1083`)
```python
else:
    # Normal flow: normalize + train + classify
    normalized_introns, u12_reference, u2_reference, normalizer = normalize_scores(
        scored_introns, config, reporter, logger
    )
```

**Normalization in `normalize_scores()` (lines 653-654):**
```python
normalizer = ScoreNormalizer()
normalizer.fit(reference_introns, dataset_type='reference')  # ← REFERENCE DATA
```

### Critical Difference: What Data Normalizer is Fitted On

| Mode | Normalizer Fitted On | Z-Score Formula |
|------|---------------------|-----------------|
| **Pretrained** | Experimental introns | `z = raw / median(\|experimental_raw\|)` |
| **Fixed C / Full** | Reference introns (U12 + U2) | `z = raw / median(\|reference_raw\|)` |

### Why This Causes Different U12 Counts

The `ZeroAnchoredRobustScaler` (`scoring/normalizer.py:81-123`) computes scales as:
```python
# For each feature (five, bp, three)
abs_values = np.abs(feature_values)
winsorized = np.clip(abs_values, 0, np.quantile(abs_values, 0.995))
scale = np.median(winsorized)  # ← This is VERY different for experimental vs reference!
```

**Experimental data distribution:**
- ~99.5% U2-type introns (weak U12 signals, strong U2 signals)
- ~0.5% U12-type introns (strong U12 signals, weak U2 signals)
- Median dominated by U2 characteristics

**Reference data distribution:**
- Curated mix of U12 and U2 introns
- Balanced or imbalanced depending on reference sets
- Median reflects both U12 and U2 signal strengths

**Result:** Different scales → Different z-scores → Different SVM predictions → Different U12 counts!

### Example Scenario

Imagine branch point scores:

**Reference introns:**
- U12 references: bp_raw_score ~ [8, 10, 12] (strong U12 signal)
- U2 references: bp_raw_score ~ [-5, -3, -2] (strong U2 signal)
- `median(|reference_bp|) ≈ 5.0`

**Experimental introns (99.5% U2):**
- Most: bp_raw_score ~ [-5, -3, -2] (U2-like)
- Few: bp_raw_score ~ [8, 10] (U12-like)
- `median(|experimental_bp|) ≈ 3.0` (dominated by U2 majority)

**For a marginal U12 with bp_raw_score = 6:**
- Pretrained normalization: `z = 6 / 3.0 = 2.0`
- Normal normalization: `z = 6 / 5.0 = 1.2`

The SVM sees completely different input values!

## Design Intent vs. Actual Use Case

### The Pretrained Model Path Was Designed For Cross-Species Adaptation

From `cli/main.py:676-683`:
```python
"""Classify introns using a pretrained model with cross-species domain adaptation.

This function enables applying a trained SVM to new species without curated references.
It implements unsupervised domain adaptation by fitting the normalizer on unlabeled
experimental data, which is statistically valid because:
- 99.5% of introns are U2-type → robust estimators learn species-specific U2 baseline
- No label leakage (normalization uses only marginal feature distribution)
- Corrects covariate shift from species-specific sequence characteristics
"""
```

**Expected use case:**
1. Train on Species A (e.g., human) with curated references
2. Save model
3. Apply to Species B (e.g., mouse) WITHOUT curated references
4. Fit normalizer on unlabeled Species B data to correct for species differences

### The Problem: Same-Species Testing

When testing on the **same species** (train on human, test on human):
1. Using `--pretrained_model`: Applies cross-species adaptation (WRONG for same-species!)
2. Using `-C` or full optimization: Uses reference normalization (CORRECT)

**This explains the discrepancy!**

## Solution Options

### Option 1: Document the Behavior (Quick Fix)

Add clear documentation that:
- `--pretrained_model` is for **cross-species** application only
- For same-species testing with a saved model, use `-C` with the optimized value
- The two modes will give different results by design

### Option 2: Add a Flag to Control Normalization Strategy (Medium Fix)

Add `--normalization_mode` argument:
- `cross_species`: Fit on experimental (current pretrained behavior)
- `same_species`: Fit on reference (current normal behavior)
- Auto-detect based on whether using pretrained model

### Option 3: Save Normalizer with Model (Larger Refactor)

Modify model saving to include the fitted normalizer:
```python
# When saving model
model_bundle = {
    'ensemble': result.ensemble,
    'threshold': config.scoring.threshold,
    'normalizer': normalizer  # ← Add fitted normalizer
}
```

Then during pretrained inference:
```python
# Load and use saved normalizer (same-species mode)
normalizer = model_data['normalizer']
normalized_introns = list(normalizer.transform(introns, dataset_type='experimental'))
```

### Option 4: Hybrid Approach (Most Flexible)

Save both ensemble and normalizer, but allow overriding:
- `--pretrained_model --cross_species`: Fit new normalizer on experimental
- `--pretrained_model` (default): Use saved normalizer from training

## Recommended Action

For immediate testing consistency:
1. **Avoid using `--pretrained_model` for same-species validation**
2. **Use `-C <optimized_value>` instead** when testing with a known good C parameter
3. This ensures reference normalization is used consistently

For long-term solution:
- Implement **Option 4** (save normalizer, allow override for cross-species)
- Clearly document the two use cases in help text and README

## Code Locations

**Pretrained normalization:**
- `cli/main.py:669-759` - `classify_with_pretrained_model()`
- Line 732: `normalizer.fit(introns, dataset_type='unlabeled')`

**Normal normalization:**
- `cli/main.py:618-666` - `normalize_scores()`
- Line 654: `normalizer.fit(reference_introns, dataset_type='reference')`

**Normalizer implementation:**
- `scoring/normalizer.py:160-413` - `ScoreNormalizer` class
- `scoring/normalizer.py:36-158` - `ZeroAnchoredRobustScaler` class

**Pipeline routing:**
- `cli/main.py:1060-1083` - If/else choosing pretrained vs normal path

---

**Last Updated:** 2025-11-10
**Analysis by:** Claude Code Investigation
