# Feature Augmentation Implementation Guide

**Date:** 2025-11-13
**Branch:** claude/false-positive-analysis-011CUt9YD4yNb4AzXF2TKAU9
**Goal:** Reduce false positives by augmenting features from 3 to 5

---

## Problem Statement

**Current Issue:** 8 false positive U12-type introns in *C. elegans* at >90% confidence (expected: ~0)

**Root Cause:** "One-end-strong" pattern - introns with a single high z-score but weak/negative scores elsewhere.

**Example extreme case:**
```
Intron: F49A5.10.1-intron_2
  5' SS z-score:  -0.15  (negative!)
  BPS z-score:     0.27  (weak)
  3' SS z-score:  11.45  (extremely high)
  → Classified as 100% U12-type (FALSE POSITIVE)
```

---

## Solution: Feature Augmentation (3→5 features)

### Current Features (3D)
```python
features = [z_5prime, z_bps, z_3prime]
```

### Proposed Features (5D)
```python
features = [
    z_5prime,              # Original
    z_bps,                 # Original
    z_3prime,              # Original
    abs(z_5prime - z_3prime),  # NEW: Asymmetry (penalize one-end-strong)
    z_5prime + z_3prime        # NEW: Sum_ends (reward consistency)
]
```

### Why This Works

**The 3D hyperplane approach:**
- Decision: `w1·z_5' + w2·z_BPS + w3·z_3' + b ≥ 0`
- Problem: A very high z_3' (like 11.45) can compensate for negative/weak other features
- The model *could* learn to penalize this implicitly, but it's fragile

**The 5D approach:**
- Explicitly represents the asymmetry property
- If model learns `w4 = -2.0` (negative weight for asymmetry):
  - Decision gets: `-2.0 × 11.6 = -23.2` penalty
  - This large penalty prevents false positive classification
- Still a linear model (interpretable, fast)
- Makes the problem more linearly separable

---

## Implementation Changes

### Overview

**3 files need modification:**

1. `classification/predictor.py` - Feature extraction during prediction
2. `classification/trainer.py` - Feature extraction during training
3. `classification/optimizer.py` - Feature extraction during hyperparameter optimization

---

## File 1: classification/predictor.py

### Location: Lines 368-372

**BEFORE:**
```python
features.append([
    intron.scores.five_z_score,
    intron.scores.bp_z_score,
    intron.scores.three_z_score
])
```

**AFTER:**
```python
# Extract z-scores for clarity
z5 = intron.scores.five_z_score
zbp = intron.scores.bp_z_score
z3 = intron.scores.three_z_score

features.append([
    z5,
    zbp,
    z3,
    abs(z5 - z3),  # Asymmetry: penalize one-end-strong
    z5 + z3         # Sum_ends: reward consistency
])
```

### Update Docstring: Line 345

**BEFORE:**
```python
Features: [five_z_score, bp_z_score, three_z_score]
```

**AFTER:**
```python
Features: [five_z_score, bp_z_score, three_z_score, asymmetry, sum_ends]
```

### Update Return Type: Line 354

**BEFORE:**
```python
Returns:
    Feature matrix (n_introns, 3)
```

**AFTER:**
```python
Returns:
    Feature matrix (n_introns, 5)
```

---

## File 2: classification/trainer.py

### Location 1: Lines ~285-290 (u12_features)

Find the code that builds u12_features (approximately line 285):

**BEFORE:**
```python
for intron in u12_introns:
    # Validation checks...

    u12_features.append([
        intron.scores.five_z_score,
        intron.scores.bp_z_score,
        intron.scores.three_z_score
    ])
```

**AFTER:**
```python
for intron in u12_introns:
    # Validation checks...

    z5 = intron.scores.five_z_score
    zbp = intron.scores.bp_z_score
    z3 = intron.scores.three_z_score

    u12_features.append([
        z5,
        zbp,
        z3,
        abs(z5 - z3),  # Asymmetry
        z5 + z3         # Sum_ends
    ])
```

### Location 2: Lines 296-300 (u2_features)

**BEFORE:**
```python
u2_features.append([
    intron.scores.five_z_score,
    intron.scores.bp_z_score,
    intron.scores.three_z_score
])
```

**AFTER:**
```python
z5 = intron.scores.five_z_score
zbp = intron.scores.bp_z_score
z3 = intron.scores.three_z_score

u2_features.append([
    z5,
    zbp,
    z3,
    abs(z5 - z3),  # Asymmetry
    z5 + z3         # Sum_ends
])
```

---

## File 3: classification/optimizer.py

### Location 1: Lines 423-427 (u12_features)

**BEFORE:**
```python
u12_features.append([
    intron.scores.five_z_score,
    intron.scores.bp_z_score,
    intron.scores.three_z_score
])
```

**AFTER:**
```python
z5 = intron.scores.five_z_score
zbp = intron.scores.bp_z_score
z3 = intron.scores.three_z_score

u12_features.append([
    z5,
    zbp,
    z3,
    abs(z5 - z3),  # Asymmetry
    z5 + z3         # Sum_ends
])
```

### Location 2: Lines 439-443 (u2_features)

**BEFORE:**
```python
u2_features.append([
    intron.scores.five_z_score,
    intron.scores.bp_z_score,
    intron.scores.three_z_score
])
```

**AFTER:**
```python
z5 = intron.scores.five_z_score
zbp = intron.scores.bp_z_score
z3 = intron.scores.three_z_score

u2_features.append([
    z5,
    zbp,
    z3,
    abs(z5 - z3),  # Asymmetry
    z5 + z3         # Sum_ends
])
```

---

## Complete Patch Summary

### Files Modified: 3
1. `classification/predictor.py` - 1 location + 2 docstring updates
2. `classification/trainer.py` - 2 locations (u12 and u2 feature lists)
3. `classification/optimizer.py` - 2 locations (u12 and u2 feature lists)

### Lines Changed: ~15-20 total

### Breaking Changes
- **Models incompatibility:** Old trained models (3 features) are incompatible with new code (5 features)
- **Solution:** Retrain models after implementing changes
- **Version:** Consider bumping to v1.6.0 to indicate algorithm change

---

## Testing Protocol

### 1. Baseline Test (H. sapiens Chr19)

**Purpose:** Ensure no regression on true U12 detection

```bash
# Run with 5-feature model
intronIC -g ../intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
         -a ../intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz \
         -n homo_sapiens_5feat

# Expected: ~30 U12-type introns (same as before)
awk '($2!="." && $2>0)' homo_sapiens_5feat.meta.iic | wc -l
```

**Success criteria:** 28-32 U12s found (within normal variance)

### 2. False Positive Test (C. elegans)

**Purpose:** Verify reduction in false positives

```bash
# Run on C. elegans
intronIC -g celegans_genome.fa.gz \
         -a celegans_annotation.gff3.gz \
         -n caenorhabditis_elegans_5feat

# Count U12s
awk '($2!="." && $2>0)' caenorhabditis_elegans_5feat.meta.iic | wc -l
```

**Success criteria:** 0-2 U12s (down from 8)

### 3. Reference Cross-Validation

**Purpose:** Ensure training quality

```bash
# Check metrics in training log
grep "F1 score" caenorhabditis_elegans_5feat.training.log
```

**Success criteria:** F1 ≥ 0.98 (should be similar to 3-feature version)

### 4. Comparative Analysis

```bash
# Run analysis script
python3 compare_true_vs_false_u12s.py
```

**Expected improvements:**
- False positive asymmetry: Reduced by 50%
- False positive sum_ends: Better separation from true U12s

---

## Expected Outcomes

| Metric | Before (3 feat) | After (5 feat) | Change |
|--------|-----------------|----------------|---------|
| **C. elegans FPs** | 8 | 0-2 | -75% to -100% |
| **H. sapiens U12s** | ~30 | ~30 | No change |
| **Reference F1** | 0.99+ | 0.98+ | Minimal |
| **Training time** | Baseline | +5-10% | Slight increase |
| **Prediction time** | Baseline | Negligible | No change |

---

## Model Weight Inspection (Optional)

After training, inspect learned weights to understand feature importance:

```python
def inspect_model_weights(ensemble):
    """Inspect learned SVM weights for 5-feature model."""
    feature_names = ['5\' SS', 'BPS', '3\' SS', 'Asymmetry', 'Sum_ends']

    print("\n" + "="*70)
    print("LEARNED MODEL WEIGHTS")
    print("="*70)

    for i, model_wrapper in enumerate(ensemble.models):
        # Get base SVM from calibrated wrapper
        base_model = model_wrapper.model.calibrated_classifiers_[0].estimator

        # For Pipeline, extract SVM step
        if hasattr(base_model, 'named_steps'):
            svm = base_model.named_steps['svm']
        else:
            svm = base_model

        weights = svm.coef_[0]  # Shape: (5,)
        intercept = svm.intercept_[0]

        print(f"\nModel {i+1}:")
        print(f"  Intercept:    {intercept:>8.3f}")
        for name, w in zip(feature_names, weights):
            print(f"  {name:>12}: {w:>8.3f}")

        # Highlight dominant features
        max_abs_w = max(abs(weights))
        dominant_idx = [j for j, w in enumerate(weights) if abs(w) == max_abs_w][0]
        print(f"\n  → Dominant: {feature_names[dominant_idx]} (|w|={max_abs_w:.3f})")

        # Check asymmetry weight (should be negative)
        asym_weight = weights[3]
        if asym_weight < 0:
            print(f"  ✓ Asymmetry weight is negative ({asym_weight:.3f}) - penalizes one-end-strong")
        else:
            print(f"  ⚠ Asymmetry weight is positive ({asym_weight:.3f}) - unexpected!")

# Usage:
# After training in trainer.py or in a diagnostic script:
# inspect_model_weights(ensemble)
```

**Expected observations:**
- `w_asymmetry` should be **negative** (penalizes high asymmetry)
- `w_sum_ends` may be **positive** (rewards both ends strong)
- `w_BPS` often largest in magnitude (BPS is strongest U12 signal)

---

## Rollback Plan

If results are worse than expected:

1. **Revert code changes:**
   ```bash
   git checkout HEAD~1 classification/predictor.py
   git checkout HEAD~1 classification/trainer.py
   git checkout HEAD~1 classification/optimizer.py
   ```

2. **Alternative approaches:**
   - Try only asymmetry feature (4 features total)
   - Try different feature combinations (e.g., min(z_5', z_3'))
   - Implement post-filter instead (see FALSE_POSITIVE_REDUCTION_PROPOSAL.md Phase 3)

---

## Documentation Updates

After successful implementation:

1. **README.md:** Update classification method section to mention 5 features
2. **CHANGELOG.md:** Add entry for v1.6.0 with feature augmentation
3. **Output format docs:** Note that old 3-feature models are incompatible

---

## Questions & Answers

### Q: Will this slow down the classifier?

**A:** Negligibly. The extra computations (`abs()` and `+`) are trivial. Training time may increase 5-10% due to higher dimensionality, but prediction time is essentially unchanged.

### Q: Why not use a non-linear kernel instead?

**A:** Non-linear kernels (RBF, polynomial) are slower, harder to interpret, and more prone to overfitting with imbalanced data. Linear models with good features are often better for this problem.

### Q: What if the 5-feature model has lower F1 on reference data?

**A:** A small drop (0.99 → 0.98) is acceptable if false positives are significantly reduced. The reference set is optimized for the old 3-feature space. Real-world performance (fewer FPs) is more important.

### Q: Can I use this with pretrained models?

**A:** No. Pretrained models expect 3 features. You must retrain after implementing 5-feature extraction. Consider creating new pretrained models and versioning them clearly (e.g., `v1.6_5feat`).

---

## References

- **Analysis script:** `analyze_false_positives.py`
- **Comparison script:** `compare_true_vs_false_u12s.py`
- **Full proposal:** `FALSE_POSITIVE_REDUCTION_PROPOSAL.md`
- **C. elegans results:** `archive/test_outputs/celegans/caenorhabditis_elegans.score_info.iic`
- **H. sapiens results:** `homo_sapiens.score_info.iic`

---

## Implementation Checklist

- [ ] Modify `classification/predictor.py` (1 location + docstrings)
- [ ] Modify `classification/trainer.py` (2 locations)
- [ ] Modify `classification/optimizer.py` (2 locations)
- [ ] Test on H. sapiens Chr19 (expect ~30 U12s)
- [ ] Test on C. elegans (expect 0-2 U12s, down from 8)
- [ ] Check reference cross-validation F1 (expect ≥0.98)
- [ ] Run comparative analysis
- [ ] Inspect model weights (optional)
- [ ] Update documentation
- [ ] Create new pretrained models (optional)
- [ ] Version bump to v1.6.0

---

**Status:** Ready for implementation
**Estimated time:** 1-2 hours for code changes + testing
**Expected impact:** 75-100% reduction in false positives
