# Original intronIC v1.5.1 Architecture Review

**Date:** 2025-11-17
**Purpose:** Understand why original v1.5.1 had ~1 FP on C. elegans vs current ~130 FPs
**Source:** `archive/v1.5.1/intronIC/intronIC.py` (6,112 lines)

---

## Executive Summary

The original intronIC v1.5.1 used a **remarkably simple architecture**:
- **Normalization:** StandardScaler (mean-centering + variance scaling)
- **Features:** Just 3 z-scores (five, BP, three) - **NO balance features**
- **No clipping, no saturation, no augmentation**
- **Result:** ~1 FP on C. elegans (0.001% FP rate)

Our current architecture with zero-anchored scaling + balance features achieves ~130 FPs (0.12% FP rate) - a **130x regression** in cross-species performance.

---

## Original v1.5.1 Architecture

### 1. Score Normalization (line 3731)

**Algorithm:** `sklearn.preprocessing.StandardScaler`
```python
score_scaler = preprocessing.StandardScaler().fit(scale_vector)
```

**StandardScaler performs:**
```python
z = (raw_score - mean(raw_scores)) / std(raw_scores)
```

**Properties:**
- **Centering:** YES (subtracts mean)
- **Scaling:** YES (divides by std)
- **Fitted on:** Reference introns only
- **Distribution:** Transforms to standard normal (mean=0, std=1)

### 2. Features (lines 5158-5172)

**Exactly 3 features:**
```python
score_label_table = {
    'five': 'five_z_score',
    'bp': 'bp_z_score',
    'three': 'three_z_score'
}
```

**NO balance features** - just the 3 z-scores directly from StandardScaler.

### 3. SVM Training (lines 5345-5428)

- Model: sklearn.svm.SVC
- Subsample U2s to match U12 count
- 80/20 train/test split
- 5-fold CV for hyperparameter C
- Ensemble via U2 subsampling

**No preprocessing in pipeline** - model receives z-scores directly.

---

## Current Architecture (v2.0.0)

### 1. Score Normalization

**ZeroAnchoredRobustScaler:**
```python
z = raw_score / median(|raw_scores|)
```

**Key difference:** NO centering (preserves semantic zero)

### 2. Features

**7+ features:**
- five_z_score, bp_z_score, three_z_score
- min_all, neg_absdiff_5_bp, neg_absdiff_5_3, neg_absdiff_bp_3

**Balance features added** to penalize imbalance.

---

## Why Original Was Better for Cross-Species

### Hypothesis: StandardScaler Centering Removes Composition Bias

**Original (with centering):**
```python
# Human training: mean_3=5.0, std_3=2.0
# C. elegans inference: raw_3=14.46
z_3 = (14.46 - 5.0) / 2.0 = +4.73�
```

**Centering removes the global shift** caused by composition bias. Feature becomes relative to training distribution.

**Current (no centering):**
```python
# Human training: median(|3'|)=1.26
# C. elegans inference: raw_3=14.46
z_3 = 14.46 / 1.26 = +11.45
```

**No centering = unbounded growth**. Composition bias directly inflates z-score.

### C. elegans Performance Comparison

| Architecture | False Positives | FP Rate | Change |
|--------------|----------------|---------|--------|
| Original v1.5.1 | ~1 | 0.001% | baseline |
| Current (no clip) | 130 | 0.12% | **130x worse** |

---

## Critical Mistake: Removing Centering

**Original reasoning:** "Preserve semantic zero (score=0 means U12HU2)"

**Problem:**
- Semantic zero only matters for **interpretation**, not for SVM
- SVM learns a hyperplane in feature space - doesn't care about absolute zero
- Without centering, features are **composition-dependent** and incomparable across species

**StandardScaler centering:**
- Creates **relative features** (relative to training distribution)
- **Accidentally solved** cross-species composition bias
- Standard practice in ML - there's a reason it's the default

---

## Recommendations (Expert-Validated)

### Architectural Principles

**Based on empirical evidence (original ~1 FP vs current 130 FPs) and expert guidance:**

1. **Don't use zero-anchoring in ML pipeline**
   - Nice for interpretability, but empirically fails cross-species
   - Centering gives better behavior in U12-free genomes

2. **Keep zero-anchored LLRs for humans, not models**
   - Use in: explanatory plots, method comparisons, exported annotations
   - Raw LLRs already encode semantic zero - that's enough

3. **Use domain-adaptive centered scaling**
   - StandardScaler OR RobustScaler (both valid, RobustScaler may handle outliers better)
   - Fit per species for domain adaptation
   - Centering is **essential**, not optional

4. **Optional: Quantile-based clipping**
   - Learn quantile q from human training
   - Recompute bounds per species (domain-adaptive)
   - Original had no clipping and worked - may not be needed

---

## Recommended Architecture

### Training (Human):
```python
Pipeline([
    ('scale', RobustScaler()),  # Or StandardScaler - both valid
    ('clip', SymmetricClipper(q=0.95, domain_adaptive=True)),  # Optional
    ('svc', LinearSVC(C=optimized))
])
```

**Features:** Start with 3 (five_z, BP_z, three_z) like original, optionally add balance features

### Inference (C. elegans, Domain-Adaptive):
```python
# Fit scaler on target species (unsupervised, no label leakage)
scaler.fit(celegans_raw_scores)

# Transform using target species parameters
z_scores = scaler.transform(celegans_raw_scores)

# Optional: Recompute clip bounds from target z-scores using learned q
if clipper:
    z_clipped = clipper.transform(z_scores)  # Recomputes bounds if domain_adaptive=True

# Classify with pretrained human SVM
predictions = svm.predict(z_clipped)
```

---

## Implementation Plan

### Phase 1: Test RobustScaler (Simplest Change)

**Hypothesis:** RobustScaler with centering will match original performance

```python
from sklearn.preprocessing import RobustScaler

Pipeline([
    ('scale', RobustScaler()),  # with_centering=True (default)
    ('svc', LinearSVC(C=optimized))
])
```

**Expected C. elegans:** <10 FPs (based on original ~1 FP with StandardScaler)

**Why RobustScaler?**
- More robust to outliers than StandardScaler (uses median/IQR instead of mean/std)
- Still has **centering** (the critical component)
- Better for data with extreme values (like our PWM scores)

### Phase 2: Add Domain Adaptation (If Needed)

If Phase 1 works, add per-species fitting:

```python
# At inference time, fit scaler on target species
scaler = RobustScaler()
scaler.fit(target_species_raw_scores)  # Unsupervised
z_scores = scaler.transform(target_species_raw_scores)
predictions = pretrained_svm.predict(z_scores)
```

### Phase 3: Optional Clipping (If FPs Remain)

```python
Pipeline([
    ('scale', RobustScaler()),
    ('clip', SymmetricClipper(quantile=0.95, domain_adaptive=True)),
    ('svc', LinearSVC(C=optimized))
])
```

But original had **zero clipping** and worked perfectly - probably won't need this.

---

## Conclusion

**The original was correct, we overcomplicated:**
- StandardScaler/RobustScaler: Standard practice, centering handles composition bias
- Simple features: Better generalization
- No zero-anchoring in ML: Interpretability ≠ model performance

**Our redesign mistake:**
- Confused interpretability (preserve semantic zero) with model needs (comparable features)
- Zero-anchoring broke cross-species by removing centering
- Added complexity (balance features, clipping) to fix symptoms instead of root cause

**The fix:**
1. **Immediate:** Return to centered scaling (RobustScaler or StandardScaler)
2. **Optional:** Add domain adaptation (fit per species)
3. **Test:** Validate on C. elegans - expect <10 FPs

**Evidence:** 130x regression correlates exactly with removing centering.

**Expert validation:** "If centering empirically gives you better behavior in U12-free genomes, that's the more important criterion."

---

## Next Steps

1. Modify training pipeline to use RobustScaler (with centering)
2. Train new model with 3 features only
3. Test on C. elegans
4. If successful (<10 FPs), document as the correct architecture
5. Then experiment with adding balance features back (with centering)