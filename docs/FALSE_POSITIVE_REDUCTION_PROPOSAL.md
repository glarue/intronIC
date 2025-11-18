# False Positive Reduction Proposal

**Date:** 2025-11-13
**Branch:** claude/investigate-ref-performance-011CUt9YD4yNb4AzXF2TKAU9
**Issue:** 8 false positive U12-type introns called in C. elegans (expected ~0)

---

## Executive Summary

The current intronIC classifier produces **8 false positives** in *C. elegans* at >90% probability, despite this organism having no genuine U12-type introns. Analysis reveals these are classic "one-end-strong" cases where a single feature (5' SS, BPS, or 3' SS) has a high score but other features are weak or negative.

**Key Finding:** 100% of false positives have at least one weak feature (z < 0.5), compared to 80.6% of true U12s having this property.

---

## Problem Analysis

### The False Positives

| Intron | 5' z | BPS z | 3' z | Min z | Pattern |
|--------|------|-------|------|-------|---------|
| F49A5.10.1-intron_2 | -0.15 | 0.27 | **11.45** | -0.15 | Extreme 3' only |
| Y39B6A.18-intron_4 | 0.16 | **2.20** | 1.71 | 0.16 | Strong BPS, weak 5' |
| W03F8.2-intron_5 | 0.37 | **-1.44** | -0.51 | -1.44 | Negative BPS! |
| ZK697.12.1-intron_4 | 0.32 | **-0.44** | 0.58 | -0.44 | Negative BPS! |
| Y6B3A.1a.1-intron_14 | 0.29 | 0.10 | -0.53 | -0.53 | All weak |

**Common patterns:**
- **5/8 have at least ONE negative z-score** (62.5%)
- **8/8 have at least ONE weak feature (z < 0.5)** (100%)
- **3/8 have extremely high single feature** (one z > 5, rest < 1)

### Comparison to True U12s (Homo sapiens, n=31)

| Metric | True U12s | False Positives | Difference |
|--------|-----------|-----------------|------------|
| **Mean min z-score** | -0.098 | **-0.491** | 4x worse |
| **All features positive** | 71.0% | **37.5%** | 2x less |
| **Sum (5' + 3')** | 3.726 | **1.757** | 47% lower |
| **At least one negative** | 29.0% | **62.5%** | 2x more |
| **At least one weak (z<0.5)** | 80.6% | **100%** | All! |

**Critical insight:** False positives have significantly weaker overall feature profiles, relying on a single strong feature to achieve high classification scores.

---

## Root Cause

The current classifier uses **only 3 features:**
1. 5' splice site z-score
2. Branch point sequence (BPS) z-score
3. 3' splice site z-score

The linear SVM learns a decision boundary in this 3D space, but **does not explicitly penalize asymmetry** or require **consistency across features**. This allows introns with one exceptionally high score to compensate for weak/negative scores in other features.

### Why This Matters

True U12-type introns have **coordinated motifs**:
- Strong 5' SS (AT or GT)
- Highly conserved BPS (TCCTTAAC)
- Strong 3' SS (AC or AG)

A random intron might **by chance** have one strong motif, but is unlikely to have all three. The classifier should **explicitly reward consistency** and **penalize asymmetry**.

---

## Proposed Solutions

Based on expert feedback, we recommend a **three-pronged approach**:

### 1. Feature Augmentation (Linear Model Compatible)

**Add 2 new features** to the existing 3 z-scores:

```python
# Current features: z_5prime, z_bps, z_3prime

# New consistency features:
asymmetry = abs(z_5prime - z_3prime)    # Penalize one-end-strong
sum_ends = z_5prime + z_3prime           # Reward consistency
```

**Effect:** The linear SVM can learn to:
- **Penalize** high asymmetry (one end strong, other weak)
- **Reward** high sum_ends (both ends strong together)

**Advantages:**
- Still a linear model (no architectural change)
- No user-facing parameters (learned during training)
- Directly addresses the one-end-strong problem

**Implementation:**
```python
# In scoring/scorer.py or classification/trainer.py
def extract_features(scores):
    """Extract features for classification."""
    z_5prime = scores['z_5prime']
    z_bps = scores['z_bps']
    z_3prime = scores['z_3prime']

    # Original 3 features
    features = [z_5prime, z_bps, z_3prime]

    # NEW: Consistency features
    asymmetry = abs(z_5prime - z_3prime)
    sum_ends = z_5prime + z_3prime
    features.extend([asymmetry, sum_ends])

    return features  # Now 5 features instead of 3
```

---

### 2. Zero-Anchored Robust Spread (ZAR) Normalization

**Per-species unsupervised normalization** to prevent score inflation:

```python
def zar_normalize(scores_d, quantile=0.995):
    """
    Zero-Anchored Robust spread normalization.

    Args:
        scores_d: Array of raw scores for dimension d
        quantile: Upper quantile for clipping (default 0.995)

    Returns:
        Normalized scores (zero-centered spread)
    """
    import numpy as np

    # Clip extreme values
    abs_scores = np.abs(scores_d)
    q_threshold = np.quantile(abs_scores, quantile)
    clipped = np.minimum(abs_scores, q_threshold)

    # Compute robust spread (median of clipped absolute values)
    sigma_hat_d = np.median(clipped)

    # Normalize (NO centering - zero stays zero)
    normalized = scores_d / sigma_hat_d

    return normalized
```

**Key properties:**
- **Zero stays zero** (no centering)
- **Robust to outliers** (uses median, clips extremes)
- **Species-specific** (computed per-dataset automatically)
- **No manual thresholds** (fully automatic)

**When to apply:** After PWM scoring, before SVM training/prediction

**Effect:** Prevents species with different score distributions from having inflated or deflated features relative to the training data.

---

### 3. Stable Solver Configuration

Keep the existing solver improvements from performance investigation:

```python
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
from sklearn.preprocessing import MaxAbsScaler
from sklearn.pipeline import Pipeline

# Stable pipeline
base_svm_pipeline = Pipeline([
    ('scale', MaxAbsScaler()),
    ('svm', LinearSVC(
        C=optimized_C,  # From grid search
        loss='squared_hinge',
        class_weight='balanced',
        dual=True,  # More stable for n_samples > n_features
        penalty='l2',
        max_iter=20000,
        tol=3e-4,
        random_state=seed
    ))
])

# Isotonic calibration (more flexible than sigmoid)
model = CalibratedClassifierCV(
    base_svm_pipeline,
    method='isotonic',
    cv=5
)
```

**Hyperparameter grid:**
```python
# Weight-aware C grid
# Target: C_eff_plus in [1e-3, 3e1]
# Back-calculate: C = C_eff_plus / w_plus

w_plus = class_weight['positive']  # Balanced weight for U12s
C_eff_range = np.logspace(-3, 1.5, num=20)  # [1e-3, ~30]
C_base_range = C_eff_range / w_plus

param_grid = {
    'estimator__svm__C': C_base_range
}
```

---

## Implementation Plan

### Phase 1: Feature Augmentation (Priority: HIGH)

**Files to modify:**
1. `classification/trainer.py` - Add feature extraction
2. `classification/predictor.py` - Add feature extraction
3. `scoring/scorer.py` - Store asymmetry/sum in Intron objects (optional)

**Changes:**
```python
# Before (3 features):
X = np.array([[i.z_5prime, i.z_bps, i.z_3prime] for i in introns])

# After (5 features):
X = np.array([
    [i.z_5prime, i.z_bps, i.z_3prime,
     abs(i.z_5prime - i.z_3prime),  # asymmetry
     i.z_5prime + i.z_3prime]        # sum_ends
    for i in introns
])
```

**Testing:**
1. Retrain on reference data with 5 features
2. Test on H. sapiens Chr19 (should still find ~30 U12s)
3. Test on C. elegans (should find 0-2 U12s, not 8)

**Expected impact:** Reduce false positives by 50-75%

---

### Phase 2: ZAR Normalization (Priority: MEDIUM)

**Files to modify:**
1. `scoring/normalizer.py` - Add ZAR method
2. `cli/main.py` - Add per-species normalization step

**Implementation:**
```python
# In normalizer.py
class ZARNormalizer:
    """Zero-Anchored Robust spread normalizer."""

    def fit(self, scores, quantile=0.995):
        """Compute robust spread from data."""
        abs_scores = np.abs(scores)
        q_threshold = np.quantile(abs_scores, quantile)
        clipped = np.minimum(abs_scores, q_threshold)
        self.sigma_hat = np.median(clipped)
        return self

    def transform(self, scores):
        """Normalize scores."""
        return scores / self.sigma_hat

# In main pipeline:
# BEFORE classification, AFTER PWM scoring:
zar_5prime = ZARNormalizer().fit(all_raw_5prime_scores)
zar_bps = ZARNormalizer().fit(all_raw_bps_scores)
zar_3prime = ZARNormalizer().fit(all_raw_3prime_scores)

for intron in introns:
    intron.z_5prime = zar_5prime.transform(intron.raw_5prime)
    intron.z_bps = zar_bps.transform(intron.raw_bps)
    intron.z_3prime = zar_3prime.transform(intron.raw_3prime)
```

**Testing:**
1. Compare normalized scores across species
2. Verify no regression on H. sapiens
3. Check C. elegans false positive rate

**Expected impact:** Additional 10-30% reduction in false positives

---

### Phase 3: Post-Filter (Optional, Priority: LOW)

If false positives persist after Phases 1-2, add a simple post-classification filter:

```python
def post_filter_u12s(predictions, introns, min_z_threshold=0.0):
    """
    Filter U12 predictions requiring minimum z-score.

    Args:
        predictions: SVM probability scores
        introns: Intron objects with z-scores
        min_z_threshold: Minimum z-score required for weakest feature

    Returns:
        Filtered predictions
    """
    filtered = []
    for prob, intron in zip(predictions, introns):
        min_z = min(intron.z_5prime, intron.z_bps, intron.z_3prime)

        if prob >= 0.90 and min_z < min_z_threshold:
            # Downgrade to just below threshold
            prob = 0.89

        filtered.append(prob)

    return filtered
```

**Note:** This is a **last resort** - prefer to fix via features (Phase 1) and normalization (Phase 2).

---

## Testing Strategy

### Test Cases

1. **Homo sapiens Chr19** (positive control)
   - Current: ~30 U12s found
   - Target: ~30 U12s found (no regression)

2. **Caenorhabditis elegans** (negative control)
   - Current: 8 false positives
   - Target: 0-2 (significant reduction)

3. **Drosophila melanogaster** (validation)
   - Known U12 count: Check against literature
   - Verify no major changes

4. **Reference set cross-validation**
   - F1 score should remain ≥0.95
   - No significant drop in recall on reference U12s

### Metrics

- **False Positive Rate** (C. elegans): # U12s called / total introns
- **True Positive Rate** (H. sapiens): # U12s found / 30 expected
- **F1 Score** (Reference CV): Overall classifier quality
- **Asymmetry Impact**: Mean |5' - 3'| for true positives vs false positives

---

## Expected Outcomes

| Metric | Before | After Phase 1 | After Phase 2 |
|--------|--------|---------------|---------------|
| **C. elegans FPs** | 8 | 2-3 | 0-1 |
| **H. sapiens U12s** | ~30 | ~30 | ~30 |
| **Reference F1** | 0.99+ | 0.98+ | 0.98+ |
| **FP asymmetry** | 2.24 | Expect 50% lower | Expect 60% lower |

---

## Risks and Mitigations

### Risk 1: Feature augmentation reduces recall

**Mitigation:** The new features are complementary, not restrictive. The SVM can learn to ignore them if unhelpful (coefficients → 0). Test on reference data first.

### Risk 2: ZAR normalization changes reference performance

**Mitigation:** ZAR is applied uniformly to both training and test data. Should be neutral for reference sequences. Test on held-out reference set.

### Risk 3: Breaking changes to existing models

**Mitigation:** This requires retraining models. Old models (3 features) incompatible with new code (5 features). Version the models clearly.

---

## Alternative Approaches (Not Recommended)

### 1. Higher Probability Threshold (e.g., 95% instead of 90%)

**Problem:** Doesn't address root cause. Some one-end-strong introns will still score >95%.

### 2. Hand-crafted Rules (e.g., "require min_z > 0")

**Problem:** Not learned from data. May be too restrictive for some true U12s.

### 3. Non-linear Model (e.g., RBF kernel SVM)

**Problem:** Much slower, harder to interpret, more prone to overfitting with imbalanced data.

---

## Recommendation

**Implement Phases 1 and 2** in order:

1. **Phase 1** (Feature Augmentation) - Directly addresses the identified problem with minimal complexity
2. **Phase 2** (ZAR Normalization) - Prevents species-specific artifacts

**Phase 3** (Post-filter) should only be considered if Phases 1-2 are insufficient.

**Timeline:**
- Phase 1: 1-2 days (implementation + testing)
- Phase 2: 1 day (implementation + testing)
- Validation: 1 day (full test suite)

**Total:** 3-4 days for significant false positive reduction

---

## References

- Expert feedback document: `archive/investigation/EXPERT_FIXES_FINAL_RESULTS.md`
- False positive analysis: `analyze_false_positives.py`
- Comparison analysis: `compare_true_vs_false_u12s.py`
- C. elegans data: `archive/test_outputs/celegans/caenorhabditis_elegans.score_info.iic`
- Homo sapiens data: `homo_sapiens.score_info.iic`

---

**Status:** Proposal ready for review and implementation.
