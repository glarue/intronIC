# Expert-Recommended Feature Augmentation Implementation

**Date:** 2025-11-13
**Expert Input:** Integrated consistency features with zero-anchored robust scaling
**Goal:** Eliminate false positives while preserving true U12 detection

---

## Overview

This guide implements expert-recommended feature augmentation to address the "one-end-strong" false positive problem in C. elegans (8 FPs at >90% confidence, expected ~0).

**Key Innovation:** Add consistency features that act as **AND-gates** - requiring all splice signals (5' SS, BPS, 3' SS) to be coordinately strong, not just one.

---

## Problem Analysis

### Current State (3 features)
```python
X = [z_5prime, z_bps, z_3prime]  # 3D feature space
decision = w1*z_5' + w2*z_BPS + w3*z_3' + b
```

**Issue:** Single very high score can compensate for weak/negative others.

**Example C. elegans false positive:**
```
z_5' = -0.15  (negative!)
z_BPS = 0.27  (weak)
z_3' = 11.45  (extreme!)
→ Classified as 100% U12 (FALSE POSITIVE)
```

### Proposed Solution (6-9 features)

Add **consistency features** that explicitly encode:
1. Overall signal strength (sum features)
2. Imbalance penalty (absolute difference features)
3. AND-gate behavior (min features)

---

## Feature Design

### Core Consistency Features (Recommended Minimal Set)

```python
def augment_features(z_5prime, z_bps, z_3prime, gamma=1.0):
    """
    Augment 3D z-scores with consistency features.

    Args:
        z_5prime, z_bps, z_3prime: Original z-scores
        gamma: Asymmetry penalty weight (tune via CV)

    Returns:
        6D feature vector
    """
    # Original features
    s5 = z_5prime
    sbp = z_bps
    s3 = z_3prime

    # Consistency features
    sum_53 = s5 + s3                      # Both-ends strength
    absdiff_53 = abs(s5 - s3) * gamma    # Imbalance penalty (weighted)
    min_53 = 0.5 * (sum_53 - abs(s5 - s3))  # AND-gate: min(s5, s3)

    return [s5, sbp, s3, sum_53, absdiff_53, min_53]
```

**Why these features work:**

1. **sum_53** rewards overall signal strength
   - True U12s: both ends strong → high sum
   - False positives: one end only → lower sum

2. **absdiff_53 * gamma** penalizes imbalance
   - True U12s: balanced → low |s5 - s3|
   - False positives: extreme asymmetry → high penalty
   - gamma weight amplifies this penalty (tunable)

3. **min_53** acts as AND-gate
   - Mathematical identity: `min(a, b) = 0.5*(a + b - |a - b|)`
   - Requires BOTH ends to be strong
   - Model can learn: `w_sum - w_absdiff ≈ w_min`

### Extended Features (Optional, More Comprehensive)

```python
def augment_features_extended(z_5prime, z_bps, z_3prime, gamma=1.0):
    """Extended feature set including BPS consistency."""
    s5 = z_5prime
    sbp = z_bps
    s3 = z_3prime

    # Core consistency (5' & 3')
    sum_53 = s5 + s3
    absdiff_53 = abs(s5 - s3) * gamma
    min_53 = 0.5 * (sum_53 - abs(s5 - s3))

    # BPS consistency (optional)
    min_5_bp = min(s5, sbp)               # 5' AND BPS
    min_3_bp = min(s3, sbp)               # 3' AND BPS
    min_all = min(s5, sbp, s3)            # ALL three
    sum_all = s5 + sbp + s3               # Total signal

    return [
        s5, sbp, s3,           # Original (3)
        sum_53, absdiff_53, min_53,  # 5'-3' consistency (3)
        min_5_bp, min_3_bp, min_all  # BPS consistency (3)
        # Total: 9 features
    ]
```

**Recommendation:** Start with **6-feature minimal set**, then test extended 9-feature set if needed.

---

## Zero-Anchored Robust (ZAR) Scaling

**Critical:** Use ZAR scaling instead of standard z-score normalization.

### Why ZAR?

1. **Zero-anchored:** Zero stays zero (0 = "U12≈U2"), unlike z-scores which center data
2. **Robust:** Uses median and quantile clipping, not mean/std
3. **Species-agnostic:** Computed per-dataset, prevents cross-species artifacts

### Implementation

```python
import numpy as np

def zar_scale(x, q=0.995, eps=1e-9):
    """
    Compute robust spread for zero-anchored scaling.

    Args:
        x: Array of scores
        q: Upper quantile for clipping (default 0.995)
        eps: Minimum scale to avoid division by zero

    Returns:
        Robust scale value
    """
    a = np.abs(x)
    cap = np.quantile(a, q) if len(a) > 100 else a.max()
    clipped = np.clip(a, 0, cap)
    return max(np.median(clipped), eps)

def zar_transform(X):
    """
    Zero-anchored robust scaling for feature matrix.

    Args:
        X: Feature matrix (n_samples, n_features)

    Returns:
        X_scaled: Scaled features
        scales: Scale values per feature (for inverse transform)
    """
    scales = np.array([zar_scale(X[:, j]) for j in range(X.shape[1])])
    X_scaled = X / scales
    return X_scaled, scales
```

### When to Apply

Apply ZAR scaling **AFTER** feature augmentation, **BEFORE** SVM training:

```python
# 1. Augment features
X_aug = augment_features(z_5prime, z_bps, z_3prime, gamma=gamma)

# 2. Apply ZAR scaling
X_scaled, scales = zar_transform(X_aug)

# 3. Train SVM on scaled features
model.fit(X_scaled, y)
```

---

## Complete Training Pipeline

```python
import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import StratifiedKFold, GridSearchCV

# ============================================================================
# Feature Augmentation Transformer
# ============================================================================

class FeatureAugmenter(BaseEstimator, TransformerMixin):
    """
    Augment 3D z-scores with consistency features.

    Parameters:
        gamma: Asymmetry penalty weight (default 1.0)
        mode: 'minimal' (6 feat) or 'extended' (9 feat)
    """
    def __init__(self, gamma=1.0, mode='minimal'):
        self.gamma = gamma
        self.mode = mode

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        """Transform (n, 3) → (n, 6) or (n, 9)"""
        s5, sbp, s3 = X[:, 0], X[:, 1], X[:, 2]

        # Core consistency features
        sum_53 = s5 + s3
        absdiff_53 = np.abs(s5 - s3) * self.gamma
        min_53 = 0.5 * (sum_53 - np.abs(s5 - s3))

        if self.mode == 'minimal':
            return np.column_stack([s5, sbp, s3, sum_53, absdiff_53, min_53])

        elif self.mode == 'extended':
            min_5_bp = np.minimum(s5, sbp)
            min_3_bp = np.minimum(s3, sbp)
            min_all = np.minimum(np.minimum(s5, sbp), s3)
            return np.column_stack([
                s5, sbp, s3,
                sum_53, absdiff_53, min_53,
                min_5_bp, min_3_bp, min_all
            ])

# ============================================================================
# ZAR Scaling Transformer
# ============================================================================

class ZARScaler(BaseEstimator, TransformerMixin):
    """Zero-Anchored Robust scaler."""

    def __init__(self, q=0.995, eps=1e-9):
        self.q = q
        self.eps = eps
        self.scales_ = None

    def fit(self, X, y=None):
        """Compute robust scales per feature."""
        self.scales_ = np.array([
            self._robust_scale(X[:, j])
            for j in range(X.shape[1])
        ])
        return self

    def transform(self, X):
        """Apply scaling."""
        if self.scales_ is None:
            raise ValueError("Must call fit() before transform()")
        return X / self.scales_

    def _robust_scale(self, x):
        """Compute robust spread for one feature."""
        a = np.abs(x)
        cap = np.quantile(a, self.q) if len(a) > 100 else np.max(a)
        clipped = np.clip(a, 0, cap)
        return max(np.median(clipped), self.eps)

# ============================================================================
# Training Pipeline
# ============================================================================

# Cross-validation strategy
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# Build pipeline
pipeline = Pipeline([
    ('augment', FeatureAugmenter(gamma=1.0, mode='minimal')),
    ('zar', ZARScaler(q=0.995)),
    ('svm', LinearSVC(
        class_weight='balanced',
        loss='squared_hinge',
        dual=True,  # More stable for n_samples > n_features
        penalty='l2',
        tol=3e-4,
        max_iter=20000,
        random_state=42
    ))
])

# Calibration wrapper (isotonic is more flexible than sigmoid)
calibrated = CalibratedClassifierCV(
    pipeline,
    method='isotonic',
    cv=5,
    ensemble='auto'
)

# ============================================================================
# Hyperparameter Grid Search
# ============================================================================

# Weight-aware C grid
# For ~1% positive class: effective C_pos in [1e-3, 30]
# Actual C ≈ C_eff / w_pos where w_pos ≈ 50 for 1% class
C_grid = np.logspace(-4.6, 0.1, num=15)  # [~2.5e-5, ~1.25]

# Gamma values for asymmetry penalty
gamma_grid = [1.0, 2.0, 4.0, 8.0]

# Calibration method
calib_methods = ['sigmoid', 'isotonic']

param_grid = {
    'estimator__svm__C': C_grid,
    'estimator__augment__gamma': gamma_grid,
    'method': calib_methods
}

# Grid search with log-loss scoring (optimizes probability quality)
grid_search = GridSearchCV(
    calibrated,
    param_grid=param_grid,
    scoring='neg_log_loss',  # Optimize probability calibration
    cv=cv,
    n_jobs=-1,
    refit=True,
    verbose=1
)

# Train
grid_search.fit(X_train, y_train)

# Best model and parameters
best_model = grid_search.best_estimator_
best_params = grid_search.best_params_

print(f"Best C: {best_params['estimator__svm__C']:.6f}")
print(f"Best gamma: {best_params['estimator__augment__gamma']}")
print(f"Best calibration: {best_params['method']}")

# ============================================================================
# Prediction
# ============================================================================

# Predict probabilities
proba = best_model.predict_proba(X_test)[:, 1]  # P(U12)

# Apply threshold (default 90%)
y_pred = (proba >= 0.90).astype(int)
```

---

## Implementation in intronIC

### Files to Modify

**Same 3 files as before, but with enhanced feature extraction:**

1. `classification/predictor.py`
2. `classification/trainer.py`
3. `classification/optimizer.py`

### Code Changes

#### 1. Add Feature Augmentation Function

**New file:** `classification/features.py`

```python
"""Feature augmentation for consistency-based classification."""

import numpy as np
from typing import Tuple

def augment_features(
    z_5prime: np.ndarray,
    z_bps: np.ndarray,
    z_3prime: np.ndarray,
    gamma: float = 1.0,
    mode: str = 'minimal'
) -> np.ndarray:
    """
    Augment 3D z-scores with consistency features.

    Adds features that reward coordinated strong signals across all splice
    regions and penalize "one-end-strong" patterns common in false positives.

    Args:
        z_5prime: 5' splice site z-scores (n,)
        z_bps: Branch point z-scores (n,)
        z_3prime: 3' splice site z-scores (n,)
        gamma: Asymmetry penalty weight (higher = stricter)
        mode: 'minimal' (6 features) or 'extended' (9 features)

    Returns:
        Augmented feature matrix (n, 6) or (n, 9)

    Features (minimal mode):
        [z_5', z_BPS, z_3', sum_53, absdiff_53*gamma, min_53]

    Features (extended mode):
        [z_5', z_BPS, z_3', sum_53, absdiff_53*gamma, min_53,
         min_5_bp, min_3_bp, min_all]
    """
    s5 = z_5prime
    sbp = z_bps
    s3 = z_3prime

    # Core consistency features
    sum_53 = s5 + s3
    absdiff_53 = np.abs(s5 - s3) * gamma
    min_53 = 0.5 * (sum_53 - np.abs(s5 - s3))

    if mode == 'minimal':
        return np.column_stack([s5, sbp, s3, sum_53, absdiff_53, min_53])

    elif mode == 'extended':
        min_5_bp = np.minimum(s5, sbp)
        min_3_bp = np.minimum(s3, sbp)
        min_all = np.minimum(np.minimum(s5, sbp), s3)
        return np.column_stack([
            s5, sbp, s3,
            sum_53, absdiff_53, min_53,
            min_5_bp, min_3_bp, min_all
        ])
    else:
        raise ValueError(f"Unknown mode: {mode}")


def zar_scale_single(x: np.ndarray, q: float = 0.995, eps: float = 1e-9) -> float:
    """
    Compute zero-anchored robust scale for one feature.

    Args:
        x: Feature values
        q: Upper quantile for clipping (default 0.995)
        eps: Minimum scale to avoid division by zero

    Returns:
        Robust scale value
    """
    a = np.abs(x)
    cap = np.quantile(a, q) if len(a) > 100 else np.max(a)
    clipped = np.clip(a, 0, cap)
    return max(np.median(clipped), eps)


def zar_transform(X: np.ndarray, q: float = 0.995) -> Tuple[np.ndarray, np.ndarray]:
    """
    Zero-anchored robust scaling for feature matrix.

    Scales each feature by its robust spread around zero (median of clipped
    absolute values). Unlike z-score normalization, this preserves zero as
    the "no signal" point.

    Args:
        X: Feature matrix (n_samples, n_features)
        q: Upper quantile for clipping (default 0.995)

    Returns:
        X_scaled: Scaled features (n_samples, n_features)
        scales: Scale values per feature (n_features,)
    """
    scales = np.array([zar_scale_single(X[:, j], q) for j in range(X.shape[1])])
    X_scaled = X / scales
    return X_scaled, scales
```

#### 2. Modify Feature Extraction in predictor.py

**Location:** `classification/predictor.py`, line ~368

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
from classification.features import augment_features

# Extract z-scores
z5_list = [i.scores.five_z_score for i in introns]
zbp_list = [i.scores.bp_z_score for i in introns]
z3_list = [i.scores.three_z_score for i in introns]

# Augment features (gamma and mode should match training)
# These should be stored in model metadata or config
gamma = 1.0  # Default, or load from model
mode = 'minimal'  # Or 'extended'

features = augment_features(
    np.array(z5_list),
    np.array(zbp_list),
    np.array(z3_list),
    gamma=gamma,
    mode=mode
)
```

**Note:** The trained model's pipeline already includes ZAR scaling, so don't apply it again during prediction.

#### 3. Modify Training in trainer.py

**Location:** `classification/trainer.py`, feature extraction

**Integration:** Use the `FeatureAugmenter` and `ZARScaler` classes in the training pipeline (see complete training pipeline above).

#### 4. Update Hyperparameter Grid in optimizer.py

**Add gamma to the grid search:**

```python
param_grid = {
    'estimator__svm__C': C_grid,
    'estimator__augment__gamma': [1.0, 2.0, 4.0, 8.0],
    'method': ['sigmoid', 'isotonic']
}
```

---

## Testing Protocol

### 1. Baseline (H. sapiens Chr19)

**Expected:** ~30 U12s (no regression)

```bash
intronIC -g test_data/Homo_sapiens.Chr19.fa.gz \
         -a test_data/Homo_sapiens.Chr19.gff3.gz \
         -n homo_sapiens_6feat
```

### 2. False Positive Test (C. elegans)

**Expected:** 0-1 U12s (down from 8)

```bash
intronIC -g celegans.fa.gz -a celegans.gff3.gz -n celegans_6feat
awk '($2>0)' celegans_6feat.meta.iic | wc -l
```

### 3. Model Inspection

```python
# After training, inspect learned weights
best_svm = grid_search.best_estimator_.calibrated_classifiers_[0].estimator.named_steps['svm']
weights = best_svm.coef_[0]

feature_names = ['5SS', 'BPS', '3SS', 'sum_53', 'absdiff_53', 'min_53']
for name, w in zip(feature_names, weights):
    print(f"{name:>12}: {w:>8.3f}")

# Check key hypotheses:
assert weights[4] < 0, "absdiff_53 should have negative weight (penalize asymmetry)"
assert weights[3] > 0 or weights[5] > 0, "sum_53 or min_53 should have positive weight"
```

---

## Expected Results

| Metric | 3-feature | 6-feature (minimal) | 9-feature (extended) |
|--------|-----------|---------------------|----------------------|
| **H. sapiens U12s** | ~30 | ~30 | ~30 |
| **C. elegans FPs** | 8 | 0-1 | 0 |
| **Reference F1** | 0.99 | 0.98-0.99 | 0.98-0.99 |
| **Training time** | 1x | 1.2x | 1.4x |

**Key improvement:** Near-complete elimination of false positives with minimal impact on true positive detection.

---

## Advantages of This Approach

1. **Linear model** - fast, interpretable, no overfitting risk
2. **Learned penalties** - gamma and weights optimized by CV, not hand-tuned
3. **Zero-preserving** - ZAR scaling maintains "no signal" = 0
4. **Species-agnostic** - ZAR computed per-dataset, handles cross-species variation
5. **Explicit AND-gates** - min features require all signals strong
6. **Tunable strictness** - gamma controls asymmetry penalty strength

---

## Implementation Checklist

- [ ] Create `classification/features.py` with augmentation functions
- [ ] Add `FeatureAugmenter` and `ZARScaler` transformers
- [ ] Update `trainer.py` to use new pipeline
- [ ] Update `optimizer.py` with gamma grid search
- [ ] Update `predictor.py` to use augmented features
- [ ] Test on H. sapiens Chr19 (~30 U12s expected)
- [ ] Test on C. elegans (0-1 U12s expected)
- [ ] Inspect model weights (verify negative absdiff, positive sum/min)
- [ ] Document gamma parameter in user-facing docs
- [ ] Version bump to 1.6.0

---

## References

- Original analysis: `FALSE_POSITIVE_REDUCTION_PROPOSAL.md`
- Basic implementation: `FEATURE_AUGMENTATION_IMPLEMENTATION.md`
- C. elegans results: `archive/test_outputs/celegans/caenorhabditis_elegans.score_info.iic`

---

**Status:** Ready for implementation
**Estimated time:** 2-3 hours (includes pipeline refactoring)
**Expected impact:** >95% false positive reduction
