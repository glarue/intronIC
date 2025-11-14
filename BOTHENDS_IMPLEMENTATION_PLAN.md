# BothEndsStrong Feature Augmentation - Implementation Plan

**Date:** 2025-11-14
**Branch:** claude/fix-scaler-centering-01C2BuWBX24F7n3cCBX1QUWu
**Goal:** Reduce false positives from 6 to ~0 in C. elegans using γ-weighted both-ends-strong features

---

## Summary

Implement pairwise "both-ends-strong" features with γ-weighted imbalance penalties to reduce "one-end-strong" false positives while maintaining linear SVM interpretability.

**Key Insight:** γ is a **feature-scaling hyperparameter** (not a loss weight) tuned via CV that changes effective regularization, making it "cheaper" for the optimizer to penalize imbalance.

---

## Architecture

### Feature Pipeline

```
Raw Intron
  ↓
Extract Base Features (3D)
  [s5, sBP, s3]  ← z-scores (LLRs)
  ↓
Scale (NO centering)
  RobustScaler(with_centering=False)
  ↓
BothEndsStrongTransformer(gamma_5_bp=γ₁, gamma_5_3=γ₂)
  Builds 7 features:
    1. s5           (original)
    2. sBP          (original)
    3. s3           (original)
    4. sum_5_bp     = s5 + sBP
    5. absdiff_5_bp = |s5 - sBP| × γ₁  ← SCALED imbalance penalty
    6. sum_5_3      = s5 + s3
    7. absdiff_5_3  = |s5 - s3| × γ₂   ← SCALED imbalance penalty
  ↓
LinearSVC(class_weight='balanced')
  ↓
CalibratedClassifierCV(method='sigmoid'|'isotonic')
  ↓
Calibrated probabilities
```

### Why This Works

**Without γ-weighting:**
- All features have similar scale after RobustScaler
- L2 regularization applies uniform penalty
- Model may ignore imbalance features

**With γ-weighting (γ > 1):**
- `absdiff_5_bp × γ` makes this feature larger
- L2 penalty now makes it "cheaper" to use this feature
- Model learns to heavily penalize 5'-BP imbalance
- Strong 5'-BP correlation is explicitly rewarded

---

## Implementation Steps

### 1. Create BothEndsStrongTransformer

**File:** `classification/transformers.py` (NEW)

```python
"""
Custom transformers for feature augmentation.
"""

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin


class BothEndsStrongTransformer(BaseEstimator, TransformerMixin):
    """
    Augment 3D features (5'SS, BPS, 3'SS z-scores) with pairwise
    both-ends-strong features.

    This transformer builds composite features that reward joint strength
    (sum) and penalize imbalance (abs difference) for correlated regions.

    The γ parameters control how strongly to weight the imbalance penalties.
    Higher γ makes the model more sensitive to one-end-strong patterns.

    Port from: Expert recommendations for reducing false positives

    Attributes:
        gamma_5_bp: Weight for 5'SS-BPS imbalance penalty (default: 1.0)
                   Suggested range for tuning: [1, 2, 4, 8]
        gamma_5_3: Weight for 5'SS-3'SS imbalance penalty (default: 1.0)
                  Suggested range for tuning: [1, 2, 4]
        include_min: Whether to include min(s5, sBP, s3) feature (default: False)
        include_total_imbalance: Whether to include total imbalance (default: False)

    Input features (3D):
        - s5:  5' splice site z-score (LLR)
        - sBP: Branch point z-score (LLR)
        - s3:  3' splice site z-score (LLR)

    Output features (7D or more):
        - s5, sBP, s3 (original, passed through)
        - sum_5_bp = s5 + sBP
        - absdiff_5_bp = |s5 - sBP| × gamma_5_bp
        - sum_5_3 = s5 + s3
        - absdiff_5_3 = |s5 - s3| × gamma_5_3
        - [optional] min_all = min(s5, sBP, s3)
        - [optional] imbalance_all = |s5-sBP| + |s5-s3| + |sBP-s3|

    Example:
        >>> transformer = BothEndsStrongTransformer(gamma_5_bp=4, gamma_5_3=2)
        >>> X = np.array([[2.0, 1.5, -0.5]])  # One intron: s5=2, sBP=1.5, s3=-0.5
        >>> X_aug = transformer.transform(X)
        >>> # X_aug = [2.0, 1.5, -0.5, 3.5, 2.0, 1.5, 5.0]
        >>> #          [s5,  sBP, s3,  sum, abs*4, sum, abs*2]
    """

    def __init__(
        self,
        gamma_5_bp: float = 1.0,
        gamma_5_3: float = 1.0,
        include_min: bool = False,
        include_total_imbalance: bool = False
    ):
        self.gamma_5_bp = gamma_5_bp
        self.gamma_5_3 = gamma_5_3
        self.include_min = include_min
        self.include_total_imbalance = include_total_imbalance

    def fit(self, X, y=None):
        """Fit transformer (no-op, stateless)."""
        return self

    def transform(self, X):
        """
        Augment 3D features to 7+D.

        Args:
            X: Array of shape (n_samples, 3) with [s5, sBP, s3]

        Returns:
            Array of shape (n_samples, 7+) with augmented features
        """
        # Ensure input is numpy array
        X = np.asarray(X)

        # Validate input shape
        if X.shape[1] != 3:
            raise ValueError(
                f"BothEndsStrongTransformer expects 3 input features "
                f"[s5, sBP, s3], got {X.shape[1]}"
            )

        # Extract base features
        s5 = X[:, 0]
        sBP = X[:, 1]
        s3 = X[:, 2]

        # Build composite features
        sum_5_bp = s5 + sBP
        absdiff_5_bp = np.abs(s5 - sBP) * self.gamma_5_bp
        sum_5_3 = s5 + s3
        absdiff_5_3 = np.abs(s5 - s3) * self.gamma_5_3

        # Stack features: originals + pairwise
        features = [
            s5[:, np.newaxis],
            sBP[:, np.newaxis],
            s3[:, np.newaxis],
            sum_5_bp[:, np.newaxis],
            absdiff_5_bp[:, np.newaxis],
            sum_5_3[:, np.newaxis],
            absdiff_5_3[:, np.newaxis]
        ]

        # Optional: min of all three
        if self.include_min:
            min_all = np.minimum(np.minimum(s5, sBP), s3)
            features.append(min_all[:, np.newaxis])

        # Optional: total imbalance
        if self.include_total_imbalance:
            imbalance_all = (
                np.abs(s5 - sBP) +
                np.abs(s5 - s3) +
                np.abs(sBP - s3)
            )
            features.append(imbalance_all[:, np.newaxis])

        return np.hstack(features)

    def get_feature_names_out(self, input_features=None):
        """Get output feature names for display."""
        names = ['s5', 'sBP', 's3', 'sum_5_bp', 'absdiff_5_bp', 'sum_5_3', 'absdiff_5_3']
        if self.include_min:
            names.append('min_all')
        if self.include_total_imbalance:
            names.append('imbalance_all')
        return np.array(names)
```

---

### 2. Update trainer.py Pipeline

**File:** `classification/trainer.py`

**Location:** Lines 178-192 (Pipeline definition)

**BEFORE:**
```python
base_pipeline = Pipeline([
    ('scale', RobustScaler(with_centering=False, with_scaling=True)),
    ('svc', LinearSVC(
        loss='squared_hinge',
        class_weight='balanced',
        dual=True,
        max_iter=20000,
        tol=3e-4,
        random_state=0
    ))
])
```

**AFTER:**
```python
from classification.transformers import BothEndsStrongTransformer

base_pipeline = Pipeline([
    ('scale', RobustScaler(with_centering=False, with_scaling=True)),
    ('augment', BothEndsStrongTransformer(
        gamma_5_bp=1.0,  # Will be tuned by optimizer
        gamma_5_3=1.0    # Will be tuned by optimizer
    )),
    ('svc', LinearSVC(
        loss='squared_hinge',
        class_weight='balanced',
        dual=True,
        max_iter=20000,
        tol=3e-4,
        random_state=0
    ))
])
```

**Note:** Feature extraction (lines 275-305) stays the same - still extracts 3 base features.

---

### 3. Update optimizer.py for γ Tuning

**File:** `classification/optimizer.py`

**Location 1:** Lines ~60-90 (compute_effective_c_range)

Add helper for γ-aware C range:

```python
def compute_effective_c_range(
    n_total: int,
    n_positive: int,
    c_eff_pos_min: float = 1e-3,
    c_eff_pos_max: float = 30.0,
    n_points: int = 15
) -> np.ndarray:
    """
    Compute C values that correspond to desired effective penalty range
    for the positive class.

    With class_weight='balanced', sklearn sets:
        w_pos = n_total / (2 * n_positive)
        w_neg = n_total / (2 * n_negative)

    The effective penalty on positive class is: C_eff_pos = C * w_pos

    For ~1% positives (typical U12 prevalence):
        w_pos ≈ 50 (for n_total=10000, n_pos=100)
        C_eff_pos ∈ [1e-3, 30] → C ∈ [2e-5, 0.6]

    Args:
        n_total: Total number of samples
        n_positive: Number of positive samples
        c_eff_pos_min: Minimum effective C for positive class
        c_eff_pos_max: Maximum effective C for positive class
        n_points: Number of C values to try

    Returns:
        Array of C values (log-spaced)
    """
    n_negative = n_total - n_positive
    w_pos = n_total / (2 * n_positive)

    # C_eff_pos = C * w_pos → C = C_eff_pos / w_pos
    c_min = c_eff_pos_min / w_pos
    c_max = c_eff_pos_max / w_pos

    return np.logspace(np.log10(c_min), np.log10(c_max), n_points)
```

**Location 2:** Lines 466-530 (optimize_linear_svc function)

**BEFORE:**
```python
# Build pipeline
base_pipeline = Pipeline([
    ('scale', RobustScaler(with_centering=False, with_scaling=True)),
    ('svc', LinearSVC(...))
])

# Calibrated wrapper
calibrated = CalibratedClassifierCV(
    estimator=base_pipeline,
    method='sigmoid',
    cv=5
)

# Parameter grid
param_grid = {
    'estimator__svc__C': np.logspace(-4, 2, 13)
}
```

**AFTER:**
```python
from classification.transformers import BothEndsStrongTransformer

# Build pipeline with augmentation
base_pipeline = Pipeline([
    ('scale', RobustScaler(with_centering=False, with_scaling=True)),
    ('augment', BothEndsStrongTransformer()),
    ('svc', LinearSVC(
        loss='squared_hinge',
        class_weight='balanced',
        dual=True,
        max_iter=20000,
        tol=3e-4,
        random_state=0
    ))
])

# Calibrated wrapper
calibrated = CalibratedClassifierCV(
    estimator=base_pipeline,
    method='sigmoid',  # Will be tuned
    cv=5
)

# Compute weight-aware C range
n_total = len(X_train) + len(X_test)
n_positive = sum(y_train) + sum(y_test)
c_values = compute_effective_c_range(
    n_total=n_total,
    n_positive=n_positive,
    c_eff_pos_min=1e-3,
    c_eff_pos_max=30.0,
    n_points=12
)

# Parameter grid with γ tuning
param_grid = {
    'estimator__svc__C': c_values,
    'estimator__augment__gamma_5_bp': [1, 2, 4, 8],     # 5'SS-BPS correlation strongest
    'estimator__augment__gamma_5_3': [1, 2, 4],         # 5'SS-3'SS secondary
    'method': ['sigmoid', 'isotonic']                   # Calibration method
}
```

**Update scoring:** Lines 540-545

```python
# Use neg_log_loss for probabilistic calibration
grid_search = GridSearchCV(
    estimator=calibrated,
    param_grid=param_grid,
    cv=5,
    scoring='neg_log_loss',  # Optimize calibrated probabilities
    n_jobs=cv_processes,
    verbose=1
)
```

---

### 4. Update predictor.py

**File:** `classification/predictor.py`

**No changes needed!** The predictor loads the trained pipeline which includes the BothEndsStrongTransformer. It still extracts 3 base features, and the transformer in the pipeline augments them to 7.

---

### 5. Testing Protocol

#### Test 1: H. sapiens Chr19 (Baseline)

```bash
intronIC -g test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
         -a test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz \
         -n homo_sapiens_bothends

# Count U12s
awk '($2!="." && $2>0)' homo_sapiens_bothends.meta.iic | wc -l
# Expected: 28-32 (similar to baseline)

# Check log loss and F1 in training log
grep -E "(Log loss|F1 score)" homo_sapiens_bothends.training.log
```

**Success criteria:**
- U12 count: 28-32
- F1 score: ≥ 0.98
- Log loss: Lower than baseline

#### Test 2: C. elegans (False Positive Reduction)

```bash
intronIC -g celegans_genome.fa.gz \
         -a celegans_annotation.gff3.gz \
         -n caenorhabditis_elegans_bothends

# Count U12s at 90% threshold
awk '($2!="." && $2>=0.90)' caenorhabditis_elegans_bothends.meta.iic | wc -l
# Expected: 0-2 (down from 6)

# Check distribution of predictions
awk '{if (NR>1 && $2!=".") print $2}' caenorhabditis_elegans_bothends.meta.iic | \
  sort -n | tail -20
# Should show much lower probabilities for top candidates
```

**Success criteria:**
- Positive rate @ 0.90: 0-2 introns (ideally 0)
- Top probabilities: < 0.85 (weaker false positives)

#### Test 3: Model Inspection

Create diagnostic script to inspect learned weights:

```python
def inspect_bothends_weights(ensemble):
    """Inspect learned weights for BothEndsStrong model."""
    feature_names = ['s5', 'sBP', 's3', 'sum_5_bp', 'absdiff_5_bp', 'sum_5_3', 'absdiff_5_3']

    print("\n" + "="*80)
    print("BOTHENDS MODEL WEIGHTS")
    print("="*80)

    for i, model_wrapper in enumerate(ensemble.models):
        # Extract SVM from calibrated pipeline
        pipeline = model_wrapper.model.calibrated_classifiers_[0].estimator
        svm = pipeline.named_steps['svc']
        augmenter = pipeline.named_steps['augment']

        weights = svm.coef_[0]
        intercept = svm.intercept_[0]

        print(f"\nModel {i+1}:")
        print(f"  Gamma_5_bp: {augmenter.gamma_5_bp}")
        print(f"  Gamma_5_3:  {augmenter.gamma_5_3}")
        print(f"  Intercept:  {intercept:>8.3f}")

        for name, w in zip(feature_names, weights):
            print(f"  {name:>15}: {w:>8.3f}")

        # Check imbalance penalty direction
        w_absdiff_5bp = weights[4]
        w_absdiff_5_3 = weights[6]

        if w_absdiff_5bp < 0:
            print(f"  ✓ absdiff_5_bp weight is negative ({w_absdiff_5bp:.3f})")
        else:
            print(f"  ⚠ absdiff_5_bp weight is positive ({w_absdiff_5bp:.3f})")

        if w_absdiff_5_3 < 0:
            print(f"  ✓ absdiff_5_3 weight is negative ({w_absdiff_5_3:.3f})")
        else:
            print(f"  ⚠ absdiff_5_3 weight is positive ({w_absdiff_5_3:.3f})")
```

**Expected observations:**
- `gamma_5_bp` likely > `gamma_5_3` (5'-BP correlation stronger)
- `w_absdiff_5bp` and `w_absdiff_5_3` should be **negative** (penalize imbalance)
- `w_sum_5_bp` may be positive (reward joint strength)
- `w_sBP` often largest in magnitude (BPS strongest signal)

---

## Expected Outcomes

| Metric | Before | After BothEnds | Change |
|--------|--------|----------------|---------|
| **C. elegans FPs @ 0.90** | 6 | 0-2 | -67% to -100% |
| **H. sapiens U12s** | ~30 | ~30 | No change |
| **Reference F1** | 0.99+ | 0.98+ | Minimal |
| **Log loss (CV)** | Baseline | Lower | Improved calibration |
| **Training time** | Baseline | +10-20% | More hyperparameters |

---

## Breaking Changes

⚠️ **Model incompatibility:** Old 3-feature models cannot be used with new code. Must retrain.

⚠️ **Pipeline structure:** New pipeline includes `BothEndsStrongTransformer` step.

⚠️ **Hyperparameter space:** Optimizer now tunes 4 parameters (C, gamma_5_bp, gamma_5_3, method) instead of just C.

---

## Implementation Checklist

- [ ] Create `classification/transformers.py` with BothEndsStrongTransformer
- [ ] Add compute_effective_c_range() to optimizer.py
- [ ] Update trainer.py pipeline to include augmentation step
- [ ] Update optimizer.py to tune γ parameters
- [ ] Test on H. sapiens Chr19 (baseline validation)
- [ ] Test on C. elegans (FP reduction)
- [ ] Inspect model weights (diagnostic)
- [ ] Document in CHANGELOG.md
- [ ] Version bump to v2.1.0

---

## Rollback Plan

If results worse than expected:

1. **Quick revert:**
   ```bash
   git checkout HEAD~1 classification/transformers.py
   git checkout HEAD~1 classification/trainer.py
   git checkout HEAD~1 classification/optimizer.py
   ```

2. **Alternative approaches:**
   - Try gamma_5_bp only (ignore gamma_5_3)
   - Try min(s5, sBP) instead of sum/absdiff
   - Implement post-filter instead of feature engineering

---

**Status:** Ready for implementation
**Estimated time:** 2-3 hours for code + testing
**Expected impact:** 67-100% reduction in C. elegans false positives
