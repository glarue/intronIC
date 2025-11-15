# False Positive Reduction - Detailed Implementation Guide
## Tiers 1 & 2 Implementation

**Goal:** Reduce C. elegans FP rate by 35-75% while maintaining recall, with species-agnostic improvements.

**C Bounds Correction:** Expert recommended `(1e-3, 1e+1)` but we'll use `(3e-4, 1e+2)` for slightly more range on upper bound.

---

## Implementation Roadmap

### Tier 1: Core Feature & Regularization Improvements

1. ✅ Add `min_all` feature (3-way AND)
2. ✅ Add `neg_absdiff_*` penalty features
3. ✅ Update default C bounds to `(3e-4, 1e+2)`
4. ✅ Add L1 penalty option to param_grid

### Tier 2: Advanced Improvements

5. ✅ Implement `OutlierClipper` transformer
6. ✅ Integrate `OutlierClipper` into pipeline
7. ✅ Fix missing `ensemble=True` in optimizer.py
8. ✅ Document sigmoid calibration preference

---

## Detailed Implementation Steps

### Step 1: Add min_all Feature

**File:** `classification/transformers.py`

**Location:** `BothEndsStrongTransformer.transform()` method

**Current feature set (include_max=False):**
```python
[s5, sBP, s3, min_5_bp, min_5_3]  # 5D
```

**New feature set (include_max=False):**
```python
[s5, sBP, s3, min_5_bp, min_5_3, min_all]  # 6D
```

**Implementation:**
```python
# After computing min_5_bp and min_5_3:

# 3-way AND: all three must be strong
# min_all = min(s5, sBP, s3)
min_all = np.minimum(np.minimum(s5, sBP), s3)
```

**Why it works:**
- U12 introns have strong 5'SS, BPS, AND 3'SS signals
- FPs often have only 2/3 strong (e.g., strong 5'SS+BPS, weak 3'SS)
- min_all = low value → likely FP
- min_all = high value → all three strong → likely true U12

**Add to features list:**
```python
features = [
    s5[:, np.newaxis],
    sBP[:, np.newaxis],
    s3[:, np.newaxis],
    min_5_bp[:, np.newaxis],
    min_5_3[:, np.newaxis],
    min_all[:, np.newaxis]  # NEW
]
```

**Update feature names:**
```python
names = [
    's5',
    'sBP',
    's3',
    'min_5_bp',
    'min_5_3',
    'min_all'  # NEW
]
```

**Update docstring:**
- Change output features from 5D/7D to 6D/8D
- Add example showing min_all computation
- Document that min_all requires ALL three signals to be strong

---

### Step 2: Add neg_absdiff Penalty Features

**File:** `classification/transformers.py`

**Location:** `BothEndsStrongTransformer.transform()` method

**Current computation (already exists for min/max):**
```python
absdiff_5_bp = np.abs(s5 - sBP)
absdiff_5_3 = np.abs(s5 - s3)
```

**New features:**
```python
# Negated absdiff = penalty for imbalance
# Model puts positive weight on these → rewards consistency
neg_absdiff_5_bp = -absdiff_5_bp
neg_absdiff_5_3 = -absdiff_5_3

# Optional: BPS-3'SS imbalance penalty
absdiff_bp_3 = np.abs(sBP - s3)
neg_absdiff_bp_3 = -absdiff_bp_3
```

**Why it works:**
- Large |s5 - sBP| → one-end-strong → FP
- neg_absdiff = large negative value → penalizes FP
- Small |s5 - sBP| → balanced signals → U12
- neg_absdiff = small negative value → doesn't penalize

**Feature set evolution:**

Without neg_absdiff_bp_3 (simpler, start here):
```python
[s5, sBP, s3, min_5_bp, min_5_3, min_all, neg_absdiff_5_bp, neg_absdiff_5_3]  # 8D
```

With neg_absdiff_bp_3 (more complete):
```python
[s5, sBP, s3, min_5_bp, min_5_3, min_all, neg_absdiff_5_bp, neg_absdiff_5_3, neg_absdiff_bp_3]  # 9D
```

**Recommendation:** Include all three neg_absdiff features for maximum FP reduction.

**Add to features list:**
```python
features = [
    s5[:, np.newaxis],
    sBP[:, np.newaxis],
    s3[:, np.newaxis],
    min_5_bp[:, np.newaxis],
    min_5_3[:, np.newaxis],
    min_all[:, np.newaxis],
    neg_absdiff_5_bp[:, np.newaxis],  # NEW
    neg_absdiff_5_3[:, np.newaxis],   # NEW
    neg_absdiff_bp_3[:, np.newaxis]   # NEW
]
```

**Update feature names:**
```python
names = [
    's5',
    'sBP',
    's3',
    'min_5_bp',
    'min_5_3',
    'min_all',
    'neg_absdiff_5_bp',  # NEW
    'neg_absdiff_5_3',   # NEW
    'neg_absdiff_bp_3'   # NEW
]
```

**Final dimensions:**
- include_max=False: 9D (base 3 + min 3 + neg_absdiff 3)
- include_max=True: 11D (+ max_5_bp + max_5_3)

**Update docstring:**
- Document neg_absdiff features and their purpose
- Add example showing imbalance penalty calculation
- Explain why negative values penalize one-end-strong patterns

---

### Step 3: Update Default C Bounds

**File:** `classification/classifier.py`

**Location:** `IntronClassifier.__init__()` default parameter

**Current:**
```python
def __init__(
    self,
    # ... other params ...
    eff_C_pos_range: tuple = (1e-3, 1e3),
    # ...
):
```

**New:**
```python
def __init__(
    self,
    # ... other params ...
    eff_C_pos_range: tuple = (3e-4, 1e+2),  # Tighter bounds reduce FPR
    # ...
):
```

**Rationale:**
- Larger C → less regularization → more complex boundary → higher FPR
- Shrinking upper bound from 1e3 to 1e2 reduces FP tendency
- Lower bound from 1e-3 to 3e-4 allows slightly more regularization
- Expert: "larger C tends to raise FPR"

**Update docstring:**
```python
"""
Args:
    ...
    eff_C_pos_range: Target effective penalty range for positive class
                     (default: (3e-4, 1e+2), tightened to reduce FPR)
    ...
"""
```

**Note:** Users can still override via config file's `c_bounds:` section.

---

### Step 4: Add L1 Penalty Option

**File:** `classification/optimizer.py`

**Location:** `_run_optimization_round()` method, param_grid definition (line ~527)

**Current param_grid:**
```python
param_grid = {
    'estimator__svc__C': C_grid,
    'estimator__augment__include_max': [False, True],
    'estimator__svc__dual': [False, True],
    'estimator__svc__intercept_scaling': [10.0, 100.0, 1000.0],
    'method': ['sigmoid', 'isotonic']
}
```

**New param_grid:**
```python
param_grid = {
    'estimator__svc__C': C_grid,
    'estimator__augment__include_max': [False, True],
    'estimator__svc__dual': [False],  # L1 requires dual=False
    'estimator__svc__penalty': ['l1', 'l2'],  # NEW
    'estimator__svc__loss': ['squared_hinge'],  # L1 requires squared_hinge
    'estimator__svc__intercept_scaling': [10.0, 100.0, 1000.0],
    'method': ['sigmoid', 'isotonic']
}
```

**CRITICAL CONSTRAINTS:**
- L1 penalty requires `dual=False`
- L1 penalty requires `loss='squared_hinge'`
- We'll simplify by forcing dual=False for all (our n_samples >> n_features anyway)

**Why it works:**
- With 9-11 features (many correlated), L1 prunes redundant ones
- Features that re-introduce one-end behavior get zeroed out
- Results in sparser, more interpretable models
- Often improves generalization

**Parameter count change:**
- Old: C × include_max × dual × intercept × method = C × 2 × 2 × 3 × 2 = 24C
- New: C × include_max × penalty × intercept × method = C × 2 × 2 × 3 × 2 = 24C
- (Same count, just swapped dual for penalty)

**Update SVMParameters dataclass:**

**File:** `classification/optimizer.py` (line ~197)

**Current:**
```python
@dataclass(frozen=True, slots=True)
class SVMParameters:
    """Optimized SVM hyperparameters for LinearSVC with BothEndsStrong features."""

    C: float
    calibration_method: str
    include_max: bool
    dual: bool
    intercept_scaling: float
    cv_score: float
    round_found: int
```

**New:**
```python
@dataclass(frozen=True, slots=True)
class SVMParameters:
    """Optimized SVM hyperparameters for LinearSVC with BothEndsStrong features."""

    C: float
    calibration_method: str
    include_max: bool
    dual: bool  # Keep for backward compatibility, but always False now
    penalty: str  # NEW: 'l1' or 'l2'
    loss: str  # NEW: 'squared_hinge' (required for L1)
    intercept_scaling: float
    cv_score: float
    round_found: int
```

**Extract from grid_search.best_params_:**
```python
best_penalty = grid_search.best_params_['estimator__svc__penalty']
best_loss = grid_search.best_params_['estimator__svc__loss']
```

**Update return statement in _run_optimization_round():**
```python
return SVMParameters(
    C=best_C,
    calibration_method=best_method,
    include_max=best_include_max,
    dual=best_dual,
    penalty=best_penalty,  # NEW
    loss=best_loss,  # NEW
    intercept_scaling=best_intercept_scaling,
    cv_score=best_score,
    round_found=round_idx
)
```

**Update all locations that create LinearSVC:**

1. `optimizer.py` - Final parameter evaluation (line ~787)
2. `trainer.py` - Ensemble training (line ~215)

Add penalty and loss parameters:
```python
LinearSVC(
    C=C,
    dual=dual,
    penalty=penalty,  # NEW
    loss=loss,  # NEW
    intercept_scaling=intercept_scaling,
    class_weight='balanced',
    max_iter=self.max_iter,
    random_state=self.random_state
)
```

---

### Step 5: Implement OutlierClipper Transformer

**File:** `classification/clipping.py` (NEW FILE)

**Purpose:** Winsorize extreme LLR outliers to prevent them from dominating the linear margin.

**Full implementation:**

```python
"""
Outlier clipping transformer for robust feature scaling.

Clips extreme values at high quantiles to prevent rare outliers from
dominating the linear SVM decision boundary and calibration.
"""

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin


class OutlierClipper(BaseEstimator, TransformerMixin):
    """
    Clip extreme outlier values at specified quantiles.

    This transformer learns symmetric clipping thresholds from training data
    and applies them to all data. Helps prevent rare huge LLR scores from
    dominating the linear model.

    Expert guidance: "Rare huge LLRs (e.g., 5' outliers) can dominate a
    linear margin and calibration."

    Attributes:
        quantile: Quantile for symmetric clipping (default: 0.999)
                 Values clipped to [-cap, +cap] where cap = quantile of |X|
        caps_: Learned clipping thresholds per feature (set during fit)

    Example:
        >>> clipper = OutlierClipper(quantile=0.999)
        >>> clipper.fit(X_train)  # Learn caps from training data
        >>> X_train_clipped = clipper.transform(X_train)
        >>> X_test_clipped = clipper.transform(X_test)  # Use same caps
    """

    def __init__(self, quantile: float = 0.999):
        """
        Initialize OutlierClipper.

        Args:
            quantile: Quantile for clipping (0.999 = clip at 99.9th percentile)
                     Higher values = less aggressive clipping
        """
        if not 0 < quantile < 1:
            raise ValueError(f"quantile must be in (0, 1), got {quantile}")
        self.quantile = quantile
        self.caps_ = None

    def fit(self, X, y=None):
        """
        Learn clipping thresholds from training data.

        Computes symmetric caps at the specified quantile of absolute values
        for each feature independently.

        Args:
            X: Training data of shape (n_samples, n_features)
            y: Target values (ignored, for sklearn compatibility)

        Returns:
            self
        """
        X = np.asarray(X)

        # Compute symmetric caps: quantile of |X| per feature
        # Shape: (n_features,)
        self.caps_ = np.quantile(np.abs(X), self.quantile, axis=0)

        return self

    def transform(self, X):
        """
        Apply learned clipping thresholds to data.

        Clips each feature to [-cap, +cap] where cap was learned during fit.

        Args:
            X: Data to clip of shape (n_samples, n_features)

        Returns:
            X_clipped: Clipped data of same shape

        Raises:
            ValueError: If transform called before fit
        """
        if self.caps_ is None:
            raise ValueError("OutlierClipper must be fit before transform")

        X = np.asarray(X)

        # Clip each feature independently
        X_clipped = np.copy(X)
        for i, cap in enumerate(self.caps_):
            X_clipped[:, i] = np.clip(X_clipped[:, i], -cap, cap)

        return X_clipped

    def get_feature_names_out(self, input_features=None):
        """
        Get output feature names (same as input).

        Args:
            input_features: Input feature names (optional)

        Returns:
            Array of output feature names
        """
        if input_features is None:
            # If not provided, use generic names
            n_features = len(self.caps_) if self.caps_ is not None else 0
            return np.array([f'x{i}' for i in range(n_features)])
        return np.asarray(input_features)
```

**Key design decisions:**

1. **Symmetric clipping:** Uses quantile of |X|, applies as [-cap, +cap]
   - Preserves zero-anchoring of LLRs
   - Treats positive and negative outliers equally

2. **Per-feature caps:** Each feature gets its own threshold
   - 5'SS, BPS, 3'SS have different dynamic ranges
   - Independent clipping respects this

3. **Fit on reference only:** Caps learned from training data
   - Prevents data leakage
   - Same caps applied to experimental introns

4. **Default quantile=0.999:** Expert recommendation
   - Clips only extreme outliers (top 0.1%)
   - Conservative, won't affect typical strong signals

---

### Step 6: Integrate OutlierClipper into Pipeline

**Files to update:**
1. `classification/optimizer.py` - Grid search pipeline
2. `classification/trainer.py` - Ensemble training pipeline
3. `classification/nested_cv.py` - Nested CV pipeline (if exists)
4. `classification/split_eval.py` - Split eval pipeline (if exists)

**Pipeline structure change:**

**Current:**
```python
Pipeline([
    ('scale', RobustScaler(with_centering=False, with_scaling=True)),
    ('augment', BothEndsStrongTransformer(include_max=...)),
    ('svc', LinearSVC(...))
])
```

**New:**
```python
Pipeline([
    ('clip', OutlierClipper(quantile=0.999)),  # NEW - clip before scaling
    ('scale', RobustScaler(with_centering=False, with_scaling=True)),
    ('augment', BothEndsStrongTransformer(include_max=...)),
    ('svc', LinearSVC(...))
])
```

**Order is critical:**
1. **Clip** raw LLRs first (prevent outliers)
2. **Scale** clipped LLRs (zero-preserving normalization)
3. **Augment** with composite features (min, neg_absdiff)
4. **SVC** for classification

**Example update for optimizer.py (line ~496):**

```python
from classification.clipping import OutlierClipper  # Add import

# In _run_optimization_round():
base_pipeline = Pipeline([
    ('clip', OutlierClipper(quantile=0.999)),  # NEW
    ('scale', RobustScaler(with_centering=False, with_scaling=True)),
    ('augment', BothEndsStrongTransformer()),
    ('svc', LinearSVC(
        dual=False,
        penalty='l2',  # Will be grid-searched
        loss='squared_hinge',
        class_weight='balanced',
        max_iter=self.max_iter,
        random_state=self.random_state
    ))
])
```

**Repeat for:**
- `trainer.py` line ~194
- Any evaluation scripts that build pipelines

**Note:** Grid search doesn't need to tune quantile - expert says 0.999 is good default.

---

### Step 7: Fix Missing ensemble=True

**File:** `classification/optimizer.py`

**Location:** Line ~787-791 in final parameter evaluation

**Current:**
```python
model = CalibratedClassifierCV(
    base_svm_pipeline,
    method=method,
    cv=5
)
```

**Fixed:**
```python
model = CalibratedClassifierCV(
    base_svm_pipeline,
    method=method,
    cv=5,
    ensemble=True  # NEW - average calibrators across folds
)
```

**Why:** sklearn ≥1.6 defaults to ensemble=True, but we should be explicit. Averaging calibrators across folds improves robustness and reduces variance in probability estimates.

**Verify other locations already have it:**
- Line ~509: Has `ensemble='auto'` ✓
- `trainer.py` line ~215: Has `ensemble='auto'` ✓

---

### Step 8: Document Sigmoid Calibration Preference

**File:** `config/training_default.yaml`

**Location:** Line ~106-114 (calibration method section)

**Current:**
```yaml
  # Probability calibration method
  # - 'sigmoid' = Platt scaling (parametric, fast)
  #               Fits sigmoid to decision function
  #               Works well when calibration data is limited
  # - 'isotonic' = Isotonic regression (non-parametric, flexible)
  #                More flexible, requires more calibration data
  #                Can overfit if calibration set is small
  # Recommendation: Test both, sigmoid often wins for small datasets
  method: ['sigmoid', 'isotonic']
```

**Enhanced:**
```yaml
  # Probability calibration method
  # - 'sigmoid' = Platt scaling (parametric, fast) **RECOMMENDED**
  #               Fits sigmoid to decision function
  #               Works well when calibration data is limited
  #               More conservative, less likely to overfit in top tail
  #               Expert: "Default to sigmoid unless isotonic clearly wins"
  # - 'isotonic' = Isotonic regression (non-parametric, flexible)
  #                More flexible, requires more calibration data
  #                Can overfit the top tail with ~80 positives/fold
  #                May produce more 0.90+ scores (higher FPR)
  # Recommendation: Grid search both, but sigmoid is safer default
  # If isotonic wins by <0.01 neg_log_loss, prefer sigmoid
  method: ['sigmoid', 'isotonic']
```

**Also add note in optimizer comments:**

```python
# Expert guidance: Prefer sigmoid unless isotonic wins by clear margin
# Isotonic can overfit top tail with limited U12s (~80 per fold)
# Sigmoid is more conservative and reduces FPR
```

---

## Summary of Changes

### Files Modified:

1. **`classification/transformers.py`**
   - Add `min_all` feature
   - Add `neg_absdiff_5_bp`, `neg_absdiff_5_3`, `neg_absdiff_bp_3` features
   - Update feature names and docstrings
   - Output: 9D (include_max=False) or 11D (include_max=True)

2. **`classification/clipping.py`** (NEW)
   - Implement `OutlierClipper` transformer
   - Quantile-based symmetric clipping
   - Per-feature caps learned from training data

3. **`classification/classifier.py`**
   - Update default `eff_C_pos_range` to `(3e-4, 1e+2)`
   - Update docstring

4. **`classification/optimizer.py`**
   - Add `penalty` and `loss` fields to `SVMParameters` dataclass
   - Update param_grid: add `penalty: ['l1', 'l2']`, `loss: ['squared_hinge']`
   - Set `dual: [False]` (required for L1)
   - Add `OutlierClipper` to pipeline
   - Fix missing `ensemble=True` at line ~791
   - Update all `best_params_` extraction code
   - Update final parameter evaluation to use new params

5. **`classification/trainer.py`**
   - Add `OutlierClipper` to pipeline
   - Update `LinearSVC` to use `penalty` and `loss` from parameters
   - Update docstrings

6. **`config/training_default.yaml`**
   - Enhance calibration method documentation
   - Add expert guidance on sigmoid preference

### Tests to Update:

1. **`tests/unit/test_classification/test_transformers.py`**
   - Update expected output dimensions: 5→9, 7→11
   - Add tests for `min_all` computation
   - Add tests for `neg_absdiff_*` features
   - Test that features are computed correctly

2. **`tests/unit/test_classification/test_clipping.py`** (NEW)
   - Test `OutlierClipper.fit()` learns correct caps
   - Test `OutlierClipper.transform()` clips correctly
   - Test symmetric clipping
   - Test per-feature independence

3. **`tests/unit/test_classification/test_optimizer.py`**
   - Update for new param_grid structure
   - Test that L1 penalty is grid-searched
   - Test SVMParameters has penalty/loss fields

4. **`tests/unit/test_classification/test_predictor.py`**
   - Update mock models to have penalty/loss attributes
   - Test prediction with 9D/11D features

### Documentation to Update:

1. **`FP_REDUCTION_IMPLEMENTATION_PLAN.md`**
   - Mark Tiers 1 & 2 as complete
   - Document actual implementation choices

2. **`CLAUDE.md`** (project documentation)
   - Update transformer description
   - Document new feature set
   - Update pipeline structure

3. **Docstrings** (throughout)
   - Update feature dimensions in examples
   - Add references to expert guidance
   - Document FP reduction motivation

---

## Testing Strategy

### Unit Tests
```bash
pytest tests/unit/test_classification/test_transformers.py -v
pytest tests/unit/test_classification/test_clipping.py -v
pytest tests/unit/test_classification/test_optimizer.py -v
```

### Integration Tests
```bash
# Run on Chr19 (baseline - should remain stable)
python -m cli.classify --train \
    -g test_data/Homo_sapiens.Chr19.fa.gz \
    -a test_data/Homo_sapiens.Chr19.gff3.gz \
    -n homo_sapiens_chr19 \
    --config config/training_fast_test.yaml

# Check that:
# - Training completes without errors
# - Feature dimensions are correct (9D/11D)
# - U12 detection rate similar to before
# - Log shows penalty and loss parameters
```

### Validation Tests
```bash
# Once C. elegans data is available:
# - Compare FP rate before/after
# - Verify recall doesn't drop
# - Check that L1 is pruning features (inspect coefficients)
```

---

## Migration & Backward Compatibility

### Breaking Changes

**Models trained with old transformer (5D/7D) are NOT compatible with new transformer (9D/11D).**

Users must retrain models after upgrade.

### Config Compatibility

Old configs will work but use new defaults:
- `eff_C_pos_range` defaults to `(3e-4, 1e+2)` (can override)
- Grid searches L1/L2 penalty (can disable via custom param_grid)
- Pipeline includes clipping (transparent to user)

### Data Compatibility

Score output files unchanged:
- min/max columns already exist
- New features are internal to model

---

## Expected Outcomes

### Performance Metrics

**C. elegans (primary target):**
- FP reduction: 35-75% (conservative estimate)
- Recall maintained: ≥95% of current
- Precision increase: 20-50%

**Other species (secondary benefit):**
- FP reduction: 10-30% (species-agnostic improvements)
- No special-casing required
- Improved generalization across phylogenetic distance

### Model Insights

**Feature importance (expected):**
- `min_all` - HIGH positive weight (requires all 3 strong)
- `neg_absdiff_5_bp` - MEDIUM positive weight (penalizes 5'SS/BPS imbalance)
- `neg_absdiff_5_3` - MEDIUM positive weight (penalizes 5'SS/3'SS imbalance)
- `s5`, `sBP`, `s3` - Still important individually
- `max_5_bp`, `max_5_3` - LIKELY ZERO (L1 prunes OR-like features)

**L1 vs L2 outcome:**
- L1 expected to win for most species (pruning helps)
- Coefficients will be sparser with L1
- Easier to interpret which signals matter

---

## Next Steps After Implementation

1. **Run comprehensive tests** on Chr19 human data
2. **Validate** FP reduction on C. elegans data
3. **Inspect** learned feature weights (which got pruned?)
4. **Document** actual FP reduction achieved
5. **Consider Tier 3** if needed (1-SE rule for model selection)

---

## Notes & Caveats

1. **OutlierClipper quantile:** 0.999 is expert recommendation
   - Could make tunable later if needed
   - Currently fixed for simplicity

2. **L1 penalty constraints:** Requires dual=False, loss=squared_hinge
   - We've hardcoded dual=False (fine since n_samples >> n_features)
   - Could add back dual search if needed for other datasets

3. **Feature count:** Now 9D-11D vs original 5D-7D
   - More features = more risk of overfitting
   - But L1 penalty mitigates this by pruning
   - Cross-validation ensures we don't overfit

4. **C bounds:** (3e-4, 1e+2) vs expert's (3e-4, 1e+1)
   - We extended upper bound slightly for flexibility
   - User requested this change
   - Can tighten further if FPR still too high

5. **Backward compatibility:** Old models won't load
   - Clear error message needed
   - Documentation should warn users to retrain
   - Consider adding version check in model loading

---

## References

- Expert guidance on FP reduction (2024)
- sklearn LinearSVC: https://scikit-learn.org/stable/modules/generated/sklearn.svm.LinearSVC.html
- sklearn CalibratedClassifierCV: https://scikit-learn.org/stable/modules/calibration.html
- Winsorization: https://en.wikipedia.org/wiki/Winsorizing
