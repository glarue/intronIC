# False Positive Reduction Implementation Plan
## Expert Recommendations for C. elegans FP Reduction

**Goal:** Reduce false positive rate in C. elegans (and other species) by enhancing the feature set and regularization strategy while maintaining species-agnostic approach.

**Source:** Expert ML guidance on handling "one-end-strong" false positives

---

## Priority Ranking

### Tier 1: High Impact, Easy Implementation ⭐⭐⭐
1. **Add min_all feature** (3-way AND)
2. **Add negative absdiff penalty features**
3. **Tighten C bounds default** (reduce FPR)
4. **Add L1 penalty option** to param_grid

### Tier 2: Medium Impact, Medium Complexity ⭐⭐
5. **Winsorization/clipping** of extreme outliers
6. **Prefer sigmoid calibration** by default
7. **Ensure ensemble=True** in all CalibratedClassifierCV calls

### Tier 3: Nice to Have, Higher Complexity ⭐
8. **1-SE rule** for model selection (prefer simpler models)

---

## Detailed Implementation Plan

### 1. Add min_all Feature (3-way AND) ⭐⭐⭐

**Why:** FPs cluster where one end is strong and third signal is modest. Requiring all three to be high cuts these.

**How:**
- File: `classification/transformers.py`
- In `BothEndsStrongTransformer.transform()`:
  ```python
  # After computing min_5_bp and min_5_3:
  min_all = np.minimum(np.minimum(s5, sBP), s3)
  ```
- Add to feature list
- Update `get_feature_names_out()` to include 'min_all'
- Update docstrings and examples

**Expected output dimension:**
- Current: 5D (include_max=False) or 7D (include_max=True)
- New: 6D (include_max=False) or 8D (include_max=True)

**Impact:** HIGH - Expert says this is critical for cutting one-end-strong FPs

---

### 2. Add Negative Absdiff Penalty Features ⭐⭐⭐

**Why:** Even with min_*, the model may not penalize one-end-strong enough. Adding negated absdiff lets model put positive weight on "consistency".

**How:**
- File: `classification/transformers.py`
- We already compute `absdiff_5_bp` and `absdiff_5_3` for min/max calculation
- Add as negative features:
  ```python
  neg_absdiff_5_bp = -absdiff_5_bp
  neg_absdiff_5_3 = -absdiff_5_3
  # Optional:
  absdiff_bp_3 = np.abs(sBP - s3)
  neg_absdiff_bp_3 = -absdiff_bp_3
  ```
- Add to feature list
- Update `get_feature_names_out()`: 'neg_absdiff_5_bp', 'neg_absdiff_5_3', ['neg_absdiff_bp_3']

**Parameters:**
- Expert mentions optional weights γ ∈ {1,2,4,8} but these can be handled by SVM coefficients
- Could add as param_grid option if needed: `estimator__augment__absdiff_weight: [1, 2, 4, 8]`
- Start simple: just add features with weight=1

**Expected output dimension:**
- With min_all + neg_absdiff_5_bp + neg_absdiff_5_3: 8D (include_max=False)
- With all + neg_absdiff_bp_3: 9D (include_max=False)
- With max features: 10D or 11D

**Key insight:** "you can drop max_* entirely" - so default include_max=False becomes even more justified

**Impact:** HIGH - Directly penalizes imbalance

---

### 3. Tighten C Bounds Default ⭐⭐⭐

**Why:** Larger C tends to raise FPR. Shrinking effective positive-class band reduces this.

**How:**
- File: `classification/classifier.py`
- Change default `eff_C_pos_range` in `__init__()`:
  ```python
  # OLD:
  eff_C_pos_range: tuple = (1e-3, 1e3)

  # NEW (expert recommendation):
  eff_C_pos_range: tuple = (3e-4, 1e+1)
  ```
- Update docstring to explain the change
- Note: C is back-solved via `C = C_eff_pos / w_pos`

**Impact:** MEDIUM-HIGH - Simple change, should reduce FPs globally

---

### 4. Add L1 Penalty Option ⭐⭐⭐

**Why:** With LLRs + min_* + neg_absdiff_* we have correlated features. L1 often zeros the ones that re-introduce one-end behavior.

**How:**
- File: `classification/optimizer.py`
- Add to param_grid (line ~527):
  ```python
  param_grid = {
      'estimator__svc__C': C_grid,
      'estimator__augment__include_max': [False, True],
      'estimator__svc__dual': [False, True],
      'estimator__svc__penalty': ['l1', 'l2'],  # NEW
      'estimator__svc__intercept_scaling': [10.0, 100.0, 1000.0],
      'method': ['sigmoid', 'isotonic']
  }
  ```

**CONSTRAINT HANDLING:**
- L1 requires `dual=False` and `loss='squared_hinge'`
- Options:
  1. Add constraint to grid: use custom ParameterGrid that filters invalid combos
  2. Let GridSearchCV handle it (it will error and use error_score=np.nan)
  3. Use sklearn's `ParameterGrid` with conditional logic

**Recommendation:** Use option 2 initially (rely on error_score=np.nan), then add proper constraint handling if it causes issues.

**Impact:** MEDIUM - May prune redundant features and reduce FPs

---

### 5. Winsorization/Clipping of Extreme Outliers ⭐⭐

**Why:** Rare huge LLRs (e.g., 5' outliers) can dominate a linear margin and calibration.

**How:**
- **Option A:** Add to scaler
  - File: `classification/scaler.py` (ZeroPreservingScaler)
  - Compute quantiles from reference data: Q_0.999 for each feature
  - Clip before scaling: `s5 = np.clip(s5, -cap5, cap5)`
  - Store caps as attributes for applying to experimental data

- **Option B:** Create separate ClippingTransformer
  - New file: `classification/clipping.py`
  - Fit on reference data to learn caps
  - Transform experimental data with same caps
  - Add to pipeline before scaler

**Recommendation:** Option B (separate transformer) is cleaner and more modular.

**Implementation:**
```python
class OutlierClipper(BaseEstimator, TransformerMixin):
    def __init__(self, quantile=0.999):
        self.quantile = quantile
        self.caps_ = None

    def fit(self, X, y=None):
        # Compute symmetric caps at quantile
        self.caps_ = np.quantile(np.abs(X), self.quantile, axis=0)
        return self

    def transform(self, X):
        X_clipped = np.copy(X)
        for i, cap in enumerate(self.caps_):
            X_clipped[:, i] = np.clip(X_clipped[:, i], -cap, cap)
        return X_clipped
```

**Impact:** MEDIUM - Prevents extreme outliers from dominating

---

### 6. Prefer Sigmoid Calibration by Default ⭐⭐

**Why:** Isotonic can overfit the top tail with ~80 positives/fold and produce more 0.90s. Platt is often more conservative.

**How:**
- File: `config/training_default.yaml`
- Reorder calibration methods to prefer sigmoid:
  ```yaml
  # OLD:
  method: ['sigmoid', 'isotonic']

  # NEW:
  method: ['sigmoid', 'isotonic']  # Order doesn't matter for GridSearchCV
  ```

- Actually, grid search picks best by CV score, not order. To "prefer" sigmoid, we'd need to:
  - **Option A:** Only search sigmoid initially: `method: ['sigmoid']`
  - **Option B:** Add 1-SE rule that prefers sigmoid when tied
  - **Option C:** Add post-selection logic to prefer sigmoid when scores are close

**Recommendation:** Document in config comments that sigmoid is recommended default, but keep both in grid. Users can override to sigmoid-only if desired.

**Impact:** LOW-MEDIUM - Mainly affects calibration quality, not core FP rate

---

### 7. Ensure ensemble=True in CalibratedClassifierCV ⭐⭐

**Why:** sklearn ≥1.6 defaults to ensemble=True, but we should be explicit. Averaging calibrators across folds improves robustness.

**How:**
- File: `classification/optimizer.py`
  - Line 509-514: Already has `ensemble='auto'` ✓
  - Line 787-791: Missing `ensemble=` parameter - ADD IT

- File: `classification/trainer.py`
  - Line 215-220: Already has `ensemble='auto'` ✓

**Fix needed:**
```python
# In optimizer.py, line 787-791:
model = CalibratedClassifierCV(
    base_svm_pipeline,
    method=method,
    cv=5,
    ensemble=True  # ADD THIS
)
```

**Impact:** LOW - We already use ensemble='auto' in most places

---

### 8. 1-SE Rule for Model Selection ⭐

**Why:** Prefers simpler models (smaller C, fewer features) when within 1 std error of best CV score. Reduces overfitting and FPs.

**How:**
- File: `classification/optimizer.py`
- Modify `_run_optimization_round()` to:
  1. Track std error of CV scores (not just mean)
  2. After finding best score, identify all configs within 1 SE
  3. Among those, select simplest:
     - Smallest C
     - Fewest features (include_max=False preferred)
     - L2 over L1 (or vice versa)

**Implementation sketch:**
```python
# After grid_search.fit():
cv_results = grid_search.cv_results_
best_idx = grid_search.best_index_
best_score = cv_results['mean_test_score'][best_idx]
std_error = cv_results['std_test_score'][best_idx]

# Find all within 1 SE
threshold = best_score + std_error  # neg_log_loss, so higher=worse, we want lower
candidates = np.where(cv_results['mean_test_score'] >= threshold)[0]

# Among candidates, pick simplest (smallest C, include_max=False, etc.)
# This requires sorting by C, then by include_max, etc.
```

**Complexity:** MEDIUM-HIGH - Requires careful implementation to avoid bugs

**Impact:** MEDIUM - Mainly helps with generalization, moderate FP reduction

---

## Implementation Order

### Phase 1: Core Feature Improvements (1-2 days)
1. ✅ Add `min_all` feature to BothEndsStrongTransformer
2. ✅ Add `neg_absdiff_*` features to BothEndsStrongTransformer
3. ✅ Update tests for new feature dimensions
4. ✅ Update documentation

### Phase 2: Regularization Tuning (1 day)
5. ✅ Tighten default C bounds
6. ✅ Add L1 penalty to param_grid
7. ✅ Fix missing ensemble=True in optimizer.py

### Phase 3: Advanced Features (2-3 days)
8. ⏸️ Implement OutlierClipper transformer
9. ⏸️ Add to pipeline and test
10. ⏸️ Update config documentation

### Phase 4: Selection Improvements (1-2 days, optional)
11. ⏸️ Implement 1-SE rule in optimizer
12. ⏸️ Add preference for sigmoid calibration

---

## Expected Impact

### Feature Set Evolution

**Before (current):**
- 3D input: [s5, sBP, s3]
- 5D output (include_max=False): [s5, sBP, s3, min_5_bp, min_5_3]
- 7D output (include_max=True): [s5, sBP, s3, min_5_bp, min_5_3, max_5_bp, max_5_3]

**After (Phase 1):**
- 3D input: [s5, sBP, s3]
- 9D output (include_max=False): [s5, sBP, s3, min_5_bp, min_5_3, min_all, neg_absdiff_5_bp, neg_absdiff_5_3, neg_absdiff_bp_3]
- 11D output (include_max=True): + [max_5_bp, max_5_3]

**After (Phase 3 - with clipping):**
- Same dimensions, but features are winsorized at Q_0.999 before scaling

### Performance Expectations

**C. elegans FP Reduction:**
- Phase 1 (features): 20-40% FP reduction
- Phase 2 (regularization): Additional 10-20% FP reduction
- Phase 3 (clipping): Additional 5-15% FP reduction
- **Total expected:** 35-75% FP reduction while maintaining recall

**Other species:**
- Improvements should be species-agnostic
- No special-casing required
- May reduce FPs globally while maintaining U12 detection

---

## Testing Strategy

### Unit Tests
- Test new transformer feature dimensions
- Test clipping transformer separately
- Test parameter grid combinations

### Integration Tests
- Run on Chr19 human data (baseline)
- Run on C. elegans data (FP reduction target)
- Compare FP rates before/after

### Regression Tests
- Ensure U12 recall doesn't drop
- Ensure human genome results remain stable
- Check that old models with fewer features still load (backward compat)

---

## Migration Notes

### Backward Compatibility

**Old models (5D or 7D features):**
- Will NOT be compatible with new transformer (9D or 11D features)
- Users must retrain models after upgrade
- Add clear warning in release notes

**Config files:**
- Old configs will work (parameters are additive)
- New configs will use tighter C bounds by default

**Data files:**
- Score output files will have same format
- Min/max columns already present

---

## Questions for Expert (Optional Follow-up)

1. **Negative absdiff weights:** Should we grid-search γ ∈ {1,2,4,8} or let SVM coefficients handle it?
2. **min_all vs min_5_bp + min_5_3:** Is min_all redundant, or does it add signal?
3. **Winsorization quantile:** Is Q_0.999 the right choice, or should we tune this?
4. **L1 vs L2:** Which is preferred when CV scores are tied?

---

## References

- Expert guidance (2024): "Feature tweaks for FP reduction"
- sklearn.svm.LinearSVC docs: https://scikit-learn.org/stable/modules/generated/sklearn.svm.LinearSVC.html
- sklearn.calibration.CalibratedClassifierCV docs: https://scikit-learn.org/stable/modules/generated/sklearn.calibration.CalibratedClassifierCV.html
- 1-SE rule: Breiman et al. (1984), "Classification and Regression Trees"
