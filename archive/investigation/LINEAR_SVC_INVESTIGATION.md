# LinearSVC Implementation Investigation

**Date:** 2025-11-07
**Branch:** `claude/investigate-ref-performance-011CUt9YD4yNb4AzXF2TKAU9`
**Status:** ⚠️ ONGOING - Poor Classification Performance

---

## Problem Summary

After implementing sklearn best practices for LinearSVC with imbalanced data, the refactored version still has poor classification performance:

**Current Results:**
- F1 Score: 0.28 (expected: 1.0)
- PR-AUC: 0.38 (expected: 1.0)
- U12s found: 0 (expected: 32)
- Convergence warnings throughout training

**Original Results (SVC):**
- F1 Score: 1.0
- PR-AUC: 1.0
- U12s found: 32
- Runtime: 2m 59s

---

## Changes Implemented

### 1. Initial LinearSVC Switch (commit 03b0a76)
```python
# From:
SVC(kernel='linear', probability=True, class_weight='balanced')

# To:
LinearSVC(class_weight='balanced', dual='auto')
+ CalibratedClassifierCV(method='sigmoid', cv=3)
```

**Result:** AttributeError - decision_function() not supported by Calibrated

ClassifierCV

### 2. Fixed decision_function() Issue (commit b238a1c)
```python
# Replaced decision_function() with log-odds calculation
relative_scores = np.log(clipped_probas / (1 - clipped_probas))
```

**Result:** Run completed but F1=0.33, 1 U12 found (still wrong)

### 3. Applied sklearn Best Practices (commit c778a1e)
```python
LinearSVC(
    C=parameters.C,
    class_weight='balanced',
    loss='squared_hinge',  # NEW
    penalty='l2',          # NEW
    dual=True,             # Changed from 'auto'
    max_iter=2000,
    random_state=seed
)
# Changed calibration cv=3 to cv=5
```

**Result:** F1=0.28, 0 U12s found, extensive convergence warnings

### 4. Increased max_iter (commit b6a93c2)
```python
LinearSVC(
    ...
    max_iter=10000,  # Increased from 2000
    tol=1e-4,        # Added explicit tolerance
    ...
)
```

**Result:** Testing now...

---

## Root Cause Analysis

### Issue 1: Convergence Failures
The model isn't converging with max_iter=2000, suggesting:
- The optimization problem is harder than expected
- The C value selected might not be optimal for LinearSVC
- Imbalanced data (387 U12 vs 20,690 U2) is challenging

### Issue 2: Optimizer/Trainer Mismatch
**Critical observation:** The optimizer optimizes **bare LinearSVC** parameters, but the trainer uses **CalibratedClassifierCV-wrapped LinearSVC**.

**Original approach:**
```python
# Optimization uses SVC with probability=True from the start
SVC(kernel='linear', probability=True, ...)
# The probability estimation is built into the model being optimized
```

**Our approach:**
```python
# Optimization uses bare LinearSVC
LinearSVC(...)
# Then trainer wraps it in calibration
CalibratedClassifierCV(LinearSVC(...), cv=5)
```

**Problem:** The optimal C value found during optimization might not be optimal for the calibrated version!

### Issue 3: dual=True May Not Be Optimal
From sklearn docs:
> "Prefer dual=False when n_samples > n_features"

Our case: ~17k training samples (80% of 21k) vs 3 features
- n_samples (17k) >> n_features (3)
- This suggests **dual=False** might actually be correct!

But sklearn LinearSVC best practices for imbalanced data specifically recommend dual=True when n_features << n_samples. This is contradictory!

---

## Hypothesis for Poor Performance

The poor classification performance might be due to:

1. **Optimizer optimizes wrong model:** We're optimizing bare LinearSVC but using calibrated version
2. **Wrong dual setting:** dual=True might not be appropriate despite best practices
3. **Convergence issues:** Model not converging even with max_iter=10000
4. **Calibration problems:** cv=5 calibration on top of imbalanced data might be failing

---

## Potential Solutions to Try

### Option A: Optimize Calibrated Model Directly
```python
# In optimizer.py, wrap LinearSVC in calibration during optimization
base_svm = LinearSVC(...)
model = CalibratedClassifierCV(base_svm, method='sigmoid', cv=5)
grid_search = GridSearchCV(model, param_grid={...})
```

**Pros:** Finds C optimal for calibrated model
**Cons:** Much slower optimization (5x more fits per C value)

### Option B: Try dual=False
```python
LinearSVC(
    dual=False,  # Use primal formulation for n_samples >> n_features
    loss='squared_hinge',
    penalty='l2',
    ...
)
```

**Pros:** May converge faster, may be more appropriate for our data
**Cons:** Contradicts sklearn best practices for imbalanced data

### Option C: Go Back to SVC
```python
# Just use the original approach that works
SVC(
    kernel='linear',
    probability=True,
    class_weight='balanced',
    cache_size=1000
)
```

**Pros:** Known to work, simpler
**Cons:** Slower (but maybe that's okay if it works)

### Option D: Use LinearSVC Without Calibration
```python
# During optimization: optimize bare LinearSVC
# During training: also use bare LinearSVC
# Skip calibration entirely
```

**Pros:** Consistent optimization/training, faster
**Cons:** No probability outputs (only decision_function)

---

## Test Results Log

| Test | Config | F1 | PR-AUC | U12s | Convergence | Runtime |
|------|--------|-----|--------|------|-------------|---------|
| Original | SVC, prob=True | 1.0 | 1.0 | 32 | ✅ | 2m 59s |
| Refactored v1 | LinearSVC, dual='auto', cv=3 | - | - | CRASH | - | - |
| Refactored v2 | LinearSVC, dual='auto', cv=3, log-odds | 0.33 | 0.42 | 1 | ⚠️ | ~40s |
| Refactored v3 | LinearSVC, dual=True, loss/penalty, cv=5 | 0.28 | 0.38 | 0 | ❌ | ~40s |
| Refactored v4 | LinearSVC, dual=True, max_iter=10k | ? | ? | ? | ? | ? |

---

## Next Steps

1. **Test with max_iter=10000** - See if convergence improves
2. **If still failing:** Try dual=False
3. **If still failing:** Optimize calibrated model directly (Option A)
4. **If still failing:** Consider reverting to SVC (Option C)

---

## Lessons Learned

1. **LinearSVC is not a drop-in replacement for SVC(kernel='linear')** - Different behavior with calibration
2. **Sklearn best practices may conflict** - dual=True vs dual=False recommendations contradict
3. **Optimizer/trainer consistency matters** - Must optimize the same model you'll use
4. **Imbalanced data is hard** - 387 vs 20,690 samples (~1.8% positive) is extremely imbalanced

---

## References

- sklearn LinearSVC docs: https://scikit-learn.org/stable/modules/generated/sklearn.svm.LinearSVC.html
- sklearn SVC docs: https://scikit-learn.org/stable/modules/generated/sklearn.svm.SVC.html
- Calibration docs: https://scikit-learn.org/stable/modules/calibration.html
- Imbalanced learning: https://scikit-learn.org/stable/modules/generated/sklearn.utils.class_weight.compute_class_weight.html
