# Performance Testing Strategy

## Problem Summary

**Current Status:**
- Original: 49.7 seconds for small dataset (97 U12, 990 U2, 10:1 imbalance)
- Refactored: ~30-40 minutes (40x slower)
- Both now have `class_weight='balanced'`, `train_test_split(0.80)`, `cache_size=1000`

**The Mystery:**
- Original has `probability=True` (should be 11x slower according to isolated tests)
- Yet original is 40x FASTER than refactored!
- Need to understand what original is actually doing

---

## Test Scripts Created

### 1. `test_training_minimal.py` - Core Comparison
**Purpose:** Minimal test with synthetic data matching your reference set size

**Tests 4 approaches:**
1. Original style: `SVC(probability=True)` + `class_weight` as attribute
2. Refactored style: `SVC(class_weight='balanced')` in constructor
3. LinearSVC with calibration (for probabilities)
4. LinearSVC alone (fastest)

**Expected runtime:** ~2-5 minutes
**Key insights:**
- Does `clone()` preserve `class_weight` when set as attribute?
- Is `probability=True` actually the bottleneck?
- How much faster is LinearSVC?

### 2. `test_high_c_behavior.py` - High C Value Analysis
**Purpose:** Test where performance diverges (your data shows issues at C≥100,000)

**Tests:** C values from 1e2 to 1e6
**Expected runtime:** 10-20 minutes
**Key insights:**
- At what C value does performance diverge?
- Does time ratio increase with C?
- How does LinearSVC perform at high C?

### 3. `test_gridsearch_behavior.py` - GridSearchCV Inspection
**Purpose:** Examine GridSearchCV's actual behavior with verbose output

**Expected runtime:** 5-10 minutes
**Key insights:**
- Are both versions testing the same C values?
- What are the actual fit times per parameter?
- Any convergence warnings or NaN scores?
- Does clone() preserve all parameters?

### 4. `test_linearsvc_with_calibration.py` - LinearSVC Integration
**Purpose:** Show how to integrate LinearSVC while maintaining probability output

**Tests 3 approaches:**
1. Optimize LinearSVC, then calibrate (RECOMMENDED)
2. Optimize pre-calibrated model
3. Original SVC with probability=True

**Expected runtime:** 3-5 minutes
**Key insight:** Best way to integrate LinearSVC into refactored code

---

## Recommended Testing Sequence

### Phase 1: Quick Diagnosis (Run First)
```bash
python test_training_minimal.py > results_minimal.log 2>&1
```

**What to look for:**
- Time comparison: Original vs Refactored vs LinearSVC
- Clone behavior: Does `class_weight` get preserved?
- If LinearSVC is 10x+ faster → That's our solution

**Decision point:**
- If LinearSVC is dramatically faster → Skip to Phase 3 (implement it)
- If Original and Refactored are similar → Investigate other factors
- If Original is faster → Proceed to Phase 2

### Phase 2: Root Cause Analysis (If Phase 1 shows Original faster)
```bash
python test_high_c_behavior.py > results_high_c.log 2>&1
python test_gridsearch_behavior.py > results_gridsearch.log 2>&1
```

**What to look for:**
- Performance divergence at specific C values
- Clone() behavior differences
- Convergence warnings
- Actual CV timing per parameter

### Phase 3: Solution Integration (LinearSVC)
```bash
python test_linearsvc_with_calibration.py > results_linearsvc.log 2>&1
```

**What to implement:**
- Use "Approach 1" from test: Optimize then calibrate
- Modify `intronIC_refactored/classification/optimizer.py`
- Keep probability output for compatibility

---

## Expected Outcomes & Solutions

### Scenario A: LinearSVC is 10-50x faster (MOST LIKELY)
**Root cause:** `SVC(kernel='linear')` uses libsvm (optimized for non-linear kernels)
**Solution:** Switch to LinearSVC (uses liblinear, optimized for linear problems)

**Implementation:**
```python
# In optimizer.py _grid_search_round():
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV

# Option 1: Optimize LinearSVC, then calibrate (RECOMMENDED)
base_svm = LinearSVC(
    class_weight='balanced',
    dual='auto',
    max_iter=2000,
    random_state=self.random_state + round_idx
)

grid_search = GridSearchCV(base_svm, param_grid={'C': C_grid}, ...)
grid_search.fit(X_train, y_train)
best_C = grid_search.best_params_['C']

# Create final model with calibration for probabilities
final_model = LinearSVC(C=best_C, class_weight='balanced', dual='auto', ...)
calibrated_model = CalibratedClassifierCV(final_model, method='sigmoid', cv=3)
```

**Expected speedup:** 10-50x vs current refactored version

---

### Scenario B: clone() doesn't preserve class_weight attribute
**Root cause:** sklearn.base.clone() may not copy attributes set after initialization
**Solution:** Always pass parameters to constructor

**Evidence to look for:**
```python
# In test output, check:
Clone from original: class_weight = 'balanced' or 'NOT SET'?
Clone from refactored: class_weight = 'balanced' (should work)
```

**If original's clone shows 'NOT SET':**
- Original might be training without class_weight in some CV folds!
- This would make it faster but scientifically wrong
- Refactored version is correct, just need to optimize

---

### Scenario C: probability=True paradox explained
**Possible explanation:** probability=True triggers different code path in libsvm
**Look for:** Different convergence behavior, fewer iterations

**If true:** This would mean original is "cheating" via early convergence
**Solution:** Still use LinearSVC for scientifically sound performance

---

## Modern Best Practices (From Field)

### For Imbalanced Linear SVM:
1. ✅ Use `LinearSVC` (not `SVC(kernel='linear')`) - **10-100x faster**
2. ✅ Set `class_weight='balanced'` - **Handles imbalance properly**
3. ✅ Use `dual='auto'` - **Automatically chooses primal/dual**
4. ✅ For probabilities: Wrap in `CalibratedClassifierCV` - **Fast calibration**
5. ✅ Set reasonable `max_iter` - **Prevents infinite loops**

### Your Current Code:
- ✅ Has class_weight='balanced'
- ✅ Has proper cross-validation
- ❌ Uses SVC instead of LinearSVC (KEY ISSUE)
- ✅ Has probability output (need to maintain with calibration)

---

## Implementation Roadmap

### Step 1: Run Tests (User)
```bash
cd /home/user/intronIC
python test_training_minimal.py > results1.log 2>&1
# Review results1.log, check timings
```

### Step 2: If LinearSVC wins, integrate it (Claude)
- Modify `intronIC_refactored/classification/optimizer.py`
- Replace SVC with LinearSVC
- Add CalibratedClassifierCV wrapper after optimization
- Test on small dataset
- Verify probability output matches

### Step 3: Benchmark on actual data
```bash
# Small reference dataset
time python -m cli.main --genome ... --reference_u12s u12_reference_small.introns.iic.gz --reference_u2s u2_reference_small.introns.iic.gz -p 4

# Should complete in ~50 seconds if LinearSVC works!
```

### Step 4: Full validation
- Run on full reference set
- Compare outputs with original
- Verify scientific correctness
- Commit optimized version

---

## Why LinearSVC is the Right Solution

### Technical Reasons:
1. **Liblinear vs Libsvm:** Liblinear is specifically optimized for linear classification
2. **Coordinate descent:** Faster convergence for linear problems
3. **Primal vs Dual:** Automatically chooses based on n_samples vs n_features
4. **No kernel cache:** Doesn't waste memory/time on kernel computations
5. **Better scaling:** O(n_features) instead of O(n_samples²)

### Your Specific Case:
- n_samples = 1087 (large enough for primal)
- n_features = 3 (small)
- Class imbalance = 10:1 (LinearSVC handles this well)
- High C values tested (liblinear converges faster)
- Linear kernel (liblinear's specialty)

### Field Validation:
- Standard approach for large-scale linear SVM
- Used in major bioinformatics tools
- Proven for imbalanced classification
- Scikit-learn documentation recommends it

---

## Files Created

1. `test_training_minimal.py` - Quick comparison test
2. `test_high_c_behavior.py` - High C value analysis
3. `test_gridsearch_behavior.py` - GridSearchCV inspection
4. `test_linearsvc_with_calibration.py` - Integration guide
5. `PERFORMANCE_TESTING_STRATEGY.md` - This document

---

## Next Actions

**For User:**
1. Run `test_training_minimal.py` first
2. Share the output
3. We'll analyze and decide on implementation strategy

**For Claude (after results):**
1. Analyze test results
2. Implement LinearSVC solution if validated
3. Test on actual data
4. Document changes
5. Commit to branch

---

## Expected Timeline

- Tests: 10-30 minutes total
- Analysis: 10 minutes
- Implementation: 30-60 minutes
- Validation: 30 minutes
- **Total:** 2-3 hours to 40x speedup

## Success Criteria

✅ Refactored version completes in <2 minutes (vs original's 49.7s)
✅ Outputs identical probabilities (within tolerance)
✅ Scientific correctness maintained
✅ Code is cleaner and more maintainable
✅ Uses modern best practices

---

**Status:** Test scripts ready, awaiting Phase 1 results
