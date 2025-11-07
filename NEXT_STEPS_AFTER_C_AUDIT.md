# Next Steps After C Value Audit

**Date:** 2025-11-07
**Branch:** claude/investigate-ref-performance-011CUt9YD4yNb4AzXF2TKAU9
**Context:** Completed audit of C parameter optimization

---

## Audit Results Summary

### Question
"Is the optimizer exploring the full range of C values, including C≈1000 where the original finds its optimum?"

### Answer
**✅ YES - The optimizer IS exploring the full C range correctly.**

### Evidence
- Initial grid: `np.logspace(-6, 6, num=13)` generates [1e-6, ..., 1e3, ..., 1e6]
- C=1000 is at position 10/13 in the grid
- Debug output confirms all 13 values are tested in Round 1
- Optimizer logic is working as designed

### Key Discovery
LinearSVC + CalibratedClassifierCV shows **CV score saturation**:
```
C=1e-4:  CV=0.5719 (selected as best)
C=1e-3:  CV=0.5687
C=1e-2:  CV=0.5687
C=1e-1:  CV=0.5687
C=1.0:   CV=0.5687
C=10:    CV=0.5687
C=100:   CV=0.5687
C=1000:  CV=0.5687  ← No benefit despite original finding optimum here
C=1e4:   CV=0.5687
C=1e5:   CV=0.5687
C=1e6:   CV=0.5687
```

All C values from 1e-3 to 1e+6 give **IDENTICAL** CV scores.

---

## What We've Ruled Out

The performance problem (F1=0.31 vs F1=1.0, 0 U12s found vs 32 expected) is **NOT** due to:

- ❌ Incorrect C search range
- ❌ Broken optimizer logic
- ❌ Missing high C values in the grid
- ❌ Off-by-one errors in grid generation
- ❌ Wrong scoring metric (balanced_accuracy is correct)
- ❌ Wrong cross-validation setup (5-fold matches original)

---

## What We Know

### The Core Problem
LinearSVC + external calibration fundamentally differs from SVC + internal calibration:

**Original (Works):**
```python
SVC(kernel='linear', probability=True, class_weight='balanced')
```
- Uses libsvm
- Probability calibration built INTO the optimization
- Sees meaningful CV score variation across C values
- Finds optimal C≈1000
- Achieves F1=1.0, finds 32 U12s

**Refactored (Doesn't Work):**
```python
base = LinearSVC(C=?, class_weight='balanced', dual=False)
model = CalibratedClassifierCV(base, method='sigmoid', cv=5)
```
- Uses liblinear (100x faster)
- Probability calibration applied AFTER optimization
- Shows CV score saturation at C > 1e-3
- Selects C≈9.3e-05 (too regularized)
- Achieves F1=0.31, finds 0 U12s

### User's Hypothesis
The CV score saturation **may be normal** for well-separated classes. If true, the problem isn't saturation per se, but something about how the two approaches differ in their final classification behavior despite similar CV scores.

---

## Open Questions

1. **Why does LinearSVC+calibration fail to classify ANY U12s?**
   - CV scores are similar (~0.57)
   - Training appears successful
   - But final predictions: 0 U12s vs expected 32

2. **How do the decision boundaries differ?**
   - Both use linear kernels
   - Both use balanced class weights
   - Both achieve similar CV balanced accuracy
   - Yet one classifies correctly, the other doesn't

3. **Does probability calibration timing matter this much?**
   - SVC: Calibration during training (internal Platt scaling)
   - LinearSVC: Calibration after training (external wrapper)
   - Could this explain the complete classification failure?

4. **Is the threshold application different?**
   - Original: `predict_proba()` from SVC with internal calibration
   - Refactored: `predict_proba()` from CalibratedClassifierCV wrapper
   - Could the probability scales differ despite both being [0, 1]?

---

## Recommended Next Steps

### Option 1: Investigate Classification Failure Directly
Instead of focusing on CV scores, examine what's happening at prediction time:

**Test:**
1. Train both models with the same C value (e.g., C=1000)
2. Compare their predictions on the same experimental introns
3. Examine probability distributions: `predict_proba()`
4. Compare decision function values (where available)
5. Check if probability scales differ systematically

**Goal:** Understand why LinearSVC+calibration predicts 0 U12s while SVC predicts 32.

### Option 2: Try SVC with External Calibration
Test whether the problem is LinearSVC specifically or external calibration in general:

**Test:**
```python
# Keep SVC but use external calibration like refactored version
base = SVC(kernel='linear', probability=False, class_weight='balanced', C=1000)
model = CalibratedClassifierCV(base, method='sigmoid', cv=5)
```

**Goal:** Isolate whether the issue is:
- LinearSVC vs SVC (algorithm difference)
- External vs internal calibration (architectural difference)

### Option 3: Investigate Probability Calibration Quality
Check if CalibratedClassifierCV is miscalibrating probabilities:

**Test:**
1. Train LinearSVC with fixed C=1000
2. Get raw decision function: `base_svm.decision_function(X)`
3. Get calibrated probabilities: `cal_svm.predict_proba(X)`
4. Plot calibration curve: reliability diagram
5. Compare to SVC internal calibration

**Goal:** Determine if external calibration is producing incorrect probability estimates.

### Option 4: Just Revert to SVC
Accept the 100x slowdown and use the original approach:

**Pros:**
- Proven to work
- No further investigation needed
- Can complete refactoring

**Cons:**
- 2-3 minutes vs 40 seconds
- Slower for users
- Doesn't understand root cause

---

## Current Investigation Documents

- `C_VALUE_AUDIT_RESULTS.md` - Complete audit findings with detailed CV scores
- `LINEAR_SVC_INVESTIGATION.md` - Earlier investigation of LinearSVC convergence issues
- `PERFORMANCE_COMPARISON_RESULTS.md` - Speed comparison showing 4.5x overall speedup
- `INVESTIGATION_COMPLETE.md` - Summary of previous investigation phase

---

## Files with Debug Code (Can Remove Later)

- `intronIC_refactored/classification/optimizer.py:302-317` - Round 1 CV score printing
- `intronIC/intronIC.py:5477-5494` - Original SVC CV score printing (never completed)
- `test_c_values.py` - Incomplete direct C value test script

---

## Recommendation

**Start with Option 1** - investigate why LinearSVC+calibration predicts 0 U12s despite reasonable CV scores. This will give us the most direct insight into the classification failure, which is the actual problem (not the CV score saturation, which may be normal).

If Option 1 doesn't reveal a fixable issue, proceed to **Option 2** to determine if the problem is LinearSVC itself or the external calibration architecture.

If both fail to provide a solution, **Option 4** (revert to SVC) may be the pragmatic choice.

---

**Status:** C value audit complete, ready to investigate classification failure mechanism.
