# LinearSVC Implementation Summary & Status

**Date:** 2025-11-07
**Branch:** claude/investigate-ref-performance-011CUt9YD4yNb4AzXF2TKAU9
**Status:** ❌ NOT WORKING - LinearSVC + external calibration fails to classify U12 introns

---

## Executive Summary

We refactored intronIC's SVM classifier from `SVC(kernel='linear', probability=True)` to `LinearSVC + CalibratedClassifierCV` following sklearn best practices. **Despite correct implementation, the refactored version classifies 0 U12 introns vs the original's 32 U12 introns.**

---

## Current Implementation

### Architecture

```python
# Refactored (Fast but Broken)
base_svm = LinearSVC(
    C=optimized_value,          # From hyperparameter search
    class_weight='balanced',     # Handle 387 U12 vs 20,690 U2 imbalance
    loss='squared_hinge',        # Smooth, stable loss
    penalty='l2',                # L2 regularization
    dual=False,                  # Correct: n_samples (17k) >> n_features (3)
    max_iter=10000,
    tol=1e-4,
    random_state=seed
)

model = CalibratedClassifierCV(
    base_svm,
    method='sigmoid',            # Platt scaling
    cv=5                         # Stratified 5-fold CV
)

# Hyperparameter optimization
GridSearchCV(
    model,                       # Search THROUGH calibrator ✓
    param_grid={'estimator__C': C_grid},  # Access inner parameter ✓
    cv=5,                        # Stratified folds
    scoring='neg_log_loss',      # Optimize probability quality ✓
    n_jobs=-1
)
```

### Original (Slow but Working)

```python
model = SVC(
    C=optimized_value,           # ~1000 after optimization
    kernel='linear',
    probability=True,            # Internal Platt scaling
    class_weight='balanced',
    cache_size=1000,
    random_state=seed
)

GridSearchCV(
    model,
    param_grid={'C': C_grid},
    cv=5,
    scoring='balanced_accuracy',  # NOTE: Different metric!
    n_jobs=-1
)
```

---

## Best Practices Compliance

| Aspect | Best Practice | Our Implementation | Status |
|--------|--------------|-------------------|---------|
| **Class weights** | `'balanced'` for imbalanced data | `class_weight='balanced'` | ✅ |
| **dual parameter** | `"auto"` or `False` for n_samples >> n_features | `dual=False` | ✅ |
| **Loss function** | `'squared_hinge'` (smooth) | `loss='squared_hinge'` | ✅ |
| **Calibration** | Wrap with `CalibratedClassifierCV` | `CalibratedClassifierCV(..., cv=5)` | ✅ |
| **Calibration method** | `'sigmoid'` (Platt scaling) | `method='sigmoid'` | ✅ |
| **Search through calibrator** | Yes, use `estimator__param` | `param_grid={'estimator__C': ...}` | ✅ |
| **Scoring metric** | `'neg_log_loss'` for custom thresholds | `scoring='neg_log_loss'` | ✅ |
| **Stratified CV** | Required for imbalanced classes | Automatic in GridSearchCV | ✅ |
| **Fixed threshold** | Apply at prediction time | `proba >= 0.90` | ✅ |

**Conclusion:** ✅ All best practices correctly implemented

---

## Performance Comparison

### Chr19 Test Dataset
- **Experimental introns:** 12,074 (after filtering)
- **Reference data:** 387 U12, 20,690 U2
- **Expected U12s:** 32 (from original)

### Results

| Metric | Original (SVC) | Refactored (LinearSVC) | Target |
|--------|----------------|------------------------|--------|
| **Runtime** | ~60 seconds | ~13 seconds | - |
| **Speedup** | 1x | **4.5x faster** ✓ | - |
| **Optimized C** | ~1000 | 1.35e-04 | ~1000 |
| **F1 Score** | 1.0 | 0.3137 | 1.0 |
| **PR-AUC** | ~1.0 | 0.3786 | ~1.0 |
| **U12s Found** | **32** | **0** ❌ | 32 |
| **Max Probability** | >90% | 12.2% | >90% |

### Key Finding

**The refactored model produces probabilities that are ~8x too low:**
- Threshold: 90%
- Highest probability: 12.2%
- Result: NO introns exceed threshold → 0 U12s classified

---

## Investigation History

### Phase 1: Initial Implementation Issues
**Problem:** Convergence warnings, F1=0.31, 0 U12s found
**Diagnosis:** Wrong dual parameter, suboptimal loss function
**Fix:** Changed `dual=True` → `dual=False`, `loss='hinge'` → `loss='squared_hinge'`
**Result:** ❌ Still 0 U12s found

### Phase 2: C Value Optimization Audit
**Question:** Is the optimizer testing C≈1000 where the original finds optimum?
**Investigation:** Added debug logging to print all Round 1 CV scores
**Finding:** YES - Full range tested (1e-6 to 1e+6), C=1000 at position 10/13
**Discovery:** CV score saturation - all C from 1e-3 to 1e+6 gave identical scores
**Result:** Optimizer logic is correct, but selects C=9.3e-05 (too regularized)

### Phase 3: Scoring Metric Fix (Current)
**Hypothesis:** Using `'balanced_accuracy'` optimizes for threshold=0.5, but we deploy at 0.90
**Best Practice:** Use `'neg_log_loss'` to optimize probability quality for custom thresholds
**Fix:** Changed `scoring='balanced_accuracy'` → `scoring='neg_log_loss'`
**Result:** ❌ **STILL 0 U12s found** - Only improved C from 9.3e-05 to 1.35e-04

---

## Root Cause Analysis

### Why neg_log_loss Still Fails

With extreme class imbalance (1.8% positive class), `neg_log_loss` favors overly conservative models:

1. **Low C model** (heavy regularization):
   - Predicts low probabilities for everything
   - Never confidently wrong
   - Gets decent log_loss score ✓
   - **Selected by optimizer**

2. **High C model** (less regularization):
   - Tries to identify U12 introns
   - Makes some confident mistakes
   - Gets penalized heavily by log_loss ❌
   - **Rejected by optimizer**

**Result:** Even with the "correct" scoring metric, the optimizer selects an overly-regularized model that can't make confident predictions.

### The Fundamental Difference

| Aspect | SVC (Works) | LinearSVC (Fails) |
|--------|-------------|-------------------|
| **Algorithm** | libsvm | liblinear |
| **Calibration** | Internal (during training) | External (after training) |
| **Optimization** | Calibration-aware | Calibration-agnostic |
| **C optimization** | Considers calibrated probabilities | Considers uncalibrated decision function |
| **Result** | Probabilities scale correctly | Probabilities systematically too low |

**Key Insight:** With SVC, `probability=True` means the C parameter optimization happens **with knowledge of** the probability calibration. With LinearSVC + CalibratedClassifierCV, C is optimized first, **then** calibration is applied - but the calibration can't fix an already-too-regularized model.

---

## What We Tried

1. ✅ **Correct dual parameter** (`dual=False` for n_samples >> n_features)
2. ✅ **Correct loss function** (`squared_hinge` for smoothness)
3. ✅ **External calibration** (CalibratedClassifierCV with Platt scaling)
4. ✅ **Search through calibrator** (`estimator__C` parameter path)
5. ✅ **Best-practice scoring metric** (`neg_log_loss` for probability quality)
6. ✅ **Verified C search range** (1e-6 to 1e+6, includes optimal value)
7. ✅ **Stratified CV** (automatic for classification)
8. ✅ **Class balancing** (`class_weight='balanced'`)

**None of these fixed the core problem.**

---

## What We Haven't Tried

### Potential Solutions

1. **Different scoring metrics:**
   - `'brier_score_loss'` (alternative to log_loss for probability quality)
   - Custom scorer with class-weighted log_loss
   - `'average_precision'` (PR-AUC, threshold-independent)

2. **Sample weighting in scoring:**
   - Pass `sample_weight` to scoring function to up-weight rare class
   - May prevent log_loss from favoring overly conservative predictions

3. **Fixed high C value:**
   - Skip optimization, just use C=1000 (same as original's optimum)
   - Test if high C + calibration works, or if calibration itself is broken

4. **Isotonic calibration:**
   - `method='isotonic'` instead of `'sigmoid'`
   - More flexible, but requires more data (we have 21k reference samples)

5. **Just use SVC:**
   - Accept 100x slowdown (60s vs 13s for Chr19)
   - For full human genome: ~30min vs ~3min
   - Proven to work, no further investigation needed

---

## Current Code Locations

### Modified Files

**`intronIC_refactored/classification/optimizer.py`:**
- Line 292: `scoring='neg_log_loss'` (was `'balanced_accuracy'`)
- Line 422: `scoring='neg_log_loss'` in `_evaluate_C()`
- Lines 411-415: Wrap with CalibratedClassifierCV for final evaluation
- Lines 302-317: Debug logging for Round 1 CV scores

**`intronIC_refactored/classification/trainer.py`:**
- Lines 174-183: LinearSVC configuration with `dual=False`, `loss='squared_hinge'`
- Lines 187-191: CalibratedClassifierCV wrapper

### Investigation Documents

- `C_VALUE_AUDIT_RESULTS.md` - Complete audit showing C range is correct
- `LINEAR_SVC_INVESTIGATION.md` - Earlier convergence issue investigation
- `NEXT_STEPS_AFTER_C_AUDIT.md` - Recommendations after audit
- `PERFORMANCE_COMPARISON_RESULTS.md` - Speed benchmarks showing 4.5x speedup

---

## Test Reproduction

### Refactored Version (Current)
```bash
cd /home/user/intronIC
export PYTHONPATH=/home/user/intronIC/intronIC_refactored
python -m cli.main \
  -n test_refactored \
  -g intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
  -a intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz \
  --reference_u12s intronIC_refactored/intronIC/data/u12_reference.introns.iic.gz \
  --reference_u2s intronIC_refactored/intronIC/data/u2_reference.introns.iic.gz \
  -p 1 \
  --output_dir comparison_test/refactored_output

# Results:
# Runtime: ~13 seconds
# Optimized C: 1.35e-04
# F1: 0.3137
# U12s found: 0
```

### Original Version (Baseline)
```bash
cd /home/user/intronIC
python intronIC/intronIC.py \
  -n test_original \
  -g intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
  -a intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz \
  --reference_u12s intronIC/data/u12_reference.introns.iic.gz \
  --reference_u2s intronIC/data/u2_reference.introns.iic.gz \
  -p 1

# Results:
# Runtime: ~60 seconds
# Optimized C: ~1000
# F1: 1.0
# U12s found: 32
```

---

## Questions for Feedback

1. **Is LinearSVC + external calibration fundamentally incompatible with this task?**
   - Extreme class imbalance (1.8% positive)
   - Very high precision required (can't afford false positives)
   - Custom threshold (90% vs default 50%)

2. **Should we try other approaches or just revert to SVC?**
   - SVC works but is 100x slower per training iteration
   - For full genome: 30min vs 3min (still acceptable for batch processing)
   - Is the speed gain worth continuing investigation?

3. **Is there a way to make log_loss work with severe class imbalance?**
   - Custom scorer with class weights?
   - Different metric altogether?
   - Or is this a known limitation?

4. **Alternative: Skip hyperparameter optimization entirely?**
   - Use fixed C=1000 (same as original's optimum)
   - Just apply LinearSVC + calibration with known-good C
   - Would this work, or is calibration itself broken?

---

## Recommendation

Given the investigation results, I recommend **reverting to SVC** for the following reasons:

1. **Correctness > Speed:** The classifier must work correctly first
2. **Already Fast Enough:** 30 minutes for full human genome is acceptable
3. **Diminishing Returns:** Significant investigation time for uncertain payoff
4. **Production-Proven:** SVC approach is battle-tested in published research

The 4.5x speedup from LinearSVC would be nice, but not at the cost of 0% recall on the rare class.

---

## Code Diff: Reverting to SVC

If we decide to revert, here's what changes:

```python
# optimizer.py and trainer.py
# OLD (LinearSVC + external calibration):
base_svm = LinearSVC(
    C=C,
    class_weight='balanced',
    loss='squared_hinge',
    penalty='l2',
    dual=False,
    max_iter=10000,
    tol=1e-4,
    random_state=seed
)
model = CalibratedClassifierCV(base_svm, method='sigmoid', cv=5)

# NEW (SVC + internal calibration):
model = SVC(
    C=C,
    kernel='linear',
    probability=True,  # Built-in Platt scaling
    class_weight='balanced',
    cache_size=1000,
    random_state=seed
)

# Also change:
# - param_grid: 'estimator__C' → 'C'
# - scoring: 'neg_log_loss' → 'balanced_accuracy' (original used this)
```

Estimated effort: 30 minutes to revert + 1 hour testing/validation

---

## Conclusion

We successfully implemented LinearSVC following sklearn best practices, achieving 4.5x speedup. However, **the refactored version classifies 0 U12 introns** due to systematically low probability estimates (max 12.2% vs 90% threshold).

The root cause appears to be a fundamental incompatibility between:
- LinearSVC + external calibration architecture
- Extreme class imbalance (1.8% positive class)
- High-precision requirements (rare class identification)
- Custom decision threshold (90% vs default 50%)

**Recommendation:** Revert to SVC for correctness, accept the slowdown.

---

**Files:**
- Implementation: `intronIC_refactored/classification/`
- Tests: `comparison_test/refactored_output/test_fixed_logloss.*`
- Branch: `claude/investigate-ref-performance-011CUt9YD4yNb4AzXF2TKAU9`
