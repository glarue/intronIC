# Scaler Centering Fix - Implementation Plan

**Goal:** Restore cross-species performance by adopting centered scaling with probability calibration and optional refinements.

**Expected Outcome:** C. elegans FP rate <10 (currently 130), consistent probability thresholds across species.

---

## Background

**Problem:** Zero-anchored scaling (no centering) causes 130x regression in cross-species FP rate
- Original v1.5.1: ~1 FP on C. elegans (StandardScaler WITH centering)
- Current: ~130 FPs on C. elegans (ZeroAnchoredRobustScaler WITHOUT centering)

**Root Cause:** Without centering, composition bias inflates z-scores unboundedly across species.

**Solution:** Return to centered scaling + add probability calibration for consistent thresholds.

---

## Three-Phase Implementation

### Phase 1: Centered Scaling with Base Features

**Goal:** Replace zero-anchored scaling with RobustScaler (centered)

**Architecture Changes:**
```python
# OLD (current)
Pipeline([
    ('scale', ZeroAnchoredRobustScaler()),  # NO centering
    ('augment', BothEndsStrongTransformer()),
    ('svc', LinearSVC(C=optimized))
])

# NEW (Phase 1)
Pipeline([
    ('scale', RobustScaler()),  # WITH centering (default)
    ('svc', LinearSVC(C=optimized))
])
```

**Files to Modify:**

1. **classification/trainer.py**
   - Import: `from sklearn.preprocessing import RobustScaler`
   - Replace: `ZeroAnchoredRobustScaler()` → `RobustScaler()`
   - Remove: BothEndsStrongTransformer (test with 3 base features first)
   - Keep: Rest of training logic unchanged

2. **classification/optimizer.py**
   - Import: `from sklearn.preprocessing import RobustScaler`
   - Replace: `ZeroAnchoredRobustScaler()` → `RobustScaler()`
   - Remove: BothEndsStrongTransformer
   - Keep: Grid search logic unchanged

3. **classification/predictor.py**
   - Verify: Works with RobustScaler (should be transparent)
   - No changes needed (uses fitted pipeline)

**Testing:**
- Train on human reference data
- Validate on human test set (expect similar metrics)
- **Critical test:** Validate on C. elegans (expect <10 FPs)

**Success Criteria:**
- ✅ C. elegans FP count <10 (down from 130)
- ✅ Human test set F1 score ≥ original
- ✅ Training completes without errors

---

### Phase 2: Probability Calibration

**Goal:** Wrap LinearSVC with calibration for consistent probability thresholds across species

**Architecture Changes:**
```python
# NEW (Phase 2)
from sklearn.calibration import CalibratedClassifierCV
from sklearn.svm import LinearSVC

Pipeline([
    ('scale', RobustScaler()),
    ('svc', CalibratedClassifierCV(
        LinearSVC(C=optimized),
        method='sigmoid',  # Platt scaling
        cv=5  # 5-fold CV for calibration
    ))
])
```

**Why Calibration Matters:**
- LinearSVC's decision_function outputs are NOT probabilities
- Platt scaling fits a sigmoid to map margins → well-calibrated probabilities
- Makes p ≥ 0.90 threshold behave consistently across species

**Files to Modify:**

1. **classification/trainer.py**
   - Import: `from sklearn.calibration import CalibratedClassifierCV`
   - Modify: Wrap LinearSVC with CalibratedClassifierCV
   - Update: Ensure probability outputs are used (not decision_function)
   - Add: Log-loss evaluation on validation set

2. **classification/optimizer.py**
   - Modify: Grid search must handle calibrated classifier
   - Note: Calibration happens AFTER hyperparameter search
   - Process:
     1. Grid search LinearSVC to find optimal C
     2. Take best C, wrap with CalibratedClassifierCV
     3. Fit calibrated model on full training set

3. **classification/predictor.py**
   - Verify: Uses `.predict_proba()` method
   - Update: Handle CalibratedClassifierCV output format
   - Ensure: Probabilities are correctly extracted

4. **cli/main.py** (if needed)
   - Verify: Threshold logic uses calibrated probabilities
   - Update: Logging to show calibration info

**Testing:**
- Train calibrated model on human
- **Calibration check:** Plot reliability diagram (predicted p vs observed frequency)
- Validate on C. elegans (expect similar FP count as Phase 1, but better probability calibration)
- **Critical test:** Check if p ≥ 0.90 threshold is consistent across human test set vs C. elegans

**Success Criteria:**
- ✅ Calibrated probabilities match empirical frequencies
- ✅ Threshold p ≥ 0.90 behaves consistently across species
- ✅ C. elegans FP count remains low
- ✅ Human test set performance unchanged

---

### Phase 3: Advanced Refinements

**Goal:** Add margin alignment, prior shift, and scaler serialization for production robustness

#### 3A. Margin Alignment

**Purpose:** Correct residual distribution shift after scaling

**Approach:**
```python
# During training (human)
margins_train = svc.decision_function(X_train_scaled)
u2_margins = margins_train[y_train == 0]  # Negative class
mu_H_neg = np.median(u2_margins)
s_H_neg = stats.median_abs_deviation(u2_margins)
# Save: (mu_H_neg, s_H_neg) with model

# During inference (target species T)
margins_T = svc.decision_function(X_T_scaled)
mu_T = np.median(margins_T)  # Dominated by U2s
s_T = stats.median_abs_deviation(margins_T)

# Affine correction
a = s_H_neg / s_T
c = mu_H_neg - a * mu_T
margins_T_corrected = a * margins_T + c

# Apply calibrator to corrected margins
probs_T = calibrator.predict_proba_from_decision(margins_T_corrected)
```

**Files to Modify:**

1. **classification/trainer.py**
   - Add: Compute and save (mu_H_neg, s_H_neg) after training
   - Store: In model metadata

2. **classification/predictor.py**
   - Add: `apply_margin_alignment()` function
   - Add: `--margin_alignment` flag support
   - Compute: Target species margin stats
   - Apply: Affine correction before calibration

**Testing:**
- Compare C. elegans predictions with/without margin alignment
- Measure: Change in FP count
- Validate: On species with known U12 prevalence

#### 3B. Prior Shift

**Purpose:** Adjust probabilities for species with different U12 prevalence

**Bayes Formula:**
```
p(U12|x, π_T) = p(U12|x, π_H) * (π_T / π_H) / Z
where:
  π_H = human U12 prevalence (~0.005)
  π_T = target species U12 prevalence (user-provided or estimated)
  Z = normalization constant
```

**Use Cases:**
- C. elegans: π_T = 0 (U12-free) → all probabilities → 0
- Fungi: π_T = 0.0001 (very rare) → shift probabilities down
- Mammals: π_T = 0.005 (similar to human) → no shift

**Files to Modify:**

1. **classification/predictor.py**
   - Add: `apply_prior_shift(probs, pi_H, pi_T)` function
   - Add: `--u12_prior` flag to specify π_T

2. **cli/args.py**
   - Add: `--u12_prior` argument (float, default=None)
   - Add: `--margin_alignment` argument (boolean, default=False)

**Testing:**
- C. elegans with π_T=0: All predictions → U2
- Synthetic species with π_T=0.01: Check probability adjustment

#### 3C. Scaler Serialization

**Purpose:** Ensure reproducible results by reusing same scaler per species

**Approach:**
```python
# First run on species S
scaler_S = RobustScaler()
scaler_S.fit(introns_S_raw_features)
save_scaler(scaler_S, f'scalers/{species_name}_scaler.pkl')

# Subsequent runs on species S
scaler_S = load_scaler(f'scalers/{species_name}_scaler.pkl')
z_scores = scaler_S.transform(introns_S_raw_features)
```

**Files to Modify:**

1. **utils/model_io.py** (or create new file)
   - Add: `save_scaler(scaler, path)` function
   - Add: `load_scaler(path)` function
   - Store: Scaler params (center_, scale_) in JSON or pickle

2. **classification/predictor.py**
   - Add: `--scaler` argument to load pre-fitted scaler
   - Add: Auto-save scaler after first fit

3. **cli/args.py**
   - Add: `--scaler` argument (path to saved scaler)
   - Add: `--save_scaler` argument (path to save fitted scaler)

**Testing:**
- Fit scaler on C. elegans, save
- Load saved scaler, verify predictions identical
- Cross-validate: Different random samples should give same scaler params

**Success Criteria (Phase 3 Overall):**
- ✅ Margin alignment reduces FP variance across species
- ✅ Prior shift correctly handles U12-free genomes
- ✅ Scaler serialization ensures reproducibility
- ✅ All features optional (flags) and backward-compatible

---

## Implementation Order

### Step 1: Phase 1 - Centered Scaling
1. Modify `classification/optimizer.py`
2. Modify `classification/trainer.py`
3. Test train on human
4. Validate on C. elegans

**Checkpoint:** If C. elegans FPs <10, proceed to Phase 2

### Step 2: Phase 2 - Calibration
1. Modify `classification/optimizer.py` (calibration wrapper)
2. Modify `classification/trainer.py` (calibrated training)
3. Modify `classification/predictor.py` (use calibrated probs)
4. Test calibration quality (reliability diagram)
5. Validate threshold consistency

**Checkpoint:** If calibration is well-behaved, proceed to Phase 3

### Step 3: Phase 3A - Margin Alignment
1. Modify `classification/trainer.py` (save margin stats)
2. Modify `classification/predictor.py` (apply alignment)
3. Add CLI arguments
4. Test on multiple species

### Step 4: Phase 3B - Prior Shift
1. Modify `classification/predictor.py` (Bayes correction)
2. Add CLI arguments
3. Test on U12-free genomes

### Step 5: Phase 3C - Scaler Serialization
1. Create/modify `utils/model_io.py`
2. Modify `classification/predictor.py`
3. Add CLI arguments
4. Test reproducibility

---

## Files Summary

**Phase 1:**
- `classification/optimizer.py` (replace scaler)
- `classification/trainer.py` (replace scaler, remove augmentation)
- `classification/predictor.py` (verify compatibility)

**Phase 2:**
- `classification/optimizer.py` (calibration in grid search)
- `classification/trainer.py` (wrap with CalibratedClassifierCV)
- `classification/predictor.py` (use calibrated probabilities)

**Phase 3:**
- `classification/trainer.py` (margin stats)
- `classification/predictor.py` (alignment, prior shift)
- `utils/model_io.py` (scaler serialization)
- `cli/args.py` (new arguments)

---

## Testing Strategy

**After Each Phase:**
1. Train model on human reference data
2. Evaluate on human test set (F1, PR-AUC)
3. **Critical:** Validate on C. elegans (FP count)
4. Compare to previous phase results

**Final Validation (After Phase 3):**
- C. elegans: 0 expected, <10 acceptable
- D. melanogaster: Known U12s, check recall
- Human (held-out chr): Check generalization
- Multiple runs: Check reproducibility

**Regression Tests:**
- Ensure no performance degradation on human
- Verify output format compatibility
- Check memory usage unchanged

---

## Success Metrics

**Phase 1:**
- C. elegans FP rate: <10 (vs current 130) ✅
- Human F1 score: ≥ current baseline ✅

**Phase 2:**
- Calibration ECE (Expected Calibration Error): <0.05 ✅
- Threshold consistency: p ≥ 0.90 similar FP rate across species ✅

**Phase 3:**
- Margin alignment: Reduces FP variance >20% ✅
- Prior shift: C. elegans with π=0 → 0 predictions ✅
- Serialization: Identical results on reload ✅

---

## Notes

- **Backward compatibility:** Phase 3 features are optional (CLI flags)
- **Documentation:** Update CLAUDE.md after Phase 2
- **Testing data:** Keep C. elegans as primary cross-species test
- **Performance:** Monitor training time (calibration adds CV overhead)

---

**Next Step:** Begin Phase 1 implementation
