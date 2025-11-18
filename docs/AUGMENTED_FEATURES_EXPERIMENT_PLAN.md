# Augmented Features Experiment Plan

**Date:** 2025-11-17
**Objective:** Compare baseline (3D) vs augmented (7D) feature sets to determine if composite features improve classification
**Status:** Planning phase - awaiting completion of scaler serialization training

---

## Background

### Current State
- **Model A (3D baseline)**: Currently training as `homo_sapiens_scaler.model.pkl`
  - Features: five_z, bp_z, three_z
  - Uses RobustScaler with centering (median/IQR normalization)
  - Scaler serialization enabled for cross-species deployment

- **Model B (7D augmented)**: Fully implemented but disabled
  - BothEndsStrongTransformer commented out in optimizer.py and trainer.py
  - Would add: min_all, neg_absdiff_5_bp, neg_absdiff_5_3, neg_absdiff_bp_3
  - Designed to suppress "one-end-strong" false positives

### Motivation
Expert recommendation: Don't drop augmented features on theory alone. Test empirically.

**Key questions:**
1. Do composite features improve human CV metrics?
2. Do they reduce FPs in U12-free species (C. elegans)?
3. Are they genuinely useful or redundant (L1 sparsity test)?

---

## Experimental Design

### Model Variants

#### Model A: Baseline (3D)
**Features:**
- five_z: 5' splice site z-score
- bp_z: Branch point z-score
- three_z: 3' splice site z-score

**Pipeline:**
```python
Pipeline([
    ('scale', RobustScaler()),
    ('svm', LinearSVC(C=optimized, penalty='l2'))
])
```

**Status:** Training now (homo_sapiens_scaler)

---

#### Model B: Augmented (7D)
**Features:**
- five_z, bp_z, three_z (base)
- min_all = min(five_z, bp_z, three_z) - requires all signals strong
- neg_absdiff_5_bp = -|five_z - bp_z| - penalizes 5'/BP imbalance
- neg_absdiff_5_3 = -|five_z - three_z| - penalizes 5'/3' imbalance
- neg_absdiff_bp_3 = -|bp_z - three_z| - penalizes BP/3' imbalance

**Pipeline:**
```python
Pipeline([
    ('scale', RobustScaler()),
    ('augment', BothEndsStrongTransformer(include_pairwise_mins=False)),
    ('svm', LinearSVC(C=optimized, penalty='l2'))
])
```

**Status:** Ready to implement

**Rationale:** Composite features explicitly encode:
- "All ends must be strong" (min_all)
- "Signals must be consistent" (neg_absdiff_*)
- Helps linear SVM reject one-end-strong FPs

---

#### Model C: Augmented + L1 Sparsity (7D with automatic feature selection)
**Features:** Same as Model B

**Pipeline:**
```python
Pipeline([
    ('scale', RobustScaler()),
    ('augment', BothEndsStrongTransformer(include_pairwise_mins=False)),
    ('svm', LinearSVC(C=optimized, penalty='l1', dual=False))
])
```

**Status:** Requires optimizer modification

**Rationale:** L1 penalty drives redundant features to zero
- If min_all / neg_absdiff_* survive → genuinely useful
- If L1 zeros them out with no performance drop → redundant

**Note:** Requires adding L1 support to optimizer grid search

---

## Evaluation Metrics

### 1. Human CV Metrics (Training Species)
**Primary:**
- Balanced accuracy (already optimizing for this)
- PR-AUC (area under precision-recall curve)
- F1 score

**Precision at Target Recall:**
- Precision when recall = 0.95 (catch 95% of U12s)
- Precision when recall = 0.99 (catch 99% of U12s)

**Calibration Quality:**
- Log-loss (lower = better calibrated)
- Expected Calibration Error (ECE)
- Maximum Calibration Error (MCE)
- Brier score

**Comparison:** 7-fold nested CV results (already implemented)

---

### 2. C. elegans Metrics (U12-Free Species)
**Primary metric:**
- False positive count at p ≥ 0.90 threshold

**Current baseline (domain adaptation):** 130 FPs
**Expected with scaler fix (Model A):** 0-1 FPs
**Goal with Model B:** ≤ Model A (ideally 0)

**Secondary metrics:**
- FP count at p ≥ 0.95 (stricter threshold)
- FP count at p ≥ 0.99 (very strict)
- Distribution of FP probabilities

---

### 3. Other Species (If Available)
**U12-containing genomes:**
- Recall: fraction of true U12s found
- Precision: fraction of predictions that are true U12s
- F1 score

**Candidates (if test data available):**
- Mouse (Mus musculus)
- Zebrafish (Danio rerio)
- Arabidopsis thaliana (plant, different U12 distribution)

---

### 4. Feature Importance Analysis (Model C only)
**From L1 model:**
- Which features have non-zero coefficients?
- Magnitude of coefficients (proxy for importance)
- How many features zeroed out?

**Interpretation:**
- Features surviving L1 = genuinely useful
- Features zeroed without performance loss = redundant

---

## Implementation Steps

### Phase 1: Complete Baseline (In Progress)
**Status:** Model A training (homo_sapiens_scaler)

✅ **Step 1.1:** Train baseline model with scaler serialization
- [x] Implement scaler serialization (cli/main.py:1467-1477)
- [x] Add normalizer_mode CLI flag (cli/args.py:356-365)
- [x] Update classification logic (cli/main.py:1085-1128)
- [x] Start training (ETA: ~30-35 minutes)

⬜ **Step 1.2:** Validate baseline on C. elegans
```bash
pixi run intronIC classify \
  -g archive/test_outputs/celegans/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz \
  -a archive/test_outputs/celegans/Caenorhabditis_elegans.WBcel235.62.gff3.gz \
  -n caenorhabditis_elegans_baseline \
  -o celegans_test \
  --model homo_sapiens_scaler.model.pkl \
  --normalizer_mode human \
  -p 10
```

**Expected result:** 0-1 FPs (vs 130 with domain adaptation)

⬜ **Step 1.3:** Record baseline metrics
- Human CV: balanced_accuracy, PR-AUC, log-loss, ECE
- C. elegans: FP count at p ≥ 0.90

---

### Phase 2: Implement Augmented Features
**Status:** Awaiting Phase 1 completion

⬜ **Step 2.1:** Enable BothEndsStrongTransformer

**File:** `classification/optimizer.py`
- Uncomment line 45: `from classification.transformers import BothEndsStrongTransformer`
- Add transformer to pipeline construction (~line 200-250):
```python
Pipeline([
    ('scale', RobustScaler()),
    ('augment', BothEndsStrongTransformer(include_pairwise_mins=False)),
    ('svm', LinearSVC(...))
])
```

**File:** `classification/trainer.py`
- Uncomment line 32: `from classification.transformers import BothEndsStrongTransformer`
- Add transformer to pipeline (~line 100-150) - mirror optimizer changes

⬜ **Step 2.2:** Train augmented model
```bash
pixi run intronIC train \
  -n homo_sapiens_augmented \
  -o . \
  --n_optimization_rounds 3 \
  --n_cv_folds 7 \
  -p 10
```

**Expected changes:**
- Feature space: 3D → 7D
- C optimization will adjust for new feature space
- Training time: similar to baseline (~30-35 min)

⬜ **Step 2.3:** Validate augmented model on C. elegans
```bash
pixi run intronIC classify \
  -g archive/test_outputs/celegans/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz \
  -a archive/test_outputs/celegans/Caenorhabditis_elegans.WBcel235.62.gff3.gz \
  -n caenorhabditis_elegans_augmented \
  -o celegans_test \
  --model homo_sapiens_augmented.model.pkl \
  --normalizer_mode human \
  -p 10
```

**Hypothesis:** FP count ≤ baseline (composite features suppress one-end-strong FPs)

⬜ **Step 2.4:** Record augmented metrics
- Human CV: balanced_accuracy, PR-AUC, log-loss, ECE
- C. elegans: FP count at p ≥ 0.90

---

### Phase 3: Compare Baseline vs Augmented
**Status:** Awaiting Phase 2 completion

⬜ **Step 3.1:** Compare human CV metrics

**Metrics to compare:**
| Metric | Model A (3D) | Model B (7D) | Winner |
|--------|--------------|--------------|--------|
| Balanced accuracy | ? | ? | ? |
| PR-AUC | ? | ? | ? |
| F1 score | ? | ? | ? |
| Log-loss | ? | ? | ? |
| ECE | ? | ? | ? |
| Precision @ Recall=0.95 | ? | ? | ? |
| Precision @ Recall=0.99 | ? | ? | ? |

**Decision criteria:**
- If Model B improves any metric without hurting others → use B
- If Model B same as A → proceed to Phase 4 (L1 test)
- If Model A better → stick with baseline

⬜ **Step 3.2:** Compare C. elegans FP counts

| Threshold | Model A (3D) | Model B (7D) | Winner |
|-----------|--------------|--------------|--------|
| p ≥ 0.90 | ? | ? | ? |
| p ≥ 0.95 | ? | ? | ? |
| p ≥ 0.99 | ? | ? | ? |

**Decision criteria:**
- If Model B reduces FPs → strong evidence for augmented features
- If Model A = Model B → augmented features not harmful (neutral)
- If Model A < Model B → investigate why (unexpected)

⬜ **Step 3.3:** Document findings
- Create comparison plots (scatter, calibration curves)
- Analyze FP cases (which introns differ between A and B?)
- Draft summary document

---

### Phase 4: L1 Sparsity Test (Optional)
**Status:** Conditional on Phase 3 results
**Trigger:** If Phase 3 shows Model B ≈ Model A (ambiguous whether features useful)

⬜ **Step 4.1:** Add L1 penalty support to optimizer

**File:** `classification/optimizer.py`
- Add penalty parameter to grid search space:
```python
penalty_options = ['l2', 'l1']  # Add L1 option
dual_options = [True, False]     # dual=False required for L1
```
- Update optimizer to handle penalty/dual combinations
- Note: L1 + dual=True is invalid, handle this constraint

⬜ **Step 4.2:** Train L1 model
```bash
pixi run intronIC train \
  -n homo_sapiens_l1_sparse \
  -o . \
  --n_optimization_rounds 3 \
  --n_cv_folds 7 \
  -p 10
```

**Expected behavior:**
- Optimizer will compare L1 vs L2
- L1 will drive some feature weights to exactly zero
- Best penalty chosen via CV performance

⬜ **Step 4.3:** Extract feature coefficients

**From trained model:**
```python
import joblib
model_data = joblib.load('homo_sapiens_l1_sparse.model.pkl')
first_model = model_data['ensemble'].models[0].model

# Get base estimator from calibrated classifier
base = first_model.calibrated_classifiers_[0].estimator

# Extract coefficients
coefs = base.named_steps['svm'].coef_[0]
feature_names = ['five_z', 'bp_z', 'three_z', 'min_all',
                 'neg_absdiff_5_bp', 'neg_absdiff_5_3', 'neg_absdiff_bp_3']

# Analyze sparsity
for name, coef in zip(feature_names, coefs):
    status = "ACTIVE" if abs(coef) > 1e-6 else "ZEROED"
    print(f"{name:20s} {coef:+.4f}  {status}")
```

⬜ **Step 4.4:** Interpret results

**If composite features survive L1:**
- min_all, neg_absdiff_* have non-zero weights
- → Features are genuinely useful, not redundant
- → Use augmented model (Model B or C)

**If composite features zeroed by L1:**
- min_all, neg_absdiff_* driven to zero
- → Features are redundant with base features
- → Use baseline model (Model A)

**If mixed:**
- Some composite features active, others zeroed
- → Simplify transformer to only include active features
- → Retrain with reduced feature set

⬜ **Step 4.5:** Final model selection
- Choose Model A, B, or C based on all evidence
- Document decision rationale
- Update default configuration

---

## Decision Matrix

### Model Selection Criteria

| Evidence | Recommendation |
|----------|---------------|
| Model B improves human CV + reduces C. elegans FPs | **Use Model B (7D augmented)** |
| Model B improves human CV, same C. elegans FPs | **Use Model B (7D augmented)** |
| Model B same human CV + reduces C. elegans FPs | **Use Model B (7D augmented)** |
| Model B same human CV + same C. elegans FPs | **Run Phase 4 (L1 test)** |
| Model B worse human CV or worse C. elegans FPs | **Use Model A (3D baseline)** |
| L1 keeps composite features active | **Use Model C (7D with L1)** |
| L1 zeros composite features | **Use Model A (3D baseline)** |

---

## Expected Outcomes

### Best Case (Model B wins)
- Augmented features reduce C. elegans FPs from 0-1 to 0
- Human CV metrics improve (especially precision at high recall)
- Clear evidence that composite features suppress one-end-strong FPs
- **Action:** Deploy Model B as production default

### Neutral Case (Model A ≈ Model B)
- Similar performance on all metrics
- Composite features not harmful, but also not clearly helpful
- **Action:** Run L1 test (Phase 4) to determine if features redundant
- If L1 zeros them → stick with Model A (simpler)
- If L1 keeps them → use Model B (may help on other species)

### Baseline Wins (Model A better)
- Model A has better metrics than Model B
- Unexpected, but possible if composite features add noise
- **Action:** Investigate why (examine feature distributions, correlations)
- Stick with Model A unless issue identified and fixed

---

## Risk Assessment

### Low Risk
- All models use same training data and optimization procedure
- Only difference is feature engineering (3D vs 7D)
- BothEndsStrongTransformer already implemented and tested
- Scaler serialization works for any pipeline

### Potential Issues
1. **Training time:** 7D model slightly slower than 3D (more features to process)
   - Mitigation: Minimal impact (~5-10% longer)

2. **Overfitting risk:** More features = more parameters
   - Mitigation: L2/L1 regularization + 7-fold CV prevents overfitting

3. **Cross-species generalization:** Composite features tuned on human
   - Mitigation: That's exactly what we're testing with C. elegans

---

## Timeline Estimate

**Phase 1 (Baseline validation):**
- Training: 30-35 min (in progress)
- C. elegans classification: ~5 min
- Metric extraction: 5 min
- **Total:** ~45 min (mostly automated)

**Phase 2 (Augmented features):**
- Code changes: 10 min
- Training: 30-35 min
- C. elegans classification: ~5 min
- Metric extraction: 5 min
- **Total:** ~55 min

**Phase 3 (Comparison):**
- Metric comparison: 15 min
- Plot generation: 10 min
- FP case analysis: 20 min
- Documentation: 15 min
- **Total:** ~60 min

**Phase 4 (L1 test, if needed):**
- Code changes: 20 min
- Training: 30-35 min
- Coefficient analysis: 15 min
- Documentation: 10 min
- **Total:** ~80 min

**Grand Total:** 2.5-4 hours depending on whether Phase 4 needed

---

## Success Criteria

### Minimum Success
- Complete Phase 1-3 comparison
- Clear decision on Model A vs Model B based on metrics
- Document findings in comparison report

### Full Success
- Complete all 4 phases
- Definitive answer on whether composite features useful
- Deploy best model with confidence
- Publish methodology as template for future feature engineering

---

## Next Steps

1. **Immediate:** Wait for `homo_sapiens_scaler` training to complete
2. **Step 1.2:** Validate baseline on C. elegans (expect 0-1 FPs)
3. **Step 1.3:** Record baseline metrics
4. **Decision point:** Proceed to Phase 2 (implement augmented features)

---

## References

- **BothEndsStrongTransformer implementation:** `classification/transformers.py`
- **FP reduction research:** `FP_ROOT_CAUSE_ANALYSIS.md`
- **Scaler serialization plan:** `SCALER_SERIALIZATION_PLAN.md`
- **Expert recommendations:** This document, "Background" section

---

**Status:** ⬜ Ready to begin Phase 1.2 after current training completes
**Estimated Time to First Results:** ~1 hour (baseline validation)
**Risk:** Low (incremental, empirical approach)
