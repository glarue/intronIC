# intronIC Scoring & Classification Pipeline - Analysis & Recommendations

**Date:** 2025-11-03
**Analyst:** Claude (Anthropic)
**Original Code:** `/home/glarue/code/intronIC/intronIC/intronIC.py`
**Focus:** Lines 3655-5900 (Scoring, Training, Classification)

---

## Executive Summary

The original intronIC scoring and classification pipeline is **generally well-designed** with several best practices for machine learning. However, there are **two significant issues** that should be addressed in the refactored version:

1. **Post-classification z-score re-normalization** (data leakage in output, not classification)
2. **Redundant normalization in recursive scoring mode**

The core ML approach (ensemble SVM with U2 subsampling, proper train/test splits, cross-validation for hyperparameter tuning) is **sound and appropriate** for this problem domain.

---

## Table of Contents

1. [Pipeline Overview](#pipeline-overview)
2. [What's Done Right](#whats-done-right)
3. [Issues & Concerns](#issues--concerns)
4. [Detailed Analysis by Component](#detailed-analysis-by-component)
5. [Recommendations for Refactor](#recommendations-for-refactor)
6. [Test Strategy Recommendations](#test-strategy-recommendations)

---

## Pipeline Overview

### High-Level Flow

```
1. Score experimental introns with PWMs → raw scores
2. Score reference introns (U2/U12) with PWMs → raw scores

3. Fit StandardScaler on reference raw scores ONLY
4. Transform both reference and experimental raw scores → z-scores

5. Optimize SVM hyperparameter C via nested CV
   - Outer loop: Multiple training iterations with U2 subsampling
   - Inner loop: GridSearchCV with 5-fold CV
   - Refinement: 5 rounds of geometric grid search

6. Train ensemble of SVMs (10+ models)
   - Each model uses different U2 subsample (same size as U12 set)
   - 80/20 train/test split per model
   - Stratified sampling

7. Classify experimental introns
   - Apply all models
   - Average U12 probabilities (F1-weighted)
   - Assign type based on threshold (default 90%)

8. [ISSUE] Re-normalize z-scores using ALL data
9. Write outputs
```

---

## What's Done Right

### ✅ 1. Proper Train/Test Separation for Normalization

**Location:** `intronIC.py:3725-3732`

```python
# scale based on the training set only, to avoid
# test set data leak
scale_vector = get_score_vector(
    raw_refs, score_names=raw_score_names)

# make a scaler to adjust raw scores --> z-scores
score_scaler = preprocessing.StandardScaler().fit(scale_vector)

scored_refs = scale_scores(raw_refs, score_scaler)
# ...
scored_introns = scale_scores(raw_introns, score_scaler)
```

**Why This is Good:**
- Fits scaler on **reference set only** (training data)
- Transforms both reference and experimental sets with same scaler
- Prevents data leakage from experimental introns into normalization parameters
- This is a **best practice** for ML pipelines

**Grade:** A+ ✅

---

### ✅ 2. Stratified Train/Test Split

**Location:** `intronIC.py:5398-5404`

```python
(train_scores, test_scores,
 train_labels, test_labels) = train_test_split(
    input_feats,
    input_labels,
    train_size=train_fraction,  # 0.80
    stratify=input_labels,       # KEY: stratified!
    random_state=seed)
```

**Why This is Good:**
- 80/20 split is standard and appropriate
- **Stratification** ensures equal class proportions in train/test
- Critical for imbalanced datasets (U12 is ~0.5% of introns)
- Fixed random seed for reproducibility

**Grade:** A+ ✅

---

### ✅ 3. Ensemble Approach with Class Balancing

**Location:** `intronIC.py:5356-5366`

```python
if iterations:
    subset_size = len(u12_vector)  # ~500 U12s
    # ...
    u2_sets = np.random.choice(
        len(u2s), size=(iterations, subset_size), replace=False)
    u2_sets = u2s[u2_sets]
```

**Why This is Good:**
- U2 introns vastly outnumber U12 introns (~200:1 ratio)
- Subsampling U2 to match U12 count prevents majority class dominance
- Training multiple models with different subsamples creates an ensemble
- Ensemble reduces variance and improves robustness
- This is a **best practice** for extreme class imbalance

**Alternative Approaches Considered:**
- SMOTE (synthetic oversampling) - not appropriate for real genomic sequences
- Class weights only - used as backup when iterations=0, but ensemble is better

**Grade:** A ✅

---

### ✅ 4. Nested Cross-Validation for Hyperparameter Tuning

**Location:** `intronIC.py:5431-5528`

**Structure:**
```
For each optimization round (5 rounds):
    For each U2 subsample (10 iterations):
        Split into train/test (80/20, stratified)
        GridSearchCV on training set:
            5-fold cross-validation
            Test multiple C values
            Select best by balanced_accuracy
        Evaluate on held-out test set
        Record F1, PR-AUC
    Average best C values across iterations
    Refine search grid around best C
```

**Why This is Good:**
- **Nested CV** structure prevents information leakage
- Inner CV (GridSearchCV) selects hyperparameters
- Outer test set evaluates generalization
- Multiple iterations with different data splits reduce overfitting
- Geometric refinement (5 rounds) efficiently explores parameter space
- `balanced_accuracy` scorer is appropriate for imbalanced data

**Grade:** A ✅

---

### ✅ 5. Conservative Rank-1 Parameter Averaging

**Location:** `intronIC.py:5306-5322`

```python
def rank1_param_avg(models):
    """
    Returns geometric mean of rank-1 parameter values across all models
    """
    params = defaultdict(list)
    for m in models:
        for p, v in m.best_params_.items():
            best_params = rank_ones(m.cv_results_, p)  # Get ALL rank-1 values
            params[p].extend(best_params)
    refined_params = {}
    for p, v in params.items():
        best_avg_param = pystats.gmean(v)  # Geometric mean
        refined_params[p] = best_avg_param
```

**Why This is Good:**
- Sklearn's `best_params_` returns **first** rank-1 value (arbitrary)
- This code finds **all** rank-1 parameters (tied for best)
- Takes geometric mean to get conservative middle value
- For C parameter (soft margin), higher C = more conservative
- Geometric mean is appropriate for log-scale parameters

**Rationale from Comments:**
> "a lower value for C is less conservative than a larger value. Therefore, in this specific application (with a single hyperparameter) taking the average of all equally-well-performing values seems appropriate."

**Grade:** A ✅ (thoughtful, documented approach)

---

### ✅ 6. F1-Weighted Ensemble Predictions

**Location:** `intronIC.py:5651-5687`

```python
def average_svm_score_info(probabilities, labels, distances, weights=None):
    # ...
    u12_probs = np.array(probabilities)
    distances = np.array(distances)
    if weights is not None:
        u12_probs = u12_probs * weights  # Weight by F1 scores
        distances = distances * weights
    avg_u12_prob = np.mean(u12_probs)
    avg_distance = np.mean(distances)
```

**Why This is Good:**
- Models with higher F1 scores contribute more to predictions
- Better-performing models have stronger influence
- Simple and interpretable weighting scheme

**Grade:** A ✅

---

### ✅ 7. Appropriate Evaluation Metrics

**Metrics Used:**
- **F1 score** - balances precision and recall
- **Precision-Recall AUC** - better than ROC-AUC for imbalanced data
- **Balanced accuracy** - equal weight to both classes
- **Classification report** - detailed per-class metrics

**Why This is Good:**
- ROC-AUC can be misleading with severe class imbalance
- PR-AUC focuses on performance on minority class (U12)
- These are **best practices** for imbalanced classification

**Grade:** A+ ✅

---

## Issues & Concerns

### ❌ ISSUE 1: Post-Classification Z-Score Re-normalization (Data Leakage)

**Location:** `intronIC.py:5247-5251`

**Code:**
```python
# re-scale introns so z-scores are based on entire data, not just
# training set (the latter is required for the scoring process, however)
final_score_vector = get_score_vector(finalized_introns, raw_score_names)
final_z_scaler = preprocessing.StandardScaler().fit(final_score_vector)
finalized_introns = scale_scores(finalized_introns, final_z_scaler)
```

**Problem:**
1. After classification is complete, z-scores are **re-normalized** using ALL data
2. This includes experimental introns that were just classified
3. The output z-scores are **different** from the z-scores used for training/classification
4. This is a form of **data leakage** - the z-scores now contain information from the test set

**Impact:**
- **Does NOT affect classification** - this happens after prediction
- **DOES affect output** - reported z-scores are misleading
- **DOES affect interpretation** - users may think introns were classified using these z-scores
- **DOES affect reproducibility** - z-scores depend on the entire dataset

**Why This Was Done:**
The comment suggests the author wanted z-scores "based on entire data" for output. This might make z-scores more interpretable for a specific dataset, but it violates ML best practices.

**Severity:** Medium-High ⚠️
- Doesn't affect published results (classification already done)
- But violates best practices and could confuse users
- Should be removed in refactor

**Recommendation:**
- **Remove** this re-normalization step entirely
- Output z-scores should be the **same** z-scores used for classification
- If users want dataset-specific statistics, provide them separately (not as replacement z-scores)

---

### ⚠️ ISSUE 2: Redundant Normalization in Recursive Scoring

**Location:** `intronIC.py:3952-3980`

**Code:**
```python
scale_vector = get_score_vector(
    raw_refs, score_names=raw_score_names)
recursive_scaler = preprocessing.RobustScaler().fit(scale_vector)  # First normalization
scored_introns = scale_scores(raw_introns, recursive_scaler)

# ... (classification happens here)

scale_vector = get_score_vector(
    recursive_refs, score_names=raw_score_names)
recursive_scaler = preprocessing.RobustScaler().fit(scale_vector)  # Second normalization!
scored_introns = scale_scores(raw_introns, recursive_scaler)
scored_refs = scale_scores(recursive_refs, recursive_scaler)
```

**Problem:**
1. Normalizes introns twice with **different** scalers
2. Second scaler is fit on a **filtered** subset of introns (high-confidence U12s + low-confidence U2s)
3. This filtering introduces selection bias into normalization parameters
4. The final scores use a scaler trained on post-classification filtered data

**Impact:**
- Creates circular dependency: classification → filtering → re-normalization → re-classification
- Second normalization is based on biased sample (excludes mid-range scores)
- Unclear why this is needed - seems like a bug or leftover from development

**Why This Was Done:**
Unclear. May have been intended to improve training set for second round, but the approach is flawed.

**Severity:** Medium ⚠️
- Affects recursive mode only (optional feature)
- Creates subtle bias in second-round training
- Should be cleaned up in refactor

**Recommendation:**
- **Simplify:** Normalize once, train once per round
- If recursive training is needed:
  1. Build new PWMs from high-confidence introns
  2. Re-score all introns with new PWMs
  3. Normalize using **all** reference sequences (not filtered subset)
  4. Re-train SVM
  5. Re-classify

---

### ⚠️ ISSUE 3: Inconsistent Scaler Choice (StandardScaler vs RobustScaler)

**Location:**
- Standard mode: `preprocessing.StandardScaler()` (line 3731)
- Recursive mode: `preprocessing.RobustScaler()` (line 3954, 3978)

**Problem:**
- `StandardScaler` uses mean and standard deviation
- `RobustScaler` uses median and IQR (robust to outliers)
- Switching scalers between rounds changes feature distributions
- Models trained on StandardScaler features won't work optimally on RobustScaler features

**Impact:**
- Inconsistent normalization across pipeline modes
- Could affect classification performance in recursive mode
- No clear justification for the switch

**Severity:** Low-Medium ⚠️
- Only affects recursive mode
- Both scalers are reasonable choices
- But inconsistency is problematic

**Recommendation:**
- **Use StandardScaler consistently** throughout pipeline
- If outliers are a concern, handle them explicitly (e.g., winsorization)
- Document the choice and rationale

---

### 💡 MINOR: Inefficient Iteration Structure in Optimization

**Location:** `intronIC.py:5465-5475`

**Code:**
```python
for search_round in range(1, n_optimize + 1):
    round_seed = seed + search_round
    search_model, performance = train_svm(
        base_model,
        u2_vector,
        u12_vector,
        iterations=iterations,  # This trains 10+ models per round!
        cv_parameters=cv_params,
        seed=round_seed)
```

**Observation:**
- Each optimization round trains `iterations` models (default 10)
- With 5 optimization rounds, this trains **50+ models** total
- Only the last round's models are kept
- Earlier rounds are just for hyperparameter search

**Impact:**
- Computationally expensive (trains many models that are discarded)
- But necessary for robust hyperparameter optimization
- Could potentially be optimized with early stopping or fewer iterations per round

**Severity:** Low 💡 (optimization opportunity, not a bug)

**Recommendation:**
- Consider **reducing iterations** during early optimization rounds (e.g., 3 models for rounds 1-4, 10 for round 5)
- Or use **smaller U2 subsamples** during optimization (commented code at line 5456 suggests this was considered)
- Trade-off between optimization quality and runtime

---

## Detailed Analysis by Component

### Component 1: PWM Scoring

**Location:** `intronIC.py:2114-2200`

**Function:** `seq_score()`

**Approach:**
- Log-odds ratio scoring: U12 frequency / U2 frequency
- Pseudocount to avoid zeros
- Option to ignore specific positions (e.g., canonical dinucleotides)

**Assessment:** ✅ Standard and appropriate
- Log-odds is correct statistical approach
- Pseudocount is standard practice
- Implementation is clean

**Recommendation:** Keep approach, improve modularity in refactor

---

### Component 2: Branch Point Detection

**Location:** `intronIC.py:2143-2200`

**Function:** `bp_score()`

**Approach:**
- Sliding window across branch point region (-55 to -5 from 3' ss)
- Score each window with PWM
- Return max score and position

**Assessment:** ✅ Simple and effective
- Exhaustive search guarantees finding best match
- Computationally efficient for short regions

**Limitations:**
- Assumes best score = true branch point
- No probabilistic model of BP positioning
- Could potentially use HMM or profile HMM for more sophisticated modeling

**Recommendation:** Keep approach for refactor, but consider advanced methods for future work

---

### Component 3: Z-Score Normalization

**Location:** `intronIC.py:3655-3665, 3725-3732`

**Functions:** `scale_scores()`, normalization in `apply_scores()`

**Approach:**
- Fit `StandardScaler` on reference sequences only
- Transform both reference and experimental sequences

**Assessment:** ✅ Correct (except for post-classification re-normalization)
- Proper train/test separation
- Standard approach in ML pipelines

**Issue:** The post-classification re-normalization (Issue #1 above)

**Recommendation:** Remove post-classification re-normalization

---

### Component 4: SVM Hyperparameter Optimization

**Location:** `intronIC.py:5431-5528`

**Function:** `optimize_svm()`

**Approach:**
1. Start with coarse grid: C from 10^-6 to 10^6 (13 values)
2. For 5 rounds:
   - Train models with GridSearchCV
   - Find all rank-1 C values
   - Take geometric mean
   - Refine grid around best value
3. Final C is geometric mean of refined grid

**Assessment:** ✅ Sophisticated and well-designed
- Geometric scale is appropriate for C parameter
- Refinement strategy efficiently narrows search
- Averaging rank-1 parameters is conservative and thoughtful

**Alternatives:**
- Random search (less efficient for single parameter)
- Bayesian optimization (overkill for 1D search, adds complexity)
- Halving grid search (similar efficiency)

**Recommendation:** Keep approach, it's excellent

---

### Component 5: SVM Training

**Location:** `intronIC.py:5345-5428`

**Function:** `train_svm()`

**Approach:**
- Train `iterations` models (default 10)
- Each model uses different U2 subsample (size = U12 count)
- 80/20 stratified split per model
- Evaluate on test set (F1, PR-AUC)

**Assessment:** ✅ Solid ensemble approach
- Handles class imbalance well
- Multiple models reduce variance
- Stratification ensures balanced splits

**Observation:**
- Uses `np.random.choice(..., replace=False)` for subsampling
- With enough data, this ensures diversity across models
- Could potentially use bootstrap (with replacement) for more diversity

**Recommendation:** Keep approach

---

### Component 6: Ensemble Prediction

**Location:** `intronIC.py:5651-5687, 5816-5853`

**Functions:** `average_svm_score_info()`, `parallel_svm_score()`

**Approach:**
- Apply all models to each intron
- Get U12 probabilities from each model
- Weight by F1 scores (optional)
- Average probabilities
- Classify based on threshold (default 90%)

**Assessment:** ✅ Standard ensemble approach
- Soft voting (average probabilities) is appropriate
- F1-weighting is sensible
- High threshold (90%) reduces false positives

**Alternative:**
- Hard voting (majority vote) - less nuanced
- Stacking (meta-model) - overkill for this problem

**Recommendation:** Keep approach

---

### Component 7: Recursive Scoring (Optional)

**Location:** `intronIC.py:3890-4050`

**Function:** `recursive_scoring()`

**Approach:**
1. Use first-round high-confidence U12s to build new PWMs
2. Re-score all introns with new PWMs
3. Filter to high-confidence introns
4. Re-train SVM on species-specific data
5. Re-classify all introns

**Assessment:** ⚠️ Good idea, flawed implementation
- Concept is sound: adapt to species-specific patterns
- Execution has issues (redundant normalization, filtering bias)

**Issues:**
- Double normalization (Issue #2)
- Uses RobustScaler instead of StandardScaler (Issue #3)
- Filtering creates biased training set for second round

**Recommendation:**
- Simplify implementation
- Fix normalization issues
- Consider whether recursive training is worth the complexity (may not improve results much)

---

## Recommendations for Refactor

### Priority 1: Fix Critical Issues

1. **Remove post-classification re-normalization** (Issue #1)
   - Delete lines 5247-5251 equivalent in refactor
   - Output z-scores should be same as classification z-scores
   - Add tests to ensure z-scores are consistent

2. **Fix recursive scoring normalization** (Issue #2)
   - Normalize once per round, not twice
   - Don't filter training data before normalization
   - Use same scaler type throughout (StandardScaler)

3. **Make scaler choice consistent** (Issue #3)
   - Use StandardScaler throughout
   - Document why StandardScaler is chosen
   - If outliers are concern, handle explicitly

### Priority 2: Improve Code Quality

4. **Modularize components**
   - Separate PWM scoring, normalization, training, prediction
   - Each component should be independently testable
   - Clear interfaces between components

5. **Add comprehensive tests**
   - Test that normalization only uses training data
   - Test that z-scores are consistent with classification
   - Test ensemble prediction logic
   - Test hyperparameter optimization convergence

6. **Document ML rationale**
   - Why ensemble approach?
   - Why geometric mean for C?
   - Why 90% threshold?
   - Add these explanations to docstrings

### Priority 3: Performance Optimizations

7. **Optimize hyperparameter search** (optional)
   - Reduce iterations in early optimization rounds
   - Use smaller U2 subsamples during optimization
   - Consider early stopping

8. **Parallelize more components**
   - PWM scoring is already parallelized
   - Consider parallel training of ensemble models
   - Batch predictions more efficiently

### Priority 4: Advanced Features (Future)

9. **Confidence intervals on predictions**
   - Use ensemble variance to estimate uncertainty
   - Helps identify borderline cases

10. **Advanced branch point modeling**
    - Profile HMM for branch point detection
    - Positional probability model
    - May improve classification accuracy

---

## Test Strategy Recommendations

### Unit Tests

**Normalization:**
```python
def test_normalization_uses_training_only():
    """Z-score scaler should ONLY see reference sequences."""
    # Create reference and experimental introns
    # Score both with PWMs
    # Fit normalizer on references only
    # Transform both
    # Assert: experimental data never seen by scaler

def test_no_post_classification_renormalization():
    """Z-scores should not change after classification."""
    # Score and classify introns
    # Save z-scores
    # Write outputs
    # Assert: z-scores unchanged
```

**Train/Test Split:**
```python
def test_stratified_split_maintains_class_ratio():
    """Train and test sets should have same U2/U12 ratio."""
    # Create balanced dataset
    # Split with stratification
    # Assert: ratios are equal (within tolerance)

def test_no_data_leakage_across_splits():
    """Same intron should not appear in both train and test."""
    # Split data
    # Assert: no overlap between train and test sets
```

**Ensemble:**
```python
def test_ensemble_uses_different_u2_subsamples():
    """Each model should train on different U2 introns."""
    # Train ensemble
    # Extract U2 training data for each model
    # Assert: subsamples are different

def test_f1_weighted_averaging():
    """Better models should have more influence."""
    # Create models with different F1 scores
    # Make predictions
    # Assert: higher F1 model has more weight
```

**Hyperparameter Optimization:**
```python
def test_optimization_converges():
    """C parameter range should narrow each round."""
    # Run optimization with known data
    # Assert: range decreases each round

def test_rank1_param_averaging():
    """Geometric mean of rank-1 params should be used."""
    # Mock GridSearchCV results with multiple rank-1 values
    # Assert: uses geometric mean, not first value
```

### Integration Tests

```python
def test_full_pipeline_on_chr19():
    """Complete scoring and classification on chr19."""
    # Load chr19 introns
    # Run full pipeline
    # Assert: known U12s are classified correctly (>95% recall)

def test_output_zscore_consistency():
    """Z-scores in output files should match classification z-scores."""
    # Run pipeline
    # Parse output files
    # Assert: z-scores match internal scores used for classification

def test_recursive_scoring_improves_results():
    """Recursive training should find more U12s (if implemented)."""
    # Run standard pipeline
    # Run recursive pipeline
    # Assert: recursive finds same or more U12s
```

### Gold Standard Tests

```python
def test_matches_original_intronIC():
    """Classification should match original intronIC (after fixes)."""
    # Run original intronIC on chr19
    # Run refactored intronIC on chr19
    # Assert: classifications match (after accounting for fixed issues)
```

---

## Summary Table: Issues & Severity

| Issue | Location | Severity | Affects Classification? | Fix Priority |
|-------|----------|----------|------------------------|--------------|
| Post-classification re-normalization | 5247-5251 | Medium-High | No (only output) | **High** |
| Redundant normalization (recursive) | 3952-3980 | Medium | Yes (recursive mode) | **High** |
| Inconsistent scaler type | 3731, 3954 | Low-Medium | Yes (recursive mode) | **Medium** |
| Inefficient optimization iterations | 5465-5475 | Low | No (performance only) | **Low** |

---

## Overall Assessment

**Grade: B+**

**Strengths:**
- Proper train/test separation for normalization (before classification)
- Sophisticated hyperparameter optimization
- Thoughtful ensemble approach for class imbalance
- Appropriate evaluation metrics
- Good use of cross-validation

**Weaknesses:**
- Post-classification re-normalization violates ML best practices
- Recursive scoring has implementation issues
- Inconsistent scaler choice
- Could be more computationally efficient

**Bottom Line:**
The core ML approach is **sound and appropriate** for this problem. The issues identified are **fixable** and mostly affect output/interpretation rather than classification accuracy. The refactored version should address these issues while preserving the strong foundational approach.

---

## References

**ML Best Practices:**
- Proper train/test separation: [Scikit-learn docs](https://scikit-learn.org/stable/modules/cross_validation.html)
- Imbalanced classification: [He & Garcia (2009)](https://doi.org/10.1109/TKDE.2008.239)
- Ensemble methods: [Dietterich (2000)](https://doi.org/10.1007/3-540-45014-9_1)

**Original intronIC Paper:**
- Moyer et al. (2020): "Comprehensive database and evolutionary dynamics of U12-type introns"
- *Nucleic Acids Research*, 48(13):7066–7078
- https://doi.org/10.1093/nar/gkaa464

---

**Next Steps:**
1. Review this analysis with project author
2. Implement fixes in refactored code
3. Add comprehensive tests for ML pipeline
4. Validate that fixes don't change classification results (except in recursive mode)
5. Document ML approach clearly in refactored code
