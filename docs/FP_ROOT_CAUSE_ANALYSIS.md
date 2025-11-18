# False Positive Root Cause Analysis

**Date:** 2025-11-15
**Issue:** Model produces false positives on C. elegans despite FP reduction features
**Status:** 🔴 ROOT CAUSE IDENTIFIED

---

## Executive Summary

The pretrained homo_sapiens model is producing false positives because **L1 penalty (Lasso) zeroed out critical FP reduction features** due to mathematical redundancy between min features.

**Key findings:**
- ✅ FP reduction features ARE implemented correctly
- ✅ Model WAS trained with all new features
- ❌ L1 regularization ZEROED the most important features:
  - `min_all` (require all three signals strong) → **coefficient = 0**
  - All `neg_absdiff` features (imbalance penalties) → **all coefficients = 0**
  - Direct `sBP` and `s3` scores → **both coefficients = 0**
- ❌ L1 KEPT `max_5_3` which **rewards** one-end-strong patterns → **coefficient ≈ 0.4**

---

## Model Inspection Results

### Trained Model: `run_tests/hsapiens/homo_sapiens.model.pkl`

**Training date:** 2025-11-15 23:32:51
**Hyperparameters:**
- C = 0.0948
- Penalty = **L1** (Lasso - induces sparsity)
- include_max = **True**
- Calibration = isotonic

### Ensemble Coefficients (All 4 Models Show Same Pattern)

**Model 1 coefficients:**
```
s5                  : +1.249035  ← DOMINANT feature
sBP                 : +0.000000  ← ZEROED!
s3                  : +0.000000  ← ZEROED!
min_5_bp            : +0.527692
min_5_3             : +0.171555
min_all             : +0.000000  ← ZEROED! (Critical for FP rejection)
neg_absdiff_5_bp    : +0.000000  ← ZEROED!
neg_absdiff_5_3     : +0.000000  ← ZEROED!
neg_absdiff_bp_3    : +0.000000  ← ZEROED!
max_5_bp            : +0.000000
max_5_3             : +0.396311  ← BAD! Rewards one-end-strong
```

**Models 2-4:** Nearly identical pattern (weights vary ±10% but same features zeroed)

---

## Why This Causes False Positives

### Example: Top False Positive #5

**Intron:** `CaeEle-WBGene00044195@F49A5.10.1-intron_2(3)`
**Scores:** Relative=7.96, Probability=98.0%

**Z-scores:**
- s5 = **-0.154** (NEGATIVE!)
- sBP = 0.301 (weak)
- s3 = **11.451** (extremely high - classic one-end-strong)

**Feature contributions:**
```
s5:          1.249 * (-0.154) = -0.192
sBP:         0.000 * 0.301    =  0.000  ← ZEROED, no contribution
s3:          0.000 * 11.451   =  0.000  ← ZEROED, no contribution
min_5_bp:    0.528 * (-0.154) = -0.081
min_5_3:     0.172 * (-0.154) = -0.026
min_all:     0.000 * (-0.154) =  0.000  ← ZEROED! Should penalize
neg_abs_5bp: 0.000 * (-0.454) =  0.000  ← ZEROED! Should penalize
neg_abs_5_3: 0.000 * (-11.60) =  0.000  ← ZEROED! Should penalize
neg_abs_bp3: 0.000 * (-11.15) =  0.000  ← ZEROED! Should penalize
max_5_bp:    0.000 * 0.301    =  0.000
max_5_3:     0.396 * 11.451   = +4.535  ← HUGE POSITIVE! Rewards FP
Intercept:                      -0.199

Total score: -0.192 + 0 + 0 - 0.081 - 0.026 + 0 + 0 + 0 + 0 + 0 + 4.535 - 0.199 = +4.037
```

**Result:** High positive score despite:
- ❌ Negative s5 score
- ❌ Huge imbalance (s3 - s5 = 11.6 standard deviations!)
- ❌ min_all = -0.154 (should strongly reject)

**The problem:** `max_5_3` with coefficient +0.396 allows a single very high 3'SS score to override everything else.

---

## Pattern Analysis: All Top False Positives

| FP # | s5 | sBP | s3 | min_all | max_5_3 | Pattern |
|------|----|----|----|---------|---------|--------------------|
| #1 | 0.45 | 1.37 | **-1.90** | **-1.90** | 0.45 | 3'SS NEGATIVE |
| #2 | 0.28 | 1.12 | 0.41 | 0.28 | 0.41 | All weak, BP high |
| #3 | 0.39 | 0.10 | 0.74 | **0.10** | 0.74 | BP very weak |
| #4 | 0.36 | 0.30 | **-2.76** | **-2.76** | 0.36 | 3'SS VERY NEGATIVE |
| #5 | **-0.15** | 0.30 | **11.45** | **-0.15** | **11.45** | Only 3'SS strong |
| #6 | 0.32 | **-0.49** | 0.58 | **-0.49** | 0.58 | BP NEGATIVE |
| #7 | 0.36 | **-0.43** | **-2.66** | **-2.66** | 0.36 | BP and 3'SS NEG |
| #8 | 0.16 | 2.41 | 1.71 | **0.16** | 1.71 | 5'SS very weak |

**Common pattern:** 6 out of 8 FPs have **negative min_all** which should immediately reject them, but min_all coefficient is ZERO so it has no effect!

---

## Root Cause: Feature Redundancy + L1 Penalty

### Mathematical Redundancy

The min features are mathematically related:
```
min_all   = min(s5, sBP, s3)      ← Strongest constraint
min_5_bp  = min(s5, sBP)          ← Weaker constraint
min_5_3   = min(s5, s3)           ← Weaker constraint
```

**Relationship:**
- `min_all ≤ min_5_bp` (always)
- `min_all ≤ min_5_3` (always)
- `min_all` is the strictest requirement

### L1 Penalty Behavior

L1 (Lasso) regularization:
1. **Induces sparsity** - zeros out coefficients
2. **Picks one from correlated features** - eliminates redundancy
3. **Selection is based on training data** - not domain knowledge

**What happened:**
1. L1 saw min_5_bp, min_5_3, and min_all as redundant
2. L1 chose to keep min_5_bp + min_5_3 and ZERO min_all
3. This seemed reasonable on training data (all true U12s have strong signals everywhere)
4. But on real data with FPs, we need the **strictest** constraint (min_all)

### Why L1 Made Bad Choices

**Training data bias:**
- True U12s: s5 ≈ sBP ≈ s3 (all high, balanced)
- For balanced signals: min_all ≈ min_5_bp ≈ min_5_3 (all similar)
- L1 couldn't distinguish importance because they're equivalent on true positives
- Only on FALSE POSITIVES do we see the difference!

**L1 also kept max_5_3:**
- Positive coefficient on max encourages "at least one strong"
- This is EXACTLY the wrong behavior for FP rejection
- We want "ALL strong" not "at least one strong"

---

## Why This Wasn't Caught Earlier

1. **High cross-validation scores** (F1=0.9987, PR-AUC=1.0)
   - Training/validation data: curated true U12s and true U2s
   - No false positives in training data to reveal the weakness

2. **Model works on human genome**
   - Human annotations are high quality
   - Most spurious matches are already filtered by annotation quality
   - FP reduction features were "nice to have" not critical

3. **C. elegans exposed the weakness**
   - More distant from training data (nematode vs vertebrate)
   - Different splice site characteristics
   - More potential false positives with one-end-strong patterns

---

## Solutions

### Option 1: Force L2 Penalty (Recommended Short-term)

**Change:**
```python
# In config/training_default.yaml
param_grid:
  estimator__svc__penalty:
    - l2  # Remove l1 from options
```

**Effect:**
- L2 (Ridge) keeps ALL features with non-zero weights
- min_all will contribute
- neg_absdiff features will contribute
- No arbitrary feature selection

**Trade-off:**
- More features = slightly more complex model
- But only 11 features total, so negligible impact

### Option 2: Remove Redundant Features (Recommended Long-term)

**Simplify BothEndsStrong transformer:**
```python
# Keep only:
features = [
    s5, sBP, s3,           # Base features
    min_all,               # Strictest "all strong" requirement
    neg_absdiff_5_bp,      # Imbalance penalties
    neg_absdiff_5_3,
    neg_absdiff_bp_3
]
# Remove: min_5_bp, min_5_3, max_5_bp, max_5_3
```

**Benefits:**
- Removes mathematical redundancy
- Forces model to use the correct features
- Works with L1 or L2 penalty
- Cleaner, more interpretable

### Option 3: Constrain L1 Feature Selection

**Use custom feature penalties:**
```python
# In LinearSVC, penalize dropping critical features
# (Requires custom implementation or different algorithm)
```

**Effect:**
- Keep L1 sparsity benefits
- Prevent critical features from being zeroed

**Trade-off:**
- More complex implementation
- May not be worth the engineering effort

---

## Immediate Action Items

1. **Retrain homo_sapiens model with L2 penalty**
   ```bash
   # Update config to force L2
   intronIC train -n homo_sapiens -p 12 --force-l2
   ```

2. **Test on C. elegans with new model**
   ```bash
   intronIC classify -m homo_sapiens.model.pkl -g ... -a ... -n caenorhabditis_elegans
   ```

3. **Compare FP counts**
   - Current: 8 FPs predicted
   - Expected with L2: 0-2 FPs (min_all and neg_absdiff will reject them)

4. **Long-term: Refactor BothEndsStrong transformer**
   - Remove redundant min features
   - Keep only min_all
   - Consider removing max features entirely (expert: "you can drop max_* entirely")

---

## Key Takeaways

1. **Feature engineering was correct** - We identified the right features
2. **Implementation was correct** - Features are computed properly
3. **Regularization made bad choices** - L1 zeroed critical features due to redundancy
4. **Cross-validation didn't catch it** - No FPs in training data
5. **Solution is straightforward** - Switch to L2 or remove redundant features

**The irony:** We added min_all specifically to reject one-end-strong FPs, but L1 zeroed it and kept max_5_3 which REWARDS one-end-strong patterns!

---

**Analysis by:** Claude Code
**Date:** 2025-11-15
**Model inspected:** run_tests/hsapiens/homo_sapiens.model.pkl
