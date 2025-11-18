# Scoring Metric Comparison: neg_log_loss vs F_β

## What Gets Optimized During Grid Search

When you run GridSearchCV, the scoring metric determines **which hyperparameters get selected** from the grid.

---

## Current: `scoring='neg_log_loss'`

### What it optimizes
**Cross-entropy loss** (also called log loss):
```
neg_log_loss = -1/N * Σ [y_i * log(p_i) + (1-y_i) * log(1-p_i)]
```

Where `p_i` is the predicted probability for sample i.

### Behavior
- **Rewards:** Well-calibrated probabilities
- **Penalizes:** Confident wrong predictions (heavily!)
- **Example:**
  - True label: U12 (y=1)
  - Predict p=0.99 (confident, correct): loss = -log(0.99) = 0.01
  - Predict p=0.51 (barely correct): loss = -log(0.51) = 0.67
  - Predict p=0.01 (confident, wrong): loss = -log(0.01) = 4.6 ← Very bad!

### What hyperparameters it favors
- Hyperparameters that produce **calibrated probabilities**
- Models that are **confident when correct, uncertain when unsure**
- Works well with `CalibratedClassifierCV` (which you use)
- **Doesn't explicitly favor precision or recall**

### Example Selection
Imagine two hyperparameter combinations:

**Model A:**
- Precision: 95%, Recall: 85%
- Probabilities: Well-calibrated (p≈0.95 → ~95% correct)
- neg_log_loss: **0.12**

**Model B:**
- Precision: 90%, Recall: 90%
- Probabilities: Poorly calibrated (p≈0.95 → only 75% correct)
- neg_log_loss: **0.28**

**GridSearchCV picks:** Model A (better calibration, lower loss)

---

## Alternative: `scoring=F_0.5` (Precision-Focused)

### What it optimizes
**F_β score** with β=0.5 (weights precision 2× more than recall):
```
F_0.5 = 1.25 * (precision * recall) / (0.25 * precision + recall)
```

### Behavior
- **Rewards:** High precision (fewer false positives)
- **Accepts:** Slightly lower recall (missing some true U12s)
- **Ignores:** Probability calibration
- **Example:**
  - Precision=95%, Recall=85% → F_0.5 = 0.924
  - Precision=85%, Recall=95% → F_0.5 = 0.864
  - First one wins (higher precision)

### What hyperparameters it favors
- Hyperparameters that **minimize false positives**
- More conservative decision boundaries
- Potentially:
  - Tighter clipping (smaller zmax)
  - Stronger regularization (larger C in L2 penalty sense... wait, LinearSVC uses inverse)
  - Actually: **Smaller C** (more regularization, simpler model, conservative)
  - `saturate=True` (compresses extreme values)

### Example Selection
Same two models:

**Model A:**
- Precision: 95%, Recall: 85%
- F_0.5: **0.924**

**Model B:**
- Precision: 90%, Recall: 90%
- F_0.5: **0.900**

**GridSearchCV picks:** Model A (higher precision, higher F_0.5)

---

## Side-by-Side Comparison

| Aspect | neg_log_loss | F_0.5 |
|--------|--------------|-------|
| **Optimizes for** | Calibrated probabilities | High precision |
| **Cares about** | Confidence levels | Binary classification |
| **False positive bias** | Neutral | Penalizes FPs 2× more |
| **Hyperparameter selection** | Best calibration | Most conservative |
| **Works with** | Probability outputs | Hard predictions |
| **Typical result** | Balanced P/R, good calibration | Higher P, lower R |

---

## Concrete Example: Grid Search on Human Data

Imagine grid searching over:
```yaml
clip__quantile: [0.95, 0.975, 0.99]
saturate__enabled: [true, false]
C: [0.001, 0.01, 0.1]
```

### With `neg_log_loss`

**Cross-validation results:**

| quantile | saturate | C | Precision | Recall | Log Loss | F_0.5 |
|----------|----------|---|-----------|--------|----------|-------|
| 0.95 | true | 0.001 | 0.94 | 0.83 | 0.15 | 0.916 |
| 0.975 | false | 0.01 | 0.88 | 0.92 | **0.11** ← Best | 0.882 |
| 0.99 | false | 0.1 | 0.85 | 0.95 | 0.18 | 0.858 |

**Selected:** `quantile=0.975, saturate=false, C=0.01`
- Best calibrated probabilities
- Balanced precision/recall (88%/92%)

### With `F_0.5`

**Cross-validation results:**

| quantile | saturate | C | Precision | Recall | Log Loss | F_0.5 |
|----------|----------|---|-----------|--------|----------|-------|
| 0.95 | true | 0.001 | 0.94 | 0.83 | 0.15 | **0.916** ← Best |
| 0.975 | false | 0.01 | 0.88 | 0.92 | 0.11 | 0.882 |
| 0.99 | false | 0.1 | 0.85 | 0.95 | 0.18 | 0.858 |

**Selected:** `quantile=0.95, saturate=true, C=0.001`
- Highest precision (94%)
- More conservative (tighter clipping, saturation enabled)

---

## Impact on C. elegans False Positives

### Current (neg_log_loss)
1. Selects hyperparameters for **best calibration**
2. May choose looser clipping (quantile=0.975 or 0.99)
3. May skip saturation (if it hurts calibration)
4. **On C. elegans:** Allows z=11.45 to pass through → FP

### With F_0.5
1. Selects hyperparameters for **fewest false positives**
2. Likely chooses tighter clipping (quantile=0.95)
3. Likely enables saturation (further FP protection)
4. **On C. elegans:** z=11.45 gets clipped/saturated → fewer FPs

---

## The Real Question

**Does changing the scoring metric actually matter?**

The answer depends on whether the architectural changes (clipping + saturation) are **sufficient on their own**, or if we also need **different hyperparameter selection** to get the FP reduction.

### Scenario 1: Architecture is enough
- The new pipeline (clip + saturate) works well at **any** reasonable hyperparameter setting
- Both scoring metrics would pick good hyperparameters
- Stick with `neg_log_loss` (maintains calibration, simpler)

### Scenario 2: Hyperparameter selection matters
- Some hyperparameters still allow FPs (e.g., loose clipping)
- Need grid search to **actively favor conservative settings**
- Switch to `F_0.5` (guides grid search toward FP reduction)

---

## My Hypothesis

**I think Scenario 2 is more likely.** Here's why:

The grid search needs to **choose** between:
- `clip__quantile=0.95` (aggressive, clips z>~3.2)
- `clip__quantile=0.99` (loose, clips z>~5.0)

**With neg_log_loss:**
- If looser clipping gives better calibration on human data, it wins
- But looser clipping allows more extreme cross-species values through

**With F_0.5:**
- Tighter clipping reduces FPs on human data → higher F_0.5
- Grid search **learns** that aggressive clipping is better
- This generalizes to cross-species (fewer FPs in C. elegans)

---

## Recommendation

**Switch to F_0.5 (or F_0.25 for even more precision focus).**

**Reasoning:**
1. Expert specifically recommended precision-focused training
2. Grid search will learn to favor conservative hyperparameters
3. This is the **mechanism** by which the architectural changes get properly tuned
4. You can still get calibrated probabilities (that's what `CalibratedClassifierCV` does)

**Keep `class_weight='balanced'`** - that handles the imbalance in the loss function, orthogonal to the grid search metric.

**Alternative:** Try both!
1. Train with `neg_log_loss` (current)
2. Train with `F_0.5` (new)
3. Compare on C. elegans test set
4. See which gives fewer FPs

---

## Code Change

**Minimal change to `classification/optimizer.py` line 550:**

```python
# BEFORE
grid_search = GridSearchCV(
    model,
    param_grid=param_grid,
    cv=cv_splitter,
    scoring='neg_log_loss',  # ← Current
    n_jobs=self.n_jobs,
    ...
)

# AFTER (Option 1: F_0.5)
from sklearn.metrics import make_scorer, fbeta_score
scorer = make_scorer(fbeta_score, beta=0.5)

grid_search = GridSearchCV(
    model,
    param_grid=param_grid,
    cv=cv_splitter,
    scoring=scorer,  # ← F_0.5 score
    n_jobs=self.n_jobs,
    ...
)

# AFTER (Option 2: Keep current, make configurable)
# Add to config file
scoring_metric = config.get('scoring', 'neg_log_loss')
if scoring_metric == 'f_beta':
    from sklearn.metrics import make_scorer, fbeta_score
    beta = config.get('beta', 0.5)
    scorer = make_scorer(fbeta_score, beta=beta)
else:
    scorer = scoring_metric  # Use string (sklearn built-in)
```

---

## Summary

**neg_log_loss:**
- Optimizes probability calibration
- Neutral on precision/recall tradeoff
- May select looser hyperparameters if they calibrate better

**F_0.5:**
- Optimizes precision (2× weight vs recall)
- Guides grid search toward conservative hyperparameters
- Should reduce C. elegans FPs by selecting tighter clipping/saturation

**My vote:** **F_0.5** (aligns with expert recommendation, guides hyperparameter selection toward FP reduction)
