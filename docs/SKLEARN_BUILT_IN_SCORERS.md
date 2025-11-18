# Built-in sklearn Scorers for intronIC

## Relevant Built-in Options

### 1. `'neg_log_loss'` (Current)
**What it is:** Negative cross-entropy loss (negative because GridSearchCV maximizes)

**Use case:** Optimizing probability calibration

**Pros:**
- Built-in string
- Good for calibrated probabilities
- Works well with CalibratedClassifierCV

**Cons:**
- Doesn't favor precision over recall
- May select hyperparameters that allow FPs if they calibrate better

**Current performance:** Allows C. elegans FPs

---

### 2. `'balanced_accuracy'`
**What it is:** Average of recall for each class
```
balanced_acc = (recall_u12 + recall_u2) / 2
```

**Use case:** Binary classification with class imbalance

**Pros:**
- Built-in string
- Handles imbalance explicitly
- Symmetric (treats both classes equally)

**Cons:**
- Based on recall, not precision
- Doesn't favor precision over recall
- Might still allow FPs if it means better U2 recall

**Verdict:** Neutral on precision/recall tradeoff (like neg_log_loss)

---

### 3. `'average_precision'`
**What it is:** Area under Precision-Recall curve (PR-AUC)

**Formula:** Average of precision values at each recall threshold

**Use case:** Imbalanced classification, summary across all thresholds

**Pros:**
- Built-in string
- Recommended for imbalanced datasets
- Considers entire precision-recall tradeoff curve
- Better than ROC-AUC for rare positive class

**Cons:**
- Doesn't explicitly favor precision over recall
- Treats all operating points equally (we care about one decision threshold)
- Summary metric (may pick different hyperparameters than F_0.5)

**Interesting property:** Correlated with F1, but considers full PR curve

**Verdict:** Good for imbalanced data, but not precision-focused

---

### 4. `'precision'`
**What it is:** Just precision (TP / (TP + FP))

**Use case:** When you only care about precision

**Pros:**
- Built-in string
- Directly optimizes precision (our goal!)

**Cons:**
- **DANGER:** Ignores recall completely
- Can achieve 100% precision by predicting ZERO positives
- Would likely select very conservative hyperparameters (good?)
- Might miss too many real U12s (bad!)

**Verdict:** Too extreme (need to balance with recall)

---

### 5. `'f1'` (Standard F1)
**What it is:** Harmonic mean of precision and recall (β=1)
```
F1 = 2 * (precision * recall) / (precision + recall)
```

**Use case:** Balanced precision/recall

**Pros:**
- Built-in string
- Standard metric
- Balances precision and recall

**Cons:**
- Equal weight (doesn't favor precision)
- Same as current behavior (neutral on P/R tradeoff)

**Verdict:** Better than nothing, but doesn't address our precision-focus need

---

### 6. `make_scorer(fbeta_score, beta=0.5)` (F_0.5)
**What it is:** F-score with adjustable precision/recall weight
```
F_β = (1 + β²) * (precision * recall) / (β² * precision + recall)
```

**β=0.5:** Precision weighted 2× more than recall

**Use case:** When you want precision focus but still care about recall

**Pros:**
- Directly implements expert recommendation
- Explicit precision bias (tunable via β)
- Still considers recall (won't go to extremes)
- Standard metric (just needs make_scorer wrapper)

**Cons:**
- Not a built-in string (requires 2 lines of code)
- Need to choose β value

**Verdict:** Best match for "precision-focused training"

---

## Comparison for Our Use Case

| Scorer | Built-in String? | Precision Focus | Handles Imbalance | Recommended? |
|--------|------------------|-----------------|-------------------|--------------|
| `'neg_log_loss'` | ✅ Yes | ❌ No | Via calibration | ❌ Current, allows FPs |
| `'balanced_accuracy'` | ✅ Yes | ❌ No | ✅ Yes | ❌ Recall-based |
| `'average_precision'` | ✅ Yes | 🟡 Indirect | ✅ Yes | 🟡 Maybe |
| `'precision'` | ✅ Yes | ✅ Yes (too much!) | ✅ Yes | ❌ Too extreme |
| `'f1'` | ✅ Yes | ❌ No (balanced) | 🟡 Neutral | ❌ No precision focus |
| `F_0.5` (make_scorer) | ❌ No (2 lines) | ✅ Yes (2× weight) | ✅ Yes | ✅ **Best choice** |

---

## Testing 'average_precision'

Let me think through whether `'average_precision'` could work:

**What it optimizes:**
- Area under PR curve
- Essentially: "average precision across all decision thresholds"
- Good for imbalanced datasets

**Would it reduce C. elegans FPs?**
- Maybe? It optimizes for good precision-recall tradeoff across all thresholds
- But doesn't explicitly favor precision over recall at the chosen operating point
- Might select different hyperparameters than F_0.5

**Hypothesis:**
- F_0.5: "Pick hyperparameters that minimize FPs at our decision threshold"
- average_precision: "Pick hyperparameters that give good P/R tradeoff across all thresholds"

For our use case (reduce FPs at a specific threshold), **F_0.5 is more direct**.

---

## Code Comparison

### Option 1: average_precision (built-in string)
```python
grid_search = GridSearchCV(
    model,
    param_grid=param_grid,
    cv=cv_splitter,
    scoring='average_precision',  # ← One word change!
    n_jobs=self.n_jobs,
    ...
)
```

**Advantage:** Literally one word change
**Disadvantage:** Indirect precision focus

### Option 2: F_0.5 (make_scorer)
```python
from sklearn.metrics import make_scorer, fbeta_score

scorer = make_scorer(fbeta_score, beta=0.5)

grid_search = GridSearchCV(
    model,
    param_grid=param_grid,
    cv=cv_splitter,
    scoring=scorer,  # ← Two extra lines
    n_jobs=self.n_jobs,
    ...
)
```

**Advantage:** Direct precision focus (expert recommendation)
**Disadvantage:** Two extra lines (hardly "from scratch"!)

---

## Recommendation

**Use F_0.5 via make_scorer.**

**Why:**
1. Expert specifically said "precision-focused training" → F_β with β < 1 is the standard way
2. It's not "from scratch" - just a simple wrapper around built-in `fbeta_score`
3. `average_precision` is good but indirect (optimizes across all thresholds, not our operating point)
4. Two lines of code is negligible

**If you really want a built-in string:**
- Try `'average_precision'` first
- But I predict F_0.5 will work better for your specific use case

**My vote:** F_0.5 (2 lines, direct solution, expert-recommended)

---

## Alternative: Make it Configurable

```python
# In optimizer.__init__() or from config
scoring_metric = config.get('scoring', 'neg_log_loss')

if scoring_metric.startswith('f_beta_'):
    # e.g., 'f_beta_0.5'
    beta = float(scoring_metric.split('_')[-1])
    from sklearn.metrics import make_scorer, fbeta_score
    scorer = make_scorer(fbeta_score, beta=beta)
else:
    # Use built-in string
    scorer = scoring_metric

grid_search = GridSearchCV(..., scoring=scorer, ...)
```

This lets you experiment with:
- `'neg_log_loss'` (current)
- `'average_precision'` (built-in alternative)
- `'f_beta_0.5'` (precision-focused)
- `'f_beta_0.25'` (very precision-focused)

Without code changes, just config file.
