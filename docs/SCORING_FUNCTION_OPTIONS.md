# Scoring Function Options for Precision-Focused Training

## Context

Expert recommendation: "Train with a scoring function that prioritizes precision"

**Why?** We want the model to be conservative—only call U12 when very sure. This reduces false positives in cross-species applications.

## Current Situation

**Current scoring:** `'balanced_accuracy'`
```python
# From classification/optimizer.py
grid_search = GridSearchCV(
    ...,
    scoring='balanced_accuracy',  # Equal weight to both classes
    ...
)
```

**Problem:** Balanced accuracy treats precision and recall equally. We want to favor precision.

---

## Option 1: F_β Score (β < 1)

**Formula:**
```
F_β = (1 + β²) × (precision × recall) / (β² × precision + recall)
```

**β values:**
- β = 1.0: Standard F1 (equal weight to precision and recall)
- β = 0.5: Weights precision 2× more than recall (F_0.5)
- β = 0.25: Weights precision 4× more than recall (F_0.25)

**Implementation:**
```python
from sklearn.metrics import fbeta_score, make_scorer

# F_0.5: Precision weighted 2x more than recall
scorer = make_scorer(fbeta_score, beta=0.5, pos_label='u12')

grid_search = GridSearchCV(
    ...,
    scoring=scorer,
    ...
)
```

**Pros:**
- Standard sklearn metric
- Single hyperparameter (β)
- Interpretable: β=0.5 means "I care 2x more about precision"
- Smooth tradeoff curve

**Cons:**
- Still allows low precision if recall is very high
- Doesn't guarantee minimum recall threshold

**Recommendation:** β = 0.5 (2× weight to precision)

---

## Option 2: Precision at Fixed Recall

**Concept:** Find the decision threshold that achieves ≥X% recall, then report precision at that threshold.

**Example:** "Give me ≥90% recall, maximize precision subject to that constraint"

**Implementation:**
```python
from sklearn.metrics import precision_recall_curve

def precision_at_recall(y_true, y_score, min_recall=0.90):
    """
    Find highest precision while maintaining min_recall.

    Args:
        y_true: True labels
        y_score: Predicted probabilities
        min_recall: Minimum acceptable recall (e.g., 0.90 = 90%)

    Returns:
        Precision at the threshold that gives min_recall
    """
    precision, recall, thresholds = precision_recall_curve(y_true, y_score)

    # Find thresholds where recall >= min_recall
    valid_indices = recall >= min_recall

    if not any(valid_indices):
        return 0.0  # Can't achieve min_recall

    # Return highest precision among valid thresholds
    return max(precision[valid_indices])

# Make it a scorer
from sklearn.metrics import make_scorer
scorer = make_scorer(
    precision_at_recall,
    min_recall=0.90,
    needs_proba=True
)
```

**Pros:**
- Guarantees minimum recall (won't miss too many real U12s)
- Directly optimizes precision (our goal)
- Very interpretable: "90% recall, best possible precision"

**Cons:**
- Requires calibrated probabilities
- Custom implementation (not built-in sklearn)
- Less smooth optimization surface
- Need to choose min_recall threshold

**Recommendation:** min_recall = 0.90 (catch 90% of U12s, maximize precision)

---

## Option 3: Precision-Recall AUC

**Concept:** Area under precision-recall curve (like ROC AUC, but for precision vs recall)

**Implementation:**
```python
from sklearn.metrics import average_precision_score

# This is actually "average_precision", a.k.a. PR-AUC
grid_search = GridSearchCV(
    ...,
    scoring='average_precision',  # Built-in sklearn
    ...
)
```

**Pros:**
- Built-in sklearn metric
- Summarizes performance across all thresholds
- Good for imbalanced classification
- Smooth metric for optimization

**Cons:**
- Less interpretable than F_β
- Doesn't explicitly favor precision over recall
- Treats all thresholds equally (we care about one operating point)

**Recommendation:** Use as secondary metric, not primary

---

## Option 4: Matthews Correlation Coefficient (MCC)

**Formula:**
```
MCC = (TP × TN - FP × FN) / sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN))
```

**Range:** -1 (worst) to +1 (perfect)

**Implementation:**
```python
from sklearn.metrics import matthews_corrcoef, make_scorer

scorer = make_scorer(matthews_corrcoef)

grid_search = GridSearchCV(
    ...,
    scoring=scorer,
    ...
)
```

**Pros:**
- Works well with imbalanced datasets
- Single number summary
- Symmetric (treats both classes fairly)

**Cons:**
- Symmetric (doesn't favor precision)
- Less interpretable than F_β
- Requires hard predictions (not probabilities)

**Recommendation:** Not ideal for our precision-focused goal

---

## Option 5: Custom Weighted Score

**Concept:** Manually combine precision and recall with custom weights

**Implementation:**
```python
def custom_weighted_score(y_true, y_pred):
    """
    70% precision, 30% recall
    """
    from sklearn.metrics import precision_score, recall_score

    prec = precision_score(y_true, y_pred, pos_label='u12', zero_division=0)
    rec = recall_score(y_true, y_pred, pos_label='u12', zero_division=0)

    return 0.7 * prec + 0.3 * rec

scorer = make_scorer(custom_weighted_score)
```

**Pros:**
- Complete control over tradeoff
- Interpretable weights (70% precision, 30% recall)
- Can adjust based on domain requirements

**Cons:**
- Not standard (hard to compare to literature)
- Weights are arbitrary
- F_β already does this (more principled)

**Recommendation:** Only if F_β doesn't meet needs

---

## Comparison Table

| Metric | Precision Focus | Interpretability | Built-in sklearn | Needs Probabilities |
|--------|-----------------|------------------|------------------|---------------------|
| **F_0.5 score** | ⭐⭐⭐ (2× weight) | ⭐⭐⭐⭐ | ✅ Yes | ❌ No |
| **Precision @ 90% recall** | ⭐⭐⭐⭐⭐ (direct) | ⭐⭐⭐⭐⭐ | ❌ Custom | ✅ Yes |
| **PR-AUC** | ⭐⭐ (indirect) | ⭐⭐ | ✅ Yes | ✅ Yes |
| **MCC** | ⭐ (symmetric) | ⭐⭐ | ✅ Yes | ❌ No |
| **Custom 70/30** | ⭐⭐⭐ (tunable) | ⭐⭐⭐ | ❌ Custom | ❌ No |

---

## Recommendation by Use Case

### If you want simplicity and standard metric:
→ **F_0.5 score** (β=0.5, precision weighted 2×)
```python
from sklearn.metrics import make_scorer, fbeta_score
scorer = make_scorer(fbeta_score, beta=0.5)
```

### If you have a hard recall requirement:
→ **Precision @ 90% recall** (ensure we catch most U12s)
```python
# Custom implementation (see Option 2)
scorer = make_scorer(precision_at_recall, min_recall=0.90, needs_proba=True)
```

### If you want extreme precision bias:
→ **F_0.25 score** (β=0.25, precision weighted 4×)
```python
scorer = make_scorer(fbeta_score, beta=0.25)
```

---

## Current Human Model Performance

**For reference,** here's what we typically see on human test set:

```
Balanced Accuracy: ~0.95
Precision: ~0.85-0.90  (85-90% of U12 calls are correct)
Recall: ~0.90-0.95     (90-95% of real U12s are caught)
F1 score: ~0.88
```

**Question:** What tradeoff do you prefer?

**Option A (F_0.5):** Precision ~92%, Recall ~88%
- Catch slightly fewer U12s, but higher confidence in calls
- Balanced tradeoff

**Option B (Precision @ 90% recall):** Precision ~87%, Recall ≥90%
- Guarantee we catch 90% of U12s
- Maximize precision subject to that constraint

**Option C (F_0.25):** Precision ~95%, Recall ~82%
- Very conservative (few false positives)
- Risk missing ~18% of real U12s

---

## What I Need to Know

To proceed with implementation, please tell me:

1. **Which metric?** (F_0.5, Precision@Recall, or other)

2. **If F_β:** What β value? (0.5 = 2× precision weight, 0.25 = 4× weight)

3. **If Precision@Recall:** What minimum recall? (0.90 = 90%, 0.85 = 85%)

4. **Trade-off preference:**
   - More important to catch all real U12s? (higher recall, accept some FPs)
   - More important to avoid false positives? (higher precision, miss some U12s)

**My recommendation:** Start with **F_0.5** (standard, interpretable, balanced precision focus). Can adjust β later if needed.
