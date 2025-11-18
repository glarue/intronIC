# Scaling Architecture Review: Double Scaling Concern

## Background

intronIC classifies introns as U12-type (rare, ~0.5%) vs U2-type (common, ~99.5%) using:
1. Position-weight matrices (PWMs) → log-likelihood ratios (LLRs) for three regions (5'SS, BPS, 3'SS)
2. Feature normalization → z-scores
3. LinearSVC classifier with augmented features

## Current Scaling Architecture

**Two cascaded scaling steps:**

### Stage 1: ZeroAnchoredRobustScaler (Pre-processing)
**Location:** `scoring/normalizer.py`
**Purpose:** Convert raw LLRs to z-scores before classification

**Algorithm:**
```python
# Zero-anchored scaling (preserves semantic zero)
z = raw_score / median(|raw_scores|)
```

**Fitted on:** Reference introns (human U2/U12 curated set)
- Uses winsorization at 99.5th percentile
- Computes scale = median(|s|) per feature
- Does NOT center (preserves s=0 = "U12≈U2")

**Output:** Z-scores stored in intron objects (`five_z_score`, `bp_z_score`, `three_z_score`)

### Stage 2: RobustScaler (Model Pipeline)
**Location:** `classification/optimizer.py:507`, `classification/trainer.py:197`
**Purpose:** Scale features before SVM training/prediction

**Pipeline:**
```python
Pipeline([
    ('clip', OutlierClipper(quantile=0.999)),
    ('scale', RobustScaler(with_centering=False, with_scaling=True)),  # ← Second scaling
    ('augment', BothEndsStrongTransformer()),
    ('svc', LinearSVC(...))
])
```

**Algorithm:**
```python
# RobustScaler with_centering=False
z_scaled = z / IQR(z)
```

**Fitted on:** Same reference introns (during model training)
- Also fitted on z-scores from Stage 1
- Uses IQR (Q75 - Q25) instead of median(|s|)
- Also doesn't center (with_centering=False)

**Result:** Z-scores get scaled AGAIN before entering SVM

## Problem Observed

When applying human-trained model to *C. elegans* (nematode, 38% GC), we see extreme z-scores:

**Example intron:** `CaeEle-WBGene00044195@F49A5.10.1-intron_2`
- **3'SS sequence:** AACTACAGTG
- **Raw LLR:** 14.46 (P(seq|U12) / P(seq|U2) ≈ 2 million)
- **Z-score (Stage 1):** 11.45 (= 14.46 / 1.26)
- **SVM prediction:** 100% U12 (likely false positive)

**Why this happens:**
1. PWMs trained on vertebrates may not generalize to nematode composition bias
2. Stage 1 scaler fitted on human data (median(|3'SS|) ≈ 1.26)
3. C. elegans has different score distribution → extreme z-scores
4. These pass through Stage 2 scaler (quantile=0.999 is very permissive)
5. Linear model weights (~0.1-0.8) can't compensate for 11.45σ signal

**Model coefficients (Stage 2 output space):**
```
s5:                +0.84
sBP:               +0.13
s3:                +0.13
min_all:           +0.26
neg_absdiff_5_bp:  -0.016  (very weak!)
neg_absdiff_5_3:   -0.073
neg_absdiff_bp_3:  -0.090
```

**Decision function for this intron:**
- `s3 contribution: 11.45 × 0.13 = +1.49`
- `All penalties combined: -1.02`
- `Net: +0.47 → Predicts U12`

The extreme single-site score overwhelms all balance penalties.

## Questions for Expert

### Q1: Is double scaling statistically valid?
Should we be applying two sequential robust scalers?
- Stage 1: `z = s / median(|s|)`
- Stage 2: `z' = z / IQR(z)`

Or should we use **only one** scaling step?

### Q2: If double scaling is correct, which should be cross-species adapted?

**Current behavior:** Both scalers fitted on human reference data

**Option A:** Only refit Stage 1 for new species
```python
# For C. elegans scoring:
normalizer.fit(celegans_introns, dataset_type="unlabeled")  # Stage 1
# Keep Stage 2 as-is (fitted on human during training)
```

**Option B:** Refit both scalers for new species
```python
# Would require retraining entire model on C. elegans data
# (defeats purpose of pretrained model)
```

**Option C:** Remove Stage 1, use only Stage 2
```python
# Skip normalizer, pass raw LLRs to model pipeline
# Let RobustScaler handle everything
```

### Q3: Is with_centering=False correct for both scalers?

Both scalers use `with_centering=False` to preserve semantic zero (s=0 = "U12≈U2").

But after Stage 1, **z-scores already preserve zero** (zero-anchored by design). Does Stage 2 need `with_centering=False`? Or could it use standard scaling since the semantic zero is already encoded?

### Q4: Should OutlierClipper be more aggressive?

**Current:** `OutlierClipper(quantile=0.999)` - only clips extreme 0.1%
**Problem:** Allows z=11.45 to pass through unchecked

Should we:
- Make quantile tunable via grid search?
- Use tighter default (0.95, 0.975, 0.99)?
- Apply clipping to Stage 1 output instead of (or in addition to) Stage 2 input?

## Proposed Solutions

### Solution 1: Domain Adaptation (Current Capability)
Refit Stage 1 scaler on unlabeled target species data:

```python
# Already implemented in ScoreNormalizer
normalizer.fit(celegans_introns, dataset_type="unlabeled")
```

**Pros:**
- Corrects for species-specific composition bias
- Statistically valid (unsupervised, no label leakage)
- Robust estimators unaffected by rare U12s (~0.5%)

**Cons:**
- Requires unlabeled data from target species
- Adds workflow complexity

### Solution 2: Tighter Outlier Clipping
Add clipping quantile to hyperparameter grid:

```yaml
param_grid:
  estimator__clip__quantile: [0.95, 0.975, 0.99]
```

**Pros:**
- Simple, no workflow changes
- Prevents extreme outliers from dominating

**Cons:**
- Band-aid, doesn't fix root cause
- May clip legitimate strong signals

### Solution 3: Architectural Change
Remove Stage 1 scaling, let model pipeline handle everything:

**Pros:**
- Simpler architecture
- Only one scaler to reason about

**Cons:**
- Would require model retraining
- Loses explicit control over pre-SVM normalization

## Test Case Data

**Human reference 3'SS distribution (from model):**
- Scale factor: 1.26 (median of absolute values)
- Typical z-scores: -2 to +4 range

**C. elegans with human scaler:**
- Raw LLR: 14.46
- Z-score: 11.45 (extreme outlier)
- Prediction: 100% U12 (likely false positive)

**Expected behavior:**
- Balance constraint should reject imbalanced patterns
- min_all = -0.15 (negative!) should heavily penalize
- But coefficient (+0.26) × -0.15 = -0.04 (tiny penalty)
- Extreme s3 (+1.49) overwhelms all penalties

## Files Referenced

- `scoring/normalizer.py:36-158` - ZeroAnchoredRobustScaler (Stage 1)
- `scoring/normalizer.py:160-413` - ScoreNormalizer wrapper
- `classification/optimizer.py:505-509` - Pipeline with RobustScaler (Stage 2)
- `classification/clipping.py:17-174` - OutlierClipper
- `config/training_default.yaml` - Hyperparameter grid

## Questions Summary

1. **Is double scaling intended or accidental?**
2. **If intended, should cross-species applications refit Stage 1, Stage 2, or both?**
3. **Should Stage 2 use with_centering=False, or is that only needed for Stage 1?**
4. **Should OutlierClipper be tuned more aggressively to prevent extreme cross-species artifacts?**

---

**Request:** Please advise on whether this scaling architecture is statistically sound and, if so, which components should be adapted for cross-species deployment.
