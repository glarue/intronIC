# Scaling Architecture: Before vs After

## BEFORE (Current - Double Scaling)

```
┌─────────────────────────────────────────────┐
│ STAGE 1: Pre-processing (outside pipeline) │
└─────────────────────────────────────────────┘
Raw LLRs: [s5, sBP, s3]
    ↓
ZeroAnchoredRobustScaler (fitted on human)
  z = s / median(|s|)
    ↓
Z-scores stored in Intron objects
    ↓
┌─────────────────────────────────────────────┐
│ STAGE 2: Model Pipeline                    │
└─────────────────────────────────────────────┘
Input: Z-scores [z5, zBP, z3]
    ↓
OutlierClipper(quantile=0.999)  ← Clips in z-space, very permissive
  clips at 99.9th percentile
    ↓
RobustScaler(with_centering=False)  ← SECOND SCALING
  z' = z / IQR(z)
    ↓
BothEndsStrongTransformer
  → [z5', zBP', z3', min_all, neg_absdiff×3]
    ↓
LinearSVC
    ↓
Prediction

PROBLEMS:
❌ Two scalers (redundant, confusing)
❌ Outlier clipping too weak (99.9%)
❌ Cross-species: z=11.45 passes through
❌ Extreme values overwhelm balance features
```

---

## AFTER (New - Single Scaling)

```
┌─────────────────────────────────────────────┐
│ Model Pipeline (everything inside)         │
└─────────────────────────────────────────────┘
Input: Raw LLRs [s5, sBP, s3]
    ↓
ZeroAnchoredRobustScaler  ← ONLY SCALING STEP
  z = s / median(|s|)
  (fitted on human raw LLRs)
    ↓
SymmetricClipper(zmax=4.0)  ← Aggressive clipping in z-space
  z_clipped = clip(z, -zmax, +zmax)
  (hyperparameter: grid search optimal zmax)
    ↓
SaturatingTransform(enabled=True)  ← Optional compression
  z_sat = sign(z) * log(1 + |z|)
  (hyperparameter: grid search on/off)
    ↓
BothEndsStrongTransformer  ← Unchanged
  → [z_sat, ..., min_all, neg_absdiff×3]
    ↓
LinearSVC  ← Unchanged
    ↓
Prediction

BENEFITS:
✅ Single scaler (clear, simple)
✅ Aggressive outlier control (zmax=3.5-5.0)
✅ Cross-species: z=11.45 → clips to 4.0 → log(5)=1.6
✅ Balance features can now compete
✅ Precision-focused training (F_β=0.5)
```

---

## Cross-Species Example: C. elegans Intron

### BEFORE

```
Raw LLR (3'SS): 14.46
    ↓
Stage 1: 14.46 / 1.26 (human scale) = 11.45
    ↓
OutlierClipper: 11.45 < cap @ 99.9% → passes through
    ↓
Stage 2: 11.45 / 1.5 (IQR) = 7.63
    ↓
SVM contribution: 7.63 × 0.13 = +0.99
min_all penalty:  -0.15 × 0.26 = -0.04
imbalance penalty: -7.78 × 0.07 = -0.54
    ↓
Net: +0.99 - 0.04 - 0.54 = +0.41
    ↓
Result: FALSE POSITIVE (100% U12)
```

### AFTER

```
Raw LLR (3'SS): 14.46
    ↓
Scaler: 14.46 / 1.26 (human scale) = 11.45
    ↓
SymmetricClipper(zmax=4.0): 11.45 → clipped to 4.0
    ↓
SaturatingTransform: sign(4.0) * log(1+4.0) = +1.61
    ↓
SVM contribution: 1.61 × w3 ≈ +0.21
min_all penalty:  (recomputed on clipped/sat values)
imbalance penalty: (recomputed)
    ↓
Net: Penalties can now compete
    ↓
Result: CORRECT (rejected, or low probability)
```

**Key Insight:** Clipping at 4σ + saturation reduces extreme outlier from +0.99 contribution to ~+0.21, making balance penalties competitive.

---

## Grid Search Changes

### BEFORE
```yaml
param_grid:
  estimator__augment__include_pairwise_mins: [false]
  estimator__augment__include_max: [false]
  estimator__svc__penalty: ['l2']
  estimator__svc__dual: [false]
  estimator__svc__C: [0.001, 0.01, 0.1]  # 3 values (example)

# Combinations per round: 1 × 1 × 1 × 1 × 3 = 3
```

### AFTER
```yaml
param_grid:
  estimator__clip__zmax: [3.5, 4.0, 5.0]           # NEW: 3 values
  estimator__saturate__enabled: [false, true]       # NEW: 2 values
  estimator__augment__include_pairwise_mins: [false]
  estimator__augment__include_max: [false]
  estimator__svc__penalty: ['l2']
  estimator__svc__dual: [false]
  estimator__svc__C: [0.001, 0.01, 0.1]            # Keep same or reduce

# Combinations per round: 3 × 2 × 1 × 1 × 1 × 1 × 3 = 18
# If too slow, reduce C to 2 values: 3 × 2 × 3 × ... = 12
```

**Mitigation:** Use parallel CV (`--cv-processes`) or reduce C search space.

---

## Code Changes Summary

### Files to Modify

| File | Change |
|------|--------|
| `classification/clipping.py` | Rename OutlierClipper → SymmetricClipper, change to zmax-based |
| **NEW:** `classification/saturating.py` | Create SaturatingTransform class |
| `classification/optimizer.py` | Update pipeline: scaler first, new clip/saturate steps |
| `classification/trainer.py` | Same pipeline changes |
| `classification/predictor.py` | Extract raw LLRs instead of z-scores |
| `scoring/normalizer.py` | Remove scaling logic (passthrough raw scores) |
| `config/training_default.yaml` | Add clip__zmax, saturate__enabled hyperparameters |
| `utils/model_io.py` | Handle new pipeline structure |
| `tools/inspect_model.py` | Show new pipeline steps |

### New Files

| File | Purpose |
|------|---------|
| `classification/saturating.py` | SaturatingTransform class |
| `docs/SCALING_ARCHITECTURE_V2.md` | Architecture documentation |
| `tests/test_new_scaling.py` | Test new pipeline components |

### Estimated Lines of Code

- New code: ~300 lines (SymmetricClipper refactor, SaturatingTransform, tests)
- Modified code: ~500 lines (pipeline updates, predictor changes, configs)
- Documentation: ~1000 lines (architecture docs, migration guide)

---

## Migration Path for Users

### Old Models (Double Scaling)

```python
# model.pkl contains:
# - Pipeline with RobustScaler (Stage 2)
# - Expects z-score inputs

# Detection:
if 'scale' in pipeline.named_steps and isinstance(..., RobustScaler):
    version = "v1_double_scaling"
```

**Strategy:**
- Detect old models in `model_io.load_model()`
- Raise clear error: "This model uses old double-scaling architecture. Please retrain with new intronIC version."
- Provide migration guide link

### New Models (Single Scaling)

```python
# model.pkl contains:
# - Pipeline with ZeroAnchoredRobustScaler (first step)
# - SymmetricClipper, SaturatingTransform
# - Expects raw LLR inputs

# Detection:
if 'scale' in pipeline.named_steps and isinstance(..., ZeroAnchoredRobustScaler):
    version = "v2_single_scaling"
```

**No migration needed:** New models work with new code out of the box.

---

## Decision Points for User

Before implementing, please decide:

### 1. Clipping Strategy
- [ ] **Option A:** Fixed zmax=[3.5, 4.0, 5.0] (recommended)
- [ ] **Option B:** Quantile-based [0.95, 0.975, 0.99]
- [ ] **Option C:** Hybrid (fit from quantile, store fixed)

### 2. Saturating Transform
- [ ] **Option A:** Always enabled (no hyperparameter)
- [ ] **Option B:** Grid search enabled=[true, false] (recommended)

### 3. Scoring Function
- [ ] **Option A:** F_0.5 score (recommended)
- [ ] **Option B:** Precision at 90% recall
- [ ] **Option C:** Custom weighted

### 4. ScoreNormalizer
- [ ] **Option A:** Keep as passthrough (minimal changes)
- [ ] **Option B:** Deprecate, create RawScorer
- [ ] **Option C:** Remove entirely

### 5. Grid Search Size
- [ ] Keep current C search (3-5 values) → larger grid
- [ ] Reduce C search (2-3 values) → manageable grid
- [ ] Use RandomizedSearchCV instead of GridSearchCV

**Recommended defaults (based on expert feedback):**
- 1A: Fixed zmax
- 2B: Grid search saturate
- 3A: F_0.5 score
- 4A: Keep ScoreNormalizer as passthrough
- 5B: Reduce C search to keep grid manageable
