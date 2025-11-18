# Scaling Redesign Implementation Progress

## Session Summary
**Date:** 2025-01-16
**Design Decisions:**
- Clipping: Quantile-based [0.95, 0.975, 0.99]
- Saturating transform: Grid search [False, True]
- Scoring: **F_0.75** (precision ~1.56× more weight than recall)
- ScoreNormalizer: Deprecate
- Grid size: Keep full C search (no reduction)

---

## ✅ Completed Tasks

### 1. SymmetricClipper (classification/clipping.py)
- ✅ Renamed `OutlierClipper` → `SymmetricClipper`
- ✅ Updated docstrings for z-space operation
- ✅ Changed default quantile: 0.999 → 0.975
- ✅ Added backward compatibility alias
- ✅ Updated all examples and comments

### 2. SaturatingTransform (classification/saturating.py)
- ✅ Created new transformer class
- ✅ Implements `f(z) = sign(z) * log(1 + |z|)`
- ✅ Supports enable/disable via hyperparameter
- ✅ Includes inverse_transform for interpretability
- ✅ Comprehensive docstrings with examples

### 3. Optimizer Pipeline Updates (classification/optimizer.py)
- ✅ Updated imports:
  - Added: `ZeroAnchoredRobustScaler`, `SymmetricClipper`, `SaturatingTransform`
  - Added: `make_scorer`, `fbeta_score`
  - Removed: `RobustScaler`

- ✅ Updated dataclasses:
  - `SVMParameters`: Added `clip_quantile`, `saturate_enabled`, `include_pairwise_mins`
  - `OptimizationRound`: Added same new fields, kept old ones as defaults

- ✅ Updated `_run_cv_round()` method:
  - New pipeline: `ZeroAnchoredRobustScaler` → `SymmetricClipper` → `SaturatingTransform` → `BothEndsStrongTransformer` → `LinearSVC`
  - New param_grid with clip/saturate hyperparameters
  - F_0.75 scoring instead of neg_log_loss
  - Extract new best parameters from grid search
  - Return updated OptimizationRound

- ✅ Updated `optimize()` method:
  - Extract new final parameters from best round
  - Update `_evaluate_params()` call signature
  - Update print statement with new parameters
  - Return `SVMParameters` with all new fields

- ✅ Updated `_evaluate_params()` method:
  - New signature with clip_quantile, saturate_enabled, include_pairwise_mins
  - New pipeline matching `_run_cv_round()`
  - F_0.75 scoring instead of neg_log_loss
  - Updated docstrings

### 4. Trainer Pipeline Updates (classification/trainer.py)
- ✅ Updated imports:
  - Added: `ZeroAnchoredRobustScaler`, `SymmetricClipper`, `SaturatingTransform`
  - Removed: `RobustScaler`, `OutlierClipper`

- ✅ Updated `_train_single_model()` pipeline:
  - New: `ZeroAnchoredRobustScaler` → `SymmetricClipper` → `SaturatingTransform` → `BothEndsStrongTransformer` → `LinearSVC`
  - Uses parameters from SVMParameters dataclass (clip_quantile, saturate_enabled, include_pairwise_mins)
  - Fixed hyperparameters: penalty='l2', dual=False, loss='squared_hinge'

### 5. Predictor Updates (classification/predictor.py)
- ✅ Updated module docstring with new architecture notes
- ✅ Updated class docstring for SVMPredictor
- ✅ Updated `_prepare_features()` method:
  - Extracts RAW LLR scores (five_raw_score, bp_raw_score, three_raw_score)
  - No longer extracts z-scores (pipeline handles scaling)
  - Updated docstrings and error messages

- ✅ Updated `_predict_chunk()` method:
  - Extracts raw scores for pipeline input
  - **Extracts fitted scaler from pipeline**
  - **Computes z-scores and stores on introns** (for output files)
  - Uses introns with z-scores for BothEndsStrong feature computation

- ✅ Updated `_predict_chunk_worker()` function:
  - Same changes as `_predict_chunk()`
  - Extracts scaler, computes z-scores, stores on introns

- ✅ Updated `predict()` method docstring:
  - Notes that introns must have raw PWM scores
  - Pipeline handles scaling internally

**Key improvement:** Z-scores are now **extracted from the fitted pipeline** and stored on intron objects for output files (.scores.iic), ensuring output shows exactly what the model saw after scaling (before clipping/saturation)

### 6. Config File Updates (config/training_default.yaml)
- ✅ Updated param_grid with new hyperparameters:
  - Added: `estimator__clip__quantile: [0.95, 0.975, 0.99]`
  - Added: `estimator__saturate__enabled: [false, true]`
  - Fixed: `estimator__augment__include_pairwise_mins: [false]`
  - Fixed: `estimator__augment__include_max: [false]`
  - Fixed: `estimator__svc__penalty: ['l2']`
  - Simplified: `method: ['isotonic']`
  - Removed: dual, loss, intercept_scaling (now fixed parameters)

- ✅ Updated documentation:
  - Comprehensive comments for clip_quantile and saturate_enabled
  - Explained new architecture changes (2025 redesign)
  - Grid size: 6 combinations per C value (3 × 2 × 1 × 1 × 1 × 1)

---

## 🚧 In-Progress Tasks

---

## 📋 Remaining Tasks (Not Started)

### 7. ScoreNormalizer Deprecation (scoring/normalizer.py)
- ✅ Added comprehensive deprecation warning to class docstring
- ✅ Clarified: NOT used in main pipeline (deprecated for training/prediction)
- ✅ Provided guidance: Extract scaler from pipeline for new code
- ✅ Documented valid use cases:
  - Standalone z-score computation for analysis/debugging
  - Extracting z-scores from old models
  - Comparing old vs new normalization
- ✅ Example code showing how to extract scaler from pipeline

### 8. Model I/O Updates (utils/model_io.py)
- ✅ Added `pipeline_architecture: 'single_scaler_v2025'` to metadata
- ✅ Current joblib serialization works fine with new architecture:
  - SVMEnsemble with new SVMParameters fields pickles correctly
  - Pipeline with new transformers (ZeroAnchoredRobustScaler, SymmetricClipper, SaturatingTransform) pickles correctly
  - No backward compatibility needed (will only use new models)

### 9. Write tests
- Test SymmetricClipper
- Test SaturatingTransform
- Test full pipeline with raw LLRs
- Test cross-species scenario

### 10. Retrain and validate
- Retrain homo_sapiens model
- Validate on C. elegans
- Compare old vs new FP rates

### 11. Update documentation
- CLAUDE.md
- config/README.md
- Migration guide

---

## Quick Resume Guide

**To continue implementation:**

1. **Finish optimizer.py updates** (7 specific locations listed above)
2. **Apply same changes to trainer.py**
3. **Update predictor.py** to extract raw LLRs
4. **Test the pipeline** with a small dataset
5. **Full retrain and validation**

**Key files modified so far:**
- ✅ `classification/clipping.py` (SymmetricClipper)
- ✅ `classification/saturating.py` (new file)
- 🚧 `classification/optimizer.py` (partial)

**Next file to modify:**
- `classification/optimizer.py` (finish remaining 7 updates)
