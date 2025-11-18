# Scaler Centering Fix - Implementation Complete

**Date:** 2025-01-17
**Status:** ✅ Implementation complete, ready for retraining

## Problem Summary

The previous implementation had two critical issues:

### 1. Train/Test Mismatch
- **Training** (trainer.py): Extracted z-scores from introns (`five_z_score`, `bp_z_score`, `three_z_score`)
- **Prediction** (predictor.py): Extracted raw scores (`five_raw_score`, `bp_raw_score`, `three_raw_score`)
- **Result**: Pipeline with RobustScaler caused double-scaling during training, wrong data during prediction

### 2. Double-Scaling Problem
- **Stage 1**: ScoreNormalizer transformed raw LLRs → z-scores
- **Stage 2**: Pipeline's RobustScaler re-scaled z-scores → double transformation
- **Result**: Features seen by model differed between training and prediction

## Solution Implemented

Following expert guidance, we implemented a corrected architecture with **single scaling step**:

### Data Flow

```
Raw PWM Scores (LLRs)
         ↓
ScoreNormalizer (RobustScaler with centering=True)
         ↓
Z-Scores [five_z_score, bp_z_score, three_z_score]
         ↓
BothEndsStrongTransformer
         ↓
Augmented Features (7D)
    [five_z, bp_z, three_z,           # Pass-through z-scores
     min_all,                           # min(five_z, bp_z, three_z)
     neg_absdiff_5_bp,                  # -|five_z - bp_z|
     neg_absdiff_5_3,                   # -|five_z - three_z|
     neg_absdiff_bp_3]                  # -|bp_z - three_z|
         ↓
LinearSVC (L2-regularized, balanced class weights)
         ↓
CalibratedClassifierCV (sigmoid or isotonic)
         ↓
U12 Probability [0-100%]
```

## Changes Made

### 1. scoring/normalizer.py
- **Changed**: Now uses `RobustScaler(with_centering=True)` instead of `ZeroAnchoredRobustScaler`
- **Reason**: Centering removes composition bias, as recommended by expert
- **Location**: Line 294-297

### 2. classification/trainer.py
- **Removed**: RobustScaler from pipeline (lines 197-212)
- **Added**: BothEndsStrongTransformer to pipeline
- **Data**: Extracts z-scores from introns (lines 287-294)
- **Result**: Pipeline = `BothEndsStrongTransformer → LinearSVC → CalibratedClassifierCV`

### 3. classification/optimizer.py
- **Removed**: RobustScaler from all 3 pipeline locations
  - Stage 1 optimization (line 562-576)
  - Calibration evaluation (line 822-837)
  - CV evaluation (line 898-913)
- **Added**: BothEndsStrongTransformer to all 3 locations
- **Parameters**: `include_max=False`, `include_pairwise_mins=False` (7D feature space)

### 4. classification/predictor.py
- **Changed**: Extracts z-scores instead of raw scores (lines 70-79)
- **Removed**: Scaler extraction code (was lines 76-86)
- **Removed**: Z-score computation code (was lines 88-97)
- **Reason**: Introns already have z-scores from ScoreNormalizer
- **Data**: Passes z-scores directly to model (line 86)

## Verification

All architecture tests passed:

```
✓ Pipeline creation works
✓ Feature transformation: 3D z-scores → 7D augmented features
✓ No double-scaling (no 'scale' step in pipeline)
✓ End-to-end training and prediction work
```

See `test_corrected_architecture.py` for full test suite.

## Key Principles

1. **Single Scaling Step**: Scaling happens OUTSIDE pipeline via ScoreNormalizer
2. **Consistent Features**: Both trainer and predictor use z-scores
3. **No Double-Scaling**: Pipeline has NO scaler (would cause double-scaling)
4. **Domain Adaptation**: ScoreNormalizer can be refitted per-species

## Expected Improvements

Based on expert guidance and previous experiments:

1. **Reduced False Positives**: Centering + augmented features should suppress one-end-strong FPs
   - Previous: 130 FPs on C. elegans (no centering, no augmented features)
   - With centering only: 49 FPs on C. elegans
   - Expected with centering + augmented features: ~6-8 FPs (matching old v1.5.1)

2. **Better Calibration**: Single scaling step prevents distribution shift between train/test

3. **Improved Generalization**: Augmented features provide explicit penalties for imbalanced signals

## Next Steps

1. ✅ Architecture implemented
2. ✅ Tests passing
3. ⏳ **Retrain model** with corrected architecture
4. ⏳ **Validate on C. elegans** (expect ~6-8 FPs)
5. ⏳ Compare performance with previous models

## Training Command

```bash
pixi run intronIC train \
    -n homo_sapiens_corrected \
    -o . \
    --n_optimization_rounds 3 \
    --n_cv_folds 5 \
    -p 10
```

## References

- Expert workflow document (provided by user)
- SCALER_ARCHITECTURE_REVIEW.md
- AUGMENTED_FEATURES_EXPERIMENT_PLAN.md
