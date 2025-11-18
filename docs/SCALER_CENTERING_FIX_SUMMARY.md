# Scaler Centering Fix - Complete Summary

## Problem Identified

After implementing the new scaling architecture, the model trained successfully on human data (F_0.75=0.9988) but **catastrophically failed on C. elegans validation**:

- **99.87% of C. elegans introns** classified with 0% U12 probability
- Only 135/109,820 introns got 100% score
- Complete loss of discriminative power

**Note:** C. elegans has **0 U12-type introns** (species lost the U12 spliceosome), so the extreme bimodal distribution (0% or 100%) indicates the model lost ability to discriminate based on actual sequence features.

## Root Cause Analysis

### Cross-Species z-Score Disparity

**Human training data** (41% GC):
- z-score range: approximately -3 to +3
- Trained clipping thresholds: ±1.48σ, ±2.91σ, ±2.72σ

**C. elegans** (38% GC):
- z-score range: **-63.9 to +16.6** (extreme due to composition bias)
- After clipping with human thresholds: **ALL** values squashed to [-1.48, +2.91]
- SVM completely lost signal

### The Fundamental Issue

**Training-species clipping thresholds don't transfer to target species.**

The `SymmetricClipper(quantile=0.95)` was fitted on human z-scores, learning numeric boundaries like "clip at ±1.5σ". When applied to C. elegans:
- Human introns cluster around z=0
- C. elegans introns span z ∈ [-64, +16] (18× wider range)
- Clipping destroys ALL discriminative information

## Solution Implemented

### Option 1: Remove Clipping Entirely (CHOSEN)

**Rationale:**
- Zero-preservation alone is sufficient for cross-species deployment
- SVM learns appropriate boundaries from training data
- Simpler architecture, fewer hyperparameters
- Clipping intended to prevent artifacts but created worse ones

### Code Changes

#### 1. config/training_default.yaml
Removed clip_quantile from parameter grid

#### 2. classification/trainer.py
Removed SymmetricClipper from pipeline

#### 3. classification/optimizer.py
- Removed clip_quantile from SVMParameters dataclass
- Removed best_clip_quantile from OptimizationRound dataclass
- Removed SymmetricClipper from pipeline construction
- Removed clip_quantile parameter from _evaluate_params() signature
- Removed all clip_quantile logging statements

## New Architecture (Final)

```
Pipeline: raw LLRs → scale → saturate → augment → svc

Components:
- ZeroAnchoredRobustScaler: z = s / median(|s|)  [ONLY scaling step]
- SaturatingTransform: Optional log compression
- BothEndsStrongTransformer: Augments 3D → 7D with balance features
- LinearSVC: L2-regularized linear classifier
- CalibratedClassifierCV: External calibration
```

**Key Design Principle:** Trust SVM to learn appropriate boundaries from training data.

## Training Status

**Current:** Retraining model without clipping (in progress)

**Parameter Grid:**
- C: 13 values (round 1), 100 values (rounds 2-4)
- saturate_enabled: [false, true]
- Total combinations: 26 (round 1) → 200 (rounds 2-4)
- **Previously:** 78 (round 1) → 600 (rounds 2-4) with 3 clip_quantile values

## Expected Improvements

### C. elegans Validation (0 U12 introns expected)
- Should restore full z-score range utilization
- Expected: **0 high-confidence U12 predictions** (C. elegans has no U12 introns)
- Current broken model: bimodal distribution (99.87% at 0%, 135 at 100%)
- Fixed model: continuous score distribution with very few/no high scores

### All Species
- No more species-specific clipping artifacts
- Simpler deployment (one less hyperparameter)
- Faster training (3× fewer grid search combinations)

## Next Steps

1. Complete training - Wait for final model with no clipping
2. Validate on C. elegans - Should show 0 or very few U12 predictions
3. Measure FP reduction - Compare to baseline
4. Update documentation - CLAUDE.md, config/README.md, etc.

---

**Date:** 2025-11-17
**Status:** In Progress (training without clipping)
