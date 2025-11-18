# Data Flow Through Scaling Pipeline

## Current Architecture

```
┌─────────────────────────────────────────────────────────────────────────┐
│ TRAINING (Human Reference Data)                                         │
└─────────────────────────────────────────────────────────────────────────┘

PWM Scoring
    │
    ├─ 5'SS:  raw LLR₅ (e.g., -35.6)
    ├─ BPS:   raw LLRᵦ (e.g., -4.7)
    └─ 3'SS:  raw LLR₃ (e.g., -0.3)
    │
    v
┌─────────────────────────────────────────────────────────────────────────┐
│ STAGE 1: ZeroAnchoredRobustScaler (Pre-processing)                     │
│   - Fitted on: Human reference raw LLRs                                 │
│   - Algorithm: z = raw / median(|raw|)                                  │
│   - Learned scales: [σ₅=34.3, σᵦ=5.83, σ₃=1.26]  (example)            │
└─────────────────────────────────────────────────────────────────────────┘
    │
    ├─ z₅  = -35.6 / 34.3 = -1.04
    ├─ zᵦ  = -4.7 / 5.83  = -0.81
    └─ z₃  = -0.3 / 1.26  = -0.24
    │
    v
┌─────────────────────────────────────────────────────────────────────────┐
│ MODEL PIPELINE (Stage 2)                                                │
│   ┌──────────────────────────────────────────────────────────────────┐  │
│   │ OutlierClipper(quantile=0.999)                                   │  │
│   │   - Clips at 99.9th percentile of |z|                           │  │
│   │   - Fitted on: Human reference z-scores from Stage 1            │  │
│   │   - Learned caps: [cap₅=4.5, capᵦ=3.2, cap₃=3.8]  (example)   │  │
│   └──────────────────────────────────────────────────────────────────┘  │
│                              │                                           │
│   ┌──────────────────────────v───────────────────────────────────────┐  │
│   │ RobustScaler(with_centering=False)                               │  │
│   │   - Algorithm: z' = z / IQR(z)                                   │  │
│   │   - Fitted on: Human reference z-scores (after clipping)        │  │
│   │   - Learned IQRs: [IQR₅=2.1, IQRᵦ=1.8, IQR₃=1.5]  (example)   │  │
│   └──────────────────────────────────────────────────────────────────┘  │
│                              │                                           │
│   ┌──────────────────────────v───────────────────────────────────────┐  │
│   │ BothEndsStrongTransformer                                        │  │
│   │   - Augments 3D → 7D with computed features                     │  │
│   │   - Features: [z₅', zᵦ', z₃', min_all, neg_absdiff × 3]       │  │
│   └──────────────────────────────────────────────────────────────────┘  │
│                              │                                           │
│   ┌──────────────────────────v───────────────────────────────────────┐  │
│   │ LinearSVC                                                         │  │
│   │   - Learns coefficients: [w₅, wᵦ, w₃, w_min, w_diff × 3]      │  │
│   │   - Trained on: Human reference (80/20 split)                   │  │
│   └──────────────────────────────────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────────────────────┘
    │
    v
Predictions (U12 probability)


┌─────────────────────────────────────────────────────────────────────────┐
│ CROSS-SPECIES DEPLOYMENT (C. elegans, no retraining)                   │
└─────────────────────────────────────────────────────────────────────────┘

PWM Scoring (same matrices as training)
    │
    ├─ 5'SS:  raw LLR₅ (e.g., -5.3)   ← Different distribution!
    ├─ BPS:   raw LLRᵦ (e.g., 1.8)    ← C. elegans has different
    └─ 3'SS:  raw LLR₃ (e.g., 14.5)   ← composition bias (38% GC)
    │
    v
┌─────────────────────────────────────────────────────────────────────────┐
│ STAGE 1: ZeroAnchoredRobustScaler (HUMAN scales!)                      │
│   - Uses: Human-fitted scales [σ₅=34.3, σᵦ=5.83, σ₃=1.26]           │
│   - Algorithm: z = raw / σ_human                                       │
│   - Problem: C. elegans distribution ≠ human distribution             │
└─────────────────────────────────────────────────────────────────────────┘
    │
    ├─ z₅  = -5.3 / 34.3  = -0.15    ← OK
    ├─ zᵦ  =  1.8 / 5.83  =  0.30    ← OK
    └─ z₃  = 14.5 / 1.26  = 11.45    ← EXTREME! (>10σ)
    │
    v
┌─────────────────────────────────────────────────────────────────────────┐
│ MODEL PIPELINE (HUMAN-fitted components!)                               │
│   ┌──────────────────────────────────────────────────────────────────┐  │
│   │ OutlierClipper(quantile=0.999)                                   │  │
│   │   - Caps: [4.5, 3.2, 3.8] from human training                   │  │
│   │   - 11.45 > 3.8, but quantile=0.999 is too permissive!         │  │
│   │   - Value passes through without clipping                       │  │
│   └──────────────────────────────────────────────────────────────────┘  │
│                              │                                           │
│   ┌──────────────────────────v───────────────────────────────────────┐  │
│   │ RobustScaler (HUMAN IQRs!)                                       │  │
│   │   - IQRs: [2.1, 1.8, 1.5] from human training                   │  │
│   │   - z₃' = 11.45 / 1.5 = 7.63  (still extreme!)                 │  │
│   └──────────────────────────────────────────────────────────────────┘  │
│                              │                                           │
│   ┌──────────────────────────v───────────────────────────────────────┐  │
│   │ BothEndsStrongTransformer                                        │  │
│   │   - Computes: min_all = min(-0.15, 0.30, 7.63) = -0.15         │  │
│   │   - Computes: neg_absdiff_5_3 = -|(-0.15) - 7.63| = -7.78     │  │
│   └──────────────────────────────────────────────────────────────────┘  │
│                              │                                           │
│   ┌──────────────────────────v───────────────────────────────────────┐  │
│   │ LinearSVC (HUMAN coefficients!)                                  │  │
│   │   - z₃' × w₃ = 7.63 × 0.13 = +0.99  ← Dominates!              │  │
│   │   - min_all × w_min = -0.15 × 0.26 = -0.04  ← Weak penalty!   │  │
│   │   - neg_absdiff × w_diff = -7.78 × 0.07 = -0.54               │  │
│   │   - Net: +0.99 - 0.04 - 0.54 = +0.41 → Predicts U12!          │  │
│   └──────────────────────────────────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────────────────────┘
    │
    v
FALSE POSITIVE: 100% U12 probability
(min_all is negative but penalties too weak to overcome extreme z₃)


┌─────────────────────────────────────────────────────────────────────────┐
│ PROPOSED: Domain Adaptation (refit Stage 1 only)                       │
└─────────────────────────────────────────────────────────────────────────┘

PWM Scoring (same)
    │
    └─ 3'SS: raw LLR₃ = 14.5
    │
    v
┌─────────────────────────────────────────────────────────────────────────┐
│ STAGE 1: ZeroAnchoredRobustScaler (C. ELEGANS scales!)                 │
│   - Refit on: C. elegans unlabeled data                                │
│   - New scale: σ₃_celegans ≈ 5.2  (hypothetical)                      │
│   - Algorithm: z = 14.5 / 5.2 = 2.79  ← More reasonable!              │
└─────────────────────────────────────────────────────────────────────────┘
    │
    └─ z₃ = 2.79 (within normal range)
    │
    v
┌─────────────────────────────────────────────────────────────────────────┐
│ MODEL PIPELINE (HUMAN-fitted, but now receives reasonable inputs)      │
│   - OutlierClipper: 2.79 < cap (passes through, as expected)          │
│   - RobustScaler: z₃' = 2.79 / 1.5 = 1.86                             │
│   - BothEndsStrong: min_all, imbalance features                        │
│   - LinearSVC: Balanced decision (penalties can work properly)         │
└─────────────────────────────────────────────────────────────────────────┘
    │
    v
CORRECT CLASSIFICATION (balanced features, no extreme outliers)
```

## Key Questions

1. **Why two scalers?**
   - Is Stage 2 redundant if Stage 1 already normalizes?
   - Or do they serve different purposes?

2. **Cross-species adaptation:**
   - Refit Stage 1 only (domain adaptation)?
   - Refit Stage 2 also (full retraining)?
   - Or architectural change?

3. **Outlier handling:**
   - Should OutlierClipper be more aggressive?
   - Should it operate on Stage 1 output instead of Stage 2 input?

4. **with_centering=False:**
   - Needed for Stage 1 (preserve semantic zero of LLRs)
   - Still needed for Stage 2 (operating on already zero-anchored z-scores)?
