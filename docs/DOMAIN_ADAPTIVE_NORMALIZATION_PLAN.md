# Domain-Adaptive Normalization - Implementation Plan

**Date:** 2025-11-16
**Status:** EXPERT-REVIEWED DESIGN
**Context:** Post-clipping removal, implementing robust cross-species normalization

---

## Executive Summary

**Problem:** Human-trained numeric clipping thresholds (±1.5σ) fail catastrophically on species with different nucleotide composition (C. elegans: z ∈ [-63, +16]).

**Expert Solution:** Make scaling AND clipping fully domain-adaptive:
1. **Train quantile hyperparameter `q`** (e.g., 0.95) on human via CV
2. **Fit scaler + compute quantiles** on target species' data (unsupervised)
3. **Classify with human-trained SVM** (supervised component unchanged)

**Key Insight:** Only the quantile rule `q` transfers across species, not the numeric thresholds. Each species gets its own scale and clip bounds computed from that species' data.

---

## Core Principle: Dimensionless Transfer

### What Transfers Across Species (Hyperparameters)
- **Quantile rule**: `q = 0.95` (or 0.975, 0.99)
- **SVM weights**: Trained on human U12/U2 labels
- **Feature augmentation**: BothEndsStrong architecture
- **Saturation rule**: Enabled/disabled

### What Does NOT Transfer (Fitted Parameters)
- ❌ **Numeric scale**: `median(|s|)` per feature
- ❌ **Numeric clip thresholds**: `lo_j`, `hi_j` per feature
- ❌ **Raw score distributions**: Species-specific

### Implementation
```python
# TRAINING (Human)
# 1. Learn hyperparameters via CV
q = optimize_quantile([0.95, 0.975, 0.99])  # → 0.95 (example)

# 2. Fit normalizer on human data
scale_j_human = median(|s_j|)  # per feature
z_j_human = s_j / scale_j_human
lo_j_human = quantile(z_j_human, 1 - q)
hi_j_human = quantile(z_j_human, q)
z_j_clipped_human = clip(z_j_human, lo_j_human, hi_j_human)

# 3. Train SVM on normalized human data
svm.fit(z_j_clipped_human, labels_human)

# 4. Save hyperparameters + human SVM
save_model({
    'q': 0.95,  # Hyperparameter (transfers)
    'svm': svm,  # Trained on human (transfers)
    'human_normalizer': {  # For reference only
        'scale': scale_j_human,
        'lo': lo_j_human,
        'hi': hi_j_human
    }
})


# DEPLOYMENT (C. elegans)
# 1. Fit normalizer on C. elegans data (UNSUPERVISED)
scale_j_celegans = median(|s_j|)  # Different from human!
z_j_celegans = s_j / scale_j_celegans
lo_j_celegans = quantile(z_j_celegans, 1 - q)  # Use same q=0.95
hi_j_celegans = quantile(z_j_celegans, q)
z_j_clipped_celegans = clip(z_j_celegans, lo_j_celegans, hi_j_celegans)

# 2. Classify with human SVM (SUPERVISED)
predictions = svm.predict(z_j_clipped_celegans)
```

**Result:** C. elegans introns with z ∈ [-63, +16] get scaled to ~[-3, +3], then clipped at C. elegans' 0.95 quantile (not human's ±1.5).

---

## Architecture: Domain-Adaptive Normalizer

### Old Architecture (FAILED on C. elegans)
```
Training (Human):
  raw LLRs → scale (human) → clip (human ±1.5σ) → SVM

Deployment (C. elegans):
  raw LLRs → scale (human) → clip (human ±1.5σ) → SVM
                                      ↑
                            CATASTROPHIC FAILURE
                            (clips -63 to -1.5, destroys signal)
```

### New Architecture (DOMAIN-ADAPTIVE)
```
Training (Human):
  raw LLRs → ZeroAnchoredRobustScaler (fit on human)
           → SymmetricClipper (q=0.95, fit on human z-scores)
           → SaturatingTransform (optional)
           → BothEndsStrong
           → LinearSVC (fit on human labels)

  Save: {q, saturate_enabled, SVM weights}

Deployment (C. elegans):
  raw LLRs → ZeroAnchoredRobustScaler (fit on C. elegans)
           → SymmetricClipper (q=0.95, fit on C. elegans z-scores)
           → SaturatingTransform (same rule)
           → BothEndsStrong (same architecture)
           → LinearSVM (same weights from human)

  Normalizer fitted on: ALL C. elegans introns (unsupervised)
  SVM applied: Human-trained weights
```

**Key Change:** Normalizer (scale + clip) is **species-conditional**, SVM is **species-invariant**.

---

## Implementation Design

### Class: DomainAdaptiveNormalizer

```python
class DomainAdaptiveNormalizer:
    """
    Unsupervised normalizer that fits scale + clipping bounds on target domain.

    Hyperparameters (trained on human, transfer to all species):
    - q: Clipping quantile (0.95, 0.975, or 0.99)
    - symmetric: Whether to use symmetric clipping (True recommended)

    Fitted Parameters (computed per species):
    - scale_: Per-feature robust scale (median(|s|))
    - clip_lo_: Per-feature lower clipping bound (q-quantile of z)
    - clip_hi_: Per-feature upper clipping bound ((1-q)-quantile of z)
    """

    def __init__(self, quantile=0.95, symmetric=True):
        """
        Args:
            quantile: Clipping quantile (hyperparameter from human CV)
            symmetric: Use symmetric clipping (±threshold)
        """
        self.quantile = quantile
        self.symmetric = symmetric

        # Fitted parameters (set by fit())
        self.scale_ = None
        self.clip_lo_ = None
        self.clip_hi_ = None
        self.fitted_on_species_ = None
        self.fitted_on_n_samples_ = None

    def fit(self, X, species_name=None, max_samples=200000):
        """
        Fit scaler + clipping bounds on unlabeled data.

        Expert recommendation: Cap at 100k-200k introns to avoid computational issues
        while maintaining robust statistics.

        Args:
            X: Raw LLR scores (n_samples, n_features)
            species_name: For logging/metadata
            max_samples: Maximum samples to use for fitting (random sample if exceeded)

        Returns:
            self
        """
        # Step 0: Random sampling if dataset too large
        # Expert: "Sample randomly if the species has millions of introns"
        n_samples = X.shape[0]
        if n_samples > max_samples:
            indices = np.random.choice(n_samples, max_samples, replace=False)
            X_fit = X[indices]
            self.fitted_on_n_samples_total_ = n_samples
            self.fitted_on_n_samples_used_ = max_samples
        else:
            X_fit = X
            self.fitted_on_n_samples_total_ = n_samples
            self.fitted_on_n_samples_used_ = n_samples

        # Step 1: Zero-anchored robust scaling
        # Scale each feature by median(|x|) to preserve zero
        # Expert: "per-feature quantiles (q per feature) rather than global"
        self.scale_ = np.median(np.abs(X_fit), axis=0)
        self.scale_[self.scale_ == 0] = 1.0  # Avoid division by zero

        # Step 2: Transform to z-scores
        Z = X_fit / self.scale_

        # Step 3: Compute clipping thresholds at this species' PER-FEATURE quantiles
        # Expert: "per-feature quantiles rather than a single global threshold"
        if self.symmetric:
            # Symmetric: use max(|lo|, |hi|) for both bounds
            lo = np.quantile(Z, 1 - self.quantile, axis=0)  # Per-feature
            hi = np.quantile(Z, self.quantile, axis=0)      # Per-feature
            threshold = np.maximum(np.abs(lo), np.abs(hi))
            self.clip_lo_ = -threshold
            self.clip_hi_ = threshold
        else:
            # Asymmetric: use actual quantiles per feature
            self.clip_lo_ = np.quantile(Z, 1 - self.quantile, axis=0)  # Per-feature
            self.clip_hi_ = np.quantile(Z, self.quantile, axis=0)      # Per-feature

        # Metadata
        self.fitted_on_species_ = species_name
        self.fitted_on_timestamp_ = datetime.now().isoformat()

        # Compute hash for reproducibility tracking
        # Expert: "Store hash/version so you can tell if two runs used same scaling"
        self.normalizer_hash_ = self._compute_hash()

        return self

    def transform(self, X):
        """
        Scale and clip using fitted parameters.

        Args:
            X: Raw LLR scores (n_samples, n_features)

        Returns:
            Z_clipped: Normalized and clipped scores
        """
        # Scale
        Z = X / self.scale_

        # Clip
        Z_clipped = np.clip(Z, self.clip_lo_, self.clip_hi_)

        return Z_clipped

    def fit_transform(self, X, species_name=None):
        """Fit and transform in one step."""
        return self.fit(X, species_name).transform(X)

    def _compute_hash(self):
        """
        Compute hash of fitted parameters for version tracking.

        Expert: "Store hash/version so you can tell if two runs used same scaling"
        """
        import hashlib
        import json

        # Create deterministic representation of fitted parameters
        params = {
            'quantile': self.quantile,
            'symmetric': self.symmetric,
            'scale': self.scale_.tolist(),
            'clip_lo': self.clip_lo_.tolist(),
            'clip_hi': self.clip_hi_.tolist()
        }

        # Compute SHA256 hash
        params_str = json.dumps(params, sort_keys=True)
        return hashlib.sha256(params_str.encode()).hexdigest()[:16]

    def get_metadata(self):
        """
        Return fitted parameters for reproducibility.

        Expert: "Log which mode was used, with N_introns used for fitting"
        """
        return {
            'quantile': self.quantile,
            'symmetric': self.symmetric,
            'scale': self.scale_.tolist(),
            'clip_lo': self.clip_lo_.tolist(),
            'clip_hi': self.clip_hi_.tolist(),
            'fitted_on_species': self.fitted_on_species_,
            'fitted_on_n_samples_total': self.fitted_on_n_samples_total_,
            'fitted_on_n_samples_used': self.fitted_on_n_samples_used_,
            'fitted_on_timestamp': self.fitted_on_timestamp_,
            'normalizer_hash': self.normalizer_hash_
        }

    @classmethod
    def from_metadata(cls, metadata):
        """Recreate normalizer from saved metadata (for reproducibility)."""
        normalizer = cls(
            quantile=metadata['quantile'],
            symmetric=metadata['symmetric']
        )
        normalizer.scale_ = np.array(metadata['scale'])
        normalizer.clip_lo_ = np.array(metadata['clip_lo'])
        normalizer.clip_hi_ = np.array(metadata['clip_hi'])
        normalizer.fitted_on_species_ = metadata['fitted_on_species']
        normalizer.fitted_on_n_samples_ = metadata['fitted_on_n_samples']
        return normalizer
```

### Integration with Pipeline

```python
# Training (Human)
def train_model(u12_introns, u2_introns, config):
    # Extract raw LLR scores
    X_u12 = extract_scores(u12_introns)  # (387, 3)
    X_u2 = extract_scores(u2_introns)    # (20690, 3)
    X = np.vstack([X_u12, X_u2])
    y = np.array([1]*len(X_u12) + [0]*len(X_u2))

    # Grid search over quantile (hyperparameter optimization)
    best_q = None
    best_score = -np.inf

    for q in [0.95, 0.975, 0.99]:
        # Create normalizer with this quantile
        normalizer = DomainAdaptiveNormalizer(quantile=q)

        # Create pipeline
        pipeline = Pipeline([
            ('normalize', normalizer),
            ('saturate', SaturatingTransform(enabled=config.saturate)),
            ('augment', BothEndsStrongTransformer(...)),
            ('svc', LinearSVC(...))
        ])

        # Cross-validate
        scores = cross_val_score(pipeline, X, y, cv=5, scoring='f_beta_0.75')
        mean_score = np.mean(scores)

        if mean_score > best_score:
            best_score = mean_score
            best_q = q

    # Train final model with best quantile
    normalizer = DomainAdaptiveNormalizer(quantile=best_q)
    normalizer.fit(X, species_name='homo_sapiens')

    pipeline = Pipeline([
        ('normalize', normalizer),
        ('saturate', SaturatingTransform(enabled=config.saturate)),
        ('augment', BothEndsStrongTransformer(...)),
        ('svc', LinearSVC(...))
    ])
    pipeline.fit(X, y)

    # Save model + metadata
    return {
        'pipeline': pipeline,
        'hyperparameters': {
            'quantile': best_q,
            'saturate_enabled': config.saturate
        },
        'human_normalizer_metadata': normalizer.get_metadata()
    }


# Deployment (C. elegans)
def classify_species(introns, model, mode='domain_adaptive'):
    # Extract raw LLR scores
    X = extract_scores(introns)  # (109820, 3)

    if mode == 'domain_adaptive':
        # Fit NEW normalizer on C. elegans data (unsupervised)
        quantile = model['hyperparameters']['quantile']  # Use same q from human
        normalizer = DomainAdaptiveNormalizer(quantile=quantile)
        normalizer.fit(X, species_name='caenorhabditis_elegans')

        # Create pipeline with C. elegans normalizer + human SVM
        pipeline = Pipeline([
            ('normalize', normalizer),
            ('saturate', model['pipeline'].named_steps['saturate']),  # Same rule
            ('augment', model['pipeline'].named_steps['augment']),    # Same arch
            ('svc', model['pipeline'].named_steps['svc'])             # Human weights
        ])

    elif mode == 'human_fixed':
        # Use human normalizer (frozen)
        pipeline = model['pipeline']  # Use as-is

    # Classify
    predictions = pipeline.predict_proba(X)

    return predictions, normalizer.get_metadata()
```

---

## Safeguards & Fallback Logic

### Minimum Sample Size Requirement

**Expert Recommendation:** Require ≥5,000-10,000 introns for robust fitting.

```python
def should_use_domain_adaptive(n_introns, user_mode, min_introns=5000):
    """
    Decide whether to use domain-adaptive or human-fixed normalization.

    Args:
        n_introns: Number of target species introns
        user_mode: User-specified mode ('auto', 'adaptive', 'human')
        min_introns: Minimum introns for robust fitting

    Returns:
        mode: 'domain_adaptive' or 'human_fixed'
        reason: Why this mode was chosen
    """
    if user_mode == 'human':
        return 'human_fixed', 'user_requested_human_mode'

    if user_mode == 'adaptive':
        if n_introns >= min_introns:
            return 'domain_adaptive', 'user_requested_adaptive_mode'
        else:
            return 'human_fixed', f'insufficient_introns_n={n_introns}_min={min_introns}'

    # user_mode == 'auto' (default)
    if n_introns >= min_introns:
        return 'domain_adaptive', f'auto_adaptive_n={n_introns}_sufficient'
    else:
        return 'human_fixed', f'auto_fallback_n={n_introns}_insufficient'
```

### CLI Flags

```bash
# Auto mode (default): Use domain-adaptive if ≥5k introns, else human-fixed
intronIC -g genome.fa -a annotation.gff3 -n species

# Force domain-adaptive (fails if insufficient introns)
intronIC -g genome.fa -a annotation.gff3 -n species --normalizer adaptive

# Force human-fixed (always use human normalizer)
intronIC -g genome.fa -a annotation.gff3 -n species --normalizer human

# Adjust minimum introns threshold
intronIC -g genome.fa -a annotation.gff3 -n species --normalizer auto --min-introns 10000

# Load pre-fitted normalizer (for exact reproducibility)
intronIC -g genome.fa -a annotation.gff3 -n species --normalizer-params celegans.normalizer.json
```

### Logging Example

```
ℹ Extracted 109,820 introns from C. elegans genome
ℹ Normalizer mode: auto (default)
ℹ Checking sample size for domain-adaptive normalization...
✓ Sample size sufficient (109,820 ≥ 5,000)
ℹ Using domain-adaptive normalization
ℹ Fitting normalizer on C. elegans data...
  - Quantile: 0.95 (from human CV)
  - Symmetric clipping: True
  - Samples: 109,820 total, 109,820 used (below 200k cap)
✓ Normalizer fitted (hash: a3f2c91b4d8e5f67)
  - Feature 1 (5'SS): scale=2.456, clip=[-4.12, +4.12]
  - Feature 2 (BP):   scale=1.823, clip=[-3.87, +3.87]
  - Feature 3 (3'SS): scale=3.012, clip=[-4.45, +4.45]
ℹ Classifying 109,820 introns with human-trained SVM...
✓ Classification complete
ℹ Saved normalizer: caenorhabditis_elegans.normalizer.json
```

**Expert note:** "Log which mode was used, with N_introns used for fitting, so it's explicit in metadata."

---

## Reproducibility

### Metadata Schema (Enhanced)

```json
{
  "model_name": "caenorhabditis_elegans.model",
  "intronIC_version": "2.1.0",
  "pipeline_architecture": "v2_domain_adaptive",

  "normalizer": {
    "mode": "domain_adaptive",
    "decision": {
      "user_mode": "auto",
      "n_introns": 109820,
      "min_introns_threshold": 5000,
      "chosen_mode": "domain_adaptive",
      "reason": "auto_adaptive_n=109820_sufficient"
    },
    "hyperparameters": {
      "quantile": 0.95,
      "symmetric": true
    },
    "fitted_parameters": {
      "scale": [2.456, 1.823, 3.012],
      "clip_lo": [-4.12, -3.87, -4.45],
      "clip_hi": [4.12, 3.87, 4.45],
      "fitted_on_species": "caenorhabditis_elegans",
      "fitted_on_n_samples_total": 109820,
      "fitted_on_n_samples_used": 109820,
      "fitted_on_timestamp": "2025-11-16T19:23:45Z",
      "normalizer_hash": "a3f2c91b4d8e5f67"
    }
  },

  "svm": {
    "source": "homo_sapiens_reference",
    "trained_on": {
      "u12_count": 387,
      "u2_count": 20690
    }
  },

  "reproducibility": {
    "command": "intronIC -g genome.fa -a annotation.gff3 -n caenorhabditis_elegans --normalizer auto --min-introns 5000 --random-seed 42",
    "normalizer_file": "caenorhabditis_elegans.normalizer.json",
    "exact_reproduction_command": "intronIC -g genome.fa -a annotation.gff3 -n caenorhabditis_elegans --normalizer-params caenorhabditis_elegans.normalizer.json --random-seed 42"
  }
}
```

### Normalizer Serialization

```json
// caenorhabditis_elegans.normalizer.json
{
  "quantile": 0.95,
  "symmetric": true,
  "scale": [2.456, 1.823, 3.012],
  "clip_lo": [-4.12, -3.87, -4.45],
  "clip_hi": [4.12, 3.87, 4.45],
  "fitted_on_species": "caenorhabditis_elegans",
  "fitted_on_n_samples_total": 109820,
  "fitted_on_n_samples_used": 109820,
  "fitted_on_timestamp": "2025-11-16T19:23:45Z",
  "normalizer_hash": "a3f2c91b4d8e5f67",
  "feature_names": ["5prime_ss", "branch_point", "3prime_ss"]
}
```

**To exactly reproduce:**
```bash
intronIC -g genome.fa -a annotation.gff3 -n caenorhabditis_elegans \
  --normalizer-params caenorhabditis_elegans.normalizer.json \
  --random-seed 42
```

This loads the exact fitted normalizer, bypassing the fitting step entirely.

---

## Why This Doesn't Reintroduce the Old Bug

### Old Bug (Double Scaling)
```
raw LLRs → human_scaler → z_human
        → StandardScaler.fit(z_human) → z_species
        → SVM (trained on z_species)
```

**Problems:**
1. **Two scaling steps** (opaque composition)
2. **Semantic zero destroyed** (z=0 no longer means "U12≈U2")
3. **Hard to reason about** what space SVM lives in

### New Approach (Single Domain-Conditional Scaling)
```
raw LLRs → species_scaler → z_species
        → species_clipper → z_clipped
        → SVM (trained on human z_clipped)
```

**Why it's safe:**
1. **One clear scaling step** per species (zero-anchored)
2. **Semantic zero preserved** (s=0 → z=0 always)
3. **Well-defined feature space**: "sign-preserving, per-species robust-normalized, moderately clipped LLR"
4. **Unsupervised normalization** decoupled from supervised classification
5. **Only SVM weights transfer** (interpretable)

**Key Difference:**
- Old: Fit second scaler on **human z-scores**, apply to target species
- New: Fit normalizer on **target species raw scores**, use with human SVM

---

## Implementation Phases

### Phase 1: Core DomainAdaptiveNormalizer (PRIORITY)
**Goal:** Replace ZeroAnchoredRobustScaler + SymmetricClipper with unified normalizer

- [ ] Implement `DomainAdaptiveNormalizer` class
  - [ ] Zero-anchored scaling (median(|s|))
  - [ ] Quantile-based clipping (species-conditional)
  - [ ] Symmetric vs asymmetric modes
  - [ ] Metadata serialization
- [ ] Add sklearn transformer interface (`fit`, `transform`, `fit_transform`)
- [ ] Add `from_metadata()` for loading pre-fitted normalizers
- [ ] Unit tests for normalizer

**Deliverable:** Standalone normalizer that can replace existing scale+clip steps

### Phase 2: Training Integration
**Goal:** Use domain-adaptive normalizer in training pipeline

- [ ] Update `optimizer.py`:
  - [ ] Replace separate `ZeroAnchoredRobustScaler` + `SymmetricClipper` with `DomainAdaptiveNormalizer`
  - [ ] Grid search over `quantile` (0.95, 0.975, 0.99)
  - [ ] Fit normalizer on human reference data
  - [ ] Save quantile hyperparameter to model
- [ ] Update `trainer.py`:
  - [ ] Use `DomainAdaptiveNormalizer` in ensemble training
  - [ ] Save normalizer metadata with model
- [ ] Retrain human model with new normalizer

**Deliverable:** Human model trained with domain-adaptive normalizer (but fitted on human)

### Phase 3: Deployment Integration
**Goal:** Fit normalizer on target species at deployment

- [ ] Update `predictor.py`:
  - [ ] Detect normalizer mode (auto/adaptive/human)
  - [ ] Check sample size threshold
  - [ ] Fit normalizer on target species if adaptive
  - [ ] Use human normalizer if fallback
  - [ ] Log decision and fitted parameters
- [ ] Add CLI flags:
  - [ ] `--normalizer {auto,adaptive,human}`
  - [ ] `--min-introns N`
  - [ ] `--normalizer-params FILE`
- [ ] Save normalizer metadata to output

**Deliverable:** Classification adapts normalizer to target species

### Phase 4: Reproducibility & Documentation
**Goal:** Ensure exact reproducibility and clear documentation

- [ ] Enhanced metadata:
  - [ ] Normalizer decision logic
  - [ ] Fitted parameters
  - [ ] Reproducibility commands
- [ ] Normalizer serialization:
  - [ ] Save fitted normalizer to JSON
  - [ ] Load pre-fitted normalizer for exact reproduction
- [ ] Documentation:
  - [ ] Update CLAUDE.md with domain-adaptive design
  - [ ] Add normalization guide
  - [ ] Update user manual
- [ ] Validation:
  - [ ] Test C. elegans with domain-adaptive normalizer
  - [ ] Compare human vs. adaptive modes
  - [ ] Verify exact reproducibility from metadata

**Deliverable:** Fully documented, reproducible domain-adaptive normalization

---

## Expected Results

### C. elegans (Before vs. After)

**Before (Human-Fixed Clipping):**
```
Raw z-scores: [-63, +16]
Clipped to: [-1.5, +1.5]  (human thresholds)
Result: 99.87% get 0% probability (catastrophic)
```

**After (Domain-Adaptive Normalizer):**
```
Raw LLRs: extreme due to 38% GC
Scaled: z ∈ [-3, +3] (C. elegans-specific scale)
Clipped: at 0.95 quantile of C. elegans z-scores (e.g., ±4.1)
Result: Diverse scores, SVM has discriminative power
```

### Human (Validation)

**Before (No Clipping):**
```
Raw LLRs → scale → SVM
Validation F_0.75: ? (currently training)
```

**After (Domain-Adaptive, fitted on human):**
```
Raw LLRs → scale (human) → clip (human 0.95 quantile) → SVM
Validation F_0.75: Should match or exceed previous clipped model
```

---

## Open Questions (RESOLVED)

### ✓ 1. Sample Size Thresholds
**Expert Answer:** ≥5,000-10,000 introns minimum for robust fitting
**Implementation:** Default 5,000, adjustable via `--min-introns`

### ✓ 2. Default Behavior
**Expert Answer:** Auto mode (adaptive if sufficient data, else fallback)
**Implementation:** `--normalizer auto` (default), with explicit override options

### ✓ 3. Scaler Fitting Strategy
**Expert Answer:** Fit on ALL introns (unsupervised), not just confident U12s
**Implementation:** Use entire target species intron set

### ✓ 4. Clipping Revisited
**Expert Answer:** YES, bring clipping back as domain-adaptive
**Implementation:** Quantile hyperparameter transfers, numeric thresholds don't

### ✓ 5. What About Confident U12s?
**Expert Answer:** Don't use them for normalization (introduces circularity)
**Implementation:** Normalizer is purely unsupervised, fits on all introns

### ✓ 6. Two-Pass vs Single-Pass?
**Expert Answer:** Single-pass (fit normalizer once on all data)
**Implementation:** No discovery pass needed, just fit on full dataset

---

## Design Principles (Updated)

### 1. Dimensionless Transfer
**Only hyperparameters transfer across species, never numeric thresholds**
- Quantile `q` transfers
- SVM weights transfer
- Feature architecture transfers
- Fitted scales/clips do NOT transfer

### 2. Unsupervised Normalization
**Normalizer fits on unlabeled target data**
- No labels required
- No confident U12 filtering
- Purely distribution-based
- Decoupled from classification

### 3. Reproducibility First
**Every run exactly recreatable from metadata**
- Save fitted normalizer parameters
- Document decision logic
- Generate exact reproduction commands
- Serialize normalizers for reuse

### 4. Conservative Defaults
**Fall back to human normalizer if uncertain**
- Require ≥5k introns for adaptation
- Auto mode with sensible fallback
- Human mode always available

### 5. Clear Communication
**Inform user of normalization decisions**
- Log mode selection and reasoning
- Show fitted parameters
- Explain fallback if triggered

---

## Success Criteria

### Minimum Viable Product
- [ ] `DomainAdaptiveNormalizer` implemented and tested
- [ ] Human model retrained with new normalizer
- [ ] C. elegans classification succeeds (diverse scores)
- [ ] Reproducibility verified (exact results from metadata)

### Full Release
- [ ] All phases complete
- [ ] Validated on 3+ species (C. elegans, D. melanogaster, mouse)
- [ ] Documentation complete
- [ ] Auto/adaptive/human modes all working
- [ ] User guide with recommendations

---

## Future Enhancements (Post-MVP)

### Expert Recommendations for Later Iterations

#### 1. Optional Saturating Transform After Clipping
**Expert note:** "Optional saturating transform after clipping if you still see occasional insane values in some species (e.g., sign(z) * log(1 + |z|))."

**Use case:** If even after domain-adaptive clipping, some species show extreme outliers that affect SVM performance.

**Implementation:**
```python
# After clipping, optional additional compression
if saturate_enabled:
    z_clipped = np.sign(z_clipped) * np.log1p(np.abs(z_clipped))
```

**Hyperparameter:** `saturate_enabled` (True/False), learned via CV on human

**Priority:** LOW (evaluate after domain-adaptive clipping is working)

---

#### 2. Species-Specific Decision Thresholds
**Expert note:** "Once you have some 'definitely no U12' species and some 'has U12' species, you might tune a slightly more conservative decision threshold globally."

**Use case:** If certain species show systematically higher/lower SVM scores due to composition bias that persists even after normalization.

**Implementation:**
- Establish "known U12-positive" species (human, mouse, fly)
- Establish "known U12-negative" species (if any exist)
- Compute species-specific calibration curve
- Adjust threshold based on species group

**Priority:** MEDIUM (wait until we have multi-species validation data)

---

#### 3. C. elegans Regression Test
**Expert note:** "Bake in C. elegans as a negative control in your test suite: Expect ≈0 U12 calls and also assert that the feature variance doesn't collapse after normalization."

**Implementation:**
```python
def test_celegans_normalization():
    """Regression test: C. elegans should not collapse after normalization."""
    # Load C. elegans data
    celegans_introns = load_celegans_data()
    X = extract_scores(celegans_introns)

    # Fit normalizer
    normalizer = DomainAdaptiveNormalizer(quantile=0.95)
    Z = normalizer.fit_transform(X, species_name='caenorhabditis_elegans')

    # Assert: Feature variance should NOT collapse
    # (old bug: all features squeezed to same value)
    for j in range(Z.shape[1]):
        var = np.var(Z[:, j])
        assert var > 0.1, f"Feature {j} variance collapsed: {var}"

    # Assert: Should get diverse predictions, not just 0% and 100%
    predictions = model.predict_proba(Z)[:, 1]
    unique_bins = len(np.histogram(predictions, bins=20)[0])
    assert unique_bins > 5, f"Predictions too concentrated: {unique_bins} bins"

    # Assert: Expect ≈0 U12 calls (C. elegans likely has no/few U12s)
    high_confidence_u12s = np.sum(predictions > 0.90)
    assert high_confidence_u12s < 100, f"Too many U12 calls: {high_confidence_u12s}"
```

**Priority:** HIGH (add immediately after MVP)

---

#### 4. Normalizer Caching Per Species
**Optimization:** Reuse fitted normalizers across runs on the same species.

**Implementation:**
```python
# First run on species X: fit and save
normalizer.fit(X, species_name='species_x')
normalizer.save('~/.intronIC/normalizers/species_x_v1.json')

# Subsequent runs: load cached normalizer
if cached_normalizer_exists('species_x'):
    normalizer = DomainAdaptiveNormalizer.from_file('~/.intronIC/normalizers/species_x_v1.json')
    # Skip fitting, use cached parameters
```

**Benefit:** Faster runs, consistent normalization across analyses

**Priority:** MEDIUM (nice-to-have after MVP)

---

#### 5. Per-Feature Scale Function Options
**Current:** Use `median(|s|)` for all features

**Future:** Allow different robust scale estimators per feature
- IQR (interquartile range)
- MAD (median absolute deviation)
- Trimmed standard deviation

**Hyperparameter:** `scale_method` per feature, learned via CV

**Priority:** LOW (median(|s|) works well, only change if needed)

---

## Related Documents

- **FP_ROOT_CAUSE_ANALYSIS.md** - C. elegans catastrophic failure analysis
- **SCALER_CENTERING_FIX_SUMMARY.md** - Original clipping bug
- **TRAINING_COMPLETE_SUMMARY.md** - Training with old clipping approach
- **ADAPTIVE_PIPELINE_PROPOSAL.md** - Original proposal (superseded by this)

---

## Expert Review Summary

**Reviewer Feedback (2025-11-16):** ✅ **APPROVED**

### Key Validation Points

1. **"q transfers, not numeric thresholds" - ✅ Exactly the right move**
   - Treats q as dimensionless hyperparameter
   - Fixes C. elegans [–63, 16] → human ±1.48 catastrophe
   - Per-feature quantiles confirmed

2. **"Unsupervised normalization on all introns" - ✅ Aligned**
   - Avoid circularity of using confident U12s
   - Robust estimators dominated by majority class (U2s)
   - Random sampling for large datasets (cap 100k-200k)

3. **"Unified DomainAdaptiveNormalizer" - ✅ Clean architecture**
   - Single responsibility: raw LLR → zero-anchored, scaled, clipped features
   - Serializable, reusable, unit-testable
   - Clear separation: fitted params vs hyperparameters

4. **"Dimensionless transfer principle" - ✅ Proper separation**
   - Hyperparameters (rules) transfer: q, saturate, architecture
   - Fitted quantities (numbers) don't transfer: scale_j, lo_j, hi_j
   - Correct domain-adaptation boundary

5. **"Auto mode with fallback" - ✅ Nice balance**
   - Auto default: adaptive if ≥5k, else human
   - Explicit override available
   - Log mode, N_introns, normalizer hash

### Additional Expert Recommendations (Incorporated)

✓ **Per-feature quantiles** - Already in design, now explicit
✓ **Random sampling cap** - Added 100k-200k limit
✓ **Enhanced logging** - Mode, N_introns, hash tracking
✓ **Future enhancements** - Saturating transform, species thresholds, C. elegans regression test

### Expert Conclusion

**"This is very much in line with what I'd recommend. It's a principled domain-adaptive preprocessor feeding a human-trained classifier, with clear separation between what's 'universal' and what's 'species-specific.'"**

---

**STATUS:** ✅ **EXPERT-APPROVED** - Ready for implementation. Design addresses root cause while maintaining reproducibility and avoiding old bugs.
