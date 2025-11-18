# Scaler Serialization Implementation Plan

**Date:** 2025-11-17
**Objective:** Fix cross-species classification by serializing the human-trained scaler
**Root Cause:** Domain adaptation refits scaler on target species, breaking composition bias correction

---

## Problem Summary

**Current Behavior:**
- Training: Fit scaler on human reference → Save ensemble only (no scaler)
- Classify: Load ensemble → **Refit NEW scaler on experimental data** ("domain adaptation")
- Result: C. elegans gets 130 FPs (composition bias not corrected)

**Expert's Insight:**
- Domain adaptation assumption (99.5% U2 → robust baseline) is statistically sound
- BUT: Combining adapted features + human decision rule → FPs in U12-free species
- Fix: Serialize human scaler, make "human mode" default, keep adaptive as opt-in

---

## Implementation Steps

### 1. Save Normalizer with Model

**File:** `cli/main.py:1467-1475`

**Current:**
```python
model_bundle = {
    'ensemble': result.ensemble,
    'threshold': config.scoring.threshold
}
joblib.dump(model_bundle, model_path, compress=3)
```

**Change to:**
```python
model_bundle = {
    'ensemble': result.ensemble,
    'normalizer': normalizer,  # ADD: Save human-trained scaler
    'threshold': config.scoring.threshold
}
joblib.dump(model_bundle, model_path, compress=3)
```

**Note:** `normalizer` is available at line 1467 (returned from `classify_introns()` at line 1380)

---

### 2. Add CLI Flag for Normalizer Mode

**File:** `cli/args.py` (classify arguments section, ~line 345-360)

**Add new argument:**
```python
classify_parser.add_argument(
    '--normalizer_mode', '--normalizer-mode',
    choices=['human', 'adaptive', 'auto'],
    default='auto',
    help='''Normalizer mode for pretrained model classification:
      human: Use scaler from training species (recommended for U12-absent genomes)
      adaptive: Refit scaler on experimental data (experimental)
      auto: Use human if available, fall back to adaptive (default)'''
)
```

**File:** `cli/config.py` (ScoringConfig, ~line 70-90)

**Add field:**
```python
@dataclass
class ScoringConfig:
    # ... existing fields ...
    normalizer_mode: str = 'auto'  # 'human', 'adaptive', or 'auto'
```

---

### 3. Update Pretrained Classification Logic

**File:** `cli/main.py:1053-1128` (`classify_with_pretrained_model()`)

**Current (lines 1085-1110):**
```python
if isinstance(model_data, dict):
    ensemble = model_data['ensemble']
    saved_threshold = model_data.get('threshold', config.scoring.threshold)
    messenger.log_only("Loaded model bundle (dict format)")
else:
    ensemble = model_data
    saved_threshold = config.scoring.threshold
    messenger.log_only("Loaded model ensemble (legacy format)")

# Domain adaptation (THE PROBLEM!)
messenger.info("Fitting normalizer on experimental data (domain adaptation)")
normalizer = ScoreNormalizer()
normalizer.fit(introns, dataset_type='unlabeled')
```

**Change to:**
```python
# Load model components
if isinstance(model_data, dict):
    ensemble = model_data['ensemble']
    saved_normalizer = model_data.get('normalizer', None)
    saved_threshold = model_data.get('threshold', config.scoring.threshold)
    messenger.log_only("Loaded model bundle (dict format)")
else:
    ensemble = model_data
    saved_normalizer = None
    saved_threshold = config.scoring.threshold
    messenger.log_only("Loaded model ensemble (legacy format)")

# Determine normalizer mode
normalizer_mode = config.scoring.normalizer_mode
if normalizer_mode == 'auto':
    # Use human scaler if available, otherwise adaptive
    if saved_normalizer is not None:
        normalizer_mode = 'human'
        messenger.log_only("Auto mode: Using saved human scaler (recommended)")
    else:
        normalizer_mode = 'adaptive'
        messenger.log_only("Auto mode: Falling back to adaptive (no saved scaler)")

# Apply normalizer strategy
if normalizer_mode == 'human':
    if saved_normalizer is None:
        raise ValueError(
            "Normalizer mode 'human' requested but model has no saved scaler. "
            "Retrain model or use '--normalizer_mode adaptive'."
        )
    messenger.info("Using human-trained normalizer (scaler from training species)")
    messenger.log_only("This preserves composition bias correction across species")
    normalizer = saved_normalizer
else:  # adaptive
    messenger.info("Fitting normalizer on experimental data (domain adaptation)")
    messenger.log_only("Note: May cause FPs in U12-absent species (e.g., C. elegans)")
    normalizer = ScoreNormalizer()
    normalizer.fit(introns, dataset_type='unlabeled')
```

---

### 4. Update Model Metadata

**File:** `cli/main.py:1477-1496` (`generate_training_metadata()` call)

**Add to metadata:**
```python
metadata = generate_training_metadata(
    # ... existing params ...
    has_saved_normalizer=True,  # NEW: Indicate scaler is serialized
    normalizer_type='RobustScaler',  # NEW: Document scaler type
    # ... rest of params ...
)
```

---

### 5. Update Documentation/Logging

**Training completion message:**
```python
messenger.success(f"Model saved successfully with human-trained scaler")
messenger.log_only("Use '--normalizer_mode human' for cross-species classification")
```

**Classification mode selection:**
```python
if normalizer_mode == 'human':
    messenger.log_only("Normalizer: Using saved human scaler (recommended)")
elif normalizer_mode == 'adaptive':
    messenger.warning("Normalizer: Domain adaptive mode (may cause FPs in U12-absent species)")
```

---

## Testing Plan

### 1. Retrain Model with Scaler Serialization

```bash
pixi run intronIC train -n homo_sapiens_with_scaler -o . --n_optimization_rounds 3 --n_cv_folds 7 -p 10
```

**Verify:**
- Model file size increases (~few KB for scaler)
- Metadata shows `has_saved_normalizer: true`

### 2. Test Human Mode on C. elegans

```bash
pixi run intronIC classify \
  -g celegans_genome.fa.gz \
  -a celegans_annotation.gff3.gz \
  -n celegans_human_mode \
  --model homo_sapiens_with_scaler.model.pkl \
  --normalizer_mode human \
  -p 10
```

**Expected Result:** 0-1 false positives (vs current 130)

### 3. Test Adaptive Mode (Baseline)

```bash
pixi run intronIC classify \
  -g celegans_genome.fa.gz \
  -a celegans_annotation.gff3.gz \
  -n celegans_adaptive_mode \
  --model homo_sapiens_with_scaler.model.pkl \
  --normalizer_mode adaptive \
  -p 10
```

**Expected Result:** 130 false positives (same as before - confirms no regression)

### 4. Test Auto Mode

```bash
pixi run intronIC classify \
  -g celegans_genome.fa.gz \
  -a celegans_annotation.gff3.gz \
  -n celegans_auto_mode \
  --model homo_sapiens_with_scaler.model.pkl \
  # No --normalizer_mode flag (defaults to auto)
  -p 10
```

**Expected Result:** Should use human mode automatically → 0-1 FPs

### 5. Test Legacy Model Compatibility

```bash
pixi run intronIC classify \
  -g celegans_genome.fa.gz \
  -a celegans_annotation.gff3.gz \
  -n celegans_legacy_test \
  --model homo_sapiens_phase1.model.pkl \  # Old model without scaler
  -p 10
```

**Expected Result:** Auto mode falls back to adaptive → 130 FPs (with warning in log)

---

## Benefits of This Approach

### 1. **Fixes Cross-Species Classification**
- Maintains composition bias correction across species
- Expected: 0-1 FPs on C. elegans (U12-free genome)
- Preserves calibration consistency (p ≥ 0.90 means same thing across species)

### 2. **Backward Compatible**
- Old models without scaler still work (auto mode falls back to adaptive)
- Legacy behavior preserved with `--normalizer_mode adaptive`

### 3. **Flexible & Future-Proof**
- Can experiment with domain adaptation later
- Can add new normalizer strategies (e.g., marginal alignment)
- Clear separation of concerns: discrimination (SVM) vs preprocessing (scaler)

### 4. **Well-Documented**
- Clear logging messages explain mode selection
- Warnings for potentially problematic modes
- Metadata tracks scaler availability

---

## Code Locations Summary

| Task | File | Lines | Description |
|------|------|-------|-------------|
| Save normalizer | `cli/main.py` | 1467-1475 | Add to model_bundle |
| Load normalizer | `cli/main.py` | 1085-1110 | Extract from model_bundle |
| Mode selection | `cli/main.py` | 1085-1128 | Implement human/adaptive/auto logic |
| CLI flag | `cli/args.py` | ~345-360 | Add --normalizer_mode argument |
| Config field | `cli/config.py` | ~70-90 | Add normalizer_mode to ScoringConfig |
| Metadata | `cli/main.py` | 1477-1496 | Document scaler in training metadata |

---

## Next Steps

1. ✅ Review this plan
2. ⬜ Implement changes to code files
3. ⬜ Retrain model with scaler serialization
4. ⬜ Test on C. elegans (human mode)
5. ⬜ Validate results (expect 0-1 FPs)
6. ⬜ Document final architecture
7. ⬜ Update user documentation

---

**Status:** ⬜ Ready to implement
**Estimated Time:** 1-2 hours (coding + testing)
**Risk:** Low (backward compatible, well-isolated changes)
