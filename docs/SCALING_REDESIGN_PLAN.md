# Scaling Architecture Redesign - Implementation Plan

## Expert Feedback Summary

**Problem:** Current double scaling (ZeroAnchoredRobustScaler → RobustScaler) is redundant and allows extreme cross-species outliers to overwhelm balance features.

**Solution:** Single zero-anchored scaler inside model pipeline, with aggressive outlier control in z-space.

---

## New Architecture

### Pipeline Design

```python
Pipeline([
    ('scale', ZeroAnchoredRobustScaler()),     # ONLY scaling step (fitted on human raw LLRs)
    ('clip', SymmetricClipper(zmax=...)),       # Clip in z-space (hyperparameter)
    ('saturate', SaturatingTransform()),        # Optional: log(1+|z|) compression
    ('augment', BothEndsStrongTransformer()),   # Balance features (unchanged)
    ('svc', LinearSVC(...))                     # Linear classifier
])
```

### Key Changes

1. **Single Scaler**
   - Move ZeroAnchoredRobustScaler FROM pre-processing INTO model pipeline
   - Remove RobustScaler (Stage 2) entirely
   - Fit on human raw LLRs (not z-scores)
   - `with_centering=False` preserves semantic zero

2. **Aggressive Outlier Clipping**
   - Clip in z-space (AFTER scaling, not before)
   - Symmetric: `[-zmax, +zmax]`
   - Hyperparameter: `estimator__clip__zmax: [3.5, 4.0, 5.0]`
   - OR quantile-based: `estimator__clip__quantile: [0.95, 0.975, 0.99]`

3. **Optional Saturating Transform**
   - Compresses extreme values: `f(z) = sign(z) * log(1 + |z|)`
   - Makes z=11.45 vs z=4 much closer in model space
   - Hyperparameter: `estimator__saturate__enabled: [true, false]`

4. **Precision-Focused Training**
   - Optimize F_β score with β < 1 (e.g., β=0.5)
   - Biases toward "only call U12 when very sure"
   - Conservative threshold for decision boundary

5. **Cross-Species Behavior**
   - **Default:** No per-species refitting (use human-trained pipeline as-is)
   - **Optional:** Advanced mode to refit scaler only on unlabeled data
   - Clipping + saturation prevent extreme outliers from dominating

---

## Implementation Checklist

### Phase 1: New Transformer Classes

- [ ] **Create `classification/clipping.py` updates:**
  - [ ] Rename `OutlierClipper` → `SymmetricClipper`
  - [ ] Change from quantile-based (on raw) to zmax-based (on z-scores)
  - [ ] Support both `zmax` (fixed) and `quantile` (fitted) modes
  - [ ] Update docstrings and examples

- [ ] **Create `classification/saturating.py`:**
  - [ ] New `SaturatingTransform` class
  - [ ] `__init__(enabled=True)`
  - [ ] `transform()`: `sign(z) * log(1 + |z|)` when enabled, identity otherwise
  - [ ] `get_feature_names_out()`: passthrough

### Phase 2: Pipeline Restructuring

- [ ] **Update `classification/optimizer.py`:**
  - [ ] Remove RobustScaler import
  - [ ] Move ZeroAnchoredRobustScaler into pipeline (first step)
  - [ ] Change OutlierClipper → SymmetricClipper (after scaler)
  - [ ] Add SaturatingTransform step
  - [ ] Update pipeline construction in `_run_cv_round()`
  - [ ] Update param_grid to include:
    - `estimator__clip__zmax: [3.5, 4.0, 5.0]`
    - `estimator__saturate__enabled: [true, false]`

- [ ] **Update `classification/trainer.py`:**
  - [ ] Same pipeline changes as optimizer
  - [ ] Ensure pipeline construction matches

- [ ] **Update `classification/predictor.py`:**
  - [ ] Change `_prepare_features()` to extract **raw LLRs** not z-scores
  - [ ] Remove z-score extraction logic
  - [ ] Extract: `[five_raw_score, bp_raw_score, three_raw_score]`
  - [ ] Pipeline will handle scaling internally

### Phase 3: Scoring Module Changes

- [ ] **Update `scoring/normalizer.py`:**
  - [ ] **Option A:** Keep ScoreNormalizer but make it ONLY store raw scores (no scaling)
  - [ ] **Option B:** Deprecate ScoreNormalizer entirely, let pipeline handle everything
  - [ ] Decision: **Option A** (less disruptive, maintains backward compatibility)
  - [ ] Remove scaling logic from `transform()`
  - [ ] Just pass through introns with raw scores populated

- [ ] **OR: Create `scoring/raw_scorer.py`:**
  - [ ] New simplified class that only computes raw LLRs
  - [ ] No scaling, no z-scores
  - [ ] Clean separation: scoring vs normalization

### Phase 4: Training Configuration

- [ ] **Update `config/training_default.yaml`:**
  ```yaml
  param_grid:
    # Feature augmentation
    estimator__augment__include_pairwise_mins: [false]
    estimator__augment__include_max: [false]

    # Outlier control (NEW)
    estimator__clip__zmax: [3.5, 4.0, 5.0]  # Symmetric clipping threshold

    # Saturating transform (NEW)
    estimator__saturate__enabled: [false, true]  # Optional compression

    # SVM hyperparameters
    estimator__svc__penalty: ['l2']
    estimator__svc__dual: [false]
    estimator__svc__C: [...]  # Geometric search as before

    # Calibration
    calibration_method: ['isotonic']

  # Scoring metric (NEW)
  scoring: 'f_beta'  # F_0.5 score (precision-focused)
  scoring_params:
    beta: 0.5  # Weight precision 2x more than recall
  ```

- [ ] **Update grid size comments:**
  - Old: 6 combinations per C value
  - New: 6 (base) × 3 (zmax) × 2 (saturate) = 36 per C
  - Or: Reduce C search points to compensate

### Phase 5: Model I/O

- [ ] **Update `utils/model_io.py`:**
  - [ ] Ensure pipeline with new steps serializes correctly
  - [ ] Add version metadata for new architecture
  - [ ] Backward compatibility: detect old models, warn user

- [ ] **Update `tools/inspect_model.py`:**
  - [ ] Handle new pipeline structure
  - [ ] Show scaler position in pipeline
  - [ ] Show clip/saturate parameters

### Phase 6: CLI and Main

- [ ] **Update `cli/main.py`:**
  - [ ] Ensure raw scores are computed before classification
  - [ ] Remove any Stage 1 scaling calls
  - [ ] Pipeline handles all scaling internally

- [ ] **Update `cli/args.py`:**
  - [ ] Add flag for optional domain adaptation mode:
    - `--adapt-scaler` (refit scaler on unlabeled target data)
  - [ ] Update help text

### Phase 7: Testing

- [ ] **Create `tests/test_new_scaling.py`:**
  - [ ] Test SymmetricClipper with zmax and quantile modes
  - [ ] Test SaturatingTransform (enabled/disabled)
  - [ ] Test full pipeline with raw LLR input
  - [ ] Verify z-scores never exceed zmax after clipping
  - [ ] Test cross-species scenario (human model on simulated outliers)

- [ ] **Update `test_transformer.py`:**
  - [ ] Ensure BothEndsStrongTransformer still works with new pipeline

### Phase 8: Documentation

- [ ] **Create `docs/SCALING_ARCHITECTURE_V2.md`:**
  - [ ] Document new single-scaler design
  - [ ] Explain outlier control strategy
  - [ ] Cross-species deployment guidelines
  - [ ] Migration guide from old models

- [ ] **Update `CLAUDE.md`:**
  - [ ] Replace double-scaling description with new architecture
  - [ ] Update pipeline diagrams
  - [ ] Update key file locations

- [ ] **Update `config/README.md`:**
  - [ ] Document new hyperparameters (zmax, saturate)
  - [ ] Explain precision-focused scoring

### Phase 9: Retraining and Validation

- [ ] **Retrain homo_sapiens model:**
  ```bash
  intronIC train -n homo_sapiens \
    --optimizer-config config/training_default.yaml \
    -p 12
  ```

- [ ] **Validate on C. elegans:**
  ```bash
  intronIC classify -g celegans.fa -a celegans.gff \
    -n caenorhabditis_elegans \
    --model homo_sapiens.model.pkl
  ```
  - [ ] Verify z-scores stay within reasonable range
  - [ ] Check false positive rate (should be much lower)
  - [ ] Inspect predictions with `tools/analyze_predictions.py`

- [ ] **Compare to old model:**
  - [ ] Same human test set: ensure similar performance
  - [ ] C. elegans: verify reduced false positives

---

## Implementation Order

### Sprint 1: Core Infrastructure (1-2 days)
1. Create SymmetricClipper class
2. Create SaturatingTransform class
3. Update optimizer.py pipeline
4. Update trainer.py pipeline
5. Basic testing

### Sprint 2: Integration (1-2 days)
6. Update predictor.py (raw LLRs)
7. Update normalizer.py (passthrough)
8. Update config files
9. Update CLI

### Sprint 3: Testing & Validation (1-2 days)
10. Write comprehensive tests
11. Retrain homo_sapiens model
12. Validate on C. elegans
13. Compare old vs new

### Sprint 4: Documentation & Polish (1 day)
14. Write architecture docs
15. Update CLAUDE.md
16. Migration guide
17. Code review and cleanup

---

## Open Design Questions

### Q1: Clipping Hyperparameter Strategy

**Option A:** Fixed zmax values
```yaml
estimator__clip__zmax: [3.5, 4.0, 5.0]
```
Pros: Simple, interpretable, species-invariant
Cons: Doesn't adapt to human training data distribution

**Option B:** Quantile-based on human data
```yaml
estimator__clip__quantile: [0.95, 0.975, 0.99]
```
Pros: Adapts to actual training distribution
Cons: Less interpretable, species-specific

**Option C:** Hybrid (fit quantile, then use fixed zmax)
```yaml
# During training: learn zmax from quantile on human data
# Store learned value in model
# During inference: use same fixed zmax
```
Pros: Best of both worlds
Cons: More complex

**Recommendation:** Start with **Option A** (fixed zmax), simplest and most interpretable.

### Q2: Saturating Transform - Always On or Optional?

**Option A:** Always enabled (no hyperparameter)
- Simpler pipeline
- One less hyperparameter to tune
- Always provides outlier protection

**Option B:** Grid search (enabled: [true, false])
- Lets model decide if compression helps
- More flexible
- Larger grid search space

**Recommendation:** **Option B** (grid search), let cross-validation decide.

### Q3: Precision-Focused Scoring Function

**Option A:** F_β score with β = 0.5
```python
from sklearn.metrics import fbeta_score, make_scorer
scorer = make_scorer(fbeta_score, beta=0.5)
```

**Option B:** Custom scorer (precision at fixed recall threshold)
```python
def precision_at_recall(y_true, y_pred, min_recall=0.90):
    # Find threshold that gives ≥90% recall
    # Return precision at that threshold
```

**Option C:** Weighted combination
```python
score = 0.7 * precision + 0.3 * recall
```

**Recommendation:** **Option A** (F_0.5 score), standard and interpretable.

### Q4: ScoreNormalizer Refactoring

**Option A:** Keep ScoreNormalizer, make it passthrough
- Minimal code changes
- Backward compatible
- Name becomes misleading (no normalization!)

**Option B:** Deprecate ScoreNormalizer, use RawScorer
- Cleaner separation of concerns
- More code changes
- Need migration path

**Option C:** Keep ScoreNormalizer name, but move ZAR scaler logic into it
- Still does normalization (but inside pipeline)
- Confusing dual purpose

**Recommendation:** **Option A** for now (keep as passthrough), rename to `RawScorer` in future refactor.

---

## Success Criteria

### Functional
- [ ] Pipeline accepts raw LLRs as input
- [ ] Single scaler fitted on human data
- [ ] Clipping prevents z-scores > zmax
- [ ] Human model performance maintained
- [ ] C. elegans false positives reduced by >50%

### Performance
- [ ] Training time similar or faster (fewer hyperparams if needed)
- [ ] Cross-validation runs complete successfully
- [ ] Memory usage unchanged

### Code Quality
- [ ] All tests pass
- [ ] Type hints on new classes
- [ ] Comprehensive docstrings
- [ ] No regression in existing functionality

### Documentation
- [ ] Architecture clearly explained
- [ ] Migration guide for users with old models
- [ ] Cross-species deployment documented

---

## Risk Mitigation

### Risk 1: Breaking Changes
**Impact:** Users with old models can't use new code
**Mitigation:**
- Detect model version in `model_io.py`
- Provide clear error message
- Document migration path

### Risk 2: Human Performance Degradation
**Impact:** New pipeline performs worse on human test set
**Mitigation:**
- Run extensive validation before deployment
- Keep old model architecture as fallback
- A/B test on human data first

### Risk 3: Grid Search Explosion
**Impact:** 36 combinations per C value is 6x slower
**Mitigation:**
- Reduce C search points (e.g., 3 rounds instead of 5)
- Use parallel CV processes
- Consider random search instead of grid

### Risk 4: Numerical Instability
**Impact:** log(1+|z|) for very small z might cause issues
**Mitigation:**
- Add epsilon: `log(1 + |z| + 1e-10)`
- Test edge cases (z=0, z→∞)
- Validate gradient stability

---

## Next Steps

**Immediate:**
1. User reviews this plan
2. Decide on open questions (Q1-Q4)
3. Begin Sprint 1 (SymmetricClipper + SaturatingTransform)

**After initial implementation:**
4. Run validation on human data
5. Test on C. elegans
6. Iterate based on results
7. Full documentation pass

**Ready to proceed?** Which sprint should we start with, or are there design questions to resolve first?
