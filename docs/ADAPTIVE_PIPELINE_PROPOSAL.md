# Adaptive Pipeline Proposal - Species-Specific Optimization

**Date:** 2025-11-16
**Status:** PROPOSAL - Awaiting Expert Review
**Context:** Post-clipping removal, designing robust cross-species adaptation

---

## Executive Summary

**Problem:** Cross-species composition bias causes extreme z-scores that challenge model generalization. Current solution (no clipping) trusts SVM to learn boundaries, but may not be optimal for all species.

**Proposal:** Implement optional species-adaptive pipeline that learns target-species statistics when sufficient data is available, while maintaining human-trained defaults for reproducibility and small datasets.

**Key Principle:** Reproducibility first. Any adaptive behavior must be fully documentable and recreatable.

---

## Background: The Clipping Failure

### What Happened
- **Old Architecture:** Human-trained clipping thresholds (±1.5-3σ)
- **C. elegans Result:** z-scores range [-63, +16] due to 38% GC vs 41% human
- **Outcome:** 99.87% of introns got 0% U12 probability (catastrophic failure)

### Root Cause
**Training-species thresholds don't transfer to target species with different nucleotide composition.**

### Current Solution (Implemented)
- **Removed clipping entirely** from pipeline
- Pipeline: `raw LLRs → scale → saturate → augment → svc`
- Trust SVM to learn appropriate boundaries from human training data
- **Trade-off:** Simpler, but may not be optimal for distant species

---

## Proposed Architecture: Two-Mode Pipeline

### Mode 1: Fixed Human-Based (DEFAULT)
**When:**
- Target species has insufficient confident U12s (below threshold)
- User explicitly requests via `--no-adapt` flag
- First-pass scoring (before any U12s identified)

**Pipeline:**
```
raw LLRs → ZeroAnchoredRobustScaler (human-fitted) → SaturatingTransform → BothEndsStrong → SVM
```

**Characteristics:**
- Uses human-trained PWMs and scaling parameters
- Fully reproducible across runs and species
- Conservative approach for novel species
- No adaptation to target species

**Metadata Captured:**
```yaml
mode: fixed_human_based
scaler_source: human_reference
pwm_source: human_u12_u2_matrices
adaptive: false
```

---

### Mode 2: Species-Adaptive (OPTIONAL)
**When:**
- Target species has ≥N confident U12s (threshold TBD, suggest N=20-50)
- User explicitly requests via `--adaptive` flag (opt-in)
- OR: Automatic if confident U12s exceed threshold AND `--auto-adapt` enabled

**Pipeline (Two-Pass):**
```
Pass 1: Fixed human-based (as above) → Identify confident U12s
Pass 2 (if threshold met):
  raw LLRs → ZeroAnchoredRobustScaler (SPECIES-fitted) → SaturatingTransform → BothEndsStrong → SVM
```

**Species-Adaptive Components:**

1. **Scaler Fitting:**
   - Fit `ZeroAnchoredRobustScaler` on target-species U12s + U2s
   - Uses species-specific median(|s|) for each feature
   - Still preserves semantic zero (s=0 → z=0)
   - **Result:** Z-scores reflect deviation from species norms, not human norms

2. **Optional: Species-Specific PWMs** (existing `--recursive` feature)
   - Build PWMs from target-species U12s
   - More accurate motif scoring for distant species

3. **Optional: SVM Retraining** (existing `--recursive` feature)
   - Retrain SVM on species-specific data
   - Better decision boundaries for species

**Characteristics:**
- Optimized for target species
- Requires sufficient confident U12s
- More complex, more parameters to document
- Higher computational cost (two passes)

**Metadata Captured (CRITICAL FOR REPRODUCIBILITY):**
```yaml
mode: species_adaptive
adaptive: true

# Pass 1 (discovery)
pass1:
  scaler_source: human_reference
  pwm_source: human_u12_u2_matrices
  u12_candidates_found: 87
  confident_u12s_threshold: 90.0
  confident_u12s_passing: 42

# Adaptation decision
adaptation_triggered: true
adaptation_reason: confident_u12s_above_threshold
adaptation_threshold: 20
adaptation_timestamp: "2025-11-16T12:34:56Z"

# Pass 2 (species-adapted)
pass2:
  scaler_source: species_fitted
  scaler_fit_data:
    u12_count: 42
    u2_count: 1000
    median_abs_5ss: 2.456
    median_abs_bp: 1.823
    median_abs_3ss: 3.012
  pwm_source: species_specific_matrices  # if --recursive
  pwm_fit_data:  # if --recursive
    u12_count: 42
    u2_count: 1000
  svm_retrained: true  # if --recursive
  svm_training_data:
    u12_count: 42
    u2_count: 1000

# Full reproducibility
random_seed: 42
intronIC_version: "2.0.1"
pipeline_architecture: "v2_no_clipping"
```

---

## Decision Logic: When to Adapt

### Proposed Thresholds (Expert Input Needed)

| Criterion | Threshold | Rationale |
|-----------|-----------|-----------|
| **Minimum confident U12s** | 20-50? | Need enough for robust statistics |
| **Confidence threshold** | 90%? | Same as default classification cutoff |
| **Minimum U2s for U2 scaler** | 500-1000? | Need U2 statistics for balanced scaling |

### Decision Tree

```
1. Run Pass 1 (human-based)
2. Count confident U12s (score ≥ threshold)
3. IF user specified --no-adapt:
     → Use Pass 1 results (DONE)
   ELIF confident_u12s < MIN_THRESHOLD:
     → Use Pass 1 results (insufficient data)
   ELIF user specified --adaptive OR --auto-adapt:
     → Run Pass 2 (species-adapted)
   ELSE:
     → Use Pass 1 results (DEFAULT: conservative)
```

### User Control Flags (Proposed)

```bash
# Conservative default (current behavior)
intronIC -g genome.fa -a annotation.gff3 -n species

# Explicit no adaptation (even if data sufficient)
intronIC -g genome.fa -a annotation.gff3 -n species --no-adapt

# Explicit adaptation request (if data sufficient)
intronIC -g genome.fa -a annotation.gff3 -n species --adaptive

# Automatic adaptation if threshold met
intronIC -g genome.fa -a annotation.gff3 -n species --auto-adapt

# Control adaptation threshold
intronIC -g genome.fa -a annotation.gff3 -n species --auto-adapt --adapt-threshold 30

# Full recursive (existing flag, implies --adaptive)
intronIC -g genome.fa -a annotation.gff3 -n species --recursive
```

---

## Reproducibility Requirements

### Critical: Every Run Must Be Fully Recreatable

1. **Metadata File (JSON/YAML)**
   - Capture ALL parameters affecting pipeline behavior
   - Include data-dependent decisions (adaptation triggered? why?)
   - Include fitted scaler parameters (if species-adapted)
   - Include random seeds, versions, timestamps

2. **Model Metadata Enhancement**
   - Current `homo_sapiens.model.metadata.json` needs expansion
   - Add `adaptation_metadata` section
   - Add `scaler_parameters` section
   - Add `decision_log` section

3. **Reproducibility Command**
   - Generate exact command to recreate run
   - Example:
   ```bash
   # This run used species-adaptive scaling
   # To reproduce exactly:
   intronIC -g genome.fa -a annotation.gff3 -n species \
     --adaptive \
     --random-seed 42 \
     --threshold 90.0 \
     --scaler-params scaler_params.json  # Load exact fitted scaler
   ```

4. **Scaler Serialization**
   - If species-adapted, save fitted scaler separately
   - Allow loading pre-fitted scaler via `--scaler-params`
   - Ensures exact reproduction even if target introns change

---

## Implementation Phases

### Phase 1: Enhanced Metadata (PREREQUISITE)
**Goal:** Ensure all runs are fully documented

- [ ] Expand model metadata structure
- [ ] Add adaptation decision logging
- [ ] Add scaler parameter capture
- [ ] Add reproducibility command generation
- [ ] Update metadata schema version

**Deliverable:** Every run (adaptive or not) has complete provenance

### Phase 2: Species-Adaptive Scaler (CORE FEATURE)
**Goal:** Fit scaler on target species when data sufficient

- [ ] Implement species-scaler fitting logic
- [ ] Add confident U12 counting/filtering
- [ ] Add threshold-based decision logic
- [ ] Add `--adaptive`, `--no-adapt`, `--auto-adapt` flags
- [ ] Add `--adapt-threshold N` parameter
- [ ] Serialize fitted scaler to JSON
- [ ] Update predictor to use species-fitted scaler

**Deliverable:** Two-pass pipeline with species-adapted scaling

### Phase 3: Integration with Recursive Training (OPTIONAL)
**Goal:** Combine adaptive scaler with existing recursive features

- [ ] Integrate with `--recursive` flag
- [ ] Ensure PWM building uses same U12/U2 sets as scaler
- [ ] Ensure SVM retraining uses same U12/U2 sets as scaler
- [ ] Add unified "full adaptation" mode

**Deliverable:** Complete species-optimization pipeline

### Phase 4: Validation & Documentation (CRITICAL)
**Goal:** Verify adaptive pipeline works correctly

- [ ] Test on C. elegans (adaptive vs. fixed)
- [ ] Test on D. melanogaster
- [ ] Test on species with varying U12 counts
- [ ] Verify reproducibility (exact results from metadata)
- [ ] Document adaptive mode in CLAUDE.md
- [ ] Update user guide with adaptation recommendations

**Deliverable:** Validated, documented adaptive pipeline

---

## Open Questions for Expert Review

### 1. Sample Size Thresholds
**Question:** What are appropriate minimum thresholds for adaptation?
- Minimum confident U12s for scaler fitting?
- Minimum U2s for balanced scaler?
- Should thresholds differ for scaler vs. PWM vs. SVM?

**Current Guess:** 20-50 U12s, 500-1000 U2s

---

### 2. Default Behavior
**Question:** Should adaptation be opt-in or opt-out?

**Option A (Conservative - RECOMMENDED):**
- Default: Fixed human-based (current)
- User must request `--adaptive` or `--auto-adapt`
- Pros: Reproducibility, predictability, safer
- Cons: Users may not know to use adaptive mode

**Option B (Aggressive):**
- Default: Auto-adapt if threshold met
- User must request `--no-adapt` to prevent
- Pros: Automatic optimization
- Cons: Different behavior across runs, harder to reproduce

**Recommendation:** Option A (opt-in) for initial release, reconsider after validation.

---

### 3. Scaler Fitting Strategy
**Question:** How to fit species-scaler when U12s/U2s are imbalanced?

**Option A (Current Plan):**
- Fit on discovered U12s + subset of confident U2s
- Use same U2 subsampling as SVM training (balanced)

**Option B (All Data):**
- Fit on all discovered U12s + all U2s
- Reflects true species distribution (imbalanced)

**Option C (Hybrid):**
- Fit 5' and 3' scalers on all data (abundant signal)
- Fit BP scaler on balanced subset (sparser signal)

**Recommendation:** Need expert input on whether scaler should reflect balanced or imbalanced distribution.

---

### 4. Clipping Revisited
**Question:** Should adaptive mode include species-fitted clipping?

**Consideration:** If we fit scaler on target species, clipping thresholds would also be species-specific (e.g., 0.95 quantile of target species z-scores, not human z-scores).

**Option A (No Clipping - CURRENT):**
- Trust SVM boundaries, even in adaptive mode
- Simpler, fewer parameters
- Consistent with fixed mode

**Option B (Species-Adaptive Clipping):**
- Clip at 0.95-0.99 quantile of species z-scores
- Protects against species-specific extremes
- More parameters to tune and document

**Recommendation:** Start without clipping (Option A), revisit if extreme values cause issues in adaptive mode.

---

### 5. Validation Strategy
**Question:** How to validate adaptive mode is actually better?

**Metrics to Compare (Fixed vs. Adaptive):**
1. **C. elegans:** Do we get intermediate scores (not just 0% and 100%)?
2. **Known U12s:** Recall on validated U12 introns
3. **False Positive Rate:** Precision on high-confidence predictions
4. **Score Distribution:** More reasonable spread vs. bimodal
5. **Cross-Species Consistency:** Do orthologous U12s score similarly?

**Test Species:**
- C. elegans (38% GC, distant from human)
- D. melanogaster (42% GC, intermediate)
- Mouse (42% GC, close to human)

**Recommendation:** Need curated U12 sets for each species to measure recall/precision.

---

### 6. Computational Cost
**Question:** Is two-pass acceptable for large genomes?

**Estimated Cost (Human Genome):**
- Pass 1: ~30 minutes (current)
- Pass 2: ~30 minutes (refit scaler + rescore)
- Total: ~60 minutes (2× current)

**Mitigation Options:**
- Only rescore high-confidence candidates in Pass 2 (not all introns)
- Cache Pass 1 results, allow user to restart from checkpoint
- Parallelize Pass 2 more aggressively

**Recommendation:** Accept 2× cost for adaptive mode (opt-in), optimize later if needed.

---

## Design Principles

### 1. Reproducibility First
- **Every run must be exactly recreatable from metadata**
- Capture all parameters, data-dependent decisions, fitted values
- Generate reproducibility command automatically
- Version metadata schema

### 2. Conservative Defaults
- **Default behavior: Fixed human-based (current)**
- Adaptation is opt-in (`--adaptive` or `--auto-adapt`)
- Prevents surprises, maintains backward compatibility

### 3. Clear User Communication
- **Inform user when adaptation triggers**
- Log adaptation decision and rationale
- Provide statistics on adaptation quality
- Example:
  ```
  ℹ Found 42 confident U12s (threshold: 20)
  ℹ Triggering species-adaptive scaling
  ℹ Refitting scaler on species data...
  ✓ Species scaler fitted (42 U12, 1000 U2)
  ℹ Re-scoring all introns with species-adapted pipeline...
  ```

### 4. Graceful Degradation
- **If adaptation fails, fall back to fixed mode**
- Log failure reason
- Continue with human-based results
- Don't fail entire run

### 5. Expert Override
- **Allow experts to control every parameter**
- `--adapt-threshold N`
- `--scaler-params file.json` (load pre-fitted)
- `--no-adapt` (force fixed mode)
- `--adaptive` (force adaptive mode)

---

## Success Criteria

### Minimum Viable Product (MVP)
- [x] Remove clipping (DONE)
- [ ] Enhanced metadata captures full provenance
- [ ] Species-adaptive scaler fitting works
- [ ] Adaptive mode improves C. elegans results
- [ ] Reproducibility verified (exact results from metadata)

### Full Release
- [ ] All phases implemented
- [ ] Validated on 3+ species
- [ ] Documentation complete
- [ ] Expert review incorporated
- [ ] User guide with recommendations

---

## Risk Analysis

### Risk 1: Over-Fitting to Target Species
**Risk:** Adaptive scaler fits noise in small U12 set
**Mitigation:** Require minimum threshold (20-50 U12s)
**Severity:** Medium

### Risk 2: Reproducibility Failure
**Risk:** Cannot recreate exact results from metadata
**Mitigation:** Comprehensive metadata, validation tests
**Severity:** HIGH - Core requirement

### Risk 3: User Confusion
**Risk:** Users don't understand when/why to use adaptive mode
**Mitigation:** Clear documentation, sensible defaults
**Severity:** Medium

### Risk 4: Increased Complexity
**Risk:** More modes = more bugs, harder to maintain
**Mitigation:** Thorough testing, modular design
**Severity:** Medium

### Risk 5: Computational Cost
**Risk:** 2× runtime unacceptable for some users
**Mitigation:** Opt-in only, optimize Pass 2
**Severity:** Low - Only affects opt-in users

---

## Next Steps

1. **Expert Review** (CURRENT STEP)
   - Review open questions
   - Provide recommendations on thresholds, defaults, fitting strategy
   - Identify any missing considerations

2. **Finalize Design**
   - Incorporate expert feedback
   - Update proposal with decisions
   - Create detailed implementation plan

3. **Implementation**
   - Phase 1: Enhanced metadata
   - Phase 2: Adaptive scaler
   - Phase 3: Integration
   - Phase 4: Validation

4. **Validation & Documentation**
   - Test on multiple species
   - Verify reproducibility
   - Document adaptive mode
   - Update user guide

---

## Related Documents

- **SCALER_CENTERING_FIX_SUMMARY.md** - Original clipping bug fix
- **TRAINING_COMPLETE_SUMMARY.md** - Training with clipping
- **FP_ROOT_CAUSE_ANALYSIS.md** - C. elegans failure analysis
- **config/training_default.yaml** - Current configuration (no clipping)

---

## Appendix: Example Metadata Files

### Example 1: Fixed Human-Based Run
```json
{
  "model_name": "caenorhabditis_elegans.model",
  "training_date": "2025-11-16 18:45:23",
  "intronIC_version": "2.0.1",
  "pipeline_architecture": "v2_no_clipping",

  "mode": "fixed_human_based",
  "adaptive": false,

  "scaler_source": "human_reference",
  "scaler_parameters": {
    "fitted_on": "homo_sapiens_u12_u2_reference",
    "median_abs_5ss": 1.234,
    "median_abs_bp": 0.987,
    "median_abs_3ss": 1.456
  },

  "pwm_source": "human_u12_u2_matrices",

  "reference_data": {
    "u12": {
      "file": "data/u12_reference.introns.iic.gz",
      "count": 387,
      "species": "homo_sapiens"
    },
    "u2": {
      "file": "data/u2_reference.introns.iic.gz",
      "count": 20690,
      "species": "homo_sapiens"
    }
  },

  "reproducibility_command": "intronIC -g genome.fa -a annotation.gff3 -n caenorhabditis_elegans --threshold 90.0 --random-seed 42"
}
```

### Example 2: Species-Adaptive Run
```json
{
  "model_name": "caenorhabditis_elegans.model",
  "training_date": "2025-11-16 18:45:23",
  "intronIC_version": "2.0.1",
  "pipeline_architecture": "v2_no_clipping",

  "mode": "species_adaptive",
  "adaptive": true,
  "adaptation_triggered_by": "auto_adapt_flag",

  "pass1": {
    "scaler_source": "human_reference",
    "pwm_source": "human_u12_u2_matrices",
    "introns_scored": 109820,
    "u12_candidates_found": 87,
    "confident_u12s": {
      "threshold": 90.0,
      "count": 42
    }
  },

  "adaptation_decision": {
    "triggered": true,
    "reason": "confident_u12s_above_threshold",
    "threshold": 20,
    "confident_u12s": 42,
    "timestamp": "2025-11-16T18:52:10Z"
  },

  "pass2": {
    "scaler_source": "species_fitted",
    "scaler_fit_data": {
      "u12_count": 42,
      "u12_source": "pass1_confident_predictions",
      "u2_count": 1000,
      "u2_source": "pass1_confident_u2_subset",
      "median_abs_5ss": 3.456,
      "median_abs_bp": 2.123,
      "median_abs_3ss": 4.567,
      "fit_timestamp": "2025-11-16T18:52:11Z"
    },
    "pwm_source": "human_u12_u2_matrices",
    "introns_rescored": 109820
  },

  "reference_data": {
    "u12": {
      "file": "data/u12_reference.introns.iic.gz",
      "count": 387,
      "species": "homo_sapiens"
    },
    "u2": {
      "file": "data/u2_reference.introns.iic.gz",
      "count": 20690,
      "species": "homo_sapiens"
    }
  },

  "reproducibility_command": "intronIC -g genome.fa -a annotation.gff3 -n caenorhabditis_elegans --auto-adapt --adapt-threshold 20 --threshold 90.0 --random-seed 42 --scaler-params caenorhabditis_elegans.scaler.json"
}
```

---

**END OF PROPOSAL - AWAITING EXPERT INPUT**
