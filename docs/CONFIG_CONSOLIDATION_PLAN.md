# Configuration File Consolidation Plan

**Date:** 2025-11-18
**Status:** Proposal - awaiting implementation

---

## Problem

Two configuration files with overlapping parameters:
- `config/default.toml` - General user-facing config (263 lines)
- `config/training_default.yaml` - Advanced training config (220+ lines)

**Overlapping parameters:**
- `n_models`: 1 (TOML) vs 12 (YAML)
- `max_iter`: 50000 (TOML) vs 75000 (YAML optimizer), 50000 (YAML training)
- `n_cv_folds`: 5 (TOML) vs 10 (YAML)
- `n_optimization_rounds`: 5 (TOML) vs 5 (YAML)
- `random_seed/random_state`: 42 (both)

**Result:** Confusion about which file takes precedence, duplication of settings.

---

## Proposed Solution

### Format Choice: YAML ✅

**Rationale:**
- Industry standard for Python ML (PyTorch, TensorFlow, sklearn, Hugging Face)
- Better support for nested structures and lists
- Native Python support via PyYAML (already in dependencies)
- More readable for complex configurations
- Easier to add defaults/overrides with anchors

**Migration:** TOML → YAML consolidation

---

## Unified Config Structure

**Single file:** `config/config.yaml` (with profiles: `quick.yaml`, `production.yaml`)

```yaml
# ============================================================================
# intronIC Configuration
# ============================================================================
# Default configuration for intronIC (both classify and train modes)
# CLI arguments override these settings
#
# Auto-loaded from (in order of precedence):
#   1. --config PATH (explicit CLI argument)
#   2. ./.intronIC.yaml (current directory)
#   3. ~/.config/intronIC/config.yaml (XDG config dir)
#   4. ~/.intronIC.yaml (user home)
#   5. <install_dir>/config/config.yaml (built-in defaults)
#
# Generate template: intronIC --generate-config
# ============================================================================

# ============================================================================
# SCORING PARAMETERS
# ============================================================================
scoring:
  # Classification threshold (0-100%)
  threshold: 90.0

  # Feature type: "cds", "exon", or "both"
  feature_type: "both"

  # Exclude non-canonical introns (GT-AG, GC-AG, AT-AC only)
  exclude_noncanonical: false

  # PWM pseudocount (log(0) prevention)
  pseudocount: 0.0001

  # Ignore terminal dinucleotides for NC introns
  ignore_nc_dinucleotides: true

  # Enable U12 boundary correction for NC introns
  u12_boundary_correction: true

  # Scoring regions (PWM application windows)
  regions:
    five_prime:
      start: -3
      end: 9
    branch_point:
      start: -55
      end: -5
    three_prime:
      start: -6
      end: 4

# ============================================================================
# TRAINING PARAMETERS
# ============================================================================
training:
  # Evaluation mode: "nested_cv", "split", or "none"
  eval_mode: "nested_cv"

  # Nested CV folds (only if eval_mode=nested_cv)
  n_cv_folds: 7

  # Test fraction (only if eval_mode=split)
  test_fraction: 0.2

  # Ensemble configuration
  ensemble:
    # Number of models (1=single, 3-10=ensemble)
    n_models: 10

    # Subsample U2 for diversity (only if n_models > 1)
    subsample_u2: true

    # U2 subsample fraction per model
    subsample_ratio: 0.85

    # Max iterations for LinearSVC
    max_iter: 40000

    # Random seed for reproducibility
    random_state: 42

  # Use fold-averaged hyperparameters from nested CV
  # - false: Re-optimize on full dataset (better training fit)
  # - true: Use geometric mean of fold C values (better cross-species generalization)
  use_fold_averaged_params: false

  # Fixed SVM C parameter (skip optimization if set)
  # fixed_C: 50.0  # Uncomment to use fixed value

# ============================================================================
# HYPERPARAMETER OPTIMIZATION
# ============================================================================
optimizer:
  # Number of optimization rounds (geometric refinement)
  n_rounds: 5

  # Initial grid points (round 1)
  n_points_initial: 13

  # Refinement grid points (rounds 2+)
  n_points_refine: 100

  # Cross-validation folds for parameter search
  cv_folds: 9

  # Parallel jobs (-1 = all cores)
  n_jobs: -1

  # Max iterations for LinearSVC convergence
  max_iter: 60000

  # Random seed
  random_state: 42

  # Verbose output
  verbose: true

  # C parameter bounds (effective positive class range)
  c_bounds:
    eff_C_pos_range: [0.001, 1000.0]
    eff_C_neg_max: null  # Auto-calculate from class weights

# ============================================================================
# PARAMETER GRID (Advanced - for custom search space)
# ============================================================================
param_grid:
  # Calibration method
  estimator__method:
    - sigmoid
    - isotonic

  # SaturatingTransform (optional log compression)
  estimator__svc__augment__saturate_enabled:
    - false

# ============================================================================
# PERFORMANCE
# ============================================================================
performance:
  # Parallel processes for classification
  processes: 1

  # Parallel processes for CV (if different)
  # cv_processes: 8  # Uncomment to override

  # Minimum intron length (bp)
  min_intron_length: 30

# ============================================================================
# ISOFORM SELECTION
# ============================================================================
isoform:
  # Include non-longest isoforms
  allow_multiple_isoforms: false

  # Exclude overlapping introns
  exclude_overlapping: false

  # Include duplicates in output
  include_duplicates: false

# ============================================================================
# OUTPUT OPTIONS
# ============================================================================
output:
  # Remove "transcript:" and "gene:" prefixes
  clean_names: true

# ============================================================================
# EXTRACTION OPTIONS
# ============================================================================
extraction:
  # Exonic flank length (bp)
  flank_length: 50

# ============================================================================
# ADVANCED OPTIONS
# ============================================================================
advanced:
  # Global random seed
  random_seed: 42

  # Suppress non-essential output
  quiet: false

  # Enable debug logging
  debug: false
```

---

## Migration Strategy

### Phase 1: Create Unified Config ✅
1. Create `config/config.yaml` with merged structure
2. Keep all parameters from both files
3. Resolve conflicts (use YAML values for training parameters)

### Phase 2: Create Config Profiles
1. `config/profiles/quick.yaml` - Fast testing (1-2 rounds, 3 folds)
2. `config/profiles/production.yaml` - Full production (5 rounds, 10 folds, 10 models)
3. `config/profiles/minimal.yaml` - Single model, no ensemble

### Phase 3: Update Config Loading Code
1. Modify `cli/config_loader.py` to handle YAML hierarchy
2. Add auto-loading from standard paths (in order):
   - CLI `--config` argument (highest priority)
   - `./.intronIC.yaml` (project-specific)
   - `~/.config/intronIC/config.yaml` (user config)
   - `~/.intronIC.yaml` (user home)
   - `<install>/config/config.yaml` (built-in defaults)
3. Ensure CLI args override config values
4. Add `--generate-config` command to create template

### Phase 4: Update Code References
1. Search for all references to `training_default.yaml`
2. Update to use unified `config.yaml`
3. Update `--optimizer-config` to `--config` (more general)
4. Keep backward compatibility (deprecation warning)

### Phase 5: Deprecate Old Files
1. Move `default.toml` → `config/deprecated/default.toml`
2. Move `training_default.yaml` → `config/profiles/production.yaml` (repurposed)
3. Add deprecation notices if old files are found

### Phase 6: Documentation
1. Update README with new config structure
2. Add config documentation: `docs/CONFIGURATION.md`
3. Update CLI help text
4. Add migration guide for existing users

---

## Auto-Loading Priority

```
Priority (highest to lowest):
1. CLI arguments (--threshold, --n-models, etc.)
2. --config FILE (explicit config path)
3. ./.intronIC.yaml (project directory)
4. ~/.config/intronIC/config.yaml (XDG config)
5. ~/.intronIC.yaml (user home)
6. <install>/config/config.yaml (built-in defaults)
7. Hardcoded defaults in code
```

**Fallback strategy:** If no config file found, use hardcoded defaults (current behavior).

---

## Benefits

1. **Single source of truth** - No more confusion about which file to edit
2. **Clear hierarchy** - Sections logically organized (scoring, training, optimizer, etc.)
3. **Industry standard** - YAML is Python ML ecosystem standard
4. **Better defaults** - Profiles for different use cases (quick, production)
5. **Auto-loading** - Just works without specifying config path
6. **CLI overrides** - Command-line args always win (important for scripts)
7. **Backward compatible** - Keep old flags working (with deprecation warnings)

---

## Implementation Checklist

- [ ] Create `config/config.yaml` with unified structure
- [ ] Create config profiles (quick, production, minimal)
- [ ] Update `cli/config_loader.py` for auto-loading
- [ ] Add `--generate-config` command
- [ ] Update CLI argument parsing to check config hierarchy
- [ ] Search and replace references to old config files
- [ ] Add backward compatibility for `--optimizer-config`
- [ ] Move old files to deprecated/
- [ ] Update documentation (README, CLI help, CONFIGURATION.md)
- [ ] Test with various priority scenarios
- [ ] Add deprecation warnings for old config files

---

## Testing Plan

1. **Test auto-loading:**
   ```bash
   # No config file - should use defaults
   intronIC train -n test1

   # Project config - should auto-load ./.intronIC.yaml
   echo "training: {n_models: 3}" > .intronIC.yaml
   intronIC train -n test2  # Should use 3 models

   # User config - should auto-load ~/.config/intronIC/config.yaml
   mkdir -p ~/.config/intronIC
   cp config/config.yaml ~/.config/intronIC/
   intronIC train -n test3
   ```

2. **Test CLI override:**
   ```bash
   # Config says n_models=10, CLI says 1 - CLI should win
   intronIC train -n test4 --n-models 1
   ```

3. **Test explicit config:**
   ```bash
   # Should use quick profile
   intronIC train -n test5 --config config/profiles/quick.yaml
   ```

4. **Test backward compatibility:**
   ```bash
   # Old flag should still work with deprecation warning
   intronIC train -n test6 --optimizer-config config/training_default.yaml
   ```

---

## Timeline

- **Planning:** 30 minutes (DONE - this document)
- **Implementation:** 2-3 hours
  - Unified config: 30 min
  - Config loader: 1 hour
  - Code updates: 1 hour
  - Testing: 30 min
- **Documentation:** 30 minutes
- **Total:** ~3-4 hours

---

## Open Questions

1. Keep `--optimizer-config` for backward compatibility or deprecate immediately?
   - **Recommendation:** Keep with deprecation warning

2. Support both TOML and YAML during transition period?
   - **Recommendation:** No - clean break, but detect TOML and show helpful error

3. What to do with existing `training_quick.yaml` and `training_fast_test.yaml`?
   - **Recommendation:** Repurpose as profiles under `config/profiles/`

4. Should profiles inherit from base config.yaml or be standalone?
   - **Recommendation:** Standalone for simplicity, but document common patterns

---

## Notes

- Current code already has good YAML parsing in `classification/config_loader.py`
- CLI args system in `cli/args.py` and `cli/config.py` is well-structured
- Main work is consolidation and auto-loading logic
- Should be straightforward migration with clear benefits
