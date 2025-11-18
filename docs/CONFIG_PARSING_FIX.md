# YAML Config Parsing Fix (2025-01-17)

## Issue

YAML training settings were not being loaded and applied during training, causing models to use CLI defaults instead of YAML configuration values.

**User report:** "I just ran it with n_models: 7 but it only reported training 4 in the log"

## Root Causes

### Root Cause 1: Missing YAML settings extraction

The config loader (`classification/config_loader.py`) was extracting the `training` section from YAML and attaching values as attributes on the optimizer object (e.g., `training_n_models`, `training_subsample_u2`), but `cli/main.py` was:

1. Only checking for `training_n_models` from YAML
2. Not checking for other training settings (`subsample_u2`, `subsample_ratio`, `max_iter`, `random_state`)
3. Not passing these values to the `IntronClassifier` constructor

### Root Cause 2: CLI argparse default overriding YAML

The `train` subcommand had `--n_models` with `default=4` in `cli/args.py:198`, which was inconsistent with:
- `TrainingConfig` dataclass default: `n_models: int = 1` (cli/config.py:97)
- Production YAML config: `n_models: 7` (config/training_default.yaml:169)

This caused the argparse default to override the YAML value, even when YAML was explicitly loaded.

### Root Cause 3: Incorrect dataclass modification method

The code was using `config._replace()` (namedtuple method) instead of `dataclasses.replace()` for frozen dataclasses:
- `_replace()` is a method on `namedtuple` objects
- Frozen dataclasses require `dataclasses.replace()` function
- This caused `AttributeError: 'IntronICConfig' object has no attribute '_replace'`

## YAML Training Settings

The `training` section in YAML configs specifies:

```yaml
training:
  n_models: 7                # Number of ensemble models
  subsample_u2: true         # Whether to subsample U2 for diversity
  subsample_ratio: 0.85      # Fraction of U2s per model
  max_iter: 25000            # Max iterations for ensemble training
  random_state: 42           # Random seed for U2 subsampling
```

**Important distinction:**
- `optimizer.max_iter`: Used during hyperparameter optimization (grid search)
- `training.max_iter`: Used during final ensemble training

## Fix Applied

### 1. classification/config_loader.py (lines 75-101)
**Status:** Already correct - was properly extracting and attaching training settings as attributes.

```python
# Extract training settings (ensemble configuration)
training_config = config.get('training', {})

# Attach training settings as attributes for later use
if training_config:
    for key, value in training_config.items():
        setattr(optimizer, f'training_{key}', value)
```

### 2. cli/main.py (lines 1323-1360)
**Changed:** Extract ALL training settings from YAML, not just `n_models`.

```python
# Initialize YAML training settings to None (may be overridden if YAML config provided)
yaml_n_models = None
yaml_subsample_u2 = None
yaml_subsample_ratio = None
yaml_training_max_iter = None
yaml_training_random_state = None

if config.training.optimizer_config_path:
    # ... load optimizer_from_yaml ...

    # Extract training settings if present in YAML (override CLI defaults)
    yaml_n_models = getattr(optimizer_from_yaml, 'training_n_models', None)
    yaml_subsample_u2 = getattr(optimizer_from_yaml, 'training_subsample_u2', None)
    yaml_subsample_ratio = getattr(optimizer_from_yaml, 'training_subsample_ratio', None)
    yaml_training_max_iter = getattr(optimizer_from_yaml, 'training_max_iter', None)
    yaml_training_random_state = getattr(optimizer_from_yaml, 'training_random_state', None)

    # Log each setting if found
    if yaml_n_models is not None:
        config = config._replace(training=config.training._replace(n_models=yaml_n_models))
        messenger.log_only(f"  Ensemble models: {yaml_n_models} (from YAML)")
    if yaml_subsample_u2 is not None:
        messenger.log_only(f"  Subsample U2: {yaml_subsample_u2} (from YAML)")
    # ... etc for other settings
```

### 3. cli/main.py (lines 1375-1413)
**Changed:** Pass all training settings to `IntronClassifier` constructor.

```python
# Use YAML training settings if available, otherwise use CLI defaults
ensemble_subsample_u2 = yaml_subsample_u2 if yaml_subsample_u2 is not None else True
ensemble_subsample_ratio = yaml_subsample_ratio if yaml_subsample_ratio is not None else 0.8
ensemble_max_iter = yaml_training_max_iter if yaml_training_max_iter is not None else optimizer_max_iter
ensemble_random_state = yaml_training_random_state if yaml_training_random_state is not None else config.training.seed

classifier = IntronClassifier(
    n_optimization_rounds=optimizer_n_rounds,
    classification_threshold=config.scoring.threshold,
    n_ensemble_models=config.training.n_models,
    subsample_u2=ensemble_subsample_u2,          # Now from YAML!
    subsample_ratio=ensemble_subsample_ratio,     # Now from YAML!
    fixed_c=config.training.fixed_C,
    optimize_c=(config.training.fixed_C is None),
    random_state=ensemble_random_state,           # Now from YAML!
    cv_processes=optimizer_cv_processes,
    classification_processes=config.performance.processes,
    max_iter=ensemble_max_iter,                   # Now from YAML training.max_iter!
    eval_mode=config.training.eval_mode,
    n_cv_folds=optimizer_cv_folds,
    test_fraction=config.training.test_fraction,
    param_grid_override=param_grid_override,
    n_points_initial=optimizer_n_points_initial,
    eff_C_pos_range=eff_C_pos_range,
    eff_C_neg_max=eff_C_neg_max
)
```

## Verification

Tested YAML loader successfully extracts all training settings:

```bash
$ pixi run python -c "from classification.config_loader import load_optimizer_config; \
  optimizer = load_optimizer_config('config/training_default.yaml'); \
  print(f'n_models: {getattr(optimizer, \"training_n_models\", \"NOT FOUND\")}'); \
  print(f'subsample_u2: {getattr(optimizer, \"training_subsample_u2\", \"NOT FOUND\")}'); \
  print(f'subsample_ratio: {getattr(optimizer, \"training_subsample_ratio\", \"NOT FOUND\")}'); \
  print(f'max_iter: {getattr(optimizer, \"training_max_iter\", \"NOT FOUND\")}'); \
  print(f'random_state: {getattr(optimizer, \"training_random_state\", \"NOT FOUND\")}')"

Training settings from YAML:
  n_models: 7
  subsample_u2: True
  subsample_ratio: 0.85
  max_iter: 25000
  random_state: 42
```

✓ All values correctly extracted and will be applied during training

## Fixes Applied

### Fix 1: Extract all YAML training settings (cli/main.py)
**Status:** ✅ Complete

Lines 1323-1413 now extract and apply all 5 training settings from YAML:
- `n_models`
- `subsample_u2`
- `subsample_ratio`
- `max_iter` (training, not optimizer)
- `random_state`

### Fix 2: Align argparse default with TrainingConfig (cli/args.py)
**Status:** ✅ Complete

Changed line 198: `default=4` → `default=1`

This ensures:
- CLI default matches `TrainingConfig` dataclass default (1)
- YAML values properly override (no conflict with argparse default of 4)
- Production config (`training_default.yaml` with `n_models: 7`) works as expected

### Fix 3: Use correct dataclass modification method (cli/main.py)
**Status:** ✅ Complete

Changed all `config._replace()` calls to use `dataclasses.replace()`:

**Line 18:** Added import `from dataclasses import replace`

**Lines 1338-1341:** Auto-load config update:
```python
# Before: config._replace(training=config.training._replace(...))
# After:
config = replace(
    config,
    training=replace(config.training, optimizer_config_path=default_config)
)
```

**Lines 1370-1373:** n_models update from YAML:
```python
# Before: config._replace(training=config.training._replace(n_models=...))
# After:
config = replace(
    config,
    training=replace(config.training, n_models=yaml_n_models)
)
```

This fixes the `AttributeError: 'IntronICConfig' object has no attribute '_replace'` error.

## Files Modified

1. `cli/main.py` - Extract and apply all YAML training settings (lines 1323-1413)
2. `cli/main.py` - Fixed dataclass modification method (lines 18, 1338-1341, 1370-1373)
3. `cli/args.py` - Changed `--n_models` default from 4 to 1 (line 198)

## Impact

Training runs using `--optimizer-config` will now correctly apply all settings from the YAML `training` section:

- **n_models**: Number of ensemble models
- **subsample_u2**: Whether to subsample U2 for diversity
- **subsample_ratio**: Fraction of U2s to use per model
- **max_iter**: Max iterations during ensemble training (separate from optimization max_iter)
- **random_state**: Random seed for reproducible U2 subsampling

## Usage

### Auto-loading (NEW - 2025-01-17)

`config/training_default.yaml` is now **auto-loaded** during training if it exists and no explicit `--optimizer-config` is provided:

```bash
# Simply run train - config/training_default.yaml is loaded automatically
intronIC train -n homo_sapiens

# Expected log output:
# Auto-loading default configuration: .../config/training_default.yaml
# Loading optimizer configuration from: .../config/training_default.yaml
#   Ensemble models: 7 (from YAML)
#   Subsample U2: True (from YAML)
#   Subsample ratio: 0.85 (from YAML)
#   Training max_iter: 25000 (from YAML)
#   Training random_state: 42 (from YAML)
# Loaded custom optimizer configuration:
#   Optimization rounds: 5
#   ...
```

### Explicit config (still supported)

```bash
# Use a different config file
intronIC train -n homo_sapiens --optimizer-config config/training_quick.yaml
```

### Disabling auto-load

To prevent auto-loading and use CLI defaults only, rename or move `config/training_default.yaml`.

## Related Files

- CLI improvements summary: `CLI_IMPROVEMENTS_SUMMARY.md`
- Default training config: `config/training_default.yaml`
- Config loader: `classification/config_loader.py`
- Classifier API: `classification/classifier.py` (IntronClassifier.__init__)
