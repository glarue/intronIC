# intronIC Configuration Files

This directory contains configuration files for intronIC.

## Configuration System

intronIC uses a unified YAML configuration system that handles all settings - both general parameters (scoring, output, performance) and training/optimization parameters.

**File:** `config.yaml` (default configuration)
**Format:** YAML
**Purpose:** All intronIC settings including scoring, training, optimization, output, and performance

---

## Configuration File (`config.yaml`)

The default `config.yaml` contains comprehensive settings with inline documentation for all intronIC parameters.

**Settings include:**
- **Scoring:** threshold, PWM scoring regions, pseudocount, boundary correction
- **Training:** ensemble configuration, evaluation mode, CV folds
- **Optimization:** hyperparameter grid search, C bounds, feature selection
- **Performance:** parallel processes, min intron length, streaming mode
- **Isoforms:** isoform selection, duplicate handling
- **Output:** naming conventions, file formats
- **Advanced:** random seed, debug mode

### Usage:

**Generate a template:**
```bash
intronIC --generate-config > my_config.yaml
```

**Use a custom config:**
```bash
# Specify config file location (--config must come BEFORE subcommand)
intronIC --config path/to/my_config.yaml \
  classify -g genome.fa -a annotation.gff -n species

# Or place in one of the auto-detected locations (in order of precedence):
# 1. --config PATH (explicit CLI argument)
# 2. ./.intronIC.yaml (current directory)
# 3. ~/.config/intronIC/config.yaml (XDG config dir)
# 4. ~/.intronIC.yaml (user home)
# 5. <install_dir>/config/config.yaml (built-in defaults)
```

**Example use case:** Set project-wide defaults for threshold, parallel processes, training parameters, and output naming without repeating CLI args.

---

## Common Configuration Tasks

### Quick Classification (use pretrained model)
```bash
# Uses default settings and pretrained model
intronIC classify -g genome.fa -a annotation.gff -n species
```

### Classification with Custom Threshold
Create a config file:
```yaml
scoring:
  threshold: 95.0  # More stringent than default 90%
```

Then run:
```bash
intronIC --config my_config.yaml \
  classify -g genome.fa -a annotation.gff -n species
```

### Training a New Model
For production training with full optimization:
```bash
# Uses all default training settings from config.yaml
intronIC train -n homo_sapiens -p 12
```

For faster testing during development:
```yaml
# fast_test.yaml
training:
  n_cv_folds: 3
  ensemble:
    n_models: 5
    max_iter: 10000

optimizer:
  n_rounds: 3
  n_points_initial: 7
  n_points_refine: 30
  cv_folds: 3
  max_iter: 10000
```

Then run:
```bash
intronIC --config fast_test.yaml train -n homo_sapiens -p 8
```

### Parallel Processing
```yaml
performance:
  processes: 8  # Use 8 cores for classification
```

Or via CLI:
```bash
intronIC classify -g genome.fa -a annotation.gff -n species -p 8
```

---

## Creating Custom Configurations

To create a custom configuration:

1. **Generate a template** with `intronIC --generate-config > my_config.yaml`
2. **Modify sections** based on your needs:
   - Adjust `scoring.threshold` for more/less stringent classification
   - Modify `training.ensemble.n_models` for faster/more robust training
   - Reduce `optimizer.n_rounds` and grid parameters for faster optimization
   - Set `performance.processes` for your CPU count
3. **Test with Chr19** data first to estimate runtime
4. **Document** what you changed and why in comments

### Example: Custom Config for Fast Testing

```yaml
# Fast testing configuration
# Reduced parameters for quick iteration during development

scoring:
  threshold: 90.0

training:
  n_cv_folds: 3  # Reduced from default 7
  ensemble:
    n_models: 5  # Reduced from default 16
    max_iter: 10000  # Reduced from default 50000

optimizer:
  n_rounds: 3  # Reduced from default 7
  n_points_initial: 7  # Reduced from default 20
  n_points_refine: 30  # Reduced from default 100
  cv_folds: 3  # Reduced from default 7
  max_iter: 10000  # Reduced from default 65000

performance:
  processes: 8
```

### Example: Custom Config for Production

```yaml
# Production configuration
# Optimized for robust, publication-quality results

scoring:
  threshold: 90.0

training:
  eval_mode: nested_cv
  n_cv_folds: 7
  ensemble:
    n_models: 16
    subsample_u2: true
    subsample_ratio: 0.8
    max_iter: 50000

optimizer:
  n_rounds: 7
  n_points_initial: 20
  n_points_refine: 100
  cv_folds: 7
  max_iter: 65000
  n_jobs: -1  # Use all available cores

performance:
  processes: 12
```

---

## Configuration Sections

The config file supports these main sections:

### `scoring`
Classification and PWM scoring parameters:
- `threshold`: U12 probability threshold (0-100%)
- `feature_type`: Extract from CDS, exons, or both
- `exclude_noncanonical`: Filter non-canonical introns
- `u12_boundary_correction`: Auto-correct annotation errors
- `regions`: Scoring window coordinates (5', BP, 3')

### `training`
Model training configuration:
- `eval_mode`: Evaluation strategy (nested_cv, split, none)
- `n_cv_folds`: Cross-validation folds for outer loop
- `ensemble`: Number of models, subsampling strategy
- `use_fold_averaged_params`: Use CV-averaged hyperparameters

### `optimizer`
Hyperparameter optimization settings:
- `n_rounds`: Geometric refinement rounds
- `n_points_initial`/`n_points_refine`: Grid density
- `cv_folds`: Cross-validation folds for inner loop
- `max_iter`: Maximum SVM iterations
- `scoring_metric`: Optimization metric (balanced_accuracy, f_beta)
- `penalty_options`: L1/L2 regularization
- `class_weight_multipliers`: Precision/recall tradeoff
- `feature_transform`: Which composite features to include

### `performance`
Runtime and resource settings:
- `processes`: Parallel processes for classification
- `cv_processes`: Parallel processes for CV
- `min_intron_length`: Minimum intron size

### `isoform`
Transcript selection rules:
- `allow_multiple_isoforms`: Include all vs longest only
- `exclude_overlapping`: Filter overlapping introns
- `include_duplicates`: Keep duplicate coordinates

### `output`
Output file formatting:
- `clean_names`: Remove "transcript:" prefixes

### `extraction`
Sequence extraction parameters:
- `flank_length`: Exonic flank size

### `advanced`
Advanced settings:
- `random_seed`: Reproducibility seed
- `quiet`: Suppress output
- `debug`: Enable debug logging

---

## Quick Reference

### I want to...

**Change classification threshold (e.g., 90% → 95%)**
```yaml
scoring:
  threshold: 95.0
```

**Speed up training for testing**
```yaml
training:
  n_cv_folds: 3
  ensemble:
    n_models: 5

optimizer:
  n_rounds: 3
  cv_folds: 3
```

**Use more CPU cores**
```yaml
performance:
  processes: 16
```

**Include all isoforms instead of longest only**
```yaml
isoform:
  allow_multiple_isoforms: true
```

**Enable debug logging**
```yaml
advanced:
  debug: true
```

---

## Notes

- **Default behavior:** If no `--config` is specified, intronIC uses built-in defaults from `config/config.yaml`

- **CLI overrides:** Command-line arguments always override config file values

- **Convergence:** If you see convergence warnings, increase `optimizer.max_iter` or `training.ensemble.max_iter`

- **Parallelization:** Using `-p` (or `performance.processes`) significantly speeds up classification. For training, `optimizer.n_jobs: -1` uses all cores for grid search
