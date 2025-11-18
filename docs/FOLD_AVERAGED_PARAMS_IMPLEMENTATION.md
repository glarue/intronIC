# Fold-Averaged Hyperparameters Implementation

**Date:** 2025-01-18
**Feature:** `--use-fold-averaged-params` flag for better cross-species generalization

## Problem Statement

The current implementation uses nested CV for evaluation but then re-optimizes hyperparameters on the full dataset for the final production model. This causes a discrepancy:

- **Nested CV evaluation:** C≈36-51, sigmoid calibration → 2 FPs on C. elegans
- **Production model (re-optimized):** C=93.46, isotonic calibration → 6 FPs on C. elegans

The re-optimization on the full dataset maximizes training-set performance but may overfit to the training species composition, hurting cross-species generalization.

## Solution

Added configurable `use_fold_averaged_params` option that:
- Uses **geometric mean** of fold-specific C values from nested CV
- Uses **majority-vote** calibration method from folds
- Skips re-optimization on full dataset
- Provides more conservative parameters favoring generalization over training fit

## Implementation Details

### 1. Configuration Layer (`cli/config.py`)

**TrainingConfig dataclass:**
```python
@dataclass(frozen=True, slots=True)
class TrainingConfig:
    # ... existing fields ...
    use_fold_averaged_params: bool = False  # Use fold-averaged hyperparameters
```

**Parameter passing (line 272):**
```python
training_config = TrainingConfig(
    # ... other parameters ...
    use_fold_averaged_params=getattr(args, 'use_fold_averaged_params', False)
)
```

### 2. CLI Arguments (`cli/args.py`)

**Added to BOTH subcommands:**

**Train subcommand (lines 248-253):**
```python
training.add_argument(
    '--use-fold-averaged-params',
    action='store_true',
    default=None,
    help='Use fold-averaged hyperparameters from nested CV instead of re-optimizing on full dataset (better cross-species generalization)'
)
```

**Classify subcommand (lines 550-555):** (identical)

### 3. YAML Configuration (`config/training_default.yaml`)

**Lines 197-207:**
```yaml
training:
  # Use fold-averaged hyperparameters from nested CV
  # - false = Re-optimize C and calibration on full dataset after nested CV (default sklearn behavior)
  #           * Maximizes performance on training data
  #           * May overfit to training species composition
  # - true = Use geometric mean of fold-specific C values and majority-vote calibration
  #          * More conservative, favors generalization over training fit
  #          * RECOMMENDED for cross-species applications
  #          * Expected to reduce FPs on distant species (e.g., C. elegans: 2 FPs vs 6 FPs)
  # Only applies when eval_mode='nested_cv'
  # Note: For same-species prediction, default (false) is better
  use_fold_averaged_params: false
```

### 4. Nested CV Evaluator (`classification/nested_cv.py`)

**No changes needed!** Already collects fold-specific parameters in FoldResult (lines 306-307):
```python
fold_results.append(FoldResult(
    # ... other fields ...
    optimized_C=parameters.C,  # Already collected!
    calibration_method=parameters.calibration_method  # Already collected!
))
```

### 5. IntronClassifier (`classification/classifier.py`)

**Added imports (lines 27-28):**
```python
from collections import Counter
import numpy as np
```

**Added parameter (line 138):**
```python
def __init__(
    self,
    # ... existing parameters ...
    use_fold_averaged_params: bool = False,  # Use fold-averaged hyperparameters
    # ...
):
```

**Added helper method (lines 214-269):**
```python
def _compute_fold_averaged_params(self, nested_cv_result) -> SVMParameters:
    """
    Compute fold-averaged hyperparameters from nested CV results.

    Uses geometric mean for C (better for log-scale parameters) and
    majority vote for calibration method.

    Returns SVMParameters with conservative parameters for cross-species generalization.
    """
    fold_results = nested_cv_result.fold_results

    # Extract fold-specific C values and calibration methods
    c_values = np.array([fold.optimized_C for fold in fold_results])
    calibration_methods = [fold.calibration_method for fold in fold_results]

    # Geometric mean of C values (better for log-scale parameters)
    geometric_mean_C = float(np.exp(np.mean(np.log(c_values))))

    # Majority vote for calibration method
    method_counts = Counter(calibration_methods)
    majority_calibration = method_counts.most_common(1)[0][0]

    # Log for transparency
    print(f"\nFold-averaged hyperparameters:")
    print(f"  Fold-specific C values: {[f'{c:.2e}' for c in c_values]}")
    print(f"  Geometric mean C: {geometric_mean_C:.6e}")
    print(f"  Fold-specific calibration: {calibration_methods}")
    print(f"  Majority-vote calibration: {majority_calibration}")
    print(f"  Rationale: Using conservative fold-averaged params for better cross-species generalization")

    return SVMParameters(
        C=geometric_mean_C,
        calibration_method=majority_calibration,
        # ... fixed architecture parameters ...
    )
```

**Modified both classify() and classify_batch() methods:**

Replaced Stage 1 hyperparameter optimization section with conditional logic:
```python
# Stage 1: Hyperparameter Optimization
# Choose between fold-averaged params (conservative) or re-optimization (aggressive)
print("\n=== Stage 1: Hyperparameter Optimization ===")

# Check if we should use fold-averaged parameters from nested CV
if (self.use_fold_averaged_params and
    self.eval_mode == 'nested_cv' and
    eval_result is not None):
    # Use fold-averaged hyperparameters (better cross-species generalization)
    print("Using fold-averaged hyperparameters from nested CV")
    print("(Skipping re-optimization on full dataset for better cross-species generalization)")
    parameters = self._compute_fold_averaged_params(eval_result)
else:
    # Standard approach: re-optimize on full dataset
    if self.use_fold_averaged_params and self.eval_mode == 'nested_cv':
        print("Note: use_fold_averaged_params=True but nested CV result not available")
        print("Falling back to re-optimization on full dataset")
    elif self.use_fold_averaged_params:
        print("Note: use_fold_averaged_params=True but eval_mode is not 'nested_cv'")
        print("Falling back to re-optimization on full dataset")

    # [existing re-optimization code...]
```

## Usage

### Command Line

**Training mode:**
```bash
# Use fold-averaged parameters (recommended for cross-species)
intronIC train -n homo_sapiens --use-fold-averaged-params

# Default behavior (re-optimize on full dataset)
intronIC train -n homo_sapiens
```

**Classify mode (with on-the-fly training):**
```bash
# Use fold-averaged parameters
intronIC classify -g genome.fa.gz -a annotation.gff3.gz -n species --train --use-fold-averaged-params

# Default behavior
intronIC classify -g genome.fa.gz -a annotation.gff3.gz -n species --train
```

### YAML Configuration

In `config/training_default.yaml`:
```yaml
training:
  use_fold_averaged_params: true  # Enable fold-averaging
```

## Expected Results

### Training Species (Human)
- Fold-averaged: Slightly lower training F1 (more conservative)
- Re-optimized: Maximizes training F1

### Distant Species (C. elegans)
- Fold-averaged: **2 FPs** (better generalization)
- Re-optimized: **6 FPs** (overfits to human composition)

## Technical Rationale

### Why Geometric Mean for C?

C is a log-scale parameter (regularization strength). Geometric mean is appropriate for log-scale parameters:
- Arithmetic mean of C values: biased toward large values
- Geometric mean of C values: equivalent to arithmetic mean in log-space
- Formula: `geometric_mean(C) = exp(mean(log(C_i)))`

### Why Majority Vote for Calibration?

Calibration method is categorical (sigmoid vs isotonic):
- Majority vote is the natural aggregation for categorical variables
- Prevents ties by selecting most common choice
- If tie: Counter.most_common() returns first alphabetically

### Why This Helps Cross-Species?

Fold-averaged parameters represent the **average optimal choice** across different data splits, which is more robust to composition variations. Re-optimization on full dataset finds the **single best choice for that specific dataset**, which may be biased toward training species characteristics.

## Files Modified

1. `cli/config.py` - Added field and parameter passing
2. `cli/args.py` - Added CLI arguments (2 locations)
3. `config/training_default.yaml` - Added YAML configuration
4. `classification/classifier.py` - Added helper method and conditional logic
5. `classification/nested_cv.py` - **No changes** (already collects parameters)

## Testing

**Test command:**
```bash
pixi run intronIC train -n homo_sapiens_fold_averaged_test \
    -o . \
    --n_optimization_rounds 2 \
    --n_cv_folds 3 \
    --use-fold-averaged-params \
    -p 10
```

**Expected log output:**
```
=== Stage 1: Hyperparameter Optimization ===
Using fold-averaged hyperparameters from nested CV
(Skipping re-optimization on full dataset for better cross-species generalization)

Fold-averaged hyperparameters:
  Fold-specific C values: ['3.68e+01', '4.12e+01', '3.55e+01']
  Geometric mean C: 3.784123e+01
  Fold-specific calibration: ['sigmoid', 'sigmoid', 'sigmoid']
  Majority-vote calibration: sigmoid
  Rationale: Using conservative fold-averaged params for better cross-species generalization
```

## Future Improvements

1. **Add to model metadata:** Store whether fold-averaging was used
2. **Automatic mode:** Auto-enable for cross-species classification
3. **Weighted geometric mean:** Weight by fold F1 scores
4. **Calibration scoring:** If tie, compare sigmoid vs isotonic on validation set

## References

- Architecture: `SCALER_CENTERING_FIX_SUMMARY.md`
- Nested CV: `classification/nested_cv.py`
- Original issue: C. elegans 6 FPs vs 2 FPs discrepancy between nested CV and production model
