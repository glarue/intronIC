# Feature Selection Implementation Plan

## Goal
Enable flexible feature selection in intronIC to allow training/classification with subsets of the three base motifs (5'SS, BP, 3'SS) instead of requiring all three.

**Use cases:**
- Train using only 5'SS scores (1D model)
- Train using 5'SS + BP (2D model)
- Train using 5'SS + 3'SS (2D model, skip BP)
- Current default: All three (3D→7D with composite features)

---

## Architecture Overview

### Current Pipeline Flow
```
1. Scoring → [five_raw, bp_raw, three_raw]
2. Normalization → [five_z, bp_z, three_z]
3. BothEndsStrongTransformer → [five_z, bp_z, three_z, composite_features...]
4. LinearSVC → Classification
```

### Required Changes
```
1. Config: Specify which features to use
2. Scoring: Only score selected motifs (skip disabled)
3. Normalization: Only normalize selected features
4. Transformer: Handle variable dimensions (1D, 2D, or 3D)
5. Training/Prediction: Build feature matrix from selected features only
```

---

## Complexity Assessment

### **Complexity Level: MODERATE-HIGH** ⚠️

**Why moderate-high:**
- Touches ~10 core files across multiple modules
- Changes frozen dataclass configs (requires careful migration)
- Transformer logic needs complete rewrite for flexibility
- Feature names/indices scattered throughout codebase
- Backward compatibility concerns with existing models
- Testing matrix grows significantly (1D, 2D combinations, 3D)

**Estimated effort:** 4-8 hours of focused work + 2-4 hours testing

---

## Detailed Implementation Plan

## Phase 1: Configuration Changes

### 1.1 Add Feature Selection to Config
**Files:** `src/intronIC/cli/config.py`, `config/config.yaml`

**Changes:**

#### config.yaml
```yaml
scoring:
  # Feature selection: which motifs to score and use in classification
  # Options: 'five' (5'SS), 'bp' (branch point), 'three' (3'SS)
  # Default: all three for full 3D feature space
  # Examples:
  #   - ['five'] → 1D model (5'SS only)
  #   - ['five', 'bp'] → 2D model (no 3'SS)
  #   - ['five', 'bp', 'three'] → 3D model (current default)
  enabled_features: ['five', 'bp', 'three']

  # Scoring regions (only used if feature is enabled)
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
```

#### config.py
```python
@dataclass(frozen=True, slots=True)
class ScoringConfig:
    """Configuration for scoring and classification."""
    # Feature selection
    enabled_features: tuple[str, ...] = ('five', 'bp', 'three')

    # Validate enabled_features in __post_init__
    def __post_init__(self):
        valid = {'five', 'bp', 'three'}
        for feat in self.enabled_features:
            if feat not in valid:
                raise ValueError(f"Invalid feature: {feat}. Must be one of {valid}")
        if not self.enabled_features:
            raise ValueError("At least one feature must be enabled")
```

**Complexity:** LOW
**Risks:** Frozen dataclass requires workaround for validation
**Alternative:** Use `@validator` from pydantic if we switch config system

---

## Phase 2: Scoring Changes

### 2.1 Conditional Scoring in PWM Scorer
**File:** `src/intronIC/scoring/scorer.py`

**Current:** Always scores all three regions
**Change:** Only score enabled features

```python
class IntronScorer:
    def __init__(self, pwm_sets, enabled_features=('five', 'bp', 'three')):
        self.pwm_sets = pwm_sets
        self.enabled_features = set(enabled_features)

    def score_intron(self, intron, ...):
        """Score only enabled features, set others to None."""
        # Score 5'SS (if enabled)
        if 'five' in self.enabled_features:
            five_score = self._score_region(...)
        else:
            five_score = None

        # Score BP (if enabled)
        if 'bp' in self.enabled_features:
            bp_score = self._score_bp_region(...)
        else:
            bp_score = None

        # Score 3'SS (if enabled)
        if 'three' in self.enabled_features:
            three_score = self._score_region(...)
        else:
            three_score = None
```

**Complexity:** LOW-MODERATE
**Impact:** ~50 lines changed in scorer.py
**Testing:** Verify None values handled correctly downstream

---

## Phase 3: Normalization Changes

### 3.1 ScoreNormalizer with Dynamic Features
**File:** `src/intronIC/scoring/normalizer.py`

**Current:** Always expects 3 scores
**Change:** Handle variable number of enabled features

```python
class ScoreNormalizer:
    def __init__(self, enabled_features=('five', 'bp', 'three')):
        self.enabled_features = tuple(enabled_features)
        self.n_features = len(enabled_features)
        self.feature_to_idx = {f: i for i, f in enumerate(enabled_features)}

    def fit(self, introns):
        """Fit using only enabled features."""
        # Build score matrix with only enabled features
        score_matrix = []
        for intron in introns:
            row = []
            for feat in self.enabled_features:
                if feat == 'five':
                    row.append(intron.scores.five_raw_score)
                elif feat == 'bp':
                    row.append(intron.scores.bp_raw_score)
                elif feat == 'three':
                    row.append(intron.scores.three_raw_score)

            # Skip if any enabled feature is None
            if None not in row:
                score_matrix.append(row)

        # Fit RobustScaler on (n_samples, n_features)
        self.scaler_.fit(np.array(score_matrix))

    def transform(self, introns):
        """Transform using only enabled features."""
        # Similar logic: build matrix, transform, assign back
        # Set z-scores to None for disabled features
```

**Complexity:** MODERATE
**Impact:** ~100 lines changed in normalizer.py
**Risks:** Need to handle None values carefully

---

## Phase 4: Feature Transformer Changes

### 4.1 Flexible BothEndsStrongTransformer
**File:** `src/intronIC/classification/transformers.py`

**Current:** Hardcoded 3D input → 7D/9D/11D output
**Change:** Handle 1D/2D/3D input dynamically

```python
class BothEndsStrongTransformer(BaseEstimator, TransformerMixin):
    def __init__(
        self,
        enabled_features=('five', 'bp', 'three'),
        composite_features=None,
        gamma_imbalance=1.0
    ):
        """
        Args:
            enabled_features: Which base features are available
            composite_features: Which composite features to compute
                Valid only if corresponding base features enabled:
                - 'min_5_bp', 'absdiff_5_bp', 'max_5_bp': need 'five' + 'bp'
                - 'min_5_3', 'absdiff_5_3', 'max_5_3': need 'five' + 'three'
                - 'absdiff_bp_3': needs 'bp' + 'three'
                - 'min_all': needs all three
        """
        self.enabled_features = tuple(enabled_features)
        self.composite_features = composite_features or []
        self.gamma_imbalance = gamma_imbalance

        # Validate composite features
        self._validate_composite_features()

    def _validate_composite_features(self):
        """Check that composite features match enabled base features."""
        enabled_set = set(self.enabled_features)

        for comp in self.composite_features:
            if comp in ['min_5_bp', 'absdiff_5_bp', 'max_5_bp']:
                if not ({'five', 'bp'} <= enabled_set):
                    raise ValueError(f"{comp} requires 'five' and 'bp'")
            elif comp in ['min_5_3', 'absdiff_5_3', 'max_5_3']:
                if not ({'five', 'three'} <= enabled_set):
                    raise ValueError(f"{comp} requires 'five' and 'three'")
            elif comp == 'absdiff_bp_3':
                if not ({'bp', 'three'} <= enabled_set):
                    raise ValueError(f"{comp} requires 'bp' and 'three'")
            elif comp == 'min_all':
                if len(enabled_set) < 3:
                    raise ValueError("min_all requires all three features")

    def transform(self, X):
        """
        Transform variable-dimension input to augmented features.

        Args:
            X: Array of shape (n_samples, n_features)
               where n_features = len(enabled_features)

        Returns:
            Array of shape (n_samples, n_base + n_composite)
        """
        X = np.asarray(X)
        n_expected = len(self.enabled_features)

        if X.shape[1] != n_expected:
            raise ValueError(
                f"Expected {n_expected} features {self.enabled_features}, "
                f"got {X.shape[1]}"
            )

        # Start with base features
        features = [X]

        # Extract base features by name
        feature_map = {}
        for i, name in enumerate(self.enabled_features):
            feature_map[name] = X[:, i]

        # Compute composite features
        for comp in self.composite_features:
            if comp == 'min_5_bp':
                val = np.minimum(feature_map['five'], feature_map['bp'])
                features.append(val.reshape(-1, 1))
            elif comp == 'absdiff_5_bp':
                val = np.abs(feature_map['five'] - feature_map['bp'])
                val *= self.gamma_imbalance
                features.append(val.reshape(-1, 1))
            # ... etc for all composite features

        return np.hstack(features)

    def get_feature_names(self):
        """Get names of output features."""
        names = list(self.enabled_features)
        names.extend(self.composite_features)
        return names
```

**Complexity:** HIGH
**Impact:** Complete rewrite of transformer (~200 lines)
**Testing:** Need to test all combinations (7 possible base feature sets)

---

## Phase 5: Training & Prediction Changes

### 5.1 Update Trainer
**File:** `src/intronIC/classification/trainer.py`

**Changes:**
- Pass `enabled_features` to normalizer
- Pass `enabled_features` and `composite_features` to transformer
- Update feature matrix construction
- Log which features are being used

```python
def train_ensemble(..., enabled_features=('five', 'bp', 'three')):
    # Create normalizer with feature selection
    normalizer = ScoreNormalizer(enabled_features=enabled_features)

    # Create transformer with feature selection
    transformer = BothEndsStrongTransformer(
        enabled_features=enabled_features,
        composite_features=composite_features
    )

    # Log configuration
    messenger.info(f"Training with features: {enabled_features}")
    messenger.info(f"Final feature space: {transformer.get_feature_names()}")
```

**Complexity:** LOW
**Impact:** ~30 lines changed

### 5.2 Update Predictor
**File:** `src/intronIC/classification/predictor.py`

**Changes:**
- Extract `enabled_features` from saved model metadata
- Build feature matrix using only enabled features
- Handle models trained with different feature sets

**Complexity:** LOW-MODERATE
**Impact:** ~50 lines changed

---

## Phase 6: Model Persistence Changes

### 6.1 Save Feature Configuration in Model
**File:** `src/intronIC/utils/model_io.py`

**Change:** Include feature configuration in model metadata

```python
def save_model(...):
    metadata = {
        'species_name': species_name,
        'trained_date': datetime.now().isoformat(),
        'intronIC_version': '2.0.0',
        'pipeline_architecture': 'single_scaler_v2025',
        'enabled_features': list(enabled_features),  # NEW
        'composite_features': list(composite_features),  # NEW
        'metrics': metrics,
    }
```

**Complexity:** LOW
**Impact:** ~5 lines added
**Note:** Provides backward compatibility detection

---

## Phase 7: CLI & Config Loading

### 7.1 Update CLI Argument Parsing
**File:** `src/intronIC/cli/main.py`

**Changes:**
- Add `--enabled-features` CLI flag
- Load from YAML config if present
- Validate feature combinations
- Pass to scoring/training/prediction

```python
parser.add_argument(
    '--enabled-features',
    nargs='+',
    choices=['five', 'bp', 'three'],
    default=None,
    help='Which motifs to score and use for classification (default: all three)'
)
```

**Complexity:** LOW
**Impact:** ~20 lines added

---

## Phase 8: Backward Compatibility

### 8.1 Detect and Handle Old Models
**File:** `src/intronIC/utils/model_io.py`

**Change:** Assume 3-feature models if metadata missing

```python
def load_model(model_path):
    model = joblib.load(model_path)
    metadata = load_model_metadata(model_path)

    # Backward compatibility: assume all three features if not specified
    enabled_features = metadata.get('enabled_features', ['five', 'bp', 'three'])
    composite_features = metadata.get('composite_features', [
        'min_all', 'absdiff_5_bp', 'absdiff_5_3', 'absdiff_bp_3'
    ])

    return model, enabled_features, composite_features
```

**Complexity:** LOW
**Impact:** ~10 lines added

---

## Phase 9: Testing & Validation

### 9.1 Unit Tests
**New file:** `tests/unit/test_feature_selection.py`

Test cases:
- 1D models (each individual feature)
- 2D models (all pairwise combinations)
- 3D model (current default)
- Composite feature validation
- Config validation
- Transformer dimensionality handling
- Normalizer with missing features

**Complexity:** MODERATE
**Effort:** ~4 hours to write comprehensive tests

### 9.2 Integration Tests
**Test scenarios:**
- Train 1D model on reference data
- Train 2D model (5' + BP)
- Train 2D model (5' + 3', skip BP)
- Verify predictions match expected dimensions
- Test backward compatibility with old models

**Complexity:** MODERATE
**Effort:** ~2 hours

---

## File Modification Summary

### Core Changes (Must Modify)
1. ✅ `src/intronIC/cli/config.py` - Add enabled_features to ScoringConfig
2. ✅ `config/config.yaml` - Add enabled_features section
3. ✅ `src/intronIC/scoring/scorer.py` - Conditional scoring
4. ✅ `src/intronIC/scoring/normalizer.py` - Variable-dimension normalization
5. ✅ `src/intronIC/classification/transformers.py` - Complete rewrite for flexibility
6. ✅ `src/intronIC/classification/trainer.py` - Pass feature config
7. ✅ `src/intronIC/classification/predictor.py` - Extract feature config from model
8. ✅ `src/intronIC/utils/model_io.py` - Save/load feature metadata
9. ✅ `src/intronIC/cli/main.py` - CLI args + config loading

### Affected Files (Update Feature Names/Indices)
10. `src/intronIC/classification/optimizer.py` - Feature matrix construction
11. `src/intronIC/classification/classifier.py` - Model application
12. `src/intronIC/classification/model_inspector.py` - Feature name display
13. `src/intronIC/visualization/plots.py` - Plot feature distributions
14. `src/intronIC/file_io/writers.py` - Output file headers

### Testing Files (New)
15. `tests/unit/test_feature_selection.py`
16. `tests/integration/test_feature_selection_training.py`

**Total files:** 16 files
**Total estimated changes:** ~600-800 lines of code

---

## Risk Assessment

### High Risks ⚠️
1. **Breaking existing models** - Old models won't have feature metadata
   - *Mitigation:* Default to 3-feature assumption, add warnings

2. **Transformer complexity explosion** - Many edge cases
   - *Mitigation:* Comprehensive unit tests for all combinations

3. **Feature matrix misalignment** - Easy to get indices wrong
   - *Mitigation:* Use named feature dictionaries, not positional indexing

### Medium Risks ⚠️
4. **Config frozen dataclass** - Hard to add validation
   - *Mitigation:* Consider switching to pydantic for validation

5. **Composite feature logic errors** - min/max/absdiff combinations
   - *Mitigation:* Test each composite feature independently

### Low Risks ✅
6. **CLI argument parsing** - Well-understood pattern
7. **YAML config loading** - Already have infrastructure

---

## Alternative Approaches

### Option A: Simple Feature Masking (Easier)
Instead of fully dynamic dimensions, always work with 3D space but mask unused features:
- Score all three regions, set disabled ones to 0
- Normalize all three, zero-out disabled features
- Use feature importance weighting instead of removal

**Pros:** Simpler, less code changes
**Cons:** Wasteful scoring, less clean conceptually

### Option B: Multiple Specialized Transformers (More complex)
Create separate transformers for 1D, 2D, 3D cases:
- `SingleFeatureTransformer` (1D → 1D, no composites)
- `TwoFeatureTransformer` (2D → 3D-5D with pairwise composites)
- `ThreeFeatureTransformer` (3D → 4D-11D, current logic)

**Pros:** Clearer separation, easier to test
**Cons:** Code duplication, more maintenance

---

## Recommendation

### Should we implement this?

**Pros:**
- ✅ Enables experimentation with minimal feature sets
- ✅ Potentially faster scoring (skip unused motifs)
- ✅ Could reveal which motifs are most important
- ✅ Matches original intronIC capability

**Cons:**
- ❌ Moderate-high complexity (~600-800 LOC)
- ❌ Significant testing burden (7 feature combinations)
- ❌ Risk of breaking existing functionality
- ❌ Current 3D model works well (99.5%+ accuracy)
- ❌ Limited practical benefit (all three motifs are informative)

### My recommendation: **DEFER UNLESS NEEDED**

**Rationale:**
1. The current 3-feature model achieves excellent performance
2. All three motifs (5'SS, BP, 3'SS) are biologically meaningful
3. The implementation complexity is not trivial
4. Testing matrix grows significantly
5. No clear user demand for this feature yet

**When to revisit:**
- User explicitly requests single-feature training
- Research question requires feature ablation studies
- Performance issues with BP scoring in certain species
- After v2.0 release when architecture is stable

---

## Implementation Time Estimate

If proceeding:
- **Phase 1-2 (Config + Scoring):** 2 hours
- **Phase 3 (Normalization):** 2 hours
- **Phase 4 (Transformer):** 3 hours
- **Phase 5-8 (Training/Prediction/CLI):** 2 hours
- **Phase 9 (Testing):** 6 hours
- **Total:** **15 hours** (2 full work days)

Plus debugging and iteration: **+25% = ~19 hours total**

---

## Questions for Decision

1. **Is there a specific use case?** (e.g., species where BP motif is uninformative)
2. **Research vs production?** (Feature ablation study vs operational need)
3. **Priority vs other work?** (Is this blocking anything?)
4. **Acceptable complexity?** (16 files, ~800 LOC, 19 hours effort)

---

## Conclusion

This is a **feasible but non-trivial** feature addition. The implementation is well-scoped and risks are manageable, but the effort-to-benefit ratio may not justify it unless there's a specific need.

**Recommendation:** Add to backlog, implement post-v2.0 if demand materializes.
