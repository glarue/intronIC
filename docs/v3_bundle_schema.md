# FI v3 model bundle schema

intronIC v2.4 introduces support for v3 model bundles — a single
self-describing pickle that pairs the trained ensemble with the
configuration and training metadata that produced it.

The v3 schema is the recommended format for new models. Existing v2.3
bundles continue to load unchanged.

## Top-level structure

```python
{
    "version":          "v3",
    "model_id":         str,        # e.g. "fi_v3_C200_g0.001_ezf0.75_seeds42-100-200"
    "config":           dict,       # see below
    "training":         dict,       # see below
    "holdout_metrics":  dict,       # {seed_id: {"f1@50": ..., "precision@50": ..., ...}}
    "v23_baseline":     dict,       # same shape as a per-seed entry
    "seeds":            dict,       # {seed_id: {"models": [42 sub-models]}}
}
```

### `config`

Identifies the model architecture and tells the loader which input
features to expect.

```python
{
    "C":                  float,    # SVC penalty (e.g. 200.0)
    "gamma":              float,    # RBF gamma (e.g. 0.001)
    "kernel":             "rbf",
    "calibration_method": "isotonic" | "sigmoid",
    "easy_fraction":      float,    # subsample ratio used at training (display only)
    "ensemble_size":      int,      # sub-models per seed (typically 42)
    "seeds":              [int, ...],
    "input_features":     [str, ...],   # IntronScores attribute names; first 3 are always
                                        # the base z-scores (five_z_score, bp_z_score, three_z_score)
    "feature_subtypes":   [str, ...],   # e.g. ["forceinclude_eb"] — informational
    # ... other keys are tolerated and preserved
}
```

The loader uses `config["input_features"][3:]` as `extra_features` on
every constructed `SVMModel`. Production reads each extra feature by
attribute name from `IntronScores` (via `getattr(scores, name, None)`),
so any feature you list there must be a real attribute or `@property`
on `IntronScores`.

### `training`

Records what the SVM was trained against. Used to compute
`training_prior` and to log provenance.

```python
{
    "corpus_rows":      int,    # full training corpus size (pre-subsample)
    "u12_positives":    int,    # corpus-wide U12 count (includes any held-out positives)
    "test_positives":   int,    # held-out positives (subtracted from u12_positives
                                #   when computing training_prior)
    "n_train":          int,    # rows actually used for training (per seed)
    "n_test":           int,
    "n_species":        int,
    "n_clades":         int,
    "trained_date":     str,    # ISO date
    "notes":            str,    # free-form provenance
}
```

### `seeds`

Maps each random seed to its trained sub-models. Each sub-model is a
fitted estimator with `predict_proba(X) -> ndarray` (typically a
`sklearn.calibration.CalibratedClassifierCV`).

```python
{
    42:  {"models": [<sub_model>, <sub_model>, ...]},
    100: {"models": [...]},
    200: {"models": [...]},
}
```

## Loading

`intronIC.utils.model_io.load_model()` returns the raw pickled object.
Pass it through `normalize_model_bundle()` to get the runtime shape:

```python
from intronIC.utils.model_io import load_model, normalize_model_bundle

raw = load_model(path)            # {"version": "v3", ...}  (or any legacy format)
runtime = normalize_model_bundle(raw)
ensemble = runtime["ensemble"]    # SVMEnsemble of 126 SVMModel wrappers
```

`normalize_model_bundle()` is idempotent: v2.3 dicts and bare
`SVMEnsemble` objects are returned unchanged. Both `cli/main.py` load
points already pipe through it, so adding new v3 features generally
doesn't require changes to consumers.

## Runtime fields synthesized from a v3 bundle

| Runtime key | v3 source | Notes |
|---|---|---|
| `ensemble` | `seeds` × `config` | `SVMEnsemble[SVMModel]`, one wrapper per sub-model |
| `normalizer` | (none) | `None` — v3 features are already z-scored, so production's adaptive mode is correct |
| `threshold` | `V3_DEFAULT_THRESHOLD` | 50.0 (matches the [0,100] score scale) |
| `training_prior` | `(u12_positives - test_positives) / n_train` | Used by `--species_prior` adjustment |
| `human_negative_stats` | (none) | `None` — currently unused by the classify pipeline |
| `_v3_bundle` | the original bundle | Preserved for downstream provenance/logging |

## Authoring a v3 bundle

The minimum a producer must populate is:

- `version` set to `"v3"`
- `config["seeds"]` listing the seed ids
- `config["input_features"]` listing the IntronScores attribute names in
  the same order they were stacked at training time, with the three base
  z-scores first
- `config["C"]`, `config["gamma"]`, `config["kernel"]`,
  `config["calibration_method"]`
- `seeds[seed_id]["models"]` — the trained sub-models for that seed
- `training["n_train"]` and `training["u12_positives"]` (so
  `training_prior` is meaningful)

Everything else is optional but encouraged for provenance.

## Compatibility

| Format | Status | Read at | Notes |
|---|---|---|---|
| bare `SVMEnsemble` | supported | both load points | Legacy backward-compat |
| v2.3 dict | supported | both load points | Pass-through |
| v3 dict | supported (v2.4+) | both load points | Translated by `normalize_model_bundle` |
