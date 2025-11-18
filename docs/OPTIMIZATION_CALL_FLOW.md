# SVM Optimization Call Flow Comparison

## ORIGINAL intronIC.py

### Entry Point: main() → line 3785-3795
```python
model = svm.SVC(probability=True, kernel='linear', cache_size=1000)
model, model_performance = optimize_svm(
    model,                    # base_model WITH probability=True
    ref_u12_vector,          # U12 z-scores
    ref_u2_vector,           # U2 z-scores
    n_optimize=OPTIMIZE_N,   # default: 5
    iterations=SVM_ITER,     # default: 0
    seed=SEED,
    cv_jobs=CV_PROCESSES     # from -p argument
)
```

**Key**: Base model created ONCE with probability=True, passed to optimize_svm

---

### optimize_svm() → line 5431-5523

#### Parameters:
- `base_model`: The SVC(probability=True...) instance
- `iterations`: Default 0
- `cv_jobs`: Parallelization level

#### Setup (lines 5446-5463):
```python
# Initial C grid: 10^-6 to 10^6 (13 points)
log_intervals = np.logspace(-6, 6, num=13)
initial_parameters = {'C': log_intervals}

cv_params = {
    'cv': 5,
    'scoring': scorer,           # 'balanced_accuracy'
    'n_jobs': cv_jobs,           # e.g., 4
    'param_grid': initial_parameters
}
```

#### Optimization Loop (lines 5465-5503):
```python
for search_round in range(1, n_optimize + 1):  # 1-5
    round_seed = seed + search_round

    # CRITICAL CALL: Passes cv_params to train_svm
    search_model, performance = train_svm(
        base_model,              # SVC(probability=True...)
        u2_vector,               # All U2 z-scores
        u12_vector,              # All U12 z-scores
        iterations=iterations,   # 0
        cv_parameters=cv_params, # GridSearchCV params
        seed=round_seed
    )

    # Get best C from this round
    best_first_params = rank1_param_avg(search_model)

    # Refine grid around best C for next round
    for p, v in best_first_params.items():
        best_index = index_of_nearest(refined_params[p], v)
        low_bound = refined_params[p][max(best_index - 1, 0)]
        high_bound = refined_params[p][min(best_index + 1, len(...) - 1)]
        p_range = np.geomspace(low_bound, high_bound, 100)  # 100 points!
        refined_params[p] = p_range

    cv_params['param_grid'] = refined_params
```

**Key**: Each round refines to 100 C values around best from previous round

---

### train_svm() → line 5345-5428

#### When cv_parameters is NOT None (during optimization):

```python
# Line 5353-5354: Convert to numpy arrays
u12s = np.array(u12_vector)
u2s = np.array(u2_vector)

# Line 5356-5372: Check iterations (0 in our case)
if iterations:
    # Would subsample U2s multiple times - NOT TAKEN
    ...
else:
    # Line 5369: CRITICAL - Sets class_weight here!
    base_model.class_weight = 'balanced'
    u2_sets = [u2s]  # Single set with all U2s
    zero_labels = np.zeros(len(u2s))

# Line 5380: Create combined labels
input_labels = np.concatenate([zero_labels, np.ones(len(u12s))])

# Line 5381-5406: Single iteration (iterations=0)
for iter_n, u2_subset in enumerate(u2_sets, start=1):

    # Line 5388-5393: Create GridSearchCV
    model = GridSearchCV(
        base_model,           # SVC(probability=True, class_weight='balanced', ...)
        error_score=np.nan,
        **cv_parameters       # cv=5, n_jobs=4, param_grid={'C': [...]}, scoring=...
    )

    # Line 5397: Concatenate U2 and U12 features
    input_feats = np.concatenate([u2_subset, u12s])

    # Line 5398-5404: CRITICAL - train_test_split BEFORE GridSearchCV
    (train_scores, test_scores,
     train_labels, test_labels) = train_test_split(
        input_feats,
        input_labels,
        train_size=0.80,          # 80% for GridSearchCV
        stratify=input_labels,
        random_state=seed
    )

    # Line 5406: Fit GridSearchCV on 80% of data
    model.fit(train_scores, train_labels)

    # Lines 5407-5422: Calculate metrics on test set
    predict_labels = model.predict(test_scores)
    subset_f1 = f1_score(test_labels, predict_labels)
    ...
```

**CRITICAL FINDINGS**:
1. GridSearchCV receives `SVC(probability=True, class_weight='balanced', cache_size=1000)`
2. train_test_split reduces dataset to 80% BEFORE GridSearchCV.fit()
3. With 1087 samples: 870 for GridSearchCV, 217 for final test

---

## REFACTORED intronIC_refactored/

### Entry Point: cli/main.py → normalize_scores() → classify_introns()

```python
# Line ~615 in main.py
ensemble = classifier.train_and_classify(
    u12_ref_introns,
    u2_ref_introns,
    experimental_introns,
    ...
)
```

---

### SVMClassifier.train_and_classify() → classification/classifier.py

```python
def train_and_classify(self, u12_refs, u2_refs, experimental):
    # 1. Optimize hyperparameters
    parameters = self.optimizer.optimize(u12_refs, u2_refs)

    # 2. Train ensemble
    ensemble = self.trainer.train_ensemble(u12_refs, u2_refs, parameters)

    # 3. Classify
    classified = self.predictor.predict(ensemble, experimental)

    return ensemble
```

---

### SVMOptimizer.optimize() → classification/optimizer.py

#### Lines 94-142:
```python
def optimize(self, u12_introns, u2_introns, initial_range=(1e-6, 1e6)):
    # Extract features
    X, y = self._prepare_training_data(u12_introns, u2_introns)

    # Initialize search range
    current_grid = self._create_initial_grid(initial_range)  # 13 points

    # Run geometric refinement
    for round_idx in range(self.n_rounds):  # 0-4
        print(f"Optimization round {round_idx + 1}/{self.n_rounds}...")

        round_result = self._grid_search_round(
            X, y, current_grid, round_idx
        )
        self.rounds_.append(round_result)

        # Prepare next round's grid
        if round_idx < self.n_rounds - 1:
            current_grid = self._refine_grid(current_grid, round_result.best_C)

    # Final parameter
    final_C = gmean(self.rounds_[-1].rank_one_Cs)
    ...
```

---

### SVMOptimizer._grid_search_round() → classification/optimizer.py

#### Lines 220-284 (CURRENT VERSION WITH MY CHANGES):
```python
def _grid_search_round(self, X, y, C_grid, round_idx):
    # Line 243-247: Create base model (NO probability parameter!)
    base_svm = SVC(
        kernel='linear',
        cache_size=1000,
        random_state=self.random_state + round_idx
    )

    # Lines 253-258: train_test_split (MY RECENT ADDITION)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y,
        train_size=0.80,
        stratify=y,
        random_state=self.random_state + round_idx
    )

    # Lines 260-268: Create GridSearchCV
    grid_search = GridSearchCV(
        base_svm,
        param_grid={'C': C_grid},
        cv=self.cv_folds,               # 5
        scoring='balanced_accuracy',
        n_jobs=self.n_jobs,              # e.g., 4
        error_score=np.nan,
        verbose=2
    )

    # Line 271: Fit on 80% subset
    grid_search.fit(X_train, y_train)
    ...
```

**CRITICAL DIFFERENCES FROM ORIGINAL**:
1. NO `probability=True` in SVC
2. NO `class_weight='balanced'` in SVC constructor
3. Uses train_test_split (I just added this)
4. Grid size: 13 on round 1, then ???

---

### _refine_grid() → classification/optimizer.py line 286

```python
def _refine_grid(self, current_grid, best_C):
    best_index = np.argmin(np.abs(current_grid - best_C))

    if best_index == 0:
        low_bound = current_grid[0]
        high_bound = current_grid[1]
    elif best_index == len(current_grid) - 1:
        low_bound = current_grid[-2]
        high_bound = current_grid[-1]
    else:
        low_bound = current_grid[best_index - 1]
        high_bound = current_grid[best_index + 1]

    # Line 307: CRITICAL - How many points in refined grid???
    return np.geomspace(low_bound, high_bound, num=self.n_points_refine)
```

**CHECK**: What is `self.n_points_refine`?

Line 67 in __init__:
```python
def __init__(
    self,
    n_rounds: int = 5,
    n_points_initial: int = 13,
    n_points_refine: int = 100,  # <--- HERE!
    ...
):
```

So rounds 2-5 each test **100 C values**!

---

## SUMMARY OF KEY DIFFERENCES

| Aspect | Original | Refactored (before my changes) | Refactored (after my changes) |
|--------|----------|-------------------------------|-------------------------------|
| Base SVC | `SVC(probability=True, kernel='linear', cache_size=1000)` | `SVC(kernel='linear', cache_size=1000, random_state=...)` | Same |
| class_weight | Set via `base_model.class_weight='balanced'` before GridSearchCV | NOT SET AT ALL! | Still NOT SET |
| train_test_split | Yes, 80/20 before GridSearchCV.fit() | NO | YES (I just added) |
| Grid size round 1 | 13 C values | 13 C values | 13 C values |
| Grid size rounds 2-5 | 100 C values each | 100 C values each | 100 C values each |
| Total CV fits | 5 * (13 + 100*4) = 2,065 fits | Same | Same |

---

## THE SMOKING GUN

**The refactored version is MISSING `class_weight='balanced'`!**

Looking at the original train_svm line 5369:
```python
base_model.class_weight = 'balanced'
```

This is set BEFORE GridSearchCV is created, so every SVC instance used by GridSearchCV has class_weight='balanced'.

In the refactored version, we create a fresh SVC in _grid_search_round() WITHOUT setting class_weight, and we NEVER set it afterward!

This means:
- Original: GridSearchCV uses SVC instances with balanced class weights
- Refactored: GridSearchCV uses SVC instances with NO class weighting (treats 990 U2s and 97 U12s equally)

**With severely imbalanced classes (10:1 ratio), NOT using class_weight='balanced' will cause SVM optimization to behave completely differently and likely much slower!**

