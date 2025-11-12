#!/usr/bin/env python3
"""
Test LinearSVC with CalibratedClassifierCV for probability output.

This shows how to integrate LinearSVC into the refactored code while
maintaining predict_proba() compatibility with the original.
"""

import time
import numpy as np
from sklearn.svm import SVC, LinearSVC
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import GridSearchCV, train_test_split

# Create reference data
np.random.seed(42)
n_u2 = 990
n_u12 = 97
n_features = 3

u2_features = np.random.randn(n_u2, n_features) - 0.5
u12_features = np.random.randn(n_u12, n_features) + 1.5

X_full = np.vstack([u2_features, u12_features])
y_full = np.array([0] * n_u2 + [1] * n_u12)

X_train, X_test, y_train, y_test = train_test_split(
    X_full, y_full,
    train_size=0.80,
    stratify=y_full,
    random_state=42
)

print("=" * 80)
print("LinearSVC WITH CALIBRATION FOR PROBABILITIES")
print("=" * 80)
print(f"Dataset: {len(y_train)} training samples")
print()

# Small grid for quick testing
C_grid = [1e-2, 1.0, 1e2, 1e4, 1e5]

print("=" * 80)
print("APPROACH 1: Calibrate AFTER optimization (RECOMMENDED)")
print("=" * 80)
print("Steps:")
print("1. Use GridSearchCV to find optimal C with LinearSVC (fast)")
print("2. Create new LinearSVC with optimal C")
print("3. Wrap in CalibratedClassifierCV for probabilities")
print()

# Step 1: Optimize LinearSVC without calibration (fast!)
print("[Step 1] Optimizing LinearSVC (no calibration for speed)...")
base_svm = LinearSVC(
    class_weight='balanced',
    dual='auto',
    max_iter=2000,
    random_state=42
)

grid = GridSearchCV(
    base_svm,
    param_grid={'C': C_grid},
    cv=5,
    scoring='balanced_accuracy',
    n_jobs=4,
    verbose=1
)

start = time.time()
grid.fit(X_train, y_train)
opt_time = time.time() - start

best_C = grid.best_params_['C']
print(f"\nOptimization time: {opt_time:.2f}s")
print(f"Best C: {best_C:.2e}")
print(f"Best CV score: {grid.best_score_:.4f}")

# Step 2 & 3: Create final model with optimal C and add calibration
print(f"\n[Step 2-3] Creating final model with C={best_C:.2e} and calibration...")
final_svm = LinearSVC(
    C=best_C,
    class_weight='balanced',
    dual='auto',
    max_iter=2000,
    random_state=42
)

# Wrap in calibration - this is FAST and gives us predict_proba()
calibrated_model = CalibratedClassifierCV(
    final_svm,
    method='sigmoid',  # 'sigmoid' is faster than 'isotonic'
    cv=3  # Use 3 folds for calibration (faster than 5)
)

start = time.time()
calibrated_model.fit(X_train, y_train)
calib_time = time.time() - start

print(f"Calibration time: {calib_time:.2f}s")
print(f"Total time: {opt_time + calib_time:.2f}s")

# Test probability output
print(f"\n[Testing] Probability predictions:")
test_sample = X_test[:5]
proba = calibrated_model.predict_proba(test_sample)
pred = calibrated_model.predict(test_sample)
print(f"Predictions: {pred}")
print(f"Probabilities (U2, U12):")
print(proba)

print("\n" + "=" * 80)
print("APPROACH 2: Calibrate DURING optimization (SLOWER)")
print("=" * 80)
print("Use GridSearchCV with already-calibrated model")
print("(Slower because calibration happens for every C value)")
print()

# Pre-wrap in calibration
precalib_base = CalibratedClassifierCV(
    LinearSVC(class_weight='balanced', dual='auto', max_iter=2000, random_state=42),
    method='sigmoid',
    cv=3
)

# Note: Parameter names are now nested!
grid2 = GridSearchCV(
    precalib_base,
    param_grid={'estimator__C': C_grid},  # Note: sklearn 0.24+ uses 'estimator__' prefix
    cv=5,
    scoring='balanced_accuracy',
    n_jobs=4,
    verbose=1
)

print("Optimizing pre-calibrated LinearSVC...")
start = time.time()
grid2.fit(X_train, y_train)
precalib_time = time.time() - start

print(f"\nTime: {precalib_time:.2f}s")
print(f"Best C: {grid2.best_params_['estimator__C']:.2e}")
print(f"Best CV score: {grid2.best_score_:.4f}")

print("\n" + "=" * 80)
print("COMPARISON WITH SVC")
print("=" * 80)

# For comparison: SVC with probability=True
svc_model = SVC(
    kernel='linear',
    probability=True,
    class_weight='balanced',
    cache_size=1000,
    random_state=42
)

grid3 = GridSearchCV(
    svc_model,
    param_grid={'C': C_grid},
    cv=5,
    scoring='balanced_accuracy',
    n_jobs=4,
    verbose=1
)

print("Optimizing SVC with probability=True...")
start = time.time()
grid3.fit(X_train, y_train)
svc_time = time.time() - start

print(f"\nTime: {svc_time:.2f}s")
print(f"Best C: {grid3.best_params_['C']:.2e}")
print(f"Best CV score: {grid3.best_score_:.4f}")

# Compare probabilities
print("\n[Comparing probability outputs]")
proba_svc = grid3.best_estimator_.predict_proba(test_sample)
proba_linear = calibrated_model.predict_proba(test_sample)

print("\nSVC probabilities:")
print(proba_svc)
print("\nLinearSVC+calibration probabilities:")
print(proba_linear)
print("\nProbability differences (should be small):")
print(np.abs(proba_svc - proba_linear))

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print(f"Approach 1 (optimize then calibrate):  {opt_time + calib_time:6.2f}s")
print(f"Approach 2 (optimize pre-calibrated):  {precalib_time:6.2f}s")
print(f"SVC with probability=True:             {svc_time:6.2f}s")
print()
print(f"Speedup Approach 1 vs SVC: {svc_time/(opt_time + calib_time):.2f}x")
print(f"Speedup Approach 2 vs SVC: {svc_time/precalib_time:.2f}x")
print()
print("RECOMMENDATION:")
print("Use Approach 1 (optimize then calibrate) for best performance.")
print("This is the pattern to use in optimizer.py!")
