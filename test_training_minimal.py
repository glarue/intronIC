#!/usr/bin/env python3
"""
Minimal test to understand original vs refactored SVM training behavior.

This test isolates JUST the training step with identical data to see:
1. What parameters are actually used during GridSearchCV
2. How clone() behaves with different initialization methods
3. Actual timing differences with minimal overhead
"""

import time
import numpy as np
from sklearn.svm import SVC, LinearSVC
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.base import clone
from sklearn.calibration import CalibratedClassifierCV

# Create reference data matching your small dataset
np.random.seed(42)
n_u2 = 990
n_u12 = 97
n_features = 3

# Generate synthetic data with similar properties
u2_features = np.random.randn(n_u2, n_features) - 0.5  # Centered around -0.5
u12_features = np.random.randn(n_u12, n_features) + 1.5  # Centered around +1.5

X_full = np.vstack([u2_features, u12_features])
y_full = np.array([0] * n_u2 + [1] * n_u12)

# 80% split for training (matching both implementations)
X_train, X_test, y_train, y_test = train_test_split(
    X_full, y_full,
    train_size=0.80,
    stratify=y_full,
    random_state=42
)

print("=" * 80)
print("MINIMAL SVM TRAINING TEST")
print("=" * 80)
print(f"Dataset: {len(y_full)} samples ({n_u2} U2, {n_u12} U12)")
print(f"Training: {len(y_train)} samples ({np.sum(y_train==0)} U2, {np.sum(y_train==1)} U12)")
print(f"Testing: {len(y_test)} samples ({np.sum(y_test==0)} U2, {np.sum(y_test==1)} U12)")
print()

# Small C grid for fast testing (just 5 values spanning range)
C_grid_small = [1e-4, 1e-2, 1.0, 100.0, 10000.0]

print("=" * 80)
print("TEST 1: Original Style - probability=True, class_weight set as attribute")
print("=" * 80)

# Exactly replicate original approach
model1 = SVC(probability=True, kernel='linear', cache_size=1000, random_state=42)
model1.class_weight = 'balanced'

# Check what clone sees
clone1 = clone(model1)
print(f"After clone:")
print(f"  - probability: {getattr(clone1, 'probability', 'NOT SET')}")
print(f"  - class_weight: {getattr(clone1, 'class_weight', 'NOT SET')}")
print(f"  - kernel: {getattr(clone1, 'kernel', 'NOT SET')}")
print()

grid1 = GridSearchCV(
    model1,
    param_grid={'C': C_grid_small},
    cv=5,
    scoring='balanced_accuracy',
    n_jobs=4,
    verbose=1
)

print(f"Starting GridSearchCV with {len(C_grid_small)} C values, 5 folds, 4 jobs...")
start1 = time.time()
grid1.fit(X_train, y_train)
elapsed1 = time.time() - start1

print(f"\nElapsed: {elapsed1:.2f} seconds")
print(f"Best C: {grid1.best_params_['C']:.2e}")
print(f"Best score: {grid1.best_score_:.4f}")
print()

print("=" * 80)
print("TEST 2: Refactored Style - NO probability, class_weight in constructor")
print("=" * 80)

model2 = SVC(
    kernel='linear',
    class_weight='balanced',
    cache_size=1000,
    random_state=42
)

# Check what clone sees
clone2 = clone(model2)
print(f"After clone:")
print(f"  - probability: {getattr(clone2, 'probability', 'NOT SET')}")
print(f"  - class_weight: {getattr(clone2, 'class_weight', 'NOT SET')}")
print(f"  - kernel: {getattr(clone2, 'kernel', 'NOT SET')}")
print()

grid2 = GridSearchCV(
    model2,
    param_grid={'C': C_grid_small},
    cv=5,
    scoring='balanced_accuracy',
    n_jobs=4,
    verbose=1
)

print(f"Starting GridSearchCV with {len(C_grid_small)} C values, 5 folds, 4 jobs...")
start2 = time.time()
grid2.fit(X_train, y_train)
elapsed2 = time.time() - start2

print(f"\nElapsed: {elapsed2:.2f} seconds")
print(f"Best C: {grid2.best_params_['C']:.2e}")
print(f"Best score: {grid2.best_score_:.4f}")
print()

print("=" * 80)
print("TEST 3: LinearSVC (modern approach for linear kernels)")
print("=" * 80)

# LinearSVC doesn't support probability directly, so wrap in calibration
base_svm = LinearSVC(
    class_weight='balanced',
    dual='auto',
    max_iter=2000,
    random_state=42
)

# Wrap in calibration to get predict_proba (if needed)
# Note: We use method='sigmoid' which is faster than 'isotonic'
model3 = CalibratedClassifierCV(base_svm, method='sigmoid', cv=3)

grid3 = GridSearchCV(
    model3,
    param_grid={'base_estimator__C': C_grid_small},  # Note the nested parameter
    cv=5,
    scoring='balanced_accuracy',
    n_jobs=4,
    verbose=1
)

print(f"Starting GridSearchCV with {len(C_grid_small)} C values, 5 folds, 4 jobs...")
start3 = time.time()
grid3.fit(X_train, y_train)
elapsed3 = time.time() - start3

print(f"\nElapsed: {elapsed3:.2f} seconds")
print(f"Best C: {grid3.best_params_['base_estimator__C']:.2e}")
print(f"Best score: {grid3.best_score_:.4f}")
print()

print("=" * 80)
print("TEST 4: LinearSVC without calibration (fastest)")
print("=" * 80)

model4 = LinearSVC(
    class_weight='balanced',
    dual='auto',
    max_iter=2000,
    random_state=42
)

grid4 = GridSearchCV(
    model4,
    param_grid={'C': C_grid_small},
    cv=5,
    scoring='balanced_accuracy',
    n_jobs=4,
    verbose=1
)

print(f"Starting GridSearchCV with {len(C_grid_small)} C values, 5 folds, 4 jobs...")
start4 = time.time()
grid4.fit(X_train, y_train)
elapsed4 = time.time() - start4

print(f"\nElapsed: {elapsed4:.2f} seconds")
print(f"Best C: {grid4.best_params_['C']:.2e}")
print(f"Best score: {grid4.best_score_:.4f}")
print()

print("=" * 80)
print("SUMMARY")
print("=" * 80)
print(f"Test 1 (Original style, probability=True):        {elapsed1:6.2f}s")
print(f"Test 2 (Refactored style, no probability):        {elapsed2:6.2f}s")
print(f"Test 3 (LinearSVC + calibration):                 {elapsed3:6.2f}s")
print(f"Test 4 (LinearSVC alone):                         {elapsed4:6.2f}s")
print()
print(f"Speedup Test1 vs Test2: {elapsed1/elapsed2:5.2f}x {'(Original FASTER)' if elapsed1 < elapsed2 else '(Refactored FASTER)'}")
print(f"Speedup Test1 vs Test3: {elapsed1/elapsed3:5.2f}x")
print(f"Speedup Test1 vs Test4: {elapsed1/elapsed4:5.2f}x")
print(f"Speedup Test2 vs Test4: {elapsed2/elapsed4:5.2f}x")
print()

# Test if predict_proba works
print("=" * 80)
print("PROBABILITY PREDICTION TEST")
print("=" * 80)
print("Test 1 (probability=True):", "predict_proba" in dir(grid1.best_estimator_))
print("Test 2 (no probability):", "predict_proba" in dir(grid2.best_estimator_))
print("Test 3 (LinearSVC + calib):", "predict_proba" in dir(grid3.best_estimator_))
print("Test 4 (LinearSVC alone):", "predict_proba" in dir(grid4.best_estimator_))

if "predict_proba" in dir(grid1.best_estimator_):
    proba1 = grid1.best_estimator_.predict_proba(X_test[:5])
    print(f"\nTest 1 probabilities sample:\n{proba1}")

if "predict_proba" in dir(grid2.best_estimator_):
    proba2 = grid2.best_estimator_.predict_proba(X_test[:5])
    print(f"\nTest 2 probabilities sample:\n{proba2}")
else:
    print("\nTest 2: No predict_proba available")

if "predict_proba" in dir(grid3.best_estimator_):
    proba3 = grid3.best_estimator_.predict_proba(X_test[:5])
    print(f"\nTest 3 probabilities sample:\n{proba3}")

print()
print("=" * 80)
print("Key Insights:")
print("=" * 80)
print("1. If Test1 and Test2 have similar times, probability=True is NOT the issue")
print("2. If Test4 is much faster, LinearSVC is the right solution")
print("3. If Test3 is slower than Test4, calibration overhead is significant")
print("4. Check clone() behavior above - does it preserve class_weight?")
