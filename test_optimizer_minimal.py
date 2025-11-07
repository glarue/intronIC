#!/usr/bin/env python
"""Minimal test to compare GridSearchCV performance between original and refactored."""

import sys
sys.path.insert(0, '/home/user/intronIC')

# Set multiprocessing start method first
from multiprocessing import get_all_start_methods, set_start_method
fork_types = get_all_start_methods()
if 'forkserver' in fork_types:
    set_start_method('forkserver')
elif 'spawn' in fork_types:
    set_start_method('spawn')

import numpy as np
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV, train_test_split
import time

# Create dataset matching real size (~21k total)
np.random.seed(42)
n_u2 = 20690  # Match real U2 count
n_u12 = 387   # Match real U12 count

# Generate fake z-scores (U12s have positive, U2s have negative)
u12_data = np.random.randn(n_u12, 3) + 2.0  # Mean shifted to +2
u2_data = np.random.randn(n_u2, 3) - 1.0     # Mean shifted to -1

X = np.vstack([u12_data, u2_data])
y = np.array([1]*n_u12 + [0]*n_u2)

print(f"Dataset: {len(X)} samples ({n_u12} U12, {n_u2} U2)")
print(f"Testing with n_jobs=4")
print()

# Test with train_test_split (like original)
print("=" * 60)
print("Test 1: With train_test_split (80/20) - Original's approach")
print("=" * 60)

X_train, X_test, y_train, y_test = train_test_split(
    X, y, train_size=0.80, stratify=y, random_state=42
)
print(f"Training on {len(X_train)} samples (80%)")

C_grid = np.logspace(-2, 3, 6)  # Just 6 C values for speed
print(f"Testing {len(C_grid)} C values: {C_grid}")

svm = SVC(kernel='linear', class_weight='balanced', cache_size=1000, random_state=42)
grid = GridSearchCV(
    svm,
    param_grid={'C': C_grid},
    cv=5,
    scoring='balanced_accuracy',
    n_jobs=4,
    error_score=np.nan
)

start = time.time()
grid.fit(X_train, y_train)
elapsed1 = time.time() - start

print(f"✓ Completed in {elapsed1:.2f} seconds")
print(f"  Best C: {grid.best_params_['C']:.6e}")
print(f"  Best score: {grid.best_score_:.4f}")
print()

# Test without train_test_split (like our buggy version)
print("=" * 60)
print("Test 2: WITHOUT train_test_split - Refactored's old approach")
print("=" * 60)
print(f"Training on {len(X)} samples (100%)")

svm2 = SVC(kernel='linear', class_weight='balanced', cache_size=1000, random_state=42)
grid2 = GridSearchCV(
    svm2,
    param_grid={'C': C_grid},
    cv=5,
    scoring='balanced_accuracy',
    n_jobs=4,
    error_score=np.nan
)

start = time.time()
grid2.fit(X, y)
elapsed2 = time.time() - start

print(f"✓ Completed in {elapsed2:.2f} seconds")
print(f"  Best C: {grid2.best_params_['C']:.6e}")
print(f"  Best score: {grid2.best_score_:.4f}")
print()

# Summary
print("=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"With train_test_split (80%):     {elapsed1:.2f}s")
print(f"Without train_test_split (100%): {elapsed2:.2f}s")
print(f"Speedup from 80/20 split:        {elapsed2/elapsed1:.2f}x")
print()

if elapsed1 < 3.0:
    print("✓ Performance looks good!")
else:
    print("✗ Still too slow - there may be another issue")
