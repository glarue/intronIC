#!/usr/bin/env python
"""Minimal test of GridSearchCV parallelization."""

from multiprocessing import get_all_start_methods, set_start_method
import os

# Set start method first (like original)
fork_types = get_all_start_methods()
if 'forkserver' in fork_types:
    set_start_method('forkserver')
elif 'spawn' in fork_types:
    set_start_method('spawn')

# Set environment variables
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

import numpy as np
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV
import time

# Create fake data
np.random.seed(42)
X = np.random.randn(1000, 3)
y = np.random.randint(0, 2, 1000)

print("Testing GridSearchCV with n_jobs=4...")
print(f"Data: {X.shape[0]} samples, {X.shape[1]} features")

# Test with verbose output
svm = SVC(kernel='linear', class_weight='balanced')
grid = GridSearchCV(
    svm,
    param_grid={'C': np.logspace(-6, 6, 13)},
    cv=5,
    n_jobs=4,
    verbose=2,
    scoring='balanced_accuracy'
)

start = time.time()
grid.fit(X, y)
elapsed = time.time() - start

print(f"\nCompleted in {elapsed:.2f} seconds")
print(f"Best C: {grid.best_params_['C']:.6e}")
print(f"Best score: {grid.best_score_:.4f}")
