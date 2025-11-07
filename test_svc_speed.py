#!/usr/bin/env python3
"""Test SVC performance with different parameter combinations."""

import time
import numpy as np
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV

# Create small dataset matching our reference size
np.random.seed(42)
n_samples = 870  # 80% of 1087
n_features = 3

X = np.random.randn(n_samples, n_features)
y = np.array([0] * 792 + [1] * 78)  # 90% U2, 10% U12

C_grid = [1e-6, 1e-4, 1e-2, 1.0, 100.0]  # Small grid for quick test

print("Testing SVC with different parameters...")
print(f"Dataset: {n_samples} samples, {n_features} features")
print(f"C grid: {C_grid}\n")

# Test 1: Original-style (probability=True, no class_weight in constructor)
print("=" * 60)
print("Test 1: Original style")
print("SVC(probability=True, kernel='linear', cache_size=1000)")
print("=" * 60)
model1 = SVC(probability=True, kernel='linear', cache_size=1000)
model1.class_weight = 'balanced'

grid1 = GridSearchCV(
    model1,
    param_grid={'C': C_grid},
    cv=5,
    scoring='balanced_accuracy',
    n_jobs=4,
    verbose=2
)

start = time.time()
grid1.fit(X, y)
elapsed1 = time.time() - start
print(f"\nElapsed: {elapsed1:.2f} seconds\n")

# Test 2: Refactored style (no probability, class_weight in constructor)
print("=" * 60)
print("Test 2: Refactored style")
print("SVC(kernel='linear', class_weight='balanced', cache_size=1000)")
print("=" * 60)
model2 = SVC(kernel='linear', class_weight='balanced', cache_size=1000, random_state=42)

grid2 = GridSearchCV(
    model2,
    param_grid={'C': C_grid},
    cv=5,
    scoring='balanced_accuracy',
    n_jobs=4,
    verbose=2
)

start = time.time()
grid2.fit(X, y)
elapsed2 = time.time() - start
print(f"\nElapsed: {elapsed2:.2f} seconds\n")

# Test 3: Like original but WITHOUT probability
print("=" * 60)
print("Test 3: Original WITHOUT probability=True")
print("SVC(kernel='linear', cache_size=1000)")
print("=" * 60)
model3 = SVC(kernel='linear', cache_size=1000)
model3.class_weight = 'balanced'

grid3 = GridSearchCV(
    model3,
    param_grid={'C': C_grid},
    cv=5,
    scoring='balanced_accuracy',
    n_jobs=4,
    verbose=2
)

start = time.time()
grid3.fit(X, y)
elapsed3 = time.time() - start
print(f"\nElapsed: {elapsed3:.2f} seconds\n")

print("=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"Test 1 (Original, probability=True):     {elapsed1:.2f}s")
print(f"Test 2 (Refactored, no probability):      {elapsed2:.2f}s")
print(f"Test 3 (Original, no probability):        {elapsed3:.2f}s")
print(f"\nDifference Test1 vs Test2: {elapsed2 - elapsed1:+.2f}s ({100*(elapsed2/elapsed1-1):+.1f}%)")
print(f"Difference Test1 vs Test3: {elapsed3 - elapsed1:+.2f}s ({100*(elapsed3/elapsed1-1):+.1f}%)")
