#!/usr/bin/env python3
"""
Examine GridSearchCV behavior in detail.

This test captures verbose output and CV results to understand:
1. Are both versions actually testing the same C values?
2. Are convergence warnings appearing?
3. What are the actual fit times per parameter?
4. Are there any NaN scores or failed fits?
"""

import numpy as np
import pandas as pd
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.base import clone

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
print("GRIDSEARCH BEHAVIOR ANALYSIS")
print("=" * 80)
print(f"Dataset: {len(y_train)} samples ({np.sum(y_train==0)} U2, {np.sum(y_train==1)} U12)")
print()

# Test with a range that includes problematic high C values
# But keep it small for reasonable runtime
C_grid = [1e-2, 1e0, 1e2, 1e4, 1e5]

print("=" * 80)
print("TEST 1: Original Style")
print("=" * 80)

model1 = SVC(probability=True, kernel='linear', cache_size=1000, random_state=42)
model1.class_weight = 'balanced'

print("\nModel configuration BEFORE GridSearchCV:")
print(f"  kernel: {model1.kernel}")
print(f"  probability: {model1.probability}")
print(f"  class_weight: {model1.class_weight}")
print(f"  cache_size: {model1.cache_size}")

print("\nChecking clone behavior:")
cloned = clone(model1)
print(f"  Cloned kernel: {cloned.kernel}")
print(f"  Cloned probability: {cloned.probability}")
print(f"  Cloned class_weight: {getattr(cloned, 'class_weight', 'NOT PRESERVED')}")
print(f"  Cloned cache_size: {cloned.cache_size}")

grid1 = GridSearchCV(
    model1,
    param_grid={'C': C_grid},
    cv=5,
    scoring='balanced_accuracy',
    n_jobs=4,
    verbose=2,
    error_score='raise'  # Fail loudly if something goes wrong
)

print(f"\nRunning GridSearchCV with C values: {C_grid}")
print("=" * 80)
grid1.fit(X_train, y_train)
print("=" * 80)

print("\nGridSearchCV Results:")
results1_df = pd.DataFrame(grid1.cv_results_)
print(results1_df[['param_C', 'mean_test_score', 'std_test_score', 'mean_fit_time', 'std_fit_time']])

print(f"\nBest parameters: {grid1.best_params_}")
print(f"Best score: {grid1.best_score_:.4f}")

# Check for any NaN scores
nan_count = results1_df['mean_test_score'].isna().sum()
print(f"Failed fits (NaN scores): {nan_count}")

print("\n" + "=" * 80)
print("TEST 2: Refactored Style")
print("=" * 80)

model2 = SVC(kernel='linear', class_weight='balanced', cache_size=1000, random_state=42)

print("\nModel configuration BEFORE GridSearchCV:")
print(f"  kernel: {model2.kernel}")
print(f"  probability: {model2.probability}")
print(f"  class_weight: {model2.class_weight}")
print(f"  cache_size: {model2.cache_size}")

print("\nChecking clone behavior:")
cloned2 = clone(model2)
print(f"  Cloned kernel: {cloned2.kernel}")
print(f"  Cloned probability: {cloned2.probability}")
print(f"  Cloned class_weight: {cloned2.class_weight}")
print(f"  Cloned cache_size: {cloned2.cache_size}")

grid2 = GridSearchCV(
    model2,
    param_grid={'C': C_grid},
    cv=5,
    scoring='balanced_accuracy',
    n_jobs=4,
    verbose=2,
    error_score='raise'
)

print(f"\nRunning GridSearchCV with C values: {C_grid}")
print("=" * 80)
grid2.fit(X_train, y_train)
print("=" * 80)

print("\nGridSearchCV Results:")
results2_df = pd.DataFrame(grid2.cv_results_)
print(results2_df[['param_C', 'mean_test_score', 'std_test_score', 'mean_fit_time', 'std_fit_time']])

print(f"\nBest parameters: {grid2.best_params_}")
print(f"Best score: {grid2.best_score_:.4f}")

# Check for any NaN scores
nan_count2 = results2_df['mean_test_score'].isna().sum()
print(f"Failed fits (NaN scores): {nan_count2}")

print("\n" + "=" * 80)
print("COMPARISON")
print("=" * 80)

comparison = pd.DataFrame({
    'C': C_grid,
    'Original_time': results1_df['mean_fit_time'].values,
    'Refactored_time': results2_df['mean_fit_time'].values,
    'Original_score': results1_df['mean_test_score'].values,
    'Refactored_score': results2_df['mean_test_score'].values,
})
comparison['Time_ratio'] = comparison['Refactored_time'] / comparison['Original_time']
comparison['Score_diff'] = comparison['Refactored_score'] - comparison['Original_score']

print(comparison.to_string(index=False))

print("\n" + "=" * 80)
print("KEY OBSERVATIONS")
print("=" * 80)
print(f"1. Clone preserves class_weight: Original={getattr(cloned, 'class_weight', 'NO')}, Refactored={cloned2.class_weight}")
print(f"2. Scores are similar: Max diff = {comparison['Score_diff'].abs().max():.6f}")
print(f"3. Time ratio at C=1e5: {comparison[comparison['C']==1e5]['Time_ratio'].values[0]:.2f}x")
print(f"4. Time increases with C: {comparison['Time_ratio'].values}")
print()
print("If time_ratio increases dramatically at high C, that's the smoking gun!")
