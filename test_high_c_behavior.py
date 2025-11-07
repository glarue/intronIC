#!/usr/bin/env python3
"""
Test SVM training behavior at HIGH C values.

The performance data shows:
- C <= 10,000: Both versions perform similarly
- C = 100,000: Refactored is 17-30s per fold (slow!)
- C = 1,000,000: Refactored is 1.7-1.9 minutes per fold (very slow!)

This test focuses on high C values to see where the divergence happens.
"""

import time
import numpy as np
from sklearn.svm import SVC, LinearSVC
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.calibration import CalibratedClassifierCV

# Create reference data matching your small dataset
np.random.seed(42)
n_u2 = 990
n_u12 = 97
n_features = 3

# Generate synthetic data
u2_features = np.random.randn(n_u2, n_features) - 0.5
u12_features = np.random.randn(n_u12, n_features) + 1.5

X_full = np.vstack([u2_features, u12_features])
y_full = np.array([0] * n_u2 + [1] * n_u12)

# 80% split
X_train, X_test, y_train, y_test = train_test_split(
    X_full, y_full,
    train_size=0.80,
    stratify=y_full,
    random_state=42
)

print("=" * 80)
print("HIGH C VALUE PERFORMANCE TEST")
print("=" * 80)
print(f"Dataset: {len(y_train)} training samples ({np.sum(y_train==0)} U2, {np.sum(y_train==1)} U12)")
print()

# Test progressively higher C values
C_values = [1e2, 1e3, 1e4, 1e5, 1e6]
results = []

for C in C_values:
    print("=" * 80)
    print(f"Testing C = {C:.0e}")
    print("=" * 80)

    # Test 1: Original style (SVC with probability=True, class_weight as attribute)
    print("\n[1] Original style: SVC(probability=True) + class_weight attribute")
    model1 = SVC(probability=True, kernel='linear', cache_size=1000, random_state=42, C=C)
    model1.class_weight = 'balanced'

    start = time.time()
    scores1 = cross_val_score(model1, X_train, y_train, cv=5, scoring='balanced_accuracy', n_jobs=1)
    elapsed1 = time.time() - start
    print(f"   Time: {elapsed1:.2f}s (per fold: {elapsed1/5:.2f}s)")
    print(f"   Score: {np.mean(scores1):.4f} ± {np.std(scores1):.4f}")

    # Test 2: Refactored style (SVC with class_weight in constructor)
    print("\n[2] Refactored style: SVC(class_weight='balanced')")
    model2 = SVC(kernel='linear', class_weight='balanced', cache_size=1000, random_state=42, C=C)

    start = time.time()
    scores2 = cross_val_score(model2, X_train, y_train, cv=5, scoring='balanced_accuracy', n_jobs=1)
    elapsed2 = time.time() - start
    print(f"   Time: {elapsed2:.2f}s (per fold: {elapsed2/5:.2f}s)")
    print(f"   Score: {np.mean(scores2):.4f} ± {np.std(scores2):.4f}")

    # Test 3: SVC without probability (to isolate probability overhead)
    print("\n[3] SVC without probability: SVC(kernel='linear') + class_weight attribute")
    model3 = SVC(kernel='linear', cache_size=1000, random_state=42, C=C)
    model3.class_weight = 'balanced'

    start = time.time()
    scores3 = cross_val_score(model3, X_train, y_train, cv=5, scoring='balanced_accuracy', n_jobs=1)
    elapsed3 = time.time() - start
    print(f"   Time: {elapsed3:.2f}s (per fold: {elapsed3/5:.2f}s)")
    print(f"   Score: {np.mean(scores3):.4f} ± {np.std(scores3):.4f}")

    # Test 4: LinearSVC (modern approach)
    print("\n[4] LinearSVC: LinearSVC(class_weight='balanced', dual='auto')")
    model4 = LinearSVC(class_weight='balanced', dual='auto', max_iter=2000, random_state=42, C=C)

    start = time.time()
    scores4 = cross_val_score(model4, X_train, y_train, cv=5, scoring='balanced_accuracy', n_jobs=1)
    elapsed4 = time.time() - start
    print(f"   Time: {elapsed4:.2f}s (per fold: {elapsed4/5:.2f}s)")
    print(f"   Score: {np.mean(scores4):.4f} ± {np.std(scores4):.4f}")

    results.append({
        'C': C,
        'orig': elapsed1,
        'refact': elapsed2,
        'no_prob': elapsed3,
        'linear': elapsed4
    })

    print(f"\n   Timing comparison:")
    print(f"   - Refactored vs Original: {elapsed2/elapsed1:.2f}x {'SLOWER' if elapsed2 > elapsed1 else 'FASTER'}")
    print(f"   - Refactored vs LinearSVC: {elapsed2/elapsed4:.2f}x slower")
    print(f"   - No-prob vs Original: {elapsed3/elapsed1:.2f}x {'SLOWER' if elapsed3 > elapsed1 else 'FASTER'}")
    print()

# Summary table
print("\n" + "=" * 80)
print("SUMMARY TABLE")
print("=" * 80)
print(f"{'C Value':>12} | {'Original':>10} | {'Refactored':>10} | {'No Prob':>10} | {'LinearSVC':>10} | {'Refact/Orig':>12}")
print("-" * 80)

for r in results:
    ratio = r['refact'] / r['orig']
    print(f"{r['C']:>12.0e} | {r['orig']:>9.2f}s | {r['refact']:>9.2f}s | {r['no_prob']:>9.2f}s | {r['linear']:>9.2f}s | {ratio:>11.2f}x")

print()
print("=" * 80)
print("KEY INSIGHTS")
print("=" * 80)
print("1. If Original ~= No-Prob: probability=True is NOT causing slowdown")
print("2. If Original < Refactored at high C: Something else is different")
print("3. If LinearSVC much faster: liblinear is the right solution")
print("4. Watch for performance divergence as C increases")
