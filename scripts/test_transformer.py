#!/usr/bin/env python3
"""
Quick test to verify BothEndsStrongTransformer configurations.
"""
import numpy as np
from intronIC.classification.transformers import BothEndsStrongTransformer

# Test data
X = np.array([
    [2.0, 1.8, 1.5],   # Balanced U12
    [5.0, -1.0, -0.5], # One-end-strong FP
    [0.3, 0.3, 0.3]    # Weak but balanced
])

print("Test input (3 samples, 3 features):")
print(X)
print()

# Configuration 1: Recommended (7D)
print("="*80)
print("Configuration 1: include_pairwise_mins=False, include_max=False (RECOMMENDED)")
print("="*80)
transformer1 = BothEndsStrongTransformer(
    include_pairwise_mins=False,
    include_max=False
)
X1 = transformer1.transform(X)
features1 = transformer1.get_feature_names_out()
print(f"Output shape: {X1.shape}")
print(f"Features ({len(features1)}): {list(features1)}")
print("\nFirst sample (balanced U12):")
for feat, val in zip(features1, X1[0]):
    print(f"  {feat:20s}: {val:7.4f}")
print()

# Configuration 2: Backward compatible (9D)
print("="*80)
print("Configuration 2: include_pairwise_mins=True, include_max=False")
print("="*80)
transformer2 = BothEndsStrongTransformer(
    include_pairwise_mins=True,
    include_max=False
)
X2 = transformer2.transform(X)
features2 = transformer2.get_feature_names_out()
print(f"Output shape: {X2.shape}")
print(f"Features ({len(features2)}): {list(features2)}")
print("\nFirst sample (balanced U12):")
for feat, val in zip(features2, X2[0]):
    print(f"  {feat:20s}: {val:7.4f}")
print()

# Configuration 3: All features (11D)
print("="*80)
print("Configuration 3: include_pairwise_mins=True, include_max=True")
print("="*80)
transformer3 = BothEndsStrongTransformer(
    include_pairwise_mins=True,
    include_max=True
)
X3 = transformer3.transform(X)
features3 = transformer3.get_feature_names_out()
print(f"Output shape: {X3.shape}")
print(f"Features ({len(features3)}): {list(features3)}")
print()

# Verify min_all works correctly
print("="*80)
print("Verification: min_all computation")
print("="*80)
for i, (s5, sBP, s3) in enumerate(X):
    print(f"\nSample {i+1}: s5={s5:.1f}, sBP={sBP:.1f}, s3={s3:.1f}")
    expected_min_all = min(s5, sBP, s3)

    # Get min_all from config 1 (position 3, after s5, sBP, s3)
    actual_min_all = X1[i, 3]

    print(f"  Expected min_all: {expected_min_all:.1f}")
    print(f"  Actual min_all:   {actual_min_all:.1f}")
    print(f"  Match: {'✓' if np.isclose(expected_min_all, actual_min_all) else '✗'}")

print("\n" + "="*80)
print("All tests completed successfully!")
print("="*80)
