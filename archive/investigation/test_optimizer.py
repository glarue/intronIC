#!/usr/bin/env python
"""Test optimizer parallelization directly."""

import sys
sys.path.insert(0, '/home/user/intronIC')

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

from intronIC_refactored.classification.optimizer import SVMOptimizer
from intronIC_refactored.core.intron import Intron, IntronScores
import numpy as np

# Create fake introns with z-scores
def make_fake_intron(idx, is_u12=False):
    intron = Intron(
        region=f"chr1",
        start=idx*100,
        stop=idx*100+50,
        strand="+",
        parent=f"transcript_{idx}"
    )
    # Fake z-scores: U12s have positive values, U2s have negative
    if is_u12:
        z_scores = (2.0 + np.random.rand(), 3.0 + np.random.rand(), 2.5 + np.random.rand())
    else:
        z_scores = (-2.0 + np.random.rand(), -3.0 + np.random.rand(), -2.5 + np.random.rand())

    intron.scores = IntronScores(
        five_z_score=z_scores[0],
        bp_z_score=z_scores[1],
        three_z_score=z_scores[2]
    )
    return intron

print("Creating fake training data...")
u12s = [make_fake_intron(i, is_u12=True) for i in range(50)]
u2s = [make_fake_intron(i, is_u12=False) for i in range(200)]

print(f"Created {len(u12s)} U12 and {len(u2s)} U2 introns")

print("\nTesting with n_jobs=4...")
optimizer = SVMOptimizer(n_rounds=2, n_jobs=4)
result = optimizer.optimize(u12s, u2s)
print(f"Result: C={result.C:.6e}, score={result.cv_score:.4f}")

print("\nDone!")
