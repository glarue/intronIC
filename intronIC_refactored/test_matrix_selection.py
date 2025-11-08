#!/usr/bin/env python3
"""
Test that matrix selection works correctly.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from scoring.pwm import PWMLoader

# Load matrices
matrix_file = Path(__file__).parent.parent / "intronIC" / "data" / "scoring_matrices.fasta.iic"
pwm_sets = PWMLoader.load_from_file(matrix_file)

print("=" * 70)
print("Matrix Selection Test")
print("=" * 70)
print()

# Check what matrices are available
five_set = pwm_sets['five']
print("Available 5' matrices:")
for key, pwm in five_set.matrices.items():
    print(f"  {key}: {pwm.name}")
print()

# Test selection for different dinucleotides
test_cases = [
    ('u12', 'gtag'),  # Should select u12_gtag
    ('u12', 'atac'),  # Should select u12_atac
    ('u2', 'gtag'),   # Should select u2_gtag
    ('u2', 'gcag'),   # Should select u2_gcag
    ('u12', 'gcag'),  # Should fall back to u12_gtag (most common)
]

print("Matrix selection tests:")
for intron_type, dnts in test_cases:
    selected = five_set.select_best(intron_type, dnts)
    print(f"  {intron_type}, {dnts:5s} → {selected.name}")
print()

# Test scoring with correct matrices
print("Scoring test with GT-AG intron:")
test_seq = "TCAGTATCCTTC"
seq_start_pos = -3

u12_pwm = five_set.select_best('u12', 'gtag')
u2_pwm = five_set.select_best('u2', 'gtag')

u12_score = u12_pwm.score_sequence(test_seq, seq_start_position=seq_start_pos)
u2_score = u2_pwm.score_sequence(test_seq, seq_start_position=seq_start_pos)

import math
log_ratio = math.log2(u12_score / u2_score)

print(f"  U12 GTAG score: {u12_score:.6e}")
print(f"  U2 GTAG score: {u2_score:.6e}")
print(f"  Log ratio: {log_ratio:.4f}")
print(f"  Expected: 18.2163")

if abs(log_ratio - 18.2163) < 0.01:
    print("  ✓ Scoring works correctly!")
else:
    print(f"  ✗ Scoring mismatch! Diff: {abs(log_ratio - 18.2163)}")

print("=" * 70)
