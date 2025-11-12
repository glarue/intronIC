#!/usr/bin/env python3
"""
Test that log ratio calculation matches original.
"""

import sys
from pathlib import Path
import math

sys.path.insert(0, str(Path(__file__).parent))

from scoring.pwm import PWMLoader

# Load matrices
matrix_file = Path(__file__).parent.parent / "intronIC" / "data" / "scoring_matrices.fasta.iic"
pwm_sets = PWMLoader.load_from_file(matrix_file)

# Get PWMs
u12_five_pwm = pwm_sets['five'].u12_canonical
u2_five_pwm = pwm_sets['five'].u2_canonical

# Test sequence and position
test_seq = "TCAGTATCCTTC"
seq_start_position = -3

print("=" * 70)
print("Log Ratio Calculation Test")
print("=" * 70)
print(f"Sequence: {test_seq}")
print(f"Start position: {seq_start_position}")
print()

# Score with U12
u12_score = u12_five_pwm.score_sequence(test_seq, seq_start_position=seq_start_position)
print(f"U12 score: {u12_score}")
print(f"U12 score (sci): {u12_score:.15e}")

# Score with U2
u2_score = u2_five_pwm.score_sequence(test_seq, seq_start_position=seq_start_position)
print(f"U2 score: {u2_score}")
print(f"U2 score (sci): {u2_score:.15e}")

# Calculate log ratio
log_ratio = math.log2(u12_score / u2_score)
print()
print(f"Log ratio (U12/U2): {log_ratio}")
print(f"Expected: 18.216262900453128")
print()

diff = abs(log_ratio - 18.216262900453128)
if diff < 0.01:
    print(f"✓ TEST PASSED! Difference: {diff:.10f}")
else:
    print(f"✗ TEST FAILED! Difference: {diff}")

print("=" * 70)
