#!/usr/bin/env python3
"""
Test that PWM scoring fix produces correct scores.

Tests the 5' sequence from ENST00000359435_5 which should score +18.22
according to the original intronIC output.
"""

import sys
from pathlib import Path

# Add parent directory to path to import modules
sys.path.insert(0, str(Path(__file__).parent))

from intronIC.scoring.pwm import PWMLoader
import math

# Load the PWM matrices
matrix_file = Path(__file__).parent.parent.parent.parent / "data" / "scoring_matrices.fasta.iic"
pwm_sets = PWMLoader.load_from_file(matrix_file)

# Get the U12 5' PWM
u12_five_pwm = pwm_sets['five'].u12_canonical

# Test sequence from ENST00000359435_5
test_seq = "TCAGTATCCTTC"
seq_start_position = -3  # 5' region starts at position -3

# Score the sequence
score = u12_five_pwm.score_sequence(test_seq, seq_start_position=seq_start_position)

# Calculate log ratio (original output shows raw score, not log ratio)
# Actually, the raw score IS the log ratio in the original code
# No wait - looking at the original output:
# "5prime raw: 18.216262900453128 z: 5.100993754895024"
# This is the log ratio, not the PWM score itself

# Let me also test the U2 score to calculate the log ratio
u2_five_pwm = pwm_sets['five'].u2_canonical
u2_score = u2_five_pwm.score_sequence(test_seq, seq_start_position=seq_start_position)

# Calculate log ratio
if u2_score > 0 and score > 0:
    log_ratio = math.log2(score / u2_score)
else:
    log_ratio = None

print("=" * 60)
print("PWM Scoring Test Results")
print("=" * 60)
print(f"Test sequence: {test_seq}")
print(f"Position: {seq_start_position}")
print(f"PWM start_index: {u12_five_pwm.start_index}")
print(f"PWM length: {u12_five_pwm.length}")
print()
print("U12 5' PWM score:", score)
print("U2 5' PWM score:", u2_score)
if log_ratio is not None:
    print("Log ratio (U12/U2):", log_ratio)
print()
print("Expected log ratio: 18.216262900453128")
print()
if log_ratio is not None:
    if abs(log_ratio - 18.216262900453128) < 0.01:
        print("✓ TEST PASSED - Score matches expected value!")
    else:
        print("✗ TEST FAILED - Score does not match!")
        print(f"  Difference: {log_ratio - 18.216262900453128}")
else:
    print("✗ TEST FAILED - Could not calculate log ratio")
print("=" * 60)
