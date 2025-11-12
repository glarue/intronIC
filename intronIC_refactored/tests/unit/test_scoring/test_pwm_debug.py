#!/usr/bin/env python3
"""
Debug PWM scoring to see exactly what's happening.
"""

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from scoring.pwm import PWMLoader, BASE_TO_INDEX

# Load the PWM matrices
matrix_file = Path(__file__).parent.parent / "intronIC" / "data" / "scoring_matrices.fasta.iic"
pwm_sets = PWMLoader.load_from_file(matrix_file)

# Get the U12 5' PWM
u12_five_pwm = pwm_sets['five'].u12_canonical

# Test sequence from ENST00000359435_5
test_seq = "TCAGTATCCTTC"
seq_start_position = -3  # 5' region starts at position -3

print("=" * 70)
print("PWM Scoring Debug")
print("=" * 70)
print(f"Sequence: {test_seq}")
print(f"seq_start_position: {seq_start_position}")
print(f"PWM start_index: {u12_five_pwm.start_index}")
print(f"PWM length: {u12_five_pwm.length}")
print()
print("Position mapping:")
print("i  base  logical_pos  matrix_index  freq")
print("-" * 70)

score = None
for i, base in enumerate(test_seq):
    logical_position = seq_start_position + i
    matrix_index = logical_position - u12_five_pwm.start_index

    # Get frequency
    if matrix_index < 0 or matrix_index >= u12_five_pwm.length:
        freq = "SKIP (out of bounds)"
        continue
    else:
        base_index = BASE_TO_INDEX[base]
        freq = u12_five_pwm.matrix[base_index, matrix_index]
        if freq == 0.0:
            freq = u12_five_pwm.pseudocount

        if score is None:
            score = freq
        else:
            score *= freq

    print(f"{i:2d}  {base:4s}  {logical_position:11d}  {matrix_index:12d}  {freq}")

print()
print(f"Final score: {score}")
print(f"Final score (scientific): {score:.15e}")
print()
print("Now let's check what the matrix values are at each position:")
print()

# Check a few key positions
for i in range(len(test_seq)):
    logical_position = seq_start_position + i
    matrix_index = logical_position - u12_five_pwm.start_index
    base = test_seq[i]

    if 0 <= matrix_index < u12_five_pwm.length:
        print(f"Position {logical_position}: {base} → ", end="")
        for b in ['A', 'C', 'G', 'T']:
            freq = u12_five_pwm.matrix[BASE_TO_INDEX[b], matrix_index]
            print(f"{b}:{freq:.4f} ", end="")
        print()

print("=" * 70)
