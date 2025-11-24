#!/usr/bin/env python3
"""
Test that matrix conversion from dict to numpy preserves values.
"""

import sys
from pathlib import Path
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).parent))

from intronIC.scoring.pwm import PWMLoader, BASE_TO_INDEX

# Load matrices using refactored loader
matrix_file = Path(__file__).parent.parent.parent.parent / "data" / "scoring_matrices.fasta.iic"
pwm_sets = PWMLoader.load_from_file(matrix_file)

# Also load using original-style parsing to compare
def load_original_style(filepath):
    """Load matrices the way original intronIC does."""
    matrices = {}

    with open(filepath, 'r') as f:
        current_name = None
        current_matrix = None
        current_start_index = 0
        bases_order = None
        row_num = None

        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            if line.startswith('>'):
                # Save previous
                if current_name and current_matrix:
                    matrices[current_name] = {
                        'matrix': dict(current_matrix),
                        'start_index': current_start_index
                    }

                # Parse header
                header_parts = line[1:].split()
                current_name = header_parts[0]
                start_index = 0
                for part in header_parts:
                    if 'start=' in part:
                        start_index = int(part.split('=')[1])
                        break
                current_start_index = start_index
                current_matrix = defaultdict(dict)
                bases_order = None
                row_num = 0

            elif bases_order is None:
                bases_order = [b for b in line.upper().split() if b in 'ACGT']
            else:
                freqs = [float(f) for f in line.split()]
                position = row_num + current_start_index  # This is the key difference!
                for base, freq in zip(bases_order, freqs):
                    current_matrix[base][position] = freq
                row_num += 1

        # Save last
        if current_name and current_matrix:
            matrices[current_name] = {
                'matrix': dict(current_matrix),
                'start_index': current_start_index
            }

    return matrices

orig_matrices = load_original_style(matrix_file)

# Compare U12 five matrix
print("=" * 70)
print("Matrix Conversion Test")
print("=" * 70)

# Find the U12 ATAC five matrix in original
u12_five_orig = None
for name, data in orig_matrices.items():
    if 'u12' in name.lower() and 'atac' in name.lower() and 'five' in name.lower():
        u12_five_orig = data
        print(f"Found original matrix: {name}")
        break

if not u12_five_orig:
    print("ERROR: Could not find U12 five matrix in original")
    sys.exit(1)

# Get refactored PWM
u12_five_pwm = pwm_sets['five'].u12_canonical

print(f"\nOriginal start_index: {u12_five_orig['start_index']}")
print(f"Refactored start_index: {u12_five_pwm.start_index}")
print(f"Refactored length: {u12_five_pwm.length}")

# Compare a few positions
test_positions = [-3, -2, -1, 0, 1, 2, 3, 4, 5]
print("\nComparing frequencies at key positions:")
print("Pos  Base  Original    Refactored  Match?")
print("-" * 70)

all_match = True
for pos in test_positions:
    for base in ['A', 'C', 'G', 'T']:
        # Original lookup
        orig_freq = u12_five_orig['matrix'].get(base, {}).get(pos, None)

        # Refactored lookup
        matrix_idx = pos - u12_five_pwm.start_index
        if 0 <= matrix_idx < u12_five_pwm.length:
            base_idx = BASE_TO_INDEX[base]
            ref_freq = u12_five_pwm.matrix[base_idx, matrix_idx]
        else:
            ref_freq = None

        if orig_freq is not None and ref_freq is not None:
            match = "✓" if abs(orig_freq - ref_freq) < 1e-10 else "✗"
        else:
            match = "✓" if orig_freq == ref_freq else "✗"

        if match == "✗":
            all_match = False
            print(f"{pos:3d}  {base}     {orig_freq}  {ref_freq}  {match}")

if all_match:
    print("\n✓ All positions match!")
else:
    print("\n✗ Some positions don't match!")

# Now test scoring with both
print("\n" + "=" * 70)
print("Scoring Comparison")
print("=" * 70)

test_seq = "TCAGTATCCTTC"
start_pos = -3

print(f"Sequence: {test_seq}")
print(f"Start position: {start_pos}")
print()

# Score with refactored
ref_score = u12_five_pwm.score_sequence(test_seq, seq_start_position=start_pos)
print(f"Refactored score: {ref_score}")

# Score with original-style
orig_score = None
for i, base in enumerate(test_seq):
    pos = start_pos + i
    freq = u12_five_orig['matrix'].get(base, {}).get(pos, 0.0001)
    if freq == 0.0:
        freq = 0.0001
    if orig_score is None:
        orig_score = freq
    else:
        orig_score *= freq

print(f"Original-style score: {orig_score}")
print()

if abs(orig_score - ref_score) < 1e-10:
    print("✓ Scores match!")
else:
    print(f"✗ Scores don't match! Difference: {abs(orig_score - ref_score)}")

print("=" * 70)
