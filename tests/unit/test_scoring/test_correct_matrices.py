#!/usr/bin/env python3
"""
Test with correct matrix selection (GT-AG intron should use gtag matrices).
"""

import sys
from pathlib import Path
import math
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).parent))

from intronIC.scoring.pwm import PWMLoader

# Load matrices
matrix_file = Path(__file__).parent.parent.parent.parent / "data" / "scoring_matrices.fasta.iic"

# Parse matrices directly to access all of them
matrices_raw = PWMLoader._parse_matrix_file(matrix_file)

# Find U12 GTAG and U2 GTAG five matrices
u12_gtag_five = None
u2_gtag_five = None

for name_tuple, data in matrices_raw.items():
    name_str = '_'.join(name_tuple)
    if 'u12' in name_tuple and 'gtag' in name_tuple and 'five' in name_tuple:
        print(f"Found U12 GTAG five: {name_tuple}")
        u12_gtag_five = PWMLoader._dict_to_pwm(
            name=name_str,
            matrix_dict=data['matrix'],
            start_index=data['start_index']
        )
    elif 'u2' in name_tuple and 'gtag' in name_tuple and 'five' in name_tuple:
        print(f"Found U2 GTAG five: {name_tuple}")
        u2_gtag_five = PWMLoader._dict_to_pwm(
            name=name_str,
            matrix_dict=data['matrix'],
            start_index=data['start_index']
        )

if not u12_gtag_five or not u2_gtag_five:
    print("ERROR: Could not find both matrices")
    sys.exit(1)

# Test sequence and position
test_seq = "TCAGTATCCTTC"
seq_start_position = -3

print("\n" + "=" * 70)
print("Correct Matrix Selection Test")
print("=" * 70)
print(f"Sequence: {test_seq}")
print(f"Start position: {seq_start_position}")
print()

# Score with U12 GTAG
u12_score = u12_gtag_five.score_sequence(test_seq, seq_start_position=seq_start_position)
print(f"U12 GTAG score: {u12_score}")
print(f"U12 GTAG score (sci): {u12_score:.15e}")

# Score with U2 GTAG
u2_score = u2_gtag_five.score_sequence(test_seq, seq_start_position=seq_start_position)
print(f"U2 GTAG score: {u2_score}")
print(f"U2 GTAG score (sci): {u2_score:.15e}")

# Calculate log ratio
log_ratio = math.log2(u12_score / u2_score)
print()
print(f"Log ratio (U12_GTAG / U2_GTAG): {log_ratio}")
print(f"Expected: 18.216262900453128")
print()

diff = abs(log_ratio - 18.216262900453128)
if diff < 0.01:
    print(f"✓ TEST PASSED! Difference: {diff:.10f}")
else:
    print(f"✗ TEST FAILED! Difference: {diff}")

print("=" * 70)
