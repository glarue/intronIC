"""
Position Weight Matrix (PWM) scoring for intron sequences.

This module implements PWM-based scoring for splice sites using log-odds ratios.
The matrices are trained on known U2 and U12 intron sequences.

Port from: intronIC.py:2114-2142 (seq_score), 1180-1264 (load_external_matrix)

Design:
- Immutable PWM objects (frozen dataclasses)
- Efficient numpy-based scoring
- Support for ignore_positions (to skip canonical dinucleotides)
- Pseudocount handling for ambiguous bases
"""

from dataclasses import dataclass
from typing import Optional, Set, Dict
from pathlib import Path
from collections import defaultdict
import numpy as np


# Base ordering for matrix indexing: ACGT
BASE_ORDER = ['A', 'C', 'G', 'T']
BASE_TO_INDEX = {base: idx for idx, base in enumerate(BASE_ORDER)}


@dataclass(frozen=True, slots=True)
class PWM:
    """
    Position Weight Matrix for sequence scoring.

    The matrix stores base frequencies at each position in a motif.
    Scoring multiplies the frequencies for each base in the sequence
    to get a log-odds ratio.

    Attributes:
        name: Identifier for this PWM (e.g., "u12_atac_five")
        matrix: Numpy array of shape (4, length) for ACGT frequencies
        length: int
        pseudocount: Value to use for zero frequencies (default: 0.0001)
        start_index: Optional offset for position numbering (default: 0)
    """
    name: str
    matrix: np.ndarray  # Shape: (4, length) for ACGT
    length: int
    pseudocount: float = 0.0001
    start_index: int = 0

    def __post_init__(self):
        """Validate PWM structure."""
        # Verify matrix has 4 rows (one for each base: ACGT)
        if self.matrix.shape[0] != 4:
            raise ValueError(
                f"PWM matrix must have 4 rows (ACGT), got {self.matrix.shape[0]}"
            )

        # Verify matrix length matches declared length
        if self.matrix.shape[1] != self.length:
            raise ValueError(
                f"PWM matrix has {self.matrix.shape[1]} positions, "
                f"but length={self.length}"
            )

    def score_sequence(
        self,
        seq: str,
        ignore_positions: Optional[Set[int]] = None
    ) -> float:
        """
        Score sequence using PWM (product of base frequencies).

        This implements the log-odds scoring from the original intronIC,
        where the score is the product of base frequencies at each position.

        Port from: intronIC.py:2114-2140 (seq_score)

        Args:
            seq: Sequence to score (must match PWM length)
            ignore_positions: Set of PWM positions to ignore (set freq to 1.0)
                             Used to skip canonical dinucleotides in scoring

        Returns:
            Score as product of frequencies (float)

        Raises:
            ValueError: If sequence length doesn't match PWM length
            ValueError: If sequence contains lowercase (should be uppercase)

        Example:
            >>> pwm = PWM("test", matrix, length=4)
            >>> pwm.score_sequence("ACGT")
            0.125
            >>> # Ignore positions 0 and 1 (canonical GT)
            >>> pwm.score_sequence("GTAAGT", ignore_positions={0, 1})
            0.5
        """
        # Validate sequence is not empty
        if len(seq) == 0:
            raise ValueError(
                "Cannot score empty sequence. "
                "This indicates a bug in sequence extraction or filtering."
            )

        # Check for lowercase (warn user to use uppercase)
        if seq != seq.upper():
            raise ValueError(
                f"Sequence must be uppercase. Got: {seq}. "
                "Please convert to uppercase before scoring."
            )

        # Initialize score
        # Port from: intronIC.py:2125
        score = None

        # Iterate through sequence positions
        # Port from: intronIC.py:2126-2138
        for i, base in enumerate(seq):
            # Calculate position in PWM (accounting for start_index)
            pwm_position = i + self.start_index

            # Skip positions outside PWM length (for flexible scoring)
            if pwm_position >= self.length:
                continue

            # Check if this position should be ignored
            # Port from: intronIC.py:2127-2128
            if ignore_positions is not None and pwm_position in ignore_positions:
                base_freq = 1.0
            else:
                # Look up base frequency in matrix
                try:
                    base_index = BASE_TO_INDEX[base]
                    # Use pwm_position to index into matrix, not i
                    base_freq = self.matrix[base_index, pwm_position]

                    # Apply pseudocount if frequency is zero
                    # Port from: intronIC.py:2134
                    if base_freq == 0.0:
                        base_freq = self.pseudocount

                except KeyError:
                    # Ambiguous base (N, etc.) - use pseudocount
                    # Port from: intronIC.py:2132-2134
                    base_freq = self.pseudocount

            # Multiply into score
            # Port from: intronIC.py:2135-2138
            if score is None:
                score = base_freq
            else:
                score *= base_freq

        return score


@dataclass(frozen=True, slots=True)
class PWMSet:
    """
    U2 and U12 PWMs for a single splice site region.

    Each region (five, bp, three) has separate PWMs for:
    - U2 canonical (e.g., GT-AG)
    - U2 non-canonical (e.g., GC-AG)
    - U12 canonical (e.g., AT-AC)
    - U12 non-canonical (rare)

    Attributes:
        u2_canonical: U2-type PWM for canonical boundaries
        u2_noncanonical: Optional U2-type PWM for non-canonical boundaries
        u12_canonical: U12-type PWM for canonical boundaries
        u12_noncanonical: Optional U12-type PWM for non-canonical boundaries
    """
    u2_canonical: PWM
    u2_noncanonical: Optional[PWM]
    u12_canonical: PWM
    u12_noncanonical: Optional[PWM]


class PWMLoader:
    """
    Load PWMs from scoring_matrices.fasta.iic format.

    The file format is:
        >matrix_name  start={n}  (n={count})
        A    C    G    T
        0.25 0.25 0.25 0.25
        0.30 0.20 0.30 0.20
        ...

    Matrix names follow pattern: {type}_{boundaries}_{region}
    Examples:
        - u12_atac_five
        - u12_atac_bp
        - u2_gtag_five
        - u2_gcag_three

    Port from: intronIC.py:1180-1264 (load_external_matrix)
    """

    @staticmethod
    def load_from_file(filepath: Path) -> Dict[str, PWMSet]:
        """
        Load all PWMs from scoring_matrices.fasta.iic file.

        Args:
            filepath: Path to scoring_matrices.fasta.iic

        Returns:
            Dictionary mapping region to PWMSet:
            {
                'five': PWMSet(...),
                'bp': PWMSet(...),
                'three': PWMSet(...)
            }

        Raises:
            FileNotFoundError: If file doesn't exist
            ValueError: If file format is invalid
        """
        if not filepath.exists():
            raise FileNotFoundError(f"Matrix file not found: {filepath}")

        # Parse all matrices from file
        # Port from: intronIC.py:1242-1264
        matrices = PWMLoader._parse_matrix_file(filepath)

        # Group matrices by region (five, bp, three)
        pwm_sets = PWMLoader._group_into_pwm_sets(matrices)

        return pwm_sets

    @staticmethod
    def _parse_matrix_file(filepath: Path) -> Dict[tuple, Dict[str, Dict[int, float]]]:
        """
        Parse matrix file into raw frequency dictionaries.

        Port from: intronIC.py:1242-1264

        Returns:
            Dictionary mapping (type, boundary, region) tuples to matrix data:
            {
                ('u12', 'atac', 'five'): {
                    'A': {0: 0.25, 1: 0.30, ...},
                    'C': {0: 0.25, 1: 0.20, ...},
                    ...
                },
                start_index: -3  # Stored in special key
            }
        """
        matrices = {}

        # Read file and parse FASTA-like format
        current_name = None
        current_matrix = None
        current_start_index = 0
        bases_order = None

        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()

                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    continue

                # Header line: >matrix_name  start={n}
                if line.startswith('>'):
                    # Save previous matrix if exists
                    if current_name is not None and current_matrix is not None:
                        parsed_name = PWMLoader._parse_matrix_name(current_name)
                        matrices[parsed_name] = {
                            'matrix': current_matrix,
                            'start_index': current_start_index
                        }

                    # Parse new header
                    # Port from: intronIC.py:1249-1254
                    header_parts = line[1:].split()  # Remove '>'
                    current_name = header_parts[0]

                    # Extract start_index if present
                    start_index = 0
                    for part in header_parts:
                        if 'start=' in part:
                            start_index = int(part.split('=')[1])
                            break
                    current_start_index = start_index

                    # Initialize new matrix
                    current_matrix = defaultdict(dict)
                    bases_order = None

                # Base order line: A  C  G  T
                elif bases_order is None:
                    # Port from: intronIC.py:1258
                    bases_order = [b for b in line.upper().split() if b in 'ACGT']

                # Frequency line: 0.25  0.25  0.25  0.25
                else:
                    # Port from: intronIC.py:1259-1262
                    freqs = [float(f) for f in line.split()]
                    position = len(current_matrix['A']) + current_start_index

                    for base, freq in zip(bases_order, freqs):
                        current_matrix[base][position] = freq

        # Save last matrix
        if current_name is not None and current_matrix is not None:
            parsed_name = PWMLoader._parse_matrix_name(current_name)
            matrices[parsed_name] = {
                'matrix': current_matrix,
                'start_index': current_start_index
            }

        return matrices

    @staticmethod
    def _parse_matrix_name(matrix_name: str) -> tuple:
        """
        Parse matrix name into (type, boundary, region) tuple.

        Port from: intronIC.py:1203-1240 (__name_parser)

        Examples:
            "u12_atac_five" -> ('u12', 'atac', 'five')
            "u2_gtag_bp" -> ('u2', 'gtag', 'bp')
            "u12_atac_three_v2" -> ('u12', 'atac', 'three', 'v2')

        Returns:
            Tuple of parsed name components
        """
        # Port from: intronIC.py:1212-1240
        subtypes = {
            "u12": ["u12", "12", "minor"],
            "u2": ["u2", "major"]
        }
        regions = {
            "five": ["five", "5"],
            "bp": ["bp", "branch-point", "branchpoint"],
            "three": ["three", "3"]
        }
        boundaries = {
            "atac": ["at-ac", "atac"],
            "gtag": ["gt-ag", "gtag"],
            "gcag": ["gc-ag", "gcag"]
        }

        name_bits = []
        name_lower = matrix_name.lower()

        # Extract each component
        for category in [subtypes, boundaries, regions]:
            for key, patterns in category.items():
                if any(pattern in name_lower for pattern in patterns):
                    name_bits.append(key)
                    break

        return tuple(name_bits)

    @staticmethod
    def _group_into_pwm_sets(
        matrices: Dict[tuple, Dict]
    ) -> Dict[str, PWMSet]:
        """
        Group parsed matrices into PWMSet objects by region.

        Args:
            matrices: Dictionary from _parse_matrix_file

        Returns:
            Dictionary mapping region name to PWMSet:
            {'five': PWMSet, 'bp': PWMSet, 'three': PWMSet}
        """
        # Group by region
        by_region = defaultdict(dict)

        for name_tuple, data in matrices.items():
            # Extract components
            if len(name_tuple) < 3:
                continue  # Skip malformed names

            intron_type = name_tuple[0]  # 'u2' or 'u12'
            boundary = name_tuple[1] if len(name_tuple) > 1 else 'gtag'
            region = name_tuple[2] if len(name_tuple) > 2 else None

            if region not in ['five', 'bp', 'three']:
                continue  # Skip unknown regions

            # Convert to PWM object
            pwm = PWMLoader._dict_to_pwm(
                name='_'.join(name_tuple),
                matrix_dict=data['matrix'],
                start_index=data['start_index']
            )

            # Store in appropriate category
            # Canonical: GTAG (U2), ATAC (U12)
            # Non-canonical: GCAG, etc.
            is_canonical = (boundary == 'gtag' and intron_type == 'u2') or \
                          (boundary == 'atac' and intron_type == 'u12')

            key = f"{intron_type}_{'canonical' if is_canonical else 'noncanonical'}"
            by_region[region][key] = pwm

        # Build PWMSets
        pwm_sets = {}
        for region in ['five', 'bp', 'three']:
            if region not in by_region:
                continue

            pwms = by_region[region]

            pwm_sets[region] = PWMSet(
                u2_canonical=pwms.get('u2_canonical'),
                u2_noncanonical=pwms.get('u2_noncanonical'),
                u12_canonical=pwms.get('u12_canonical'),
                u12_noncanonical=pwms.get('u12_noncanonical')
            )

        return pwm_sets

    @staticmethod
    def _dict_to_pwm(
        name: str,
        matrix_dict: Dict[str, Dict[int, float]],
        start_index: int
    ) -> PWM:
        """
        Convert parsed dictionary to PWM object with numpy array.

        Args:
            name: PWM name
            matrix_dict: {'A': {0: 0.25, 1: 0.30}, 'C': {...}, ...}
            start_index: Starting position index

        Returns:
            PWM object with numpy matrix
        """
        # Determine matrix length
        # Port from: intronIC.py uses position keys to determine length
        all_positions = set()
        for base_freqs in matrix_dict.values():
            all_positions.update(base_freqs.keys())

        if not all_positions:
            raise ValueError(f"Empty matrix for {name}")

        min_pos = min(all_positions)
        max_pos = max(all_positions)
        length = max_pos - min_pos + 1

        # Build numpy array (4 x length)
        matrix = np.zeros((4, length))

        for base_idx, base in enumerate(BASE_ORDER):
            if base in matrix_dict:
                for pos, freq in matrix_dict[base].items():
                    # Convert absolute position to array index
                    array_idx = pos - min_pos
                    matrix[base_idx, array_idx] = freq

        return PWM(
            name=name,
            matrix=matrix,
            length=length,
            start_index=start_index
        )
