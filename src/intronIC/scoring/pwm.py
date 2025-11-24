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
from typing import Optional, Set, Dict, Tuple
from pathlib import Path
from collections import defaultdict
import numpy as np
import re


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
        seq_start_position: int = 0,
        ignore_positions: Optional[Set[int]] = None
    ) -> float:
        """
        Score sequence using PWM (product of base frequencies).

        This implements the log-odds scoring from the original intronIC,
        where the score is the product of base frequencies at each position.

        Port from: intronIC.py:2114-2140 (seq_score)

        Args:
            seq: Sequence to score
            seq_start_position: Logical position of first base in sequence
                               (e.g., -3 for 5' region starting at -3)
                               Corresponds to enumerate(seq, start=start_index)
            ignore_positions: Set of logical positions to ignore (set freq to 1.0)
                             Used to skip canonical dinucleotides in scoring

        Returns:
            Score as product of frequencies (float)

        Raises:
            ValueError: If sequence is empty
            ValueError: If sequence contains lowercase (should be uppercase)

        Example:
            >>> # Score 5' region starting at position -3
            >>> pwm = PWM("u12_five", matrix, length=40, start_index=-20)
            >>> pwm.score_sequence("TCAGTATCCTTC", seq_start_position=-3)
            0.00384  # U12 score
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
        # The original uses: for i, e in enumerate(seq, start=start_index)
        # where start_index is the logical position of the first base
        for i, base in enumerate(seq):
            # Calculate logical position of this base
            # For example: seq_start_position=-3, i=0 → logical_pos=-3
            logical_position = seq_start_position + i

            # Convert logical position to matrix array index
            # For example: logical_pos=-3, start_index=-20 → matrix_index=17
            matrix_index = logical_position - self.start_index

            # Skip positions outside PWM length (for flexible scoring)
            if matrix_index < 0 or matrix_index >= self.length:
                continue

            # Check if this position should be ignored
            # Port from: intronIC.py:2127-2128
            # Note: ignore_positions uses logical positions, not matrix indices
            if ignore_positions is not None and logical_position in ignore_positions:
                base_freq = 1.0
            else:
                # Look up base frequency in matrix
                try:
                    base_index = BASE_TO_INDEX[base]
                    base_freq = self.matrix[base_index, matrix_index]

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

        # If all positions were skipped (score still None), use pseudocount
        # This can happen when sequence position doesn't overlap with PWM coverage
        if score is None:
            score = self.pseudocount

        return score


@dataclass(frozen=True, slots=True)
class PWMSet:
    """
    U2 and U12 PWMs for a single splice site region, organized by dinucleotide type.

    Each region (five, bp, three) has separate PWMs for different intron subtypes:
    - U2 GT-AG (most common U2)
    - U2 GC-AG (non-canonical U2)
    - U12 GT-AG (most common U12 - ~70% of U12 introns)
    - U12 AT-AC (classic U12 marker - ~25% of U12 introns)

    Note: U12 introns are identified by motif characteristics, not just dinucleotides.
    Most U12 introns have GT-AG boundaries but U12-specific motifs throughout.

    Attributes:
        matrices: Dictionary mapping (subtype, boundary[, version]) to PWM
                 e.g., {('u12', 'gtag'): PWM, ('u12', 'gtag', 'A10'): PWM, ...}
                 Keys may have 2 or 3 elements depending on whether version exists
    """
    matrices: Dict[Tuple[str, ...], PWM]

    def select_best(self, intron_type: str, dinucleotides: str) -> PWM:
        """
        Select the best PWM for an intron based on type and dinucleotides.

        If multiple versions exist (e.g., vA9, vA10), returns the first match.
        For scoring with all versions, use get_all_versions() instead.

        Args:
            intron_type: 'u12' or 'u2'
            dinucleotides: e.g., 'gtag', 'atac', 'gcag'

        Returns:
            Best matching PWM

        Example:
            >>> pwm_set.select_best('u12', 'gtag')  # Returns U12 GT-AG matrix
            >>> pwm_set.select_best('u12', 'atac')  # Returns U12 AT-AC matrix
            >>> pwm_set.select_best('u2', 'gtag')   # Returns U2 GT-AG matrix
        """
        # First try exact match without version
        key = (intron_type, dinucleotides)
        if key in self.matrices:
            return self.matrices[key]

        # Try finding any version of this type+dinucleotide combination
        for matrix_key, pwm in self.matrices.items():
            if len(matrix_key) >= 2 and matrix_key[0] == intron_type and matrix_key[1] == dinucleotides:
                return pwm

        # Fall back to most common type for this subtype
        if intron_type == 'u12':
            # Try GT-AG first (most common U12), then AT-AC
            for fallback in ['gtag', 'atac']:
                # Try without version first
                key = ('u12', fallback)
                if key in self.matrices:
                    return self.matrices[key]
                # Try with any version
                for matrix_key, pwm in self.matrices.items():
                    if len(matrix_key) >= 2 and matrix_key[0] == 'u12' and matrix_key[1] == fallback:
                        return pwm
        else:  # u2
            # Try GT-AG first (most common U2), then GC-AG
            for fallback in ['gtag', 'gcag']:
                # Try without version first
                key = ('u2', fallback)
                if key in self.matrices:
                    return self.matrices[key]
                # Try with any version
                for matrix_key, pwm in self.matrices.items():
                    if len(matrix_key) >= 2 and matrix_key[0] == 'u2' and matrix_key[1] == fallback:
                        return pwm

        # Last resort: return any matrix of the right type
        for key, pwm in self.matrices.items():
            if key[0] == intron_type:
                return pwm

        raise ValueError(f"No PWM found for {intron_type}")

    def get_all_versions(self, intron_type: str, dinucleotides: str) -> list[PWM]:
        """
        Get all PWM versions for an intron type and dinucleotides.

        This returns all matrices (including different versions like vA9, vA10)
        that match the given type and dinucleotides. The scorer can then try
        each version and select the highest-scoring one.

        Port from: intronIC.py:2915-2942 (multi_matrix_score loops through all matrices)

        Args:
            intron_type: 'u12' or 'u2'
            dinucleotides: e.g., 'gtag', 'atac', 'gcag'

        Returns:
            List of all matching PWMs (may include multiple versions)

        Example:
            >>> pwm_set.get_all_versions('u12', 'gtag')
            [PWM(name='u12_gtag_bp_vA10'), PWM(name='u12_gtag_bp_vA9')]
        """
        matching_pwms = []

        # Find all matrices matching (intron_type, dinucleotides)
        # Keys may be (type, dnt) or (type, dnt, version)
        for key, pwm in self.matrices.items():
            if len(key) >= 2 and key[0] == intron_type and key[1] == dinucleotides:
                matching_pwms.append(pwm)

        # If no exact matches, try fallback dinucleotides
        if not matching_pwms:
            if intron_type == 'u12':
                fallbacks = ['gtag', 'atac']
            else:  # u2
                fallbacks = ['gtag', 'gcag']

            for fallback in fallbacks:
                for key, pwm in self.matrices.items():
                    if len(key) >= 2 and key[0] == intron_type and key[1] == fallback:
                        matching_pwms.append(pwm)
                if matching_pwms:
                    break

        # Last resort: return any matrix of the right type
        if not matching_pwms:
            for key, pwm in self.matrices.items():
                if key[0] == intron_type:
                    matching_pwms.append(pwm)

        if not matching_pwms:
            raise ValueError(f"No PWM found for {intron_type}")

        return matching_pwms

    @property
    def u2_gtag(self) -> Optional[PWM]:
        """Convenience property for U2 GT-AG matrix."""
        return self.matrices.get(('u2', 'gtag'))

    @property
    def u2_gcag(self) -> Optional[PWM]:
        """Convenience property for U2 GC-AG matrix."""
        return self.matrices.get(('u2', 'gcag'))

    @property
    def u12_gtag(self) -> Optional[PWM]:
        """Convenience property for U12 GT-AG matrix."""
        return self.matrices.get(('u12', 'gtag'))

    @property
    def u12_atac(self) -> Optional[PWM]:
        """Convenience property for U12 AT-AC matrix."""
        return self.matrices.get(('u12', 'atac'))


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
    def load_from_file(filepath: Path, pseudocount: float = 0.0001) -> Dict[str, PWMSet]:
        """
        Load all PWMs from file (supports legacy .iic, YAML, and JSON formats).

        Args:
            filepath: Path to PWM file (.iic, .yaml, .yml, or .json)
            pseudocount: Pseudocount value to use for PWM scoring (default: 0.0001)

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

        Note:
            Automatically detects format based on file extension:
            - .json: JSON format (grouped)
            - .yaml or .yml: YAML format (flat)
            - Otherwise: Legacy .iic format
        """
        if not filepath.exists():
            raise FileNotFoundError(f"Matrix file not found: {filepath}")

        # Auto-detect format based on extension
        if filepath.suffix.lower() == '.json':
            matrices = PWMLoader._parse_json_file(filepath)
        elif filepath.suffix.lower() in ('.yaml', '.yml'):
            matrices = PWMLoader._parse_yaml_file(filepath)
        else:
            # Legacy .iic format
            matrices = PWMLoader._parse_matrix_file(filepath)

        # Group matrices by region (five, bp, three)
        pwm_sets = PWMLoader._group_into_pwm_sets(matrices, pseudocount)

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
    def _parse_yaml_file(filepath: Path) -> Dict[tuple, Dict[str, Dict[int, float]]]:
        """
        Parse YAML format PWM file into raw frequency dictionaries.

        Args:
            filepath: Path to YAML PWM file

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
        import yaml

        with open(filepath, 'r') as f:
            data = yaml.safe_load(f)

        matrices = {}

        # Process each matrix
        for matrix_name, matrix_data in data.get('matrices', {}).items():
            # Parse matrix name into components
            parsed_name = PWMLoader._parse_matrix_name(matrix_name)

            # Convert YAML matrix format to internal format
            # YAML format: [[A, C, G, T], [A, C, G, T], ...]
            # Internal format: {'A': {0: val, 1: val}, 'C': {0: val}, ...}

            bases = matrix_data.get('bases', ['A', 'C', 'G', 'T'])
            yaml_matrix = matrix_data.get('matrix', [])
            start_index = matrix_data.get('start_index')
            if start_index is None:
                start_index = 0

            # Build internal matrix structure
            internal_matrix = {base: {} for base in bases}

            for pos_idx, row in enumerate(yaml_matrix):
                for base_idx, base in enumerate(bases):
                    internal_matrix[base][pos_idx] = row[base_idx]

            # Store with metadata
            matrices[parsed_name] = {
                'matrix': internal_matrix,
                'start_index': start_index
            }

        return matrices

    @staticmethod
    def _parse_json_file(filepath: Path) -> Dict[tuple, Dict[str, Dict[int, float]]]:
        """
        Parse JSON format PWM file into raw frequency dictionaries.

        JSON format uses grouped matrices with metadata:
        {
            "format_version": "1.0",
            "matrix_groups": [
                {
                    "description": ["comment1", "comment2", ...],
                    "matrices": {
                        "u12_atac_five": {
                            "bases": ["A", "C", "G", "T"],
                            "start_index": -20,
                            "sample_size": 114,
                            "matrix": [[0.26, 0.21, ...], ...]
                        },
                        ...
                    }
                },
                ...
            ]
        }

        Args:
            filepath: Path to JSON PWM file

        Returns:
            Dictionary mapping (type, boundary, region) tuples to matrix data
        """
        import json

        with open(filepath, 'r') as f:
            data = json.load(f)

        matrices = {}

        # Flatten grouped matrices
        for group in data.get('matrix_groups', []):
            for matrix_name, matrix_data in group.get('matrices', {}).items():
                # Parse matrix name into components
                parsed_name = PWMLoader._parse_matrix_name(matrix_name)

                # Convert JSON matrix format to internal format
                # JSON format: [[A, C, G, T], [A, C, G, T], ...]
                # Internal format: {'A': {0: val, 1: val}, 'C': {0: val}, ...}

                bases = matrix_data.get('bases', ['A', 'C', 'G', 'T'])
                json_matrix = matrix_data.get('matrix', [])
                start_index = matrix_data.get('start_index')
                if start_index is None:
                    start_index = 0

                # Build internal matrix structure
                internal_matrix = {base: {} for base in bases}

                for pos_idx, row in enumerate(json_matrix):
                    for base_idx, base in enumerate(bases):
                        internal_matrix[base][pos_idx] = row[base_idx]

                # Store with metadata
                matrices[parsed_name] = {
                    'matrix': internal_matrix,
                    'start_index': start_index
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

        # Extract version tag if present (e.g., 'vA9', 'vA10')
        # Port from: intronIC.py:1236-1238
        # Pattern: non-letter + 'v' + optional '.' + version string
        matrix_version = re.findall(r'[^A-Za-z]v\.?([^_\s]+)', matrix_name)
        if matrix_version:
            name_bits.append(matrix_version[0])

        return tuple(name_bits)

    @staticmethod
    def _group_into_pwm_sets(
        matrices: Dict[tuple, Dict],
        pseudocount: float = 0.0001
    ) -> Dict[str, PWMSet]:
        """
        Group parsed matrices into PWMSet objects by region.

        Args:
            matrices: Dictionary from _parse_matrix_file
            pseudocount: Pseudocount value for PWM scoring

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
                start_index=data['start_index'],
                pseudocount=pseudocount
            )

            # Store by (intron_type, boundary, version?) key
            # If version tag exists (e.g., 'A9', 'A10'), include it in key
            # Otherwise use (intron_type, boundary) only
            # Port from: intronIC.py:2916-2917 (preserves full matrix_key tuple)
            version = name_tuple[3] if len(name_tuple) > 3 else None
            if version:
                # e.g., ('u12', 'gtag', 'A10'), ('u12', 'gtag', 'A9')
                key = (intron_type, boundary, version)
            else:
                # e.g., ('u12', 'gtag'), ('u12', 'atac')
                key = (intron_type, boundary)
            by_region[region][key] = pwm

        # Build PWMSets
        pwm_sets = {}
        for region in ['five', 'bp', 'three']:
            if region not in by_region:
                continue

            pwm_sets[region] = PWMSet(matrices=by_region[region])

        return pwm_sets

    @staticmethod
    def _dict_to_pwm(
        name: str,
        matrix_dict: Dict[str, Dict[int, float]],
        start_index: int,
        pseudocount: float = 0.0001
    ) -> PWM:
        """
        Convert parsed dictionary to PWM object with numpy array.

        Args:
            name: PWM name
            matrix_dict: {'A': {0: 0.25, 1: 0.30}, 'C': {...}, ...}
            start_index: Starting position index
            pseudocount: Pseudocount value for PWM scoring

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
            pseudocount=pseudocount,
            start_index=start_index
        )
