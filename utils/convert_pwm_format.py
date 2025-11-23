#!/usr/bin/env python3
"""
Bidirectional converter between legacy .iic PWM format and YAML format.

Usage:
    # Convert legacy to YAML
    python convert_pwm_format.py scoring_matrices.fasta.iic scoring_matrices.yaml

    # Convert YAML to legacy
    python convert_pwm_format.py scoring_matrices.yaml scoring_matrices.fasta.iic
"""

import sys
from pathlib import Path
from typing import Dict, List, Any
import yaml


def parse_legacy_pwm(filepath: Path) -> Dict[str, Any]:
    """
    Parse legacy .iic PWM format.

    Returns:
        Dictionary with metadata and matrices
    """
    result = {
        'metadata': {
            'format_version': '1.0',
            'sources': [],
            'notes': []
        },
        'matrices': {}
    }

    current_name = None
    current_matrix = []
    current_start_index = None
    current_sample_size = None
    in_header = True

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()

            # Skip empty lines
            if not line:
                continue

            # Parse comment lines for metadata
            if line.startswith('#'):
                comment = line[1:].strip()
                if 'External sources:' in comment or 'Internal sources:' in comment:
                    result['metadata']['sources'].append(comment)
                else:
                    result['metadata']['notes'].append(comment)
                continue

            # Header line: >matrix_name
            if line.startswith('>'):
                # Save previous matrix if exists
                if current_name is not None and current_matrix:
                    result['matrices'][current_name] = {
                        'start_index': current_start_index,
                        'sample_size': current_sample_size,
                        'bases': ['A', 'C', 'G', 'T'],
                        'matrix': current_matrix
                    }

                # Parse new header
                header_parts = line[1:].split()
                current_name = header_parts[0]

                # Extract metadata from header
                current_start_index = None
                current_sample_size = None

                for part in header_parts:
                    if 'start=' in part:
                        current_start_index = int(part.split('=')[1])
                    elif part.startswith('(n=') and part.endswith(')'):
                        current_sample_size = int(part[3:-1])

                # Reset matrix
                current_matrix = []
                in_header = True
                continue

            # Base header line (A C G T)
            if line in ('A\tC\tG\tT', 'A C G T'):
                in_header = False
                continue

            # Matrix data rows
            if not in_header and current_name is not None:
                # Parse frequency row
                freqs = [float(x) for x in line.split()]
                if len(freqs) == 4:  # Valid row
                    current_matrix.append(freqs)

        # Save last matrix
        if current_name is not None and current_matrix:
            result['matrices'][current_name] = {
                'start_index': current_start_index,
                'sample_size': current_sample_size,
                'bases': ['A', 'C', 'G', 'T'],
                'matrix': current_matrix
            }

    return result


def write_yaml_pwm(data: Dict[str, Any], filepath: Path):
    """Write PWM data to YAML format."""

    # Create clean output structure
    output = {
        'format_version': data['metadata'].get('format_version', '1.0'),
        'metadata': {
            'description': 'Position Weight Matrices for intronIC U12/U2 intron classification',
            'sources': data['metadata'].get('sources', []),
            'notes': data['metadata'].get('notes', [])
        },
        'matrices': {}
    }

    # Convert matrices
    for name, matrix_data in data['matrices'].items():
        matrix_output = {
            'bases': matrix_data.get('bases', ['A', 'C', 'G', 'T']),
            'matrix': matrix_data['matrix']
        }

        # Only include start_index if it's not None
        if matrix_data.get('start_index') is not None:
            matrix_output['start_index'] = matrix_data['start_index']

        # Only include sample_size if it's not None
        if matrix_data.get('sample_size') is not None:
            matrix_output['sample_size'] = matrix_data['sample_size']

        output['matrices'][name] = matrix_output

    # Write with custom formatting for readability
    with open(filepath, 'w') as f:
        f.write(f"# intronIC Position Weight Matrices (YAML format)\n")
        f.write(f"# Generated from: {Path(filepath).stem}\n\n")
        yaml.dump(output, f, default_flow_style=False, sort_keys=False, width=120)


def write_legacy_pwm(data: Dict[str, Any], filepath: Path):
    """Write PWM data to legacy .iic format."""

    with open(filepath, 'w') as f:
        # Write metadata comments
        for source in data['metadata'].get('sources', []):
            f.write(f"# {source}\n")

        for note in data['metadata'].get('notes', []):
            f.write(f"# {note}\n")

        # Write matrices
        for name, matrix_data in data['matrices'].items():
            # Header line
            header_parts = [f">{name}"]

            if matrix_data.get('start_index') is not None:
                header_parts.append(f"start={matrix_data['start_index']}")

            if matrix_data.get('sample_size') is not None:
                header_parts.append(f"(n={matrix_data['sample_size']})")

            f.write('\t'.join(header_parts) + '\n')

            # Base header
            bases = matrix_data.get('bases', ['A', 'C', 'G', 'T'])
            f.write('\t'.join(bases) + '\n')

            # Matrix rows
            for row in matrix_data['matrix']:
                f.write('\t'.join(str(x) for x in row) + '\n')


def convert(input_file: Path, output_file: Path):
    """
    Automatically detect format and convert.

    Args:
        input_file: Input PWM file (.iic or .yaml/.yml)
        output_file: Output PWM file (.iic or .yaml/.yml)
    """
    input_ext = input_file.suffix.lower()
    output_ext = output_file.suffix.lower()

    print(f"Converting {input_file} → {output_file}")

    # Determine input format
    if input_ext in ('.yaml', '.yml'):
        # Load YAML
        with open(input_file, 'r') as f:
            data = yaml.safe_load(f)
        print(f"  Loaded YAML format ({len(data['matrices'])} matrices)")
    else:
        # Assume legacy .iic format
        data = parse_legacy_pwm(input_file)
        print(f"  Loaded legacy .iic format ({len(data['matrices'])} matrices)")

    # Write output format
    if output_ext in ('.yaml', '.yml'):
        write_yaml_pwm(data, output_file)
        print(f"  ✓ Written YAML format")
    else:
        write_legacy_pwm(data, output_file)
        print(f"  ✓ Written legacy .iic format")

    print(f"\nMatrices converted:")
    for name, matrix in data['matrices'].items():
        matrix_shape = len(matrix['matrix']) if matrix['matrix'] else 0
        print(f"  - {name}: {matrix_shape} positions")


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        print("\nError: Requires exactly 2 arguments")
        sys.exit(1)

    input_file = Path(sys.argv[1])
    output_file = Path(sys.argv[2])

    if not input_file.exists():
        print(f"Error: Input file not found: {input_file}")
        sys.exit(1)

    convert(input_file, output_file)


if __name__ == '__main__':
    main()
