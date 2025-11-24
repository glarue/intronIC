#!/usr/bin/env python3
"""
Converter between legacy .iic PWM format, YAML, and JSON formats.

Usage:
    # Convert legacy to JSON (recommended)
    python convert_pwm_format.py scoring_matrices.fasta.iic intronIC_scoring_PWMs.json

    # Merge multiple legacy files into one JSON
    python convert_pwm_format.py --merge file1.iic file2.iic output.json

    # Convert JSON to legacy
    python convert_pwm_format.py intronIC_scoring_PWMs.json scoring_matrices.iic

    # Convert legacy to YAML (legacy support)
    python convert_pwm_format.py scoring_matrices.fasta.iic scoring_matrices.yaml
"""

import sys
import json
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
import yaml


def parse_legacy_pwm_grouped(filepath: Path) -> Dict[str, Any]:
    """
    Parse legacy .iic PWM format, preserving comment-to-matrix group associations.

    Returns:
        Dictionary with format_version and matrix_groups
    """
    result = {
        'format_version': '1.0',
        'matrix_groups': []
    }

    current_group_comments = []
    current_group_matrices = {}

    current_name = None
    current_matrix = []
    current_start_index = None
    current_sample_size = None
    current_metadata_str = None
    in_header = True

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()

            # Skip empty lines
            if not line:
                continue

            # Comment lines - start new group if we have matrices already
            if line.startswith('#'):
                comment = line[1:].strip()

                # If we already have matrices, save current matrix then group
                if current_group_matrices or (current_name is not None and current_matrix):
                    # Save current matrix if exists
                    if current_name is not None and current_matrix:
                        matrix_data = {
                            'bases': ['A', 'C', 'G', 'T'],
                            'matrix': current_matrix
                        }
                        if current_start_index is not None:
                            matrix_data['start_index'] = current_start_index
                        if current_sample_size is not None:
                            matrix_data['sample_size'] = current_sample_size
                        if current_metadata_str is not None:
                            matrix_data['metadata'] = current_metadata_str

                        current_group_matrices[current_name] = matrix_data
                        current_name = None
                        current_matrix = []

                    # Save current group and start new one
                    result['matrix_groups'].append({
                        'description': current_group_comments,
                        'matrices': current_group_matrices
                    })
                    current_group_comments = []
                    current_group_matrices = {}

                current_group_comments.append(comment)
                continue

            # Header line: >matrix_name
            if line.startswith('>'):
                # Save previous matrix if exists
                if current_name is not None and current_matrix:
                    matrix_data = {
                        'bases': ['A', 'C', 'G', 'T'],
                        'matrix': current_matrix
                    }
                    if current_start_index is not None:
                        matrix_data['start_index'] = current_start_index
                    if current_sample_size is not None:
                        matrix_data['sample_size'] = current_sample_size
                    if current_metadata_str is not None:
                        matrix_data['metadata'] = current_metadata_str

                    current_group_matrices[current_name] = matrix_data

                # Parse new header
                header_parts = line[1:].split()
                current_name = header_parts[0]

                # Extract metadata from header
                current_start_index = None
                current_sample_size = None
                current_metadata_str = None

                for part in header_parts[1:]:
                    if 'start=' in part:
                        current_start_index = int(part.split('=')[1])
                    elif part.startswith('(n=') and part.endswith(')'):
                        current_sample_size = int(part[3:-1])
                    elif part.startswith('(') and part.endswith(')'):
                        # Other metadata like (mercer2015)
                        current_metadata_str = part[1:-1]

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
            matrix_data = {
                'bases': ['A', 'C', 'G', 'T'],
                'matrix': current_matrix
            }
            if current_start_index is not None:
                matrix_data['start_index'] = current_start_index
            if current_sample_size is not None:
                matrix_data['sample_size'] = current_sample_size
            if current_metadata_str is not None:
                matrix_data['metadata'] = current_metadata_str

            current_group_matrices[current_name] = matrix_data

        # Save last group
        if current_group_matrices:
            result['matrix_groups'].append({
                'description': current_group_comments,
                'matrices': current_group_matrices
            })

    return result


def parse_legacy_pwm_flat(filepath: Path) -> Dict[str, Any]:
    """
    Parse legacy .iic PWM format without preserving groups (for backward compatibility).

    Returns:
        Dictionary with metadata and matrices (flat structure)
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
                freqs = [float(x) for x in line.split()]
                if len(freqs) == 4:
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


def write_json_pwm(data: Dict[str, Any], filepath: Path):
    """Write PWM data to JSON format with custom formatting for visual alignment."""

    # Custom JSON encoder for pretty-printing with aligned columns
    def format_matrix_row(row: List[float]) -> str:
        """Format a matrix row with fixed-width decimals for alignment."""
        return '[' + ', '.join(f'{val:.10f}' for val in row) + ']'

    with open(filepath, 'w') as f:
        f.write('{\n')
        f.write('  "format_version": "1.0",\n')
        f.write('  "matrix_groups": [\n')

        for group_idx, group in enumerate(data['matrix_groups']):
            is_last_group = group_idx == len(data['matrix_groups']) - 1

            f.write('    {\n')

            # Write description (comments)
            f.write('      "description": [\n')
            for comment_idx, comment in enumerate(group['description']):
                is_last_comment = comment_idx == len(group['description']) - 1
                comma = '' if is_last_comment else ','
                # Escape special JSON characters
                escaped = comment.replace('\\', '\\\\').replace('"', '\\"')
                f.write(f'        "{escaped}"{comma}\n')
            f.write('      ],\n')

            # Write matrices
            f.write('      "matrices": {\n')
            matrix_items = list(group['matrices'].items())
            for matrix_idx, (name, matrix_data) in enumerate(matrix_items):
                is_last_matrix = matrix_idx == len(matrix_items) - 1

                f.write(f'        "{name}": {{\n')
                f.write(f'          "bases": {json.dumps(matrix_data["bases"])},\n')

                # Optional fields
                if 'start_index' in matrix_data:
                    f.write(f'          "start_index": {matrix_data["start_index"]},\n')
                if 'sample_size' in matrix_data:
                    f.write(f'          "sample_size": {matrix_data["sample_size"]},\n')
                if 'metadata' in matrix_data:
                    f.write(f'          "metadata": "{matrix_data["metadata"]}",\n')

                # Write matrix with aligned columns
                f.write('          "matrix": [\n')
                for row_idx, row in enumerate(matrix_data['matrix']):
                    is_last_row = row_idx == len(matrix_data['matrix']) - 1
                    comma = '' if is_last_row else ','
                    f.write(f'            {format_matrix_row(row)}{comma}\n')
                f.write('          ]\n')

                matrix_comma = '' if is_last_matrix else ','
                f.write(f'        }}{matrix_comma}\n')
            f.write('      }\n')

            group_comma = '' if is_last_group else ','
            f.write(f'    }}{group_comma}\n')

        f.write('  ]\n')
        f.write('}\n')


def write_yaml_pwm(data: Dict[str, Any], filepath: Path):
    """Write PWM data to YAML format (legacy support)."""

    # Convert grouped format to flat format for YAML
    if 'matrix_groups' in data:
        flat_data = {
            'metadata': {
                'format_version': data.get('format_version', '1.0'),
                'notes': []
            },
            'matrices': {}
        }

        # Flatten groups
        for group in data['matrix_groups']:
            flat_data['metadata']['notes'].extend(group['description'])
            flat_data['matrices'].update(group['matrices'])

        data = flat_data

    # Create output structure
    output = {
        'format_version': data.get('format_version', data.get('metadata', {}).get('format_version', '1.0')),
        'metadata': {
            'description': 'Position Weight Matrices for intronIC U12/U2 intron classification',
            'sources': data.get('metadata', {}).get('sources', []),
            'notes': data.get('metadata', {}).get('notes', [])
        },
        'matrices': {}
    }

    # Convert matrices
    matrices = data.get('matrices', {})
    for name, matrix_data in matrices.items():
        matrix_output = {
            'bases': matrix_data.get('bases', ['A', 'C', 'G', 'T']),
            'matrix': matrix_data['matrix']
        }

        if matrix_data.get('start_index') is not None:
            matrix_output['start_index'] = matrix_data['start_index']

        if matrix_data.get('sample_size') is not None:
            matrix_output['sample_size'] = matrix_data['sample_size']

        output['matrices'][name] = matrix_output

    # Write with custom formatting
    class CompactDumper(yaml.SafeDumper):
        """Custom YAML dumper that uses flow style for lists."""
        pass

    def represent_list(dumper, data):
        """Represent lists in flow style (e.g., [A, C, G, T])."""
        return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)

    CompactDumper.add_representer(list, represent_list)

    with open(filepath, 'w') as f:
        f.write(f"# intronIC Position Weight Matrices (YAML format)\n")
        f.write(f"# Generated from: {filepath.stem}\n\n")
        yaml.dump(output, f, Dumper=CompactDumper, default_flow_style=False, sort_keys=False, width=120)


def write_legacy_pwm(data: Dict[str, Any], filepath: Path):
    """Write PWM data to legacy .iic format."""

    with open(filepath, 'w') as f:
        # Handle grouped format
        if 'matrix_groups' in data:
            for group_idx, group in enumerate(data['matrix_groups']):
                # Write group comments
                for comment in group['description']:
                    f.write(f"# {comment}\n")

                # Write matrices in this group
                for name, matrix_data in group['matrices'].items():
                    # Header line
                    header_parts = [f">{name}"]

                    if matrix_data.get('start_index') is not None:
                        header_parts.append(f"start={matrix_data['start_index']}")

                    if matrix_data.get('metadata') is not None:
                        header_parts.append(f"({matrix_data['metadata']})")

                    if matrix_data.get('sample_size') is not None:
                        header_parts.append(f"(n={matrix_data['sample_size']})")

                    f.write('\t'.join(header_parts) + '\n')

                    # Base header
                    bases = matrix_data.get('bases', ['A', 'C', 'G', 'T'])
                    f.write('\t'.join(bases) + '\n')

                    # Matrix rows
                    for row in matrix_data['matrix']:
                        f.write('\t'.join(str(x) for x in row) + '\n')

        # Handle flat format (legacy)
        else:
            # Write metadata comments
            metadata = data.get('metadata', {})
            for source in metadata.get('sources', []):
                f.write(f"# {source}\n")

            for note in metadata.get('notes', []):
                f.write(f"# {note}\n")

            # Write matrices
            for name, matrix_data in data.get('matrices', {}).items():
                header_parts = [f">{name}"]

                if matrix_data.get('start_index') is not None:
                    header_parts.append(f"start={matrix_data['start_index']}")

                if matrix_data.get('sample_size') is not None:
                    header_parts.append(f"(n={matrix_data['sample_size']})")

                f.write('\t'.join(header_parts) + '\n')

                bases = matrix_data.get('bases', ['A', 'C', 'G', 'T'])
                f.write('\t'.join(bases) + '\n')

                for row in matrix_data['matrix']:
                    f.write('\t'.join(str(x) for x in row) + '\n')


def merge_legacy_files(input_files: List[Path], output_file: Path):
    """Merge multiple legacy .iic files into one JSON file."""

    merged_data = {
        'format_version': '1.0',
        'matrix_groups': []
    }

    for input_file in input_files:
        print(f"  Reading {input_file}...")
        grouped_data = parse_legacy_pwm_grouped(input_file)
        merged_data['matrix_groups'].extend(grouped_data['matrix_groups'])

    print(f"\nWriting merged output to {output_file}...")
    write_json_pwm(merged_data, output_file)

    total_matrices = sum(len(g['matrices']) for g in merged_data['matrix_groups'])
    print(f"  ✓ Merged {len(merged_data['matrix_groups'])} groups, {total_matrices} matrices")


def convert(input_file: Path, output_file: Path):
    """
    Automatically detect format and convert.

    Args:
        input_file: Input PWM file (.iic, .json, .yaml/.yml)
        output_file: Output PWM file (.iic, .json, .yaml/.yml)
    """
    input_ext = input_file.suffix.lower()
    output_ext = output_file.suffix.lower()

    print(f"Converting {input_file} → {output_file}")

    # Determine input format
    if input_ext == '.json':
        with open(input_file, 'r') as f:
            data = json.load(f)
        total_matrices = sum(len(g['matrices']) for g in data.get('matrix_groups', []))
        print(f"  Loaded JSON format ({len(data.get('matrix_groups', []))} groups, {total_matrices} matrices)")

    elif input_ext in ('.yaml', '.yml'):
        with open(input_file, 'r') as f:
            data = yaml.safe_load(f)
        print(f"  Loaded YAML format ({len(data.get('matrices', {}))} matrices)")

    else:
        # Legacy .iic format - use grouped parsing for JSON output, flat for others
        if output_ext == '.json':
            data = parse_legacy_pwm_grouped(input_file)
            total_matrices = sum(len(g['matrices']) for g in data['matrix_groups'])
            print(f"  Loaded legacy .iic format ({len(data['matrix_groups'])} groups, {total_matrices} matrices)")
        else:
            data = parse_legacy_pwm_flat(input_file)
            print(f"  Loaded legacy .iic format ({len(data['matrices'])} matrices)")

    # Write output format
    if output_ext == '.json':
        write_json_pwm(data, output_file)
        print(f"  ✓ Written JSON format")
    elif output_ext in ('.yaml', '.yml'):
        write_yaml_pwm(data, output_file)
        print(f"  ✓ Written YAML format")
    else:
        write_legacy_pwm(data, output_file)
        print(f"  ✓ Written legacy .iic format")

    # Print summary
    if 'matrix_groups' in data:
        print(f"\nMatrix groups:")
        for group_idx, group in enumerate(data['matrix_groups'], 1):
            print(f"  Group {group_idx}: {len(group['matrices'])} matrices")
            for name, matrix in group['matrices'].items():
                matrix_shape = len(matrix['matrix'])
                print(f"    - {name}: {matrix_shape} positions")
    else:
        print(f"\nMatrices converted:")
        for name, matrix in data.get('matrices', {}).items():
            matrix_shape = len(matrix['matrix'])
            print(f"  - {name}: {matrix_shape} positions")


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        print("\nError: Requires at least 2 arguments")
        sys.exit(1)

    # Check for --merge flag
    if sys.argv[1] == '--merge':
        if len(sys.argv) < 4:
            print("Error: --merge requires at least 2 input files and 1 output file")
            sys.exit(1)

        input_files = [Path(arg) for arg in sys.argv[2:-1]]
        output_file = Path(sys.argv[-1])

        for input_file in input_files:
            if not input_file.exists():
                print(f"Error: Input file not found: {input_file}")
                sys.exit(1)

        print(f"Merging {len(input_files)} files → {output_file}")
        merge_legacy_files(input_files, output_file)

    else:
        input_file = Path(sys.argv[1])
        output_file = Path(sys.argv[2])

        if not input_file.exists():
            print(f"Error: Input file not found: {input_file}")
            sys.exit(1)

        convert(input_file, output_file)


if __name__ == '__main__':
    main()
