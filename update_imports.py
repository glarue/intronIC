#!/usr/bin/env python3
"""
Script to update all imports from flat layout to src layout.

Changes:
  from intronIC.cli.X import Y          -> from intronIC.cli.X import Y
  from intronIC.extraction.X import Y   -> from intronIC.extraction.X import Y
  import intronIC.cli.X                 -> import intronIC.cli.X
  etc.
"""

import re
import sys
from pathlib import Path

# Modules that moved to src/intronIC/
MODULES = [
    'intronIC.cli',
    'core',
    'extraction',
    'file_io',
    'scoring',
    'classification',
    'visualization',
    'utils',
    'data',
]

def update_imports_in_file(filepath: Path, dry_run: bool = False) -> tuple[int, list[str]]:
    """Update imports in a single file. Returns (num_changes, changed_lines)."""

    try:
        content = filepath.read_text()
    except Exception as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)
        return 0, []

    original_content = content
    changes = []

    # Pattern 1: from MODULE.X import Y
    # Matches: from cli.main import run
    # Result:  from intronIC.cli.main import run
    for module in MODULES:
        pattern1 = rf'^(\s*)from ({module})(\.|\s+import)'
        replacement1 = rf'\1from intronIC.\2\3'

        new_content, count = re.subn(pattern1, replacement1, content, flags=re.MULTILINE)
        if count > 0:
            changes.append(f"  Pattern 1 (from {module}.*): {count} replacements")
            content = new_content

    # Pattern 2: import MODULE or import MODULE.X
    # Matches: import cli.main
    # Result:  import intronIC.cli.main
    for module in MODULES:
        pattern2 = rf'^(\s*)import ({module})(\.|$|\s)'
        replacement2 = rf'\1import intronIC.\2\3'

        new_content, count = re.subn(pattern2, replacement2, content, flags=re.MULTILINE)
        if count > 0:
            changes.append(f"  Pattern 2 (import {module}.*): {count} replacements")
            content = new_content

    # Pattern 3: importlib.import_module("MODULE") or similar string references
    for module in MODULES:
        pattern3 = rf'''(['"])({module})(\.|\1)'''

        def replace_if_import_context(match):
            """Only replace if it looks like an import string."""
            full_match = match.group(0)
            # Check if preceded by import-related keywords
            start_pos = content.find(full_match)
            if start_pos > 20:
                context = content[start_pos-20:start_pos]
                if 'import' in context.lower() or 'module' in context.lower():
                    return match.group(1) + 'intronIC.' + match.group(2) + match.group(3)
            return full_match

        new_content = re.sub(pattern3, replace_if_import_context, content)
        if new_content != content:
            changes.append(f"  Pattern 3 (string '{module}'): updated")
            content = new_content

    if content != original_content:
        if not dry_run:
            try:
                filepath.write_text(content)
            except Exception as e:
                print(f"Error writing {filepath}: {e}", file=sys.stderr)
                return 0, []
        return len(changes), changes

    return 0, []


def main():
    """Update imports in all Python files."""

    dry_run = '--dry-run' in sys.argv

    if dry_run:
        print("DRY RUN MODE - No files will be modified\n")

    # Paths to search
    search_paths = [
        Path('src/intronIC'),
        Path('tests'),
        Path('scripts'),
        Path('.'),  # Root level (debug_cycles.py, etc.)
    ]

    total_files = 0
    total_changes = 0

    for search_path in search_paths:
        if not search_path.exists():
            continue

        # Get all .py files
        if search_path == Path('.'):
            # Only root-level .py files, not recursive
            py_files = [f for f in search_path.glob('*.py') if f.is_file()]
        else:
            py_files = list(search_path.rglob('*.py'))

        for filepath in py_files:
            # Skip __pycache__ and .pixi directories
            if '__pycache__' in str(filepath) or '.pixi' in str(filepath):
                continue

            num_changes, changes = update_imports_in_file(filepath, dry_run=dry_run)

            if num_changes > 0:
                total_files += 1
                total_changes += num_changes
                print(f"\n{filepath}:")
                for change in changes:
                    print(change)

    print(f"\n{'[DRY RUN] ' if dry_run else ''}Summary:")
    print(f"  Files modified: {total_files}")
    print(f"  Total changes: {total_changes}")

    if dry_run:
        print("\nRun without --dry-run to apply changes")


if __name__ == '__main__':
    main()
