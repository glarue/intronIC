#!/usr/bin/env python3
"""
Check that version is consistent across configuration files.

This ensures pyproject.toml and pixi.toml have the same version.
"""

import sys
from pathlib import Path

try:
    import tomli
except ImportError:
    import tomllib as tomli  # Python 3.11+

def check_versions():
    """Check version consistency across files."""
    repo_root = Path(__file__).parent.parent

    # Read pyproject.toml
    pyproject_path = repo_root / "pyproject.toml"
    with open(pyproject_path, "rb") as f:
        pyproject = tomli.load(f)
        pyproject_version = pyproject["project"]["version"]

    # Read pixi.toml
    pixi_path = repo_root / "pixi.toml"
    with open(pixi_path, "rb") as f:
        pixi = tomli.load(f)
        pixi_version = pixi["workspace"]["version"]

    # Check consistency
    if pyproject_version != pixi_version:
        print("❌ ERROR: Version mismatch!")
        print(f"  pyproject.toml: {pyproject_version}")
        print(f"  pixi.toml:      {pixi_version}")
        print()
        print("Please update both files to have the same version.")
        return False

    print(f"✓ Version consistent: {pyproject_version}")
    return True

if __name__ == "__main__":
    success = check_versions()
    sys.exit(0 if success else 1)
