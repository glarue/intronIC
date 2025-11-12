"""
Pytest configuration for intronIC tests.

Ensures proper imports by adding the project root to Python path.
Provides shared fixtures for common paths.
"""

import pytest
import sys
from pathlib import Path

# Add the project root to the Python path
project_root = Path(__file__).parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))


# Path fixtures
@pytest.fixture(scope="session")
def project_root():
    """Return path to project root directory."""
    return Path(__file__).parent.parent


@pytest.fixture(scope="session")
def data_dir(project_root):
    """Return path to data directory."""
    return project_root / "data"


@pytest.fixture(scope="session")
def test_data_dir(data_dir):
    """Return path to test data directory."""
    return data_dir / "test_data"


@pytest.fixture(scope="session")
def matrix_file(data_dir):
    """Return path to scoring matrices file."""
    return data_dir / "scoring_matrices.fasta.iic"
