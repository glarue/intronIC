"""
Pytest configuration for intronIC tests.

Ensures proper imports by adding the project root to Python path.
"""

import sys
from pathlib import Path

# Add the project root to the Python path
project_root = Path(__file__).parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))
