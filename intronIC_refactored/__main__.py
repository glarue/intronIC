"""
Entry point for running intronIC as a module.

Usage:
    python -m intronIC_refactored [arguments]
"""

# Add package directory to sys.path for relative imports
import sys
from pathlib import Path
package_dir = Path(__file__).parent
if str(package_dir) not in sys.path:
    sys.path.insert(0, str(package_dir))

# CRITICAL: Set multiprocessing start method BEFORE any imports
# Must match original's approach (intronIC.py:5042-5046)
# IMPORTANT: Only set when running as main script, not when spawned by multiprocessing
if __name__ == '__main__':
    from multiprocessing import get_all_start_methods, get_context, set_start_method
    fork_types = get_all_start_methods()

    # Only set if not already set (avoid RuntimeError in child processes)
    try:
        if 'forkserver' in fork_types:
            set_start_method('forkserver', force=False)
        elif 'spawn' in fork_types:
            set_start_method('spawn', force=False)
    except RuntimeError:
        # Context already set (e.g., by parent process or previous import)
        pass

# NOTE: Intentionally NOT setting BLAS thread environment variables here
# to test if they're interfering with joblib's worker spawning.
# Joblib's loky backend should handle thread management automatically.

from cli.main import main

if __name__ == '__main__':
    sys.exit(main())
