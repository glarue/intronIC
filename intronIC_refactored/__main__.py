"""
Entry point for running intronIC as a module.

Usage:
    python -m intronIC_refactored [arguments]
"""

# CRITICAL: Set multiprocessing start method BEFORE any imports
# Must match original's approach (intronIC.py:5042-5046)
from multiprocessing import get_all_start_methods, set_start_method
fork_types = get_all_start_methods()
if 'forkserver' in fork_types:
    set_start_method('forkserver')
elif 'spawn' in fork_types:
    set_start_method('spawn')

# NOTE: Intentionally NOT setting BLAS thread environment variables here
# to test if they're interfering with joblib's worker spawning.
# Joblib's loky backend should handle thread management automatically.

import sys
from cli.main import main

if __name__ == '__main__':
    sys.exit(main())
