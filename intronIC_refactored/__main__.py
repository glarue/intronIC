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

# Set BLAS thread environment variables to prevent thread oversubscription
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'

import sys
from cli.main import main

if __name__ == '__main__':
    sys.exit(main())
