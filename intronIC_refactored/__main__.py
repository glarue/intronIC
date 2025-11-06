"""
Entry point for running intronIC as a module.

Usage:
    python -m intronIC_refactored [arguments]
"""

# CRITICAL: Set multiprocessing start method BEFORE any imports
# This prevents thread oversubscription when using parallelization (n_jobs > 1).
#
# Problem: Default 'fork' method copies parent process memory including initialized
# BLAS state with multi-threading. When GridSearchCV spawns workers, each inherits
# multi-threaded BLAS, causing n_jobs × BLAS_threads competing threads.
#
# Solution: Use 'forkserver' or 'spawn' which create fresh processes without
# inheriting BLAS state, enabling true parallelization.
#
# See: https://scikit-learn.org/stable/computing/parallelism.html
import multiprocessing as mp

# Must be called before any other code that uses multiprocessing
# Note: This is the same approach used in the original intronIC.py:5042-5046
try:
    available_methods = mp.get_all_start_methods()
    if 'forkserver' in available_methods:
        mp.set_start_method('forkserver', force=True)
    elif 'spawn' in available_methods:
        mp.set_start_method('spawn', force=True)
except RuntimeError:
    # Already set, ignore
    pass

# Secondary defense: Set thread environment variables
# (Though forkserver/spawn should handle this, these provide extra safety)
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'

import sys
from cli.main import main

if __name__ == '__main__':
    sys.exit(main())
