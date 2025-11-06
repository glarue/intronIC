"""
Entry point for running intronIC as a module.

Usage:
    python -m intronIC_refactored [arguments]
"""

# CRITICAL: Set threading environment variables BEFORE any imports
# This prevents thread oversubscription when using parallelization (n_jobs > 1)
# GridSearchCV and multiprocessing Pool spawn worker processes that use BLAS,
# and multi-threaded BLAS (OpenBLAS, MKL) causes n_jobs × BLAS_threads contention.
# These settings force single-threaded BLAS in each worker for true parallelization.
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'

import sys
from cli.main import main

if __name__ == '__main__':
    sys.exit(main())
