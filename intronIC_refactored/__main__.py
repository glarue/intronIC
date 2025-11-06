"""
Entry point for running intronIC as a module.

Usage:
    python -m intronIC_refactored [arguments]
"""

# Set BLAS thread environment variables BEFORE any imports
# This prevents thread oversubscription: n_jobs × BLAS_threads competing threads
# Must be set before numpy/sklearn initialization
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'

# Let joblib/sklearn handle multiprocessing backend selection
# Joblib's 'loky' backend (default in modern sklearn) handles process spawning correctly

import sys
from cli.main import main

if __name__ == '__main__':
    sys.exit(main())
