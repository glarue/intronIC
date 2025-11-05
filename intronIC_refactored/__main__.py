"""
Entry point for running intronIC as a module.

Usage:
    python -m intronIC_refactored [arguments]
"""

import sys
from cli.main import main

if __name__ == '__main__':
    sys.exit(main())
