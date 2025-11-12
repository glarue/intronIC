"""
intronIC - Refactored version.

A bioinformatics tool for extracting and classifying introns from genomic data.
"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("intronIC")
except PackageNotFoundError:
    # Package not installed yet (e.g., during development)
    __version__ = "2.0.0-dev (not installed)"
