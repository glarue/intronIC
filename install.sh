#!/usr/bin/env bash
#
# Simple installation script for intronIC
# Usage: ./install.sh [--dev]
#

set -e  # Exit on error

echo "======================================"
echo "intronIC Installation Script"
echo "======================================"
echo ""

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "ERROR: python3 not found. Please install Python 3.8 or later."
    exit 1
fi

# Check Python version
PYTHON_VERSION=$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')
echo "Found Python $PYTHON_VERSION"

# Check if pip is available
if ! python3 -m pip --version &> /dev/null; then
    echo "ERROR: pip not found. Please install pip."
    exit 1
fi

# Determine installation mode
if [ "$1" == "--dev" ]; then
    echo "Installing in DEVELOPMENT mode..."
    python3 -m pip install -e ".[dev]"
else
    echo "Installing in USER mode..."
    echo "(Use './install.sh --dev' for development mode)"
    python3 -m pip install .
fi

echo ""
echo "======================================"
echo "Installation Complete!"
echo "======================================"
echo ""
echo "Verify installation:"
echo "  intronIC --version"
echo "  intronIC --help"
echo ""
echo "Run a test:"
echo "  intronIC -g genome.fa -a annotation.gff3 -n species_name"
echo ""
echo "For more information, see:"
echo "  - INSTALL.md for detailed installation options"
echo "  - README.md for usage guide"
echo "  - intronIC --help for all command-line options"
echo ""
