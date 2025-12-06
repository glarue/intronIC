#!/bin/bash
#
# Publish intronIC to PyPI
#
# Usage:
#   ./scripts/publish_pypi.sh          # Publish to TestPyPI (default)
#   ./scripts/publish_pypi.sh --prod   # Publish to production PyPI
#
# Prerequisites:
#   - PyPI API token configured in ~/.pypirc or as TWINE_PASSWORD env var
#   - pixi environment with build and twine installed
#
# Steps performed:
#   1. Check version consistency (pyproject.toml vs pixi.toml)
#   2. Clean previous builds
#   3. Build source distribution and wheel
#   4. Upload to PyPI (test or production)
#

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Determine script directory and repo root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$REPO_ROOT"

echo "========================================"
echo "  intronIC PyPI Publishing Script"
echo "========================================"
echo ""

# Parse arguments
USE_PROD=false
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --prod|--production)
            USE_PROD=true
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --prod, --production  Upload to production PyPI (default: TestPyPI)"
            echo "  --dry-run             Build but don't upload"
            echo "  -h, --help            Show this help message"
            exit 0
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            exit 1
            ;;
    esac
done

# Step 1: Check version consistency
echo -e "${YELLOW}Step 1: Checking version consistency...${NC}"
if ! pixi run python scripts/check_version.py; then
    echo -e "${RED}Version mismatch detected. Please fix before publishing.${NC}"
    exit 1
fi
echo ""

# Get version from pyproject.toml
VERSION=$(grep '^version = ' pyproject.toml | head -1 | sed 's/version = "\(.*\)"/\1/')
echo -e "Version to publish: ${GREEN}${VERSION}${NC}"
echo ""

# Step 2: Clean previous builds
echo -e "${YELLOW}Step 2: Cleaning previous builds...${NC}"
rm -rf dist/ build/ *.egg-info src/*.egg-info
echo "✓ Cleaned dist/, build/, and egg-info directories"
echo ""

# Step 3: Build
echo -e "${YELLOW}Step 3: Building source distribution and wheel...${NC}"
pixi run python -m build
echo ""

# List built files
echo "Built files:"
ls -la dist/
echo ""

# Step 4: Upload
# Note: twine is in the dev environment
if [ "$DRY_RUN" = true ]; then
    echo -e "${YELLOW}Dry run mode - skipping upload${NC}"
    echo ""
    echo "To upload manually:"
    if [ "$USE_PROD" = true ]; then
        echo "  pixi run -e dev twine upload dist/*"
    else
        echo "  pixi run -e dev twine upload --repository testpypi dist/*"
    fi
else
    if [ "$USE_PROD" = true ]; then
        echo -e "${YELLOW}Step 4: Uploading to ${RED}PRODUCTION PyPI${NC}..."
        echo ""
        echo -e "${RED}⚠️  WARNING: You are about to upload to PRODUCTION PyPI!${NC}"
        echo -e "${RED}   This cannot be undone for this version number.${NC}"
        echo ""
        read -p "Are you sure? Type 'yes' to continue: " confirm
        if [ "$confirm" != "yes" ]; then
            echo "Aborted."
            exit 1
        fi
        echo ""
        pixi run -e dev twine upload dist/*
    else
        echo -e "${YELLOW}Step 4: Uploading to TestPyPI...${NC}"
        pixi run -e dev twine upload --repository testpypi dist/*
    fi
fi

echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  Done!${NC}"
echo -e "${GREEN}========================================${NC}"

if [ "$DRY_RUN" = false ]; then
    if [ "$USE_PROD" = true ]; then
        echo ""
        echo "Package published to PyPI:"
        echo "  https://pypi.org/project/intronIC/${VERSION}/"
        echo ""
        echo "Install with:"
        echo "  pip install intronIC==${VERSION}"
    else
        echo ""
        echo "Package published to TestPyPI:"
        echo "  https://test.pypi.org/project/intronIC/${VERSION}/"
        echo ""
        echo "Install from TestPyPI with:"
        echo "  pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ intronIC==${VERSION}"
        echo ""
        echo "When ready for production, run:"
        echo "  ./scripts/publish_pypi.sh --prod"
    fi
fi
