#!/bin/bash
# Bump version across all project files
#
# Usage: ./scripts/bump_version.sh 2.0.2

set -e

if [ -z "$1" ]; then
    echo "Usage: $0 <new_version>"
    echo "Example: $0 2.0.2"
    exit 1
fi

NEW_VERSION=$1

echo "Bumping version to $NEW_VERSION..."

# Update pyproject.toml
sed -i "s/^version = \".*\"/version = \"$NEW_VERSION\"/" pyproject.toml
echo "✓ Updated pyproject.toml"

# Update pixi.toml
sed -i "s/^version = \".*\"/version = \"$NEW_VERSION\"/" pixi.toml
echo "✓ Updated pixi.toml"

echo ""
echo "Version bumped to $NEW_VERSION in all files."
echo ""
echo "Next steps:"
echo "  1. Review changes: git diff"
echo "  2. Commit: git commit -am 'chore: bump version to $NEW_VERSION'"
echo "  3. Tag: git tag v$NEW_VERSION"
echo "  4. Push: git push origin main v$NEW_VERSION"
echo "  5. Build: rm -rf dist/ && python -m build"
echo "  6. Upload: python -m twine upload dist/*"
