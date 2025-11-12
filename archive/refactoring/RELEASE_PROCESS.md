# intronIC Release Process

## Versioning

intronIC uses **Semantic Versioning 2.0.0**:

```
MAJOR.MINOR.PATCH

Examples:
- 2.0.0  # Major release (breaking changes)
- 2.1.0  # Minor release (new features, backward compatible)
- 2.0.1  # Patch release (bug fixes)
```

**Version Number in Files:**
- `pyproject.toml` → `[project] version = "2.0.0"`
- `pixi.toml` → `[workspace] version = "2.0.0"`
- `__init__.py` → Reads from package metadata automatically
- `cli/args.py` → Reads from package metadata for `--version`

## Pre-Release Checklist

Before creating a new release:

### 1. Ensure All Tests Pass
```bash
# Run full test suite
pytest

# Run with chr19 data
pytest -m requires_chr19

# Check coverage
pytest --cov=. --cov-report=term
```

### 2. Update Version Number
```bash
# Update pyproject.toml
sed -i 's/version = ".*"/version = "2.0.0"/' pyproject.toml

# Update pixi.toml
sed -i 's/version = ".*"/version = "2.0.0"/' pixi.toml

# Check consistency
python scripts/check_version.py
```

### 3. Update CHANGELOG
Add new section to `CHANGELOG.md`:

```markdown
## [2.0.0] - 2024-11-11

### Added
- New features...

### Changed
- Modified behavior...

### Fixed
- Bug fixes...

### Breaking Changes
- API changes...
```

### 4. Update README if Needed
- New features
- Changed command-line options
- Updated examples

### 5. Commit Version Bump
```bash
git add pyproject.toml pixi.toml CHANGELOG.md README.md
git commit -m "Bump version to 2.0.0"
```

## Creating a Release

### Step 1: Tag the Release
```bash
# Create annotated tag
git tag -a v2.0.0 -m "Release version 2.0.0

Major refactoring of intronIC with improved:
- Modular architecture
- Test coverage
- Documentation
- Performance

See CHANGELOG.md for full details."

# Verify tag
git tag -l -n9 v2.0.0
```

### Step 2: Build Distribution
```bash
# Install build tools (one-time)
pip install --upgrade build twine

# Clean previous builds
rm -rf dist/ build/ *.egg-info

# Build source and wheel distributions
python -m build

# Verify contents
tar -tzf dist/intronIC-2.0.0.tar.gz | head -20
unzip -l dist/intronIC-2.0.0-py3-none-any.whl | head -20
```

### Step 3: Test the Build
```bash
# Create clean test environment
python -m venv test_release_env
source test_release_env/bin/activate

# Install from local wheel
pip install dist/intronIC-2.0.0-py3-none-any.whl

# Test basic functionality
intronIC --version
intronIC --help

# Test with real data (if available)
intronIC -g test_data/*.fa.gz -a test_data/*.gff3.gz -n test -s

# Clean up
deactivate
rm -rf test_release_env
```

### Step 4: Upload to Test PyPI (Optional)
```bash
# Configure Test PyPI (one-time)
# Create account at https://test.pypi.org/
# Get API token from https://test.pypi.org/manage/account/token/

# Upload to Test PyPI
twine upload --repository testpypi dist/*

# Test installation from Test PyPI
pip install --index-url https://test.pypi.org/simple/ \
    --extra-index-url https://pypi.org/simple/ \
    intronIC

# Verify
intronIC --version
```

### Step 5: Upload to PyPI
```bash
# Configure PyPI (one-time)
# Create account at https://pypi.org/
# Get API token from https://pypi.org/manage/account/token/

# Upload to PyPI
twine upload dist/*

# Verify on PyPI
# Visit: https://pypi.org/project/intronIC/

# Test installation
pip install intronIC
intronIC --version
```

### Step 6: Push Tags and Code
```bash
# Push code
git push origin main

# Push tags
git push origin v2.0.0

# Or push all tags
git push --tags
```

### Step 7: Create GitHub Release
1. Go to https://github.com/glarue/intronIC/releases/new
2. Select tag: `v2.0.0`
3. Title: `intronIC v2.0.0`
4. Description: Copy from CHANGELOG.md
5. Attach files: `dist/intronIC-2.0.0.tar.gz`
6. Publish release

## Post-Release

### Update Development Version
After release, update version to next dev version:

```bash
# pyproject.toml
version = "2.0.1-dev"

# pixi.toml
version = "2.0.1-dev"

git commit -am "Bump version to 2.0.1-dev"
git push
```

### Announce Release
- Update project README
- Announce on relevant mailing lists/forums
- Update documentation website (if any)

## Hotfix Releases

For urgent bug fixes:

```bash
# Create hotfix branch from tag
git checkout -b hotfix-2.0.1 v2.0.0

# Make fixes
git commit -am "Fix critical bug"

# Update version to 2.0.1
# ... update files ...

# Tag and release
git tag -a v2.0.1 -m "Hotfix release 2.0.1"

# Merge back to main
git checkout main
git merge hotfix-2.0.1

# Push everything
git push origin main hotfix-2.0.1 v2.0.1
```

## Version History

| Version | Date | Type | Description |
|---------|------|------|-------------|
| 2.0.0 | 2024-11-11 | Major | Complete refactoring |
| 1.5.1 | 2020-06-15 | Patch | Bug fixes (original) |

## Rollback Procedure

If a release has issues:

```bash
# Delete tag locally
git tag -d v2.0.0

# Delete tag remotely
git push origin :refs/tags/v2.0.0

# If uploaded to PyPI, cannot delete
# Must release new version with fixes
```

## Troubleshooting

### Build Fails
```bash
# Clean everything
rm -rf dist/ build/ *.egg-info
pip install --upgrade build

# Try again
python -m build
```

### Upload Fails
```bash
# Check credentials
twine check dist/*

# Re-upload specific files
twine upload dist/intronIC-2.0.0.tar.gz
```

### Version Mismatch
```bash
# Check version consistency
python scripts/check_version.py

# Fix mismatches
vim pyproject.toml pixi.toml
```

## References

- [Semantic Versioning](https://semver.org/)
- [Python Packaging Guide](https://packaging.python.org/)
- [PyPI Upload Guide](https://packaging.python.org/tutorials/packaging-projects/)
- [Twine Documentation](https://twine.readthedocs.io/)
