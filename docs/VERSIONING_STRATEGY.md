# intronIC Versioning Strategy

## Current Situation

### Original intronIC
- Uses **versioneer** (old tool, last updated 2018)
- Version comes from git tags
- Generates `_version.py` file
- Version in setup.cfg: `1.5.1`

### Refactored intronIC
- Static version in `pyproject.toml`: `1.5.1`
- Static version in `pixi.toml`: `1.5.1`
- Static version in `__init__.py`: `2.0.0-dev`
- **No automatic versioning from git tags**

## Problem

Three different version strings in the codebase! This leads to:
- Confusion about actual version
- Manual updates required in multiple places
- Risk of inconsistency
- Not following DRY principle

## Modern Solutions (2024)

Versioneer is **outdated**. Modern alternatives:

### 1. **setuptools-scm** (Recommended ✓)
- Most popular modern replacement for versioneer
- Actively maintained (2024)
- Works with pyproject.toml
- Derives version from git tags automatically
- Supports PEP 440 version schemes
- Works with hatchling, setuptools, etc.

### 2. **hatch-vcs**
- Hatchling's version plugin
- Similar to setuptools-scm
- Integrates with hatchling backend

### 3. **Static version** (Simplest)
- Manually update version in pyproject.toml
- Simple but error-prone
- Requires discipline in release process

### 4. **Bump2version / bump-my-version**
- CLI tool to bump versions
- Updates all files at once
- Still requires manual triggering

## Recommendation: setuptools-scm

**Why?**
- Industry standard in 2024
- Automatic versioning from git tags
- Supports development versions (e.g., `2.0.0.dev42+g1234abc`)
- Single source of truth (git tags)
- Works with PyPI publication
- Minimal configuration

**How it works:**
1. You tag a release: `git tag v2.0.0`
2. setuptools-scm reads the tag
3. Version is automatically set to `2.0.0`
4. Between tags: version is `2.0.0.dev42+g1234abc` (dev version)
5. Published wheels have correct version embedded

## Implementation Plan

### Option A: setuptools-scm (Recommended)

**1. Update pyproject.toml:**
```toml
[build-system]
requires = ["hatchling", "hatch-vcs"]
build-backend = "hatchling.build"

[project]
name = "intronIC"
dynamic = ["version"]  # Version comes from git
description = "..."
# ... rest of metadata

[tool.hatch.version]
source = "vcs"

[tool.hatch.build.hooks.vcs]
version-file = "_version.py"
```

**2. Add .git_archival.txt** (for source distributions):
```
ref-names: $Format:%D$
```

**3. Tag current state:**
```bash
git tag v2.0.0
```

**4. Version available in code:**
```python
from importlib.metadata import version
__version__ = version("intronIC")
```

**5. Development versions:**
- `v2.0.0` tag → version is `2.0.0`
- 5 commits after tag → version is `2.0.1.dev5+gabc1234`
- Dirty working tree → version is `2.0.1.dev5+gabc1234.d20241111`

### Option B: Simple Static Version (Simpler but manual)

**1. Single source in pyproject.toml:**
```toml
[project]
version = "2.0.0"
```

**2. Remove from other files:**
- Remove version from `pixi.toml` (not needed)
- Update `__init__.py` to read from package metadata

**3. Create version.py:**
```python
# version.py
from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("intronIC")
except PackageNotFoundError:
    __version__ = "unknown"
```

**4. Use in __init__.py:**
```python
from .version import __version__
```

**5. Bump version manually:**
Edit `pyproject.toml` before each release.

### Option C: Hybrid (Static with validation)

Keep static version but add consistency check:

**1. pyproject.toml** (source of truth):
```toml
[project]
version = "2.0.0"
```

**2. __init__.py:**
```python
from importlib.metadata import version
__version__ = version("intronIC")
```

**3. pixi.toml** (synced):
```toml
version = "2.0.0"  # Keep in sync manually
```

## Comparison

| Feature | setuptools-scm | Static | Hybrid |
|---------|----------------|--------|--------|
| Auto from git tags | ✓ | ✗ | ✗ |
| Dev versions | ✓ | ✗ | ✗ |
| Simplicity | Medium | High | Medium |
| Error-prone | Low | High | Medium |
| PyPI-ready | ✓ | ✓ | ✓ |
| Maintenance | Low | High | Medium |

## Versioning Scheme

Following **Semantic Versioning 2.0.0**:

```
MAJOR.MINOR.PATCH[-PRERELEASE][+BUILD]

Examples:
- 2.0.0         # Stable release
- 2.0.0-alpha.1 # Alpha pre-release
- 2.0.0-beta.2  # Beta pre-release
- 2.0.0-rc.1    # Release candidate
- 2.0.1         # Patch release
- 2.1.0         # Minor release (new features)
- 3.0.0         # Major release (breaking changes)
```

**Recommended versions for intronIC:**
- Current original: `1.5.1`
- First refactored release: `2.0.0` (major refactoring = breaking change)
- Development: `2.0.0.dev0+...` (if using setuptools-scm)

## Release Workflow

### With setuptools-scm:
```bash
# 1. Ensure all changes committed
git add .
git commit -m "Prepare for v2.0.0 release"

# 2. Tag the release
git tag -a v2.0.0 -m "Release version 2.0.0"

# 3. Build distributions
python -m build

# 4. Upload to PyPI
twine upload dist/*

# 5. Push tag
git push origin v2.0.0
```

### With static version:
```bash
# 1. Update version in pyproject.toml
sed -i 's/version = ".*"/version = "2.0.0"/' pyproject.toml

# 2. Update pixi.toml if needed
sed -i 's/version = ".*"/version = "2.0.0"/' pixi.toml

# 3. Commit version bump
git add pyproject.toml pixi.toml
git commit -m "Bump version to 2.0.0"

# 4. Tag the release
git tag -a v2.0.0 -m "Release version 2.0.0"

# 5. Build and upload
python -m build
twine upload dist/*

# 6. Push everything
git push && git push --tags
```

## Recommended Choice: Static Version (for now)

**Reasoning:**
1. **Simplicity**: Easy to understand and maintain
2. **Explicit**: Clear what version is being released
3. **No dependencies**: Doesn't require setuptools-scm
4. **Works everywhere**: pip, pixi, both work identically
5. **Can migrate later**: Easy to switch to setuptools-scm later if needed

**Implementation:**
1. Keep `version = "2.0.0"` in `pyproject.toml` (single source of truth)
2. Update `__init__.py` to read from package metadata
3. Update `pixi.toml` to match (for display purposes)
4. Add version command: `intronIC --version`
5. Add pre-commit hook to check version consistency

## Making Version Available

**In CLI:**
```python
# cli/main.py
from importlib.metadata import version

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--version',
        action='version',
        version=f'intronIC {version("intronIC")}'
    )
```

**In code:**
```python
# __init__.py
from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("intronIC")
except PackageNotFoundError:
    # Package not installed
    __version__ = "unknown (not installed)"
```

## Version Consistency Check

Create `scripts/check_version.py`:
```python
#!/usr/bin/env python3
"""Check that version is consistent across files."""

import sys
import tomli

# Read pyproject.toml
with open("pyproject.toml", "rb") as f:
    pyproject = tomli.load(f)
    pyproject_version = pyproject["project"]["version"]

# Read pixi.toml
with open("pixi.toml", "rb") as f:
    pixi = tomli.load(f)
    pixi_version = pixi["workspace"]["version"]

# Check consistency
if pyproject_version != pixi_version:
    print(f"ERROR: Version mismatch!")
    print(f"  pyproject.toml: {pyproject_version}")
    print(f"  pixi.toml: {pixi_version}")
    sys.exit(1)

print(f"✓ Version consistent: {pyproject_version}")
```

## Implementation Steps

1. **Decide on version number**: `2.0.0` for first refactored release
2. **Update pyproject.toml**: Set version to `2.0.0`
3. **Update pixi.toml**: Set version to `2.0.0`
4. **Update __init__.py**: Read version from package metadata
5. **Add --version to CLI**: Show version on command line
6. **Add version check script**: Ensure consistency
7. **Document in RELEASE_PROCESS.md**: How to make releases
8. **Test**: Verify `intronIC --version` works

## Future Migration to setuptools-scm

If automatic versioning is desired later:

```bash
# 1. Install setuptools-scm
pip install setuptools-scm

# 2. Update pyproject.toml
# (change version to dynamic)

# 3. Tag current state
git tag v2.0.0

# 4. Test
pip install -e .
python -c "import intronIC; print(intronIC.__version__)"
```

## Summary

**Current Recommendation: Static versioning in pyproject.toml**

- Simple, explicit, works everywhere
- Manual but requires discipline
- Easy to migrate to setuptools-scm later if desired
- Version number: Start with `2.0.0` for refactored release

This is **PyPI-ready** - static versions work perfectly fine for publication.
