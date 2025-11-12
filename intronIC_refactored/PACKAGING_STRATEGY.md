# intronIC Packaging Strategy

## Current Situation Analysis

### What We Have
1. **pyproject.toml**: Modern Python packaging with hatchling backend
   - Defines package metadata, dependencies, and entry point
   - Entry point: `intronIC = "cli.main:main"`
   - Ready for pip installation

2. **pixi.toml**: Conda-based environment management
   - Handles both conda and PyPI dependencies
   - Includes development tasks
   - Cross-platform support

3. **__main__.py**: Module execution support
   - Allows `python -m intronIC_refactored`
   - Handles multiprocessing setup

### What Users Want
- Simple installation: `pip install intronIC`
- Simple execution: `intronIC [args]` or `python -m intronIC [args]`
- Easy to get help: `intronIC --help`
- Optional: Development installation for contributors

## Recommended Approach: **Hybrid Strategy**

Support BOTH pip and pixi to maximize accessibility:

### For End Users (Simple Pip Installation)
```bash
# Install from PyPI (when published)
pip install intronIC

# Or install from GitHub
pip install git+https://github.com/glarue/intronIC.git

# Or install locally
git clone https://github.com/glarue/intronIC.git
cd intronIC
pip install .

# Run
intronIC -g genome.fa -a annotation.gff3 -n species_name
```

### For Developers/Contributors (Pixi for Reproducibility)
```bash
# Clone repository
git clone https://github.com/glarue/intronIC.git
cd intronIC/intronIC_refactored

# Install pixi (if not already installed)
curl -fsSL https://pixi.sh/install.sh | bash

# Create environment and install in dev mode
pixi install
pixi run install  # Runs: pip install -e .

# Run
pixi run intronIC -g genome.fa -a annotation.gff3 -n species_name
# Or after install: intronIC [args]
```

## Why Hybrid is Best

### Advantages of Keeping Pip Support:
1. **Accessibility**: Most bioinformaticians know pip, not all know conda/pixi
2. **Simplicity**: One command installation
3. **PyPI Distribution**: Can publish to PyPI for `pip install intronIC`
4. **CI/CD**: Works with standard Python testing infrastructure
5. **Compatibility**: Works in any Python environment (venv, conda, system)
6. **Lightweight**: No need to install conda/pixi for basic usage

### Advantages of Keeping Pixi:
1. **Reproducibility**: Lock file ensures exact dependency versions
2. **Cross-platform**: Handles platform-specific binaries
3. **Speed**: Conda solver often faster than pip for scientific packages
4. **Development**: Great for contributor environment setup
5. **Non-Python Dependencies**: Can handle system libraries if needed
6. **Environment Isolation**: Built-in virtual environment management

### Why Not Choose One?
- **Pip-only**: Loses reproducibility guarantees, harder cross-platform support
- **Pixi-only**: Alienates users unfamiliar with conda ecosystem, can't publish to PyPI

## Implementation Plan

### 1. Fix Package Structure
Current issue: Modules are in root directory, not in a package directory.

**Current structure:**
```
intronIC_refactored/
├── cli/
├── core/
├── extraction/
├── scoring/
├── classification/
├── file_io/
├── utils/
├── intronIC/data/  # Data files
└── __main__.py
```

**Should be:**
```
intronIC_refactored/
├── intronIC/  # Main package directory
│   ├── __init__.py
│   ├── __main__.py
│   ├── cli/
│   ├── core/
│   ├── extraction/
│   ├── scoring/
│   ├── classification/
│   ├── file_io/
│   ├── utils/
│   └── data/  # Data files included in package
├── tests/
├── pyproject.toml
├── pixi.toml
└── README.md
```

**Why?**
- Standard Python packaging expects a top-level package directory
- Makes imports cleaner: `from intronIC.core import Intron`
- Prevents namespace collisions
- Matches original intronIC structure

### 2. Update pyproject.toml
```toml
[project.scripts]
intronIC = "intronIC.cli.main:main"

[tool.hatch.build.targets.wheel]
packages = ["intronIC"]

[tool.hatch.build.targets.wheel.force-include]
"intronIC/data" = "intronIC/data"
```

### 3. Update pixi.toml
```toml
[tasks]
# Install in development mode
install = "pip install -e ."

# Run intronIC (works after install)
intronIC = "intronIC"

# Or run directly without install
run = "python -m intronIC"
```

### 4. Create Simple Installation Scripts

**install.sh** (for Unix):
```bash
#!/bin/bash
set -e
echo "Installing intronIC..."
pip install .
echo "Installation complete! Run 'intronIC --help' to get started."
```

**install.bat** (for Windows):
```batch
@echo off
echo Installing intronIC...
pip install .
echo Installation complete! Run 'intronIC --help' to get started.
```

### 5. Update README with Clear Installation Instructions

See separate README update below.

## Migration Steps

1. **Restructure directories** (use git mv to preserve history)
2. **Update pyproject.toml** with correct package paths
3. **Update pixi.toml** with correct paths
4. **Update all imports** in source files
5. **Update test imports**
6. **Test pip installation**: `pip install .`
7. **Test command works**: `intronIC --help`
8. **Test pixi workflow**
9. **Update documentation**

## Testing the Package

### Test Pip Installation
```bash
# Create clean venv
python -m venv test_venv
source test_venv/bin/activate

# Install
pip install .

# Test
intronIC --help
intronIC --version

# Test with real data
intronIC -g test_data/genome.fa -a test_data/annotation.gff3 -n test

# Clean up
deactivate
rm -rf test_venv
```

### Test Pixi Installation
```bash
# Clean environment
pixi clean
rm -rf .pixi

# Reinstall
pixi install

# Test
pixi run intronIC --help
pixi run test-small  # Run built-in test
```

## Future: PyPI Publication

Once stable, publish to PyPI:

```bash
# Build distributions
pip install build
python -m build

# Upload to TestPyPI first
pip install twine
twine upload --repository testpypi dist/*

# Test installation from TestPyPI
pip install --index-url https://test.pypi.org/simple/ intronIC

# If all good, upload to real PyPI
twine upload dist/*
```

Then users can simply:
```bash
pip install intronIC
```

## Recommendation

**Implement the hybrid approach** with these priorities:

1. **High Priority**: Fix package structure and make pip installation work smoothly
2. **High Priority**: Update README with clear installation instructions for both methods
3. **Medium Priority**: Keep pixi.toml updated for development workflow
4. **Low Priority**: Publish to PyPI (after thorough testing)

This gives users maximum flexibility while maintaining reproducibility for development.
