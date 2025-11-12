# intronIC Installation Guide

## Quick Start

### For End Users (Recommended)

The simplest way to install intronIC:

```bash
# Install with pip
pip install git+https://github.com/glarue/intronIC.git#subdirectory=intronIC_refactored

# Or clone and install locally
git clone https://github.com/glarue/intronIC.git
cd intronIC/intronIC_refactored
pip install .

# Test installation
intronIC --help
intronIC --version
```

### For Developers

For reproducible development environment:

```bash
# Install pixi (one-time setup)
curl -fsSL https://pixi.sh/install.sh | bash

# Clone and setup
git clone https://github.com/glarue/intronIC.git
cd intronIC/intronIC_refactored

# Create environment and install dependencies
pixi install

# Install intronIC in development mode
pixi run install

# Now you can use either:
intronIC --help           # Direct command (after install)
pixi run intronIC --help  # Through pixi
python -m intronIC --help # As module
```

---

## Detailed Installation Options

### Option 1: Pip Install (Simple, Universal)

**Requirements:**
- Python 3.8 - 3.12
- pip

**Steps:**

1. **Create a virtual environment** (recommended):
   ```bash
   python -m venv intronic-env
   source intronic-env/bin/activate  # On Windows: intronic-env\Scripts\activate
   ```

2. **Install intronIC**:
   ```bash
   pip install git+https://github.com/glarue/intronIC.git#subdirectory=intronIC_refactored
   ```

3. **Verify installation**:
   ```bash
   intronIC --version
   intronIC --help
   ```

4. **Run intronIC**:
   ```bash
   intronIC -g genome.fa -a annotation.gff3 -n species_name
   ```

**Pros:**
- Simple, one-command installation
- Works in any Python environment (venv, conda, system)
- Familiar to most Python users
- Lightweight

**Cons:**
- Need to manage Python version yourself
- Dependency resolution might differ across platforms

---

### Option 2: Pixi Install (Reproducible, Cross-Platform)

**Requirements:**
- [Pixi](https://pixi.sh/) package manager

**Steps:**

1. **Install pixi** (one-time setup):
   ```bash
   # Linux/macOS
   curl -fsSL https://pixi.sh/install.sh | bash

   # Windows PowerShell
   iwr -useb https://pixi.sh/install.ps1 | iex
   ```

2. **Clone repository**:
   ```bash
   git clone https://github.com/glarue/intronIC.git
   cd intronIC/intronIC_refactored
   ```

3. **Install dependencies and package**:
   ```bash
   pixi install          # Creates locked environment
   pixi run install      # Installs intronIC in dev mode
   ```

4. **Run intronIC**:
   ```bash
   # Option A: Direct command (if installed in development mode)
   intronIC -g genome.fa -a annotation.gff3 -n species_name

   # Option B: Through pixi
   pixi run intronIC -g genome.fa -a annotation.gff3 -n species_name

   # Option C: As Python module
   pixi run python -m intronIC -g genome.fa -a annotation.gff3 -n species_name
   ```

**Pros:**
- Guaranteed reproducibility (locked dependencies)
- Cross-platform binary handling
- Environment management built-in
- Great for development
- Includes dev tools (pytest, black, mypy)

**Cons:**
- Requires installing pixi
- Slightly more complex setup
- Larger disk footprint

---

### Option 3: Local Development Install

For contributing to intronIC development:

```bash
# Clone repository
git clone https://github.com/glarue/intronIC.git
cd intronIC/intronIC_refactored

# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install in editable mode with development dependencies
pip install -e ".[dev]"

# Run tests
pytest

# Format code
black .

# Type check
mypy cli core extraction scoring classification file_io utils
```

**Development dependencies include:**
- pytest (testing)
- pytest-cov (coverage)
- black (formatting)
- ruff (linting)
- mypy (type checking)
- ipython (interactive shell)

---

## Usage Examples

### Basic Classification

```bash
intronIC \
  -g genome.fa.gz \
  -a annotation.gff3.gz \
  -n homo_sapiens \
  -p 8
```

### Extract Sequences Only

```bash
intronIC \
  -g genome.fa.gz \
  -a annotation.gff3.gz \
  -n species_name \
  --sequences_only
```

### With Pretrained Model

```bash
intronIC \
  -g genome.fa.gz \
  -a annotation.gff3.gz \
  -n species_name \
  --pretrained_model homo_sapiens.model.pkl
```

### Recursive Training

```bash
intronIC \
  -g genome.fa.gz \
  -a annotation.gff3.gz \
  -n distant_species \
  --recursive
```

---

## Verifying Installation

After installation, verify it works:

```bash
# Check version
intronIC --version

# View help
intronIC --help

# Run quick test (if you have test data)
intronIC \
  -g test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
  -a test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz \
  -n test_run \
  -s
```

Expected output:
- No errors
- Creates output files: `test_run.seqs.iic`, `test_run.meta.iic`, etc.
- Completes in 1-2 minutes

---

## Troubleshooting

### ImportError: No module named 'intronIC'

**Problem**: Package not installed properly.

**Solution**:
```bash
pip install .  # or pip install -e . for development
```

### Command 'intronIC' not found

**Problem**: Entry point not in PATH or not installed.

**Solutions**:
1. Ensure virtual environment is activated
2. Use full path: `python -m intronIC`
3. Reinstall: `pip install --force-reinstall .`

### ModuleNotFoundError: networkx/scikit-learn/etc.

**Problem**: Dependencies not installed.

**Solution**:
```bash
pip install --upgrade pip
pip install .  # Reinstall with all dependencies
```

### Tests fail after installation

**Problem**: Test dependencies not installed.

**Solution**:
```bash
pip install ".[test]"  # Install with test dependencies
pytest
```

### Multiprocessing errors on macOS

**Problem**: macOS default multiprocessing method incompatibility.

**Solution**: intronIC automatically handles this. If issues persist:
```bash
export OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES
intronic [args]
```

---

## Uninstalling

```bash
pip uninstall intronIC
```

For pixi:
```bash
# Remove environment
pixi clean

# Or remove entire .pixi directory
rm -rf .pixi
```

---

## System Requirements

### Minimum
- Python 3.8+
- 4 GB RAM
- 2 CPU cores

### Recommended
- Python 3.10+
- 8+ GB RAM
- 8+ CPU cores (for parallel processing)
- 10 GB disk space (for large genomes)

### Dependencies
Automatically installed via pip/pixi:
- numpy < 2.0
- scipy >= 1.5.0
- scikit-learn >= 0.22, < 2.0
- biogl >= 0.1.0
- matplotlib >= 3.3.0
- networkx >= 2.5.1
- rich >= 10.0

---

## Getting Help

- **Documentation**: See [README.md](README.md)
- **Issues**: https://github.com/glarue/intronIC/issues
- **Citation**: Moyer et al. (2020) NAR 48(13):7066–7078

---

## Platform-Specific Notes

### Linux
- Works out of the box
- Tested on Ubuntu 20.04+, CentOS 7+

### macOS
- Works with both Intel and Apple Silicon
- May need Xcode Command Line Tools: `xcode-select --install`

### Windows
- Works with Python 3.8+
- Use Git Bash or PowerShell
- Some paths may need Windows-style separators (`\` instead of `/`)

---

## Alternative: Docker (Future)

For complete environment isolation (coming soon):

```bash
docker pull glarue/intronic:latest
docker run -v $(pwd):/data glarue/intronic \
  -g /data/genome.fa \
  -a /data/annotation.gff3 \
  -n species_name
```
