![intronIC_logo](https://user-images.githubusercontent.com/6827531/82829967-62872480-9e69-11ea-94e9-fa7306c7df1b.png)

# intronIC - (intron <ins>I</ins>nterrogator and <ins>C</ins>lassifier)

**Version 2.0.0** - Refactored Edition with Corrected Architecture

`intronIC` is a bioinformatics tool for extracting and classifying intron sequences as **U12-type (minor)** or **U2-type (major)** using a support vector machine (SVM) trained on position-weight matrix (PWM) scores. It can be used with a genome and annotation file, or with pre-extracted intron sequences. Alternatively, `intronIC` can extract all annotated intron sequences without classification (using the `extract` subcommand).

---

## About This Refactored Version

This refactored version maintains **100% algorithmic fidelity** and **CLI compatibility** with the original intronIC while providing a modernized, maintainable codebase:

### Key Improvements

- **Corrected ML Architecture (v2.0)**: Fixed double-scaling issue and train/test mismatch
  - Single scaling step via RobustScaler with centering (removes composition bias)
  - Configurable augmented features with 5D standard (`absdiff_bp_3`, `absdiff_5_bp`) or custom feature sets
  - Two-stage optimization (C via balanced_accuracy, calibration via log-loss)
  - L1/L2 penalty search with class weight multiplier optimization
- **Modular Architecture**: Organized into logical packages (extraction, scoring, classification, output) instead of a single 6,000+-line file
- **Enhanced Code Quality**: Type hints throughout, immutable data structures, better error handling
- **Bug Fixes**: Corrected data leakage in z-score normalization, fixed type_id assignment logic
- **Better Testing**: Structured for unit and integration testing
- **Modern Tooling**: Support for `pixi` and `uv` package managers
- **Enhanced Logging**: Clearer progress reporting with section markers and detailed training logs
- **Improved Documentation**: Comprehensive inline documentation and external guides

### What's Preserved

- **Same Classification Algorithm**: Linear SVM with balanced class weights
- **Same Feature Extraction**: PWM scoring of 5' splice site, branch point, and 3' splice site
- **Same Output Formats**: All `.iic` files maintain compatibility (with minor enhancements)
- **Same Performance**: Comparable runtime and memory usage to original
- **Validated Accuracy**: Identical classification results on test data

For detailed information about the refactoring approach, see [REFACTOR_SUMMARY.md](REFACTOR_SUMMARY.md).

---

## Scientific Background

### Minor (U12-type) vs Major (U2-type) Introns

Most eukaryotic introns (~99.5%) are spliced by the **major (U2-type) spliceosome** and typically have:
- 5' splice site: **GT**
- 3' splice site: **AG**
- Branch point: **A** within a loose consensus

A small fraction (~0.5%) are spliced by the **minor (U12-type) spliceosome** and typically have:
- 5' splice site: **AT** (AT-AC type) or **GT** (GT-AG type)
- 3' splice site: **AC** (AT-AC type) or **AG** (GT-AG type)
- Branch point: Highly conserved **TCCTTAAC** motif

### Classification Approach

intronIC uses a **three-step scoring and classification pipeline**:

1. **PWM Scoring**: Apply position-weight matrices to three key regions (5' splice site, branch point, 3' splice site) to calculate raw log-odds scores
2. **Normalization**: Convert raw scores to z-scores using parameters fit on reference sequences only (prevents data leakage)
3. **SVM Classification**: Train an ensemble of linear SVMs on reference U12/U2 introns, output probability scores (0-100%)

The output probability represents the classifier's confidence that an intron is U12-type. By default, introns with scores **>90%** are considered high-confidence U12-type predictions.

### ML Pipeline Architecture

intronIC uses a **single scaling step** architecture to prevent double-scaling and ensure train/test consistency:

```
Raw PWM Scores (LLRs)
         ↓
ScoreNormalizer (EXTERNAL to pipeline)
  - RobustScaler(with_centering=True)
  - Fitted on reference introns only
  - Transforms: raw LLRs → z-scores
  - Removes composition bias via centering
         ↓
Z-Scores [five_z_score, bp_z_score, three_z_score]
         ↓
ML Pipeline (NO scaler inside)
  ├─ BothEndsStrongTransformer
  │  └─ Augments 3D → 5D features (standard config):
  │     • Pass-through: five_z, bp_z, three_z
  │     • absdiff_bp_3 = |bp_z - three_z| (BP/3' imbalance penalty)
  │     • absdiff_5_bp = |five_z - bp_z| (5'/BP imbalance penalty)
  │  └─ Or custom 4D-7D with different features:
  │     • min_all, absdiff_5_3, min_5_bp, max_5_bp, etc.
  ├─ LinearSVC
  │  └─ L1 or L2 penalty (grid-searched), balanced class weights
  └─ CalibratedClassifierCV
     └─ External calibration (sigmoid or isotonic)
         ↓
U12 Probability (0-100%)
```

**Key Design Principles:**

1. **Single Scaling Step**: Scaling happens ONLY in ScoreNormalizer (external to pipeline). The pipeline receives pre-scaled z-scores and does NOT re-scale them. This prevents double-scaling.

2. **Train/Test Consistency**: Both training and prediction extract z-scores from introns and pass them to the pipeline, ensuring identical data transformations.

3. **Domain Adaptation**: ScoreNormalizer can be refitted per-species (adaptive mode) or reused from training species (human mode) for cross-species classification.

4. **Feature Engineering**: BothEndsStrongTransformer adds configurable composite features. The standard 5D configuration adds `absdiff_bp_3` and `absdiff_5_bp` (BP/3' and 5'/BP imbalance penalties) based on L1 regularization analysis. See `config/config.yaml` for all available features.

5. **Hyperparameter Optimization**:
   - **Grid search over**: C parameter, L1/L2 penalty, class weight multipliers
   - **Stage 1**: Optimize C using balanced_accuracy (discrimination quality)
   - **Stage 2**: Select calibration method (sigmoid vs isotonic) using log-loss (probability quality)

6. **YAML Configuration**: All optimizer settings are configurable via `config/config.yaml` including feature selection, penalty options, class weight multipliers, and CV parameters.

This architecture was validated on C. elegans, achieving **1 false positives** (1/109,830) vs 130 with uncentered scaling.

---

## Installation

### Method 1: Using `pixi` (Recommended)

[Pixi](https://pixi.sh/) is a modern package manager that handles all dependencies automatically:

```bash
# Install pixi (if not already installed)
curl -fsSL https://pixi.sh/install.sh | bash

# Clone the repository
git clone https://github.com/glarue/intronIC.git
cd intronIC/intronIC_refactored

# Install dependencies and set up environment
pixi install

# Run on test data
pixi run test-small
```

### Method 2: Using `pip`

Traditional Python installation:

```bash
# Clone the repository
git clone https://github.com/glarue/intronIC.git
cd intronIC/intronIC_refactored

# Install in development mode
pip install -e .

# intronIC is now callable from command line
intronIC -h
```

### Method 3: Using `uv` (Fast Alternative)

[uv](https://github.com/astral-sh/uv) is a fast Python package installer:

```bash
# Install uv (if not already installed)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clone and set up
git clone https://github.com/glarue/intronIC.git
cd intronIC/intronIC_refactored

uv venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
uv pip install -e .

# Run intronIC
intronIC -h
```

---

## Dependencies

`intronIC` requires Python 3.8-3.12 and the following packages:

* **[numpy](https://numpy.org/)** `>=1.19.0, <2.0` - Numerical operations
* **[scipy](https://scipy.org/)** `>=1.5.0` - Scientific computing
* **[scikit-learn](https://scikit-learn.org/)** `>=0.22, <2.0` - SVM classifier
* **[matplotlib](https://matplotlib.org/)** `>=3.3.0` - Plotting
* **[networkx](https://networkx.org/)** `>=2.5.1` - Graph operations for annotation parsing
* **[rich](https://rich.readthedocs.io/)** `>=10.0` - Terminal progress bars
* **[biogl](https://github.com/glarue/biogl)** `>=0.1.0` - Bioinformatics utilities

All dependencies are automatically installed by `pixi`, `uv`, or `pip`.

`intronIC` was developed on Linux and has been tested on macOS and Windows.

---

## Quick Start

### Installation (30 seconds)

```bash
# Clone and install
git clone https://github.com/glarue/intronIC.git
cd intronIC
pip install -e .

# Verify installation
intronIC --version
```

### Basic Commands

```bash
# Classify introns (train on-the-fly)
intronIC -g genome.fa.gz -a annotation.gff3.gz -n species_name -p 8

# Use pretrained model (faster)
intronIC -g genome.fa.gz -a annotation.gff3.gz -n species_name --model model.pkl -p 8

# Extract sequences only (no classification)
intronIC extract -g genome.fa.gz -a annotation.gff3.gz -n species_name -p 8

# Train a model (no genome needed)
intronIC train -n my_model -p 8
```

### Test Run (Human Chr19)

```bash
# With test data included in repository
intronIC -g test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
         -a test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz \
         -n homo_sapiens_chr19 -p 4
```

Expected results:
- **~29,000** introns extracted
- **~30** U12-type introns (score ≥90%)
- **~8** AT-AC type U12 introns
- Output files: `homo_sapiens_chr19.*.iic`

---

## Usage

### Commands

intronIC supports three subcommands:

| Command | Description |
|---------|-------------|
| (default) | Classify introns from genome + annotation |
| `train` | Train a model on reference data (no genome needed) |
| `extract` | Extract sequences only (no classification) |

### Default Mode: Classify Introns

The default mode extracts introns and classifies them as U12 or U2 type:

```bash
# Basic usage (trains model on-the-fly)
intronIC -g genome.fa.gz -a annotation.gff3.gz -n species_name

# With pretrained model (faster, recommended)
intronIC -g genome.fa.gz -a annotation.gff3.gz -n species_name \
  --model homo_sapiens.model.pkl -p 8

# Memory-efficient streaming mode
intronIC -g genome.fa.gz -a annotation.gff3.gz -n species_name \
  --model homo_sapiens.model.pkl --streaming -p 8
```

### Train Subcommand

Train a classifier model without needing a genome:

```bash
# Basic training with built-in references
intronIC train -n my_model -p 8

# With custom configuration
intronIC --config config/config.yaml train -n my_model -p 12

# With custom reference sequences
intronIC train -n my_model -p 8 \
  --reference_u12s custom_u12.iic \
  --reference_u2s custom_u2.iic

# Quick training (skip nested CV evaluation)
intronIC train -n my_model --eval_mode none -p 8
```

Output: `my_model.model.pkl` - use with `--model` for classification.

### Extract Subcommand

Extract intron sequences without classification:

```bash
# Extract from annotation (streaming mode by default)
intronIC extract -g genome.fa.gz -a annotation.gff3.gz -n species_name

# Extract from BED file
intronIC extract -g genome.fa.gz -b introns.bed -n species_name

# With custom flank length
intronIC extract -g genome.fa.gz -a annotation.gff3.gz -n species_name --flank-len 20
```

Output: `.introns.iic`, `.meta.iic`, `.bed.iic` files (no classification scores).

### Required Arguments

| Argument | Short | Description |
|----------|-------|-------------|
| `--genome` | `-g` | Genome FASTA file (gzip supported) |
| `--annotation` | `-a` | GFF3/GTF annotation file (gzip supported) |
| `--species-name` | `-n` | Species name / output prefix |

Alternative inputs (instead of `-a`):
- `-b FILE` - BED file with intron coordinates
- `-q FILE` - Pre-extracted sequences file

### Common Options

| Argument | Short | Description | Default |
|----------|-------|-------------|---------|
| `--processes` | `-p` | Number of CPU cores | 1 |
| `--threshold` | `-t` | U12 probability threshold (0-100) | 90 |
| `--model` | | Pretrained model file | None |
| `--config` | | YAML configuration file | Auto-discovered |
| `--streaming` | | Memory-efficient mode | False |
| `--sequences-only` | `-s` | Skip classification | False |
| `--feature-type` | `-f` | Feature type: `cds`, `exon`, or `both` | both |
| `--allow-multiple-isoforms` | `-i` | Include all isoforms | False (longest only) |
| `--exclude-overlapping` | `-v` | Exclude overlapping introns | False |
| `--no-nc` | | Exclude non-canonical introns | False |
| `--recursive` | | Recursive training | False |

### Usage Examples

**1. Basic classification:**
```bash
intronIC -g genome.fa.gz -a annotation.gff3.gz -n my_species -p 8
```

**2. With pretrained model (recommended for speed):**
```bash
intronIC -g genome.fa.gz -a annotation.gff3.gz -n my_species \
  --model homo_sapiens.model.pkl -p 8
```

**3. Streaming mode for large genomes:**
```bash
intronIC -g genome.fa.gz -a annotation.gff3.gz -n my_species \
  --model homo_sapiens.model.pkl --streaming -p 8
```

**4. Extract sequences only:**
```bash
intronIC extract -g genome.fa.gz -a annotation.gff3.gz -n my_species -p 8
```

**5. Train custom model:**
```bash
intronIC train -n my_trained_model -p 8 --config config/config.yaml
```

**6. Stricter threshold (95%):**
```bash
intronIC -g genome.fa.gz -a annotation.gff3.gz -n my_species -t 95 -p 8
```

**7. Include all isoforms, CDS only:**
```bash
intronIC -g genome.fa.gz -a annotation.gff3.gz -n my_species -i -f cds -p 8
```

**8. Classify from BED coordinates:**
```bash
intronIC -g genome.fa.gz -b intron_coordinates.bed -n my_species \
  --model homo_sapiens.model.pkl
```

---

## Output Files

All output files are tab-delimited with the `.iic` extension and named `{species_name}.{type}.iic`.

### Main Output Files

**1. `.meta.iic` - Comprehensive Metadata**

Contains detailed information for each intron:
- Intron name/label with tags
- Relative score (distance from threshold)
- Terminal dinucleotides (e.g., GT-AG, AT-AC)
- Motif schematic (showing branch point context)
- Branch point region sequence
- Intron length
- Parent transcript ID
- Grandparent gene ID
- Intron index and total family size
- Fractional position in transcript
- Exon phase
- Type ID (u2 or u12)
- Attributes (longest_isoform, corrected, etc.)

**2. `.bed.iic` - BED-format Coordinates**

Standard BED format with scores:
- Chromosome
- Start (0-based, BED standard)
- Stop (1-based)
- Label (`intron_name;probability`)
- SVM score (0-100, integer)
- Strand

**3. `.seqs.iic` - Sequences**

Intron sequences with flanking regions:
- Intron name
- 5' flanking sequence (exonic)
- Intron sequence
- 3' flanking sequence (exonic)
- SVM score (if classification performed)

**4. `.scores.iic` (or `.score_info.iic`) - Detailed Scoring**

Per-intron breakdown of all scores:
- Name and scores (relative, SVM, decision_distance)
- 5' splice site: sequence, raw score, z-score
- Branch point: sequences (U12 and U2 versions), raw score, z-score
- 3' splice site: sequence, raw score, z-score

**5. Mapping Files**

- `.dupe_map.iic` - Maps duplicate introns to their representative
- `.overlap_map.iic` - Maps overlapping intron coordinates

**6. Visualization Files (`.png`)**

- `*_scatter.png` - 2D scatter plot of classified introns with marginal distributions
- `*_training_scatter.png` - Scatter plot of training data
- `*_training_hexplot.png` - Hexbin density plot of reference introns
- `*_pr_curve.png` - Precision-Recall AUC curves for model evaluation

**7. Log Files**

- `.log` - Main log file with pipeline progress and summary statistics
- `.training.log` - Detailed training log (when models are trained, not with `--model`)

### Identifying U12-type Introns

U12-type introns are identified by their **relative score > 0** (equivalent to SVM score > threshold):

```bash
# Extract U12-type introns from meta file
awk '($2!="." && $2>0)' species_name.meta.iic

# Count U12-type introns
awk '($2!="." && $2>0)' species_name.meta.iic | wc -l

# Get U12-type intron names
awk '($2!="." && $2>0) {print $1}' species_name.meta.iic

# Filter by higher confidence (relative score > 10)
awk '($2!="." && $2>10)' species_name.meta.iic
```

### Understanding the Scores

**SVM Score (0-100):**
- Probability that the intron is U12-type
- 50 = equal probability of U2 or U12
- >90 = high confidence U12 (default threshold)
- <10 = high confidence U2

**Relative Score:**
- Distance from the threshold
- Calculated as: `svm_score - threshold`
- Positive values = above threshold (U12-type at chosen confidence)
- Negative values = below threshold (U2-type)
- Makes filtering easier: just check if > 0

**Type ID (u2 or u12):**
- Binary classification based on raw classifier decision (50% boundary)
- Independent of the user-chosen threshold
- Used for organizing output and statistics

**Decision Distance:**
- Log-odds ratio: `log(probability / (1 - probability))`
- 0 = equal probability (50%)
- Positive = favors U12
- Negative = favors U2
- Useful for understanding classifier confidence

---

## A Note on the `-n` (Name) Argument

By default, `intronIC` expects species names in **binomial format** (genus, species) separated by a non-alphanumeric character:

* `homo_sapiens` ✅
* `homo.sapiens` ✅
* `homo-sapiens` ✅

`intronIC` formats the name internally into a tag for intron IDs (e.g., `HomSap`), using only the first two elements.

**Output files** are named using the full argument supplied to `-n`, so:
- `homo_sapiens` → files named `homo_sapiens.*`
- `homo_sapiens.v2` → files named `homo_sapiens.v2.*`
- Intron IDs in both would use `HomSap` tag

To use the name argument exactly as provided without any parsing, add the `--na` flag:
```bash
intronIC -g genome.fa.gz -a annotation.gff3.gz -n "My Custom Name" --na
```

---

## Resource Usage

### Memory

Memory usage scales with the number of annotated introns in the genome:

* **Small genomes** (<50,000 introns): <1 GB
* **Typical genomes** (50,000-200,000 introns): 1-3 GB
* **Large genomes** (>200,000 introns): 3-5 GB
* **Human genome** (Ensembl 95, ~1 million introns): ~5 GB
* **Streaming mode** (with `--model --streaming`): ~0.5 GB regardless of genome size

Most modern computers should handle even large genomes without issue. For memory-constrained environments, use streaming mode with a pretrained model.

### Runtime

Runtime depends on genome size, annotation density, and whether models are pre-trained:

| Genome | Introns | Train Mode | Pretrained (`--model`) |
|--------|---------|------------|------------------------|
| Chr19 (test) | ~29,000 | 5-15 min | 1-3 min |
| Small genome | ~50,000 | 10-30 min | 2-5 min |
| Human (full) | ~1,000,000 | 20-40 min | 5-10 min |

**Tips for faster runs:**
- Use `--model` with a pretrained model to skip training (fastest)
- Use `-p N` for parallel processing (recommended: 4-8 cores)
- Use `--streaming` with `--model` for large genomes with memory constraints
- Use small reference sets for testing (`--reference_u12s`, `--reference_u2s`)
- Extract sequences first with `extract` subcommand, then classify separately if iterating on parameters

---

## Advanced Usage

### Using Pretrained Models

For cross-species classification using a model trained on another species:

```bash
# Use a specific trained model file
intronIC -g genome.fa.gz -a annotation.gff3.gz -n species_name \
  --model /path/to/trained_species.model.pkl
```

This is the recommended approach for:
- Classifying species without curated U12 references
- Applying a human-trained model to other vertebrates
- Fast classification when training data is unavailable

The pretrained model contains:
- Trained SVM ensemble with optimized hyperparameters
- Frozen scaler from training species (for cross-species normalization)
- Model metadata (training parameters, feature configuration)

### Streaming Mode

For large genomes with memory constraints, streaming mode processes introns per-chromosome:

```bash
# Memory-efficient streaming with pretrained model
intronIC -g genome.fa.gz -a annotation.gff3.gz -n species_name \
  --model trained.model.pkl --streaming -p 8
```

Streaming mode provides ~90% memory savings by:
- Processing one chromosome at a time
- Writing results immediately (not accumulating in memory)
- Using the frozen scaler from the pretrained model

**Requirements**: Streaming mode requires `--model` (pretrained model with frozen scaler).

### Configuration Files

intronIC uses YAML configuration files for advanced parameter tuning:

```bash
# Use custom configuration
intronIC -g genome.fa.gz -a annotation.gff3.gz -n species_name \
  --config config/profiles/production.yaml
```

Configuration files are auto-discovered from (in priority order):
1. `--config PATH` (explicit CLI argument)
2. `./.intronIC.yaml` (current directory)
3. `~/.config/intronIC/config.yaml` (XDG config)
4. Built-in defaults

Key configurable parameters include:
- **Feature selection**: Choose which augmented features to use (5D standard or custom)
- **Penalty options**: L1, L2, or both for regularization search
- **Class weight multipliers**: Fine-tune precision/recall tradeoff
- **CV parameters**: Number of folds, optimization rounds
- **Ensemble settings**: Number of models, subsampling ratio

See `config/config.yaml` for full documentation of all options.

### Recursive Training

For species distant from the training data, recursive training can improve accuracy:

```bash
intronIC -g genome.fa.gz -a annotation.gff3.gz -n species_name --recursive
```

This performs two passes:
1. Initial classification to identify high-confidence U12-type introns
2. Build species-specific PWMs and retrain models
3. Re-classify all introns with the updated models

### Custom Reference Sequences

For specialized analyses, you can provide custom reference sequences:

```bash
intronIC -g genome.fa.gz -a annotation.gff3.gz -n species_name \
  --reference_u12s my_u12_introns.iic \
  --reference_u2s my_u2_introns.iic
```

Reference files should follow the `.iic` format (tab-delimited: name, 5'_flank, intron_seq, 3'_flank).

### Two-Stage Workflow

For large genomes or parameter tuning, you can separate extraction from classification:

**Stage 1: Extract sequences only**
```bash
intronIC extract -g genome.fa.gz -a annotation.gff3.gz -n species_name
# Produces: species_name.introns.iic (and .meta.iic, .bed.iic)
```

**Stage 2: Classify extracted sequences**
```bash
intronIC -q species_name.introns.iic -n species_name -t 95
# Much faster for testing different thresholds or references
```

---

## Troubleshooting

### Common Issues

**"No U12-type introns found"**
- Normal for some small genomes or chromosomes
- Try lowering threshold: `-t 80`
- Check that annotation contains sufficient introns
- Consider using `--recursive` for distant species

**"Out of memory" errors**
- Use a machine with more RAM for very large genomes
- Try processing chromosomes separately using BED input
- Reduce parallelization: `-p 1` or `-p 2`

**"No introns extracted"**
- Check that genome and annotation use matching chromosome names
- Verify annotation format (GFF3 or GTF)
- Try different feature type: `-f cds` or `-f exon`
- Check annotation file is not corrupted

**Slow performance**
- Use parallel processing: `-p 4` or `-p 8`
- Use `--model` to skip model training
- Use smaller reference sets for testing
- Consider extracting sequences first (`extract` subcommand), then classify separately

**Classification results differ from original intronIC**
- Minor differences can occur due to:
  - Random seed in cross-validation
  - sklearn version differences
  - Floating-point precision
- Major differences are unexpected; please file an issue

### Getting Help

* **Documentation**: See the [original wiki](https://github.com/glarue/intronIC/wiki) for detailed guides
* **Issues**: Report bugs at [GitHub Issues](https://github.com/glarue/intronIC/issues)
* **Questions**: Open a discussion at [GitHub Discussions](https://github.com/glarue/intronIC/discussions)

For refactoring-specific questions, see [REFACTOR_SUMMARY.md](REFACTOR_SUMMARY.md).

---

## Testing the Installation

To verify your installation works correctly:

```bash
# With pixi
cd intronIC_refactored
pixi run test-small  # Should complete in 2-5 minutes

# With pip
intronIC -g ../intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
         -a ../intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz \
         -n test_run
```

Expected output:
- Several `.iic` files named `test_run.*` or `homo_sapiens_test_small.*`
- A `.log` file with classification summary
- PNG plots showing score distributions
- Console output showing ~30 U12-type introns found

---

## Project Structure

The refactored codebase is organized into logical modules:

```
intronIC_refactored/
├── cli/                 # Command-line interface and orchestration
│   ├── main.py          # Pipeline entry point
│   ├── config.py        # Configuration management
│   └── reporter.py      # Progress reporting
├── core/                # Core data structures
│   ├── intron.py        # Intron class and related types
│   └── reference.py     # Reference sequence management
├── extraction/          # Intron extraction from various sources
│   ├── annotation.py    # GFF3/GTF parsing
│   ├── bed.py           # BED file parsing
│   ├── sequences.py     # Sequence file parsing
│   └── filter.py        # Quality control and filtering
├── scoring/             # PWM scoring and normalization
│   ├── pwm.py           # Position-weight matrix operations
│   ├── scorer.py        # Score calculation
│   └── normalizer.py    # Z-score normalization
├── classification/      # SVM training and prediction
│   ├── trainer.py       # Model training with nested CV
│   ├── predictor.py     # Ensemble prediction
│   ├── nested_cv.py     # Nested cross-validation
│   └── split_eval.py    # Evaluation utilities
├── output/              # Output file generation
│   ├── writers.py       # All output writers
│   └── formatter.py     # Formatting utilities
├── visualization/       # Plotting functions
│   └── plots.py         # All visualization code
├── utils/               # Utility modules
│   ├── genome.py        # Genome file handling
│   ├── logging_utils.py # Enhanced logging
│   └── sequences.py     # Sequence utilities
├── __main__.py          # Module entry point
├── pixi.toml            # Pixi configuration
├── pyproject.toml       # Python project metadata
└── README.md            # This file
```

---

## Citing intronIC

If you use `intronIC` in your research, please cite:

**Devlin C Moyer, Graham E Larue, Courtney E Hershberger, Scott W Roy, Richard A Padgett.** *Comprehensive database and evolutionary dynamics of U12-type introns.* **Nucleic Acids Research,** Volume 48, Issue 13, 27 July 2020, Pages 7066–7078. <https://doi.org/10.1093/nar/gkaa464>

---

## About intronIC

`intronIC` was created to provide a customizable, open-source method for identifying minor (U12-type) spliceosomal introns from genomic data. U12-type introns are rare (~0.5% of introns) but functionally important, and contain distinct splicing motifs that make them amenable to computational identification.

### Why intronIC?

**Earlier U12 databases** (U12DB, SpliceRack, ERISdb) were valuable resources but:
- Static by design (not updated with new genome releases)
- Based on older genome annotations
- Limited to pre-selected species
- Used heuristic classification criteria

**intronIC addresses these limitations:**
- Works with any genome + annotation
- Uses the well-established SVM classification approach
- Produces interpretable probability scores
- Allows customization of training data and parameters
- Provides extensive metadata for downstream analysis
- Regularly updated with algorithm improvements

### Classification Method

intronIC's approach combines sequence motif analysis with machine learning:

1. **Position-Weight Matrices (PWMs)**: Capture sequence preferences at three key regions
   - 5' splice site (donor): Recognizes GT/AT at intron start
   - Branch point: Identifies TCCTTAAC-like motifs in U12-type introns
   - 3' splice site (acceptor): Recognizes AG/AC at intron end

2. **Z-Score Normalization**: Converts raw PWM scores to standardized features
   - Fit on reference sequences only (prevents data leakage)
   - Accounts for different score ranges across regions

3. **Linear SVM Classifier**: Learns decision boundary in 3D feature space
   - Trained on curated U12-type and U2-type reference sets
   - Balanced class weights handle imbalanced data (~0.5% expected U12-type)
   - Probability calibration provides confidence estimates

4. **Ensemble Averaging**: Reduces variance through multiple models
   - Each model trained on different U2 subsamples
   - F1-weighted voting combines predictions
   - Produces robust, reliable probabilities

This approach avoids arbitrary score thresholds and provides probabilistic classifications that researchers can interpret based on their specific needs (e.g., high-confidence predictions for experimental validation vs. comprehensive catalogs).

### The Refactoring

This refactored version maintains complete algorithmic fidelity to the original while dramatically improving code organization and maintainability. The original 6,093-line monolithic file has been restructured into 15+ focused modules, each with a single responsibility.

Key improvements include:
- Fixed data leakage bug in z-score normalization
- Corrected type_id assignment logic
- Added comprehensive type hints
- Immutable data structures for thread safety
- Better logging and error messages
- Structured for testing and extension

For complete details, see [REFACTOR_SUMMARY.md](REFACTOR_SUMMARY.md).

---

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Make your changes with tests
4. Submit a pull request

For major changes, please open an issue first to discuss the proposed changes.

---

## License

`intronIC` is released under the [GNU General Public License v3.0](LICENSE).

---

## Acknowledgments

Developed by **Graham E. Larue** with contributions from the Roy Lab and Padgett Lab.

Reference database curation: **Devlin C. Moyer, Courtney E. Hershberger**

Special thanks to the bioinformatics community for tools and libraries that make this work possible.

---

**For more detailed documentation, algorithm descriptions, and examples, visit the [intronIC wiki](https://github.com/glarue/intronIC/wiki).**
