![intronIC_logo](https://user-images.githubusercontent.com/6827531/82829967-62872480-9e69-11ea-94e9-fa7306c7df1b.png)

# intronIC Refactored - (intron <ins>I</ins>nterrogator and <ins>C</ins>lassifier)

**Version 1.5.1** - Refactored with improved performance

`intronIC` is a program that can be used to classify intron sequences as minor (U12-type) or major (U2-type), using a genome and annotation or the sequences themselves. Alternatively, `intronIC` can be used to simply extract all intron sequences without classification (using `-s`).

## What's New in the Refactored Version

This refactored version maintains **full CLI compatibility** with the original while providing:

- **5-10x faster SVM training** via `SVC(probability=False)` + `CalibratedClassifierCV`
- **Grid-searched calibration** (sigmoid vs isotonic methods)
- **Better probability scoring** using `neg_log_loss` metric
- **Modular architecture** for easier maintenance and testing
- **Type hints** throughout the codebase
- **Modern package management** support (pixi, uv)

### Performance Comparison

| Operation | Original | Refactored | Speedup |
|-----------|----------|------------|---------|
| Full optimization (Chr19) | 25-40 min | 2-5 min | 5-10x |
| Small ref optimization | 5-8 min | 1-2 min | 3-5x |
| Single CV fold | 20-120s | 1-15s | 10-20x |

---

## Quick Start

### Using pixi (Recommended)

```bash
# Install pixi
curl -fsSL https://pixi.sh/install.sh | bash

# Navigate to refactored directory
cd intronIC_refactored

# Install and run on test data (fast)
pixi install
pixi run test-small
```

### Using uv (Fast Alternative)

```bash
# Install uv
curl -LsSf https://astral.sh/uv/install.sh | sh

# Setup and run
cd intronIC_refactored
uv venv
source .venv/bin/activate
uv pip install -e .

python -m intronIC_refactored -g ../intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
  -a ../intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz -n test_run
```

### Using pip (Traditional)

```bash
cd intronIC_refactored
pip install -e .

intronIC -g genome.fa.gz -a annotation.gff3.gz -n species_name
```

---

## Documentation

**For complete setup and usage instructions, see [SETUP.md](SETUP.md)**

The SETUP.md file includes:
- Detailed installation options (pixi, uv, pip)
- Running methods and examples
- Common use cases with small/full reference sets
- Performance tuning tips
- Troubleshooting guide
- Testing instructions

For general information about intronIC's features and algorithms, see the [original wiki](https://github.com/glarue/intronIC/wiki).

---

## Key Features

- **Three input modes**: Genome + annotation, BED file, or pre-extracted sequences
- **Binary classification**: U2-type (major, ~99.5%) vs U12-type (minor, ~0.5%)
- **SVM-based scoring**: Probability scores for each intron (0-100%)
- **Parallel processing**: Multi-core support for faster runtime
- **Comprehensive output**: Multiple file formats with detailed metadata
- **Custom training**: Support for species-specific reference data

---

## Common Usage Examples

### Classify all introns (fast test with small references)

```bash
pixi run test-small
# Or:
python -m intronIC_refactored \
  -g ../intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
  -a ../intronIC/test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz \
  -n homo_sapiens \
  --reference_u12s ../intronIC/data/u12_reference_small.introns.iic.gz \
  --reference_u2s ../intronIC/data/u2_reference_small.introns.iic.gz \
  -p 4
```

### Production run with full references

```bash
pixi run test-full
# Or:
intronIC -g genome.fa.gz -a annotation.gff3.gz -n species_name -p 8
```

### Extract sequences only (no classification)

```bash
pixi run test-seqs-only
# Or:
intronIC -g genome.fa.gz -a annotation.gff3.gz -n species_name -s
```

---

## CLI Arguments (Fully Compatible)

The refactored version maintains 100% CLI compatibility with the original. Key arguments:

**Required:**
- `-n species_name` - Name in binomial form (e.g., homo_sapiens)
- `-g genome.fa.gz -a annotation.gff3.gz` - Genome + annotation, **OR**
- `-g genome.fa.gz -b coordinates.bed` - Genome + BED file, **OR**
- `-q sequences.iic` - Pre-extracted sequences

**Common options:**
- `-p N` - Parallel processes (significantly reduces runtime)
- `-s` - Extract sequences only (skip classification)
- `--reference_u12s` / `--reference_u2s` - Custom reference sets
- `-t 90` - Classification threshold (default: 90%)
- `-i` - Include multiple isoforms (default: longest only)
- `--no_nc` - Exclude non-canonical introns

---

## Output Files

All output files use the `.iic` extension and contain tab-delimited data:

- **`.meta.iic`** - Comprehensive metadata for each intron
- **`.bed.iic`** - BED-format coordinates with scores
- **`.seqs.iic`** - Intron sequences with flanking exonic regions
- **`.scores.iic`** - Detailed scoring breakdown (PWM scores, z-scores)
- **`.dupe_map.iic`** - Maps duplicate introns to representatives
- **`.overlap_map.iic`** - Maps overlapping intron coordinates
- **`.png`** - Visualization plots (scatter, hexbin, PR curves)

### Identifying U12-type Introns

U12-type introns have **relative scores > 0** (equivalent to probability > 90% by default):

```bash
# Filter U12-type introns from meta file
awk '($2!="." && $2>0)' species_name.meta.iic

# Count U12-type introns
awk '($2!="." && $2>0)' species_name.meta.iic | wc -l
```

---

## Dependencies

- Python >=3.8, <3.13
- numpy >=1.19.0, <2.0
- scipy >=1.5.0
- scikit-learn >=0.22, <2.0
- biogl >=0.1.0
- matplotlib >=3.3.0
- networkx >=2.5.1
- rich >=10.0

All dependencies are automatically installed via pixi or uv/pip.

---

## Performance Notes

### Memory Usage
- Human genome (Ensembl 95): ~5 GB
- Most genomes: <2 GB

### Runtime
- **Refactored version (Chr19)**: 2-5 min with small refs, 5-15 min with full refs
- **Original version (Chr19)**: 25-40 min
- Scales with `-p` parallel processes

---

## Cite

If you find this tool useful, please cite:

Devlin C Moyer, Graham E Larue, Courtney E Hershberger, Scott W Roy, Richard A Padgett, Comprehensive database and evolutionary dynamics of U12-type introns, Nucleic Acids Research, Volume 48, Issue 13, 27 July 2020, Pages 7066–7078, <https://doi.org/10.1093/nar/gkaa464>

---

## About

`intronIC` was written to provide a customizable, open-source method for identifying minor (U12-type) spliceosomal introns from annotated intron sequences. Minor introns usually represent ~0.5% (at most) of a given genome's introns, and contain distinct splicing motifs which make them amenable to bioinformatic identification.

`intronIC` uses a support-vector machine (SVM) classification method trained on position-weight matrix (PWM) scores from three key regions:
- 5' splice site (donor)
- Branch point
- 3' splice site (acceptor)

This produces an easy-to-interpret "probability of being U12-type" score for each intron, avoiding the heuristic fuzziness of other classification methods.

The refactored version significantly improves training performance while maintaining the same high-quality classification accuracy, making it practical to run intronIC on large datasets or to iterate quickly during parameter tuning.

For more details about the scientific background, algorithms, and output formats, see the [original wiki](https://github.com/glarue/intronIC/wiki).
