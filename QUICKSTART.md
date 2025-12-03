# intronIC Quick Start Guide

## Installation (30 seconds)

```bash
# Clone repository
git clone https://github.com/glarue/intronIC.git
cd intronIC/intronIC_refactored

# Install
pip install .

# Verify
intronIC --version
```

## Basic Usage

### 1. View Help
```bash
intronIC --help
```

### 2. Extract and Classify Introns
```bash
intronIC \
  --genome genome.fa.gz \
  --annotation annotation.gff3.gz \
  --species-name homo_sapiens \
  --processes 8
```

**Short form:**
```bash
intronIC -g genome.fa.gz -a annotation.gff3.gz -n homo_sapiens -p 8
```

### 3. Extract Sequences Only (No Classification)
```bash
intronIC -g genome.fa.gz -a annotation.gff3.gz -n species_name --sequences-only
```

### 4. Use Pretrained Model
```bash
intronIC -g genome.fa.gz -a annotation.gff3.gz -n species_name \
  --model homo_sapiens.model.pkl
```

## Key Options

| Option | Short | Description |
|--------|-------|-------------|
| `--genome` | `-g` | Genome FASTA file (required) |
| `--annotation` | `-a` | GFF3/GTF annotation (required) |
| `--species-name` | `-n` | Output prefix (required) |
| `--processes` | `-p` | Number of CPU cores (default: 1) |
| `--sequences-only` | `-s` | Skip classification |
| `--threshold` | `-t` | U12 probability threshold (default: 90) |
| `--model` | | Use existing model |
| `--recursive` | | Recursive training for distant species |

## Output Files

After running, you'll get:

- `{species}.meta.iic` - Complete metadata for all introns
- `{species}.bed.iic` - BED format coordinates
- `{species}.seqs.iic` - Intron sequences
- `{species}.scores.iic` - Detailed scoring information
- `{species}.log.iic` - Run log with statistics

## Examples

### Human Chr19 (Test Run)
```bash
intronIC \
  -g test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
  -a test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz \
  -n human_chr19_test \
  -p 4
```

### Extract Only, All Isoforms
```bash
intronIC \
  -g genome.fa.gz \
  -a annotation.gff3.gz \
  -n species_seqs \
  --sequences-only \
  --allow-multiple-isoforms
```

### Strict Threshold
```bash
intronIC \
  -g genome.fa.gz \
  -a annotation.gff3.gz \
  -n species_strict \
  --threshold 95 \
  -p 8
```

### Distant Species (Recursive Training)
```bash
intronIC \
  -g non_model_genome.fa.gz \
  -a non_model_annotation.gff3.gz \
  -n distant_species \
  --recursive \
  -p 8
```

## Understanding Output

### U12 Introns
Introns with **SVM score ≥ threshold** (default 90%) are classified as U12-type.

Look for entries in `.meta.iic` with:
- High relative scores (> 0)
- `type_id` = "u12"

### Statistics
Check the `.log.iic` file for:
- Total introns extracted
- Number filtered (and why)
- U12 count and percentage
- Runtime statistics

## Common Issues

**"Command not found: intronIC"**
→ Run `pip install .` again or use `python -m intronIC`

**"No module named 'networkx'"**
→ Dependencies not installed. Run `pip install .`

**Very slow performance**
→ Use `-p 8` (or more) for parallel processing

**No U12 introns found**
→ This may be normal for some species. Try `--recursive` for non-model organisms.

## Getting Help

```bash
intronIC --help              # Full usage
intronIC --version           # Check version
```

See also:
- [INSTALL.md](INSTALL.md) - Detailed installation guide
- [README.md](README.md) - Complete documentation
- [GitHub Issues](https://github.com/glarue/intronIC/issues) - Report problems

## Citation

If you use intronIC in your research, please cite:

Moyer et al. (2020). "Comprehensive database and evolutionary dynamics of U12-type introns."
*Nucleic Acids Research*, 48(13):7066–7078.
https://doi.org/10.1093/nar/gkaa464
