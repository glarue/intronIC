# intronIC Quick Start Guide

> **Note:** This guide has been integrated into the main [README.md](README.md#quick-start). Please refer to the README for the most up-to-date quick start instructions.

For quick reference, here are the essential commands:

## Installation

```bash
git clone https://github.com/glarue/intronIC.git
cd intronIC
pip install -e .
```

## Basic Commands

```bash
# Classify introns (train on-the-fly)
intronIC -g genome.fa.gz -a annotation.gff3.gz -n species_name -p 8

# Use pretrained model (faster)
intronIC -g genome.fa.gz -a annotation.gff3.gz -n species_name --model model.pkl -p 8

# Extract sequences only
intronIC extract -g genome.fa.gz -a annotation.gff3.gz -n species_name -p 8

# Train a model
intronIC train -n my_model -p 8
```

## Full Documentation

See the [README.md](README.md) for:
- Complete usage documentation
- All command-line options
- Output file descriptions
- Advanced usage examples
- Troubleshooting guide
