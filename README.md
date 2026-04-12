![intronIC_logo](https://user-images.githubusercontent.com/6827531/82829967-62872480-9e69-11ea-94e9-fa7306c7df1b.png)

# intronIC - (intron <ins>I</ins>nterrogator and <ins>C</ins>lassifier)

`intronIC` is a bioinformatics tool for extracting and classifying intron sequences as **U12-type (minor)** or **U2-type (major)** using a support vector machine trained on position-weight matrix scores.

---

## Quick Start

### Installation

```bash
pip install intronIC
```

### Basic Usage

```bash
# Classify introns (default model loaded automatically)
intronIC -g genome.fa.gz -a annotation.gff3.gz -n species_name -p 8

# Extract sequences only (no classification)
intronIC extract -g genome.fa.gz -a annotation.gff3.gz -n species_name -p 8

# Train a custom model (optional - most users don't need this)
intronIC train -n my_model -p 8
```

### Test Run

```bash
# Quick installation test using bundled test data
intronIC test -p 4

# Or show where test data is located
intronIC test --show-only
```

---

## Documentation

* **[Changelog](CHANGELOG.md)** - Release notes and version history

For complete documentation, see the **[intronIC Wiki](https://github.com/glarue/intronIC/wiki)**:

* **[Quick Start Guide](https://github.com/glarue/intronIC/wiki/Quick-start)** - Installation, dependencies, resource usage
* **[Overview](https://github.com/glarue/intronIC/wiki/Overview)** - Classification approach and scientific background
* **[Usage Info](https://github.com/glarue/intronIC/wiki/Usage-info)** - Complete CLI reference
* **[Output Files](https://github.com/glarue/intronIC/wiki/Output-files)** - File formats and interpretation
* **[Technical Details](https://github.com/glarue/intronIC/wiki/Technical-algorithm)** - Algorithm and ML architecture
* **[Example Usage](https://github.com/glarue/intronIC/wiki/Example-usage)** - Common workflows
* **[About](https://github.com/glarue/intronIC/wiki/About)** - Background and motivation

---

## What's New in v2.2

- **New 8D RBF SVM default model** trained on expanded reference data (472 U12-type + 30,155 U2-type introns)
- **Five new classification features**: branch point offset, BPS motif sharpness, polypyrimidine tract metrics, and multi-site support scoring
- **Reduced false positives**: 0 confident false calls in *C. elegans* (was 2), 1 in *Ascaris* (was 47)
- See [CHANGELOG.md](CHANGELOG.md) for full release history

---

## Key Features

- **RBF SVM classification** with probability scores (0-100%) using 8 sequence-derived features
- **Default pretrained model** loaded automatically — works for virtually all species
- **Streaming mode** (default) for ~85% memory reduction on large genomes
- **Parallel processing** for improved performance (`-p 8` recommended)
- **Fast runtimes**: ~6-10 minutes for human genome with default settings
- **Comprehensive metadata** including phase, position, parent gene/transcript

---

## Scientific Background

Most eukaryotic introns (~99.5%) are spliced by the **major (U2-type) spliceosome**, while a small fraction (~0.5%) are spliced by the **minor (U12-type) spliceosome**. U12-type introns have:

- Highly conserved **TCCTTAAC** branch point motif
- Terminal dinucleotides: **AT-AC** (~25%) or **GT-AG** (~75%)
- Functional importance and evolutionary conservation

intronIC identifies U12-type introns using:

1. **PWM Scoring**: Apply position-weight matrices to 5' splice site, branch point, and 3' splice site regions
2. **Normalization**: Convert raw scores to z-scores via robust scaling
3. **Feature Engineering**: Compute composite features (multi-site corroboration, BP position, PPT metrics, BPS motif sharpness)
4. **SVM Classification**: RBF SVM ensemble with balanced class weights outputs probability scores

For detailed algorithm description, see the [Technical Details](https://github.com/glarue/intronIC/wiki/Technical-algorithm) wiki page.

---

## Citation

If you use `intronIC` in your research, please cite:

**Devlin C Moyer, Graham E Larue, Courtney E Hershberger, Scott W Roy, Richard A Padgett.** *Comprehensive database and evolutionary dynamics of U12-type introns.* **Nucleic Acids Research,** Volume 48, Issue 13, 27 July 2020, Pages 7066–7078. <https://doi.org/10.1093/nar/gkaa464>

---

## Support

* **Documentation**: [intronIC Wiki](https://github.com/glarue/intronIC/wiki)
* **Issues**: [GitHub Issues](https://github.com/glarue/intronIC/issues)
* **Discussions**: [GitHub Discussions](https://github.com/glarue/intronIC/discussions)

---

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

```bash
git clone https://github.com/glarue/intronIC.git
cd intronIC
make install    # Set up development environment
make test       # Run tests
```

---

## License

`intronIC` is released under the [GNU General Public License v3.0](LICENSE).
