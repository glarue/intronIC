![intronIC_logo](https://user-images.githubusercontent.com/6827531/82829967-62872480-9e69-11ea-94e9-fa7306c7df1b.png)

# intronIC (intron <ins>I</ins>nterrogator and <ins>C</ins>lassifier)

Classify intron sequences as **U12-type** (minor spliceosome) or **U2-type** (major spliceosome). A 42-model RBF SVM ensemble scores each intron against position-weight matrices and outputs a calibrated probability (0-100%).

---

## Quick Start

```bash
pip install intronIC
```

```bash
# Classify introns (loads default model automatically)
intronIC -g genome.fa.gz -a annotation.gff3.gz -n species_name -p 8

# Extract sequences without classification
intronIC extract -g genome.fa.gz -a annotation.gff3.gz -n species_name -p 8

# Verify installation with bundled test data
intronIC test -p 4
```

---

## What's New in v2.3

- **42-model RBF SVM ensemble** on a streamlined 6D feature set
- **Bayesian score adjustment** suppresses false positives in species lacking a distinct U12-type intron population, using a species-level valley prior and per-intron ensemble agreement
- **Species-specific U2-type background correction** for cross-species composition bias
- **Default threshold raised to 95%** for higher-confidence calls
- See [CHANGELOG.md](CHANGELOG.md) for full release history

---

## Key Features

- **Probability scores** (0-100%) from a 42-model ensemble with isotonic calibration
- **Pretrained model** loaded automatically for cross-species analysis
- **Streaming mode** (default) reduces memory ~85% on large genomes
- **Parallel scoring** via `-p N` for linear speedup
- **Comprehensive metadata**: phase, position, parent gene/transcript

---

## How It Works

Most eukaryotic introns (~99.5%) use the **major (U2-type) spliceosome**. A small fraction (~0.5%) use the **minor (U12-type) spliceosome**, characterized by a conserved **TCCTTAAC** branch point motif and either **AT-AC** (~25%) or **GT-AG** (~75%) terminal dinucleotides.

intronIC identifies U12-type introns in five stages:

1. **PWM scoring** — score the 5' splice site, branch point, and 3' splice site against position-weight matrices
2. **Background correction** — blend species-specific nucleotide frequencies into U2-type PWMs to correct composition bias
3. **Normalization** — convert raw log-odds to z-scores via robust scaling
4. **SVM classification** — 42-model RBF SVM ensemble produces per-intron probabilities and ensemble agreement (sigma)
5. **Score adjustment** — adjust probabilities using a species-level valley prior and an ensemble disagreement penalty

See [Technical Details](https://github.com/glarue/intronIC/wiki/Technical-algorithm) for the full algorithm description.

---

## Documentation

Full documentation lives in the **[intronIC Wiki](https://github.com/glarue/intronIC/wiki)**:

- **[Quick Start](https://github.com/glarue/intronIC/wiki/Quick-start)** — Installation, dependencies, resource usage
- **[Overview](https://github.com/glarue/intronIC/wiki/Overview)** — Classification approach and scientific background
- **[Output Files](https://github.com/glarue/intronIC/wiki/Output-files)** — File formats and score interpretation
- **[Technical Details](https://github.com/glarue/intronIC/wiki/Technical-algorithm)** — Algorithm, features, score adjustment
- **[Usage Info](https://github.com/glarue/intronIC/wiki/Usage-info)** — Complete CLI reference
- **[Example Usage](https://github.com/glarue/intronIC/wiki/Example-usage)** — Common workflows
- **[Changelog](CHANGELOG.md)** — Release notes and version history

---

## Citation

If you use intronIC in your research, please cite:

> Moyer DC, Larue GE, Hershberger CE, Roy SW, Padgett RA. (2020) Comprehensive database and evolutionary dynamics of U12-type introns. *Nucleic Acids Research* 48(13):7066-7078. [doi:10.1093/nar/gkaa464](https://doi.org/10.1093/nar/gkaa464)

---

## Support

- [intronIC Wiki](https://github.com/glarue/intronIC/wiki) — Documentation
- [GitHub Issues](https://github.com/glarue/intronIC/issues) — Bug reports
- [GitHub Discussions](https://github.com/glarue/intronIC/discussions) — Questions and ideas

---

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

```bash
git clone https://github.com/glarue/intronIC.git
cd intronIC
make install    # Set up development environment
make test       # Run tests
```

---

## License

[GNU General Public License v3.0](LICENSE)
