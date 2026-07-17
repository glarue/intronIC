![intronIC_logo](https://user-images.githubusercontent.com/6827531/82829967-62872480-9e69-11ea-94e9-fa7306c7df1b.png)

# intronIC (intron <ins>I</ins>nterrogator and <ins>C</ins>lassifier)

Classify intron sequences as **U12-type** (minor spliceosome) or **U2-type** (major spliceosome). A calibrated 42-model RBF-SVM ensemble scores each intron's splice-site and branch-point motifs — background-corrected against the species' own intron pool — into a species-agnostic probability (`P_motif`), and a per-species adjudicator decides whether the genome carries a detectable U12-type population.

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

## How It Works

Most eukaryotic introns (~99.5%) are spliced by the **major (U2-type) spliceosome**; a small fraction (~0.5%) are spliced by the **minor (U12-type) spliceosome**. U12-type introns carry a conserved **TCCTTAAC** branch point motif and have either **AT-AC** (~25%) or **GT-AG** (~75%) terminal dinucleotides. intronIC scores those motifs and classifies each intron.

Given a genome FASTA and a GFF/GTF annotation, intronIC runs five stages:

1. **Extract** — pull introns from the annotation's CDS/exon coordinates (isoform handling, canonical/non-canonical terminus filters, partial-CDS handling).
2. **Score motifs** — score each intron's 5′ splice site, branch point, and 3′ splice site against position-weight matrices built from a U12-type reference set.
3. **Background-correct** — score against a background model built from the species' own intron pool, so a high score reflects genuine U12-type signal rather than that genome's base composition. These are the **raw** motif features.
4. **Classify** — feed the raw features to a calibrated 42-model RBF-SVM ensemble → a per-intron ensemble margin → **`P_motif`**, a single species-agnostic probability that the motif is U12-type (Platt-calibrated; also reported 0–100 as `adjusted_score`).
5. **Adjudicate** — at the genome level, test whether the strong U12-type calls form a real population standing above the species' own U2-type score tail (`z_excess`, plus a call-strength gate) → **`motif_category`** ∈ {`DETECTED`, `INCONCLUSIVE`, `NOT_DETECTED`, `UNASSESSABLE`}. An intron is labelled `type_id = u12` when `P_motif ≥ 0.5` and the genome is not `NOT_DETECTED`.

See [Technical Details](https://github.com/glarue/intronIC/wiki/Technical-algorithm) in the wiki for the full algorithm.

---

## Relationship to the published method (v1)

intronIC v3 keeps the PWM motif scoring of the Moyer et al. (2020) method but replaces the per-species z-normalization that made scores comparable across genomes. That normalization assumed every genome matched the human reference geometry, and failed in two ways: it under-called divergent bearers (e.g. *Amborella*, *Oryza*) and produced false positives in genomes that have lost the minor spliceosome. v3 instead scores background-corrected **raw** motif log-odds — six features, an ensemble of 42 SVMs trained on 97 species across 14 clades — and calibrates the output to `P_motif`, a species-agnostic probability. All per-species reasoning moves to one output-level adjudicator, which reports whether a genome carries a detectable U12-type population (`motif_category`), referenced to that genome's own U2-type background.

| | v1 (Moyer 2020) | v3 (current) |
|---|---|---|
| Features | 3 motif scores | 6 raw features (→9 with interactions) |
| Classifier | 1 human-trained SVM | 42-model ensemble, 97 species |
| Cross-species | per-species z-normalization | species-agnostic `P_motif` |
| Per-species logic | folded into normalization | output-level adjudicator |

For the full comparison, see [Technical Details → Relationship to the published method (v1)](https://github.com/glarue/intronIC/wiki/Technical-algorithm#relationship-to-the-published-method-v1).

---

## Key Features

- **Pretrained, cross-species model** loaded automatically — no per-run training or normalization step.
- **`P_motif`** (0–1) / **`adjusted_score`** (0–100) per intron, plus a genome-level **`motif_category`** call.
- **Streaming mode** (default) roughly halves peak memory on large genomes (~5.3 GB for full human at `-p 5`) and is bit-identical to in-memory.
- **Parallel scoring** via `-p N`, and comprehensive per-intron metadata (phase, transcript position, parent gene/transcript).

---

## Documentation

Full documentation lives in the **[intronIC Wiki](https://github.com/glarue/intronIC/wiki)**:

- **[Quick Start](https://github.com/glarue/intronIC/wiki/Quick-start)** — Installation, dependencies, resource usage
- **[Overview](https://github.com/glarue/intronIC/wiki/Overview)** — Classification approach and scientific background
- **[Output Files](https://github.com/glarue/intronIC/wiki/Output-files)** — File formats and score interpretation
- **[Technical Details](https://github.com/glarue/intronIC/wiki/Technical-algorithm)** — Algorithm, features, species adjudication
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
