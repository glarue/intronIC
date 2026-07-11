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

## What's New in v3.0.0

v3.0.0 replaces the entire z-normalization + mode-separation + continuous-discount stack
that accreted through v2.5–v2.7 with the **raw-feature `pmotif_adjudicated`** architecture.
The eval-corpus work showed that per-species z-normalization actually *hurts* cross-species
discrimination, so all per-species adaptation was moved out of scoring and into a single
output-level adjudicator.

- **Species-agnostic motif probability (`P_motif`)**: a calibrated 42-model RBF-SVM ensemble
  scores the six background-corrected **raw** motif features (`5'_raw`, `BP_raw`, `3'_raw`,
  `bp_offset`, `bp_scan_confidence`, `support2_raw`) → a per-intron ensemble margin → a
  Platt-calibrated `P_motif`. No per-species z-scaler, no second scoring pass.
- **Species adjudicator → `motif_category`**: a per-species `z_excess` (Poisson significance
  of the strong-call *count* against the genome's own U2-type tail) plus a per-genome strength
  gate decide `motif_category` ∈ {`DETECTED`, `INCONCLUSIVE`, `NOT_DETECTED`, `UNASSESSABLE`}.
- **Label rule**: `type_id = u12` iff `adjusted_score ≥ 50` (equivalently `P_motif ≥ 0.5`)
  **and** `motif_category ≠ NOT_DETECTED`. Only `NOT_DETECTED` suppresses calls;
  `INCONCLUSIVE` / `UNASSESSABLE` flag species-level ambiguity but still let strong-motif
  introns call.
- **Removed** modules and now-inert flags: the normalizer, mode-separation, and
  continuous-discount code is gone; `--load-normalizer` / `--save-normalizer` / `--no-mode-sep`
  / `--no-continuous-discount` / `--discount-*` are accepted-but-ignored no-ops. Legacy
  z/mode-sep bundles remain reproducible from the `pre-zstack-removal` git tag.
- **New outputs**: `score_info.iic` gains `P_motif`, `z_excess`, `cs_p95`, `p_gumbel_p95`,
  `motif_category` (and related) columns; a per-species `tail_model.iic.json` adjudicator
  sidecar + diagnostic figure; authoritative output-format spec under
  [`docs/output_formats.md`](docs/output_formats.md).
- **`--streaming` (default) and `--in-memory` remain bit-identical**; mode choice affects
  only the runtime/memory tradeoff. Dropping the second scoring pass makes v3 substantially
  faster than the v2.6–v2.7 two-pass architecture.

Full release history — including the superseded v2.3–v2.7 z-normalization / mode-separation /
continuous-discount architecture — is in [CHANGELOG.md](CHANGELOG.md).

---

## Key Features

- **Species-agnostic probability** (`P_motif`, 0–1; also as `adjusted_score` 0–100) from a calibrated 42-model RBF-SVM ensemble on background-corrected raw motif features, plus a per-species adjudicator (`motif_category`) that gates whether a genome carries a detectable U12-type population
- **Pretrained model** loaded automatically for cross-species analysis
- **Streaming mode** (default) roughly halves peak memory on large genomes (e.g., ~5.3 GB for full human at `-p 5`); bit-identical to in-memory
- **Parallel scoring** via `-p N` for linear speedup
- **Comprehensive metadata**: phase, position, parent gene/transcript

---

## How It Works

Most eukaryotic introns (~99.5%) are spliced by the **major (U2-type) spliceosome**; a small fraction (~0.5%) are spliced by the **minor (U12-type) spliceosome**. U12-type introns carry a conserved **TCCTTAAC** branch point motif and have either **AT-AC** (~25%) or **GT-AG** (~75%) terminal dinucleotides.

intronIC identifies U12-type introns in five stages:

1. **Extract** — pull introns from the annotation's CDS/exon coordinates (isoform handling, canonical/non-canonical terminus filters, partial-CDS handling)
2. **PWM scoring** — score the 5' splice site, branch point, and 3' splice site against position-weight matrices trained on a U12-type reference set
3. **Background correction** — correct the motif log-odds against the species' own intron pool to remove compositional bias, yielding the **raw** motif features
4. **Classification** — score the raw features through the calibrated 42-model RBF-SVM ensemble → per-intron ensemble margin → **`P_motif`** (species-agnostic, Platt-calibrated)
5. **Adjudication** — at the species level, compute `z_excess` and the strength gate → **`motif_category`**, then label `type_id = u12` where `P_motif ≥ 0.5` and the genome is not `NOT_DETECTED`

See [Technical Details](https://github.com/glarue/intronIC/wiki/Technical-algorithm) in the wiki for the full algorithm description.

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
