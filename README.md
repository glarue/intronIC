![intronIC_logo](https://user-images.githubusercontent.com/6827531/82829967-62872480-9e69-11ea-94e9-fa7306c7df1b.png)

# intronIC (intron <ins>I</ins>nterrogator and <ins>C</ins>lassifier)

Classify intron sequences as **U12-type** (minor spliceosome) or **U2-type** (major spliceosome). A two-pass RBF SVM pipeline (first-pass cluster-aware + second-pass per-species mode-separation, each a 126-model multispecies ensemble) scores each intron against position-weight matrices and outputs a calibrated probability (0-100%) along with a continuous per-intron overcall discount.

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

## What's New in v2.7

- **Continuous per-intron discount** (`adjusted_score` column): non-positive
  log-odds penalty for SVM overcalls relative to motif log-LR or for weak
  motif evidence. Empirically derived defaults; preserves panel TPs ≥99%
  while trimming the long-tail of loose-or-NA calls.
- **`adjusted_score` is the new recommended call column** — `svm_score`
  remains the raw classifier output (auditability preserved).
- **Diagnostic surface**: per-intron `raw_sum`, `svm_vs_naive`, `voting_frac`
  columns added to score_info.iic; per-species `boundary_mass` reported in
  `.modesep.json` (diagnostic only — no gate role).
- New CLI flags: `--no-continuous-discount`, `--discount-k-overcall`,
  `--discount-tau-overcall`, `--discount-k-weakmot`, `--discount-tau-motif`.
- 628 unit + 13 integration tests passing.

## What's New in v2.6

- **Mode-separation classifier (default)**: per-species recalibration places
  the U2 mode at z=0 and U12 mode at z=1 in every species. Plant recall
  jumps (AmbTri 90% → 100%, OrySat 94% → 100%) and Apostasia IPA recall
  17/21 → 20/21 without inflating false positives. See
  [Technical Details](https://github.com/glarue/intronIC/wiki/Technical-algorithm#stage-4-two-pass-classification)
  in the wiki for the architecture.
- **Three-check gate** (n_eff floor + μ_U12 location prior + multi-bandwidth
  Fisher-KDE valley depth) protects against the failure modes of per-species
  recalibration; U12-absent species fall back to first-pass scores cleanly.
- **Diagnostic JSON sidecar** (`.modesep.json`) with route, gate reason,
  μ_U2/U12, valley depth, ensemble σ on called introns, and an A/B/C/F
  quality tier per species. Per-intron `ensemble_sigma`, `first_pass_svm`,
  and `modesep_route` columns are added to `score_info.iic`.
- New CLI flags: `--no-mode-sep`, `--mode-sep-z-floor`,
  `--mode-sep-valley-min`, `--mode-sep-n-floor`,
  `--mode-sep-mu-u12-tolerance`.
- v4 cluster-aware bundles still load; pre-v2.6 behavior preserved for them.

## What's New in v2.4

- **Default model is now the v3 multispecies bundle**: 3 seeds × 42 calibrated SVMs (126 total) trained on 41,333 introns across 90 species and 14 clades. Holdout F1 = 1.000 vs the v2.3 default's 0.9975, and ~54% lower production-equivalent FPR on U12-absent species.
- **Default classification threshold lowered from 95 → 90**, made safe by the v3 model's tighter calibration. Pass `--threshold 95` to restore prior behavior.
- **`--streaming` (default) and `--in-memory` now produce bit-identical classifications**. Mode choice affects only the runtime/memory tradeoff. Reference run (v2.7, mode-sep two-pass) on Homo sapiens GRCh38.p13 + NCBI RefSeq GFF, `-p 5`, 257k scored introns: streaming ~40 min / 5.3 GB peak. In-memory was not re-measured for v2.7 (expected similar wall time at roughly 2× peak memory based on v2.4 ratios). The wall-time growth from v2.4 (~16 min) reflects the v2.6+ two-pass mode-separation architecture; both passes run the full 126-model ensemble.
- **Self-describing model bundles** carry config + training metadata alongside the weights; see [`docs/v3_bundle_schema.md`](docs/v3_bundle_schema.md).
- v2.3 model bundles continue to load unchanged; old runs reproduce by passing `--model <v2.3-bundle.pkl>`.
- See [CHANGELOG.md](CHANGELOG.md) for full release history.

## What's New in v2.3

- **42-model RBF SVM ensemble** on a streamlined 6D feature set
- **Bayesian score adjustment** suppresses false positives in species lacking a distinct U12-type intron population, using a species-level valley prior and per-intron ensemble agreement
- **Species-specific U2-type background correction** for cross-species composition bias
- **Default threshold raised to 95%** for higher-confidence calls (now lowered to 90 in v2.4)

---

## Key Features

- **Probability scores** (0-100%) from two 126-model calibrated SVM ensembles (3 seeds × 42 sub-models each, isotonic calibration) — a first-pass cluster-aware classifier and a second-pass per-species mode-separation classifier
- **Pretrained model** loaded automatically for cross-species analysis
- **Streaming mode** (default) roughly halves peak memory on large genomes (e.g., ~5.3 GB for full human at `-p 5`); bit-identical to in-memory
- **Parallel scoring** via `-p N` for linear speedup
- **Comprehensive metadata**: phase, position, parent gene/transcript

---

## How It Works

Most eukaryotic introns (~99.5%) are spliced by the **major (U2-type) spliceosome**; a small fraction (~0.5%) are spliced by the **minor (U12-type) spliceosome**. U12-type introns carry a conserved **TCCTTAAC** branch point motif and have either **AT-AC** (~25%) or **GT-AG** (~75%) terminal dinucleotides.

intronIC v2.7 identifies U12-type introns in seven stages:

1. **PWM scoring** — score the 5' splice site, branch point, and 3' splice site against position-weight matrices
2. **Background correction** — blend species-specific nucleotide frequencies into U2-type PWMs to correct composition bias
3. **Adaptive normalizer fit** — score sampled introns and fit a per-species robust z-scaler (median/IQR) for the first-pass features
4. **First-pass classification** — score every intron through the 126-model cluster-aware RBF SVM ensemble (`v4_aug`); produces `first_pass_svm` and the candidate weights used to estimate per-species U12/U2 modes
5. **Mode estimation + gate** — estimate per-species μ_U12 / μ_U2 from soft candidate weights; gate against three checks (n_eff floor, μ_U12 location prior, Fisher-discriminant KDE valley depth)
6. **Second-pass classification (mode-separation)** — on gate-pass, re-z-score motif features so U2 → 0 and U12 → 1 in every species, then score eligible introns through the 126-model `v5_modesep_aug` ensemble (`svm_score`). On gate-fail, keep first-pass scores and apply the legacy Bayesian valley-depth + ensemble-agreement adjustment.
7. **Continuous per-intron discount** — apply a non-positive log-odds penalty for SVM overcalls relative to motif log-LR; produces `adjusted_score` (the calling column).

See [Technical Details](https://github.com/glarue/intronIC/wiki/Technical-algorithm) in the wiki for the full algorithm description, including the [two-pass mode-separation architecture](https://github.com/glarue/intronIC/wiki/Technical-algorithm#stage-4-two-pass-classification) and the [v2.7 continuous discount](https://github.com/glarue/intronIC/wiki/Technical-algorithm#v27-continuous-discount-per-intron-overcall-penalty).

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
