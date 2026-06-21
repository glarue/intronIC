![intronIC_logo](https://user-images.githubusercontent.com/6827531/82829967-62872480-9e69-11ea-94e9-fa7306c7df1b.png)

# intronIC (intron <ins>I</ins>nterrogator and <ins>C</ins>lassifier)

Classify intron sequences as **U12-type** (minor spliceosome) or **U2-type** (major spliceosome). A two-pass RBF SVM pipeline (first-pass cluster-aware + second-pass per-species mode-separation, each a 42-model multispecies ensemble) scores each intron against position-weight matrices and outputs a calibrated probability (0-100%) along with a continuous per-intron overcall discount.

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

## What's New in v2.7.1

- **3× faster classification** via single-seed ensemble bundle. Default
  bundle dropped from 3 seeds × 42 = 126 models per pass to 1 × 42 = 42.
  Empirically validated bit-identical TP/FN/FP across 73-species IPA-gold
  panel + AmbTri (mode-sep recovery) + Symbio (divergent gate-fail) +
  Drosophila smoke test. Bundle file 82 MB → 27 MB.
- **Unified per-intron labels** in `score_info.iic`:
  - `type_id` ∈ {u12, u2} — single-source-of-truth binary call on
    `adjusted_score >= 50`. Supersedes prior inconsistent threshold logic.
  - `confidence` ∈ {strong, borderline} — within-call gradation
    (u12_strong: adj ≥ 90; u2_strong: first_pass_svm < 10, empirically
    derived from 73-species panel: captures 99.92% of true U2s with
    0.02% U12 loss).
  - `history` ∈ {stable, promoted, demoted} — pipeline-path signature
    (promoted = mode-sep rescue, demoted = discount/mode-sep suppression).
- **New `metrics.iic.json` counts**: `u12_strong_count`,
  `u12_borderline_count`, `u12_promoted_count`, `u2_strong_count`,
  `u2_borderline_count`, `u2_demoted_count`. `high_confidence_u12`
  retained as back-compat alias. Fixes prior misleading `u2_count`
  computed as `total − high_confidence_u12`.
- **Tier name rename**: `modesep.json` `quality_tier` strings changed
  from opaque single-character codes to self-descriptive names
  (`modesep_strong` / `modesep_standard` / `modesep_weak` /
  `first_pass_fallback`). "F" was easily misread as "the run failed".
- **New canonical doc**: [`docs/scoring_pipeline.md`](docs/scoring_pipeline.md)
  — Mermaid flowchart of all six stages + parameter reference + label rubric.
- Backfill script for already-classified species in
  `scripts/backfill_v271_labels.py`.

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
  μ_U2/U12, valley depth, ensemble σ on called introns, and a four-tier
  quality classification per species (renamed in v2.7.1 to
  `modesep_strong` / `modesep_standard` / `modesep_weak` /
  `first_pass_fallback`; previously single-character codes A/B/C/F).
  Per-intron `ensemble_sigma`, `first_pass_svm`, and `modesep_route`
  columns are added to `score_info.iic`.
- New CLI flags: `--no-mode-sep`, `--mode-sep-z-floor`,
  `--mode-sep-valley-min`, `--mode-sep-n-floor`,
  `--mode-sep-mu-u12-tolerance`.
- v4 cluster-aware bundles still load; pre-v2.6 behavior preserved for them.

## What's New in v2.4

- **Default model is now the v3 multispecies bundle**: 3 seeds × 42 calibrated SVMs (126 total) trained on 41,333 introns across 90 species and 14 clades. Holdout F1 = 1.000 vs the v2.3 default's 0.9975, and ~54% lower production-equivalent FPR on U12-absent species.
- **Default classification threshold lowered from 95 → 90**, made safe by the v3 model's tighter calibration. Pass `--threshold 95` to restore prior behavior.
- **`--streaming` (default) and `--in-memory` now produce bit-identical classifications**. Mode choice affects only the runtime/memory tradeoff. Reference run on Homo sapiens GRCh38.p13 + NCBI RefSeq GFF, `-p 5`: v2.7.0 (126 models per pass, mode-sep two-pass) ~40 min / 5.3 GB peak; v2.7.1 (42 models per pass) ~3× faster at the same parallelism, expected ~13-15 min. In-memory ratios from v2.4 suggest ~2× peak memory; not re-measured for v2.7.x. The wall-time growth from v2.4 (~16 min single-pass) reflects the v2.6+ two-pass mode-separation architecture.
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

- **Probability scores** (0-100%) from two 42-model calibrated SVM ensembles (single seed group each as of v2.7.1; isotonic calibration) — a first-pass cluster-aware classifier and a second-pass per-species mode-separation classifier
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
4. **First-pass classification** — score every intron through the 42-model cluster-aware RBF SVM ensemble (`v4_aug`); produces `first_pass_svm` and the candidate weights used to estimate per-species U12/U2 modes
5. **Mode estimation + gate** — estimate per-species μ_U12 / μ_U2 from soft candidate weights; gate against three checks (n_eff floor, μ_U12 location prior, Fisher-discriminant KDE valley depth)
6. **Second-pass classification (mode-separation)** — on gate-pass, re-z-score motif features so U2 → 0 and U12 → 1 in every species, then score eligible introns through the 42-model `v5_modesep_aug` ensemble (`svm_score`). On gate-fail, keep first-pass scores and apply the legacy Bayesian valley-depth + ensemble-agreement adjustment.
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
