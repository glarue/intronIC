# Changelog

All notable changes to intronIC will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [3.1.0] - 2026-07-27

### Added

- **`bg_fdr`** — a new report-only per-intron column in `score_info.iic`: for each strong call, the
  fraction of the calls at least that strong which the genome's *own* U2 tail already explains. It is the
  per-intron decomposition of `z_excess` (which evaluates the same arithmetic once, at the call core),
  Benjamini-Hochberg-adjusted so it is monotone in call strength. It gates nothing.

  It exists because the calling columns cannot express one real distinction. *Symbiodinium natans* and
  *Phytophthora infestans* are both `NOT_DETECTED` and both fully zeroed, and `z_excess` does not separate
  them (−1.35 vs −0.24) — but Symbiodinium's best of 109 calls is one its background expects 2.4 of (every
  `bg_fdr` = 1), while Phytophthora holds one at 6.9e-3. `P_motif` cannot separate them either: it
  saturates, reading 1.0000 at the top of both.

- **Low-k escape hatch** in the species adjudicator. Genomes with fewer than `cs_min_calls` (3) strong calls
  could not reach `DETECTED` by *any* path — the strength gate is hard-gated on the call count, and
  `z_excess` is bounded below the loss ceiling at low k by construction (k=1 ⇒ z ≤ 0; k=2 ⇒ z ≤ 1). They
  fell through to `NOT_DETECTED`, which zeroes every score. Not adjudicated as a loss: never adjudicated,
  and the evidence discarded — including for genomes carrying a complete minor spliceosome.

  `bg_fdr` is the only statistic here that survives k=1, so such a genome now resolves to **`INCONCLUSIVE`**
  (not `DETECTED`) when its strongest call clears `lowk_bg_fdr_threshold` (3e-3). Scores and calls survive
  for downstream corroboration; no population is asserted. New params `lowk_gate_enabled` /
  `lowk_bg_fdr_threshold`; set the former to `False` to reproduce 3.0.0 behaviour exactly.

  The threshold is a preserve-the-evidence trust line, not a bearer/loss separator — it sits above the
  presumed-loss floor (6.5e-4) on purpose, since the gate abstains rather than asserts. It also sits in the
  only empty band in the low-k distribution, so any value in `[2.65e-3, 6.04e-3]` selects the same set.

  Effect on the 2,785-genome WtMTA v3 corpus: **8 genomes**, all `NOT_DETECTED → INCONCLUSIVE` — four from
  bearer-prior lineages (3 *Mucor*, *Pythium*, *Seison*) and four from loss-prior ones (*Leishmania*,
  *Babesia*, *Porospora*), which is the expected shape for an abstaining gate. Nothing at k ≥ 3 moves.

### Changed

- **`score_info.iic` is now 35 columns (was 38).** `q`, `P_adj`, `P_adj_lo` and `P_adj_hi` were removed —
  all four were deterministic functions of `P_motif` + `motif_category` (`q` a per-species constant;
  `P_adj_lo`/`P_adj_hi` bit-identical to `P_adj`). `adjusted_score` is computed directly as
  `100·q_eff·P_motif`. **`type_id` moved from column 38 to 35**, so any position-based parser must be
  updated; parse by header name.
- `INCONCLUSIVE` now also covers the low-k case above, not only "`z_excess` in the empirical gap". Both mean
  *abstain, corroborate downstream*, but consumers inferring a `z_excess` range from the label will be wrong
  for those genomes.
- `ADJUDICATOR_PARAMS_VERSION` → `zexcess_gap_pgumbel_cs_lowk_2026-07-27`.

### Fixed

- The `apply_pmotif_adjudication` docstring claimed `rel_score = 100*P_motif - 90` while the code computes
  `adjusted_score - 90`. These agree only when `q_eff == 1`, so they diverged completely on `NOT_DETECTED`
  genomes. The code was right: `rel_score` must stay an exact linear transform of `adjusted_score`, because
  that identity is what makes `rel_score > 0` a self-sufficient one-column high-confidence filter.

## [3.0.0] - 2026-07-09

intronIC 3.0.0 is the consolidated, vetted release of the cross-species architecture
developed across the v2 series. It replaces the per-species z-normalization lineage
(v2.5–v2.7) with a raw-feature classifier and a single output-level species adjudicator,
and carries forward — now validated on a 97-species / 14-clade corpus — the improvements
introduced incrementally through v2.2–v2.4. Legacy z/mode-separation bundles remain
reproducible from the `pre-zstack-removal` git tag.

### Relative to the published method (v1; Moyer et al. 2020)

- **Training data.** v1 trained one SVM on human U12-type introns. v3 trains a 42-model
  RBF-SVM ensemble on 41,333 introns from 97 species across 14 clades, labeled by
  comparative genomic analysis across orthologs rather than a single reference genome.
- **Motif model.** The branch-point PWM is built from CoLa-seq empirical branch points
  (Zeng et al. 2022); the 5′SS/3′SS PWMs from a human + comparative-genomics-conserved gold standard. Each
  raw motif score is background-corrected against the species' own U2-type intron pool.
- **Features.** Three motif scores (v1) → six raw features (the three background-corrected
  motif log-odds plus a branch-point offset, a branch-point-scan confidence, and a support
  term), expanded to nine with deterministic interaction terms.
- **Cross-species comparability.** v1's per-species z-normalization under-called divergent
  bearers (e.g. *Amborella*, *Oryza*) and inflated false positives in genomes lacking the
  minor spliceosome. v3 removes it: the ensemble margin is Platt-calibrated to `P_motif`, a
  species-agnostic probability.
- **Genome-level call.** v3 reports whether a genome carries a detectable U12-type
  population (`motif_category`).

### Carried forward from the v2 series (now vetted on the full corpus)

Collected here because most v2.x releases were not published as GitHub Releases:

- CoLa-seq branch-point PWM and comparative-genomics-conserved gold-standard labeling (v2.2.0).
- Species-specific U2-type background correction, including in streaming mode (v2.3.0).
- The `support2`, `bp_offset`, and `bp_scan_confidence` features (v2.2.0–v2.3.0).
- The multispecies training bundle as the default model (v2.4.0).
- Bit-identical `--streaming` / `--in-memory` classification (v2.4.0), now a parity gate.

### Removed — the v2.5–v2.7 z-normalization stack (breaking)

Evaluation on the cross-species corpus showed per-species z-normalization *hurts*
cross-species discrimination — it manufactures false positives in U12-absent species
("z-inflation"; leave-clade-out AUC 0.916 for raw features vs 0.786 for z-normalized).
Removed:

- the adaptive/frozen normalizer and its `MIN_ADAPTIVE_INTRONS` fallback, and the Path C′
  scaler gate (v2.5.0);
- the mode-separation classifier (v2.6.0);
- the continuous per-intron discount (v2.7.0).

Deleted modules: `scoring/normalizer.py`, `scoring/mode_separation.py`,
`scoring/cluster_validation.py`, `scoring/margin_alignment.py`,
`classification/classifier.py`, `classification/mode_sep_pipeline.py`. The
`--load-normalizer` / `--save-normalizer` / `--no-continuous-discount` / `--discount-*`
flags are now accepted-but-ignored no-ops (retained so existing command lines keep
running); legacy z/mode-separation bundles are rejected at load, with a pointer to the
`pre-zstack-removal` tag.

### New — the raw-feature `pmotif_adjudicated` architecture

- **`P_motif`.** Scoring operates directly on the background-corrected raw motif log-odds;
  the per-intron ensemble margin is Platt-calibrated (`σ(2.796·margin − 1.178)`) to
  `P_motif`, a species-agnostic motif probability that can be thresholded post-hoc at any
  confidence level.
- **Species adjudicator** (`scoring/species_adjudicator.py`). `z_excess` (the Poisson
  significance of the strong-call *count* against the genome's own U2-type tail) drives a
  gap gate → `motif_category` ∈ {`DETECTED`, `INCONCLUSIVE`, `NOT_DETECTED`,
  `UNASSESSABLE`} (anchors `loss_ceiling_z = 2.60`, `bearer_floor_z = 5.50`), with a
  per-genome strength gate (`p_gumbel_p95 ≤ 0.01`, `cs_p95 ≥ 5.0` co-fallback) that
  recovers few-but-strong divergent bearers.
- **Calling rule.** `type_id = u12` iff `P_motif ≥ 0.5` (equivalently `adjusted_score ≥
  50`) AND `motif_category ≠ NOT_DETECTED`. Only `NOT_DETECTED` suppresses calls;
  `INCONCLUSIVE` / `UNASSESSABLE` flag species-level ambiguity without vetoing individual
  introns.

### Added

- `score_info.iic` columns: `P_motif`, `z_excess`, `cs_p95`, `p_gumbel_p95`,
  `cs_p95_lo`, `cs_p95_hi`, `motif_category`.
- Per-species **`tail_model.iic.json`** adjudicator sidecar + a per-species tail-model
  diagnostic figure.
- `metrics.iic.json` keys: `motif_category`, `z_excess`, `cs_p95`, `p_gumbel_p95`,
  `motif_called_u12`, `confident_u12_motif`, `normalizer_used`.
- Expanded output-format documentation (see the [Output files](https://github.com/glarue/intronIC/wiki/Output-files) wiki page).

### Figures — per-species diagnostic + scatter/hexbin

- All figure text uses `U2-type`/`U12-type` consistently (tail-model diagnostic, scatter,
  hexbin, 3D, and SVM-diagnostic panels).
- Tail-model diagnostic: the inline "U2 model ≈ …/bin here" callout is replaced by a
  horizontal dashed line at the fitted U2-type tail's expected count where the 95th-pct call
  margin lands, its value labeled on the line over the U2-type distribution (clear of the
  y-axis ticks); the 95th-pct call-margin line is a distinct accent color so it no longer
  shares the U12-type call color.
- Row-locked the scatter / 3D / histogram inputs in **both** the streaming and in-memory plot
  paths so the position array and the tier-score array are built 1:1 (never a positional zip
  of two independently NA-filtered lists), while the density hexbin still keeps every
  valid-(5′,BP) row. Bit-identical on current data (every row carries both); a defensive fix
  against a future scoring-only path emitting position-without-score rows.
- Scatter U2-type gray hexbin: colormap floored at 0.22 (`Greys`→`_U2_DENSITY_CMAP`) so the
  sparse single-count fringe cells lift off the white background instead of mapping to ~white
  on the log scale (they read as "missing" vs the standalone hexbin, though no cell is dropped —
  both use `mincnt=1`). The dense core is unchanged, so the layer still never buries the
  colored U12-type markers.
- Scatter marginal distributions are drawn as smooth filled **KDE traces** instead of bar
  histograms (`Count` → `Density`); huge genomes are stride-subsampled for the KDE and
  degenerate marginals (n<2 / zero spread) are skipped gracefully.
- **Unified color scheme** across the per-species figures: U12-type confidence reads as a
  **green→yellow→red** ramp (high→low) everywhere it is tiered — scatter/3D markers (`>90` green
  Dusty-Olive, `84–90` yellow Honey-Bronze, `≤84` red Oxidized-Iron), tail-model calls (green),
  histogram threshold (green) — and **U2-type is neutral gray** throughout (light-gray bars/fills,
  dark-gray fitted tail / KDE trace / expected-count line). Two palette colors are kept as bright
  tail-model anchor accents: **Baltic-Blue** P=0.9 call line and **Blaze-Orange** 95th-pct call
  margin (q90/exp_max stay gray). The category box maps DETECTED=green / INCONCLUSIVE=yellow /
  NOT_DETECTED=red / UNASSESSABLE=gray. Titles standardized to a neutral dark bold (removes the
  tail-model's old motif_category-colored green title; category still shown in its stats box).
  Density colormaps (inferno hexbin, floored Greys U2-type layer) unchanged.

### Fixed — `-i` / `-d` isoform-flag parity (merged from PR #17 / v2.7.x maintenance)

Both `-i` (`--allow-multiple-isoforms`) and `-d` (`--include-duplicates`) are now
respected end-to-end in **both** `--in-memory` and `--streaming` classify modes, with
streaming/in-memory parity holding across all four flag combinations (regression-guarded
by `tests/integration/test_isoform_flag.py` and the parametrized
`test_streaming_matches_in_memory_with_flags`). Reconciled onto the raw-feature pipeline:

- `extraction/filters.py::should_extract_sequences_for` no longer skips coord-duplicate
  sequence extraction unconditionally — it skips only when `include_duplicates=False`, so
  `-d` duplicates keep their sequences and reach the writer instead of being silently
  dropped as un-scoreable.
- The streaming per-contig worker filter and its worker `config_dict`s now carry and honor
  `allow_multiple_isoforms` (was hardcoded `longest_only=True`), and the streaming
  filtering summary reflects the actual isoform setting.

### Bundled model

- `data/default_pretrained.model.pkl` = the `pmotif_adjudicated` raw-feature ensemble
  (42 sub-models / 210 sub-estimators, RBF `C=200 γ=0.001` isotonic), with adjudicator
  params `ADJUDICATOR_PARAMS_VERSION=zexcess_gap_pgumbel_cs_2026-07-03` (Platt
  `2.796 / −1.178`, anchors `2.60 / 5.50`, strength gate `p_gumbel≤0.01 / cs≥5.0`).

## [2.7.1] - 2026-05-25

### Changed — unified per-intron labels (fixes misleading u2_count)

Per-intron output (`score_info.iic`) now carries three orthogonal label
axes derived after all score adjustments have been applied:

  type_id     ∈ {"u12", "u2"}              binary call on adjusted_score
  confidence  ∈ {"strong", "borderline"}    within-call gradation
  history     ∈ {"stable", "promoted", "demoted"}  pipeline path indicator

Decision logic:
  type_id = "u12" if adjusted_score >= 50 else "u2"
  confidence:
    type_id=u12 → "strong" if adjusted_score >= 90 else "borderline"
    type_id=u2  → "strong" if first_pass_svm < 10 else "borderline"
  history:
    u12 with first_pass_svm < 50  → "promoted"  (mode-sep rescued)
    u2  with first_pass_svm >= 50 → "demoted"   (discount suppressed)
    otherwise                     → "stable"

Thresholds derived empirically from 73-species IPA-gold panel:
`first_pass_svm < 10` captures 99.92% of true U2s with 0.02% U12 loss;
demoted introns have median gap of ~70 score-units (no minimum-gap
filter needed).

`metrics.iic.json` now exposes counts on each axis: `u12_count`,
`u12_strong_count`, `u12_borderline_count`, `u12_promoted_count`,
`u2_count`, `u2_strong_count`, `u2_borderline_count`, `u2_demoted_count`.
`high_confidence_u12` retained as alias for `u12_strong_count`.

**Pre-v2.7.1 behavior superseded**: `u2_count` was previously computed as
`total_classified - high_confidence_u12`, treating uncertain U2s the
same as confidently U2; per-intron `type_id` was set with inconsistent
thresholds across multiple code paths (50% in predictor.py, 90% in
legacy margin alignment and prior adjustment in main.py). All three
assignment paths superseded by the single source-of-truth function in
`scoring/labeling.py`.

Backfill script provided (`scripts/backfill_v271_labels.py`) to add the
new columns to existing score_info.iic files without re-classifying.

### Changed — single-seed ensemble bundle (3× faster classification)

Default model bundle (`default_pretrained.model.pkl`) repacked to a single
seed group of 42 sub-models per pass (down from 3 seeds × 42 = 126).
Empirically validated identical TP/FN/FP/recall across:
- 73-species IPA-gold panel (first-pass): bit-identical at N=42 vs N=126
- 8-species gold panel (second-pass, dedup-fixed): bit-identical
- AmbTri (mode-sep recovery target): 51/51 calls preserved
- Symbio (divergent, gate-fail): 27/27 calls, sigma distributions within 3%
- Drosophila smoke test: 20/20 calls preserved

Runtime: ~3× per-pass speedup measured (e.g., Symbio classification:
1h 53m → 39m 52s at -p 8). Bundle file 82 MB → 27 MB (33%). Theoretical
justification: bagged ensemble variance saturates by N≈40 under 75%-overlap
subsampling (variance floor at ρσ², not 1/N), so additional seeds provide
no measurable boundary refinement. See `docs/scoring_pipeline.md` for
details.

The 3-seed bundles for both passes (v4_aug_cluster_aware,
v5_modesep_aug) are preserved separately and can be reconstituted via
`scripts/rebundle_single_seed.py`.

### Changed — quality tier strings renamed

Mode-separation diagnostic tiers in `<species>.modesep.json` renamed from
single-character school grades to descriptive names:

  "A"  → "modesep_strong"
  "B"  → "modesep_standard"
  "C"  → "modesep_weak"
  "F"  → "first_pass_fallback"

Motivation: the opaque "F" was easily misread as "the run failed" rather
than "gate failed, the pipeline routed through first-pass + legacy
adjustment cleanly". Pre-v2.7.1 `.modesep.json` files retain the old
strings; new runs use the descriptive names.

### Added — scoring pipeline documentation

New `docs/scoring_pipeline.md` with a Mermaid flowchart covering all six
stages (z-scoring → first-pass → quality gate → mode-sep recalibration →
second-pass → continuous discount), the per-intron z_5p eligibility
filter inside the mode-sep branch, the quality-tier rubric, and a
complete parameter reference.

### Fixed — gold-merge inflation bug in evaluation scripts

Five evaluation scripts (`ensemble_size_sweep.py`, `sweep_broad_panel.py`,
`eval_panel_modesep.py`, `eval_panel_v5_1.py`,
`rescore_with_cluster_strong_u2.py`) merged the gold corpus per
alignment_id rather than per intron, inflating TP/FN/FP counts by the
per-intron IPA-partner multiplicity (typically 5-25×). Fixed all five.
Cross-comparison conclusions were unaffected (the inflation applied
symmetrically), but absolute counts in prior TSV outputs should be
considered superseded.


## [2.7.0] - 2026-05-20

### Added — Continuous per-intron discount

Applied to every modesep run (gate-pass and gate-fail), the discount is a
non-positive log-odds penalty:

  penalty_overcall = k_overcall × max(0, svm_vs_naive − τ_overcall)
  penalty_weakmot  = k_weakmot  × max(0, τ_motif       − raw_sum)

where `svm_vs_naive = logit(p_svm) − raw_sum` and
`raw_sum = 5'_raw + bp_raw + 3'_raw`. Both terms are zero in the healthy
regime; they activate only when the SVM overcalls relative to motif log-LR
or when motif evidence is weak.

**Empirically tuned defaults** (against 14-species panel + Salpingoeca):
`k_overcall = 2.0`, `τ_overcall = 0.0`, `k_weakmot = 0.0` (DISABLED by default),
`τ_motif = 10.0`. The overcall penalty alone suppresses Salpingoeca-class
overcalls (29 → 6 pre-legacy) while preserving all 1950 IPA-validated
panel TPs. The weak-motif penalty is opt-in via CLI because it loses
4 IPA-validated borderline TPs (SUDS3, ARPC5, ap4e1, mios) where the SVM
correctly leverages non-motif features (bp_offset, bp_scan_confidence,
support2) for borderline-motif U12s.

### Added — Diagnostic surface

Per-intron columns added to `score_info.iic`:
- `raw_sum` — unweighted motif log-LR sum (5'_raw + bp_raw + 3'_raw)
- `svm_vs_naive` — calibration delta = logit(p_svm) − raw_sum
- `voting_frac` — fraction of ensemble sub-models voting U12 (P > 0.5)

Per-species fields added to `.modesep.json`:
- `boundary_mass` — fraction of eligible introns whose second-pass mean P
  sits in [0.1, 0.9] (OOD-distribution diagnostic)
- `n_called_pre_discount`, `n_called_post_discount`, `continuous_discount_applied`

### Changed — Output semantics

- `svm_score` column unchanged from v2.6 (raw classifier output, preserved
  for auditability)
- `adjusted_score` column is the new recommended call score after the
  continuous discount
- Calling decision: `adjusted_score ≥ threshold` (default 90)

### CLI flags

- `--no-continuous-discount` — disable v2.7 discount (preserves pre-v2.7
  semantics where svm_score is the calling column)
- `--discount-k-overcall` / `--discount-tau-overcall` — tune the overcall
  penalty
- `--discount-k-weakmot` / `--discount-tau-motif` — tune the weak-motif penalty

### Test coverage

- 27 unit tests for `intronIC.scoring.mode_separation` (14 new for v2.7
  covering continuous discount math, voting fraction, boundary mass)
- 9 integration tests for the full classify pipeline (4 new for v2.7
  covering discount preservation of TPs, svm_score invariance under
  discount, and zero-coefficient pass-through)

### Backward compatibility

- v2.6 v5_modesep bundle continues to work; v2.7 discount applies on top
- `--no-continuous-discount` reproduces v2.6 behavior
- v4 cluster-aware bundles continue to load

### HP optimality verification

Before v2.7 release, the inherited SVM hyperparameters (C=200, γ=0.001,
ezf=0.75 — carried over from v4_aug cluster-aware) were empirically
verified against 6 neighbor configurations:

**Subsample sweep** (50K-row corpus, 14-species panel):
- γ=0.0001 marginally edged out control by +4 TPs (1934 vs 1930)
- C/2, 2C, ezf±0.15 statistically tied
- γ=0.01 lost significantly (−24 TPs)

**Full-corpus verification** (502K rows, γ=0.0001 only):
- γ=0.0001 LOST 4 TPs vs control on full corpus (1946 vs 1950)
- Bundle size penalty: 119 MB vs 78 MB (1.5×)
- Subsample → full-corpus generalization failed; wider kernel
  over-smooths fine-structure on dense data

Conclusion: current HPs are at/near the global optimum within the
tested neighborhood. Inherited HPs retained.

## [2.6.0] - 2026-05-17

### Added — Mode-separation classifier
- **New default classifier (`v5_modesep_aug`)** uses per-species mode-
  separation z-scoring to calibrate each species to its own U12/U2
  spectrum: `z = (raw − μ_U2_species) / (μ_U12_species − μ_U2_species)`.
  Collapses the cross-species variation in normalized U12/U2 separation
  (the root cause of plant U12 recall problems in vertebrate-trained
  classifiers).
- **Two-pass classify orchestration**: first pass with the cluster-aware
  ensemble produces candidate weights; mode estimation; gate decision;
  second pass with the recalibrated SVM ensemble. Implemented as a
  post-classify step on `score_info.iic`
  (`intronIC.classification.mode_sep_pipeline`).
- **Gate** combines three checks (`intronIC.scoring.mode_separation.evaluate_gate`):
  - `n_eff_candidates` ≥ 5
  - μ_U12_5p offset from the cross-species anchor (15.671 raw PWM units)
    within ±3.6 (Phase 0 empirical derivation, 61 species, 12 phyla)
  - multi-bandwidth Fisher-discriminant KDE valley depth ≥ 0.30
    (reuses `cluster_validation.compute_valley_depth` via dependency
    injection — robust 3D detector, not a moment proxy)
- **Diagnostic JSON** (`.modesep.json` sidecar) with route, gate reason,
  μ_U2/U12, valley depth, ensemble σ on called introns, and an A/B/C/F
  quality tier per species. Per-intron columns added to score_info.iic:
  `ensemble_sigma`, `first_pass_svm`, `modesep_route`.
- **CLI flags**: `--no-mode-sep`, `--mode-sep-z-floor`,
  `--mode-sep-valley-min`, `--mode-sep-n-floor`,
  `--mode-sep-mu-u12-tolerance`.

### Changed — Bundle schema
- v3 bundles now support mode-sep mode via `config.normalizer_mode = "modesep"`.
  Modesep bundles ship with an EMBEDDED first-pass ensemble
  (`first_pass_seeds` + `first_pass_config`) and a `modesep_params` block
  (z floor, valley min, n_eff floor, universal anchors, location-prior
  tolerance, first-pass model id, PWM set id). `utils/model_io.py`
  surfaces both ensembles at runtime.

### Performance
- Plant recall: AmbTri 90% → 100% (51/51), OrySat 94% → 100% (32/32),
  AraTha 96% → 98% (47/48).
- Apostasia (IPA-validated U12s, held-out from training): 17/21 → 20/21.
- 14-species panel: 1944 → 1950 TP, 9 → 3 FN, FP_strong 0 → 0
  (Pareto improvement over both v2.3 D and cluster-aware v4_aug).
- v2.7 panel net vs v2.6: TP=1950 (unchanged), FN=3 (unchanged),
  FP_strong=0 (unchanged), FP_any=216 (was 217). All IPA-validated TPs
  preserved; per-intron overcall defense added without recall regression.
- U12-absent species (CaeEle, SacCer, SchPom, TetThe, ChlRei, AscSuu)
  gate-fail cleanly to first-pass scores; zero FP_strong leakage.

### Gate-fail defense-in-depth
- When the modesep gate refuses to recalibrate (n_eff floor, location-prior,
  or no-valley), the legacy valley-based log-odds score adjustment
  (`compute_adjusted_scores_batch`) is applied to the first-pass scores.
  This suppresses spurious first-pass calls in noisy or non-bimodal
  species (e.g., Salpingoeca rosetta: 29 first-pass calls → 0 after
  valley-depth=0.035 discount). Matches pre-v2.6 no-valley behavior.

### Backwards compatibility
- v4 cluster-aware bundles still load via `--model <path>`.
- Path C' scaler gate is removed (`scaler_gate.py`, `--no-scaler-gate`).
  v4 bundles now run adaptive-only at classify time; mode-separation
  replaces the scaler-gate's call-asymmetry routing.

## [2.5.0] - 2026-05-14

### Added
- **Path C' species-level scaler gate** (`intronIC.classification.scaler_gate`).
  At classify time, each species is scored under BOTH the per-species
  adaptive RobustScaler and the bundled frozen multispecies scaler.
  Call-set asymmetry between the two modes decides which route the
  species takes: `adaptive` (well-agreed; vertebrate-class), `frozen`
  (adaptive-demotion of borderline U12s; apostasia-class), or `strict`
  (adaptive over-calls into a U12-absent void; tetrahymena / ChlRei-class).
  Threshold defaults (`only_f_thr=0.10`, `only_a_thr=0.50`,
  `n_both_max=5`, `n_total_min=5`) were calibrated on a 45-species
  pressure test (12 dev + 33 validation) with an asymmetric `n_total`
  floor: frozen requires `n_total >= 5`, strict has no floor (so
  ChlRei-class species with `n_total=2, only_a=100%, n_both=0` route
  strict correctly). Gate is ON by default; disable via `--no-scaler-gate`.
- **`--no-scaler-gate` CLI flag** for forcing the pre-2.5.0 adaptive-only
  behavior when desired (debugging, reproducing v2.4.x output).
- **New `score_info.iic` columns**: `svm_score_adaptive`, `svm_score_frozen`,
  `scaler_used` — preserves per-mode scores and the route the gate
  picked, so downstream tools can audit the decision.
- **v4_aug default model** (`default_pretrained.model.pkl`). Same
  corpus, same hyperparameters as the v2.4.x v3 bundle, but trained
  with scaler-strategy=augment: each row appears in two forms
  (adaptive-z and frozen-z) during training, doubling the effective
  corpus to ~1.005M rows. The model learns a single decision boundary
  that is robust to both scaling regimes — which is what makes the
  gate's route-switching safe at inference. v2.4.2 model archived as
  `default_pretrained.model.pkl.bak_v2.4.2`.

### Changed
- **Default-on scaler gate**: Production behavior now differs from
  v2.4.x for species that hit a non-adaptive route. The vast majority
  of well-represented species (vertebrates, drosophila, plants with
  modest U12 counts) route adaptive and see identical behavior to
  v2.4.2. Divergent species (apostasia-class, tetrahymena-class,
  U12-absent species with adaptive-collapse) now route frozen or
  strict and produce materially different call sets.
- **`scaler_used`-aware downstream tools**: scripts that consume
  `score_info.iic` should treat the three new columns as optional;
  they are NA for pre-2.5.0 outputs and for the streaming-classify
  path (where the gate is not yet wired — adaptive-only is enforced
  with a warning).

### Limitations
- The scaler gate is not yet wired into the streaming-classify path
  (`--streaming`). Streaming-mode users see a warning at run start and
  get adaptive-only behavior. For divergent species, use in-memory
  mode. Wiring the gate into streaming requires either a two-pass
  approach or a post-streaming re-emit; planned for v2.5.x.

## [2.4.2] - 2026-05-10

### Added
- **Bundled multispecies fallback scaler** (`v3_fallback_normalizer.pkl`).
  v2.4.0 shipped without a saved scaler in the v3 model bundle, which
  prevented single-intron / tiny-annotation scoring (regression from
  v1 / early v2). The fallback is a frozen RobustScaler fit on 476,848
  raw scores from the 90-species v3 training corpus; loader hooks it
  into the runtime model dict via `_v3_to_runtime`. The classify
  pipeline (both streaming and in-memory) falls through to it when the
  adaptive fit pool is empty or below `MIN_ADAPTIVE_INTRONS`. Validated
  cross-species: chr19 fallback restores valley detection (depth 0.79
  vs adaptive 0.08); Drosophila fallback calls 17 U12s vs adaptive 20
  (gold ~18); Tetrahymena fallback correctly calls 0 (adaptive 8 raw,
  both reach 0 post-valley adjustment). Single-intron scoring via `-q`
  now works again.
- **`--load-normalizer` honored in streaming mode** in addition to
  in-memory. When set, overrides whatever scaler the bundle ships and
  skips the adaptive pre-pass.

### Changed
- **Valley-depth detection now uses Fisher's linear discriminant in 3D**
  (5'z, BPz, 3'z) instead of the naive 2D (5'z + BPz) projection. Fisher's
  direction Σ⁻¹(μ_U12 − μ_U2) — with shrinkage — gives the maximum-
  separation axis for the multi-bandwidth KDE valley search, which makes
  the depth metric more robust on lineages where the 3' motif carries
  meaningful signal (e.g. dinoflagellates, some basal eukaryotes) without
  reducing sensitivity on canonical U12-rich species. No false positives
  introduced on the U12-absent panel.
- **MIN_ADAPTIVE_INTRONS = 200** (was 30 in dev). 30 is just-stable for
  the median; the IQR estimate had ~20% standard-error noise there.
  At 200 IQR noise drops to ~7%. Realistic inputs are almost always
  either single-intron / handful or full-genome, so the band where the
  threshold matters is narrow — but in that band the bundled
  multispecies fallback gives cleaner z-scores than a noisy per-input
  RobustScaler fit.
- **`intronIC test` smoke fixture remains chr19**, now run with
  `--load-normalizer` pointing at the bundled `v3_fallback_normalizer.pkl`.
  This sidesteps the chr19+adaptive valley failure (chr19's narrow
  intron pool compresses the U2 IQR) without needing a larger fixture
  genome — valley fires cleanly at depth ≈ 0.79.

## [2.4.1] - 2026-05-10

### Changed
- **v3 multispecies training corpus is now bundled** as the default
  reference set for `intronIC train`. Files: `u12_reference_multispecies.introns.iic.gz`
  (10,003 U12 introns) and `u2_reference_multispecies.introns.iic.gz`
  (31,330 U2 introns) — the actual 41,333-row training pool the v3
  ensemble was fit on (90 species, 14 clades, post-singleton-decay
  filter). Sequences are bundled with 50 bp flanks per side, oriented
  so the bases adjacent to each splice site are preserved (inside-out
  trim).
- **Legacy v2.3 reference sets renamed** to `u12_reference_human.introns.iic.gz`
  and `u2_reference_human.introns.iic.gz`. Still bundled and selectable
  via `--reference-u12s` / `--reference-u2s` for backwards-compatible
  retraining.

## [2.4.0] - 2026-05-09

### Changed
- **Default bundled model upgraded to v3 multispecies** (3 seeds × 42 = 126
  calibrated SVMs trained on 41,333 introns from 90 species across 14
  clades, with 7 evaluation-only holdout species). Holdout F1 = 1.000 vs
  v2.3 default 0.9975 on the 5-species recall set; production-equivalent
  (valley-adjusted) FPR on U12-absent species is ~54% lower than v2.3
  (12 vs 26 across ~330k scored introns at threshold 90).
- **Default classification threshold lowered from 95 → 90.** The v3
  multispecies model's adjusted-score distribution is tightly calibrated
  (Brier ≈ 4×10⁻⁶) and dropping to 90 picks up real conserved U12s
  without measurable FPR penalty (+1 absolute FP across 330k scored U12-absent
  introns). Pass `--threshold 95` to restore the prior behavior.
- **Streaming classify path refactored** to share its extract+filter+score
  logic with the new adaptive-fit worker (see Added). v2.3-format bundles
  with a saved normalizer continue to behave identically; the refactor is
  covered by the existing `-p 1` vs `-p 3` equivalence tests.
- **`--streaming` and `--in-memory` now produce bit-identical
  classifications** on the same input. Previously the two paths diverged
  on score values (BG correction and KDE valley depth depended on
  upstream iteration order) and on which introns reached `score_info.iic`
  (in-memory included [d]-tagged duplicates and AMBIGUOUS-omitted introns
  with bogus 0.0 raws; streaming did not). Mode choice now affects only
  runtime/memory tradeoff. Locked in by
  `tests/integration/test_streaming_equivalence.py::TestStreamingMatchesInMemory`.

### Added
- **v3 multispecies model bundle support** (`--model <bundle.pkl>`).
  v3 bundles ship the full ensemble plus config and training metadata
  in a single self-describing dict. Loading is automatic —
  `normalize_model_bundle()` detects v3 by the `version` key and
  translates into the runtime shape the rest of the pipeline already
  expects.
- **Adaptive normalizer fitting in streaming mode.** When a model bundle
  has no saved scaler (v3 multispecies default), `intronIC classify
  --streaming` runs a lightweight per-contig pre-pass that scores all
  introns with the (BG-corrected, if enabled) PWMs, then fits a
  `RobustScaler` on the pooled raw 5'/BP/3' distribution. The classify
  pass uses that scaler. Cost: one extra genome pass — comparable to
  the BG-correction pre-pass that already runs by default. Bundles with
  a saved normalizer (v2.3 format) skip this step.
- **`IntronScores.support2` derived feature** — second-largest of the
  clipped-at-zero z-scores (5', BP, 3'). Used by v3 multispecies as a
  6th feature; a `@property` on the dataclass so it stays derivable
  from the existing z-scores at zero memory cost.
- `ScoreNormalizer.fit_from_array()` — fit directly from a raw-score
  matrix (skips the intron-object iteration step). Used by the
  streaming adaptive-fit pre-pass.
- `apply_scaler_to_scored_batch()` — splits the scale-then-classify
  step out of `score_and_normalize_batch()` so the streaming worker can
  score once and apply a scaler that's resolved later in the same run.
- v3 bundle schema spec at `docs/v3_bundle_schema.md`.

### Fixed
- `SpeciesBackground.build_corrected_pwm_sets` gated on `n_introns`
  (`len(_intron_seqs)`), which is only populated by the in-memory
  `accumulate()` path. Streaming uses `merge_worker_counts()` and left
  `_intron_seqs` empty, so the `min_introns` guard always tripped and
  the streaming path silently returned the human U2 PWMs unchanged
  for every run with BG correction enabled. Switched to `n_accumulated`
  (which has the right fallback over `_RegionAccumulator` counts).
- `_build_final_pwm_sets` iterated `acc.get_all_subtypes()` in dict
  insertion order. Multiple non-canonical 5' dnts collapse to the
  same `('u2', 'gtag')` slot via the `FIVE_DNT_TO_SUBTYPE` fallback,
  and last-writer-wins made the final PWM order-dependent. Sort
  subtypes before iteration in both the 5'/3' and BPS loops.
- `compute_valley_depth` subsampled `u2_proj` with `np.random.choice`
  under a fixed seed, so the indices were deterministic but the points
  at those indices inherited the upstream iteration order. Sort
  `u2_proj` (and `u12_proj` for symmetry) before sampling so the KDE
  input is data-determined, not order-determined.
- In-memory `score_introns` ran on omitted introns and produced bogus
  0.0 raw scores plus arbitrary z-scores for them; streaming was
  already filtering on `omitted == NONE`. Filter the in-memory
  `introns_for_scoring` list to match. Omitted introns still flow
  through `merge_scored_and_omitted_introns` into `meta.iic` /
  `bed.iic` with their omitted reason.
- In-memory `apply_species_background`, `classify_with_pretrained_model`,
  and `write_outputs` now uniformly skip coordinate-duplicate isoforms
  when `include_duplicates=False`, matching streaming. In-memory
  `score_info.iic` no longer includes `[d]`-tagged duplicate rows.

### Notes
- Existing v2.3 model bundles continue to load unchanged (the legacy
  `{"ensemble": ..., "normalizer": ...}` dict is pass-through).

## [2.3.0] - 2026-04-23

### Changed
- **Default model upgraded from 8D to 6D RBF SVM** with 42-model ensemble (was 15-model, 8D)
  - Kernel: RBF (C=35.11, gamma=0.01), isotonic calibration
  - 6 features: s5z, BPz, s3z, support2, bp_offset, bp_scan_confidence
  - Dropped ppt_t_weighted and ppt_longest_run (minimal discriminative value in RBF kernel)
  - U2 subsample ratio: 75% per model (was 80%)
- **Default threshold raised from 90% to 95%** for higher-confidence U12 calls
- `rel_score` now computed as `adjusted_score - threshold` (was `svm_score - threshold`)

### Added
- **Bayesian score adjustment** combining two independent signals:
  - Species-level valley prior from 2D KDE bimodality (5'z, BPz) — discounts probabilities in species lacking distinct U12 clusters
  - Per-intron ensemble sigma penalty — penalizes introns with high model disagreement
  - Formula: `logit(p_adj) = logit(p_svm) + log(pi_species / 0.5) - k_sigma * sigma`
  - Configurable via `scoring.score_adjustment` in YAML config
- **Species-specific U2 background correction** now works in streaming classification mode (was disabled, forced fallback to non-streaming)
- New `score_info.iic` columns: `adjusted_score` (post-adjustment probability) and `ensemble_sigma` (per-intron model agreement)
- Valley depth and ensemble sigma reported in metrics JSON output
- Shared background correction helpers (`_create_background_accumulator`, `_finalize_background_correction`) eliminate three-way code duplication
- `ScoreAdjustmentConfig` dataclass for score adjustment parameters

### Fixed
- Streaming classify mode now supports species-specific background correction (previously gated off)
- Metrics JSON `high_confidence_u12` count now reflects adjusted scores, not raw SVM scores

## [2.2.0] - 2026-04-11

### Changed
- **Default model upgraded from 4D linear SVM to 8D RBF SVM**
  - Kernel: RBF (C=65.88, gamma=0.01), 15-model ensemble with isotonic calibration
  - 8 features: s5z, BPz, s3z, support2, bp_offset, ppt_t_weighted, ppt_longest_run, bp_scan_confidence
  - Training data expanded to 472 U12 + 30,155 U2 reference introns (was 387 U12 + 20,690 U2)
  - U12 gold standard derived from IPA-conserved introns with CoLa-seq empirical branch points
  - U2 training set merges original human introns with IPA-conserved hard negatives
  - Cross-species false positive reduction: 0 in C. elegans (was 2), 1 in Ascaris (was 47)
- 3'SS scoring window narrowed to core-only (-6 to +3) to cleanly separate PPT signal from acceptor motif
- BPS search region no longer overlaps with 3'SS scoring region (was 1bp overlap)

### Added
- New scoring features computed for all introns and written to score_info output:
  - `support2`: second-largest positive z-score — encodes "at least two sites support U12"
  - `bp_offset`: branch point adenosine position relative to 3'SS (definitional, using PWM reference position)
  - `ppt_t_weighted`: T-weighted pyrimidine score (T=1.0, C=0.5, purine=0) over fixed 20nt window anchored at 3'SS boundary
  - `ppt_longest_run`: longest uninterrupted C/T run in the same 20nt window
  - `bp_scan_confidence`: BPS motif sharpness — log2(best/mean) of U12 PWM scores across the BP search window
  - Per-site absolute fit scores: `fit_u12_5'`, `fit_u12_bp`, `fit_u12_3'`, `min_fit_bp_3`
- New composite features in BothEndsStrongTransformer: `support2`, `support_rest`, `concentration`
- `reference_offset` field on PWM for definitional branch point adenosine position
- Extra features pipeline: `extra_feature_names` parameter threaded through classifier, optimizer, trainer, predictor, and transformer for configurable feature sets beyond the 3 base z-scores

### Fixed
- BED input mode crash when introns lack metadata (missing `IntronMetadata` initialization)
- BP search region properly excluded from 3'SS scoring region boundary

## [2.1.3] - 2026-03-10

### Added
- New `NO_SEQUENCE` omission reason (`[x]` tag) for introns on genome regions missing from the FASTA (e.g. organellar genes referencing contigs not included in nuclear assemblies)

### Fixed
- Pipeline no longer crashes when introns reference contigs absent from the genome FASTA
  - `SequenceExtractor` now yields these introns marked as `omitted_no_sequence` instead of silently dropping them
  - `SequenceWriter` skips introns without sequence data instead of raising `ValueError`
  - Affected species include those with mitochondrial/plastid gene annotations (e.g. *Manihot esculenta*, *Protopterus annectens*, *Spinacia oleracea*)

### Dependencies
- Bump pillow from 12.0.0 to 12.1.1

## [2.1.2] - 2026-02-15

### Fixed
- Log files no longer contain ANSI color codes
- Disabled line wrapping in log file output

## [2.1.0] - 2024-12-15

### Added
- New `intronIC test` subcommand for quick installation validation
  - Runs classification on bundled Chr19 test data
  - Displays intron counts and classification results
  - `--show-only` flag to show test data location without running
  - `--output-dir` flag to specify output location

### Changed
- **Default model switched to isotonic calibration** for improved cross-species performance
  - Better classification accuracy in C. elegans and other non-human species
  - Previous sigmoid-calibrated model preserved as `default_pretrained.model.sigmoid.pkl`
- Removed Python version upper bounds (now compatible with Python 3.14+)
- Removed numpy <2.0 constraint (verified compatibility in commit 0525903)

### Fixed
- `intronIC test` now displays correct intron counts from metrics file
  - Reads from `.metrics.iic.json` instead of looking for non-existent JSON
  - Shows `total_scored` (classified introns) instead of `total_introns_generated`
  - Shows `high_confidence_u12` to match classification results table
  - Added thousand separators for better readability
- Fixed UnboundLocalError in test command by initializing variables before loop
- Completely suppressed sklearn InconsistentVersionWarning
  - All `joblib.load()` calls replaced with `load_model()` wrapper
  - Warning filter now catches all sklearn warnings, not just UserWarning
  - Applies to pretrained model, normalizer, and ensemble loading paths
  - Suppresses version mismatch warnings for LinearSVC across sklearn 1.7.2 → 1.8.0

## [2.0.10] - 2024-12-15

### Changed
- Default model switched to isotonic calibration
- Previous sigmoid model saved as `default_pretrained.model.sigmoid.pkl`

### Fixed
- Sklearn version warnings suppressed in model loading

## [2.0.9] - 2024-12-15

### Fixed
- `intronIC test` command now reads correct metrics from `.metrics.iic.json`
- Displays accurate intron counts with thousand separators

## [2.0.8] - 2024-12-15

### Fixed
- UnboundLocalError in test command variable initialization
- Initial attempt at suppressing sklearn warnings

## [2.0.7] - 2024-12-14

### Added
- `intronIC test` subcommand for installation testing
- Bundled Chr19 test data in PyPI package (~20MB)

### Changed
- Updated README with new test command

## [2.0.6] - 2024-12-08

### Changed
- Remove Python version upper bounds
- Remove numpy <2.0 constraint
- Update pixi.toml to match pyproject.toml constraints

## [2.0.0] - 2024-11-26

### Changed
- Complete rewrite with modular architecture
- src/ layout for better packaging isolation
- Streaming mode for memory-efficient processing
- Domain adaptation with adaptive/frozen normalizer modes
- Improved CLI with classify/train/extract subcommands

### Added
- Per-contig streaming classification
- JSON metrics output (.metrics.iic.json)
- Rich console output with progress tracking
- Comprehensive logging system
- Model metadata tracking

### Fixed
- Memory optimization for large genomes
- Coordinate system handling
- Feature extraction pipeline

## [1.5.1] - Earlier

Previous monolithic version. See git history for details.

---

[2.3.0]: https://github.com/glarue/intronIC/compare/v2.2.0...v2.3.0
[2.2.0]: https://github.com/glarue/intronIC/compare/v2.1.3...v2.2.0
[2.1.3]: https://github.com/glarue/intronIC/compare/v2.1.2...v2.1.3
[2.1.2]: https://github.com/glarue/intronIC/compare/v2.1.0...v2.1.2
[2.1.0]: https://github.com/glarue/intronIC/compare/v2.0.10...v2.1.0
[2.0.10]: https://github.com/glarue/intronIC/compare/v2.0.9...v2.0.10
[2.0.9]: https://github.com/glarue/intronIC/compare/v2.0.8...v2.0.9
[2.0.8]: https://github.com/glarue/intronIC/compare/v2.0.7...v2.0.8
[2.0.7]: https://github.com/glarue/intronIC/compare/v2.0.6...v2.0.7
[2.0.6]: https://github.com/glarue/intronIC/compare/v2.0.0...v2.0.6
[2.0.0]: https://github.com/glarue/intronIC/compare/v1.5.1...v2.0.0
