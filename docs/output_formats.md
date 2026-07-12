# intronIC v3 Output File Formats — authoritative spec

**Applies to:** intronIC **v3.0.0** (`pmotif_adjudicated` raw-feature architecture).
**Scope:** every file the **classify** pipeline writes for one species. Each schema below
was verified against the writer source (`file:line` cited) and cross-checked against a live
production run (`bos_taurus-ARS-UCD1.2`). Where a docstring disagreed with the code, the
**code is authoritative** and the discrepancy is flagged.

> Companion: `docs/generated_files_manifest.tsv` (machine-readable index of every file).
> Pipeline narrative: `docs/scoring_pipeline.md`. Model bundle: `docs/v3_bundle_schema.md`.
> Adjudicator: `docs/adjudicator_design.md`.

---

## 0. Conventions

- **Run name / key.** All files are prefixed by the `-n` run name (`<base_filename>`). In the
  WtMTA corpus that is `<key>.full`, where `<key>` = the species directory basename. Files land
  in the `-o` output directory (in the corpus: `v3_runs/<key>/`).
- **Two execution paths, one result.** Production is always **streaming**
  (`classify_streaming_per_contig`, `cli/main.py:3740`); an in-memory path exists and produces
  **bit-identical** per-intron files (enforced by `tests/integration/test_streaming_equivalence.py`).
  The mapping files (`.dupe_map.iic`/`.overlap_map.iic`) are written **only** by the streaming path.
- **Delimiter** is a literal **TAB** in every `.iic` table; row terminator `\n`.
- **Two missing-value sentinels** (this matters for parsing):
  - Columns written by the base writers use the string **`NA`**.
  - Numeric columns appended by the adjudicator use lowercase **`nan`** (`%.6g` formatting).
- **Scored-only vs full inventory.** `score_info.iic` contains **only scored introns**
  (`ScoreWriter.write_intron` early-returns on omitted introns, `writers.py:1272-1274`).
  `meta.iic`, `bed.iic`, and `introns.iic` contain the **full inventory** including omitted
  introns (rows whose `type_id`/scores are `NA`). Join all files on the **`name`** column.
- **Authoritative call source.** `score_info.iic` is the source of truth for the per-intron
  call; `meta.iic` and `bed.iic` are **synced from it** post-adjudication
  (`_sync_calls_to_meta_and_bed`, `main.py:311-334`). Do **not** derive the U12/U2 call from
  `introns.iic`'s numeric score (it is the pre-adjudication `svm_score`).

---

## 1. The call, in one place (read this before using any file)

A per-intron call and a per-species population verdict are **two different things**:

| Question | Field | File(s) | Rule |
|---|---|---|---|
| Is this intron U12 or U2? | **`type_id`** | score_info (auth.), meta | `u12` iff **`P_motif ≥ 0.5` AND `motif_category ≠ NOT_DETECTED`**, else `u2` |
| How U12-like is the motif? | **`P_motif`** ∈ [0,1] | score_info | species-agnostic, Platt-calibrated motif probability (col 27) |
| High-confidence U12? | `rel_score > 0` | score_info, meta | ≡ `adjusted_score > 90` ≡ `P_motif > 0.9` (strict; `rel_score` is centered on the 90% threshold) |
| Does this genome have a U12 population? | **`motif_category`** | metrics, score_info, tail_model | `DETECTED` / `INCONCLUSIVE` / `NOT_DETECTED` / `UNASSESSABLE` (per-species, same on every row) |

`motif_category` semantics: `DETECTED` = a U12-motif population is present (corroborate
downstream — **not** a confirmed bearer). `NOT_DETECTED` = calls consistent with the U2
background; it is the **only** category that suppresses calls (forces `u12_count → 0`).
`INCONCLUSIVE`/`UNASSESSABLE` flag species-level ambiguity but **still let** strong-motif
introns call, for downstream corroboration.

**`rel_score` sign convention (important):** signed and centered on the **90%** U12-likelihood
threshold, range **[−90, +10]**. `rel_score < 0` means "below the 90% HC threshold", **not**
"anti-U12". A U2 call and a sub-HC U12 call can both be negative.

---

## 2. `score_info.iic` — the authoritative per-intron matrix (38 columns)

Tab-delimited, **one header row**, uncompressed, **scored introns only**. Written in **two
stages**: `ScoreWriter` (`writers.py:1211-1401`) writes the header + **columns 1-26**; then
`apply_pmotif_adjudication` (`scoring/species_adjudicator.py:615-730`) re-reads the file with
pandas, appends **columns 27-38**, and **overwrites** `rel_score` (col 2) and `adjusted_score`
(col 25). Cols 1-26 use `NA`; adjudicator cols use `nan`.

Scope key: **[I]** = per-intron; **[S]** = per-species constant (identical on every row — store
once per species if normalizing).

| # | Column | Meaning | Type / range | NA | Scope |
|---|--------|---------|--------------|----|----|
| 1 | `name` | Intron ID `{SpAbbr}-{gene}@{transcript}_{idx}(family)[tags]` (see §8) | string | — | I |
| 2 | `rel_score` | *Adjudicator-overwritten* `= adjusted_score − 90`. `>0` ⇔ HC U12 | float [−90,+10] | nan | I |
| 3 | `svm_score` | **First-pass** ensemble U12 prob ×100 (isotonic). **NOT the call column** | float 0-100 | NA | I |
| 4 | `5'_seq` | 5′SS extraction-window sequence | DNA | NA | I |
| 5 | `5'_raw` | 5′SS raw motif log-odds (U12/U2), background-corrected | float, signed | NA | I |
| 6 | `bp_seq` | Branch-point seq under the U12 BPS PWM | DNA | NA | I |
| 7 | `bp_seq_u2` | Branch-point seq under the U2 BPS model | DNA | NA | I |
| 8 | `bp_raw` | Branch-point raw motif log-odds, bg-corrected | float, signed | NA | I |
| 9 | `3'_seq` | 3′SS extraction-window sequence | DNA | NA | I |
| 10 | `3'_raw` | 3′SS raw motif log-odds, bg-corrected | float, signed | NA | I |
| 11 | `decision_dist` | First-pass SVM decision distance; `>0` ⇒ first-pass u12 | float, signed | NA | I |
| 12 | `bp_offset` | Distance (nt) from branch point to 3′SS | int (negative) | NA | I |
| 13 | `bp_scan_confidence` | Confidence of the BPS scan pick | float | NA | I |
| 14 | `ppt_ct` | Polypyrimidine-tract score (fraction pyrimidine) | float ~0-1 | NA | I |
| 15 | `ppt_raw` | PPT raw log-odds | float, signed | NA | I |
| 16 | `core_3'_raw` | Core-3′ raw log-odds. **`NA` on every row in the shipped bundle** (unpopulated) | float/NA | NA | I |
| 17 | `fit_u12` | Absolute goodness-of-fit to the U12 model | float, signed | NA | I |
| 18 | `fit_u2` | Absolute goodness-of-fit to the U2 model | float, signed | NA | I |
| 19 | `fit_u12_5'` | U12 fit, 5′SS component | float | NA | I |
| 20 | `fit_u12_bp` | U12 fit, branch-point component | float | NA | I |
| 21 | `fit_u12_3'` | U12 fit, 3′SS component | float | NA | I |
| 22 | `min_fit_bp_3` | `min(fit_u12_bp, fit_u12_3')` (robust conjunctive fit) | float | NA | I |
| 23 | `ppt_longest_run` | Longest pyrimidine run in the PPT window | int | NA | I |
| 24 | `ppt_t_weighted` | T-weighted PPT score | float | NA | I |
| 25 | `adjusted_score` | *Adjudicator-overwritten* `= 100·P_adj`. The 0-100 calling scale | float 0-100 | nan | I |
| 26 | `ensemble_sigma` | SD of per-intron prob across the 42 ensemble models (uncertainty) | float ≥0 | NA | I |
| 27 | **`P_motif`** | **Authoritative motif probability** `σ(2.796·margin − 1.178)`. Species-agnostic | float [0,1] | nan | I |
| 28 | `z_excess` | **Primary species stat**: Poisson significance of strong-call count vs the genome's U2 tail | float ≥0 | nan | **S** |
| 29 | `cs_p95` | 95th pct of un-clipped call margins (call strength); strength-gate co-fallback | float | nan | **S** |
| 30 | `p_gumbel_p95` | **Primary strength driver**: per-genome Gumbel tail-prob at `cs_p95`; `≤0.01` ⇒ DETECTED | float [0,1] | nan | **S** |
| 31 | `cs_p95_lo` | Bootstrap CI lower bound of `cs_p95` — **annotation, not a gate** | float | nan | **S** |
| 32 | `cs_p95_hi` | Bootstrap CI upper bound of `cs_p95` — annotation | float | nan | **S** |
| 33 | **`motif_category`** | **Authoritative species gate** (enum, never nan) | DETECTED/INCONCLUSIVE/NOT_DETECTED/UNASSESSABLE | — | **S** |
| 34 | `q` | **LEGACY/vestigial** back-compat gate ∈ {0,1}: 0 iff NOT_DETECTED | float {0,1} | nan | **S** |
| 35 | `P_adj` | **LEGACY** `= q·P_motif`. Feeds adjusted_score/rel_score/type_id (still load-bearing) | float [0,1] | nan | I |
| 36 | `P_adj_lo` | **LEGACY** `= q_lo·P_motif` (CI collapsed to point) | float [0,1] | nan | I |
| 37 | `P_adj_hi` | **LEGACY** `= q_hi·P_motif` (== P_adj) | float [0,1] | nan | I |
| 38 | **`type_id`** | **Resolved call** `u12`/`u2` (never NA here) | string | — | I |

Column origins: cols 2-26 first-pass values `writers.py:1276-1358`; header `writers.py:1215-1242`;
adjudicator overwrites/appends `species_adjudicator.py:685-704`; enum `species_adjudicator.py:126-135`.

**DB guidance.** Key off **`P_motif` (27) + `motif_category` (33) + `type_id` (38)**. Columns
**34-37 (`q`, `P_adj`, `P_adj_lo/hi`) are explicitly legacy** (`species_adjudicator.py:634-636,695`;
`q` is degenerate 0/1, the `P_adj` CI has collapsed to a point — species uncertainty now lives in
`motif_category`, not a CI). Note the subtlety: `adjusted_score`/`rel_score`/`type_id` are *derived
from* `P_adj`, so the chain is still functionally load-bearing even though the four raw columns are
flagged vestigial. `svm_score` (3) is a **different calibration** than `P_motif` and is **not** the
call — don't threshold it.

---

## 3. `meta.iic` — per-intron metadata (16 columns)

Tab-delimited, **has a header**, full inventory (omitted rows carry `NA`). Written by `MetaWriter`
(`writers.py:835-1014`; header `:861-878`, row `:978-995`). `type_id`/`rel_score` are **synced from
`score_info.iic`** after adjudication.

| # | Column | Meaning | Type | NA |
|---|--------|---------|------|----|
| 1 | `name` | Same intron ID as score_info col 1 | string | — |
| 2 | `rel_score` | Signed, 90%-centered [−90,+10]; `>0` = HC U12 (synced from score_info) | float 4dp | NA |
| 3 | `dnts` | Terminal dinucleotides, e.g. `GT-AG`, `AT-AC` | string | NA |
| 4 | `motif_schematic` | `{exon[-3:]}\|{5'seq}...{bp}/{bp_u2}...{3'seq}\|{exon[:3]}` (`writers.py:586-644`) | string | NA |
| 5 | `bp_context` | BP-region seq with branch point wrapped `[...]` + 3′ display (`writers.py:670-713`) | string | NA |
| 6 | `bp_offset` | Branch-point → 3′SS distance (nt, negative) | int | NA |
| 7 | `length` | Intron length (nt) — **never NA** | int | — |
| 8 | `parent` | Transcript ID | string | NA |
| 9 | `grandparent` | Gene ID | string | NA |
| 10 | `index` | Intron ordinal within transcript | int | NA |
| 11 | `family_size` | Introns in transcript | int | NA |
| 12 | `frac_pos` | Fractional position along the gene | float 4dp | NA |
| 13 | `phase` | Coding phase {0,1,2} | int | NA |
| 14 | `type_id` | `u12`/`u2` (synced from score_info); `NA` for omitted | string | NA |
| 15 | `feature` | Which feature defined the intron: `cds` or `exon` (**not** the U12/U2 type) | string | NA |
| 16 | `attributes` | Comma-joined flags (`noncanonical`, `not_longest_isoform`, `omitted_*`, …) | string | NA |

*Docstring discrepancy:* the class docstring lists 14 columns (omits `bp_offset`, `attributes`) and
says `type_id` writes `'.'`; the code writes **16 columns** with `NA`. Code is authoritative.

---

## 4. `bed.iic` — genomic coordinates (7 columns, no header)

Written by `BEDWriter` (`writers.py:721-827`). **No header.** Full inventory.

| # | Field | Meaning | NA |
|---|-------|---------|----|
| 1 | chrom | Contig/sequence ID (e.g. `NC_037328.1`) | — |
| 2 | start | **0-based** BED start (`intron.start − 1`, `writers.py:774`) | — |
| 3 | stop | End (1-based inclusive == 0-based exclusive → standard BED half-open) | — |
| 4 | name | Same intron ID as score_info/meta col `name` | — |
| 5 | score | **`adjusted_score`** (adjudicated 0-100 calling scale, synced from score_info post-adjudication, `main.py:336`). The 50-split **is** the call (`≥ 50` ⟺ `u12`); `NA` if unscored | NA |
| 6 | strand | `+` / `-` (genomic; coords are never strand-flipped) | — |
| 7 | attributes | Same flag list as meta col 16 | NA |

*Docstring discrepancy:* docstring says 6 columns and `'.'` for missing; code writes **7 columns**
(adds `attributes`) with `NA`. Code is authoritative.

---

## 5. `introns.iic` — sequences (5 columns, no header)

Written by `SequenceWriter` (`writers.py:1022-1173`). **No header.** Full inventory (every intron
that has sequence — a superset of score_info). **Shipped gzipped as `.introns.iic.gz`** — but the
compression is done by the **orchestrator**, not intronIC (intronIC writes plain `.introns.iic`).

| # | Field | Meaning |
|---|-------|---------|
| 1 | name | Same intron ID |
| 2 | score | `svm_score` rounded 2dp (first-pass, 0-100), `NA` if unscored — **not** the adjudicated call |
| 3 | upstream_flank | Upstream exonic flank (length = `config.extraction.flank_len`; production run used **50**) |
| 4 | seq | Intron sequence |
| 5 | downstream_flank | Downstream exonic flank |

This is the input format for `-q`/`--sequence-file` re-scoring mode.

---

## 6. `metrics.iic.json` — per-species aggregate (26 keys)

Assembled by `_finalize_classification_metrics` (`cli/main.py:544-652`, the single source of truth
for both paths); counts tallied post-adjudication from the finalized `meta.iic`
(`count_calls_from_meta`, `summarize_boundaries_from_meta`, `writers.py:51-162`). **All 26 keys
always present**; adjudicator stats become JSON `null` when the U2 tail can't be referenced.

| Key | Type | Meaning |
|---|---|---|
| `total_introns` | int | Introns **written** to output (superset; includes unscored/omitted `type_id==NA` rows) |
| `threshold` | float | High-confidence / `rel_score`-centering threshold (90.0) — **not** the call boundary; the U12 **call** is `adjusted_score ≥ 50` (§1) |
| `total_genes` | int | Genes processed |
| `total_introns_generated` | int | All introns generated from CDS across isoforms, pre-dedup/filter (largest total) |
| `total_scored` | int | Introns scored+classified (= `u12_count + u2_count`) |
| `pretrained` | bool | Always true (classify path) |
| `model_path` | str | Absolute path of the loaded model `.pkl` |
| `streaming_mode` | str | `per_contig` or `in_memory` (only intentional cross-path difference) |
| `normalizer_used` | str/null | `"none (raw features)"` (legacy/no-op field) |
| `u12_count` | int | Scored U12 calls (`type_id==u12`; P_motif≥0.5 ∧ ≠NOT_DETECTED) |
| `u2_count` | int | Scored U2 calls |
| `high_confidence_u12` | int | Strong U12 subset (`rel_score>0` ≡ P_motif≥0.9); ≤ u12_count |
| `u12_percentage` | float | `u12_count/total_scored·100` |
| `high_confidence_percentage` | float | `high_confidence_u12/total_scored·100` |
| `u12_boundaries` | dict | Terminal-dinucleotide tally over U12 calls (top-20). **`AT-AC` = the diagnostic minor subtype** |
| `u2_boundaries` | dict | Same over U2 calls |
| `high_confidence_u12_boundaries` | dict | Dinucleotide tally over the HC U12 subset |
| `high_confidence_u12_by_feature` | dict | HC U12 split by `cds`/`exon`; sums to `high_confidence_u12` |
| `motif_category` | str/null | Species gate (see §1) |
| `z_excess` | float/null | Primary population statistic |
| `cs_p95` | float/null | Call-strength (95th pct of margins) |
| `p_gumbel_p95` | float/null | Strength-gate driver (Gumbel outlier prob) |
| `cs_p95_lo` | float/null | CI annotation on cs_p95 |
| `motif_called_u12` | int/null | **Ungated** count (P_motif≥0.5 ignoring motif_category); == u12_count unless NOT_DETECTED |
| `feature_type` | str | Extraction features: `cds`/`exon`/`both` |
| `confident_u12_motif` | int | `high_confidence_u12` iff DETECTED else 0 (the "safe to headline" HC count) |

Totals nest: `total_introns_generated ⊋ total_introns ⊋ total_scored`; `u12_count + u2_count ==
total_scored`; `high_confidence_u12 ⊆ u12_count`. **Zero-AT-AC flag:** a `DETECTED` genome with many
U12s but `u12_boundaries['AT-AC'] == 0` is the hallmark of an **incomplete annotation** (short minor
introns dropped by the gene models) — track it downstream.

---

## 7. `tail_model.iic.json` — adjudicator sidecar

Assembled by `build_tail_model` (`species_adjudicator.py:556-612`), written at `:709-713`. Only
emitted by `pmotif_adjudicated` bundles. **Always-present** keys: `params_version`, `platt_a`,
`platt_c`, `call_threshold` (0.90), `u2_threshold` (0.50), `evt_tail_pct` (90.0), `evt_tail_frac`
(0.10), `n_u2`, `n_call` (count of `P_motif ≥ 0.90`, **ungated** — equals metrics `high_confidence_u12`
except on `NOT_DETECTED` genomes, where HC is suppressed to 0 but `n_call` is not), `z_excess`, `cs_p95`, `p_gumbel_p95`,
`motif_category`, `call_margins` (the per-call ensemble margins — the strength signal Platt saturation
discards), `assessable`, `reason`. **Conditional** (only when `n_u2 ≥ 200`): `u2_hist_edges`,
`u2_hist_counts`, `med`, `mad`, `q_tail`, `q90`, `lam`, `exp_max` (the fitted U2 EVT tail).

Gate logic (`classify_motif_category`, `species_adjudicator.py:138-168`), evaluated in order:
non-finite `z_excess` → `UNASSESSABLE`; `z_excess ≥ 5.50` → `DETECTED`; strength gate
(`n_calls ≥ 3` ∧ (`p_gumbel_p95 ≤ 0.01` **OR** `cs_p95 ≥ 5.0`)) → `DETECTED`; `z_excess ≤ 2.60` →
`NOT_DETECTED`; else `INCONCLUSIVE`. Anchors: `loss_ceiling_z=2.60`, `bearer_floor_z=5.50`.

---

## 8. Mapping files, log, plots

- **`dupe_map.iic`** / **`overlap_map.iic`** — `MappingWriter` (`writers.py:1432-1466`). Two columns,
  no header: `representative_intron_id  \t  excluded_intron_id`, one row per excluded member. **IDs
  are raw `intron_id`s** (`rna-XM_…:intron_start_end`), **not** the formatted `name`. `dupe_map`
  written iff duplicates exist (normally present); `overlap_map` iff overlaps exist (normally absent).
  Streaming path only (`main.py:4304-4323`).
- **`iic.log`** (`<key>.full.iic.log`) — freeform human-readable text (Rich box tables), the logger's
  only sink (`setup_logging`, `main.py:1498-1542`). Captures the command+config banner, stage
  messages, the `pmotif_adjudicated:` adjudication summary line, and three tables (Intron Filtering
  Summary, Classification Results, Top-20 boundaries). **Not machine-parseable** — use the JSON/TSV
  files for data.
- **Plots** (`visualization/plots.py`), plotted in background-corrected raw-motif-log-odds space:
  - Always: `.plot.hex.iic.png` (5′×BP density), `.plot.scatter.iic.png` (U12 tiers + marginals).
  - Conditional: `.plot.scatter3d.iic.png` (5′×BP×3′, iff 3′ scores), `.plot.score_histogram.iic.png`
    (calling-score histogram), `.plot.tail_model.iic.png` (adjudicator diagnostic, iff sidecar + U2
    population), `.plot.scatter_raw.iic.png` (raw-motif companion scatter, iff `motif_category ==
    NOT_DETECTED` — colours introns by per-intron `100·P_motif` so the motif-strong introns the species
    call suppressed stay visible; subtitle keeps the true `0 U12-type` count).
  - *Training-only* suffixes (`.AUC.iic.png`, `.ref_hex.iic.png`, `.plot.training_scatter.iic.png`,
    `.plot.decision_surface.iic.png`) are **never** emitted by a classify run.

### Name string format (`generate_intron_name`, `writers.py:395-498`)
`{SpAbbr}-{grandparent}@{parent}_{index}({family_size}){tags}` — e.g.
`BosTau-gene-CRYZL1@rna-NM_001035037.3_11(12)`: `BosTau` = 3+3 species abbreviation; `gene-CRYZL1`
= gene (grandparent); `rna-NM_001035037.3` = transcript (parent); `_11` = intron index; `(12)` =
introns in transcript. Optional trailing `;[...]` tags encode omission reason / dynamic flags.

---

## 9. Not intronIC's own outputs (external orchestrator)

These appear in a corpus run dir but are written by the **WtMTA orchestrator**, not intronIC — a
consumer should not attribute them to the intronIC spec:

- **`.introns.iic.gz`** — intronIC writes plain `.introns.iic`; the orchestrator gzips it post-run.
- **`.runner.log`** — the orchestrator's per-task log (distinct from intronIC's `.iic.log`).
- **`*.gene_contig_count`** — a genome-dir sidecar the orchestrator caches for its mode heuristic.
- **`.modesep.json`** — legacy; only written if mode-separation runs, which it never does in the
  deployed pmotif bundle.
