# raw_gated scoring — progress & plan

Single source of truth for the **raw-score + continuous-gate** scoring architecture: a candidate replacement
for the z-normalizer + mode-separation + continuous-discount stack. This doc tracks what's built, what's
validated, the honest verdict, and the remaining path to supplant.

- **Branch:** `raw-gated-scoring` (off `v6-graduated-tail`). Uncommitted as of this writing.
- **Analysis provenance / science:** `/mnt/data/u12/ipa/conservation_corpus/eval_corpus/` —
  `FINDINGS.md` (§6–§6h, the empirical case), `PROPOSED_ARCHITECTURE.md` (design + panel), `WHY_RAW_NOT_Z.md`
  (the conceptual "why", with the global-vs-per-species control). Cross-session record: the
  `eval-corpus-framework` memory.
- **Status (2026-06-26):** implemented end-to-end, runs in the real pipeline, and on a 6-species head-to-head
  **matches the heavily-tuned production v6 on loss-FP control with bearer recall preserved** — at a fraction
  of the complexity. Not yet supplant-ready (small panel, post-hoc τ, build-script bundle, output not fully
  synced).

## 0. MAJOR UPDATE (2026-06-27): `P_motif`-as-default supersedes the per-species layer

A broad-panel investigation (chasing an interpretable per-intron probability) found that **the per-species
adjustment — gate, discount, AND the empirical-Bayes prior shift — is largely unnecessary and, for
bearers, actively harmful.** The raw-motif SVM score *alone*, calibrated to a probability `P_motif`, is the
robust default for essentially every real genome.

- **`P_motif` = σ(a·margin + b)**, the Platt-calibrated **ensemble decision_function margin** of the
  raw-feature SVM (use the *margin*, not the isotonic prob, which saturates at 0/1). Species-agnostic.
- **Why the per-species shift is wrong (the z-norm mistake, redux):** every real genome has U12 prevalence
  **< 1%** (only *Physarum* exceeds it; human ~0.33% is merely top-decile), so a posterior referenced to the
  balanced-training prior (π=0.243) deflates *every* genome — there is no neutral genome. Worse, the
  training U2 set carries hard high-margin negatives, so the training likelihood **understates** how clean a
  real bearer's U2 background is; empirically the motif is *more* specific in the wild than in training.
  `P_motif` alone recalls **623/630 human conserved U12s with ~741 total calls** (≈ the true U12 count); the
  prior shift dropped that to 565. The prior shift's "interpretable probability" goal is real, but it is
  delivered by `P_motif` itself + the **two-number split** (`P_motif` species-agnostic; an optional
  per-species posterior only where a genome is *confirmed* depleted).
- **Broad-panel validation (the decision-grade evidence):**
  - **FP-inflation, 19 snRNA-confirmed true losses** (every U12-absent clade: green algae, ciliates,
    Apicomplexa, Dikarya fungi, nematode, kinetoplastid, microsporidian, diplomonad):
    `P_motif≥0.9` gives **median 3, max 15, mean 4.6 FPs/genome — no inflation anywhere** (the irreducible
    motif-strong floor). Script `eval_corpus/loss_fp_one.py`, table `lossfp2.tbl`.
  - **Bearer precision** (`eval_corpus/bearer_precision.py`): the called set is **85–100% canonical-U12-5′SS
    -concordant down to the bottom**; non-conserved calls are real U12s the conservation set missed, not
    junk (human 89%, drosophila/aphanomyces 100%, oryza/amborella 92–94%, gut fungi 85–96%).
  - **ascaris** (`role=loss`, was `SNRNA_NOT_SEARCHED`): fresh cmsearch → **all 4 minor snRNAs absent**
    (definitive loss); `P_motif` = **7 calls** (loss-like). Manifest updated to LOSS/AGREE_loss.
  - **symbiodinium** (the conflict case: ambiguous snRNA + 840k introns, historically FP-prone): intronIC
    re-run in flight to get its `P_motif` count vs v6 (sister strain Smic_CassKB8 has u11+u4atac @4.4e-5 →
    divergent bearer). [pending]
- **Implied architecture:** **`P_motif` is the call everywhere; a *gated loss-detector* flips on suppression
  only for confirmed-depleted genomes** (signal: high-`P_motif` rate ≲2e-4 AND one-end-strong shape vs
  bearers' all-ends-strong/U12-5′SS, + snRNA-absence backstop). This retires the uniform prior shift and
  collapses the z/mode-sep/discount/gate stack to one motif score + a thin, evidence-gated loss layer.
- **Open:** build + validate the gated loss-detector (fire on all 19 losses, none of the bearers — the
  loss-vs-very-poor-bearer boundary at count ~15 vs aph 22 needs the shape signal); fill 2 old-schema gaps
  (symbiodinium [running], ustilago); the margin/Platt base calibration should be made leave-clade-out (v1
  used in-sample templates). Artifacts: `eval_corpus/{prior_shift_v2_margin.py, threshold_tradeoff.py,
  loss_fp_one.py, bearer_precision.py, pmotif_validation.png, panel_recall_cost.png}`.

### 0a. The species-level adjudicator — robust depth-beyond-U2-tail (2026-06-27)

Converged architecture (with the user): **per-intron `P_motif` (granular call) + a per-species
separation statistic that adjudicates the *motif-strong FPs that don't form a coherent cluster*** (loss /
recent-loss). The species signal is a confidence layer that *complements*, never overrides, per-intron
`P_motif`.

- **The existing `cluster_validation.compute_valley_depth` (median_depth / gap_fraction) is size-blind**:
  it requires the U12 mode to have density *mass*, so it FALSE-NEGATIVES on U12-poor bearers —
  **aphanomyces (8 IPA U12s) gives median_depth 0.00**, identical to losses, even though the scatter
  plainly shows a separated U12 cloud at 5′≈13 (`eval_corpus/separation_scatters.png`). 22 U12 points
  among 56k U2 don't dent the density. A formal mass-based clustering (GMM/silhouette) would re-break the
  same way — so it offers no advantage here.
- **The fix (size-robust, principled, NOT special-cased):** `depth_tail` = (robust median call-core −
  U2's own 99.9th pct) / MAD_U2 — center/tail-based not density-based, and the *median* call-core peels the
  U2-edge fringe of mis-calls (the user's point). Bearer call-clouds sit **deep** beyond their U2 tail,
  ancient-loss FP-clouds sit **shallow/tail-adjacent**. Boundary set (`eval_corpus/robust_sep_one.py`):
  all 4 losses `depth_tail` **≤2.10**, all 5 bearers **≥3.39** (incl. divergent neocal/amborella and
  U12-poor aphanomyces 4.69 / drosophila 6.09) — clean gap ~2.7, **and aphanomyces flips to clearly
  separated** where the density-valley failed.
- **Symbio resolves to LOSS at the single-genome level — no snRNA needed:** `depth_tail` **1.30** (lowest of
  all, below the ancient losses). Its 151 motif-strong `P_motif` calls hug its *own* (long, U12-like) U2
  tail rather than forming a deep coherent cluster — exactly the **recent-loss "long U12-like U2 tail"**
  signature (user's hypothesis). v6 also gave symbio 0 calls; symbio has no IPA-conserved U12 and (this
  assembly) no snRNA, so loss/recent-loss is the right read; the only bearer hint is a sister strain's
  snRNA. (Note: a recent-loss genome divergent enough is single-genome-irreducible vs a divergent bearer —
  that residual is functional, not motif.)
- **Principled rule (replacing the empirical ~2.7 cut):** an **excess-over-U2-tail outlier test** — fit the
  U2 margin tail, extrapolate to the call region, ask whether there are significantly *more* high-margin
  introns than U2's own tail predicts (`log_excess` / Poisson `z_excess`). Loss: calls ARE the U2 tail
  (excess≈0). Bearer: a separate population U2 can't explain (large excess). Threshold-free (significance).
- **Validation IN PROGRESS:** `depth_tail` + the principled excess across the full 37-genome panel (all 19
  losses + ascaris + symbio + 14 bearers) to confirm the gap holds and pick the rule. Artifacts:
  `eval_corpus/{robust_sep_one.py, valley_detect_one.py, separation_scatters.py, separation_scatters.png}`.
- **Principled outlier rule (the threshold replacement):** the separation cut is an **EVT test on the MAD
  distance** — is the `P_motif` call-core beyond the genome's OWN U2 *expected maximum* (Gumbel, adaptive to
  N)? Loss: call-core ≈ U2 expected-max (it IS the tail) → `p_gumbel ≈ 0.5–1.0`. Bearer: beyond it → a
  separate population → `p_gumbel ≪ α`. Sanity (5 cases): bearers p≈3e-4–9e-3, losses/symbio p≈0.5–1.0.
  Decision = `p_gumbel < α` (~0.01), not an empirical depth cut. (`robust_sep_one.py`: `excess_z`/`p_gumbel`.)
- **Training-confound check (resolved-in-principle):** the panel is *not* meaningfully confounded by
  training membership — raw-motif features are fixed-PWM + genome-intrinsic U2 background (not learned), the
  SVM sees only the 6 universal features (never genome identity → no per-genome memorization), and the
  separation statistic references each genome's *own* U2 tail. The only residual is ordinary classifier
  generalization, already bounded by leave-clade-out AUC 0.916. **PLANNED before finalization (deferred,
  user 2026-06-27):** re-score the separation panel with **leave-clade-out (Stage-1 OOF) margins** so every
  genome's `P_motif` is strictly out-of-sample — bulletproof confirmation that `depth_tail`/`p_gumbel`
  separation isn't inflated by in-sample optimism. Low-risk but required for the final write-up.
- **FULL PANEL RESULT (37 genomes, 2026-06-27):** BOTH statistics separate ALL 13 bearers from ALL 17
  computable losses with a clean gap — `excess_z` bearers +0.55..+3.88 vs losses −2.25..+0.30; `depth_tail`
  bearers ≥3.39 vs losses ≤2.91. Hard cases all correct: aphanomyces (density-valley FN) → bearer +2.35;
  symbio → lowest (−2.25) → loss. 8 intron-poor losses → "too few" (no U12 population) → trivially loss.
  Figure `eval_corpus/panel_separation.png`, data `/tmp/fullpanel_results.txt`.
- **AIRTIGHT PASS (deferred; both items together):** (1) leave-clade-out OOF margins [above]; (2) replace the
  exponential (ξ=0, MoM, fixed-90th-pct) U2-tail with a **proper POT-GPD** (shape ξ via MLE/PWM,
  stability-selected threshold) → honest `p_gumbel`, decision = `p < α`. The 3 borderline losses (chlamy/
  aspergillus/schizo at excess_z +0.19/+0.20/+0.30) are the exponential-tail underestimate; the GPD either
  re-centers them toward 0 or exposes a real motif-strong-FP floor. **User (2026-06-27): worth doing but
  won't change the verdict — the margin is large.**
- **AIRTIGHT PASS PART A — DONE (2026-06-27, `eval_corpus/airtight_oof.py`): in-sample optimism does NOT
  inflate the depth_tail separation.** Re-scored all 24 computable panel genomes with an ensemble that
  EXCLUDES each genome's own training clade (strictly OOF; diaspora genomes not in the 41,257-row v23 corpus
  — oomycetes, ciliates, green algae — are already OOF via the full ensemble), with the Platt re-fit on
  pooled leave-clade-out OOF margins. The like-with-like signal is **`delta = dt_oof − dt_full` (same serial
  ensemble, only the clade-holdout differs): small for ALL genomes — vertebrate bearers +2.4..+3.3 (more
  separated OOF), plant bearers −0.7..−1.25, gut fungi ~0; losses ascaris/caeno ~0, aspergillus −1.65,
  toxo/vitrella −0.3** — i.e. holding out a genome's clade barely moves its depth_tail, so the bearer/loss
  separation is NOT a training-overlap artifact. Caveat: the 15-model *serial proxy* ensemble (the real
  126-SVC bundle's joblib pool segfaults in this exec context) has a tighter U2 scale (depth_tail ~3–5× the
  42-model bundle) and noisier loss-FP suppression, so its leave-clade-out Firth-q gives **22/24 correct**
  with 2 *proxy-scale* boundary cases (aphanomyces q=0.39, caenorhabditis q=0.58) whose **`delta`≈0 (NOT an
  OOF effect)** — both separate cleanly in the canonical 42-model bundle (caeno dt=2.17 loss, aphan dt=4.69
  bearer). **Verdict: Part A confirms the separation is real, not in-sample-inflated.** Corpus: 41,261
  introns (10,044 U12 / 31,213 U2, 13 clades).
- **AIRTIGHT PASS PART B — DONE (2026-06-27, `eval_corpus/gpd_tail_analysis.py`): keep the exponential
  secondary; the borderline residual is a real motif-strong-FP floor.** Fit a proper POT-GPD to the U2-margin
  tails. The motivating premise (heavier-than-exponential tail) is **wrong in sign**: shape ξ is robustly
  **NEGATIVE** (~−0.10..−0.16; bootstrap CIs exclude 0 at high exceedance counts — aspergillus ξ=−0.160,
  schizo ξ=−0.126) — the U2 tail is *lighter*-than-exponential / bounded, so the exponential is already mildly
  conservative. The GPD does NOT re-center the 3 borderline losses toward 0: it shifts the *whole* panel up by
  a near-constant ~+1.5 via genpareto MLE's small-sample negative-ξ bias (reproduced ~+1.0 on **genuinely
  exponential** data), an estimator artifact, not signal. The +0.2/+0.3 borderline residual (chlamy/asperg/
  schizo) survives both models and is **independently corroborated by the non-parametric `depth_tail`**
  (chlamydomonas ranks *above* the clear losses on depth) → a genuine **motif-strong-FP floor**, not a
  tail-model deficiency; no tail model should absorb it. Bearer/loss ordering unchanged. → no code change to
  the (labelled, diagnostic-only) secondary.
- **Low-N fallback — SETTLED (2026-06-27).** When a species has fewer than `min_u2` (200) U2 introns it
  cannot reference its own tail. **Default = option A:** `LOW_N → q_eff=1` (P_adj=P_motif, no suppression —
  matches "P_motif everywhere; suppress only confirmed-depleted"; correct for the dominant `-q`
  curated-candidate case, and production genome runs never hit it). **Option B (opt-in):**
  `--adjudicator-min-u2 <N>` (config `scoring.adjudicator_min_u2`) lowers the threshold so a small genome
  **self-adjudicates against its OWN (noisier) U2 tail** — recovering tiny-loss-genome FP suppression.
  Verified end-to-end (147-intron `-q`: default → LOW_N q=1.0; `--adjudicator-min-u2 50` → self-adjudicates,
  depth_tail=3.31 q=0.76). A **pooled GLOBAL U2 reference was built and REJECTED**: a single global tail
  can't span the per-genome q99.9 spread (−5.35..−1.69), so it is dominated by long-U2-tail genomes and
  **over-suppresses clean-background ones** (suppressed 37 real human U12s → q=0.07, even with correct
  features). Option B avoids this because each genome references its OWN U2.

### 0b. CLOSED scoring design — two per-intron numbers + a hierarchical uncertainty (2026-06-27)

The scoring is now design-complete (validation/refinements deferred per §0a). **Two per-intron outputs:**

- **`P_motif`** = σ(Platt(ensemble margin)) — pure sequence-level motif congruence, **species-agnostic**.
  A textbook-`GTATCC` U12 scores ~0.99 regardless of how many siblings its genome has. The granular call.
- **`P_adj` = `q` · `P_motif`** — per-intron posterior that the intron is a *functional* U12, where
  `q = P(this genome is a real U12 bearer)` from the species adjudicator. Clean because
  `P(U12 | motif, not-a-bearer) ≈ 0` (no spliceosome → no functional U12s). The **low-numbers uncertainty
  flows through `q`'s bootstrap CI** into `P_adj`.

**`q` calibration + bootstrap** (`eval_corpus/q_bootstrap.py`): `q = σ(2.20·excess_z − 1.31)` (logistic on
the panel labels, regularised → smooth; q=0.5 at excess_z=+0.59, in the loss/bearer gap). Per-genome CI by
**bootstrapping the call-core** (resample the `P_motif≥0.9` calls → distribution of excess_z → of `q`).
Results: bearers `q` 0.77–1.00 (drosophila 18 calls→1.00, aphanomyces 22→0.98, saprolegnia 34→0.77,
human→0.87); losses `q` 0.04–0.29 (tetrahymena→0.04, chlamydomonas 3 calls→0.29) — so a loss's motif-strong
FP gets `P_adj≈0.04` despite `P_motif≈1.0`, while bearer U12s keep high `P_adj`.

**Key (refines "low-N ⇒ uncertain"):** uncertainty is driven by **depth (distance past the genome's U2
expected-max), not raw N**. Down-sampling human's deep U12 core to **k=2** still gives `q≈0.99` confident —
two *genuinely deep* U12s are individually too far beyond what `N_u2` draws can produce by chance to be FPs
(the EVT test already prices this in). Low-numbers uncertainty surfaces **near the boundary** (moderate-depth
/ divergent low-N bearers → wide CI → *undetermined*), never as "loss". So an arbitrarily-low-N bearer is
captured as **confident** (deep motifs) or **undetermined-wide-CI** (borderline) — **count is never a
threshold**; it only modulates `q`'s CI. Matches the user's intent exactly.

### 0c. IMPLEMENTATION PLAN (wire the closed design into the real pipeline)

Replaces the §5 supplant steps (those targeted the older species_gate). Ordered:

1. **Margin-capable bundle.** The raw ensemble must expose the **decision_function margin** (currently the
   bundle's `P` is isotonic-saturated). Store, in the bundle: the raw SVC ensemble (margin), the **Platt**
   (margin→`P_motif`) params `(a,c)`, the **`q`-calibration** `(qa,qb)`, and the **EVT/adjudicator settings**
   (U2 def = `P_motif<0.5`, tail-fit recipe, bootstrap B). Stamp `scoring_mode:"pmotif_adjudicated"`.
2. **Adjudicator module** (`scoring/species_adjudicator.py`, extends/replaces `species_gate.py`): pure fns
   — `evt_params(margins,P_motif)` → (med,MAD,tail) ; `excess_z(call_core, params)` ; `q_and_ci(calls)`
   (calibrated `q` + bootstrap CI). Plus a file-side `apply_pmotif_adjudication(score_info)` post-process.
3. **Calibration in `main_train`** (or the build script first): fit Platt (margin→`P_motif`, **OOF**), fit
   `q` (excess_z→`P(bearer)`) on the labelled adjudicator panel, freeze EVT settings → bundle. The `q`-fit is
   a *new* labelled calibration input (the snRNA/IPA panel) — document its provenance.
4. **Inference wiring** (`cli/main.py`, the existing `scoring_mode` dispatch): per-run, compute `P_motif`
   from margins, run the adjudicator → `q`+CI, write `P_adj=q·P_motif` + CI. Must run identically on the
   **streaming and in-memory** paths (the v2.7.1 parity hazard — add the parity test).
5. **Output schema** (`file_io/writers.py`): `score_info` + `meta` gain `P_motif, P_adj, q, P_adj_lo,
   P_adj_hi`; `type_id` from a `P_adj` threshold (≥0.5 called / ≥0.9 high-conf — decide), with the CI
   surfaced; `metrics.iic.json` gains the species `q`+CI. Update `docs/scoring_pipeline.md` schema.
6. **Gates:** chr19 smoke + streaming==in-memory parity; then the §0a airtight pass (POT-GPD tail +
   leave-clade-out OOF margins) before finalising.
7. **Supplant PR:** flip the bundled default to the `pmotif_adjudicated` model, DELETE the z stack
   (`normalizer.py` adaptive-z, `mode_separation.py`+`mode_sep_pipeline.py`, `prior_adjustment.py`).

**Implementation status (2026-06-27, post-committee — see §0d):** `scoring/species_adjudicator.py` is built and
unit-tested (`tests/unit/test_scoring/test_species_adjudicator.py`, 22 tests). Done:
- **Ship-blocker #1 (depth_tail swap):** `depth_tail = (call-core − q99.9_U2)/MAD_U2` is the PRIMARY,
  size-invariant q-driver. The old size-aware `excess_z`/`p_gumbel` are retained as a labelled **secondary**
  diagnostic (`secondary_available` flips false when the exponential tail can't be fit — `depth_tail` does not
  depend on it). The new fit supersedes the excess_z fit (2.20/−1.31).
- **Ship-blocker #2 (DONE):** the q-fit is now a **separation-safe Firth-penalized** logistic (the panel is
  cleanly separated so a plain MLE diverges; Firth = Jeffreys-prior penalty): **`q = σ(3.64·depth_tail −
  10.86)`** (q=0.5 at depth_tail=2.98, inside the gap [2.91, 3.39]). **Leave-clade-out validated**
  (`eval_corpus/q_firth_leaveclade.py`, 13 clades): **26/27 genomes classify correctly out-of-sample**; the
  lone OOF miss is **chlamydomonas** (deepest loss, 3 calls) tipping to q≈0.60 *only* when its whole
  Chlorophyta clade is held out — and with 3 calls its bootstrap CI is wide → it surfaces as `UNDETERMINED`,
  not a confident bearer FP (exactly the undetermined band's job). Symbiodinium recoded **CONFLICT** and
  excluded from the loss class. The **undetermined band** = "bootstrap q-CI straddles 0.5 → `UNDETERMINED`".
- **Ship-blocker #3 (operational guards):** `AdjStatus` codes
  (`ADJUDICATED`/`UNDETERMINED`/`LOW_N`/`DEGENERATE_TAIL`/`SCHEMA_FAIL`), the degenerate-branch contract (MAD≈0
  and non-finite paths return a status, never a silent NaN), shape/NaN input guards, and a version pin
  (`ADJUDICATOR_PARAMS_VERSION = "depth_tail_firth_2026-06-27"`, `AdjudicatorParams.params_version`).

All three ship-blockers are now closed in `scoring/species_adjudicator.py`.

**Pipeline integration — steps 1, 3, 4 DONE (2026-06-27):**
- **Step 1 (margin-capable bundle):** `scripts/build_pmotif_adjudicated_bundle.py` re-stamps a raw-feature
  ensemble into a `scoring_mode:"pmotif_adjudicated"` bundle carrying `adjudicator_params`
  (Platt + q + EVT settings, version-pinned) + the raw SVC ensemble (margin via `decision_function`).
- **Step 3 (inference wiring):** `apply_pmotif_adjudication(score_info, ensemble_models, params)` (file-side,
  mirrors `apply_raw_gated_postprocess` but uses the **margin**, not the saturated `svm_score`) computes
  `P_motif` → adjudicate → `P_adj = q_eff·P_motif`. Dispatched from the **shared**
  `_run_post_classification_pipeline` (`cli/main.py`), so streaming and in-memory are identical *by
  construction* — verified bit-identical on chr19 (`type_id`/`P_motif`/`P_adj`/`adjusted_score`/`rel_score`).
  Inconclusive species (`LOW_N`/`DEGENERATE_TAIL`/`SCHEMA_FAIL`) fall back to `q_eff=1` (`P_adj=P_motif`, no
  silent suppression).
- **Step 4 (output schema):** `score_info.iic` gains `P_motif, q, P_adj, P_adj_lo, P_adj_hi`; `type_id` from
  `P_adj≥0.5`; `adjusted_score=100·P_adj`, `rel_score=100·P_adj−90` (existing conventions preserved).
- **meta.iic/bed.iic sync (DONE):** the file-side post-process modes write their calls *after* meta/bed are
  already on disk, so the shared dispatch now calls `_sync_calls_to_meta_and_bed` (reuses
  `_sync_meta_from_score_info` for meta `type_id`+`rel_score`; rewrites bed col 4 ← `adjusted_score`). This
  fixes a real correctness bug — stale meta `type_id` silently corrupted the `metrics.iic.json` boundary
  tables whenever a call flipped. Applied to **both** `pmotif_adjudicated` and `raw_gated`. Verified on
  chr19: meta/bed full-row **identical** streaming vs in-memory; meta `type_id`/`rel_score` and bed col 4
  match `score_info` 1.0000.
- **End-to-end verified on all three input modes:** genome+annotation (`--in-memory` and `--streaming`,
  bit-identical, chr19 → q=0.993 bearer, 39 U12 calls); `-q` pre-extracted sequences full (bit-identical to
  the genome path); `-q` low-N subset (120 introns → `LOW_N` → `q_eff=1` fallback, `P_adj=P_motif`, no crash).
  Unit-tested in `test_species_adjudicator.py` (file-side glue + low-N fallback) — 24 tests green.

**Step 2 (`main_train` emission) — PARTIAL + a key calibration finding (2026-06-27):** a `--scoring-mode
{zscore,raw_gated,pmotif_adjudicated}` train flag (`cli/args.py`, `cli/config.py`) now routes `main_train`
through `_train_and_save_raw_bundle` (`cli/main.py`): for the raw modes it **skips z-normalization +
`IntronClassifier`** (both z-only) and trains the raw-feature ensemble directly (mirrors
`build_raw_gated_bundle.py`), stamping `scoring_mode` + the mode's post-process params. Verified: emits a
valid, loadable, runnable bundle (e.g. `intronIC train --scoring-mode pmotif_adjudicated --n-models 42`).

**BUT — the calibration (Platt + q) is ensemble-specific and does NOT transfer to a re-trained ensemble.**
Empirically, a `main_train` 42-model ensemble on the *same* built-in v23 reference (10003 U12 / 31330 U2)
gives chr19 `depth_tail≈2.7` / `q≈0.30` (UNDETERMINED) vs the canonical eval ensemble's `depth_tail≈4.4` /
`q≈0.993` — the frozen `DEFAULT_PLATT_*` / `DEFAULT_Q_*` constants are tied to one specific trained
ensemble's margin/`depth_tail` scale and don't survive re-training (different feature provenance + SVM
non-determinism shift the scale). So `main_train` emits the **ensemble + a hard `messenger.warning` + a
`calibration_provenance` stamp**, but a freshly-trained pmotif bundle is **NOT production-calibrated**.
Consequence for the supplant plan: a shippable `pmotif_adjudicated` default is fundamentally an **eval_corpus
artifact** — train ensemble → **fit Platt OOF on the reference** → **fit q on the species panel for THAT
ensemble** → assemble. `main_train` provides the ensemble-training + z-bypass plumbing; the q-fit
specifically needs the dev snRNA/IPA panel and cannot be done from `intronIC train` alone.

**Airtight pass §0a — both parts DONE** (Part A leave-clade-out OOF margins; Part B POT-GPD diagnostic; see
§0a). **Low-N fallback — SETTLED** (option A default + `--adjudicator-min-u2` opt-in; pooled global-q
rejected; see §0a "Low-N fallback — SETTLED"). Committed pmotif parity test DONE.

**Calibration pipeline — DONE (2026-06-27, `eval_corpus/calibrate_pmotif_bundle.py`): the shipped constants
are confirmed self-consistent for the canonical raw_gated_42 ensemble, and the derivation is reproducible for
any ensemble.** The pipeline (1) re-derives **Platt** from the ensemble's margins over the v23 corpus →
`σ(2.7958·margin − 1.1778)`, **reproducing the shipped `DEFAULT_PLATT_*` (2.796/−1.178) exactly** (|Δ|=0.000,
in-sample ECE 4e-4); (2) re-computes per-genome `depth_tail` (this ensemble) and Firth-fits **q**,
leave-clade-out **24/25** (the lone miss is chlamydomonas — the known §0a borderline loss); (3) assembles the
bundle with the re-derived Platt + the **shipped full-resolution q** (`DEFAULT_Q_*` 3.64/−10.86). Important:
the CAP-subsampled q re-fit (4.31/−11.9, boundary 2.76) is a **method cross-check only** — *inference computes
`depth_tail` at full resolution*, so the bundled q must be the full-res fit; bundling the subsampled re-fit
would be a scale mismatch. The **shippable bundle** = re-stamp of raw_gated_42 + these verified constants
(`scripts/build_pmotif_adjudicated_bundle.py`; version `depth_tail_firth_2026-06-27`); end-to-end chr19
q=0.993, 39 calls.

**SUPPLANT STEP 1 — DONE (2026-06-28, the reversible default flip; z stack KEPT).** The verified pmotif
bundle (re-stamp of `raw_gated_42` + the version-pinned calibrated constants, 2.0 MB) is now placed at the
bundled-default path `src/intronIC/data/default_pretrained.model.pkl` (replacing the 27 MB v3 modesep bundle,
preserved in git history + backed up out-of-repo). It stamps `scoring_mode:"pmotif_adjudicated"`, so the
shared dispatch's pmotif branch fires by default and returns before the modesep check — **default runs now
score via pmotif** (chr19, no `--model`: status ADJUDICATED, q=0.993, 39 U12 calls). The z stack
(`normalizer.py` adaptive-z, `mode_separation.py`+`mode_sep_pipeline.py`, `prior_adjustment.py`) is LEFT
INTACT — reversible (revert the bundle to restore modesep), and modesep stays fully testable against its own
explicit `v5_modesep` bundle (those tests pin it, not the default, so they're unaffected). Validated under the
new default: input-modes + CLI (36), unit suite (723), modesep tests skip cleanly (dev-only bundle absent),
streaming==in-memory parity [the canonical integration gate]. **Step 2 (delete the z stack) remains** — a
separate PR once this default has soaked. The §5 formal gates / §0a airtight pass are complete (§0a).

**SUPPLANT STEP 2a — DONE (2026-06-28, the post-SVM gate deletion; normalizer rip-out deferred to 2b).**
Recon revealed the z stack is two distinct-risk pieces: (i) the post-SVM mode-sep + continuous-discount +
cluster-validation GATE, cleanly removable without touching the scoring path; (ii) the adaptive-z
NORMALIZER, woven through the scoring/streaming hot path and currently running harmlessly under raw
bundles. 2a does (i). Tagged `pre-zstack-removal` (old zscore/modesep + v3/v6 bundles reproducible there),
then: extracted `file_io/meta_sync.sync_meta_from_score_info` out of the doomed `mode_sep_pipeline` (the
pmotif path's only dependency on it); removed the legacy dispatch branch from
`_run_post_classification_pipeline` (raw_gated + pmotif return early); deleted the orphaned helpers
`_apply_post_classification_adjustment` / `_write_pi_species_adjusted_scores` and the dead
`_apply_margin_alignment` wrapper; deleted the four gate modules `scoring/mode_separation.py`,
`scoring/cluster_validation.py`, `scoring/margin_alignment.py`, `classification/mode_sep_pipeline.py`.
Added a fail-fast guard `model_io.assert_scoreable_bundle` (called at both classify load sites): a legacy
zscore/modesep bundle now errors AT LOAD with a message pointing to the `pre-zstack-removal` tag (verified
on the real 27 MB v3 bundle). Flipped the default `--scoring-mode` to `pmotif_adjudicated` (zscore kept as
a deprecated/unscoreable choice). `main.py` 7100→6690 lines; 9 z-stack test files removed +
`tests/unit/test_bundle_guard.py` added. VERIFIED: chr19 byte-identical to a pre-surgery pmotif golden
(score_info/meta/bed/introns, BOTH `--streaming` and `--in-memory`); full unit suite 669 passed;
streaming==in-memory parity gate green. KEPT (NOT part of the gate): `scoring/prior_adjustment.py` still
backs the live `--species-prior` Bayesian adjustment; `scoring/normalizer.py` (adaptive-z) — runs
harmlessly under raw bundles and the in-memory cv-array filter still keys on `*_z` presence.

**SUPPLANT STEP 2b — the deeper normalizer/z-feature rip-out (multi-pass; IN PROGRESS).** 2b removes
z-normalization from the scoring/streaming hot path entirely and deletes `normalizer.py` + its closure.
Recon (the `ScoreNormalizer` closure spans 9 src files) split it into four verifiable sub-passes:

- **2b-1 — DONE (commit 3cedf7a).** Remove the dead zscore *training* + bundle loader (inference output
  unchanged, chr19 byte-identical). `utils/model_io.py`: deleted the v3/modesep translator
  (`_build_v3_ensemble` / `_v3_to_runtime` / `_load_fallback_normalizer`); `normalize_model_bundle` is now a
  pass-through (a v3/zscore bundle passes through and is rejected at `assert_scoreable_bundle`). `main_train`'s
  zscore `else`-branch and `main_classify`'s no-pretrained "normalize + train + classify" `else`-branch now
  `raise` (both used `ScoreNormalizer` + `IntronClassifier`). `--scoring-mode` drops the `zscore` choice.
  `test_model_io` rewritten for the pass-through loader. `classify_introns` / `IntronClassifier` /
  `normalize_scores` are now DEAD (deleted in 2b-4 with `normalizer.py`, which they import).

- **2b-2 — DONE (commit 7befa5d), the risky core, PROVEN.** Ripped z-normalization out of both classify hot
  paths. **Controlled probe:** with the adaptive `RobustScaler` fully skipped, EVERY scoring column
  (`svm_score`/`P_motif`/`P_adj`/`q`/`type_id`/`adjusted_score`/`rel_score`/raw/CIs) is byte-identical to the
  pre-change golden ⇒ z-norm is irrelevant to the calls (it only produced unused diagnostic columns).
  `predictor.py`: guard the three BothEndsStrong `min_5_bp`/`min_5_3` blocks (computed unconditionally from
  z-scores → `None+None` crash on raw bundles → `None` when z absent). In-memory predicts on the raw-scored
  introns (no `ScoreNormalizer.fit/transform`); the streaming worker drops `apply_scaler_to_scored_batch` and
  classifies the raw batch directly; both the cv-filter and the streaming accumulator key on `svm_score` and
  carry raw motif scores. chr19 both modes: scoring columns identical to golden **and** streaming==in-memory
  parity perfect (0 mismatches, all columns). The `5'_z`/`bp_z`/`3'_z`/`min_5_bp`/`min_5_3`/
  `svm_score_adaptive`/`svm_score_frozen` columns are now `NA` in both paths.

- **2b-3 — REMAINING.** Drop the now-`NA` z/min/adaptive columns from `ScoreWriter` (header + data) and
  CONVERT the 2D/3D scatter plots in `visualization/plots.py` (4 functions read `5'_z`/`bp_z`/`3'_z`) to
  **raw-feature space** `[5'_raw, bp_raw, 3'_raw]` (the actual classifier space; user decision). Re-anchor the
  chr19 golden (z columns gone), re-verify parity.

- **2b-4 — REMAINING.** Delete `scoring/normalizer.py` (`ScoreNormalizer` / `ZeroAnchoredRobustScaler` /
  `DatasetType`), the now-dead streaming adaptive-scaler pre-pass + in-memory scaler resolution, the dead
  `classify_introns` / `IntronClassifier` / `normalize_scores`, `trainer.Z_BASE_FEATURES` + the optimizer
  z-default, the `scoring/__init__` `ScoreNormalizer`/`DatasetType` exports, `clipping.py` docstrings, and the
  remaining z-stack tests (`test_normalizer`, `test_new_scaling_architecture`, the z bits of `test_scorer` /
  `test_classifier` / `test_edge_cases`). Also clear the stale always-None `cluster_validation_result` /
  `mode_separation` metrics plumbing in `main.py`, and update `CLAUDE.md` + `DEV_PATHS_MANIFEST.md` layout.

**Committed parity test (DONE):** `tests/integration/test_pmotif_adjudicated_parity.py` builds a tiny
1-model `pmotif_adjudicated` bundle on the fly (`intronIC train --scoring-mode pmotif_adjudicated`; no
committed binary, no dev-panel dependency — parity is independent of calibration) and asserts streaming vs
in-memory `score_info` + `meta` are byte-identical and the `P_motif/q/P_adj/P_adj_lo/P_adj_hi` columns are
present. ~100s, skips gracefully without the dev env / chr19 data.

**Clarified supplant path (2026-06-27):** a shippable default does NOT require the from-scratch calibration
pipeline — the canonical eval ensemble (`raw_gated_42.model.pkl`) + the version-pinned DEFAULT constants is
*already* a calibrated, validated pmotif bundle (re-stamp via `scripts/build_pmotif_adjudicated_bundle.py` →
chr19 q=0.993, clean 37-genome panel). The step-7 supplant PR is therefore: place the re-stamped canonical
bundle at the bundled-default path, flip the default scoring_mode, delete the z stack — gated only on the
§0a airtight pass (leave-clade-out OOF margins), which the user has said "won't change the verdict."

### 0d. Committee review + meta-review record (2026-06-27) — auditability

The §0a–§0c statistical core (technical write-up: `eval_corpus/STATISTICAL_METHODS.md`) was put through an
**8-agent review**: a 5-agent expert committee, then a 3-agent meta-review (statistical adjudicator,
deployment steelman, citation fact-checker) to separate legitimate findings from incongruous ones and to
confirm method claims by direct literature research (user's instruction). Full transcripts saved at
`eval_corpus/reviews/{01_EVT_methodology, 02_calibration_bayes, 03_standard_methods_comparison,
04_genomics_domain, 05_adversarial_redteam}.md`.

**Consolidated verdict: the method is sound and materially better than the z-norm stack it replaces.** No
finding overturns the architecture. Three *cheap* ship-blockers; everything heavier was judged non-blocking
or rejected.

**Three ship-blockers (committee-endorsed, all low-cost):**
1. **Swap the `q`-driver `excess_z` → `depth_tail`.** `excess_z`'s expected-max carries a `ln(N_u2)` term
   (the Gumbel location grows with U2 count), making it a *size-aware significance* test (`Q_sig`: "beyond
   THIS genome's U2 chance, look-elsewhere-corrected"). The authors actually want `Q_pop` ("a real U12
   *population* regardless of genome size"): **`depth_tail` = (median_call_core − q99.9_U2)/MAD_u2**, which is
   **size-invariant** (no `ln N`). It already exists in `eval_corpus/robust_sep_one.py` (three lines above
   `excess_z`) and separates the panel *better* (bearer/loss gap +0.48 vs excess_z +0.25). Keep `p_gumbel`/
   `excess_z` as a **labelled secondary** (the significance view), not the q-driver. **This is the fix for the
   "#2 N-dependence" objection below.**
2. **Separation-safe `q` calibration.** The 37-genome panel is cleanly separated → plain logistic
   coefficients diverge. Use **Firth / log-F(1,1) penalized** (or weakly-informative Bayesian) logistic;
   add **leave-clade-out (clade-weighted)** validation and an explicit **undetermined band** (wide-CI →
   "undetermined", never silently "loss"). Re-code *Symbiodinium*'s label: it is a **CONFLICT** case (this
   assembly snRNA-absent, sister strain positive), currently entered as a hard loss in the q-fit — exclude or
   down-weight it rather than letting it anchor the loss class.
3. **Operational guards.** Status codes (`ADJUDICATED` / `LOW_N` / `SCHEMA_FAIL` / `DEGENERATE_TAIL` / …), a
   **degenerate-branch contract** (the `len(exc) ≤ 20` NaN path in `_u2_expected_max` / its `depth_tail`
   analog must return a status, not a silent NaN), NaN guards throughout, and a **version pin** on the
   calibration constants in the bundle.

**Rejected / deferred as over-engineering (do NOT block on these):** a full POT-GPD shape-ξ pipeline
(diagnostic-only — margin is large; §0a airtight pass); hierarchical Bayesian `P(bearer)`; a 3-component
Olthof "hybrid-class" mixture (see citation check); a Markov/HMM tail rebuild. The exponential (ξ=0) tail is
an acceptable approximation given the gap; GPD stays a one-time diagnostic.

**Citation fact-check (direct-research mandate). Toolbox confirmed solid; two corrections, one vindication:**
- **CORRECTION — Cheng-Hall mis-attribution.** The `n^−3/5` modality/mode-detectability rate was cited to
  "Cheng-Hall"; direct check re-attributes the relevant detectability result to **Hall / Hall & York** (the
  calibration-of-the-dip-test and bandwidth-detectability line). Re-cite accordingly; the *claim* (KDE/dip
  modality tests are size-blind at tiny π₁) stands and is in fact the reason valley-depth false-negatives
  aphanomyces.
- **CORRECTION — MoM-vs-ξ=0 framing.** The real statistical objection to the current tail is **the
  exponential *shape* assumption (ξ=0)**, not the method-of-moments *estimator* of λ. PBdH (Pickands–
  Balkema–de Haan) gives GPD(σ,ξ) for exceedances; ξ=0 is the special exponential case. Frame the
  approximation as "we assume ξ=0," and `depth_tail` (which is non-parametric in the tail shape — it only
  uses q99.9 + MAD) sidesteps the shape question entirely, another reason it is the better q-driver.
- **VINDICATED — user's Olthof skepticism.** Direct check of Olthof et al. (the "third hybrid intron class"):
  the paper proposes a **continuum**, not concrete discrete clusters; the hybrid class is **hypothesized**,
  from **one unreplicated paper**. The committee's suggestion to model a 3-component mixture is therefore
  **rejected** — we do not bake a contested taxonomy into the adjudicator.
- **CONFIRMED — the no-gains biology (backs objection #2's resolution).** Verified verbatim: U12 introns are
  ancient/single-origin and essentially never gained de-novo (**Physarum** the lone exception = expansion +
  transformed motifs, not de-novo minor-spliceosome gain); loss runs via **deletion ≫ conversion**
  (Lin/Roy ratios ~2:1, 9:1, 5:1). So the "false bearer" failure corner (a small cluster of motif-strong
  introns in a true non-bearer) is biologically near-empty: small-N, motif-strong, *not* a population ⇒
  ancient-loss relics, which are *degrading* (the panel's 6 small-N losses all had `n_call=0`). This
  **de-escalates** the worry that drove objection #2.

**Objection #2 in plain terms (the ELI5, recorded):** the old `excess_z` cutoff is "is this genome's
high-scoring cluster beyond what *this genome's own U2 pile* could throw up by chance?" Bigger genome ⇒ more
U2 draws ⇒ a higher bar (`ln N_u2`) ⇒ a *small* real bearer with few deep U12s can be pushed under the bar
purely because its genome is large — an artifact of asking a size-aware significance question. The fix is to
ask the **size-free** question instead ("is there a real U12 *population* sitting in territory the U2 bulk
doesn't reach?") = `depth_tail`. The biology backstops it: because U12s are never gained, the only way to get
a small motif-strong non-population is ancient-loss relics, which sit *shallow* (tail-adjacent), exactly where
`depth_tail` puts them.

**FINAL NOTE (user, 2026-06-27) — WGD / gene-family duplication weighs on #2.** There IS a mechanism that
*increases* a genome's U12 count without de-novo gains: **whole-genome or gene-family duplication** (e.g.
**salmonids** carry many U12s following a recent WGD). This *strengthens*, not weakens, the #2 resolution:
duplication raises the U12 count by copying *existing, conserved-motif* U12s, so it inflates the **deep,
coherent population** (high `depth_tail`) — never the shallow motif-strong-FP floor. And because WGD scales
the whole intron complement, it tends to raise `N_u2` and the U12 count *together*, so the size-free
`depth_tail` reads it correctly while a size-aware `excess_z` would partly cancel the U12 gain against the
larger U2 denominator. Net: the small-N false-bearer corner stays empty and the genuine-bearer corner
(including duplication-expanded bearers) is exactly what `depth_tail` is built to catch.

### 0e. Independent statistical audit of the SETTLED design + response (2026-06-27)

A fresh independent audit (web-research-backed) of the *settled* design — given the now-known data shape
(clean separation, tight gap, chlamydomonas the irreducible borderline, ensemble-specific calibration).
**Overall: SOUND, ship-able; the architecture and the classification (LCO 24–26/27) are robust because clean
separation pins them. Every real gap is in the *uncertainty* layer (the `P_adj ± CI`), not the point calls.**

- **#1 — degenerate q-CI (CONFIRMED bug → FIXED, commit f4829cb).** The audit verified, and I reproduced, that
  the plain bootstrap-of-median gives a width-0 CI on tied/few calls, silently defeating the UNDETERMINED band
  in the low-N regime it exists for (3 tied borderline calls → q=0.81, CI width 0.00). Fixed with a
  **count-aware smoothed bootstrap** (`ci_smooth_floor_frac`): borderline 3-call CI width 0.00→0.44, count-aware
  (3→0.44, 15→0.20, 60→0.10), chr19 unchanged, deterministic (parity preserved). Honest limit: 3 *deep* tied
  calls stay confident-bearer — that's the biologically-irreducible case (3 deep relics in a loss vs 3 deep
  low-N-bearer U12s are single-genome-indistinguishable; needs conservation/snRNA, not a CI fix).
- **#2 — Firth biases q toward 0.5 → use FLIC (CHECKED → no-op here; NO change).** Empirically the Firth fit's
  mean predicted q (0.467) already equals the panel event rate (0.464); FLIC's intercept correction shifts the
  boundary only +0.016 (2.988→3.003, inside the boundary bootstrap SD 0.14) and barely moves any q. The audit's
  bias result is a *rare-event/imbalance* phenomenon; the balanced ~50/50 panel with the boundary in the gap
  has no meaningful bias. The "loss q up to 0.43" is chlamydomonas, which *should* be ~0.43 (borderline);
  clearly-loss genomes sit at q 0.001–0.1. So the shipped Firth constants stand. **Contingency:** if border-
  genome sourcing skews the panel toward imbalance, apply FLIC in the calibration step.
- **HIERARCHICAL Bayesian two-groups model (the named "cleanest" alternative) — DEFERRED (real downsides).**
  It propagates call-core + boundary + U2-reference uncertainty coherently and gives the EIG stopping rule for
  free, but: a new heavyweight production dependency (PyMC/Stan/numpyro), MCMC/VI determinism vs the
  streaming==in-memory bit-identical requirement, bundling a posterior vs 4 version-pinned constants, and a full
  re-validation. And the audit's own finding: it **does not change the classification**, and the cheap fixes
  capture ~80% of the uncertainty-honesty benefit. Deferred for discussion; worth it only if `P_adj ± CI`
  becomes a load-bearing quantitative downstream output.
- **Lower-priority/noted (not blocking):** `Q99.9(U2)` is the 0.2th-order statistic at the min_u2 floor
  (high-variance, unpropagated reference) — a GPD-based or shrunk high-quantile would stabilize it; the rejected
  pooled-global skipped **partial pooling / empirical-Bayes shrinkage** of each genome's tail toward a global
  hyper-distribution (the established middle that would fix both the floor instability *and* avoid the
  over-suppression) — worth prototyping; the `LOW_N→q_eff=1` default is already a safe operating point. Marginal
  algebra of `P_adj=q·P_motif` is clean; the call double-use (calls define `depth_tail` and receive `P_adj`) is
  real but median-bounded and second-order (disclose, don't re-engineer). Full audit transcript context in the
  session record; reviewed `species_adjudicator.py` + the eval fit scripts.

## 1. The case (why)
- **raw motif features > z** for discrimination, leave-clade-out, with the entire gain on the z-inflated
  loss-species FP class (real-classifier AUC 0.916 vs 0.786; `FINDINGS §6c`).
- **The per-species z-normalizer + mode-sep z-anchor MANUFACTURE the loss FPs** (verified: mode-sep inflates
  motif-empty introns 12–26×, the discount then claws back ~99.9%); the per-species U2 *background*
  correction inside `5'_raw`, by contrast, correctly *deflates* them (`FINDINGS §6e/§6h`). The raw log-odds
  already carry the conserved absolute motif signal, so no per-species feature normalization is needed.
- **Divergent-bearer recall** that mode-sep provides survives as a **gated per-species output threshold**
  (not per-intron z re-anchoring), which cannot inflate motif-empty introns and so needs no discount.

## 2. What's implemented (branch `raw-gated-scoring`)
All changes are **backward-compatible**: `base_features` defaults to the z triple, so existing
modesep/z bundles are byte-identical. **112 classification/scoring/edge unit tests green.**

| area | file | change |
|---|---|---|
| feature | `core/intron.py` | `support2_raw` property (raw analog of `support2`) |
| feature | `classification/trainer.py` | `_extract_feature_vector(base_names=…)` (bundle-driven); `Z_/RAW_BASE_FEATURES`; `SVMTrainer(base_features=…)` threaded to the trained model via `replace` |
| feature | `classification/optimizer.py` | `SVMParameters.base_features` field |
| feature | `classification/predictor.py` | all 3 prediction sites read `base_features` from the bundle |
| gate | `scoring/species_gate.py` | **NEW** — `SpeciesGateParams`, `apply_species_gate` (pure, 11 unit tests), `apply_raw_gated_postprocess` (file-side post-classification) |
| inference | `cli/main.py` | **one branch** in the shared post-classification dispatch: `scoring_mode == "raw_gated"` → run the gate instead of mode-sep + discount |
| bundle | `scripts/build_raw_gated_bundle.py` | builds a raw bundle (raw ensemble via the threaded trainer; `N_MODELS` env) |
| tests | `tests/unit/test_scoring/test_species_gate.py` | gate + `support2_raw` |

**The gate** (`scoring/species_gate.py`): per-species `τ` interpolates strict→relaxed by a product of two
continuous signals from the raw classifier output — `core` (# all-ends-strong high-P introns) × `gap` (median
`support2_raw` of the top-K by P). `g = ramp(core)·ramp(gap)`; `τ = τ_high·(1−g) + τ_low·g`; call iff
`P ≥ τ`. Product ⇒ robust to a few stray strong introns in a loss species.

## 3. Validation (proxy → real classifier → end-to-end)
- **Proxy** (sklearn ensemble on existing score_info): raw+global ≫ z on the z-inflated FP class
  (`PROPOSED_ARCHITECTURE.md`).
- **Real classifier, leave-clade-out** (Stage 1): AUC 0.916 vs 0.786; fidelity gate passed
  (`FINDINGS §6c`).
- **End-to-end, real pipeline** (Stage 2-lite + head-to-head): the raw_gated bundle runs through real
  intronIC; the gate computes per-species and replaces mode-sep+discount.

**Head-to-head (same genome, same `type_id==u12` call definition; 42-model bundle):**

LOSS — FP burden (lower=better), total over tetrahymena/volvox/chlamydomonas/coccomyxa:

| repo default | production v6 (tuned) | raw_gated @τ0.99 | raw_gated @τ0.999 |
|---|---|---|---|
| 40 | **6** | 14 | **6** |

BEARER — conservation recall / total calls @τ0.999: saprolegnia **8/8, 25** (vs default 43); aphanomyces
**7/8, 17** (vs default 49; the 1 miss is the weak-motif tail).

→ **raw_gated's operating curve passes through production v6's loss-FP point (6) with bearer recall preserved
and fewer total calls than the default** — i.e. it matches the tuned production's FP/recall tradeoff with no
discount/graduated-tail/z-normalizer. The earlier "14 vs 6 gap" was an uncalibrated `τ_high` (0.99), not the
architecture. Reproduce: `eval_corpus/headtohead.py` (RAW_TAG=raw42).

### 3a. Principled `τ_high` calibration (plan step 1 — decision-grade)
`τ_high` is the strict threshold applied when `g=0` (no U12 mode), so it **alone** sets the loss-species FP
floor. Calibrated against **two independent negative sources** that the head-to-head test panel does not
include (`eval_corpus/tau_calibration.py`, figure `tau_calibration.png`):

- **(A) deployed raw42 model on the 4 loss GENOMES** (~395k real out-of-distribution true-negatives):
  per-genome FPR.
- **(B) leave-clade-out OOF on the 41k training corpus** (13 clades; serial sklearn ensemble matching the
  bundle's C=200/γ=0.001/rbf/isotonic): out-of-sample U2 FPR + per-clade spread.

| τ_high | (A) loss FP (4 gen) / per-gen FPR | (B) OOF U2 FPR / per-clade spread | OOF U12 recall |
|---|---|---|---|
| 0.99–0.997 | 14 / 3.5e-5 | 6.4e-5 / 0..2.8e-4 | 0.98 |
| **0.999** | **6 / 1.5e-5** | **3.2e-5 / 0..1.4e-4** | **0.96** |
| 0.9999 | 4 / 1.0e-5 | 3.2e-5 / 0..1.4e-4 | 0.93 |

**The two negatives concur:** at 0.999 the OOF U2 FPR (3.2e-5) *predicts* the loss-genome FP
(3.2e-5 × ~100k introns ≈ 3/genome ✓); the curve flattens in the tail; the per-clade spread is tight (no
single clade falls off a cliff at 0.999); and OOF U12 recall is still 0.96. So **τ_high=0.999 ≈ a per-genome
FPR of ~10⁻⁵ — the irreducible motif-strong FP floor — and it is *not* panel-overfit** (confirmed on negatives
the panel excludes). Looser (0.997) triples loss FP for overcall-band recall only; stricter (0.9999) saves 2
FP but erodes OOF recall to 0.93. Caveat: 0.999 sits on the *falling edge* of the absolute loss-FP count
(14→6→4), so the FP **count** has small-number (Poisson) sensitivity there even though the **rate** is flat.

### 3b. Divergent-bearer recall audit (what the conservative point actually drops)
For sap/aph, the introns CURRENT calls but raw_gated@0.999 does not (`eval_corpus/dropped_recall_audit.py`,
`dropped_{saprolegnia,aphanomyces}.tsv`), cross-referenced to IPA conservation class + raw motif vector:

- **Saprolegnia (18 dropped, 0 conserved-U12 lost):** 13/18 IPA-flagged non-U12 (`hard_negative` ×8,
  `U12_fp_suspect` ×3, `focal_not_called` ×2); the 5 unlabeled are all **one-end-strong** (a single U12-shaped
  motif end, the other two weak/negative — the z-inflation FP signature). Current gave nearly all svm 95–99.9.
  **Clean: every drop is an overcall.**
- **Aphanomyces (33 dropped, 1 conserved-U12 lost):** 32 follow the same pattern. The 1 real loss,
  `H257_11044`, is a **genuinely divergent U12** — textbook U12 branch point (`bp_raw=11.31`) + 12 U12
  orthologs, but a 5′SS one substitution off consensus (`GTAACCTT` vs `GTATCCTT`, +4 T→A) → `5'_raw=-1.59` →
  raw-SVM **P=0.90**. It is lost at *any* τ_high (P<0.95) and irrecoverable by relaxation (aph is U12-poor,
  `g=0`). At the single-genome level it is **motif-indistinguishable** from the dropped FPs — only
  cross-species conservation reveals it. Recovering it requires admitting the FP class *or* a conservation
  signal (the known single-genome floor).

**Net:** at the conservative (g=0) point we lose ≈1 divergent U12 per U12-poor bearer — specifically the
degraded-5′SS kind — and in exchange trim the entire one-end-strong overcall band. Consistent with the stated
"rather be conservative on calls than raise FPs" preference. raw_gated also *recovered* 1 aph intron current
missed (`H257_15513`, canonical `GTATCCT` 5′SS + moderate BP).

## 4. Honest verdict + caveats
**Verdict:** supplant is **plausible and well-supported** — the simpler architecture matches tuned-production
FP control with recall held, and `τ_high=0.999` is now **principled, not panel-overfit** (§3a: two independent
negative sources concur on a ~10⁻⁵ per-genome FPR; per-clade-stable; recall-preserving). The audited bearer
losses are ≈1 degraded-5′SS divergent U12 per U12-poor species, irreducible at the single-genome level (§3b).
**Not yet fully earned**, because:
- **The gate's *relaxation* arm (`g>0`) is unexercised by this panel.** Both bearers compute `g=0` (U12-poor:
  top-100 median `support2_raw` ≈1.3–1.6 < gap floor), so their recall is **pure τ_high specificity**, not
  gated relaxation. Validating the novel relaxation needs a **U12-rich-but-shifted** species (the
  Amborella/Oryza class) where conserved U12s sit in the τ_low..τ_high band. ← next calibration target.
- Small panel (4 loss + 2 bearer); independent-negative calibration (§3a) widens it on the negative side only.
- No production-v6 *bearer* recall available → can show raw_gated **holds its own**, not that it **beats**
  production.
- Bundle built via `scripts/build_raw_gated_bundle.py` on the 41k corpus (not `main_train`'s curated ref);
  `meta.iic`/`bed.iic` still carry first-pass calls (only `score_info` `type_id` is rewritten).

## 5. Plan — remaining steps to supplant (in order)
1. **Broader panel + principled `τ` calibration.** ✅ **τ_high calibration DONE (§3a):** 0.999 confirmed
   principled + per-clade-stable on two independent negative sources; recall audit DONE (§3b). **Remaining
   sub-item:** add a **U12-rich-but-shifted bearer** (Amborella/Oryza class) to exercise the gate's `g>0`
   relaxation arm, which the current U12-poor bearers (`g=0`) do not test.
2. **`main_train` integration.** Add a `raw_gated` training mode (skip the normalizer fit, train on raw
   features, stamp `scoring_mode` + gate params) so the bundle is produced reproducibly by the real entry
   point (main.py:5826 → `classify_introns`). Replaces the build-script.
3. **Output sync.** `apply_raw_gated_postprocess` (and the dispatch) update `meta.iic` + `bed.iic` to the gate
   calls, not just `score_info`. Add a streaming==in-memory parity test for the raw_gated path.
4. **Supplant PR.** Once (1)–(3) hold: make the bundled default a raw bundle, flip the default `scoring_mode`,
   and DELETE the z stack — `normalizer.py` adaptive-z, `mode_separation.py` + `mode_sep_pipeline.py`,
   `prior_adjustment.py` (discount). The flag-gated dual path is transitional scaffolding only.

## 6. Artifacts
- Bundles: `eval_corpus/raw_gated.model.pkl` (9-model proof), `eval_corpus/raw_gated_42.model.pkl` (42-model).
- Scripts: `scripts/build_raw_gated_bundle.py`; `eval_corpus/{headtohead.py, raw_gated_v2.py, raw_gated_pipeline.py}`.
- Calibration/audit (§3a/§3b): `eval_corpus/{tau_calibration.py, tau_calibration_plot.py, dropped_recall_audit.py}`;
  outputs `eval_corpus/{tau_calibration.png, tau_calib.log, dropped_saprolegnia.tsv, dropped_aphanomyces.tsv}`.
- Runs: `eval_corpus/stage2_runs/{*_v3ctrl (current arm), *_raw42 (raw arm)}/`.
