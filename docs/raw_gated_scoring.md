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

Deferred (airtight pass, §0a): leave-clade-out OOF *margins* (distinct from the q leave-clade-out — this
re-scores `P_motif` itself out-of-sample); POT-GPD as a one-time tail diagnostic. **Still ahead:** the
eval_corpus calibration pipeline (Platt-OOF + q-panel for the bundled ensemble — the real prerequisite for a
shippable default), a committed pmotif integration/parity test (needs a lightweight bundled fixture), the §5
formal gates, and the step-7 supplant PR (flip default + delete the z stack). The low-N/discrete-input
fallback refinement (a bundled global-`q` or snRNA backstop instead of `q_eff=1`) is the other open item.

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
