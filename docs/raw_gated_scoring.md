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

## 4. Honest verdict + caveats
**Verdict:** supplant is **plausible and well-supported** — the simpler architecture matches tuned-production
FP control with recall held. **Not yet earned**, because:
- `τ_high=0.999` was picked post-hoc to land on production's FP count. It is *recall-preserving* (recall is
  identical at 0.985 and 0.999), but the default `τ` must come from a **principled calibration on independent
  negatives**, not the test losses.
- Small panel (4 loss + 2 bearer).
- No production-v6 *bearer* recall available → can show raw_gated **holds its own**, not that it **beats**
  production.
- Bundle built via `scripts/build_raw_gated_bundle.py` on the 41k corpus (not `main_train`'s curated ref);
  `meta.iic`/`bed.iic` still carry first-pass calls (only `score_info` `type_id` is rewritten).

## 5. Plan — remaining steps to supplant (in order)
1. **Broader panel + principled `τ` calibration.** Add more loss + bearer species; set `τ_high` from
   independent negatives (e.g. a target FPR on training U2 / held-out loss species), confirm ~0.999 is robust
   and not panel-overfit. Re-run `headtohead.py`. *(This is the decision-grade test.)*
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
- Runs: `eval_corpus/stage2_runs/{*_v3ctrl (current arm), *_raw42 (raw arm)}/`.
