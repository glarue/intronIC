# Raw-SVM scaler choice + hyperparameter search — design & plan

Status: **DONE — searched 2026-06-28; verdict = keep current (no rebuild).** Owner: Graham. Created
2026-06-28 after the z-stack removal (`docs/raw_gated_scoring.md`) + the 6-species gut-check.

> ## RESULTS (2026-06-28) — TL;DR: both questions answered, **no bundle change warranted**
>
> Code step 1 (make the scaler swappable + bundle-stamped) is **implemented and shipped** (see §3.1
> below); steps 3–6 (rebuild / re-calibrate / re-stamp) were **evaluated and deliberately NOT done** —
> the search shows the shipped config is already optimal.
>
> **Search**: `eval_corpus/scaler_hpsearch.py` — raw 6-feature single-SVC (isotonic, cv=3) replica,
> leave-CLADE-out on the conservation-anchored eval labels (13 train clades / 17 eval clades; pos=2589,
> neg=982 = ipa + snRNA hard-negatives). Coarse grid `{none,standard,robust} × C{20..1000} × γ{3e-4..1e-2,scale}`
> (90 cfgs) + a refinement `standard × C{200,300,500,700} × γ{2e-4..1e-3}` (20 cfgs). Logs:
> `eval_corpus/scaler_hpsearch.{coarse,refine}.log`.
>
> **(1) Scaler — keep StandardScaler.** RobustScaler never wins (and is catastrophic with γ='scale':
> snRNA@5%→0.28, a RobustScaler-doesn't-unit-variance × fixed-γ mismatch); `none` occasionally edges
> TPR@5% but loses AUC + halves ipa@5%. The fat motif tails do *not* degrade StandardScaler materially.
>
> **(2) C / γ — keep C=200, γ=0.001.** The inherited z-SVM values sit on the raw space's Pareto frontier:
> baseline (standard,200,0.001) = **AUC 0.922, TPR@5% 0.34, ipa@5% 0.09** — tied-best on every *stable*
> metric. AUC spread across the whole refinement is 0.919–0.923 (≤0.004 = single-SVC noise). The only
> metric where baseline looked low (snRNA@5% 0.92) is **noise**: it oscillates non-monotonically across
> adjacent γ (C=200: γ5e-4→0.92, 7e-4→0.98, 1e-3→0.92), i.e. the small snRNA negative set makes the
> single-FPR-point interpolation jumpy — not a real γ effect.
>
> **Decisive tie-break — holdout fidelity is *identical*:** baseline vs the best candidate
> (standard,500,3e-4) both give recall F1@50=**0.999** and consensus_fp FP=**18/454 (4.0%)**.
>
> **Caveat (now closed — see below):** the grid above is a single-SVC cv=3 *proxy*, not the shipped
> 42-model ensemble. So it was confirmed on the production object.
>
> **Ensemble confirmation (`eval_corpus/ensemble_confirm.py`, log alongside).** Each leave-clade-out fold
> trains the REAL pipeline — `SVMTrainer` ensemble of **N_MODELS=42** (shipped count) with stratified U2
> subsampling + isotonic `CalibratedClassifierCV(cv=5)` — and scores the held-out clade with the
> production `ensemble_margin` (mean `decision_function` over all sub-estimators; AUC/TPR are rank-based on
> the margin, no per-fold Platt needed). Baseline vs the two refinement co-leaders:
>
> | config | scaler/C/γ | AUC | TPR@5% | snRNA@5% | ipa@5% |
> |---|---|---|---|---|---|
> | **baseline** | standard/200/0.001 | **0.931** | **0.36** | 0.98 | **0.07** |
> | candA | standard/500/3e-4 | 0.930 | 0.35 | 0.98 | 0.06 |
> | candB | standard/700/2e-4 | 0.930 | 0.35 | 0.98 | 0.05 |
>
> The shipped config is the **outright best of the three** (candidates are marginally *worse*: ΔAUC
> −0.001, ΔTPR@5% −0.01). And the proxy's noisy snRNA@5% (0.92–0.98) collapses to a flat **0.98 for all
> three** in the ensemble — confirming the bias-variance prediction that ensembling erases the
> single-FPR-point jitter (candA/candB's apparent snRNA edge was proxy noise). The ensemble also beats any
> single SVC outright (AUC 0.931 vs ~0.922). **No config beats the shipped one on the production object.**
>
> **Net:** the search **retroactively validates** the inherited hyperparameters on the raw 6-feature
> space (they were never re-searched before; now they have been). The shipped
> `default_pretrained.model.pkl` stays; no re-Platt / re-q / version bump needed. The scaler is now an
> explicit, bundle-stamped, swappable parameter (`SVMParameters.scaler`) so a future re-search is a
> one-line build-env change — but none is indicated today.

---

## 0. TL;DR (and a correction)

While discussing whether to add a global feature scaler for the RBF-SVM, we discovered the
bundled model **already has one**: the sklearn pipeline is `['transform', 'scale', 'svc']`
with a fitted **global `StandardScaler`** (step `'scale'`), so the RBF sees standardized
features. The earlier "we don't do a global scaler" premise was wrong. So the open work is
NOT "add a scaler" — it is:

1. **Is `StandardScaler` the right scaler**, vs `RobustScaler`, given our fat-tailed,
   signal-in-the-tail features? (StandardScaler's mean/std are pulled by the heavy tails.)
2. **Are `C=200`, `γ=0.001` actually optimal for the raw feature space?** They were
   *inherited verbatim from the v3 z-SVM* (the bundle is literally
   `pmotif_adjudicated_v23corpus_C200_g0.001`) and **never re-searched** on the raw 6-feature
   space. The grid-search machinery exists (`SVMOptimizer` / nested CV) but
   `_train_and_save_raw_bundle` hardcodes these.

Settle both with **one nested-CV search over {scaler ∈ Standard/Robust} × C × γ** on the raw
corpus. THEN — the critical part — **re-calibrate the pmotif Platt + q constants**, because
they are fit to the *current* ensemble margins (§4).

---

## 1. Verified current state

Bundle `src/intronIC/data/default_pretrained.model.pkl` (`scoring_mode: pmotif_adjudicated`):

- **Pipeline** (`classification/trainer.py:461-478`, fit by `SVMTrainer`):
  `BothEndsStrongTransformer` → **`StandardScaler()`** → `SVC(kernel=rbf, C=200, gamma=0.001)`,
  wrapped in `CalibratedClassifierCV(method=isotonic, cv=5, ensemble='auto')`.
- The StandardScaler is **global** (fit on training data inside each CV fold; averaged across
  folds/seeds) — NOT per-species. This is the "good half of z" (per-feature balancing for the
  kernel) WITHOUT the per-species re-anchoring we removed.
- **6 features** = base `[5'_raw, bp_raw, 3'_raw]` + extra `[bp_offset, bp_scan_confidence,
  support2_raw]`. `support2_raw` = 2nd-largest of the clipped motif log-odds.
- The bundle's scaler `scale_` ≈ `[25.3, 5.4, 2.58, 3.2, 1.18, 4.86]`, `mean_` ≈
  `[-21.1, 0.0, -0.46, 1.79, 0.54, 2.59]` — i.e. the features really are on disparate scales
  and the scaler is doing real per-feature work.
- **Output columns stay raw** (`5'_raw` etc.) — the scaler lives *inside* the estimator, so it
  never touches the interpretable `score_info.iic` columns.
- Provenance: re-stamp of `eval_corpus/raw_gated_42.model.pkl`
  (`raw_gated_v23corpus_C200_g0.001_isotonic`); `C/γ` inherited from the z-SVM, never
  re-searched on raw.

### Feature-scale facts (human, ~40k introns)
| feature | std | IQR | range | note |
|---|---|---|---|---|
| 5'_raw | 12.4 | 16.8 | [−76, +21] | heavy negative tail; U12 lobe is the positive tail |
| bp_offset | 10.5 | 16.0 | [−46, −9] | integer |
| bp_scan_confidence | 3.6 | 4.6 | [0, 26] | |
| bp_raw | 4.7 | 6.1 | [−29, +14] | |
| 3'_raw | 2.1 | 3.0 | [−7, +7] | |

~6× spread → a single RBF γ cannot balance the features; the in-pipeline scaler does. The
heavy tails are exactly why **RobustScaler may beat StandardScaler** (median/IQR ignore the
−76-type tails that inflate StandardScaler's mean/std).

---

## 2. The questions to settle

1. **Scaler**: `StandardScaler` (current) vs `RobustScaler` (median/IQR, robust to the fat
   tails). Both linear/affine per-feature — faithful geometry. (Reject MinMax/MaxAbs:
   outlier-fragile; reject `Normalizer`: per-sample L2, nonlinear across features, wrong tool.)
2. **C / γ**: re-tune on the raw + chosen-scaler space (inherited z-values are likely
   suboptimal but may be fine — the search tells us). Also reconsider `calibration_method`
   (isotonic vs sigmoid) while we're at it.

Hard constraint: the scaler stays **global / frozen** (fit once on the training corpus),
**NEVER per-species** — per-species refit is the z-inflation we removed (see
`docs/raw_gated_scoring.md`, eval-corpus z-vs-raw findings).

---

## 3. Implementation plan (staged)

### 3.1 Implemented (step 1 — shipped 2026-06-28)

The scaler is now an explicit, swappable, bundle-stamped parameter; steps 2–6 ran as the *evaluation* in
the RESULTS box above (verdict: keep current → no rebuild).

- `svm_factory.make_scaler_step(name)` → `StandardScaler` / `RobustScaler` / `None` (+ `SCALER_CHOICES`).
- `SVMParameters.scaler: str = "standard"` (frozen-dataclass field, **appended last** + custom
  `__getstate__`/`__setstate__` mirroring `SVMModel`, so existing pickled bundles back-fill
  `scaler='standard'`). **This pickle-compat shim was required**: a mid-list insert into the `slots=True`
  positional pickle silently shifted every later field and dropped `calibration_isotonic_logloss`,
  breaking inference on the shipped bundle — caught by `test_streaming_equivalence`. The shim also
  hardens *all* future field additions to `SVMParameters`.
- `SVMTrainer._train_single_model` builds the `'scale'` step from `parameters.scaler` (drops it for
  `none`); `SVMOptimizer(scaler=...)` threads it through all 3 search pipelines + stamps the returned
  `SVMParameters`. Scaler lives *inside* the estimator, so `score_info.iic` columns stay raw and the
  predictor/adjudicator need no change (margin = `decision_function` through the scaler).
- `scripts/build_raw_gated_bundle.py` reads `SCALER` / `SVM_C` / `SVM_GAMMA` env (defaults =
  standard/200/0.001) and stamps them into `model_id` + `provenance`.
- Tests: `tests/unit/test_classification/test_trainer.py::TestScalerSwap` (6). Default `standard` keeps
  the bundle byte-identical (`test_streaming_equivalence` 15/15; full unit suite 608 pass).

### 3.2 Original staged plan (for reference)

1. **Make the scaler swappable in `SVMTrainer`.** Add a `scaler` knob (e.g.
   `scaler ∈ {"standard", "robust", "none"}`) threaded onto `SVMParameters` (frozen dataclass)
   + stamped onto the bundle, mirroring how `base_features` was made swappable in 2b-prep.
   Pipeline picks `StandardScaler()` / `RobustScaler()` / drops the `'scale'` step accordingly.
   Default stays `standard` so existing behaviour is byte-identical.
2. **Search harness** (in `eval_corpus/`, extend `path1_calibrated_test.py` /
   `stage1_real_classifier_test.py` / `q_firth_leaveclade.py` machinery): nested **leave-clade-out**
   CV over the grid `{scaler} × {C} × {γ}` (+ maybe `calibration_method`), scored with the
   **FP-weighted** metric we've used (TPR@low-FPR on the conservation-anchored eval labels),
   plus the bundle's held-out recall set + `consensus_fp` set. Faithful replica = the real
   6-feature set + isotonic + the 41k-row `multispecies_v23` / `stage1_trainmatrix` corpus.
3. **Build the winning raw ensemble** via `scripts/build_raw_gated_bundle.py` with the chosen
   `scaler / C / γ` (N_MODELS=42, 3 seeds, as the canonical bundle).
4. **Re-calibrate the pmotif layer** (§4) — Platt + q — and bump
   `ADJUDICATOR_PARAMS_VERSION`.
5. **Re-stamp** the pmotif bundle (`scripts/build_pmotif_adjudicated_bundle.py`) and swap it in
   at `src/intronIC/data/default_pretrained.model.pkl`.
6. **Validate** (§5).

---

## 4. CRITICAL downstream dependency — re-calibration (do not skip)

The pmotif scoring layer is fit to the *current* ensemble's margin distribution:

- `P_motif = σ(PLATT_A·margin + PLATT_C)` with `PLATT_A=2.796, PLATT_C=-1.178`
  (`scoring/species_adjudicator.py:45-46`).
- `q = σ(Q_A·depth_tail + Q_B)` with `Q_A=3.64, Q_B=-10.86` (Firth fit; `:47-48`).
- `depth_tail` itself is computed from the ensemble call-margins vs the species U2 margin tail.

**Any change to the scaler / C / γ shifts the ensemble margins**, so all of the above become
stale. The full chain is therefore:

  search → new raw ensemble → **re-fit Platt** (`eval_corpus/calibrate_pmotif_bundle.py`) →
  **re-fit q** leave-clade-out (`eval_corpus/q_firth_leaveclade.py`) → bump
  `ADJUDICATOR_PARAMS_VERSION` → re-stamp bundle → validate.

A new SVM with stale Platt/q would silently mis-score. This coupling is the main reason this is
a multi-step, carefully-tested change and not a one-line scaler swap.

---

## 5. Test / validation plan

- **Search-level**: held-out leave-clade-out FP-weighted metric vs the current
  (Standard/200/0.001) baseline — accept only if ≥ baseline on the conservation eval AND no
  regression on the loss-species FP control.
- **End-to-end panel** (the gut-check baseline to beat / match): re-run the 6 species and check
  U12 calls stay sane — human ~700–800, danio rich, nematostella rich, arabidopsis ~290–323,
  drosophila ~20, **c. elegans 0** (loss control). Use the corrected metrics (post the
  count-consolidation fix; `scripts/patch_metrics_counts.py` if comparing old files).
- **chr19 smoke** + **streaming == in-memory parity** (`tests/integration/test_streaming_equivalence.py`)
  — must stay bit-identical between paths.
- **Metrics consistency asserts** (already in `_finalize_classification_metrics`): HC ≤ u12,
  boundaries ⊆ counts.
- **Labeling invariant** (verified post-fix, keep as a check): `type_id==u12 ⇔ adjusted_score≥50`,
  `rel_score>0 ⇔ adjusted_score≥90`, `meta == score_info` type_id.

---

## 6. Notes / context for a fresh session

- This builds on the **z-stack removal** (supplant 2a + 2b, all committed on branch
  `raw-gated-scoring`, not pushed; tag `pre-zstack-removal` preserves legacy bundles) and the
  **count-consolidation** fixes. `docs/raw_gated_scoring.md` is the authoritative architecture
  doc; this file is the scaler/HP-search follow-up.
- The "+0.006 AUC from adding StandardScaler" in the eval-corpus `control_scaling_test` was a
  bare single-SVC *gaining* the scaler the production pipeline already has — i.e. production
  already realizes it. The remaining gain (if any) is Robust-vs-Standard + re-tuned C/γ.
- Why no scaler choice was ever made deliberately: the z-SVM's features were pre-scaled
  per-species (so the in-pipeline StandardScaler was near-redundant there); porting to raw kept
  the StandardScaler + the z-era C/γ without re-searching.
- Build/calibration scripts: `scripts/build_raw_gated_bundle.py`,
  `scripts/build_pmotif_adjudicated_bundle.py`;
  `eval_corpus/{calibrate_pmotif_bundle,q_firth_leaveclade,path1_calibrated_test,stage1_real_classifier_test}.py`.
- Training corpus: `multispecies_v23` / `eval_corpus/stage1_trainmatrix.npz`.
