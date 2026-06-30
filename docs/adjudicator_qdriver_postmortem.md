# Post-mortem: the `depth_tail → q → adjusted_score` adjudicator (and why it's being replaced)

**Status:** the method described here is what currently ships in `scoring/species_adjudicator.py`
(`ADJUDICATOR_PARAMS_VERSION = "depth_tail_firth_2026-06-27"`). This document records *why it is being
superseded* by a two-number design, so the failure modes are not re-introduced. Written 2026-06-29 from
the bake-off + leave-clade-out + two-groups analyses in `eval_corpus/` (`tail_stat_bakeoff.py`,
`refit_q_2feature.py`, `anchor_lco.py`, `twogroups_posterior_demo.py`). Cross-ref the
`adjudicator-qdriver-bakeoff` memory.

---

## 1. What the method was

Per-intron motif probability `P_motif = σ(2.796·margin − 1.178)` (calibrated, species-agnostic). On top of
it, a per-species overlay:

- `depth_tail = (median(call-margins) − U2 q99.9) / MAD_U2` — a size-invariant separation statistic
  (call set = `P_motif ≥ 0.9`; U2 ref = `P_motif < 0.5`).
- `q = σ(3.64·depth_tail − 10.86)` — a **Firth-penalized logistic** fit to bearer/loss species labels;
  `q = 0.5` lands at `depth_tail = 2.98`, inside the panel gap `[2.91, 3.39]`.
- `P_adj = q · P_motif`; `adjusted_score = 100·P_adj`; `type_id = u12` iff `adjusted_score ≥ 50`;
  **high-confidence U12 (HC)** iff `adjusted_score ≥ 90` — the number reported per species.

Because for a confirmed bearer `P_motif ≈ 1` on its calls, **HC ⟺ q ≥ 0.9 ⟺ depth_tail ≥ 3.59**.

## 2. The failure

**Genuine bearers reported 0 high-confidence U12.** On the labeled panel, the depth_tail driver put only
**35/43 bearers** over the HC line. The casualties were exactly the band-dwellers:

- *Neocallimastix* (167 strong calls, q 0.67 → **HC 0**), *Anaeromyces* (119 calls, q 0.73 → **HC 0**),
  *Polychytrium* (q 0.80 → HC 0) — all called bearers (q > 0.5), none confident.
- *Chytridium lagenaria* (real bearer, u11+u12 present): `depth_tail = 1.96`, **below** the loss
  *Chlamydomonas* (2.91) — depth_tail mis-ranks it as more loss-like than a loss.

The reported deliverable (`#(adjusted_score ≥ 90)`) was wrong for the cases that mattered most.

## 3. Why it failed — four intertwined issues

### 3a. The evaluation gap (root cause)
Every development check thresholded on **`P_motif ≥ 0.9`** (motif probability) or evaluated the species
fit at **`q > 0.5`** (the call boundary). **Nobody computed the actual reported metric** —
`#(adjusted_score ≥ 90)` per species — against per-species ground truth. The whole apparatus (FP-inflation
control, bearer-precision, 37-genome panel separation, leave-clade-out, the 8-agent committee, the
independent §0e audit) validated *separation* and the *call*, one full threshold below where the report
lives. Review 02 caught the *precursor* (q's magnitude is unidentified under separation; bearers span
0.77–1.00) but framed it as an uncertainty-layer defensibility issue and never traced it to HC.
**Lesson: validate the exact quantity you publish, against truth — not a proxy threshold.**

### 3b. Threshold-baking under separation (the q is not a calibrated probability)
The panel is cleanly separated, so a plain logistic MLE diverges; Firth keeps it finite. Under separation
the fit splits in two:
- the **0.5 *location*** (the ratio −Q_B/Q_A) is well-identified — it sits robustly in the gap;
- the **slope** (Q_A) is **penalty-determined, not data-determined** — and the slope fixes where *every
  other* threshold lands relative to 0.5.

So the fit legitimately anchors **exactly one** threshold (0.5) and lets the placement of all others —
**including the reported 0.9** — be an artifact of the penalty. That is why bearers were fine at q>0.5 and
vanished at q≥0.9: 0.5 was anchored in the gap; 0.9 was read off an arbitrary slope (landing at
depth_tail 3.59, above the weakest bearers at 3.39).
**General principle: a classifier fit to separated binary labels can honestly anchor one threshold.
Whichever you tune (50% here) becomes the boundary; every other threshold is a slope-dependent,
under-determined byproduct. Building explicitly around 90% would have the identical defect, mirror-imaged
(it would relocate the arbitrariness to 50%). Do not read a second threshold off a fitted-on-separated
probability.**

### 3c. `depth_tail` is count-blind
`depth_tail` is a *central-tendency* statistic (median call-core). It was chosen for size-invariance — but
size-invariance was achieved by **discarding the count** of strong calls. That fails *numerous-but-shallow*
divergent bearers: *Chytridium* has 33 strong calls (a real population) but a shallow call-core, so its
median-based depth (1.96) reads loss-like. The bake-off (`tail_stat_bakeoff.py`, 43 cross-clade bearers +
14 losses + the 840k *Symbiodinium* probe) is decisive:

| statistic | bearer/loss separation | divergent bearers | huge-loss (symbio) |
|---|---|---|---|
| `depth_tail` (median) | ❌ OVERLAP (Chytridium 1.96 < chlamy 2.91) | ❌ | ✓ |
| `excess_z` (older, size-aware EVT max) | ❌ OVERLAP (Chytridium −0.29) | ❌ | ✓ |
| `n_call` (raw count) | ❌ OVERLAP | ✓ | ❌ inflates (181 calls → false bearer) |
| `density` (count/N) | ❌ OVERLAP (small-loss spike) | ✓ | ✓ |
| **`z_excess`** (Poisson sig. of strong-call **count** vs the genome's own U2 tail) | ✅ **AUC 1.000, gap +1.4** | ✅ | ✅ (−3.2, no inflation) |

`z_excess` is the "is there a real U12 *population*" statistic done right: a count, but referenced to the
genome's own U2-tail expectation, so it is robust at both size extremes where raw `n_call` and `density`
fail. **Lesson: "does this genome have a population of X" is a count question; a central-tendency
statistic answers the wrong question.**

### 3d. A species signal cannot be folded into a calibrated per-intron probability
`P_adj = q · P_motif` multiplies a calibrated per-intron probability by a per-species near-binary gate and
thresholds the product. The product is not a calibrated probability of any nameable event, so thresholding
it at 90% (or 75%) is not a principled certainty choice. The deeper reason (not a fixable bug):

- **Calibration needs a populated transition (overlap).** "q = 0.7 is calibrated" requires a population of
  genomes where ~70% are bearers. U12 status is **near-binary** (ancient, single-origin, no de-novo gain),
  so species sit at the two ends — there is no continuum of "60%-bearer species" to calibrate against.
- **Each species is n = 1.** A novel/divergent genome's bearer-ness can't be a calibrated frequency — there
  is one of it. Its honest uncertainty is *epistemic*, not a frequency.
- The rigorous combination (a two-groups / local-fdr posterior using the species **prevalence** as the
  prior) *exists*, but its honest output is **interval-valued**: `twogroups_posterior_demo.py` shows the
  per-intron posterior band stays tight for a confident bearer (human, 0.89 ± 0.005), balloons for a
  band-dweller (Chytridium, 0.53 ± 0.14), and collapses to 0 for losses. It also reduces to ≈`P_motif`
  for confident bearers — and §0 already found the formal prior-shift *hurts* recall (623→565). So folding
  the species number in buys little where it works and fabricates precision where it doesn't (the q-ceiling
  band was exactly that fabrication).

**Lesson: keep a non-calibratable signal out of a calibrated quantity. Report the two separately.**

## 4. What replaces it (the two-number design)

1. **Per-intron: report `P_motif`** — calibrated on the per-intron continuum (LCO ECE ~7e-4), so a user
   thresholds it **post-hoc** at any certainty (0.75 / 0.90 / 0.99), each principled. This *preserves* the
   cross-species-comparable threshold v1 wanted (P_motif is species-agnostic) — and improves on it
   (v1's comparability came from z-normalization, which broke for divergent species like Amborella/Oryza).
2. **Per-species: an empirical gate** — `z_excess` (the population statistic) placed against an
   **empirical bearer/loss gap** (`loss_ceiling`, `bearer_floor`; supervised, frozen). Tier =
   confident-bearer / **suspect** / confident-loss, with the continuous `z_excess` ± bootstrap CI reported
   alongside. No fitted probability slope ⇒ no "one anchored threshold, rest arbitrary" pathology — a
   separated classifier is asked for the *one* thing it can give (a boundary), nothing more.
3. **Conservatism via out-of-support abstention** (motif-only): a confident call requires being inside the
   calibration support; the gap and beyond-the-extremes abstain to **suspect**. A future bearer more
   extreme than anything sampled lands in the gap and abstains by construction, rather than mislabeling.
4. **Scope discipline:** the pipeline is **motif-only at runtime** — no phylogeny, no snRNA. Those signals
   belong to the *database layer* downstream. Taxonomy is used only at *calibration time* (coverage audit /
   data-generation targeting), never in the classifier.

Validation so far: `z_excess`-gate leave-clade-out (`anchor_lco.py`) = **0 cross-errors** (held-out
genomes never mis-cross; worst case is an honest "suspect"). The binding open work is *coverage* of the
divergent-bearer tips that pin `bearer_floor` (early fungi, oomycetes, basal metazoans), not the method.

## 4a. Frozen calibration (2026-06-30)

Validated end-to-end on cross-clade data (the 43-genome panel + 68 freshly-run divergent bearers
[early fungi / oomycetes / basal metazoans] + 9 snRNA-confirmed losses; analysis in
`eval_corpus/{corroborate_v2.py, refit_anchors_expanded.py, batch_zexcess_68.py}`).

**snRNA bearer/loss rule (calibration-time only; never at runtime):**
- `significant` = best cmsearch E ≤ **0.01** (Infernal `--incE` inclusion default — NOT a stricter cut;
  divergent snRNAs sit in the [8e-14, 1.1e-3] void and a 1e-6 cut wrongly drops them).
- **BEARER** = `n_sig ≥ 3` of {u11,u12,u4atac,u6atac}. This is *self-consistently defining-aware*: a loss
  retains only the conserved u4atac/u6atac (caps at 2/4 at 0.01), so reaching 3 *requires* a defining
  snRNA (u11 or u12) — **without** explicitly demanding the u11+u12 pair, which would wrongly reject the
  flagship gut fungi (neocallimastix/anaeromyces/piromyces have strong u12 but *undetectable, divergent*
  u11: u12 ≤2.5e-7, u11 0.1–2.2). Don't require the pair.
- **DIVERGENT_BEARER (rescue)** = `n_sig ≥ 1` AND multiple borderline hits AND **a clear strong-U12 motif
  population (high z_excess)**. Both required — snRNA borderline *and* motif corroboration. Without the
  population requirement, a 4-weak-hits Mortierellomycotina (`mortierella sp.`) would false-promote at 0.1.
- **LOSS** otherwise (incl. 2-conserved-only, like aspergillus coremiiformis which keeps u4atac/u6atac but
  lost u11+u12). **UNKNOWN/REVIEW** if not cmsearched, or a threshold-edge case with no population.
- **High-z bearer corroboration is asymmetric:** the motif population is primary; *any* snRNA corroborates;
  the only red flag is searched-and-found-ZERO. Result: **0 suspicious of 96 bearers (z≥4)**, 82 corroborated,
  14 obvious bearers (corals/sponges/Hydra/Adiantum/AMF, z≥100) pending cmsearch.

**Frozen species-gate anchors (z_excess):** `loss_ceiling = 2.60` (aspergillus coremiiformis; all loss-anchors
confirmed loss at 0.01, incl. the textbook 0/4 ciliates/apicomplexan/green-algae), `bearer_floor = 4.00`
(*Mycotypha*, snRNA-confirmed, lowest motif-detectable corroborated bearer). Gap **[2.60, 4.00] CLEAN**;
**leave-clade-out = 0 cross-errors / 87**. Adding the 68 divergent genomes did NOT crash the floor — the low
ones were snRNA-confirmed losses (clade→bearer was too coarse; e.g. Mortierellomycotina are real U12 losses)
or **motif-undetectable confirmed bearers** (`plasmopara halstedii`, `pythium insidiosum` — oomycetes, ≥3 snRNA,
z<2.6) that define the irreducible motif-detection floor and are correctly deferred to the snRNA/database layer.

**Output = two numbers:** per-intron `P_motif` (calibrated, post-hoc-thresholdable at any certainty) + the
species `motif_category` gate {`DETECTED` / `BORDERLINE` / `NOT_DETECTED` / `UNASSESSABLE`} with `z_excess ± CI`.
The category names the **motif evidence** for a U12 population, NOT biological bearer/loss truth — `DETECTED` ==
'a U12-motif population, by motif, corroborate downstream'; `NOT_DETECTED` is **not** a loss call (motif-silent
divergent bearers — plasmopara/pythium — land here). The corroborated species verdict (`u12_status` =
U12_POSITIVE / U12_NEGATIVE / CONFLICT = motif × snRNA) is a separate database-layer concept. Out-of-support
abstention is motif-only; taxonomy is calibration-time bookkeeping; snRNA/phylogeny live downstream.

**Open TODOs:** cmsearch the 14 obvious bearers (won't move the floor — all z≥6); vet the 2 in-gap species
(`mucor_saturninus` 2.70 / `parasitella` 2.72; could nudge loss_ceiling→2.72, gap stays clean); route
`mortierella sp.` → REVIEW.

**IMPLEMENTED 2026-06-30** in `scoring/species_adjudicator.py` (`ADJUDICATOR_PARAMS_VERSION =
zexcess_gap_2026-06-30`): `compute_z_excess` + `classify_motif_category` (frozen anchors) are the PRIMARY gate,
emitting the `MotifCategory` enum {`DETECTED` / `BORDERLINE` / `NOT_DETECTED` / `UNASSESSABLE`}; `_effective_q`
is now a BINARY gate (suppress only `NOT_DETECTED`, no q-deflation); score_info emits `z_excess` +
`motif_category`; `type_id` = u12 iff `P_motif>=0.5` AND not `NOT_DETECTED`; `rel_score` = `100*P_motif-90`
(so `high_confidence_u12` = strong-motif calls in a non-`NOT_DETECTED` genome, qualified by `motif_category`);
the `depth_tail->q->P_adj` chain is retained as labelled secondary/back-compat. The metrics summary surfaces
`motif_category` + `z_excess` and a `confident_u12_motif` count (= strong calls where `motif_category==DETECTED`).
Validated: adjudicator + metrics/writers unit tests + streaming/in-memory & pmotif PARITY all green
(bit-identical across paths). Consumer (`get_iic_stats.py`) migrated to read `motif_category` +
`confident_u12_motif` (legacy keys kept as fallback). Remaining (separate repo): migrate the WtMTA pipeline.

## 5. Generalizable lessons (the point of this document)

1. **Grade the deliverable, not a proxy.** Compute the exact published metric against ground truth. A green
   board of separation/precision/AUC checks said nothing about the HC count we actually reported.
2. **A classifier on separated binary labels anchors exactly one threshold.** Don't turn it into a
   probability and read other thresholds off the (unidentified) slope. Tuning 50% vs 90% is the same trap.
3. **"Does this have a population of X" is a count question.** A central-tendency statistic (median depth)
   is the wrong tool; a background-referenced count (`z_excess`) is the right one.
4. **Don't fold a non-calibratable signal into a calibrated number.** A near-binary, n=1-per-unit species
   signal cannot become a calibrated per-intron probability; report the two values separately and let
   thresholding happen post-hoc on the calibrated one.
5. **Keep runtime scope tight.** Calibration may use rich external evidence (snRNA, taxonomy); the runtime
   classifier should depend only on its own input. Conservatism comes from *abstention within support*, not
   from importing metadata.
