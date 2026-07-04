# The species adjudicator — current design (authoritative)

**Status: SHIPPED.** `ADJUDICATOR_PARAMS_VERSION = "zexcess_gap_pgumbel_cs_2026-07-03"`.
Code: `src/intronIC/scoring/species_adjudicator.py`. Bundle: `src/intronIC/data/default_pretrained.model.pkl`.

This is the single current source of truth for the **output-level species adjudicator** — the layer that sits
on top of the raw-feature `P_motif` SVM ensemble and makes the per-species U12-motif call. It consolidates and
supersedes the design narrative in `raw_gated_scoring.md §0a–§0e` (the `depth_tail → q → P_adj` design, now
historical) and folds in the two refinements that shipped after it. For the *why* behind each change, the two
rationale docs remain authoritative on their own subject and are linked at each step:

- **`adjudicator_qdriver_postmortem.md`** — why `depth_tail → q → P_adj` was abandoned and the two-number
  `z_excess` gap gate replaced it (the four generalizable lessons).
- **`adjudicator_strength_evt.md`** — the per-genome Gumbel strength gate (`p_gumbel_p95`) that recovers
  divergent low-count bearers, with the absolute `cs_p95` retained as a null-free co-fallback.

---

## 1. What it produces — two numbers, reported separately

The adjudicator deliberately emits **two independent quantities** and never fuses them into one thresholded
score (that fusion was the core mistake of the superseded design — see the postmortem, lesson 4):

1. **`P_motif`** — the per-intron, **species-agnostic**, Platt-calibrated motif probability
   `P_motif = σ(2.796·margin − 1.178)`, where `margin` is the raw-feature SVM ensemble's
   `decision_function` (use the margin, not the isotonic prob, which saturates at 0/1). A textbook `GTATCC`
   U12 scores ~0.99 regardless of how many siblings its genome has. Calibrated on the per-intron continuum
   (leave-clade-out ECE ~7e-4), so a user thresholds it **post-hoc** at any certainty (0.75 / 0.90 / 0.99),
   each principled. This *is* the cross-species-comparable score v1 wanted from z-normalization — without
   z-normalization (which broke for divergent species like Amborella/Oryza).

2. **`motif_category`** — a per-species gate over the genome's `P_motif` calls: one of
   **`DETECTED` / `INCONCLUSIVE` / `NOT_DETECTED` / `UNASSESSABLE`**, reported with the continuous
   `z_excess` (± bootstrap CI) alongside. It names the **motif evidence for a U12 population**, *not*
   biological bearer/loss truth — see §5.

## 2. `z_excess` — the population statistic

The gate's primary driver is **`z_excess`**: the Poisson significance of the *count* of strong calls
(`P_motif ≥ 0.9`) against what the genome's **own U2 background tail** predicts. Fit a Gumbel EVT null to the
genome's U2 margin tail (`_u2_reference`: `med`, `mad`, `z_expmax`, `evt_beta`), extrapolate to the call
region, and ask: are there significantly *more* high-margin introns than U2's own tail would throw up by
chance? A loss's motif-strong FPs **are** the U2 tail (`z_excess ≈ 0`); a bearer has a separate population U2
can't explain (large `z_excess`).

`z_excess` is "does this genome have a U12 *population*" done right. The postmortem's bake-off (§3c there) is
why it beat every alternative: `depth_tail`/`excess_z` (central-tendency) overlap on numerous-but-shallow
divergent bearers; raw `n_call` inflates on huge genomes; `density` spikes on tiny losses. Referencing the
count to the genome's *own* U2 tail is robust at both size extremes (`z_excess` AUC 1.000, gap +1.4 on the
bake-off panel; Symbiodinium's 181 calls → −3.2, no false inflation).

## 3. The gate — `classify_motif_category` (exact logic)

```
classify_motif_category(z_excess, cs_p95, p_gumbel_p95, n_calls, params):
  if not finite(z_excess):                    -> UNASSESSABLE   # no per-genome U2 null (LOW_N / unfit tail)
  if z_excess >= bearer_floor_z (5.50):       -> DETECTED       # count gate: a clear U12 population
  if cs_gate_enabled and n_calls >= cs_min_calls (3) and (
        (finite(p_gumbel_p95) and p_gumbel_p95 <= p_gumbel_threshold (0.01))   # PRIMARY strength gate
     or (finite(cs_p95)       and cs_p95       >= cs_point_threshold (5.0)) ): # CO-FALLBACK (null-free)
                                              -> DETECTED       # strength gate: few-but-genuinely-strong calls
  if z_excess <= loss_ceiling_z (2.60):       -> NOT_DETECTED   # calls consistent with the U2 background
  else:                                       -> INCONCLUSIVE   # the gap [2.60, 5.50): abstain, corroborate
```

**Ordering matters:** the strength gate is checked **before** the loss ceiling, so a genome whose `z_excess`
sits *below* 2.60 (count-wise loss territory) but which has a few genuinely strong calls is **rescued to
DETECTED**. This is the divergent-bearer recovery path (worked example: `achlya_hypogyna`, an oomycete with
`z_excess = 2.13` but `cs_p95 = 7.3` → DETECTED; snRNA-corroborated bearer, all four minor snRNAs at E ≤ 1e-9).

### The strength gate (`p_gumbel_p95`) — why it exists

`z_excess` is a **count** statistic → blind to a divergent bearer with only a *few* genuinely strong U12
calls (low count → low `z_excess`, reads as loss). The strength gate fixes this by evaluating the **same
per-genome Gumbel U2 null at the call upper tail**: `p_gumbel_p95 = P(the genome's own U2 background produces
a max ≥ the call cs_p95)`, where `cs_p95` is the 95th percentile of the *un-clipped* ensemble margins on the
calls. `p_gumbel_p95 ≤ 0.01` = "these strong calls are a ~1% outlier against this genome's own U2 noise" →
DETECTED. It is composition-adaptive (judged against each genome's own `med`/`mad`/`z_expmax`), which is why
it replaced the earlier tuned constant `cs_p95 ≥ 5.0` as the primary — that 5.0 was a per-genome ~0.5% outlier
test *all along* (`cs_p95 = 5.0 ↔ p_gumbel_p95 ≈ 0.005`), now stated as one and adaptive where a global
constant would drift. `cs_p95 ≥ 5.0` is retained as a null-free **co-fallback** (OR) for marginal EVT fits;
the two agree on the panel, so this is **no-regression by construction** (worst loss `p_gumbel_p95 = 0.044`,
well above 0.01 → no new loss FP). Upper-tail evaluation is the unlock: recovery climbs monotone
median 12/46 → p95 36/46 and saturates at p95. Full evidence: `adjudicator_strength_evt.md`.

## 4. How the gate turns into calls

Downstream labeling uses a **binary** effective gate (`_effective_q`): suppress **only** `NOT_DETECTED`.

- `q_eff = 0` iff `motif_category == NOT_DETECTED`, else `q_eff = 1`.
- `P_adj = q_eff · P_motif`; `adjusted_score = 100·P_adj`; `rel_score = 100·P_adj − 90`.
- **`type_id = u12`** iff `P_motif ≥ 0.5` **AND** `motif_category != NOT_DETECTED`.
- **High-confidence (HC)** = `rel_score > 0` ⇔ `P_motif > 0.9` in a non-`NOT_DETECTED` genome.
- `confident_u12_motif` (metrics) = strong-motif calls (`rel_score > 0`) in a **DETECTED** genome.

Consequences worth internalizing:
- **`INCONCLUSIVE` and `UNASSESSABLE` do NOT suppress calls.** They flag the *species-level* motif evidence as
  ambiguous, but their strong-motif introns still get `type_id = u12` / HC (`q_eff = 1`). The category is the
  confidence qualifier, not a per-intron veto — corroborate INCONCLUSIVE genomes downstream (snRNA/phylogeny).
- Only `NOT_DETECTED` zeroes the calls (`P_adj = 0`, `type_id = u2`, `rel_score = −90`).

> **Column note.** The `q`, `P_adj`, `P_adj_lo`, `P_adj_hi`, `adjusted_score` columns are **legacy** — the
> re-derived (binary-`q`) remains of the superseded `q·P_motif` chain, kept so existing consumers don't break.
> `P_adj_lo == P_adj_hi == P_adj` now (the species uncertainty lives in `motif_category`, not a per-intron CI).
> **New consumers should read `P_motif` + `motif_category`.**

## 5. What the categories mean (and do not)

`motif_category` names the **motif evidence for a U12 population**, not biological bearer/loss truth:

- **`DETECTED`** = "a U12-motif population is present, by motif — corroborate downstream." *Not* a proof of
  functional U12.
- **`NOT_DETECTED`** = "the strong calls are consistent with this genome's U2 background." **This is NOT a loss
  call.** Motif-silent divergent bearers (e.g. the oomycetes `plasmopara halstedii`, `pythium insidiosum` —
  ≥3 snRNAs but `z_excess < 2.6`) land here and are correctly deferred to the snRNA/database layer.
- **`INCONCLUSIVE`** = "in the gap [2.60, 5.50) and not strength-rescued" — abstain within support; corroborate.
- **`UNASSESSABLE`** = "cannot reference a per-genome U2 null" (LOW_N / degenerate tail); see §6.

The corroborated species verdict (`u12_status` = U12_POSITIVE / U12_NEGATIVE / CONFLICT = motif × snRNA) is a
separate **database-layer** concept. **Runtime scope is motif-only**: no phylogeny, no snRNA. Taxonomy and
snRNA are used only at *calibration time* (coverage audit, anchor sourcing), never by the runtime classifier.
Conservatism comes from **abstention within support** (the gap and the `UNASSESSABLE` band), not from importing
metadata. A future bearer more extreme than anything sampled lands in the gap and abstains by construction
rather than being mislabeled.

## 6. `-q` / low-N fallback

The per-genome test presumes the genome *has* a U2 background to reference:

- **Full-complement genome / `-q`** (≥ `min_u2 = 200` U2 introns) → EVT tail fits → normal adjudication
  (135/135 real genomes fit their tail).
- **Curated / tiny `-q`** (< 200 U2) → `LOW_N` → `motif_category = UNASSESSABLE`, `q_eff = 1`: per-intron
  `P_motif` is reported as-is, no bearer/loss call — the honest answer (no per-genome null ⇒ no principled
  species call; an absolute strength threshold on a hand-picked set is meaningless). **No new `-q` fragility**
  from the strength gate: `p_gumbel_p95` is NaN exactly when `z_excess` is NaN (same EVT block), so both
  strength paths are behind the `z_excess`-finite guard and short-circuit to `UNASSESSABLE`.
- **Opt-in** `--adjudicator-min-u2 <N>` (config `scoring.adjudicator_min_u2`) lets a known-representative
  small genome self-adjudicate against its own (noisier) U2 tail.

A well-powered genome (≥ `min_u2` U2) with **zero** strong calls short-circuits to `NOT_DETECTED` (the
clearest loss), before any tail fit.

## 7. Frozen calibration + provenance

Supervised, frozen anchors on `z_excess` (calibration-time snRNA/IPA panel: the 43-genome cross-clade panel +
68 freshly-run divergent bearers [early fungi / oomycetes / basal metazoans] + snRNA-confirmed losses):

- **`loss_ceiling_z = 2.60`** — *aspergillus coremiiformis*, the highest-`z` snRNA-confirmed loss (keeps only
  the conserved u4atac/u6atac; lost u11+u12). All loss anchors confirmed loss at the Infernal `--incE ≤ 0.01`
  inclusion rule, incl. textbook 0/4 ciliates / apicomplexan / green algae.
- **`bearer_floor_z = 5.50`** — a **TRUST threshold, not a separator** (widened from the original 4.00). The
  near-floor zone is genuinely mixed: a probable-loss (Monocercomonoides `z ≈ 4.24`, snRNA-confirmed loss)
  sits *adjacent* to a probable-bearer (Blyttiomyces `z ≈ 4.54`), so `z_excess` alone cannot separate there.
  Above 5.50, a `DETECTED`-by-count call is trustworthy without corroboration; between the ceiling and the
  floor is the `INCONCLUSIVE` gap (corroborate); below is loss unless strength-rescued.
- Platt `(a, c) = (2.796, −1.178)`; strength gate `p_gumbel_threshold = 0.01`, `cs_point_threshold = 5.0`,
  `cs_min_calls = 3`, `cs_p95_pct = 95`; `min_u2 = 200`, `min_evt_excesses = 20`.

**snRNA bearer/loss rule (calibration-time only, never at runtime):** `significant` = best cmsearch
E ≤ 0.01; **BEARER** = ≥3 of {u11,u12,u4atac,u6atac} significant (self-consistently defining-aware — a loss
caps at 2/4, so reaching 3 requires a defining snRNA, *without* demanding the u11+u12 pair, which would wrongly
reject gut fungi with divergent/undetectable u11). See the postmortem §4a. The clade-conditional reliability of
this rule (which snRNAs are trustworthy per clade) is separately calibrated in
`/mnt/data/u12/species_index/SNRNA_CLADE_SENSITIVITY.md` (U12/U6atac robust; U11/U4atac fragile).

Leave-clade-out validation: **0 cross-errors** — held-out genomes never mis-cross the gap; the worst case is an
honest `INCONCLUSIVE`. The binding open work is *coverage* of the divergent-bearer tips that pin `bearer_floor`
(early fungi, oomycetes, basal metazoans), not the method.

## 8. Design evolution (so the failure modes aren't reintroduced)

| stage | driver | outcome | why replaced / kept |
|---|---|---|---|
| v1 | per-species z-normalization | cross-species comparability | z-inflated loss FPs; broke on Amborella/Oryza. Supplanted (`raw_gated_scoring.md`). |
| 2026-06-27 | `depth_tail → q → P_adj` | settled design, `depth_tail_firth_2026-06-27` | reported **0 HC U12 for genuine band-dweller bearers**; count-blind; baked a 2nd threshold off a separated-fit slope. Postmortem. |
| 2026-06-30 | `z_excess` gap gate, two numbers | `zexcess_gap_2026-06-30`; anchors 2.60 / **4.00**; cats DETECTED/BORDERLINE/NOT_DETECTED/UNASSESSABLE | population count done right; abstain in the gap. |
| 2026-06-30c | widen floor 4.00 → **5.50**; rename **BORDERLINE → INCONCLUSIVE** | trust threshold (near-floor zone is mixed, not separable) | commits `9cb9f5c`, `bc8d9de`. |
| 2026-07-03 | + `cs_p95 ≥ 5.0` strength gate | recover few-but-strong divergent bearers | commits `8b8571a`, `5a63ff3`. |
| 2026-07-03b | + `p_gumbel_p95 ≤ 0.01` PRIMARY (cs_p95 co-fallback) | `zexcess_gap_pgumbel_cs_2026-07-03` (**current**) | dissolves the tuned-constant brittleness; per-genome, adaptive. `adjudicator_strength_evt.md`. |

**The lessons that pin this design (postmortem §5):** grade the deliverable not a proxy; a classifier on
separated binary labels anchors exactly one threshold (don't read others off the slope); "does this have a
population" is a **count** question; don't fold a non-calibratable n=1 species signal into a calibrated
per-intron number — report the two separately; keep runtime scope tight (calibration may use snRNA/taxonomy,
the runtime classifier depends only on its own input).

## 9. Code / bundle / tests

- **Code:** `scoring/species_adjudicator.py` — `MotifCategory`, `AdjudicatorParams`, `classify_motif_category`,
  `compute_z_excess` / `compute_p_gumbel` / `_u2_reference`, `_assess`, `apply_pmotif_adjudication`
  (file-side; the shared `_run_post_classification_pipeline` in `cli/main.py` dispatches it so streaming ==
  in-memory by construction).
- **Bundle:** `data/default_pretrained.model.pkl` carries `adjudicator_params` (version-pinned). Re-stamped —
  not retrained — for each gate change; calibrated Platt/anchors preserved. Backups: `*.bak_*_2026-*`.
- **Tests:** `tests/unit/test_scoring/test_species_adjudicator.py`; parity `tests/integration/
  test_pmotif_adjudicated_parity.py` + `test_streaming_equivalence.py`; chr19 end-to-end smoke.

## 10. Related docs

- `raw_gated_scoring.md` — the raw-feature architecture + z-stack supplant log (adjudicator §0a–§0e now historical → this doc).
- `adjudicator_qdriver_postmortem.md`, `adjudicator_strength_evt.md` — the two rationale docs (see §top).
- `adjudicator_call_strength_plan.md` — the `cs_p95` gate plan (superseded by `adjudicator_strength_evt.md`).
- `training_set_construction_protocol.md` — difficulty-stratified phylo-allocated training-set protocol (the sanctioned lever for the coverage gaps that pin `bearer_floor`).
- `/mnt/data/u12/species_index/SNRNA_CLADE_SENSITIVITY.md` — clade-conditional snRNA reliability for the corroboration layer.
