# Call-strength integration into the species adjudicator — IMPLEMENTED (single point gate)

> **SUPERSEDED as the PRIMARY driver (2026-07-03b).** `cs_p95 >= 5.0` is now the null-free **co-fallback**;
> the primary strength driver is the per-genome Gumbel test `p_gumbel_p95 <= 0.01` (a derived, adaptive
> restatement of the same 5.0 outlier bar). See [`adjudicator_strength_evt.md`](adjudicator_strength_evt.md)
> and the consolidated [`adjudicator_design.md`](adjudicator_design.md). This doc remains accurate for the
> `cs_p95` statistic itself and its stress tests.

**Status: IMPLEMENTED 2026-07-03.** The adjudicator now gates on the POINT `cs_p95 >= 5.0` (single track).
The earlier three-track OR design (this doc's prior draft) was **rejected by the stress tests** — a single
point-p95 gate dominates. Evidence of record: `eval_corpus/STRENGTH_STRESS_FINDINGS.md`. Finding + validation:
`docs/call_strength_recovery.md`; memory `call-strength-divergent-recovery`.

## 1. What it adds and why
`z_excess` (a COUNT statistic) recovers 0 of the count-blind divergent bearers — a genome with a *few*
genuinely strong U12 calls (low count → low z) is missed. The un-clipped ensemble MARGIN retains call
STRENGTH that both `z_excess` and the Platt-saturated `P_motif` discard. `cs_p95` = the 95th percentile of the
call margins (never max): a loss's calls top out ~bounded, a bearer's real U12s reach much higher; p95
discounts a lone relic. Gating on `cs_p95 >= 5.0` recovers the divergent bearers (Trichinella/Trichuris,
oomycetes, Mucor, mites) at 0 loss FP.

## 2. Design — a SINGLE point-p95 gate (three-track rejected)
`motif_category = DETECTED` if `z_excess >= bearer_floor_z` (existing) **OR** (`cs_gate_enabled` and
`cs_p95 >= cs_point_threshold` and `n_calls >= cs_min_calls`). No lower-bound (`cs_p95_lo`) track.

**Why not the CI lower bound (the crux).** A `cs_p95_lo` gate self-calibrates to call count — but that is
exactly wrong: divergent bearers ARE few-call-but-strong, so the CI penalizes the signal we recover.
Trichinella nativa has 5 strong calls (point 6.24) yet `lo@10` = 2.24 == the relic-loss Micromonas (2.08):
the lower bound cannot separate "few but strong" (bearer) from "few and lucky" (relic loss) and suppresses the
bearers (lo@10 recovers 8/14 vs point's 10/14 at the same 0 FP; 4 of the 10 strongest bearers — all Trichinella
+ Pythium — fire point but fail lo). Relic-loss robustness instead comes from three things that DON'T penalize
strong-but-sparse signal: **p95** (ignores a lone above-95th relic), **`cs_min_calls >= 3`** (kills degenerate
k<3 where p95==max, e.g. Micromonas), and the **fat threshold margin** (5.0 sits ~1.36 above the k≥3 loss
ceiling ~3.64). `cs_p95_lo`/`cs_p95_hi` are still computed + logged as a per-genome confidence annotation.

## 3. Calibration (final panel of record)
58 bearers / 38 losses / 25 clades — index-labeled (losses by PRIOR incl composition-FP losses), snRNA-
adjudicated (Oxytricha→bearer; Blastocystis×3→bearer/uncertain-gain; Hypsibius/Heterostelium/Emiliania→excluded
UNCERTAIN; Paramacrobiotus snRNA searched + written back). **point `cs_p95 >= 5.0`: 10/14 z-missed divergent
bearers @ 0 loss FP, TPR 91%, FLAT across [4.25, 5.25].** 5.5 halves recovery (5/14); 4.25 is the low end of
the flat band. 5.0 chosen: near the top of the empty [4.11→5.37] separation band → max injection margin while
keeping all 10. `cs_min_calls=3`, `cs_p95_pct=95`. Injection stress: a single spurious strong call must exceed
~5.5 (and only lifts low-k losses); two coordinated strong calls defeat any peak gate (inherent → downstream
snRNA/IPA corroboration on DETECTED). Bonus: `cs_p95 >= 5.0` also cleanly adjudicates CONFLICT_motif_blind
(validated divergent bearers all ≥5.37, uncertain ones all ≤4.11).

## 4. Code state (DONE)
`src/intronIC/scoring/species_adjudicator.py`:
- `classify_motif_category(z_excess, cs_p95, n_calls, params)` gates on the POINT `cs_p95 >= cs_point_threshold`.
- `AdjudicatorParams.cs_point_threshold = 5.0` (`DEFAULT_CS_POINT_THRESHOLD`); `cs_min_calls=3`, `cs_p95_pct=95`,
  `cs_gate_enabled=True`. `cs_p95_threshold=4.50` kept as a LEGACY (unused) field so old-stamped bundles load.
- `cs_p95_lo`/`cs_p95_hi` still computed in the shared bootstrap (q-CI parity preserved) + written to
  score_info + metrics as confidence annotations.
- `ADJUDICATOR_PARAMS_VERSION = "zexcess_gap_cs_point_2026-07-03"` (z_excess stays PRIMARY; cs_point additive).
- Tests updated (`tests/unit/test_scoring/test_species_adjudicator.py`): gate test uses `cs_point_threshold`.

**Verified:** 33/33 adjudicator unit tests; 621 full unit suite; streaming↔in-memory parity (canonical gate);
pmotif-adjudicated parity; chr19 end-to-end (cs_p95 columns present, DETECTED via z-gate).

## 5. Remaining (release step)
**Bundle re-stamp.** The deployed `default_pretrained.model.pkl` stamps `adjudicator_params` from the pre-
call-strength era (no `cs_*` keys) → they fall back to code defaults, so the gate is ALREADY live at 5.0 and
correct functionally. For clean provenance, re-stamp with `cs_point_threshold=5.0`, `cs_min_calls=3`,
`cs_p95_pct=95.0`, `cs_ci_lo_pct=10.0`, `cs_gate_enabled=True`, `params_version=zexcess_gap_cs_point_2026-07-03`
per `docs/v3_bundle_schema.md`'s bundle-change release pattern.

## 6. Honest residuals
- Two coordinated injected strong calls defeat any peak gate — mitigate via downstream corroboration on DETECTED.
- Genuine few-strong-among-many bearers (Blastocystis-atcc: 2 calls @ ~6.95 of 93) are missed by point p95 — a
  count-of-strong / top-k statistic would catch them (future option; using p95 per GL). Oxytricha (real ciliate
  bearer, point 4.04) is likewise missed — inseparable from the loss ceiling on this axis; caught downstream.
