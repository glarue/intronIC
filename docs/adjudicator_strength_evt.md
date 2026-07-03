# Adjudicator strength gate — per-genome Gumbel outlier test (`p_gumbel_p95`)

**Status: IMPLEMENTED 2026-07-03b** (`ADJUDICATOR_PARAMS_VERSION = zexcess_gap_pgumbel_cs_2026-07-03`).
The species adjudicator's call-strength gate is now driven by a **per-genome anomaly test on the call upper
tail** — the Gumbel tail-probability that a genome's own U2 background produces a maximum as strong as the
call `cs_p95` — with the previous absolute `cs_p95 ≥ 5.0` retained as a null-free **co-fallback**.

Supersedes the tuned-constant framing in `docs/adjudicator_call_strength_plan.md` (still accurate for the
`cs_p95` statistic and stress tests; this doc replaces its §2 gate rule). Finding + full evidence:
`eval_corpus/{evt_strength_test,evt_sweep,analyze_evt}.py`.

## 1. Why (the problem with a tuned constant)

`z_excess` is a **count** statistic → blind to divergent bearers that have a *few* genuinely strong U12
calls (low count → low z). The call-strength gate `cs_p95 ≥ 5.0` recovers them, but 5.0 was a **global
constant tuned to the eval panel** — it only works because the deployed genomes happen to share a similar
U2 margin scale, and would silently mis-calibrate on a genome with an unusual U2 background. That fragility
was the objection this change resolves.

## 2. The statistic

The adjudicator already fits, per genome, a Gumbel EVT null of the U2 margin tail (`z_expmax`, `evt_beta`
in `_u2_reference`) — the same null that drives `z_excess`. The shipped `excess_z`/`p_gumbel` evaluated that
null at the **call median** (`call_core`), which is central-tendency blind: a divergent bearer's median call
is *buried in the U2 background*. The fix is one line — evaluate the **same null at the call upper tail**:

```
p_gumbel_p95 = compute_p_gumbel(cs_p95, ref)      # Gumbel tail-prob at the p95 of the call margins
```

`p_gumbel_p95` = P(the genome's own U2 background produces a max ≥ the call `cs_p95`) — a **one-sided
outlier p-value referenced to each genome's own U2 null**. It is composition-adaptive: the same evidence is
judged against each genome's own `med`/`mad`/`z_expmax`, not a global margin constant.

## 3. The gate (`classify_motif_category`)

```
if not finite(z_excess):                 -> UNASSESSABLE        # cannot reference the genome's own tail
if z_excess >= bearer_floor_z (5.5):     -> DETECTED            # count gate (unchanged)
if n_calls >= cs_min_calls (3) and (
     p_gumbel_p95 <= p_gumbel_threshold (0.01)    # PRIMARY: per-genome Gumbel outlier test
     or cs_p95 >= cs_point_threshold (5.0) ):     # CO-FALLBACK: null-free absolute
                                         -> DETECTED
if z_excess <= loss_ceiling_z (2.6):     -> NOT_DETECTED
else:                                    -> INCONCLUSIVE
```

- **Primary = `p_gumbel_p95 ≤ 0.01`** — the per-genome ~1% outlier test.
- **Co-fallback = `cs_p95 ≥ 5.0`** (OR) — null-free, stays stable at marginal EVT fits where the per-genome
  null is noisy. The two agree on the panel, so the OR recovers the same set at 0 FP; the fallback only adds
  robustness, never new false positives (worst loss `p_gumbel_p95 = 0.044`, well above 0.01).
- Both strength paths are **behind the `z_excess`-finite guard**, so a genome that cannot fit its own U2
  tail stays `UNASSESSABLE` — the gate never fabricates a call without a per-genome null (see §5).

## 4. Evidence (46 z-missed bearers / 26 losses with ≥3 calls, deployed model)

**Upper tail is the unlock** — same Gumbel null at each call percentile:

| evaluate at | recovery @ 0-FP | AUC |
|---|---|---|
| median (`excess_z`, shipped secondary) | 12/46 | 0.72 |
| p80 | 23/46 | 0.81 |
| p90 | 31/46 | 0.87 |
| **p95** | **36/46** | 0.89 |
| p99 / max | 36/46 | 0.92 / 0.92 |

Recovery **triples** median→p95 and **saturates at p95** (p99/max add nothing but `max` is fragile to a lone
relic) — so p95 is the principled knee. Per-genome: the Trichinella `excess_z` is *negative* through the
median (reads as loss) and only crosses into significance in the top ~10%.

**Equivalent to the tuned gate, but derived.** Concordance of recoveries: the `p_gumbel_p95 ≤ 0.01` and
`cs_p95 ≥ 5.0` gates pick the **same genomes** (differ on one borderline Blastocystis at the very ceiling).
Correspondence: `excess_z_p95 = 0.87·cs_p95 − 2.21` (r = 0.94); **`cs_p95 = 5.0 ↔ p_gumbel_p95 ≈ 0.005`**
(genomes near 5.0 span p = 2.6e-3 … 1.2e-2). So the tuned 5.0 was a per-genome ~0.5–1% outlier test all along
— now stated as one, and adaptive where the constant would drift.

## 5. `-q` / low-N fallback (important)

The per-genome test presumes a genome's U2 background exists. In `-q` (pre-extracted sequence) mode:

- **Full-complement `-q`** (thousands of U2) → EVT tail fits → normal adjudication. (135/135 real genomes fit.)
- **Curated / tiny `-q`** (`< min_u2 = 200` U2) → `LOW_N` → `UNASSESSABLE` species call, per-intron `P_motif`
  reported as-is, no bearer/loss call — the honest answer (no per-genome null ⇒ no principled *or* meaningful
  call; an absolute strength threshold on a hand-picked set is meaningless anyway).
- Opt-in `adjudicator_min_u2` (CLI/config) lets a known-representative small set self-adjudicate against its
  own noisier tail.

Adoption introduces **no new `-q` fragility**: `p_gumbel_p95` is NaN exactly when `z_excess` is NaN (same EVT
block), and the shipped gate already short-circuits to `UNASSESSABLE` on non-finite `z_excess`.

## 6. Safety / no-regression

The gate is an **OR that retains `cs_p95 ≥ 5.0`**, so every genome DETECTED under the prior gate is still
DETECTED — the change can only *add* the `p_gumbel` path, and that path adds DETECTED only for `p ≤ 0.01`
while the worst loss reaches `p = 0.044` → no new loss FP. Verified: chr19 end-to-end (column present,
`p_gumbel_p95 = 3.2e-4`, DETECTED), adjudicator unit suite, streaming↔in-memory parity, bundle re-stamp
(all calibrated Platt/q/anchor/`cs_point` values preserved).

## 7. Code / bundle

- `src/intronIC/scoring/species_adjudicator.py`: `p_gumbel_p95 = compute_p_gumbel(cs_p95, ref)` in `_assess`;
  `classify_motif_category(z_excess, cs_p95, p_gumbel_p95, n_calls, params)`; `AdjudicatorParams
  .p_gumbel_threshold = 0.01`; `AdjudicatorResult.p_gumbel_p95`; written to `score_info.iic` + `metrics`.
- `src/intronIC/cli/main.py`: `p_gumbel_p95` plumbed through `PostClassResult` + the metrics assembler
  (both classify paths).
- Bundle: `default_pretrained.model.pkl` re-stamped (adds `p_gumbel_threshold`, bumps `params_version`;
  gate was already live at 0.01 via `from_dict` fallback). Backup: `*.bak_pre_pgumbel_2026-07-03`.
