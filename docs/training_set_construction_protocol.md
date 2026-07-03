# Training-set construction protocol (difficulty-stratified, phylo-allocated)

**Status: DRAFT proposed 2026-07-03.** A formal, reusable protocol for building the U12/U2
classifier's training corpus so the decision boundary is set by the cases we care about
(U12/U2-confusable motifs, composition-FP-prone clades) with phylogenetic representation controlled
*independently* of example difficulty — so no over-sequenced clade (vertebrates, grasses) drowns the
signal, and no dilution softens the boundary.

This protocol is grounded in two experiments that **closed** the obvious corpus levers, and it exists
so future corpus work is a disciplined, falsifiable procedure rather than intuition:

- **Corpus augmentation is closed** — injecting divergent U12 positives (Trichinella + oomycete +
  Diptera) regressed both axes (−341 clean-bearer recall / +91 loss FP). Divergent U12 and residual
  loss-FP overlap in feature space; adding divergent positives admits more loss FPs. See
  `docs/V6_CORPUS_AUGMENT_PLAN.md` (in the training-data dir) and memory `training-corpus-review`.
- **Phylo-balanced negatives are closed** — swapping the deployed stage1 per-clade-capped hard
  negatives for a phylo-balanced 100k pool (positives held fixed) held the loss ceiling (4.04→4.10,
  0/85 FP) but **regressed divergent-bearer recovery 36→29/50**, concentrated on the nematode
  Trichinella (the out-of-clade cases with 0 training positives). Mechanism: ~2× call inflation =
  looser boundary from diluted hard-negative density. See memory `phylobal-negatives-abtest`;
  scripts `eval_corpus/{train_phylobal_ab,eval_focused,analyze_focused}.py`.

---

## 1. First principles

1. **The boundary is made of boundary cases.** For the RBF-SVM ensemble, only support vectors —
   boundary-adjacent examples — define the decision surface; interior mass is inert to it. Corpus
   construction should therefore **maximize the density and quality of boundary cases** and treat
   interior examples as calibration ballast, not signal. (Empirically: phylobal's 100k mostly-easy
   negatives gave a *looser* boundary than stage1's 31k concentrated ones.)

2. **Representation ⟂ difficulty.** *Which clades* are represented and *which introns within a clade*
   are hard are orthogonal knobs. **Balance the allocation across clades; hard-mine within each
   clade.** Balancing away difficulty — filling a balanced per-clade budget with random/easy examples
   — is the only thing that genuinely "bastardizes balanced," and it is exactly the failure mode the
   phylobal experiment exhibited.

3. **Error costs are asymmetric → allocation is risk-weighted, not flat.** Composition FPs concentrate
   in loss-prior clades (Dikarya, kinetoplastids, microsporidia, red algae, AT-rich nematode
   Chromadorea). Those clades receive a *larger* hard-negative budget. "Imbalance toward the hard loss
   clades" is a **feature** — it is what holds the loss ceiling and preserves the strong-U12 margin the
   call-strength gate reads.

---

## 2. Definitions

- **phylo-unit** — a clade/phylum at the corpus's working rank (the ~14-clade axis). The unit of
  allocation budgeting.
- **class** — U12 (positive) or U2 (negative).
- **difficulty tier** — `boundary` or `interior`, defined per class by **model-independent** criteria
  (§3 step 2). Never defined by the incumbent model's own scores (anti-circularity, §5).
- **cell** — a `(phylo-unit × class × tier)` bucket. Construction allocates and fills cells.

---

## 3. Construction procedure

1. **Partition** all candidate introns into cells `(phylo-unit × class × tier)`.

2. **Define "boundary" by model-independent criteria** (not the incumbent SVM margin):
   - *Boundary negatives* = U2 introns matching U12-conserved positions (U12 donor/branch consensus —
     the `easy_consensus` construction), **plus** U2 from U12-absent / loss-prior species
     (`easy_biology`). Motif- and biology-defined.
   - *Boundary positives* = conservation-anchored U12 (IPA orthology evidence: ≥3 species / ≥2 phyla).
   - Everything else is `interior`.

3. **Deduplicate by orthology/isoform BEFORE allocation.** Collapse a U12 intron conserved across N
   species to **one representative per ortholog group** (reciprologs/orthogroups from the IPA stack),
   and isoform-collapse within a species. This removes vertebrate pseudoreplication *at the source*
   (positives are otherwise K-fold-replicated vertebrate orthologs, e.g. HomSap 536) rather than
   papering over it with post-hoc reweighting.

4. **Set per-clade budgets** `B_c = B_base × risk_c`:
   - `B_base` — the balanced floor giving each clade roughly equal voice (representation).
   - `risk_c ≥ 1` — a multiplier boosting the *negative* budget of loss-prior clades (the ceiling).

5. **Fill budgets boundary-first.** Draw from the `boundary` tier until the budget or the tier is
   exhausted; only then add `interior` examples, capped as ballast. This maximizes hard-case density
   per clade *under* a representation-balanced allocation — the two principles composed.

6. **Positives: represent + dedupe, do NOT divergent-augment.** Divergent-U12 recall is explicitly
   **out of scope** for the corpus (§ closed lever). The corpus's job is a clean, tight, representative
   boundary; divergent recovery is the output gate's job (`species_adjudicator` call-strength gate).

---

## 4. Mandatory validation gate

Every candidate corpus is validated on a **held-out** panel never used in training —
snRNA-corroborated bearers (including divergent) + loss-prior-clade genomes (including composition-FP
ones). This is exactly the A/B harness in `eval_corpus/{train_*,eval_focused,analyze_focused}.py`.

- **Two axes, both measured at the FIXED deployed operating point** (operating-point-invariant — you
  may **not** move the gate to manufacture recovery):
  1. **Divergent-bearer recovery at 0 loss-FP** — # z-missed bearers whose `cs_p95` clears the loss
     ceiling.
  2. **Loss ceiling** — max loss `cs_p95` (k ≥ 3 calls).
- **No-regression rule** — the candidate must not regress *either* axis vs the incumbent. Lifting
  recovery while raising the ceiling (or vice versa) **fails**.
- **Scale note** — cross-ensemble comparison uses each ensemble's own Platt call threshold
  (margin@P=0.9) so `cs_p95` is comparable; verify the thresholds coincide (they did here: dep/ctl/trt
  ≈ 1.21) or normalize before comparing.
- Any **accepted** change triggers a full Platt/q/cs_p95 re-calibration + re-validation.

---

## 5. Guardrails

- **Anti-circularity** — boundary membership comes from stable motif/phylo/conservation criteria, never
  the current model's scores. Mining hard cases with the incumbent margin entrenches its blind spots
  and cannot discover the cases it already gets wrong (see memory `divergent-bearer-recall`).
- **Closed lever** — no divergent-positive augmentation to chase recall. That path is closed by
  experiment; recall is an output-gate objective.
- **Re-calibration is not optional** — any composition change invalidates the Platt/q/cs_p95
  calibration.

---

## 6. How the deployed corpus maps onto this protocol

The deployed `stage1` corpus already instantiates most of the protocol, which is why it wins:

| Protocol element | stage1 (deployed) | gap / future |
|---|---|---|
| boundary negatives, model-independent | ✓ `easy_consensus` + `easy_biology` | — |
| risk-weighted loss-clade negatives | ✓ Nematode 7148 / Dikarya 3846 | — |
| boundary-first fill | ✓ (both neg components are hard-mined) | — |
| phylo-balanced allocation | partial (per-clade caps 2k/5k) | equalize `B_base` across clades |
| positives deduped by orthology | ✗ vertebrate pseudoreplication remains | **the main open item** |
| no divergent augmentation | ✓ | — |

The one clean, low-risk future improvement the protocol surfaces is **orthology-dedupe of the
positives** (step 3) to remove vertebrate over-weighting — additive to, not a replacement of, the
hard-negative density that already works. Everything else the deployed corpus already does well; the
protocol mainly makes *why* explicit and gives future changes a falsifiable gate.

---

## References

- `docs/raw_gated_scoring.md` — the raw pmotif_adjudicated architecture the corpus feeds.
- `docs/adjudicator_call_strength_plan.md` — the output-level call-strength gate (where divergent
  recall lives).
- memory: `phylobal-negatives-abtest`, `training-corpus-review`, `call-strength-divergent-recovery`,
  `divergent-bearer-recall`.
- `eval_corpus/{train_phylobal_ab,eval_focused,analyze_focused}.py` — the validation-gate harness.
