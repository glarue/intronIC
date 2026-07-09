# Unified training/holdout index — `training_data_index.tsv` + `holdout_index.tsv`

> **In this repo the two TSVs are committed gzipped** (`training_data_index.tsv.gz`,
> `holdout_index.tsv.gz`) — the raw `holdout_index.tsv` is 123 MB and exceeds GitHub's
> 100 MB/file limit. `gunzip -k *.gz` to read them. They are a **versioned artifact** of
> the deployed v3 model's data; the canonical build (below) runs against the training
> corpus at `/mnt/data/u12/ipa/training_data/multispecies_v23/` (external to this repo —
> the multi-GB corpus, per-species `.iic` runs, and IPA conservation matrix are not
> vendored here), where `build_training_index.py` writes the raw TSVs.

**What it is.** One row per intron (identical 44-column schema in both files), carrying
(a) identity + provenance, (b) the exact raw features and label, (c) intron sequence +
flanks, and (d) per-intron IPA conservation:

- **`training_data_index.tsv`** — the **deployed model's training rows** (`is_holdout==""`),
  41,261 rows. Retrain-reproducible (see below).
- **`holdout_index.tsv`** — the **held-out evaluation rows** (`is_holdout!=""`),
  202,642 rows across three tags: `recall` (AraTha/DroMel/GalGal/NemVec/RhiIrr — the
  5-species / **790-U12-positive** recall eval, F1=1.0), `protist` (ToxGon loss control),
  `consensus_fp` (18-species U2 false-positive control). Lets the holdout evaluation be
  reproduced. `used_in_training=0` on every holdout row.

It serves two purposes:

1. **Documentation / audit** of the data and where every example comes from — including
   the full intron identity (`parent_gene`, `parent_transcript`, `intron_index` /
   `n_introns_in_tx`, parsed from the intron name).
2. **Retrain-reproducibility** — the training file's `five_raw…support2_raw` + `y` are the
   *exact* matrix the deployed ensemble was fit on, so retraining from it reproduces the
   model state.

Built by `build_training_index.py` (in this directory; emits both raw files into the
external corpus dir it reads from). Rerun (requires the external corpus present):
`/mnt/data/u12/intronIC_v2_env/bin/python build_training_index.py`, then re-gzip the two
outputs into this directory.

## Reproducibility guarantee (how identity was recovered)

The deployed ensemble (`…/conservation_corpus/eval_corpus/raw_gated_42.model.pkl`,
calibrated + adjudicator-stamped into the bundled `default_pretrained.model.pkl`) was
fit on `eval_corpus/stage1_trainmatrix.npz`, which stores only `Xraw/Xz/y/clade` — **no
per-row identity**. The builder (`eval_corpus/stage1_real_classifier_test.py`, Part A)
assembles that matrix from `multispecies_v23/.reorg_backup/multispecies_corpus.tsv`
(non-holdout rows) + a raw-feature join from each species' `score_info.iic`, in corpus
order, discarding names at cache time.

`build_training_index.py` **replays that assembly verbatim while keeping the
`intron_name`/`species` aligned to each row**, then **asserts** the reconstructed
`Xz/Xraw/y/clade` equal the cached npz element-for-element. That assert PASSES
(41,261 rows), which *proves* the emitted names are the exact per-row identity of the
deployed training matrix.

### To reproduce the model from this index
```python
import numpy as np, csv, gzip
from sklearn.svm import SVC
from sklearn.calibration import CalibratedClassifierCV
rows = [r for r in csv.DictReader(gzip.open("training_data_index.tsv.gz", "rt"), delimiter="\t")
        if r["used_in_training"] == "1"]                     # post-NaN-filter fit set
Xraw = np.array([[float(r[c]) for c in
        ("five_raw","bp_raw","three_raw","bp_offset","bp_scan_confidence","support2_raw")]
        for r in rows])
y = np.array([int(r["y"]) for r in rows])
# The deployed ensemble is a 42-model bagged RBF-SVM (see eval_corpus/raw_gated_v2.py);
# a single base fit that matches the fidelity gate is:
clf = CalibratedClassifierCV(SVC(C=200, gamma=0.001, kernel="rbf",
        class_weight="balanced"), method="isotonic", cv=3).fit(Xraw, y)
```
For the full 42-model bagging + Platt/adjudicator stamping, see
`eval_corpus/raw_gated_v2.py` → `calibrate_pmotif_bundle.py`.

## Row set
- **41,261 rows** = the exact assembled deployed matrix (u12 = 10,045, u2 = 31,216).
- `used_in_training == 1` → 41,257 rows (4 dropped by the builder's NaN filter; that
  subset is the actual `.fit()` input).
- Note: this differs slightly from the "41,333 / 10,003 U12" figure in
  `docs/v27_dev_guide.md` — that counts the `multispecies_v23` z-corpus; the deployed
  raw model was fit on the marginally different `eval_corpus/stage1` matrix. This index
  is faithful to **what the model actually saw** (assert-verified).

## Columns
44 columns, identical in both files.

| group | columns |
|---|---|
| identity/provenance | `intron_name` (canonical), `species` (abbrev), `species_full`, `clade`, `phylum`, `genus`, `parent_gene`, `parent_transcript`, `intron_index`, `n_introns_in_tx` (the last four parsed from the intron name: `{sp}-{gene}@{transcript}_{index}({total})`), `chr`, `start_0based`, `end`, `strand`, `position_id`, `source_iic` (originating `.introns.iic` run) |
| label + set membership | `label` (u12/u2), `y` (1/0), `subtype` (provenance class: `positive` / `nc_positive` / `anchor_tp` / `easy_consensus` / `hard_neg_consensus` / …), `gold_tier`, `is_holdout` (`""` = training; `recall` / `protist` / `consensus_fp` in the holdout file), `used_in_training` (1 = in the deployed `.fit()`; 0 for holdout + the 4 NaN-dropped rows) |
| features (reproducibility) | **`five_raw`, `bp_raw`, `three_raw`, `bp_offset`, `bp_scan_confidence`, `support2_raw`** (= the deployed `Xraw`), plus reference `five_z/bp_z/three_z`, `svm_score`, `n_species_in_group`, `n_phyla_in_group` |
| sequence | `intron_len`, `was_truncated`, `up_flank` (50 bp), `intron_seq_trunc`, `down_flank` (50 bp) |
| conservation | `conserved_as_u12` (`1` = passes rule, `0` = assessed-but-fails, `""` = no IPA conservation record), `n_species_conserved_u12`, `n_indep_genera_conserved_u12` (K), `conservation_self_adaptive` (focal score), `conserved_species` (comma-sep partner **binomials**) |

## Sequence scheme
`intron_seq_trunc` is passed through from `corpus_sequences.v6.1.tsv` unchanged. Introns
≤ ~200 bp are full; longer introns are stored as **`first100 + "..." + last100`
(= 203 bp)** — the corpus's own truncation, which already uses the `...` indicator.
`up_flank`/`down_flank` are 50 bp. (If a 250 + `...` + 250 = 500 bp representation is
wanted later, re-extract the intron bodies from the per-species genome FASTAs by
`chr/start_0based/end/strand` — flanks can be reused as-is.)

## Conservation rule (locked; from `audit_conservation_matrix.py::truth_count`)
Source: `intron_conservation_long.tsv` (4.3 M focal×partner rows), grouped by focal
`intron_key`. A partner counts if `partner_adaptive > 90`.
- `n_species_conserved_u12` = distinct such partner species.
- `n_indep_genera_conserved_u12` = **K** = distinct such partner *genera*, excluding the
  focal's own genus (guards against inflation by close relatives).
- `conserved_as_u12` = `(self_adaptive≥90 & K≥1) | (≥80 & K≥2) | (≥75 & K≥3)`.
- `conserved_species` lists the qualifying partner species as **full binomials**
  (underscore form, e.g. `gallus_gallus,homo_sapiens,latimeria_chalumnae`), sorted.
  Fixed-panel abbreviations from the source matrix (`GalGal`, `XenTro`, …) are mapped to
  binomials via `eval_corpus/abbrev_score_info.tsv` (`to_binom` in the builder). Full
  per-partner detail (scores, alignment ids) is in `intron_conservation_long.tsv`.

Coverage: 3,346 of 10,045 U12 positives have ≥1 qualifying partner; 3,310 pass the
graded rule. U2 negatives carry no conservation (blank / 0). The conservation *scores*
are the v2.6-era `adaptive` (post-mode-sep) values that the IPA matrix was built with —
they are the conservation evidence, independent of the v3 raw-feature classifier.

## Note on `label_source`
Deliberately not included: in `gold_corpus.tsv` it is a near-constant `conservation`
tag (records that gold labels came from the IPA-conservation pipeline), so it carries no
per-intron information. The informative label provenance is the `subtype` column
(`positive` / `nc_positive` / `anchor_tp` / `easy_consensus` / `hard_neg_consensus` / …).
