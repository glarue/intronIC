#!/usr/bin/env python3
"""
Unified index of the intronIC v3 (pmotif_adjudicated) TRAINING + HOLDOUT data.

Emits two same-schema tables in this directory:
  - training_data_index.tsv : the deployed model's training rows (is_holdout == "")
  - holdout_index.tsv        : the held-out evaluation rows (is_holdout != "")

Each row = one intron with (a) identity + provenance (incl. parsed parent gene /
transcript / ordinal), (b) the exact raw features + label, (c) sequence + flanks,
(d) per-intron IPA conservation.

Reproducibility: the training assembly replays Part A of
eval_corpus/stage1_real_classifier_test.py VERBATIM and ASSERTS the reconstructed
Xz/Xraw/y/clade equal the cached stage1_trainmatrix.npz element-for-element, so the
emitted raw features + y are exactly what the deployed ensemble was fit on. The holdout
rows are assembled by the identical path (matching the builder's load_holdout()), so the
5-species/790-positive recall eval (F1=1.0) is reproducible from holdout_index.tsv.

See REBUILD.md for the column dictionary, the conservation rule, and the retrain recipe.
"""
import csv, gzip, os, re, sys
import numpy as np
from collections import defaultdict

sys.path.insert(0, "/mnt/data/u12/ipa/conservation_corpus/scripts")
import iic_io
cn = iic_io.canonical_name
csv.field_size_limit(1 << 24)

TD    = "/mnt/data/u12/ipa/training_data"
M     = f"{TD}/multispecies_v23"
EVAL  = "/mnt/data/u12/ipa/conservation_corpus/eval_corpus"
CORP  = f"{M}/.reorg_backup/multispecies_corpus.tsv"
CACHE = f"{EVAL}/stage1_trainmatrix.npz"
SEQS  = f"{M}/corpus_sequences.v6.1.tsv"
CONS  = f"{M}/intron_conservation_long.tsv"
OUT_TRAIN = f"{M}/training_index/training_data_index.tsv"
OUT_HOLD  = f"{M}/training_index/holdout_index.tsv"

# ---------- builder helpers (copied VERBATIM from stage1_real_classifier_test.py) ----------
def _f(x):
    try:
        v = float(x); return v if v == v else np.nan
    except (TypeError, ValueError):
        return np.nan
def support2(a, b, c):
    return np.sort(np.clip(np.column_stack([a, b, c]), 0, None), axis=1)[:, 1]

ab2dir, ab2clade = {}, {}
for r in csv.DictReader(open(f"{M}/panel_lineage.tsv"), delimiter="\t"):
    ab2dir[r["abbrev"]] = r["dir_name"]; ab2clade[r["abbrev"]] = r["clade"]
GROUP_OVERRIDE = {
    "AchHyp": "Oomycetes", "AphAst": "Oomycetes", "SapPar": "Oomycetes", "ThrCla": "Oomycetes",
    "ChlRei": "GreenAlgae", "CocSub": "GreenAlgae", "HaeLac": "GreenAlgae",
    "VolAfr": "GreenAlgae", "VolCar": "GreenAlgae", "VolRet": "GreenAlgae", "KleNit": "Charophytes",
    "AdiCap": "Basal_plants", "AdiNel": "Basal_plants", "CerRic": "Basal_plants",
    "DipCom": "Basal_plants", "SelMoe": "Basal_plants", "CerPur": "Basal_plants",
    "BasRan": "EarlyFungi", "SynFus": "EarlyFungi", "SynPlu": "EarlyFungi", "SynPse": "EarlyFungi",
    "ParSed": "EarlyFungi", "SpiPun": "EarlyFungi", "TetThe": "Ciliates",
}
CLADE_REMAP = {"Bryophytes_basal_plants": "Basal_plants", "Chytrid_Zoopago": "EarlyFungi"}
def clade_of(ab, base):
    return GROUP_OVERRIDE.get(ab) or CLADE_REMAP.get(base, base)

ab2si_eval, ab2bino = {}, {}
for r in csv.DictReader(open(f"{EVAL}/abbrev_score_info.tsv"), delimiter="\t"):
    ab2si_eval[r["abbrev"]] = r["score_info"]; ab2bino[r["abbrev"]] = r["binomial"]

# conserved_species is emitted as full binomials (underscore form). Partner ids in the
# conservation matrix mix fixed-panel abbrevs (GalGal, XenTro, ...) with binomials;
# map the abbrevs, pass binomials through unchanged.
AB2BINOM = {a: b.strip().replace(" ", "_") for a, b in ab2bino.items()}
_ABBR = re.compile(r"[A-Z][a-z]{2,}[A-Z][a-z]{2,}\d*$")
def to_binom(sp):
    return AB2BINOM.get(sp, sp) if _ABBR.fullmatch(sp) else sp

def si_for(ab):
    d = ab2dir.get(ab)
    if d:
        p = f"{TD}/v23_runs/{d}/{d}.score_info.iic"
        if os.path.exists(p): return p
    return ab2si_eval.get(ab)

RAWCOLS = ["5'_raw", "bp_raw", "3'_raw"]
def raw_for_species(ab, names):
    p = si_for(ab); out = {}
    if not p or not os.path.exists(p): return out
    op = gzip.open(p, "rt") if p.endswith(".gz") else open(p)
    hdr = op.readline().rstrip("\n").split("\t")
    if not all(c in hdr for c in RAWCOLS): op.close(); return out
    ix = [hdr.index(c) for c in RAWCOLS]; ni = hdr.index("name")
    for line in op:
        f = line.rstrip("\n").split("\t"); c = iic_io.canonical_name(f[ni])
        if c in names: out[c] = [_f(f[j]) for j in ix]
    op.close(); return out

# parse the canonical intron name -> parent gene / transcript / ordinal
NAME_RE = re.compile(r"^(?P<sp>[^-]+)-(?P<gene>.+?)@(?P<tx>.+?)_(?P<idx>\d+)\((?P<tot>\d+)\)$")
def parse_name(nm):
    m = NAME_RE.match(nm)
    if not m:
        return "", "", "", ""
    gene, tx = m.group("gene"), m.group("tx")
    for p in ("gene-", "gene:"):
        if gene.startswith(p): gene = gene[len(p):]; break
    for p in ("rna-", "rna:", "mrna-", "transcript-", "transcript:"):
        if tx.startswith(p): tx = tx[len(p):]; break
    return gene, tx, m.group("idx"), m.group("tot")

# ---------- assembly (shared by both sets); mirrors the deployed builder ----------
def assemble(keep):
    """keep(row)->bool selects corpus rows. Returns arrays + aligned names/species/prov."""
    rows = []; by_sp = defaultdict(set); prov = {}
    for r in csv.DictReader(open(CORP), delimiter="\t"):
        if not keep(r):
            continue
        nm = cn(r["intron_name"]); sp = r["species"]
        rows.append((nm, sp, clade_of(sp, ab2clade.get(sp, "?")), 1 if r["label"] == "u12" else 0,
                     _f(r["five_z"]), _f(r["bp_z"]), _f(r["three_z"]), _f(r["bp_offset"]), _f(r["bp_scan_confidence"])))
        by_sp[sp].add(nm); prov[nm] = r
    rawmap = {sp: raw_for_species(sp, names) for sp, names in by_sp.items()}
    Xz, Xraw, y, clade, names, specs = [], [], [], [], [], []
    miss = 0
    for nm, sp, cl, yy, z5, zb, z3, boff, bconf in rows:
        rw = rawmap.get(sp, {}).get(nm)
        if rw is None:
            miss += 1; continue
        Xz.append([z5, zb, z3, boff, bconf]); Xraw.append(rw + [boff, bconf])
        y.append(yy); clade.append(cl); names.append(nm); specs.append(sp)
    Xz = np.array(Xz, float); Xraw = np.array(Xraw, float); y = np.array(y); clade = np.array(clade)
    Xz = np.column_stack([Xz, support2(Xz[:, 0], Xz[:, 1], Xz[:, 2])])
    Xraw = np.column_stack([Xraw, support2(Xraw[:, 0], Xraw[:, 1], Xraw[:, 2])])
    return Xz, Xraw, y, clade, names, specs, prov, miss

print("Assembling TRAINING set (is_holdout == '')...", flush=True)
tXz, tXraw, ty, tclade, tnames, tspecs, tprov, tmiss = assemble(lambda r: r["is_holdout"].strip() == "")
print(f"  {tXraw.shape}  raw-join misses: {tmiss}", flush=True)

# FIDELITY ASSERT vs the cached deployed matrix (reproducibility proof)
d = np.load(CACHE, allow_pickle=True)
assert tXz.shape == d["Xz"].shape and tXraw.shape == d["Xraw"].shape, "shape mismatch vs npz"
assert np.allclose(tXz, d["Xz"], equal_nan=True) and np.allclose(tXraw, d["Xraw"], equal_nan=True), "feature mismatch vs npz"
assert np.array_equal(ty, d["y"]) and np.array_equal(tclade, d["clade"]), "label/clade mismatch vs npz"
print(f"  IDENTITY ASSERT PASSED: {len(tnames)} rows align 1:1 to stage1_trainmatrix.npz "
      f"(u12={int(ty.sum())} u2={int((ty==0).sum())})", flush=True)

print("Assembling HOLDOUT set (is_holdout != '')...", flush=True)
hXz, hXraw, hy, hclade, hnames, hspecs, hprov, hmiss = assemble(lambda r: r["is_holdout"].strip() != "")
print(f"  {hXraw.shape}  raw-join misses: {hmiss}  "
      f"(u12={int(hy.sum())} u2={int((hy==0).sum())})", flush=True)

# ---------- sequences + conservation, streamed ONCE for the union ----------
allnames = set(tnames) | set(hnames)
print(f"Joining sequences for {len(allnames)} introns...", flush=True)
seqmap = {}
for r in csv.DictReader(open(SEQS), delimiter="\t"):
    k = cn(r["intron_name"])
    if k in allnames:
        seqmap[k] = r

print("Aggregating IPA conservation (streaming the 4.3M-row long table)...", flush=True)
cons = {}
for r in csv.DictReader(open(CONS), delimiter="\t"):
    k = cn(r["intron_key"])
    if k not in allnames:
        continue
    e = cons.get(k)
    if e is None:
        e = cons[k] = {"sc": _f(r["self_adaptive"]), "genus": r["genus"], "partners": []}
    e["partners"].append((r["partner_species"], r["partner_genus"], _f(r["partner_adaptive"])))
print(f"  focal introns with conservation data: {len(cons)}", flush=True)

def conserved_fields(key):
    e = cons.get(key)
    if e is None:
        # no IPA conservation record for this intron: conserved_as_u12 left blank
        # ("" = not assessed, distinct from "0" = assessed-but-fails-rule)
        return "", 0, 0, "", "NA"
    sc = e["sc"]; fg = e["genus"]
    strong = [(psp, pgen) for (psp, pgen, pa) in e["partners"] if pa == pa and pa > 90]
    sp_set = sorted({to_binom(psp) for psp, pgen in strong})
    K = len({pgen for psp, pgen in strong if pgen and pgen != fg})
    truth = bool(sc == sc and (
        (sc >= 90 and K >= 1) or (sc >= 80 and K >= 2) or (sc >= 75 and K >= 3)))
    sc_str = "NA" if sc != sc else f"{sc:.2f}"
    return ",".join(sp_set), len(sp_set), K, ("1" if truth else "0"), sc_str

# ---------- emit ----------
RAWNAMES = ["five_raw", "bp_raw", "three_raw", "bp_offset", "bp_scan_confidence", "support2_raw"]
HEADER = (
    ["intron_name", "species", "species_full", "clade", "phylum", "genus",
     "parent_gene", "parent_transcript", "intron_index", "n_introns_in_tx",
     "chr", "start_0based", "end", "strand", "position_id",
     "label", "y", "subtype", "gold_tier", "source_iic", "is_holdout", "used_in_training",
     "n_species_in_group", "n_phyla_in_group"]
    + RAWNAMES
    + ["five_z", "bp_z", "three_z", "svm_score",
       "intron_len", "was_truncated", "up_flank", "intron_seq_trunc", "down_flank",
       "conserved_as_u12", "n_species_conserved_u12", "n_indep_genera_conserved_u12",
       "conservation_self_adaptive", "conserved_species"]
)

def g(dct, *keys):
    for k in keys:
        v = dct.get(k, "")
        if v not in ("", None):
            return v
    return ""

def emit(path, names, specs, clade, Xz, Xraw, y, prov):
    used = ~(np.isnan(Xz).any(1) | np.isnan(Xraw).any(1))
    with open(path, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t"); w.writerow(HEADER)
        for i, nm in enumerate(names):
            pr = prov.get(nm, {}); sq = seqmap.get(nm, {})
            pg, ptx, pidx, ptot = parse_name(nm)
            cf_sp, cf_nsp, cf_K, cf_truth, cf_sc = conserved_fields(nm)
            hold = pr.get("is_holdout", "").strip()
            row = [
                nm, specs[i], g(sq, "species_full"), clade[i],
                g(sq, "phylum"), g(sq, "genus"),
                pg, ptx, pidx, ptot,
                g(pr, "chr"), g(pr, "start_0based"), g(pr, "end"), g(pr, "strand"), g(pr, "position_id"),
                g(pr, "label"), int(y[i]), g(pr, "subtype"), g(sq, "gold_tier"), g(sq, "source_iic"),
                hold, int(used[i]) if not hold else 0,
                g(pr, "n_species_in_group"), g(pr, "n_phyla_in_group"),
            ]
            row += [f"{Xraw[i, j]:.6g}" if Xraw[i, j] == Xraw[i, j] else "NA" for j in range(6)]
            row += [f"{Xz[i, j]:.6g}" if Xz[i, j] == Xz[i, j] else "NA" for j in range(3)]
            row += [g(pr, "svm_score"),
                    g(sq, "intron_len"), g(sq, "was_capped"), g(sq, "up_flank"),
                    g(sq, "intron_seq"), g(sq, "down_flank"),
                    cf_truth, cf_nsp, cf_K, cf_sc, cf_sp]
            w.writerow(row)
    print(f"  wrote {len(names)} rows -> {path}", flush=True)

print("Writing indices...", flush=True)
emit(OUT_TRAIN, tnames, tspecs, tclade, tXz, tXraw, ty, tprov)
emit(OUT_HOLD,  hnames, hspecs, hclade, hXz, hXraw, hy, hprov)
print("DONE.", flush=True)
