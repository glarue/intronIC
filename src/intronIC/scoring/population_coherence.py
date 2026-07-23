"""Per-genome CONFIDENCE signals on the U12-like (strong-call) population — complementary to, and
independent of, ``z_excess`` (``species_adjudicator.py``).

MOTIVATION
----------
``z_excess`` is a COUNT statistic: it asks whether the strong-call count exceeds what the genome's own
U2 tail predicts under an exponential null. It says nothing about whether that null is *trustworthy* or
whether the strong-call set is a coherent population. On genomes whose splicing signal is poorly defined
(e.g. dinoflagellates: the U2-type logos are themselves a near-structureless smear), a "U12-like" set is
far more likely to be the noisy edge of that smear than a distinct minor-spliceosome population — yet it
can still clear ``z_excess`` (Amoebophrya A120: z=4.27, over the loss ceiling). This module quantifies the
"is there even a distinct, coherent population here, and is the genome's splicing signal well-defined
enough to trust the count-based gate?" intuition.

TWO ORTHOGONAL, SELF-CONTAINED, CLADE-BLIND SIGNALS (corpus correlation ~-0.16)
------------------------------------------------------------------------------
1. ``background_logo_ic`` (SEQUENCE-side PRIOR) — total information content (bits) of the U2-type
   background splice logos (5'SS + BP + 3'SS). Low IC == the genome's *major*-spliceosome signal is a
   smear, so the prior that a distinct minor population can be found above it is low. NOT selection-biased
   (the background is the bulk, not chosen for motif strength). Available even with zero strong calls.
2. ``terminus_eff`` (PUTATIVE-side) — effective number of terminal dinucleotides (2**Miller-Madow entropy)
   in the strong-call set. A genuine U12 population — canonical (GT-AG/AT-AC) OR a divergent non-canonical
   one — collapses onto a FEW coherent termini; a compositional smear spreads ~uniformly across many. This
   is definition-free (does not privilege canonical termini), so it does not penalize divergent bearers.

CALIBRATION (see docs/background_coherence_confidence_investigation.md; artifacts in
eval_corpus/population_coherence_2026-07/). Thresholds set on a 2026-07 corpus of 1715 snRNA-confirmed
bearers (busco>=85, n_strong>=20) + 37 reliable divergent losses (snRNA-absent, reliable U12/U6atac axis):
  - background_logo_ic: every bearer AND every loss >= 15.56; all three dino suspects <= 14.55. FLOOR=15.0
    flags the dinos and 0/1752 real genomes (incl. the most divergent oomycetes, Saprolegnia ~15.6).
  - terminus_eff: max over all bearers = 9.48 (Macrostomum, a real divergent-terminus flatworm); A120=54.
    CEIL=12.0 flags A120's smear with a wide margin and 0/1715 bearers.

REPORT-ONLY BY DESIGN. This is a corroboration/confidence axis, NOT a gate: it does not change ``type_id``
or ``motif_category``. It flags *how much to trust* a detection (and whether the population is coherent),
not whether an intron is U12. It does not discriminate bearer-vs-loss (a real loss has a crisp GT-AG
background too); it discriminates well-defined-splicing from structureless-splicing.

The score-side alternative — testing/replacing the exponential U2-tail (heavier-tailed model) — was
evaluated and is a NULL result: A120's margin tail is bounded/light (cv_excess=0.76 < 1), and cv_excess is
~constant corpus-wide (0.84-0.96). Tail SHAPE does not separate these cases; sequence coherence does.

FIRST-CLASS WIRE-IN POINT (not done here to keep the v3.0.0 release behavior unchanged):
``species_adjudicator.apply_pmotif_adjudication`` already holds the score_info ``df`` (columns ``5'_seq``,
``bp_seq``, ``3'_seq``) and the computed ``p_motif``; ``dnts`` is derivable there or read from the meta
writer. Calling :func:`compute_from_records` there and stamping the fields into the tail-model sidecar /
metrics would promote this to a per-genome reported field alongside ``motif_category``.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Optional, Sequence
import numpy as np

# --- calibrated constants (corpus 2026-07; see module docstring + investigation doc) ---
BG_ICSW_FLOOR: float = 36.5        #: PRIMARY: background_ic_sw (composition-relative, structure-weighted)
                                   #: below this => splicing signal degenerate. Corpus-calibrated 2026-07 on
                                   #: 1872 clean well-annotated bearers (corroborated, BUSCO>=80, canonical
                                   #: >=0.92; full-background/uncapped IC): their IC_sw min=37.52, median 44.0.
                                   #: The dinoflagellate substrates top out ~35.9, so 36.5 sits centered in the
                                   #: empty ~35.9-37.52 gap: 0/1872 clean-bearer false flags, all 10 dinos caught.
BG_LOGO_IC_FLOOR: float = 15.0     #: LEGACY (absolute 2-H IC) — superseded by BG_ICSW_FLOOR; kept for the
                                   #: still-reported legacy background_logo_ic. Do not use for the flag.
TERMINUS_EFF_CEIL: float = 12.0    #: putative-set effective #termini above this => terminus smear
MIN_STRONG_FOR_TERMINUS: int = 20  #: fewer strong calls than this => terminus coherence not assessed
MIN_U2_FOR_BACKGROUND: int = 1000  #: fewer U2 background introns than this => background IC not assessed

CANONICAL_MINOR_TERMINI = frozenset({"GT-AG", "AT-AC"})  #: canonical minor-spliceosome termini

# --- annotation-quality QC via terminal dinucleotides (canonical splice termini) ---
# Empirically grounded (2026-07 corpus scan of 1873 snRNA-corroborated bearers, BUSCO>=80): canonical
# terminus fraction is >=0.90 in 99.6% of REAL bearers across EVERY divergent clade (oomycetes median 0.994,
# Mucoromycota 0.989, nematodes 0.997). Genuine biological terminal drift is effectively unobserved; sub-0.85
# genomes are annotation pollution (polyploid grasses, shared-pipeline artifacts). So low canonical fraction
# is a safe (but SOFT, not hard) annotation-quality signal — see docs/background_ic_annotation_qc_*.md.
QC_CANONICAL_TERMINI = frozenset({"GT-AG", "GC-AG", "AT-AC"})  #: spliceosomal-canonical termini for QC
GTAG_CLEAN: float = 0.92   #: canonical terminus fraction >= this => annotation "clean"
GTAG_POOR: float = 0.85    #: canonical terminus fraction < this => annotation "poor" (pollution likely)

# Structure weights for the composition-relative background IC (structure_weighted_ic): the expected
# per-position splice information, learned as the mean background KL profile of 5 clean reference bearers
# (H. sapiens, A. thaliana, S. bicolor, M. musculus, D. rerio), normalized to mean 1.0. Layout matches the
# score_info windows 5'_seq(12) + bp_seq(12) + 3'_seq(10); peaks at the donor core, the branch-A, and the
# acceptor AG, ~0 at uninformative flanks. Re-derive with eval_corpus/population_coherence_2026-07 tooling.
STRUCTURE_WEIGHTS = (
    # 5'SS (12): donor core (+3/+4) dominates
    0.479, 0.874, 1.935, 3.973, 2.915, 1.137, 0.892, 1.234, 0.189, 0.058, 0.081, 0.088,
    # BP (12): branch-A region (+7 / +9)
    0.045, 0.097, 0.080, 0.155, 0.327, 0.521, 1.166, 2.776, 0.136, 2.860, 0.303, 0.103,
    # 3'SS (10): acceptor AG (+4/+5)
    0.459, 0.647, 0.127, 2.069, 3.618, 3.964, 0.587, 0.039, 0.015, 0.052,
)


def _mm_entropy(counts: np.ndarray) -> float:
    """Shannon entropy (bits) of a category-count vector with Miller-Madow small-sample bias correction."""
    counts = np.asarray(counts, float)
    N = counts.sum()
    if N <= 1:
        return 0.0
    p = counts[counts > 0] / N
    H = -(p * np.log2(p)).sum()
    K = int((counts > 0).sum())
    return float(H + (K - 1) / (2 * N * np.log(2)))


def logo_information_content(seqs: Sequence[str], length: int) -> float:
    """Total bias-corrected information content (bits) of aligned, equal-length motif sequences.

    Per column: IC = max(0, 2 - H_MM), summed over ``length`` columns (max = 2*length). Non-ACGT symbols
    are ignored per column. Returns NaN with fewer than 8 usable sequences. Computed over ALL sequences (no
    subsampling) so the result is ORDER-INDEPENDENT — required for streaming/in-memory parity, since the two
    classify paths may emit score_info rows in different order (per-column counts are commutative).
    """
    seqs = [s for s in seqs if isinstance(s, str) and len(s) >= length]
    n = len(seqs)
    if n < 8:
        return float("nan")
    arr = np.frombuffer("".join(s[:length] for s in seqs).encode(), dtype="S1").reshape(n, length)
    total = 0.0
    for j in range(length):
        col = arr[:, j]
        cnt = np.array([(col == b).sum() for b in (b"A", b"C", b"G", b"T")], float)
        if cnt.sum() < 1:
            continue
        total += max(0.0, 2.0 - _mm_entropy(cnt))
    return float(total)


def _col_counts(seqs: Sequence[str], length: int):
    """length x 4 ACGT count matrix over aligned equal-length sequences (None if < 8 usable). Over ALL
    sequences (no subsampling) — order-independent, for streaming/in-memory parity (see logo_information_content)."""
    seqs = [s[:length] for s in seqs if isinstance(s, str) and len(s) >= length]
    if len(seqs) < 8:
        return None
    n = len(seqs)
    arr = np.frombuffer("".join(seqs).encode(), dtype="S1").reshape(n, length)
    return np.stack([(arr == b).sum(0) for b in (b"A", b"C", b"G", b"T")], 1).astype(float)


def structure_weighted_ic(seq5: Sequence[str], seqbp: Sequence[str], seq3: Sequence[str],
                          weights: Sequence[float] = STRUCTURE_WEIGHTS) -> float:
    """Composition-relative, structure-weighted background information content (bits).

    Two corrections over :func:`logo_information_content` (absolute ``2 - H`` vs uniform), motivated by our
    priors about intron motifs (see docs/background_ic_annotation_qc_investigation.md):
      1. **Composition-relative**: each column's information is KL(p_col || q), where q is the genome's OWN
         background base composition (pooled over the motif windows) — not uniform 0.25. Removes GC/AT-skew
         inflation and matches the pipeline's background-correction convention (``background.py``).
      2. **Structure-weighted**: columns are weighted by ``STRUCTURE_WEIGHTS`` (expected per-position splice
         information), so a genuinely-divergent-but-real bearer with a weak *extended* consensus but intact
         *anchors* is not dragged to the noise floor by its uninformative positions.
    Windows: 5'SS(12) + BP(12) + 3'SS(10) = 34 columns. NaN if any window has < 8 usable sequences.
    """
    counts = [(_col_counts(seq5, 12), 12), (_col_counts(seqbp, 12), 12), (_col_counts(seq3, 10), 10)]
    if any(c is None for c, _ in counts):
        return float("nan")
    pooled = np.sum([c.sum(0) for c, _ in counts], 0)
    q = pooled / pooled.sum()
    kl = []
    for c, L in counts:
        p = c / np.maximum(c.sum(1, keepdims=True), 1)
        for j in range(L):
            pj = p[j]; nz = pj > 0
            kl.append(max(0.0, float((pj[nz] * np.log2((pj[nz] + 1e-9) / (q[nz] + 1e-9))).sum())))
    w = np.asarray(weights, float)
    return float((w * np.asarray(kl)).sum())


def canonical_terminus_fraction(dnts: Sequence[str]) -> float:
    """Fraction of introns with a spliceosomal-canonical terminus (GT-AG / GC-AG / AT-AC). Annotation-QC:
    <~0.85 signals annotation pollution (junk introns with random boundaries), empirically not real drift."""
    d = [x for x in dnts if isinstance(x, str) and x]
    if not d:
        return float("nan")
    return float(np.mean([x in QC_CANONICAL_TERMINI for x in d]))


def annotation_quality(canonical_fraction: float) -> str:
    """clean / suspect / poor / unassessed from the canonical terminus fraction (thresholds GTAG_CLEAN/POOR)."""
    if not np.isfinite(canonical_fraction):
        return "unassessed"
    if canonical_fraction >= GTAG_CLEAN:
        return "clean"
    if canonical_fraction < GTAG_POOR:
        return "poor"
    return "suspect"


def substrate_confidence_weight(canonical_fraction: float) -> float:
    """SOFT [0,1] confidence multiplier from annotation quality: 1.0 at canonical>=0.92, ramping to 0 at 0.60
    (~40% junk). Soft by design (never a hard suppressor) — a real bearer on a polluted assembly is still a
    real bearer; this only down-weights *confidence* in its count-based adjudication."""
    if not np.isfinite(canonical_fraction):
        return 1.0
    return float(np.clip((canonical_fraction - 0.60) / (GTAG_CLEAN - 0.60), 0.0, 1.0))


def effective_termini(dnts: Sequence[str]) -> float:
    """Effective number of terminal dinucleotides = 2**(Miller-Madow entropy of the dnt distribution)."""
    vals, counts = np.unique([d for d in dnts if isinstance(d, str) and d], return_counts=True)
    if counts.sum() < 1:
        return float("nan")
    return float(2 ** _mm_entropy(counts))


def top2_terminus_fraction(dnts: Sequence[str]) -> float:
    """Fraction of the set explained by its two commonest terminal dinucleotides (n-robust coherence)."""
    vals, counts = np.unique([d for d in dnts if isinstance(d, str) and d], return_counts=True)
    tot = counts.sum()
    if tot < 1:
        return float("nan")
    return float(np.sort(counts)[-2:].sum() / tot)


def canonical_log2fc(put_dnts: Sequence[str], bg_dnts: Sequence[str]) -> float:
    """log2 fold-change of the canonical-minor terminus fraction (putative vs the genome's OWN background).

    Background-calibrated: strongly negative == the strong-call set is depleted in canonical termini
    relative to the genome's own intron pool. Note: this axis alone would false-flag a genuinely
    non-canonical divergent bearer; use it together with :func:`effective_termini` (coherence), which does
    not. Reported as a diagnostic, not a flag.
    """
    put = [d for d in put_dnts if isinstance(d, str) and d]
    bg = [d for d in bg_dnts if isinstance(d, str) and d]
    if not put or not bg:
        return float("nan")
    cfp = np.mean([d in CANONICAL_MINOR_TERMINI for d in put])
    cfb = np.mean([d in CANONICAL_MINOR_TERMINI for d in bg])
    eps = 0.5 / len(put)
    return float(np.log2((cfp + eps) / (cfb + eps)))


@dataclass
class PopulationConfidence:
    """Report-only per-genome confidence signals on the U12-like population (see module docstring)."""
    n_strong: int
    n_u2: int
    background_logo_ic: float = float("nan")   #: SEQUENCE-side prior: total bg splice-logo IC (bits, ABSOLUTE)
    background_5ss_ic: float = float("nan")     #: bg 5'SS-only IC (bits), the single cleanest component
    background_ic_sw: float = float("nan")      #: PRIMARY bg IC: composition-relative + structure-weighted (bits)
    terminus_eff: float = float("nan")          #: PUTATIVE-side: effective #terminal dinucleotides
    terminus_top2_frac: float = float("nan")    #: fraction explained by the 2 commonest termini
    canonical_log2fc: float = float("nan")      #: diagnostic: canonical-terminus depletion vs own background
    canonical_terminus_fraction: float = float("nan")  #: annotation QC: GT-AG/GC-AG/AT-AC fraction (all introns)
    annotation_quality: str = "unassessed"      #: clean / suspect / poor (from canonical_terminus_fraction)
    substrate_confidence: float = 1.0           #: SOFT [0,1] confidence multiplier from annotation quality
    splicing_signal_low: bool = False           #: background_ic_sw < BG_ICSW_FLOOR (degenerate splicing signal)
    terminus_smear: bool = False                #: terminus_eff > TERMINUS_EFF_CEIL
    confidence: str = "unassessed"              #: coherent / low_signal / terminus_smear / low_signal+smear

    def as_dict(self) -> dict:
        return asdict(self)

    def substrate_quality(self) -> dict:
        """The compact per-genome REPORT block for metrics.iic.json (see module docstring / MVP doc)."""
        return {
            "canonical_terminus_fraction": self.canonical_terminus_fraction,
            "annotation_quality": self.annotation_quality,
            "background_structure_ic": self.background_ic_sw,
            "splicing_signal": ("degenerate" if self.splicing_signal_low else "well_defined"),
            "population_confidence": self.confidence,
            "substrate_confidence": self.substrate_confidence,
        }


def _classify(bg_ic_sw: float, term_eff: float, n_strong: int, n_u2: int) -> tuple[bool, bool, str]:
    # PRIMARY splicing-signal flag on the composition-relative, structure-weighted IC (BG_ICSW_FLOOR).
    low = bool(np.isfinite(bg_ic_sw) and bg_ic_sw < BG_ICSW_FLOOR)
    smear = bool(np.isfinite(term_eff) and term_eff > TERMINUS_EFF_CEIL)
    bg_ok = n_u2 >= MIN_U2_FOR_BACKGROUND and np.isfinite(bg_ic_sw)
    term_ok = n_strong >= MIN_STRONG_FOR_TERMINUS and np.isfinite(term_eff)
    if not bg_ok and not term_ok:
        return low, smear, "unassessed"
    if low and smear:
        conf = "low_signal+smear"
    elif smear:
        conf = "terminus_smear"
    elif low:
        conf = "low_signal"
    else:
        conf = "coherent"
    return low, smear, conf


def compute_from_records(p_motif: np.ndarray,
                         seq5: Sequence[str], seqbp: Sequence[str], seq3: Sequence[str],
                         dnts: Sequence[str],
                         call_threshold: float = 0.9,
                         u2_threshold: float = 0.5) -> PopulationConfidence:
    """Compute the confidence signals from per-intron arrays (pipeline-friendly).

    All array-likes are aligned per intron. ``p_motif`` splits the set: putative == P_motif >= call_threshold
    (the strong calls / z_excess numerator), background == P_motif < u2_threshold (the U2-type bulk).
    ``seq5``/``seqbp``/``seq3`` are the aligned 5'SS/BP/3'SS motif strings (score_info ``5'_seq``/``bp_seq``/
    ``3'_seq``); ``dnts`` the terminal-dinucleotide strings (meta ``dnts``, e.g. "GT-AG").
    """
    p = np.asarray(p_motif, float)
    seq5 = np.asarray(seq5, dtype=object); seqbp = np.asarray(seqbp, dtype=object)
    seq3 = np.asarray(seq3, dtype=object); dnts = np.asarray(dnts, dtype=object)
    put = p >= call_threshold
    bg = p < u2_threshold
    n_strong = int(put.sum()); n_u2 = int(bg.sum())

    res = PopulationConfidence(n_strong=n_strong, n_u2=n_u2)
    # SEQUENCE-side prior (background) — computable whenever there is an ample U2 background.
    if n_u2 >= MIN_U2_FOR_BACKGROUND:
        bg5, bgbp, bg3 = list(seq5[bg]), list(seqbp[bg]), list(seq3[bg])
        ic5 = logo_information_content(bg5, 12)
        icbp = logo_information_content(bgbp, 12)
        ic3 = logo_information_content(bg3, 10)
        res.background_5ss_ic = ic5
        if np.isfinite(ic5) and np.isfinite(icbp) and np.isfinite(ic3):
            res.background_logo_ic = ic5 + icbp + ic3
        # PRIMARY background IC: composition-relative + structure-weighted (see structure_weighted_ic).
        res.background_ic_sw = structure_weighted_ic(bg5, bgbp, bg3)
    # ANNOTATION QC — canonical terminus fraction over the SCORABLE introns (finite P_motif). Restricting to
    # scorable is essential: unscorable introns (N-containing termini from assembly gaps, e.g. "GT-NN"/"NN-AG")
    # are NOT annotation over-prediction and would spuriously tank the fraction on gappy assemblies (tupaia:
    # 0.43 over all meta rows vs 0.94 over the classified u2 set). Genuine pollution (junk gene models with
    # real non-canonical termini, e.g. sugarcane 0.63) survives this restriction.
    scorable_dnts = list(dnts[np.isfinite(p)])
    res.canonical_terminus_fraction = canonical_terminus_fraction(scorable_dnts)
    res.annotation_quality = annotation_quality(res.canonical_terminus_fraction)
    res.substrate_confidence = substrate_confidence_weight(res.canonical_terminus_fraction)
    # PUTATIVE-side coherence — needs enough strong calls for a stable estimate.
    if n_strong >= MIN_STRONG_FOR_TERMINUS:
        res.terminus_eff = effective_termini(list(dnts[put]))
        res.terminus_top2_frac = top2_terminus_fraction(list(dnts[put]))
        res.canonical_log2fc = canonical_log2fc(list(dnts[put]), list(dnts[bg]))

    res.splicing_signal_low, res.terminus_smear, res.confidence = _classify(
        res.background_ic_sw, res.terminus_eff, n_strong, n_u2)
    return res


def compute_from_iic(score_info_path: str, meta_path: Optional[str] = None) -> PopulationConfidence:
    """Post-hoc, REPORT-ONLY convenience: compute the signals from a finished genome's output files.

    Reads ``5'_seq``/``bp_seq``/``3'_seq``/``P_motif`` from ``score_info.iic`` and ``dnts`` from the sibling
    ``meta.iic`` (auto-located if ``meta_path`` omitted). Works on any existing v3 run with no re-classify.
    """
    import pandas as pd
    s = pd.read_csv(score_info_path, sep="\t", dtype=str, keep_default_na=False,
                    usecols=["name", "P_motif", "5'_seq", "bp_seq", "3'_seq"])
    p = pd.to_numeric(s["P_motif"], errors="coerce").to_numpy()
    if meta_path is None:
        meta_path = score_info_path.replace(".score_info.iic", ".meta.iic")
    try:
        m = pd.read_csv(meta_path, sep="\t", dtype=str, keep_default_na=False, usecols=["name", "dnts"])
        dnts = s.merge(m, on="name", how="left")["dnts"].to_numpy()
    except (FileNotFoundError, ValueError):
        dnts = np.array([""] * len(s), dtype=object)   # background-only signal still available
    return compute_from_records(p, s["5'_seq"].to_numpy(), s["bp_seq"].to_numpy(),
                                s["3'_seq"].to_numpy(), dnts)
