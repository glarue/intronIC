"""Tests for the report-only per-genome population-coherence confidence signals.

These signals (background-logo IC + putative-set terminus coherence) are corpus-calibrated
(2026-07) to flag genomes whose splicing signal is too poorly defined to trust the count-based
z_excess gate, or whose strong-call set is a terminus smear rather than a coherent population.
See docs/background_coherence_confidence_investigation.md.
"""
import numpy as np
import pytest

from intronIC.scoring import population_coherence as pc


# ---- metric primitives -------------------------------------------------------------------

def test_logo_ic_identical_sequences_is_maximal():
    # a perfectly conserved motif carries 2 bits/column
    seqs = ["ACGTACGTACGT"] * 500
    assert pc.logo_information_content(seqs, 12) == pytest.approx(2.0 * 12, abs=1e-6)


def test_logo_ic_random_sequences_is_near_zero():
    rng = np.random.default_rng(0)
    seqs = ["".join(rng.choice(list("ACGT"), 12)) for _ in range(4000)]
    # bias-corrected IC of uniform-random columns sits just above 0
    assert pc.logo_information_content(seqs, 12) < 0.5


def test_logo_ic_too_few_sequences_is_nan():
    assert np.isnan(pc.logo_information_content(["ACGTACGTACGT"] * 3, 12))


def test_effective_termini_single_vs_uniform():
    assert pc.effective_termini(["GT-AG"] * 300) == pytest.approx(1.0, abs=1e-6)
    # 8 termini in equal proportion -> ~8 effective termini
    eight = [d for d in ("GT-AG", "AT-AC", "GT-CC", "CT-GC", "GA-AG", "GT-AT", "CT-AC", "AT-GC")]
    eff = pc.effective_termini(eight * 40)
    assert 7.0 < eff < 8.5


def test_top2_fraction():
    dnts = ["GT-AG"] * 70 + ["AT-AC"] * 25 + ["CT-GC"] * 5
    assert pc.top2_terminus_fraction(dnts) == pytest.approx(0.95, abs=1e-6)


def test_canonical_log2fc_depletion_is_negative():
    # putative set avoids the canonical termini the background is rich in -> strongly negative
    put = ["GT-CC", "CT-GC", "GA-AG", "GT-AT"] * 25
    bg = ["GT-AG"] * 900 + ["AT-AC"] * 100
    assert pc.canonical_log2fc(put, bg) < -2.0


# ---- end-to-end classification on synthetic genomes --------------------------------------

def _bg(n, coherent):
    """Synthesize a U2 background: coherent -> crisp GT-AG logos; incoherent -> random smear."""
    rng = np.random.default_rng(1)
    if coherent:
        # conserved at the columns STRUCTURE_WEIGHTS peaks on: donor +2..+5, branch-A (BP +6..+9),
        # acceptor around 3'SS +3..+6 — so IC_sw clears the real-calibrated BG_ICSW_FLOOR.
        s5 = ["GTATCC" + "".join(rng.choice(list("ACGT"), 6)) for _ in range(n)]           # crisp donor core
        sbp = ["".join(rng.choice(list("ACGT"), 5)) + "CTAACG" + "".join(rng.choice(list("ACGT"), 1)) for _ in range(n)]
        s3 = ["".join(rng.choice(list("ACGT"), 3)) + "TCAGGT" + "".join(rng.choice(list("ACGT"), 1)) for _ in range(n)]  # acceptor AG at +3..+5
        dn = ["GT-AG"] * n
    else:
        s5 = ["".join(rng.choice(list("ACGT"), 12)) for _ in range(n)]              # structureless
        sbp = ["".join(rng.choice(list("ACGT"), 12)) for _ in range(n)]
        s3 = ["".join(rng.choice(list("ACGT"), 10)) for _ in range(n)]
        dn = ["".join(rng.choice(list("ACGT"), 2)) + "-" + "".join(rng.choice(list("ACGT"), 2)) for _ in range(n)]
    return s5, sbp, s3, dn


def _make_genome(n_bg=3000, n_strong=120, bg_coherent=True, strong_smear=False):
    b5, bbp, b3, bdn = _bg(n_bg, bg_coherent)
    rng = np.random.default_rng(2)
    if strong_smear:
        p5 = ["".join(rng.choice(list("ACGT"), 12)) for _ in range(n_strong)]
        pbp = ["".join(rng.choice(list("ACGT"), 12)) for _ in range(n_strong)]
        p3 = ["".join(rng.choice(list("ACGT"), 10)) for _ in range(n_strong)]
        pdn = ["".join(rng.choice(list("ACGT"), 2)) + "-" + "".join(rng.choice(list("ACGT"), 2)) for _ in range(n_strong)]
    else:
        p5 = ["GTATCCTTT" + "".join(rng.choice(list("ACGT"), 3)) for _ in range(n_strong)]  # U12-like donor
        pbp = ["".join(rng.choice(list("ACGT"), 5)) + "TTTGAAT" for _ in range(n_strong)]
        p3 = ["".join(rng.choice(list("ACGT"), 8)) + "AC" for _ in range(n_strong)]
        pdn = ["AT-AC"] * n_strong
    p_motif = np.concatenate([np.full(n_strong, 0.99), np.full(n_bg, 0.01)])
    seq5 = np.array(p5 + b5, dtype=object); seqbp = np.array(pbp + bbp, dtype=object)
    seq3 = np.array(p3 + b3, dtype=object); dnts = np.array(pdn + bdn, dtype=object)
    return p_motif, seq5, seqbp, seq3, dnts


def test_coherent_genome_is_unflagged():
    r = pc.compute_from_records(*_make_genome(bg_coherent=True, strong_smear=False))
    assert r.background_logo_ic > pc.BG_LOGO_IC_FLOOR
    assert r.terminus_eff < pc.TERMINUS_EFF_CEIL
    assert not r.splicing_signal_low and not r.terminus_smear
    assert r.confidence == "coherent"


def test_fuzzy_background_flags_low_signal():
    r = pc.compute_from_records(*_make_genome(bg_coherent=False, strong_smear=False))
    assert r.background_logo_ic < pc.BG_LOGO_IC_FLOOR
    assert r.splicing_signal_low
    assert "low_signal" in r.confidence


def test_terminus_smear_flags_smear():
    r = pc.compute_from_records(*_make_genome(bg_coherent=True, strong_smear=True))
    assert r.terminus_eff > pc.TERMINUS_EFF_CEIL
    assert r.terminus_smear
    assert "smear" in r.confidence


def test_dino_like_genome_flags_both():
    # fuzzy background AND smeared strong calls == the A120 signature
    r = pc.compute_from_records(*_make_genome(bg_coherent=False, strong_smear=True))
    assert r.splicing_signal_low and r.terminus_smear
    assert r.confidence == "low_signal+smear"


def test_terminus_not_assessed_below_min_strong():
    # too few strong calls -> terminus unassessed, but background prior still classifies
    p_motif, s5, sbp, s3, dn = _make_genome(n_strong=5, bg_coherent=True)
    r = pc.compute_from_records(p_motif, s5, sbp, s3, dn)
    assert np.isnan(r.terminus_eff)
    assert not r.terminus_smear
    assert r.confidence == "coherent"   # from the background axis alone


# ---- improved metrics: structure-weighted IC + annotation QC -----------------------------

def test_classify_flag_keys_on_icsw_floor():
    # splicing_signal_low is driven by background_ic_sw vs BG_ICSW_FLOOR (not the legacy absolute IC)
    low, _, conf = pc._classify(pc.BG_ICSW_FLOOR - 5, float("nan"), 0, pc.MIN_U2_FOR_BACKGROUND)
    assert low and conf == "low_signal"
    high, _, conf2 = pc._classify(pc.BG_ICSW_FLOOR + 5, float("nan"), 0, pc.MIN_U2_FOR_BACKGROUND)
    assert not high and conf2 == "coherent"


def test_structure_weighted_ic_orders_structured_above_noise():
    b5c, bbpc, b3c, _ = _bg(4000, coherent=True)
    b5n, bbpn, b3n, _ = _bg(4000, coherent=False)
    ic_clean = pc.structure_weighted_ic(b5c, bbpc, b3c)
    ic_noise = pc.structure_weighted_ic(b5n, bbpn, b3n)
    assert ic_clean > ic_noise
    assert np.isnan(pc.structure_weighted_ic(["GTAAGT"] * 3, ["A"] * 3, ["A"] * 3))  # too few


def test_canonical_terminus_fraction_and_quality():
    clean = ["GT-AG"] * 970 + ["GC-AG"] * 20 + ["AT-AC"] * 10
    assert pc.canonical_terminus_fraction(clean) == pytest.approx(1.0)
    assert pc.annotation_quality(pc.canonical_terminus_fraction(clean)) == "clean"
    polluted = ["GT-AG"] * 620 + ["".join(__import__("random").Random(i).choices("ACGT", k=2)) + "-" +
                "".join(__import__("random").Random(i + 1).choices("ACGT", k=2)) for i in range(380)]
    cf = pc.canonical_terminus_fraction(polluted)
    assert cf < 0.85
    assert pc.annotation_quality(cf) == "poor"
    assert pc.annotation_quality(0.88) == "suspect"


def test_substrate_confidence_weight_is_soft_ramp():
    assert pc.substrate_confidence_weight(0.98) == pytest.approx(1.0)   # clean -> full confidence
    assert pc.substrate_confidence_weight(0.62) < 0.15                  # sugarcane-like -> near zero
    assert 0.0 <= pc.substrate_confidence_weight(0.80) <= 1.0           # monotone in-between
    assert pc.substrate_confidence_weight(float("nan")) == 1.0          # unknown -> no penalty


def test_substrate_quality_block_shape():
    r = pc.compute_from_records(*_make_genome(bg_coherent=True, strong_smear=False))
    q = r.substrate_quality()
    assert set(q) == {"canonical_terminus_fraction", "annotation_quality", "background_structure_ic",
                      "splicing_signal", "population_confidence", "substrate_confidence"}
    assert q["splicing_signal"] in ("well_defined", "degenerate")
