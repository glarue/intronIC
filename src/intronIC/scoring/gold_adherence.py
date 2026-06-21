"""Absolute U12-adherence feature (gold expectation).

Ported verbatim (logic) from the validated corpus tool
``multispecies_v23/gold_adherence.py`` (ADHERENCE_DECISION.md, 2026-06-09). This
is the *scoring-side* copy: it must compute byte-identical adh5/adhbp to the
corpus builder (``add_adherence_features.py``) so that the features the model is
trained on match the features it is scored with. See RETRAIN_PLAN.md ("Feature
parity train vs score: the #1 risk").

The gold PWMs are FULL human-PWM span (5' -20..+19, BP -9..+2), so scored
sub-coords are selected by *position label*, not array index — the feature is
window-agnostic. Give it any 5' window plus the coordinate of its first base
(``five_start``) and it scores the configured ``core`` labels by matching
position labels to the PWM. The intronIC scorer's 5' region (``five_coords``,
default ``(-3, 9)``) feeds it directly with ``five_start = five_coords[0]``.

    adh5  = sum over the selected 5' core labels of log2( gold_U12_freq / 0.25 )
    adhbp = sum over the BP positions of log2( gold_U12_freq / 0.25 )

Dinucleotide handling (mirrors intronIC base scoring, code-traced):
  - canonical set = {(GT,AG), (GC,AG), (AT,AC)}      (extraction/sequences.py)
  - PWM class via select_best: atac iff canonical AT-AC; else gtag (incl. NC)
  - ignore_nc_dnts=True: for NON-canonical introns, skip 5' positions {0,+1}

This is a FEATURE for the SVM, not a hard gate; base scoring stays on the human
PWMs.
"""
import json
import math
from pathlib import Path
from typing import Optional

CANONICAL_PAIRS = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
DEFAULT_CORE = ('0', '+1', '+2', '+3', '+4', '+5')   # 5' core sub-coords (selectable/shiftable)
DNT_5P_LABELS = {'0', '+1'}                           # the 5' dinucleotide (skipped if non-canonical)
LOG0 = 1e-6

# Packaged gold-expectation PWM (shipped alongside the model bundle).
DEFAULT_GOLD_PWM = Path(__file__).resolve().parent.parent / "data" / "gold_u12_expectation_pwms.json"


class GoldAdherence:
    """Score an intron's absolute adherence to the corpus gold-U12 expectation."""

    def __init__(self, gold_json_path=None, core=DEFAULT_CORE, ignore_nc_dnts=True):
        if gold_json_path is None:
            gold_json_path = DEFAULT_GOLD_PWM
        g = json.load(open(gold_json_path))
        if not isinstance(g.get("five_prime"), dict) or "gtag" not in g["five_prime"]:
            raise ValueError(
                "gold json is not dinucleotide-split; rebuild with build_gold_u12_pwms.py"
            )
        self.path = str(gold_json_path)
        self.five = g["five_prime"]          # {"gtag":[N {base:freq}], "atac":[N ...]}
        self.bp = g["branch_point"]
        self.P5 = g["positions_5p"]           # full span, e.g. -20..+19
        self.BPK = g["positions_bp"]          # -9..+2
        self.core = list(core)
        self.ignore_nc_dnts = ignore_nc_dnts
        self._p5idx = {lab: i for i, lab in enumerate(self.P5)}
        missing = [c for c in self.core if c not in self._p5idx]
        if missing:
            raise ValueError(
                f"core labels {missing} not in gold PWM span {self.P5[0]}..{self.P5[-1]}"
            )

    @classmethod
    def from_default(cls, **kwargs) -> Optional["GoldAdherence"]:
        """Load the packaged gold PWM; return None (with no error) if absent.

        Lets the scorer enrich output with adherence whenever the bundle ships
        the PWM, while degrading gracefully for older bundles / sequence-only
        contexts that have no gold expectation file.
        """
        if not DEFAULT_GOLD_PWM.exists():
            return None
        return cls(DEFAULT_GOLD_PWM, **kwargs)

    @staticmethod
    def _lab_to_int(lab):
        return int(lab.replace("+", ""))

    @staticmethod
    def classify(dnt5, dnt3):
        """(pwm_class, is_canonical) from the splice dinucleotides, mirroring select_best."""
        d5, d3 = dnt5.upper(), dnt3.upper()
        is_canonical = (d5, d3) in CANONICAL_PAIRS
        cls = "atac" if (d5, d3) == ("AT", "AC") else "gtag"
        return cls, is_canonical

    def adh5(self, five_seq, five_start, cls, is_canonical):
        """Score the core labels. ``five_seq`` begins at coordinate ``five_start``
        (e.g. -20 for the introns 40nt window, -3 for the scorer's 5' region)."""
        s = five_seq.upper()
        M = self.five[cls]
        skip = (not is_canonical) and self.ignore_nc_dnts
        total = 0.0
        for lab in self.core:
            if skip and lab in DNT_5P_LABELS:
                continue
            pi = self._p5idx[lab]
            si = self._lab_to_int(lab) - five_start   # offset into five_seq
            if 0 <= si < len(s) and s[si] in "ACGT":
                total += math.log2(max(M[pi].get(s[si], LOG0), LOG0) / 0.25)
        return total

    def adhbp(self, bp_seq, cls):
        s = bp_seq.upper()
        M = self.bp[cls]
        return sum(
            math.log2(max(M[i].get(s[i], LOG0), LOG0) / 0.25)
            for i in range(min(len(s), len(self.BPK))) if s[i] in "ACGT"
        )

    def compute(self, five_seq, five_start, bp_seq, dnt5, dnt3):
        """Return dict(adh5, adhbp, dnt_class, is_canonical) for one intron.

        five_seq: any 5' window; five_start: coord of its first base.
        bp_seq:   the located branch-point window (None -> adhbp is None).
        dnt5/dnt3: the 5'/3' splice dinucleotides (intron[:2], intron[-2:])."""
        cls, is_canon = self.classify(dnt5, dnt3)
        adhbp = self.adhbp(bp_seq, cls) if bp_seq else None
        return {
            "adh5": self.adh5(five_seq, five_start, cls, is_canon),
            "adhbp": adhbp,
            "dnt_class": cls,
            "is_canonical": is_canon,
        }


# --- global adh z-scale constants: MUST match the offline _adh_feats that produced the graduated_tail features ---
ADH5_MED, ADH5_IQR = -16.5067, 30.8676
ADHBP_MED, ADHBP_IQR = 3.4674, 8.8282
_DEFAULT_GA = None


def adh_features(df, valid, want5=True, ga=None):
    """Globally z-scaled (adh5z, adhbp_z) for the `valid` rows of a score_info frame.

    Convention-aware 5' window (len>=12 => -3-anchored GT@idx3; len 10 => 0-anchored GT@idx0) + 3' AG/AC slot.
    Mirrors the offline `_adh_feats` (identical logic + constants) used to fit the graduated tail, so the deployed
    adh matches the fit. Returns (adh5z|None, adhbp_z)."""
    import numpy as np
    global _DEFAULT_GA
    if ga is None:
        if _DEFAULT_GA is None:
            _DEFAULT_GA = GoldAdherence.from_default()
        ga = _DEFAULT_GA
    s5 = df["5'_seq"].astype(str).values[valid]
    bp = df["bp_seq"].astype(str).values[valid]
    s3 = df["3'_seq"].astype(str).values[valid]
    a5 = np.full(len(bp), ADH5_MED, float)
    abp = np.full(len(bp), ADHBP_MED, float)
    for i in range(len(bp)):
        x5, xb, x3 = s5[i].upper(), bp[i].upper(), s3[i].upper()
        if not xb or xb == "NAN" or len(x5) < 5:
            continue
        fs = -3 if len(x5) >= 12 else 0
        gt = 3 if fs == -3 else 0
        d5 = x5[gt:gt + 2]
        c3 = [x3[-2:], x3[-6:-4] if len(x3) >= 6 else ""]
        d3 = next((c for c in c3 if c in ("AG", "AC")), c3[1] or c3[0])
        r = ga.compute(x5, fs, xb, d5, d3)
        a5[i] = r["adh5"]
        abp[i] = r["adhbp"]
    a5z = (a5 - ADH5_MED) / ADH5_IQR if want5 else None
    return a5z, (abp - ADHBP_MED) / ADHBP_IQR


if __name__ == "__main__":
    import sys
    ga = GoldAdherence(sys.argv[1] if len(sys.argv) > 1 else None)
    print(f"5' span: {ga.P5[0]}..{ga.P5[-1]} ({len(ga.P5)})  core={ga.core}")
    for cls in ("gtag", "atac"):
        c5 = "".join(max(ga.five[cls][i], key=ga.five[cls][i].get) for i in range(len(ga.P5)))
        print(f"{cls} 5' consensus: {c5}")
