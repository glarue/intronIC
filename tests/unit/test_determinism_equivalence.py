"""Step 6 validation: determinism + streaming-vs-in-memory equivalence of the full gate-gapfrac chain.

Hard requirement: same input → bit-identical output per run, AND independent of intron iteration order
(the streaming path feeds contigs via imap_unordered; in-memory feeds all at once). The species-level
stats (gap_fraction, UCL, csig) and the per-intron π_species+discount must be invariant to row order.

This exercises the actual functions end-to-end on a synthetic score_info:
  compute_valley_depth (gate stats) → evaluate_gate (route) → _confidence_shrunk_pi_species /
  compute_adjusted_score (π_species) → apply_continuous_per_intron_discount (file rewrite).

Run: PYTHONPATH=<wt>/src <env>/bin/python3 tests/unit/test_determinism_equivalence.py
"""
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from intronIC.scoring.cluster_validation import compute_valley_depth, compute_adjusted_score
from intronIC.scoring.mode_separation import evaluate_gate
from intronIC.classification.mode_sep_pipeline import apply_continuous_per_intron_discount

COLS = ["name", "svm_score", "rel_score", "adjusted_score", "ensemble_sigma",
        "5'_raw", "bp_raw", "3'_raw", "5'_z", "bp_z", "3'_z", "first_pass_svm"]


def _synth(n=600, seed=0):
    """A realistic-ish species: a U12 cluster + a U2 bulk, with motif + z columns."""
    rng = np.random.default_rng(seed)
    n12 = 80
    rows = []
    for i in range(n):
        u12 = i < n12
        if u12:
            svm = float(np.clip(rng.normal(96, 3), 0, 100))
            raw5, rawb, raw3 = rng.normal(8, 2), rng.normal(6, 2), rng.normal(2, 1)
            z5, zb, z3 = rng.normal(3, 0.5), rng.normal(2.5, 0.5), rng.normal(1.5, 0.5)
        else:
            svm = float(np.clip(rng.normal(20, 15), 0, 100))
            raw5, rawb, raw3 = rng.normal(-2, 2), rng.normal(-1, 2), rng.normal(-1, 1)
            z5, zb, z3 = rng.normal(0, 0.5), rng.normal(0, 0.5), rng.normal(0, 0.5)
        rows.append([f"intron_{i:04d}", round(svm, 4), round(svm - 90, 4), "NA",
                     round(float(rng.uniform(2, 8)), 4),
                     round(raw5, 4), round(rawb, 4), round(raw3, 4),
                     round(z5, 4), round(zb, 4), round(z3, 4), round(svm, 4)])
    return pd.DataFrame(rows, columns=COLS)


def _write(df, path):
    df.to_csv(path, sep="\t", index=False)


def _full_chain(df, d, tag):
    """Run gate stats + π_species + discount on a (possibly reordered) df; return the per-intron
    adjusted_score keyed by name (order-independent dict for comparison)."""
    # species stats from the confident cluster vs U2 bulk (Fisher axis)
    svm = df["svm_score"].to_numpy()
    z = np.column_stack([df["5'_z"], df["bp_z"], df["3'_z"]])
    conf, u2 = svm >= 90.0, svm < 90.0
    vr = compute_valley_depth(z[conf], z[u2])
    ucl = vr["gap_fraction_ucl"]
    sigma = df["ensemble_sigma"].to_numpy()
    # π_species (k_σ=0, modesep route) then discount via the real file rewriter
    pi = np.array([compute_adjusted_score(float(s), ucl, float(g), k_sigma=0.0)
                   for s, g in zip(svm, sigma)])
    df2 = df.copy()
    df2["adjusted_score"] = pi
    p = Path(d) / f"{tag}.score_info.iic"
    _write(df2, p)
    apply_continuous_per_intron_discount(p, threshold=90.0, input_column="adjusted_score")
    out = pd.read_csv(p, sep="\t")
    return (vr["gap_fraction"], vr["gap_fraction_ucl"], vr.get("centroid_sigma"),
            dict(zip(out["name"], out["adjusted_score"].round(6))))


def test_determinism_and_order_invariance():
    d = tempfile.mkdtemp()
    try:
        df = _synth()
        # run A: original order
        gfA, uclA, csA, adjA = _full_chain(df, d, "A")
        # run B: same input again (determinism)
        gfB, uclB, csB, adjB = _full_chain(df, d, "B")
        # run C: shuffled row order (streaming-vs-in-memory equivalence)
        df_shuf = df.sample(frac=1.0, random_state=12345).reset_index(drop=True)
        gfC, uclC, csC, adjC = _full_chain(df_shuf, d, "C")

        # species stats bit-identical across all three
        assert gfA == gfB == gfC, ("gap_fraction differs", gfA, gfB, gfC)
        assert uclA == uclB == uclC, ("UCL differs", uclA, uclB, uclC)
        assert csA == csB == csC, ("centroid_sigma differs", csA, csB, csC)

        # per-intron adjusted_score identical run-to-run
        assert adjA == adjB, "non-deterministic: adjusted_score differs on identical reruns"
        # per-intron adjusted_score identical regardless of row order (join by name)
        names = set(adjA) | set(adjC)
        maxd = max(abs(adjA[n] - adjC[n]) for n in names)
        assert maxd == 0.0, f"order-dependent: max |Δadjusted| = {maxd}"

        print(f"  gap_fraction={gfA:.6f}  UCL={uclA:.6f}  csig={csA:.4f}")
        print(f"  determinism: identical reruns ✓ | order-invariance: max Δ={maxd} ✓")
        print(f"  introns={len(adjA)}  called(adj>=90)={sum(1 for v in adjA.values() if v >= 90)}")
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_gate_route_order_invariant():
    """The gate decision (route + carried stats) must not depend on candidate row order."""
    rng = np.random.default_rng(3)
    u2 = rng.normal(0, 1, (1000, 3))
    u12 = rng.normal(3.0, 0.7, (50, 3))
    g1 = evaluate_gate(u12, u2, mu_u12_5p_raw=15.0, n_eff_candidates=50.0,
                       separation_5p=2.0, valley_depth_fn=compute_valley_depth)
    perm = rng.permutation(len(u12))
    g2 = evaluate_gate(u12[perm], u2, mu_u12_5p_raw=15.0, n_eff_candidates=50.0,
                       separation_5p=2.0, valley_depth_fn=compute_valley_depth)
    assert g1.passes == g2.passes and g1.reason == g2.reason
    assert abs(g1.gap_fraction - g2.gap_fraction) < 1e-12, (g1.gap_fraction, g2.gap_fraction)
    assert abs(g1.gap_fraction_ucl - g2.gap_fraction_ucl) < 1e-12
    assert abs(g1.centroid_sigma - g2.centroid_sigma) < 1e-12
    print(f"  gate route '{g1.reason}' order-invariant: gf={g1.gap_fraction:.6f} "
          f"ucl={g1.gap_fraction_ucl:.6f} csig={g1.centroid_sigma:.4f} ✓")


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"ALL {len(fns)} DETERMINISM/EQUIVALENCE CHECKS PASSED")
