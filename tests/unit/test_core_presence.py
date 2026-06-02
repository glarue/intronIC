"""Step 7: the core-presence species-level term (label-free FP suppression for U12-loss species).

The gap_fraction + centroid_sigma gate PASSES U12-loss species (TetThe/Chlamydomonas) whose first-pass
candidate cluster is geometrically real in z-space but motif-empty in raw space. The per-intron discount
can't use that species context (it's per-intron). The core-presence term — w_core keyed on
frac(candidate raw_sum > bar) — supplies it: a species with no strong-motif core gets a strong species
prior down-shift; a real bearer (however diverged) is untouched.

Run: PYTHONPATH=<wt>/src <env>/bin/python3 tests/unit/test_core_presence.py
"""
import math
import os
import numpy as np

from intronIC.scoring.mode_separation import (
    evaluate_gate, GateDecision, DEFAULT_CORE_RAW_SUM_BAR, DEFAULT_CORE_GATE_FLOOR,
)
from intronIC.scoring.cluster_validation import (
    compute_valley_depth, compute_adjusted_score, _confidence_shrunk_pi_species,
)

# Hygiene: the experimental INTRONIC_FORCE_FALLBACK diagnostic hook forces every gate-passer to
# fallback; unset it so the routing assertions below reflect the real gate, not the override.
os.environ.pop("INTRONIC_FORCE_FALLBACK", None)


def test_core_fraction_computed_and_carried():
    """The gate computes core_fraction from the candidate raw points and carries it on GateDecision."""
    rng = np.random.default_rng(0)
    u2 = rng.normal([-25, -8, -4], 3, (500, 3))
    # strong-motif candidates: raw_sum well above the bar
    u12 = rng.normal([9, 6, 3], 2, (80, 3))
    g = evaluate_gate(u12, u2, mu_u12_5p_raw=9.0, n_eff_candidates=80.0,
                      separation_5p=2.0, valley_depth_fn=compute_valley_depth)
    assert g.core_fraction is not None
    assert g.core_fraction > 0.9, g.core_fraction          # strong-motif cluster → core present
    assert "core_fraction" in GateDecision.__dataclass_fields__


def test_loss_species_suppressed_real_untouched():
    """UNIFIED core-gate contract (two mechanisms):
      PRIMARY  — a motif-empty LOSS-like cluster (core < floor) is ROUTED TO FALLBACK by the gate
                 (reason 'low_core_fraction'), so it never enters the 2nd-pass balloon; a REAL-like
                 cluster (core ≈ 1) PASSES to modesep ('ok'). (Legacy design let BOTH pass and relied
                 only on the π_species term — that's what changed.)
      SECONDARY — the π_species core term still suppresses loss / spares real on EITHER route, so a
                 loss that slips through (or any fallback-routed loss) is down-shifted there too."""
    rng = np.random.default_rng(1)
    u2 = rng.normal([-25, -8, -4], 3, (500, 3))
    loss = rng.normal([3, -1, -1], 1.5, (80, 3))    # raw_sum ≈ 1 → no core
    real = rng.normal([9, 6, 3], 2, (80, 3))         # raw_sum ≈ 18 → core present

    g_loss = evaluate_gate(loss, u2, mu_u12_5p_raw=3.0, n_eff_candidates=80.0,
                           separation_5p=2.0, valley_depth_fn=compute_valley_depth)
    g_real = evaluate_gate(real, u2, mu_u12_5p_raw=9.0, n_eff_candidates=80.0,
                           separation_5p=2.0, valley_depth_fn=compute_valley_depth)
    # PRIMARY: the gate now SEPARATES them — loss → fallback, real → modesep.
    assert g_loss.core_fraction < DEFAULT_CORE_GATE_FLOOR < g_real.core_fraction, \
        (g_loss.core_fraction, g_real.core_fraction)
    assert (not g_loss.passes) and g_loss.reason == "low_core_fraction", g_loss.reason
    assert g_real.passes and g_real.reason == "ok", g_real.reason

    # SECONDARY: the π_species core term suppresses loss / spares real on EITHER route (route-independent).
    def adj(g, core):
        return compute_adjusted_score(99.0, g.gap_fraction_ucl, 5.0, k_sigma=0.0, core_fraction=core)

    loss_with = adj(g_loss, g_loss.core_fraction)
    loss_without = adj(g_loss, None)
    real_with = adj(g_real, g_real.core_fraction)
    real_without = adj(g_real, None)

    # LOSS: the core term strongly suppresses (≥40-point drop), REAL: barely moved (<2 points)
    assert loss_without - loss_with > 40, (loss_without, loss_with)
    assert loss_with < 50, loss_with
    assert abs(real_without - real_with) < 2, (real_without, real_with)
    assert real_with > 95, real_with
    print(f"  gate:   LOSS→{g_loss.reason} (core {g_loss.core_fraction:.2f});  REAL→{g_real.reason} (core {g_real.core_fraction:.2f})")
    print(f"  π-term: LOSS svm99 {loss_without:.1f}→{loss_with:.1f};  REAL svm99 {real_without:.1f}→{real_with:.1f}")


def test_core_gate_disable_switch():
    """Back-compat escape hatch: core_gate_floor < 0 disables the routing gate — a loss-like cluster
    then PASSES (reason 'ok') exactly as in the legacy design (π_species term still suppresses it)."""
    rng = np.random.default_rng(2)
    u2 = rng.normal([-25, -8, -4], 3, (500, 3))
    loss = rng.normal([3, -1, -1], 1.5, (80, 3))
    g_off = evaluate_gate(loss, u2, mu_u12_5p_raw=3.0, n_eff_candidates=80.0,
                          separation_5p=2.0, valley_depth_fn=compute_valley_depth,
                          core_gate_floor=-1.0)
    assert g_off.passes and g_off.reason == "ok", g_off.reason
    assert g_off.core_fraction < DEFAULT_CORE_GATE_FLOOR, g_off.core_fraction


def test_core_term_inert_when_none():
    """core_fraction=None ⇒ the term is fully inert (back-compat: pre-Step-7 behavior unchanged)."""
    p_none, _ = _confidence_shrunk_pi_species(0.45, core_fraction=None)
    p_high, _ = _confidence_shrunk_pi_species(0.45, core_fraction=0.9)
    # high core ≈ no penalty ≈ the None (inert) case
    assert abs(p_none - p_high) < 0.02, (p_none, p_high)


def test_bar_is_conservative_low_valley():
    """The bar sits at the LOW end of the safe valley (mitigation for human-PWM anchoring), not gold-p5."""
    assert DEFAULT_CORE_RAW_SUM_BAR <= 10.0, DEFAULT_CORE_RAW_SUM_BAR


def test_diverged_real_with_moderate_core_survives():
    """A diverged bearer with only ~45% core (the worst-case pressure-test margin) is still trusted."""
    # core_fraction 0.44 (MucorLusit at bar=8) → w_core should keep π_species near π_train
    p, applies = _confidence_shrunk_pi_species(0.40, core_fraction=0.44)
    assert applies
    a = compute_adjusted_score(99.0, 0.40, 0.0, k_sigma=0.0, core_fraction=0.44)
    assert a > 90, a   # the worst-case diverged real bearer's calls survive


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"ALL {len(fns)} CORE-PRESENCE CHECKS PASSED")
