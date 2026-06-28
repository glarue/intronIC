#!/usr/bin/env python3
"""Build a ``pmotif_adjudicated`` model bundle from an existing raw-feature SVM ensemble.

The ``pmotif_adjudicated`` scoring mode (docs/raw_gated_scoring.md §0b/§0c) computes a per-intron motif
probability ``P_motif`` from the ensemble decision_function MARGIN (Platt-calibrated) and a per-species
bearer probability ``q`` from the depth_tail adjudicator, writing ``P_adj = q * P_motif``. It reuses the same
raw-feature ensemble as the ``raw_gated`` mode; the only difference is the post-classification layer and the
stamped scoring_mode + adjudicator params.

This is a re-stamp utility: it takes an existing raw-ensemble bundle (e.g. the one from
``build_raw_gated_bundle.py``) and writes a pmotif_adjudicated bundle. (Production / step-7 will instead
emit this directly from ``main_train``; this script is the transitional builder for end-to-end testing.)

Usage:
    build_pmotif_adjudicated_bundle.py <input_raw_bundle.pkl> <output_pmotif_bundle.pkl>
"""
import sys
import pickle
from dataclasses import asdict

from intronIC.scoring.species_adjudicator import AdjudicatorParams, ADJUDICATOR_PARAMS_VERSION

if len(sys.argv) != 3:
    sys.exit(__doc__)
IN, OUT = sys.argv[1], sys.argv[2]

src = pickle.load(open(IN, "rb"))
if not isinstance(src, dict) or "ensemble" not in src:
    sys.exit(f"input bundle {IN} has no 'ensemble' key (got keys: "
             f"{sorted(src.keys()) if isinstance(src, dict) else type(src)})")

adj_params = asdict(AdjudicatorParams())   # calibrated defaults (Platt + q + EVT settings), version-pinned

bundle = {
    "version": "pmotif_adjudicated_v1",
    "model_id": "pmotif_adjudicated_v23corpus_C200_g0.001",
    "ensemble": src["ensemble"],
    "scoring_mode": "pmotif_adjudicated",
    "input_features": src.get("input_features"),
    "adjudicator_params": adj_params,
    "threshold": 50.0,
    "provenance": {
        "rebuilt_from": IN,
        "src_version": src.get("version"),
        "src_model_id": src.get("model_id"),
        "adjudicator_params_version": ADJUDICATOR_PARAMS_VERSION,
        "note": "re-stamp of a raw-feature ensemble for pmotif_adjudicated scoring; "
                "see docs/raw_gated_scoring.md §0c step 1",
    },
}
with open(OUT, "wb") as fh:
    pickle.dump(bundle, fh)

print(f"wrote pmotif_adjudicated bundle -> {OUT}")
print(f"  ensemble: {len(bundle['ensemble'].models)} models   scoring_mode={bundle['scoring_mode']}")
print(f"  adjudicator_params_version={ADJUDICATOR_PARAMS_VERSION}")
print(f"  q = sigmoid({adj_params['q_a']}*depth_tail + {adj_params['q_b']})  "
      f"P_motif = sigmoid({adj_params['platt_a']}*margin + {adj_params['platt_c']})")
