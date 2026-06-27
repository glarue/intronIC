#!/usr/bin/env python3
"""Build a raw_gated intronIC model bundle (transitional scaffolding for the raw-score architecture).

Trains the SVM ensemble on RAW background-corrected log-odds features (via the threaded SVMTrainer with
base_features=RAW_BASE_FEATURES) on the validated v23 training matrix, and packages it in the inference
loader's v3-style format plus `scoring_mode: "raw_gated"` and the continuous-gate params. The inference
path detects `scoring_mode` and skips z-norm/mode-sep/discount, running scoring/species_gate instead.

Training data: the same v23 corpus raw features used to validate the architecture
(eval_corpus/stage1_trainmatrix.npz: columns [5'_raw, bp_raw, 3'_raw, bp_offset, bp_scan_confidence,
support2_raw] + labels). support2_raw is recomputed as a property from the raw triple.

Usage: python scripts/build_raw_gated_bundle.py [out.pkl]
"""
import sys, os, pickle
from types import SimpleNamespace
import numpy as np

from intronIC.core.intron import IntronScores
from intronIC.classification.trainer import SVMTrainer, RAW_BASE_FEATURES
from intronIC.classification.optimizer import SVMParameters
from intronIC.scoring.species_gate import SpeciesGateParams

MATRIX = "/mnt/data/u12/ipa/conservation_corpus/eval_corpus/stage1_trainmatrix.npz"
OUT = sys.argv[1] if len(sys.argv) > 1 else "/mnt/data/u12/ipa/conservation_corpus/eval_corpus/raw_gated.model.pkl"
EXTRA = ("bp_offset", "bp_scan_confidence", "support2_raw")
N_MODELS = int(os.environ.get("N_MODELS", "9"))

d = np.load(MATRIX, allow_pickle=True)
Xraw, y = d["Xraw"], d["y"]
ok = ~np.isnan(Xraw).any(1)
Xraw, y = Xraw[ok], y[ok]
print(f"training matrix: {Xraw.shape}  u12={int(y.sum())} u2={int((y==0).sum())}", flush=True)

_rng = np.random.RandomState(7)
def make_introns(mask):
    out = []
    idxs = np.where(mask)[0]
    lengths = _rng.randint(80, 1500, len(idxs))   # for the length-stratified U2 subsampling
    for j, i in enumerate(idxs):
        sc = IntronScores(
            five_raw_score=float(Xraw[i, 0]), bp_raw_score=float(Xraw[i, 1]),
            three_raw_score=float(Xraw[i, 2]), bp_offset=float(Xraw[i, 3]),
            bp_scan_confidence=float(Xraw[i, 4]),
        )
        # _subsample_u2 reads .length and .sequences (None -> GC defaults to 0.5)
        out.append(SimpleNamespace(intron_id=f"train_{i}", scores=sc,
                                   length=int(lengths[j]), sequences=None))
    return out

u12_introns = make_introns(y == 1)
u2_introns = make_introns(y == 0)

params = SVMParameters(
    C=200.0, calibration_method="isotonic", saturate_enabled=False, include_max=False,
    include_pairwise_mins=False, penalty="l2", class_weight_multiplier=1.0, loss="squared_hinge",
    kernel="rbf", gamma=0.001, extra_features=EXTRA, base_features=RAW_BASE_FEATURES,
)
trainer = SVMTrainer(
    n_models=N_MODELS, random_state=42, kernel="rbf", max_iter=20000,
    extra_feature_names=list(EXTRA), base_features=RAW_BASE_FEATURES,
)
print(f"training {N_MODELS}-model raw ensemble (C=200, gamma=0.001, rbf, isotonic)...", flush=True)
ensemble = trainer.train_ensemble(u12_introns, u2_introns, params, subsample_u2=True, subsample_ratio=0.8)

m0 = ensemble.models[0]
print("model0 params: base=", m0.parameters.base_features, " extra=", m0.parameters.extra_features, flush=True)
# sanity: the calibrated model predicts on a raw 6-vector
import numpy as _np
_p = m0.model.predict_proba(_np.array([[14.0, 8.0, 1.5, -12.0, 7.0, 8.0]]))[0, 1]
print(f"sanity predict_proba on a strong raw U12 vector: {_p:.3f}", flush=True)

bundle = {
    "version": "raw_gated_v1",
    "model_id": "raw_gated_v23corpus_C200_g0.001_isotonic",
    "ensemble": ensemble,
    "scoring_mode": "raw_gated",
    "input_features": list(RAW_BASE_FEATURES) + list(EXTRA),
    "raw_gate_params": vars(SpeciesGateParams()),   # calibrated defaults
    "threshold": 50.0,
    "provenance": {
        "training_matrix": MATRIX, "n_u12": int((y == 1).sum()), "n_u2": int((y == 0).sum()),
        "n_models": N_MODELS, "note": "transitional raw-score bundle; see eval_corpus/PROPOSED_ARCHITECTURE.md",
    },
}
with open(OUT, "wb") as fh:
    pickle.dump(bundle, fh)
print(f"\nwrote raw_gated bundle -> {OUT}")
print(f"  scoring_mode={bundle['scoring_mode']}  input_features={bundle['input_features']}")
print(f"  gate params: {bundle['raw_gate_params']}")
