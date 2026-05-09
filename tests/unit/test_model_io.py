"""Tests for `intronIC.utils.model_io.normalize_model_bundle`.

Covers all three model formats the loader has to accept:
- bare SVMEnsemble (legacy)
- v2.3 dict (current production training save format)
- v3 dict (v3 multispecies bundle)
"""
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from intronIC.utils.model_io import (
    V3_VERSION,
    _build_v3_ensemble,
    normalize_model_bundle,
)


def _make_v3_bundle(*, n_seeds=2, n_models_per_seed=3) -> dict:
    """Build a minimal v3 bundle whose sub-models are MagicMock stubs.

    Sufficient for testing the loader's wrapping / metadata translation;
    not sufficient for testing predict_proba (use the real bundle for that).
    """
    seeds = list(range(42, 42 + n_seeds))
    return {
        "version": V3_VERSION,
        "model_id": "test_v3",
        "config": {
            "C": 200.0,
            "gamma": 0.001,
            "kernel": "rbf",
            "calibration_method": "isotonic",
            "easy_fraction": 0.75,
            "ensemble_size": n_models_per_seed,
            "seeds": seeds,
            "input_features": [
                "five_z_score",
                "bp_z_score",
                "three_z_score",
                "bp_offset",
                "bp_scan_confidence",
                "support2",
            ],
        },
        "training": {
            "corpus_rows": 10_000,
            "u12_positives": 1_000,
            "test_positives": 100,
            "n_train": 9_000,
            "n_test": 1_000,
        },
        "seeds": {
            seed: {"models": [MagicMock(name=f"sub_s{seed}_m{i}")
                              for i in range(n_models_per_seed)]}
            for seed in seeds
        },
    }


class TestNormalizeModelBundle:
    """`normalize_model_bundle` should be format-agnostic."""

    def test_bare_ensemble_passes_through(self):
        sentinel = MagicMock(name="bare_SVMEnsemble")
        out = normalize_model_bundle(sentinel)
        assert out is sentinel

    def test_v23_dict_passes_through(self):
        v23_bundle = {
            "ensemble": MagicMock(name="ensemble"),
            "normalizer": MagicMock(name="normalizer"),
            "threshold": 90.0,
            "training_prior": 0.005,
            "human_negative_stats": None,
        }
        out = normalize_model_bundle(v23_bundle)
        assert out is v23_bundle  # same object — no translation needed

    def test_v3_dict_translates_to_runtime_shape(self):
        bundle = _make_v3_bundle()
        out = normalize_model_bundle(bundle)
        # Runtime shape = legacy v2.3 keys
        for key in ("ensemble", "normalizer", "threshold", "training_prior"):
            assert key in out, f"missing key: {key}"
        # Original bundle preserved for downstream introspection
        assert out.get("_v3_bundle") is bundle

    def test_v3_normalizer_is_none(self):
        out = normalize_model_bundle(_make_v3_bundle())
        assert out["normalizer"] is None  # forces adaptive at runtime

    def test_v3_threshold_default_is_50(self):
        out = normalize_model_bundle(_make_v3_bundle())
        assert out["threshold"] == 50.0


class TestV3EnsembleConstruction:
    """`_build_v3_ensemble` shape and parameter wiring."""

    def test_ensemble_size_equals_seeds_times_models(self):
        bundle = _make_v3_bundle(n_seeds=3, n_models_per_seed=4)
        ensemble, _ = _build_v3_ensemble(bundle)
        assert len(ensemble.models) == 12

    def test_extra_features_excludes_three_base_zscores(self):
        bundle = _make_v3_bundle()
        _, extras = _build_v3_ensemble(bundle)
        assert extras == ("bp_offset", "bp_scan_confidence", "support2")

    def test_each_submodel_has_correct_parameters(self):
        bundle = _make_v3_bundle()
        ensemble, _ = _build_v3_ensemble(bundle)
        params = ensemble.models[0].parameters
        assert params.kernel == "rbf"
        assert params.C == 200.0
        assert params.gamma == 0.001
        assert params.extra_features == ("bp_offset", "bp_scan_confidence", "support2")
        assert params.include_max is False  # never used in new arch
        assert ensemble.models[0].dropped_feature is None  # no feature dropout in v3 multispecies

    def test_rejects_too_few_input_features(self):
        bundle = _make_v3_bundle()
        bundle["config"]["input_features"] = ["only_one"]
        with pytest.raises(ValueError, match="input_features"):
            _build_v3_ensemble(bundle)


class TestV3TrainingPrior:
    """training_prior = train-side U12 fraction, excluding held-out positives."""

    def test_uses_train_only_count(self):
        # u12_positives=1000 corpus-wide, 100 held out, n_train=9000
        # → training_prior = 900 / 9000 = 0.1
        out = normalize_model_bundle(_make_v3_bundle())
        assert out["training_prior"] == pytest.approx(0.1)

    def test_falls_back_to_half_when_invalid(self):
        bundle = _make_v3_bundle()
        bundle["training"] = {}  # no n_train, no u12 counts
        out = normalize_model_bundle(bundle)
        assert out["training_prior"] == 0.5


_FI_V3_BUNDLE = Path(
    "/mnt/data/u12/ipa/training_data/multispecies_v23/v3_multispecies_canonical.pkl"
)


@pytest.mark.skipif(
    not _FI_V3_BUNDLE.exists(),
    reason="v3 multispecies canonical bundle not present in this environment",
)
class TestRealFIv3Bundle:
    """Integration: the real v3_multispecies_canonical.pkl loads and is well-formed."""

    def test_loads_with_expected_size_and_metadata(self):
        from intronIC.utils.model_io import load_model
        runtime = normalize_model_bundle(load_model(_FI_V3_BUNDLE))
        ensemble = runtime["ensemble"]
        assert len(ensemble.models) == 126  # 3 seeds × 42
        v3 = runtime["_v3_bundle"]
        assert v3["model_id"].startswith("v3_multispecies_")
        assert v3["config"]["seeds"] == [42, 100, 200]

    def test_predict_proba_path_works(self):
        """Wrapped sub-models are still callable for prediction."""
        import numpy as np
        from intronIC.utils.model_io import load_model
        runtime = normalize_model_bundle(load_model(_FI_V3_BUNDLE))
        # 6D feature vector matching the bundle's input_features
        X = np.array([
            [3.0, 3.0, 3.0, -25.0, 0.7, 3.0],   # strong U12-like
            [-1.0, -1.0, -1.0, -25.0, 0.5, 0.0],  # clear non-U12
        ])
        sub = runtime["ensemble"].models[0].model
        probs = sub.predict_proba(X)[:, 1]
        assert probs.shape == (2,)
        assert probs[0] > probs[1]  # strong > weak signal
