"""Tests for `intronIC.utils.model_io.normalize_model_bundle`.

After supplant 2b the loader is a pass-through: the legacy v3/modesep translation
(z-feature ensemble build + fallback normalizer) was removed with the z stack, so
every bundle passes through unchanged and a non-raw bundle is rejected later by the
load-time guard (`assert_scoreable_bundle`, covered in test_bundle_guard.py).
"""
from unittest.mock import MagicMock

from intronIC.utils.model_io import normalize_model_bundle


class TestNormalizeModelBundle:
    """`normalize_model_bundle` is now a pass-through for every input shape."""

    def test_bare_ensemble_passes_through(self):
        sentinel = MagicMock(name="bare_SVMEnsemble")
        assert normalize_model_bundle(sentinel) is sentinel

    def test_raw_bundle_dict_passes_through(self):
        raw_bundle = {
            "ensemble": MagicMock(name="ensemble"),
            "scoring_mode": "pmotif_adjudicated",
            "adjudicator_params": {},
        }
        assert normalize_model_bundle(raw_bundle) is raw_bundle

    def test_legacy_v3_dict_passes_through_unchanged(self):
        # A legacy v3 bundle is NO LONGER translated — it passes through as-is and is
        # rejected at the load-time guard (no scoring_mode). The loader must not mutate it.
        v3_bundle = {"version": "v3", "model_id": "legacy", "config": {}, "seeds": {}}
        out = normalize_model_bundle(v3_bundle)
        assert out is v3_bundle
        assert "ensemble" not in out  # no runtime translation happened
