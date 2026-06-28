"""The post-z-stack bundle guard: only raw-feature bundles are scoreable.

Replaces the removed test_v3_bundle_classify end-to-end test (v3/zscore bundles no
longer score). Asserts assert_scoreable_bundle() passes raw_gated/pmotif_adjudicated
and fails fast on every legacy shape (v3/modesep runtime, zscore, bare ensemble).
"""
import pytest

from intronIC.utils.model_io import (
    assert_scoreable_bundle,
    UnsupportedBundleError,
    SCOREABLE_SCORING_MODES,
)


@pytest.mark.parametrize("mode", ["raw_gated", "pmotif_adjudicated"])
def test_scoreable_bundle_passes(mode):
    # A minimal bundle dict carrying a scoreable scoring_mode is accepted.
    assert_scoreable_bundle({"scoring_mode": mode, "ensemble": object()})


@pytest.mark.parametrize("bundle", [
    {"scoring_mode": "zscore"},          # explicit legacy mode
    {"normalizer_mode": "modesep"},      # v3 modesep runtime dict (no scoring_mode)
    {"version": "v3"},                   # raw v3 bundle
    {},                                  # no scoring_mode at all
    object(),                            # bare (non-dict) SVMEnsemble-style object
])
def test_legacy_bundle_is_rejected(bundle):
    with pytest.raises(UnsupportedBundleError):
        assert_scoreable_bundle(bundle)


def test_error_message_points_to_recovery_path():
    with pytest.raises(UnsupportedBundleError) as exc:
        assert_scoreable_bundle({"scoring_mode": "zscore"})
    msg = str(exc.value)
    assert "pre-zstack-removal" in msg            # the reproducibility tag
    assert "pmotif_adjudicated" in msg            # the retrain instruction


def test_scoreable_modes_constant():
    assert SCOREABLE_SCORING_MODES == ("raw_gated", "pmotif_adjudicated")
