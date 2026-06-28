"""
Integration test: the pmotif_adjudicated scoring mode is streaming/in-memory bit-identical.

The pmotif adjudicator runs in the SHARED _run_post_classification_pipeline, so streaming and in-memory
must produce identical score_info/meta/bed. This test builds a tiny (1-model) pmotif_adjudicated bundle on
the fly via `intronIC train --scoring-mode pmotif_adjudicated` (no committed binary, no dev-panel
dependency — the bundle is intentionally uncalibrated, which is fine because PARITY is independent of
calibration), then classifies the bundled Chr19 human data in both modes and asserts equality + that the
new interpretable columns (P_motif, q, P_adj, P_adj_lo, P_adj_hi) are present.

Runs in ~2-3 minutes (one tiny train + two classify passes). Skips gracefully if the dev env / test data
are unavailable.
"""
import subprocess
import tempfile
from pathlib import Path

import pytest


DATA_DIR = Path(__file__).parent.parent.parent / "src" / "intronIC" / "data"
TEST_GENOME = DATA_DIR / "test_data" / "Homo_sapiens.Chr19.Ensembl_91.fa.gz"
TEST_ANNOTATION = DATA_DIR / "test_data" / "Homo_sapiens.Chr19.Ensembl_91.gff3.gz"

PMOTIF_COLUMNS = ("P_motif", "q", "P_adj", "P_adj_lo", "P_adj_hi")


@pytest.fixture(scope="module")
def intronIC_bin():
    import shutil
    dev_bin = Path("/mnt/data/u12/intronIC_v2_env/bin/intronIC")
    if dev_bin.exists():
        return str(dev_bin)
    found = shutil.which("intronIC")
    if found:
        return found
    pytest.skip("intronIC not found on PATH or in dev env")


@pytest.fixture(scope="module")
def pmotif_bundle(intronIC_bin, tmp_path_factory):
    """Train a tiny 1-model pmotif_adjudicated bundle from the built-in reference (uncalibrated; parity
    does not depend on calibration). Skips if training is unavailable in this environment."""
    out = tmp_path_factory.mktemp("pmotif_train")
    cmd = [
        intronIC_bin, "train",
        "--scoring-mode", "pmotif_adjudicated",
        "--n-models", "1",
        "-n", "pmotif_fixture",
        "-o", str(out),
        "-p", "2",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if result.returncode != 0:
        pytest.skip(
            "could not build a pmotif_adjudicated fixture bundle (built-in reference unavailable?):\n"
            f"{result.stderr[-1500:]}"
        )
    bundles = list(out.glob("*.model.pkl"))
    if not bundles:
        pytest.skip("pmotif train produced no .model.pkl")
    return str(bundles[0])


def _run_classify(intronIC_bin, bundle, output_dir, mode):
    cmd = [
        intronIC_bin, "classify",
        "-g", str(TEST_GENOME),
        "-a", str(TEST_ANNOTATION),
        "-n", "pmparity",
        "--model", bundle,
        "-o", str(output_dir),
        "-p", "2",
        f"--{mode}",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if result.returncode != 0:
        pytest.fail(
            f"intronIC classify (pmotif, {mode}) failed:\n"
            f"STDOUT:\n{result.stdout[-1500:]}\nSTDERR:\n{result.stderr[-1500:]}"
        )
    return output_dir


def _read_sorted(path):
    with open(path) as f:
        header = f.readline()
        rows = sorted(f.readlines())
    return header, rows


def _find(output_dir, glob_pat):
    matches = list(Path(output_dir).glob(glob_pat))
    assert matches, f"no {glob_pat} in {list(Path(output_dir).iterdir())}"
    return matches[0]


@pytest.mark.skipif(
    not (TEST_GENOME.exists() and TEST_ANNOTATION.exists()),
    reason="Chr19 test data not found",
)
class TestPmotifAdjudicatedParity:
    """Streaming and in-memory pmotif_adjudicated output must be bit-identical."""

    @pytest.fixture(scope="class")
    def run_results(self, intronIC_bin, pmotif_bundle):
        with tempfile.TemporaryDirectory(prefix="intronIC_pmotif_parity_") as tmpdir:
            dir_im = Path(tmpdir) / "inmem"
            dir_st = Path(tmpdir) / "stream"
            dir_im.mkdir()
            dir_st.mkdir()
            _run_classify(intronIC_bin, pmotif_bundle, dir_im, "in-memory")
            _run_classify(intronIC_bin, pmotif_bundle, dir_st, "streaming")
            yield {
                "im_score": _read_sorted(_find(dir_im, "*score_info*.iic")),
                "st_score": _read_sorted(_find(dir_st, "*score_info*.iic")),
                "im_meta": _read_sorted(_find(dir_im, "*meta*.iic")),
                "st_meta": _read_sorted(_find(dir_st, "*meta*.iic")),
            }

    def test_pmotif_columns_present(self, run_results):
        header = run_results["im_score"][0].rstrip("\n").split("\t")
        for col in PMOTIF_COLUMNS:
            assert col in header, f"{col} missing from score_info header"

    def test_score_info_identical(self, run_results):
        im_header, im_rows = run_results["im_score"]
        st_header, st_rows = run_results["st_score"]
        assert im_header == st_header
        assert im_rows == st_rows, "streaming vs in-memory score_info rows differ"

    def test_meta_identical(self, run_results):
        """meta.iic must also match (the post-process syncs type_id/rel_score into it on both paths)."""
        im_header, im_rows = run_results["im_meta"]
        st_header, st_rows = run_results["st_meta"]
        assert im_header == st_header
        assert im_rows == st_rows, "streaming vs in-memory meta rows differ"
