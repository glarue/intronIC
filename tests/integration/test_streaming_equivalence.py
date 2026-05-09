"""
Integration test: verify -p 1 and -p 3 produce identical classifications.

Exercises both --in-memory and --streaming paths against the bundled
Chr19 human test data. Streaming-with-v3 fits a fresh adaptive
RobustScaler on a per-contig pre-pass before classify; this test also
verifies that pre-pass produces identical results regardless of
parallelism.

Runs in ~3-5 minutes depending on hardware (two full pipelines per mode).
"""

import subprocess
import sys
import tempfile
from pathlib import Path

import pytest


DATA_DIR = Path(__file__).parent.parent.parent / "src" / "intronIC" / "data"
TEST_GENOME = DATA_DIR / "test_data" / "Homo_sapiens.Chr19.Ensembl_91.fa.gz"
TEST_ANNOTATION = DATA_DIR / "test_data" / "Homo_sapiens.Chr19.Ensembl_91.gff3.gz"


@pytest.fixture(scope="module")
def intronIC_bin():
    """Return the intronIC binary path, preferring the dev-installed version."""
    import shutil
    # Check the venv that has intronIC installed in dev mode
    dev_bin = Path("/mnt/data/u12/intronIC_v2_env/bin/intronIC")
    if dev_bin.exists():
        return str(dev_bin)
    # Fall back to PATH
    found = shutil.which("intronIC")
    if found:
        return found
    pytest.skip("intronIC not found on PATH or in dev env")


def _run_classify(intronIC_bin, output_dir, processes, species_name, mode="in-memory"):
    """Run intronIC classify in either --in-memory or --streaming mode."""
    cmd = [
        intronIC_bin, "classify",
        "-g", str(TEST_GENOME),
        "-a", str(TEST_ANNOTATION),
        "-n", species_name,
        "-o", str(output_dir),
        "-p", str(processes),
        f"--{mode}",
    ]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=600,
    )
    if result.returncode != 0:
        pytest.fail(
            f"intronIC classify ({mode}, p={processes}) failed:\n"
            f"STDOUT:\n{result.stdout[-2000:]}\n"
            f"STDERR:\n{result.stderr[-2000:]}"
        )
    return output_dir


def _read_score_info(output_dir, species_name):
    """Read score_info.iic and return sorted rows (excluding header)."""
    pattern = f"{species_name}.full.score_info.iic"
    matches = list(Path(output_dir).glob(f"*score_info*"))
    if not matches:
        pytest.fail(f"No score_info file found in {output_dir}: {list(Path(output_dir).iterdir())}")
    path = matches[0]
    with open(path) as f:
        header = f.readline()
        rows = sorted(f.readlines())
    return header, rows


def _read_meta(output_dir):
    """Read meta.iic and return sorted rows (excluding header)."""
    matches = list(Path(output_dir).glob("*meta*"))
    if not matches:
        pytest.fail(f"No meta file found in {output_dir}")
    path = matches[0]
    with open(path) as f:
        header = f.readline()
        rows = sorted(f.readlines())
    return header, rows


@pytest.mark.skipif(
    not TEST_GENOME.exists() or not TEST_ANNOTATION.exists(),
    reason="Chr19 test data not found",
)
class TestStreamingEquivalence:
    """Verify -p 1 and -p 3 produce identical results."""

    @pytest.fixture(scope="class")
    def run_results(self, intronIC_bin):
        """Run intronIC with -p 1 and -p 3, return output directories."""
        with tempfile.TemporaryDirectory(prefix="intronIC_equiv_") as tmpdir:
            dir_p1 = Path(tmpdir) / "p1"
            dir_p3 = Path(tmpdir) / "p3"
            dir_p1.mkdir()
            dir_p3.mkdir()

            _run_classify(intronIC_bin, dir_p1, processes=1, species_name="equiv_test")
            _run_classify(intronIC_bin, dir_p3, processes=3, species_name="equiv_test")

            # Read results before temp dir is cleaned up
            h1, score_rows_p1 = _read_score_info(dir_p1, "equiv_test")
            h3, score_rows_p3 = _read_score_info(dir_p3, "equiv_test")
            mh1, meta_rows_p1 = _read_meta(dir_p1)
            mh3, meta_rows_p3 = _read_meta(dir_p3)

            yield {
                "score_header_p1": h1,
                "score_header_p3": h3,
                "score_rows_p1": score_rows_p1,
                "score_rows_p3": score_rows_p3,
                "meta_header_p1": mh1,
                "meta_header_p3": mh3,
                "meta_rows_p1": meta_rows_p1,
                "meta_rows_p3": meta_rows_p3,
            }

    def test_score_info_headers_match(self, run_results):
        assert run_results["score_header_p1"] == run_results["score_header_p3"]

    def test_score_info_row_count_match(self, run_results):
        n1 = len(run_results["score_rows_p1"])
        n3 = len(run_results["score_rows_p3"])
        assert n1 == n3, f"Row count mismatch: p1={n1}, p3={n3}"
        assert n1 > 0, "No introns classified"

    def test_score_info_rows_identical(self, run_results):
        rows_p1 = run_results["score_rows_p1"]
        rows_p3 = run_results["score_rows_p3"]
        mismatches = []
        for i, (r1, r3) in enumerate(zip(rows_p1, rows_p3)):
            if r1 != r3:
                mismatches.append((i, r1.strip()[:120], r3.strip()[:120]))
        if mismatches:
            sample = mismatches[:5]
            msg = f"{len(mismatches)} mismatched rows (of {len(rows_p1)}):\n"
            for idx, r1, r3 in sample:
                msg += f"  [{idx}] p1: {r1}\n       p3: {r3}\n"
            pytest.fail(msg)

    def test_meta_row_count_match(self, run_results):
        n1 = len(run_results["meta_rows_p1"])
        n3 = len(run_results["meta_rows_p3"])
        assert n1 == n3, f"Meta row count mismatch: p1={n1}, p3={n3}"

    def test_meta_rows_identical(self, run_results):
        rows_p1 = run_results["meta_rows_p1"]
        rows_p3 = run_results["meta_rows_p3"]
        mismatches = []
        for i, (r1, r3) in enumerate(zip(rows_p1, rows_p3)):
            if r1 != r3:
                mismatches.append((i, r1.strip()[:120], r3.strip()[:120]))
        if mismatches:
            sample = mismatches[:5]
            msg = f"{len(mismatches)} mismatched meta rows:\n"
            for idx, r1, r3 in sample:
                msg += f"  [{idx}] p1: {r1}\n       p3: {r3}\n"
            pytest.fail(msg)


@pytest.mark.skipif(
    not TEST_GENOME.exists() or not TEST_ANNOTATION.exists(),
    reason="Chr19 test data not found",
)
class TestStreamingEquivalenceWithV3:
    """Verify streaming + v3 default: -p 1 and -p 3 produce identical results.

    Exercises the adaptive normalizer pre-pass (added in v2.4) end-to-end
    in the per-contig streaming-classify path. Identical output across
    parallelism levels confirms the fitting pass and the classify pass
    are both deterministic regardless of contig dispatch order.
    """

    @pytest.fixture(scope="class")
    def run_results(self, intronIC_bin):
        with tempfile.TemporaryDirectory(prefix="intronIC_stream_equiv_") as tmpdir:
            dir_p1 = Path(tmpdir) / "p1"
            dir_p3 = Path(tmpdir) / "p3"
            dir_p1.mkdir()
            dir_p3.mkdir()

            _run_classify(
                intronIC_bin, dir_p1, processes=1,
                species_name="stream_test", mode="streaming",
            )
            _run_classify(
                intronIC_bin, dir_p3, processes=3,
                species_name="stream_test", mode="streaming",
            )

            h1, score_rows_p1 = _read_score_info(dir_p1, "stream_test")
            h3, score_rows_p3 = _read_score_info(dir_p3, "stream_test")
            mh1, meta_rows_p1 = _read_meta(dir_p1)
            mh3, meta_rows_p3 = _read_meta(dir_p3)

            yield {
                "score_header_p1": h1,
                "score_header_p3": h3,
                "score_rows_p1": score_rows_p1,
                "score_rows_p3": score_rows_p3,
                "meta_header_p1": mh1,
                "meta_header_p3": mh3,
                "meta_rows_p1": meta_rows_p1,
                "meta_rows_p3": meta_rows_p3,
            }

    def test_score_info_headers_match(self, run_results):
        assert run_results["score_header_p1"] == run_results["score_header_p3"]

    def test_score_info_row_count_match(self, run_results):
        n1 = len(run_results["score_rows_p1"])
        n3 = len(run_results["score_rows_p3"])
        assert n1 == n3, f"streaming row count mismatch: p1={n1}, p3={n3}"
        assert n1 > 0, "streaming produced no introns"

    def test_score_info_rows_identical(self, run_results):
        rows_p1 = run_results["score_rows_p1"]
        rows_p3 = run_results["score_rows_p3"]
        for i, (r1, r3) in enumerate(zip(rows_p1, rows_p3)):
            if r1 != r3:
                pytest.fail(
                    f"streaming p1 vs p3 differ at row {i}:\n"
                    f"  p1: {r1.strip()[:160]}\n"
                    f"  p3: {r3.strip()[:160]}"
                )

    def test_meta_rows_identical(self, run_results):
        rows_p1 = run_results["meta_rows_p1"]
        rows_p3 = run_results["meta_rows_p3"]
        for i, (r1, r3) in enumerate(zip(rows_p1, rows_p3)):
            if r1 != r3:
                pytest.fail(
                    f"streaming meta p1 vs p3 differ at row {i}:\n"
                    f"  p1: {r1.strip()[:160]}\n"
                    f"  p3: {r3.strip()[:160]}"
                )
