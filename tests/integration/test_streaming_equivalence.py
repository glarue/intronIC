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


def _normalize_score_row(row: str, header: str) -> tuple:
    """Drop the species-name prefix from the label column, keep everything else.

    The species-name prefix differs between runs (e.g. ``HomSt-`` vs
    ``HomIm-``) but the rest of the row should be identical between
    classify modes on the same input.
    """
    fields = row.rstrip("\n").split("\t")
    if not fields:
        return ()
    label = fields[0]
    if "-" in label:
        label = label.split("-", 1)[1]
    fields[0] = label
    return tuple(fields)


@pytest.mark.skipif(
    not TEST_GENOME.exists() or not TEST_ANNOTATION.exists(),
    reason="Chr19 test data not found",
)
class TestStreamingMatchesInMemory:
    """Streaming and in-memory must produce equivalent classifications.

    The two paths ultimately call the same scoring + classification
    components on the same intron set; the only legitimate differences
    in output are species-name prefixes and ordering. This test asserts
    every scored intron has the same SVM/adjusted/relative scores in
    both modes, and that they classify the same locus the same way.

    Regression guard: in v2.4-development the streaming-classify worker
    inherited a pre-existing bug where coordinate-duplicate introns were
    scored (one entry per isoform) while in-memory's prefilter dropped
    them, producing ~5,900 extra rows and ~21 inflated U12 calls on
    full human. The v2.4 fix unifies the duplicate handling so this
    test passes.
    """

    @pytest.fixture(scope="class")
    def run_results(self, intronIC_bin):
        with tempfile.TemporaryDirectory(prefix="intronIC_mode_equiv_") as tmpdir:
            dir_stream = Path(tmpdir) / "streaming"
            dir_inmem = Path(tmpdir) / "inmemory"
            dir_stream.mkdir()
            dir_inmem.mkdir()

            _run_classify(
                intronIC_bin, dir_stream, processes=2,
                species_name="ModeStream", mode="streaming",
            )
            _run_classify(
                intronIC_bin, dir_inmem, processes=2,
                species_name="ModeInmem", mode="in-memory",
            )

            sh, srows = _read_score_info(dir_stream, "ModeStream")
            ih, irows = _read_score_info(dir_inmem, "ModeInmem")

            yield {
                "stream_header": sh,
                "inmem_header": ih,
                "stream_rows": srows,
                "inmem_rows": irows,
            }

    def test_score_info_headers_match(self, run_results):
        assert run_results["stream_header"] == run_results["inmem_header"]

    def test_scored_intron_set_matches(self, run_results):
        """Every locus scored in one mode must be scored in the other.

        Compares the set of intron labels (sans species prefix) that have
        a non-NA SVM score in each output. Symmetric difference must be empty.
        """
        def _scored_labels(rows):
            out = set()
            for row in rows:
                fields = row.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                if fields[2] in ("NA", ""):
                    continue
                label = fields[0]
                if "-" in label:
                    label = label.split("-", 1)[1]
                out.add(label)
            return out

        stream_set = _scored_labels(run_results["stream_rows"])
        inmem_set = _scored_labels(run_results["inmem_rows"])

        only_stream = stream_set - inmem_set
        only_inmem = inmem_set - stream_set

        if only_stream or only_inmem:
            pytest.fail(
                f"Scored-intron set differs between modes:\n"
                f"  streaming-only: {len(only_stream)} (sample: "
                f"{sorted(only_stream)[:3]})\n"
                f"  in-memory-only: {len(only_inmem)} (sample: "
                f"{sorted(only_inmem)[:3]})"
            )

    def test_per_intron_scores_match(self, run_results):
        """Every scored intron has identical SVM/adjusted/relative scores in both modes."""
        def _by_label(rows):
            out = {}
            for row in rows:
                fields = row.rstrip("\n").split("\t")
                if len(fields) < 3 or fields[2] in ("NA", ""):
                    continue
                label = fields[0]
                if "-" in label:
                    label = label.split("-", 1)[1]
                out[label] = fields
            return out

        stream = _by_label(run_results["stream_rows"])
        inmem = _by_label(run_results["inmem_rows"])

        common = sorted(set(stream) & set(inmem))
        assert common, "no intron labels in common between modes"

        # Compare columns: rel_score (1), svm_score (2), adjusted_score (-2), ensemble_sigma (-1)
        mismatches = []
        for label in common:
            s_row, i_row = stream[label], inmem[label]
            for idx, name in [(1, "rel_score"), (2, "svm_score"),
                              (-2, "adjusted_score"), (-1, "ensemble_sigma")]:
                if s_row[idx] != i_row[idx]:
                    mismatches.append((label, name, s_row[idx], i_row[idx]))
                    break  # one mismatch per intron is enough
        if mismatches:
            sample = mismatches[:5]
            msg = (
                f"{len(mismatches)} introns have differing scores between modes:\n"
                + "\n".join(
                    f"  [{lab}] {name}: streaming={s} vs in-memory={i}"
                    for lab, name, s, i in sample
                )
            )
            pytest.fail(msg)
