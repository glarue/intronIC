"""
Integration test: classify must accept all three input modes end-to-end.

Three modes are documented in `intronIC classify --help`:
  (1) -g + -a   annotation file (gff3/gtf) + genome — covered by other tests
  (2) -g + -b   BED file with intron coordinates + genome
  (3) -q        pre-extracted intron sequences in .iic format

This file covers (2) and (3). The fixtures are produced at test time by
running `intronIC extract` on the bundled chr19 data, which keeps the
test self-contained without committing binary `.iic` blobs to the repo.
"""

import subprocess
from pathlib import Path

import pytest


DATA_DIR = Path(__file__).parent.parent.parent / "src" / "intronIC" / "data"
TEST_GENOME = DATA_DIR / "test_data" / "Homo_sapiens.Chr19.Ensembl_91.fa.gz"
TEST_ANNOTATION = DATA_DIR / "test_data" / "Homo_sapiens.Chr19.Ensembl_91.gff3.gz"


@pytest.fixture(scope="module")
def intronIC_bin():
    """Return the intronIC binary path, preferring the dev-installed version."""
    import shutil
    dev_bin = Path("/mnt/data/u12/intronIC_v2_env/bin/intronIC")
    if dev_bin.exists():
        return str(dev_bin)
    found = shutil.which("intronIC")
    if found:
        return found
    pytest.skip("intronIC not found on PATH or in dev env")


@pytest.fixture(scope="module")
def extract_outputs(intronIC_bin, tmp_path_factory):
    """Run `intronIC extract` once per module to produce bed.iic + introns.iic.

    Both alternative-input tests consume these artifacts; doing extraction
    once amortizes the ~15s cost across the BED and sequence test classes.
    """
    out_dir = tmp_path_factory.mktemp("extract_chr19")
    cmd = [
        intronIC_bin, "extract",
        "-g", str(TEST_GENOME),
        "-a", str(TEST_ANNOTATION),
        "-n", "ExtChr19",
        "-o", str(out_dir),
        "-p", "2",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if result.returncode != 0:
        pytest.fail(
            "intronIC extract failed:\n"
            f"STDOUT:\n{result.stdout[-1500:]}\n"
            f"STDERR:\n{result.stderr[-1500:]}"
        )

    bed = next(out_dir.glob("*.bed.iic"), None)
    seqs = next(out_dir.glob("*.introns.iic"), None)
    if bed is None or seqs is None:
        pytest.fail(
            f"extract did not produce expected artifacts in {out_dir}: "
            f"{[p.name for p in out_dir.iterdir()]}"
        )
    return {"bed": bed, "seqs": seqs}


def _scored_count(score_info_path: Path) -> int:
    """Count rows in score_info.iic with a non-NA SVM score (column 3)."""
    with open(score_info_path) as f:
        f.readline()  # header
        return sum(
            1 for line in f
            if line.split("\t", 3)[2] not in ("NA", "")
        )


@pytest.mark.skipif(
    not TEST_GENOME.exists() or not TEST_ANNOTATION.exists(),
    reason="Chr19 test data not found",
)
class TestBedInput:
    """`classify -g <fa> -b <bed>` must complete and emit scored introns.

    Regression guard: the BED-input path bypasses the annotation parser
    and feeds coordinates directly to extraction. It exercises a different
    code path from the GFF/GTF input that the streaming-equivalence tests
    cover.
    """

    @pytest.fixture(scope="class")
    def run_results(self, intronIC_bin, extract_outputs, tmp_path_factory):
        out_dir = tmp_path_factory.mktemp("classify_bed")
        cmd = [
            intronIC_bin, "classify",
            "-g", str(TEST_GENOME),
            "-b", str(extract_outputs["bed"]),
            "-n", "BedChr19",
            "-o", str(out_dir),
            "-p", "2",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            pytest.fail(
                "intronIC classify -b failed:\n"
                f"STDOUT:\n{result.stdout[-1500:]}\n"
                f"STDERR:\n{result.stderr[-1500:]}"
            )
        return out_dir

    def test_score_info_has_scored_rows(self, run_results):
        score_path = next(Path(run_results).glob("*score_info*"))
        n = _scored_count(score_path)
        assert n > 1000, f"BED-input classify produced only {n} scored introns"

    def test_meta_file_exists(self, run_results):
        assert any(Path(run_results).glob("*meta*")), "no meta.iic written"

    def test_bed_file_exists(self, run_results):
        assert any(Path(run_results).glob("*.bed.iic")), "no bed.iic written"


@pytest.mark.skipif(
    not TEST_GENOME.exists() or not TEST_ANNOTATION.exists(),
    reason="Chr19 test data not found",
)
class TestSequenceInput:
    """`classify -q <introns.iic>` must complete and emit scored introns.

    Regression guard: the sequence-input path skips genome extraction and
    consumes pre-extracted sequences directly. It is the only mode that
    works without a genome FASTA.
    """

    @pytest.fixture(scope="class")
    def run_results(self, intronIC_bin, extract_outputs, tmp_path_factory):
        out_dir = tmp_path_factory.mktemp("classify_seq")
        cmd = [
            intronIC_bin, "classify",
            "-q", str(extract_outputs["seqs"]),
            "-n", "SeqChr19",
            "-o", str(out_dir),
            "-p", "2",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            pytest.fail(
                "intronIC classify -q failed:\n"
                f"STDOUT:\n{result.stdout[-1500:]}\n"
                f"STDERR:\n{result.stderr[-1500:]}"
            )
        return out_dir

    def test_score_info_has_scored_rows(self, run_results):
        score_path = next(Path(run_results).glob("*score_info*"))
        n = _scored_count(score_path)
        assert n > 1000, f"sequence-input classify produced only {n} scored introns"

    def test_meta_file_exists(self, run_results):
        assert any(Path(run_results).glob("*meta*")), "no meta.iic written"
