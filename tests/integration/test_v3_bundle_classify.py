"""Integration tests: end-to-end classify with the v3 multispecies bundle.

Exercises the v3 load path under various CLI flag combinations to make
sure flags that consumers actually rely on (--threshold, --no-nc,
--no-nc-ss-adjustment, --species-prior, --normalizer-mode) interact
cleanly with the v3 bundle. The unit tests in test_model_io.py cover
the bundle translation; these tests cover what happens after.

Skipped automatically if either:
  - intronIC isn't installed in the dev env, or
  - v3_multispecies_canonical.pkl isn't present (it's not committed to the repo)
"""
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).parent.parent.parent
DATA_DIR = REPO_ROOT / "src" / "intronIC" / "data"
TEST_GENOME = DATA_DIR / "test_data" / "Homo_sapiens.Chr19.Ensembl_91.fa.gz"
TEST_ANNOTATION = DATA_DIR / "test_data" / "Homo_sapiens.Chr19.Ensembl_91.gff3.gz"
V3_BUNDLE = Path(
    "/mnt/data/u12/ipa/training_data/multispecies_v23/v3_multispecies_canonical.pkl"
)

# A chr19 classify run takes ~30-60s on this machine. Each test variant
# is a fresh run, so we cap at 300s per to leave headroom for slower CI.
RUN_TIMEOUT = 300


pytestmark = pytest.mark.skipif(
    not (TEST_GENOME.exists() and TEST_ANNOTATION.exists() and V3_BUNDLE.exists()),
    reason="Chr19 test data or v3_multispecies_canonical.pkl not present",
)


@pytest.fixture(scope="module")
def intronIC_bin():
    dev_bin = Path("/mnt/data/u12/intronIC_v2_env/bin/intronIC")
    if dev_bin.exists():
        return str(dev_bin)
    found = shutil.which("intronIC")
    if found:
        return found
    pytest.skip("intronIC not found on PATH or in dev env")


def _run(intronIC_bin, output_dir, *, extra_args=(), expect_fail=False):
    """Run intronIC classify against chr19 with the v3 bundle, --in-memory."""
    cmd = [
        intronIC_bin, "classify",
        "-g", str(TEST_GENOME),
        "-a", str(TEST_ANNOTATION),
        "-n", "v3_test",
        "-o", str(output_dir),
        "--model", str(V3_BUNDLE),
        "--in-memory",
        "-p", "2",
        *list(extra_args),
    ]
    result = subprocess.run(
        cmd, capture_output=True, text=True, timeout=RUN_TIMEOUT,
    )
    if result.returncode != 0 and not expect_fail:
        pytest.fail(
            f"intronIC classify failed (rc={result.returncode}):\n"
            f"CMD: {' '.join(cmd)}\n"
            f"STDOUT (tail):\n{result.stdout[-2000:]}\n"
            f"STDERR (tail):\n{result.stderr[-2000:]}"
        )
    if result.returncode == 0 and expect_fail:
        pytest.fail(
            f"intronIC classify unexpectedly succeeded:\n"
            f"CMD: {' '.join(cmd)}\nSTDOUT:\n{result.stdout[-2000:]}"
        )
    return result


def _parse_classification_summary(stdout: str) -> dict:
    """Pull the U12/U2 counts out of the rich-formatted summary table."""
    out = {}
    m = re.search(r"U12-type \(total\)\s*│\s*(\d+(?:,\d+)*)", stdout)
    if m:
        out["u12_total"] = int(m.group(1).replace(",", ""))
    m = re.search(r"U12-type \(AT-AC\)\s*│\s*(\d+(?:,\d+)*)", stdout)
    if m:
        out["u12_at_ac"] = int(m.group(1).replace(",", ""))
    m = re.search(r"U2-type\s*│\s*(\d+(?:,\d+)*)", stdout)
    if m:
        out["u2"] = int(m.group(1).replace(",", ""))
    m = re.search(r"valley_depth=([\d.]+)", stdout)
    if m:
        out["valley_depth"] = float(m.group(1))
    return out


# ── Baseline run (shared across multiple assertions) ──────────────────


def _read_run_log(output_dir: Path) -> str:
    logs = list(output_dir.glob("*.iic.log"))
    if not logs:
        return ""
    return logs[0].read_text()


@pytest.fixture(scope="module")
def baseline_run(intronIC_bin):
    """Vanilla v3 classify on chr19 with default args (threshold=95)."""
    with tempfile.TemporaryDirectory(prefix="intronIC_v3_baseline_") as tmpdir:
        out = Path(tmpdir)
        result = _run(intronIC_bin, out)
        summary = _parse_classification_summary(result.stdout)
        run_log = _read_run_log(out)
        # Snapshot the score_info rows for the determinism test
        si_files = list(out.glob("*score_info*"))
        assert si_files, "no score_info output produced"
        with open(si_files[0]) as f:
            header = f.readline()
            rows = sorted(f.readlines())
        yield {
            "stdout": result.stdout,
            "stderr": result.stderr,
            "log": run_log,
            "summary": summary,
            "score_header": header,
            "score_rows": rows,
            "n_score_rows": len(rows),
        }


class TestBaseline:
    """Vanilla v3 classify produces reasonable output."""

    def test_run_logs_v3_model_id(self, baseline_run):
        # The v3 detection log line goes to the run's .iic.log file via
        # messenger.log_only(...) — not to stdout.
        run_log = baseline_run["log"]
        assert "v3 multispecies bundle" in run_log, (
            "v3 detection / model_id log not emitted — "
            "loader may have fallen through to the legacy branch"
        )
        assert "model_id=v3_multispecies_" in run_log

    def test_run_uses_adaptive_normalizer(self, baseline_run):
        # v3 ships normalizer=None → auto mode falls through to adaptive
        out = baseline_run["stdout"]
        assert ("Fitting normalizer on experimental data" in out
                or "adaptive" in out.lower()), (
            "v3 should fall through to adaptive normalizer (saw no log line)")

    def test_classification_summary_present(self, baseline_run):
        s = baseline_run["summary"]
        assert "u12_total" in s
        assert "u2" in s
        # Chr19 has roughly 80-100 U12 calls at threshold=95 with the v3 model.
        # Allow a generous range so this isn't brittle to model retrains.
        assert 30 <= s["u12_total"] <= 200, (
            f"unexpected U12 count on Chr19: {s['u12_total']}"
        )
        # AT-AC should be ~10-30% of U12 calls in human
        if "u12_at_ac" in s:
            assert s["u12_at_ac"] <= s["u12_total"]

    def test_valley_depth_reported(self, baseline_run):
        # Score adjustment uses valley_depth from cluster validation
        assert "valley_depth" in baseline_run["summary"]
        vd = baseline_run["summary"]["valley_depth"]
        assert 0.0 <= vd <= 1.0


# ── Determinism ───────────────────────────────────────────────────────


class TestDeterminism:
    """Two runs with identical args produce bit-identical score output."""

    def test_two_runs_identical(self, intronIC_bin, baseline_run):
        with tempfile.TemporaryDirectory(prefix="intronIC_v3_det_") as tmpdir:
            out2 = Path(tmpdir)
            _run(intronIC_bin, out2)
            si = list(out2.glob("*score_info*"))[0]
            with open(si) as f:
                header2 = f.readline()
                rows2 = sorted(f.readlines())
        assert header2 == baseline_run["score_header"]
        assert len(rows2) == baseline_run["n_score_rows"]
        assert rows2 == baseline_run["score_rows"], (
            "v3 classify is non-deterministic across runs"
        )


# ── Threshold flag ────────────────────────────────────────────────────


class TestThreshold:
    """--threshold overrides the bundle-default threshold (50)."""

    def test_threshold_lower_admits_more(self, intronIC_bin, baseline_run):
        with tempfile.TemporaryDirectory(prefix="intronIC_v3_t50_") as tmpdir:
            out = Path(tmpdir)
            r = _run(intronIC_bin, out, extra_args=("-t", "50"))
            s = _parse_classification_summary(r.stdout)
        # Lower threshold => at least as many U12 calls
        assert s["u12_total"] >= baseline_run["summary"]["u12_total"], (
            f"-t 50 should admit ≥ U12 calls than default (95): "
            f"got {s['u12_total']} vs {baseline_run['summary']['u12_total']}"
        )

    def test_threshold_higher_admits_fewer(self, intronIC_bin, baseline_run):
        with tempfile.TemporaryDirectory(prefix="intronIC_v3_t99_") as tmpdir:
            out = Path(tmpdir)
            r = _run(intronIC_bin, out, extra_args=("-t", "99"))
            s = _parse_classification_summary(r.stdout)
        assert s["u12_total"] <= baseline_run["summary"]["u12_total"], (
            f"-t 99 should admit ≤ U12 calls than default (95): "
            f"got {s['u12_total']} vs {baseline_run['summary']['u12_total']}"
        )


# ── Noncanonical / boundary correction flags ──────────────────────────


class TestNoncanonical:
    """--no-nc and --no-nc-ss-adjustment compose with the v3 bundle."""

    def test_no_nc_completes_successfully(self, intronIC_bin):
        # --no-nc excludes noncanonical introns entirely from scoring
        with tempfile.TemporaryDirectory(prefix="intronIC_v3_nonc_") as tmpdir:
            out = Path(tmpdir)
            r = _run(intronIC_bin, out, extra_args=("--no-nc",))
        s = _parse_classification_summary(r.stdout)
        assert "u12_total" in s
        # With --no-nc we still expect U12 calls (most U12s are canonical)
        assert s["u12_total"] > 0

    def test_no_nc_ss_adjustment_completes_successfully(self, intronIC_bin):
        # --no-nc-ss-adjustment disables the boundary-correction post-pass
        # for noncanonical splice sites (search_u12_boundary). Should still
        # produce a clean classification.
        with tempfile.TemporaryDirectory(prefix="intronIC_v3_nossadj_") as tmpdir:
            out = Path(tmpdir)
            r = _run(intronIC_bin, out, extra_args=("--no-nc-ss-adjustment",))
        s = _parse_classification_summary(r.stdout)
        assert "u12_total" in s
        assert s["u12_total"] > 0


# ── Species prior (uses the new training_prior field) ────────────────


class TestSpeciesPrior:
    """--species-prior triggers the v3 training_prior code path."""

    def test_species_prior_logs_prior_adjustment(self, intronIC_bin):
        # When --species-prior is passed, production runs _apply_prior_adjustment
        # which reads training_prior from the bundle. The v3 loader synthesizes
        # this from training stats; verify the adjustment actually fires.
        with tempfile.TemporaryDirectory(prefix="intronIC_v3_prior_") as tmpdir:
            out = Path(tmpdir)
            r = _run(intronIC_bin, out, extra_args=("--species-prior", "0.005"))
        # The prior-adjustment helper logs an info line when it runs
        assert ("prior adjustment" in r.stdout.lower()
                or "prior=" in r.stdout.lower()), (
            "expected prior-adjustment log when --species-prior is set"
        )


# ── Normalizer mode ───────────────────────────────────────────────────


class TestNormalizerMode:
    """v3 ships normalizer=None — `auto` falls through, `human` errors."""

    def test_auto_mode_succeeds(self, intronIC_bin):
        # Default is --normalizer-mode auto; baseline already covers this,
        # but we make it explicit here for clarity.
        with tempfile.TemporaryDirectory(prefix="intronIC_v3_auto_") as tmpdir:
            out = Path(tmpdir)
            _run(intronIC_bin, out, extra_args=("--normalizer-mode", "auto"))

    def test_adaptive_mode_succeeds(self, intronIC_bin):
        with tempfile.TemporaryDirectory(prefix="intronIC_v3_adaptive_") as tmpdir:
            out = Path(tmpdir)
            _run(intronIC_bin, out, extra_args=("--normalizer-mode", "adaptive"))

    def test_human_mode_errors_cleanly(self, intronIC_bin):
        # v3 has no saved normalizer; `--normalizer-mode human` should
        # produce a clear ValueError, not a cryptic AttributeError.
        with tempfile.TemporaryDirectory(prefix="intronIC_v3_human_") as tmpdir:
            out = Path(tmpdir)
            r = _run(intronIC_bin, out,
                     extra_args=("--normalizer-mode", "human"),
                     expect_fail=True)
        combined = (r.stdout + r.stderr).lower()
        assert "normalizer" in combined, (
            "error from --normalizer-mode human should mention 'normalizer'"
        )
