"""
Integration tests for the -i (--allow-multiple-isoforms) and
-d (--include-duplicates) CLI flags.

Regression guard for the v2.4.3 fix: prior to v2.4.3 both flags were
silently ignored on multiple code paths:

- In-memory ``classify`` hardcoded ``longest_only=True`` and
  ``include_duplicates=False`` when constructing the post-extraction
  ``IntronFilter`` (cli/main.py ~6198-6199), so neither flag had any
  effect on the in-memory output.
- The streaming per-contig worker hardcoded ``longest_only=True``
  (cli/main.py ~3556) and did not have ``allow_multiple_isoforms``
  plumbed into its worker config_dict, so ``-i`` was silently ignored
  by the streaming pipeline as well.
- ``should_extract_sequences_for`` (extraction/filters.py) skipped all
  coord-duplicate introns from sequence extraction regardless of
  ``include_duplicates``, so even with the filter fix the in-memory
  ``-d`` path could not emit duplicates (they reached the writer with
  no sequences attached).

The fixture exercises every relevant flag combination against both
``--in-memory`` and ``--streaming``, and asserts:

  1. Per-flag behavior (counts + presence of the specific test introns).
  2. Cross-mode parity (streaming bed.iic == in-memory bed.iic for
     every flag combo).

See ``tests/data/isoforms/`` for the synthetic FASTA + GFF3.

Fixture layout (1-based GFF3 coords on a single 5000bp ``chr1``):

    geneA  100-3500  +
      mRNA isoA  100-3500  +   (longest)
        intron A1: 500-799    (shared coord with B1)
        intron A2: 1200-1499  (unique to isoA)
      mRNA isoB  100-3000  +   (shorter, alt-spliced)
        intron B1: 500-799    (coord-duplicate of A1)
        intron B2: 1200-1799  (alt-spliced; only in isoB; what -i recovers)

Expected scored-intron counts (column-5 != "NA"):
    no flags  : 2   (A1, A2)
    -i        : 3   (A1, A2, B2)
    -d        : 2   (A1, A2 — B1 is omitted as ISOFORM before dup check)
    -i -d     : 4   (A1, A2, B2, plus B1 as a kept duplicate)
"""

import subprocess
import sys
import tempfile
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
FIXTURE_DIR = REPO_ROOT / "tests" / "data" / "isoforms"
FIXTURE_FASTA = FIXTURE_DIR / "synthetic.fa"
FIXTURE_GFF = FIXTURE_DIR / "synthetic.gff3"


# All four (allow-multiple-isoforms, include-duplicates) flag combinations.
# Tag is used in the temp output dir name; flags is the actual CLI form.
FLAG_COMBOS = [
    ("none", []),
    ("i", ["-i"]),
    ("d", ["-d"]),
    ("id", ["-i", "-d"]),
]

EXPECTED_SCORED_COUNT = {
    "none": 2,
    "i": 3,
    "d": 2,
    "id": 4,
}


def _intronic_cli():
    """Invoke the in-tree intronIC via the current Python interpreter."""
    return [
        sys.executable, "-c",
        "from intronIC.cli.main import main; import sys as _s; "
        "_s.argv = ['intronIC'] + _s.argv[1:]; main()",
    ]


def _run_classify(tmp_dir, species_name, mode, flags):
    """Run intronIC classify on the synthetic fixture and return the
    Path to the output bed.iic."""
    out_dir = Path(tmp_dir) / f"{mode}_{species_name}"
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = _intronic_cli() + [
        "classify",
        "-g", str(FIXTURE_FASTA),
        "-a", str(FIXTURE_GFF),
        "-n", species_name,
        "-o", str(out_dir),
        "-p", "1",
        "-f", "exon",
        f"--{mode}",
    ] + flags

    env = {"PYTHONPATH": str(REPO_ROOT / "src")}
    import os as _os
    full_env = dict(_os.environ)
    full_env.update(env)

    result = subprocess.run(
        cmd, capture_output=True, text=True, timeout=300, env=full_env,
    )
    if result.returncode != 0:
        pytest.fail(
            f"intronIC classify ({mode}, flags={flags}) failed:\n"
            f"STDOUT (tail):\n{result.stdout[-1500:]}\n"
            f"STDERR (tail):\n{result.stderr[-1500:]}"
        )

    bed_path = out_dir / f"{species_name}.bed.iic"
    if not bed_path.exists():
        pytest.fail(f"Expected bed.iic at {bed_path}, output dir contains "
                    f"{list(out_dir.iterdir())}")
    return bed_path


def _read_bed(bed_path):
    """Return list of bed.iic rows (each a list of tab-separated fields)."""
    rows = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            rows.append(line.split("\t"))
    return rows


def _scored_rows(rows):
    """Rows whose column-5 (SVM score) is not the omitted sentinel 'NA'."""
    return [r for r in rows if len(r) >= 5 and r[4] != "NA"]


def _coords(rows):
    """Set of (chrom, start, stop) from rows."""
    return {(r[0], int(r[1]), int(r[2])) for r in rows}


@pytest.fixture(scope="module")
def all_runs():
    """Run intronIC classify across (mode, flag_combo) and yield rows.

    Returned dict is keyed by (mode, tag); value is the list of bed.iic
    rows (post-split). Scoped to ``module`` so the four behaviour tests
    plus the equivalence test all share one set of subprocess runs.
    """
    runs = {}
    with tempfile.TemporaryDirectory(prefix="intronIC_iso_") as tmp_dir:
        for mode in ("in-memory", "streaming"):
            for tag, flags in FLAG_COMBOS:
                species = f"iso{tag.capitalize() or 'None'}"
                bed_path = _run_classify(tmp_dir, species, mode, flags)
                runs[(mode, tag)] = _read_bed(bed_path)
        yield runs


@pytest.mark.parametrize("mode", ["in-memory", "streaming"])
@pytest.mark.parametrize("tag", ["none", "i", "d", "id"])
def test_scored_count(all_runs, mode, tag):
    """The number of scored (non-NA) introns matches the truth table."""
    rows = all_runs[(mode, tag)]
    scored = _scored_rows(rows)
    expected = EXPECTED_SCORED_COUNT[tag]
    assert len(scored) == expected, (
        f"{mode} flags={tag}: expected {expected} scored, got {len(scored)}.\n"
        f"All rows:\n" + "\n".join("\t".join(r) for r in rows)
    )


@pytest.mark.parametrize("mode", ["in-memory", "streaming"])
def test_isoB_unique_intron_only_recovered_with_i(all_runs, mode):
    """The unique alt-spliced intron 1200-1799 (isoB-only) must be scored
    iff ``-i`` was passed."""
    target = ("chr1", 1199, 1799)  # bed.iic is 0-based half-open

    for tag, flags in FLAG_COMBOS:
        rows = all_runs[(mode, tag)]
        scored = _coords(_scored_rows(rows))
        if "-i" in flags:
            assert target in scored, (
                f"{mode} flags={tag}: isoB-only intron 1200-1799 should be "
                f"scored but isn't. Scored coords: {sorted(scored)}"
            )
        else:
            assert target not in scored, (
                f"{mode} flags={tag}: isoB-only intron 1200-1799 should NOT "
                f"be scored without -i. Scored coords: {sorted(scored)}"
            )


@pytest.mark.parametrize("mode", ["in-memory", "streaming"])
def test_duplicate_intron_emitted_only_with_i_and_d(all_runs, mode):
    """The shared intron coord 500-799 should appear twice (once per isoform)
    in the bed.iic only when both ``-i`` and ``-d`` are set.

    With ``-d`` alone, B1 is dropped as ISOFORM before the duplicate
    check ever runs, so it doesn't get a duplicate-emission slot. With
    ``-i`` alone, the duplicate is collapsed to a single representative.
    Only ``-i -d`` together emit every isoform's copy of a shared coord.
    """
    target = ("chr1", 499, 799)
    for tag, flags in FLAG_COMBOS:
        rows = all_runs[(mode, tag)]
        scored = _scored_rows(rows)
        copies_of_target = sum(
            1 for r in scored
            if (r[0], int(r[1]), int(r[2])) == target
        )
        if "-i" in flags and "-d" in flags:
            assert copies_of_target == 2, (
                f"{mode} flags={tag}: shared intron 500-799 should appear "
                f"twice (one per isoform) with -i -d; got {copies_of_target}"
            )
        else:
            assert copies_of_target == 1, (
                f"{mode} flags={tag}: shared intron 500-799 should appear "
                f"once; got {copies_of_target}"
            )


@pytest.mark.parametrize("tag", ["none", "i", "d", "id"])
def test_streaming_matches_in_memory(all_runs, tag):
    """Streaming and in-memory must produce identical bed.iic content
    (modulo the species-name prefix, which we strip).

    This extends the existing TestStreamingMatchesInMemory class to
    cover all four flag combinations, not just the no-flag case.
    """
    def _normalize(rows):
        out = []
        for r in rows:
            name = r[3]
            if "-" in name:
                name = name.split("-", 1)[1]
            new_r = list(r)
            new_r[3] = name
            out.append(tuple(new_r))
        return sorted(out)

    stream = _normalize(all_runs[("streaming", tag)])
    inmem = _normalize(all_runs[("in-memory", tag)])

    if stream != inmem:
        # Build a readable diff for the failure message.
        from difflib import unified_diff
        s_text = ["\t".join(map(str, r)) for r in stream]
        i_text = ["\t".join(map(str, r)) for r in inmem]
        diff = "\n".join(unified_diff(
            i_text, s_text, lineterm="",
            fromfile=f"in-memory[{tag}]",
            tofile=f"streaming[{tag}]",
        ))
        pytest.fail(
            f"streaming/in-memory bed.iic differ for flags={tag}:\n{diff}"
        )
