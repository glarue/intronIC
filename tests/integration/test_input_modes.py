"""
Integration test: classify must accept all three input modes end-to-end
and stay consistent with the annotation path on equivalent intron sets.

Three modes are documented in `intronIC classify --help`:
  (1) -g + -a   annotation file (gff3/gtf) + genome — covered by other tests
  (2) -g + -b   BED file with intron coordinates + genome
  (3) -q        pre-extracted intron sequences in .iic format

This file covers (2) and (3): smoke-tests, equivalence vs annotation
(when the input set matches), and the duplicate-handling contract for
each mode. Fixtures are produced at test time by running annotation
classify on the bundled chr19 data, so we don't commit any binary
`.iic` blobs.
"""

import re
import subprocess
from pathlib import Path

import pytest


DATA_DIR = Path(__file__).parent.parent.parent / "src" / "intronIC" / "data"
TEST_GENOME = DATA_DIR / "test_data" / "Homo_sapiens.Chr19.Ensembl_91.fa.gz"
TEST_ANNOTATION = DATA_DIR / "test_data" / "Homo_sapiens.Chr19.Ensembl_91.gff3.gz"

_TAG_SUFFIX_RE = re.compile(r"(;\[[^\]]+\])+$")


def _strip_tag_suffix(label: str) -> str:
    """Drop trailing ;[n], ;[i], ;[c:N], ;[d], ;[o:XX] etc.

    Allows comparing labels across modes that re-derive tags from
    metadata; the score columns themselves are independent of the tags.
    """
    return _TAG_SUFFIX_RE.sub("", label)


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
def annotation_classify(intronIC_bin, tmp_path_factory):
    """Run annotation classify once per module to produce reference outputs.

    The resulting bed.iic and introns.iic feed the round-trip equivalence
    tests; the score_info.iic is the authoritative comparison target.
    """
    out_dir = tmp_path_factory.mktemp("annotation_classify")
    cmd = [
        intronIC_bin, "classify",
        "-g", str(TEST_GENOME),
        "-a", str(TEST_ANNOTATION),
        "-n", "AnnRef",
        "-o", str(out_dir),
        "-p", "2",
        "--in-memory",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if result.returncode != 0:
        pytest.fail(
            "annotation classify failed:\n"
            f"STDOUT:\n{result.stdout[-1500:]}\n"
            f"STDERR:\n{result.stderr[-1500:]}"
        )

    bed = out_dir / "AnnRef.bed.iic"
    seqs = out_dir / "AnnRef.introns.iic"
    score_info = out_dir / "AnnRef.score_info.iic"
    for p in (bed, seqs, score_info):
        if not p.exists():
            pytest.fail(f"annotation classify did not produce {p.name}")
    return {"dir": out_dir, "bed": bed, "seqs": seqs, "score_info": score_info}


def _scored_rows(score_info_path: Path):
    """Return list of (label_no_prefix, fields) for each scored intron."""
    rows = []
    with open(score_info_path) as f:
        f.readline()  # header
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or fields[2] in ("NA", ""):
                continue
            label = fields[0]
            if "-" in label:
                label = label.split("-", 1)[1]
            rows.append((label, fields))
    return rows


def _scored_count(score_info_path: Path) -> int:
    return len(_scored_rows(score_info_path))


@pytest.mark.skipif(
    not TEST_GENOME.exists() or not TEST_ANNOTATION.exists(),
    reason="Chr19 test data not found",
)
class TestBedInputSmoke:
    """`classify -g <fa> -b <bed>` must complete and emit scored introns."""

    @pytest.fixture(scope="class")
    def run_results(self, intronIC_bin, annotation_classify, tmp_path_factory):
        out_dir = tmp_path_factory.mktemp("classify_bed_smoke")
        cmd = [
            intronIC_bin, "classify",
            "-g", str(TEST_GENOME),
            "-b", str(annotation_classify["bed"]),
            "-n", "BedSmoke",
            "-o", str(out_dir),
            "-p", "2",
            "--in-memory",
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
        assert _scored_count(score_path) > 1000

    def test_meta_and_bed_files_exist(self, run_results):
        assert any(Path(run_results).glob("*meta*"))
        assert any(Path(run_results).glob("*.bed.iic"))


@pytest.mark.skipif(
    not TEST_GENOME.exists() or not TEST_ANNOTATION.exists(),
    reason="Chr19 test data not found",
)
class TestSequenceInputSmoke:
    """`classify -q <introns.iic>` must complete and emit scored introns."""

    @pytest.fixture(scope="class")
    def run_results(self, intronIC_bin, annotation_classify, tmp_path_factory):
        out_dir = tmp_path_factory.mktemp("classify_seq_smoke")
        cmd = [
            intronIC_bin, "classify",
            "-q", str(annotation_classify["seqs"]),
            "-n", "SeqSmoke",
            "-o", str(out_dir),
            "-p", "2",
            "--in-memory",
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
        assert _scored_count(score_path) > 1000

    def test_meta_file_exists(self, run_results):
        assert any(Path(run_results).glob("*meta*"))


@pytest.mark.skipif(
    not TEST_GENOME.exists() or not TEST_ANNOTATION.exists(),
    reason="Chr19 test data not found",
)
class TestBedAnnotationEquivalence:
    """BED classify on the same intron set as annotation must produce
    bit-identical scores.

    Construction: take annotation classify's bed.iic, filter it down to
    only the rows whose corresponding score_info entry has a non-NA SVM
    score (i.e., the introns that survived prefilter + IntronFilter +
    omitted-filter and got scored). Re-run classify in BED mode against
    that filtered BED. Because the BED set matches the annotation set,
    BG correction and the adaptive normalizer fit see the same points
    in both runs, so every per-intron score column should match exactly.

    Tags in column 1 are stripped before comparison: BED-input introns
    lack transcript metadata, so a U12-corrected intron in annotation
    mode carries `[c:N]` (correction recorded) while in BED mode it
    arrives already at the corrected coordinates and does not. The
    score values themselves are independent of which tags were applied.
    """

    @pytest.fixture(scope="class")
    def filtered_bed(self, annotation_classify, tmp_path_factory):
        """Build a BED with only the introns annotation classify scored."""
        out_dir = tmp_path_factory.mktemp("scored_only_bed")
        bed_out = out_dir / "scored_only.bed"
        # Names that have a non-NA SVM score in annotation's score_info
        scored_names = set()
        with open(annotation_classify["score_info"]) as f:
            f.readline()  # header
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 3 and fields[2] not in ("NA", ""):
                    scored_names.add(fields[0])

        # Filter bed.iic and emit standard BED6 (drop omitted_reason col 7)
        with open(annotation_classify["bed"]) as bed_in, open(bed_out, "w") as bed_w:
            for line in bed_in:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 4 and fields[3] in scored_names:
                    bed_w.write("\t".join(fields[:6]) + "\n")
        return bed_out

    @pytest.fixture(scope="class")
    def run_results(self, intronIC_bin, filtered_bed, annotation_classify, tmp_path_factory):
        out_dir = tmp_path_factory.mktemp("classify_bed_equiv")
        cmd = [
            intronIC_bin, "classify",
            "-g", str(TEST_GENOME),
            "-b", str(filtered_bed),
            "-n", "BedEquiv",
            "-o", str(out_dir),
            "-p", "2",
            "--in-memory",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            pytest.fail(
                "BED equivalence classify failed:\n"
                f"STDOUT:\n{result.stdout[-1500:]}\n"
                f"STDERR:\n{result.stderr[-1500:]}"
            )
        return {
            "ann_score_info": annotation_classify["score_info"],
            "bed_score_info": next(Path(out_dir).glob("*score_info*")),
        }

    def test_scored_count_matches_annotation(self, run_results):
        ann_n = _scored_count(run_results["ann_score_info"])
        bed_n = _scored_count(run_results["bed_score_info"])
        assert ann_n == bed_n, (
            f"BED-input scored {bed_n} introns but annotation scored {ann_n}; "
            "filtered BED should produce identical count"
        )

    def test_score_columns_bit_identical(self, run_results):
        """Every per-intron score column matches between modes."""
        ann_rows = {
            _strip_tag_suffix(label): fields
            for label, fields in _scored_rows(run_results["ann_score_info"])
        }
        bed_rows = {
            _strip_tag_suffix(label): fields
            for label, fields in _scored_rows(run_results["bed_score_info"])
        }

        common = sorted(set(ann_rows) & set(bed_rows))
        only_ann = set(ann_rows) - set(bed_rows)
        only_bed = set(bed_rows) - set(ann_rows)
        assert not only_ann and not only_bed, (
            f"label set mismatch (after tag-strip): "
            f"only-annotation={list(only_ann)[:3]}, only-bed={list(only_bed)[:3]}"
        )

        mismatches = []
        # All columns except the first (label, which carries tag suffixes)
        # must match exactly.
        for label in common:
            a = ann_rows[label]
            b = bed_rows[label]
            for col in range(1, len(a)):
                if a[col] != b[col]:
                    mismatches.append((label, col, a[col], b[col]))
                    break
        if mismatches:
            sample = mismatches[:5]
            pytest.fail(
                f"{len(mismatches)} introns have differing score columns:\n"
                + "\n".join(
                    f"  [{lab}] col{c}: ann={a} vs bed={b}"
                    for lab, c, a, b in sample
                )
            )


@pytest.mark.skipif(
    not TEST_GENOME.exists() or not TEST_ANNOTATION.exists(),
    reason="Chr19 test data not found",
)
class TestBedDuplicateFiltering:
    """BED mode must dedupe coordinate-duplicate rows by default and
    score every row when -d is passed.

    Two BED rows with identical (chrom, start, stop, strand) but
    distinct names should:
      - default: produce one scored entry (the first)
      - with -d: produce both, with the second carrying the [d] tag
    Unique rows always score regardless of the flag.
    """

    @pytest.fixture(scope="class")
    def dup_bed(self, tmp_path_factory):
        """A BED with one duplicate pair plus two unique rows."""
        out_dir = tmp_path_factory.mktemp("dup_bed")
        bed = out_dir / "dups.bed"
        bed.write_text(
            "19\t8969879\t8971663\tintron_a_dup1\t0.0\t-\n"
            "19\t8969879\t8971663\tintron_a_dup2\t0.0\t-\n"
            "19\t8967189\t8969774\tintron_b_unique\t0.0\t-\n"
            "19\t8943672\t8945496\tintron_c_unique\t0.0\t-\n"
        )
        return bed

    def _run(self, intronIC_bin, dup_bed, tmp_path, name, extra_args=()):
        out_dir = tmp_path / name
        out_dir.mkdir()
        cmd = [
            intronIC_bin, "classify",
            "-g", str(TEST_GENOME),
            "-b", str(dup_bed),
            "-n", name,
            "-o", str(out_dir),
            "-p", "1",
            "--in-memory",
            *extra_args,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            pytest.fail(
                f"BED dup classify ({name}) failed:\n"
                f"STDOUT:\n{result.stdout[-1500:]}\n"
                f"STDERR:\n{result.stderr[-1500:]}"
            )
        return next(out_dir.glob("*score_info*"))

    def test_default_filters_coord_duplicates(self, intronIC_bin, dup_bed, tmp_path):
        score_path = self._run(intronIC_bin, dup_bed, tmp_path, "DupDef")
        labels = [label for label, _ in _scored_rows(score_path)]
        labels_no_tag = [_strip_tag_suffix(lbl) for lbl in labels]
        assert sorted(labels_no_tag) == sorted([
            "intron_a_dup1", "intron_b_unique", "intron_c_unique"
        ]), f"unexpected scored set: {labels}"
        # The duplicate variant must NOT appear by default
        assert not any("intron_a_dup2" in lbl for lbl in labels)

    def test_dash_d_includes_duplicates(self, intronIC_bin, dup_bed, tmp_path):
        score_path = self._run(intronIC_bin, dup_bed, tmp_path, "DupD", extra_args=("-d",))
        labels = [label for label, _ in _scored_rows(score_path)]
        # Both duplicate variants present, second tagged [d]
        assert any(lbl == "intron_a_dup1" for lbl in labels), (
            f"first dup missing: {labels}"
        )
        assert any(
            lbl.startswith("intron_a_dup2") and "[d]" in lbl for lbl in labels
        ), f"second dup missing or untagged: {labels}"
        assert any(_strip_tag_suffix(lbl) == "intron_b_unique" for lbl in labels)
        assert any(_strip_tag_suffix(lbl) == "intron_c_unique" for lbl in labels)


@pytest.mark.skipif(
    not TEST_GENOME.exists() or not TEST_ANNOTATION.exists(),
    reason="Chr19 test data not found",
)
class TestBedDuplicatesExcludedFromBackground:
    """BED-mode species background must accumulate over unique
    coordinates, not raw input rows.

    Otherwise an input with many copies of the same intron would skew
    the empirical U2 PWM toward that intron's composition. We construct
    a heavily-duplicated BED (one intron copied 1000 times alongside
    a base of unique introns) and assert that:
      - the BG-accumulator count reported in the log equals the unique
        coord count, not the total row count;
      - the corrected scores for the unique introns are identical to a
        run on the unique-only BED.
    Both conditions follow from the dedup guard in
    apply_species_background.
    """

    @pytest.fixture(scope="class")
    def beds(self, annotation_classify, tmp_path_factory):
        """Build (unique, dupe-heavy) BED pair from annotation's
        scored_only set; pick enough unique introns to exceed the
        500-intron min_introns BG threshold so BG actually runs."""
        out_dir = tmp_path_factory.mktemp("bg_dedup_bed")

        # Take the first 700 scored introns from annotation as the unique base
        scored_names = []
        with open(annotation_classify["score_info"]) as f:
            f.readline()
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 3 and fields[2] not in ("NA", ""):
                    scored_names.append(fields[0])
                if len(scored_names) >= 700:
                    break

        # Read bed.iic, keep only those rows, emit BED6
        unique_bed = out_dir / "unique.bed"
        scored_set = set(scored_names)
        first_row = None
        with open(annotation_classify["bed"]) as bed_in, open(unique_bed, "w") as bw:
            for line in bed_in:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 4 and fields[3] in scored_set:
                    bw.write("\t".join(fields[:6]) + "\n")
                    if first_row is None:
                        first_row = fields[:6]

        # Dupe-heavy: same as unique + 1000 copies of the first intron under
        # distinct names so each row passes BED parsing as a separate entry.
        # The IntronFilter / BG dedup keys on (chrom, start, stop, strand),
        # so name divergence does not save them from being deduped.
        dupe_bed = out_dir / "dupe_heavy.bed"
        with open(unique_bed) as src, open(dupe_bed, "w") as dst:
            dst.write(src.read())
            chrom, start, stop, _name, score, strand = first_row
            for i in range(1000):
                dst.write(f"{chrom}\t{start}\t{stop}\tdupe_{i}\t{score}\t{strand}\n")

        return {"unique": unique_bed, "dupe_heavy": dupe_bed}

    def _run(self, intronIC_bin, bed, tmp_path, name):
        out_dir = tmp_path / name
        out_dir.mkdir()
        cmd = [
            intronIC_bin, "classify",
            "-g", str(TEST_GENOME),
            "-b", str(bed),
            "-n", name,
            "-o", str(out_dir),
            "-p", "2",
            "--in-memory",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            pytest.fail(
                f"BED bg-dedup classify ({name}) failed:\n"
                f"STDOUT:\n{result.stdout[-1500:]}\n"
                f"STDERR:\n{result.stderr[-1500:]}"
            )
        return {
            "score_info": next(out_dir.glob("*score_info*")),
            "stdout": result.stdout,
            "stderr": result.stderr,
        }

    def test_bg_count_excludes_duplicates(self, intronIC_bin, beds, tmp_path):
        """BG accumulator log should report the unique-coord count for
        both BEDs, not the inflated row count of the dupe-heavy one."""
        unique = self._run(intronIC_bin, beds["unique"], tmp_path, "BgUnique")
        dupes = self._run(intronIC_bin, beds["dupe_heavy"], tmp_path, "BgDupes")

        # The "computing empirical U2 from N introns" line is emitted only
        # if N >= min_introns (500); we built unique=700, so it must fire.
        log_re = re.compile(r"computing empirical U2 from (\d+) introns")
        u_match = log_re.search(unique["stdout"])
        d_match = log_re.search(dupes["stdout"])
        assert u_match is not None, (
            "unique BED did not produce BG correction log line — "
            "did the unique count drop below min_introns?\n"
            + unique["stdout"][-1500:]
        )
        assert d_match is not None, (
            "dupe-heavy BED did not produce BG correction log line\n"
            + dupes["stdout"][-1500:]
        )
        assert u_match.group(1) == d_match.group(1), (
            f"BG count differs: unique={u_match.group(1)} "
            f"vs dupe-heavy={d_match.group(1)} — duplicates are leaking "
            "into the empirical U2 background"
        )

    def test_unique_intron_scores_unchanged_by_duplicates(
        self, intronIC_bin, beds, tmp_path
    ):
        """Adding 1000 dupes of one intron must not change the
        BG-corrected scores of the other ~700 unique introns."""
        unique = self._run(intronIC_bin, beds["unique"], tmp_path, "BgUnique2")
        dupes = self._run(intronIC_bin, beds["dupe_heavy"], tmp_path, "BgDupes2")

        u_rows = {
            _strip_tag_suffix(label): fields
            for label, fields in _scored_rows(unique["score_info"])
        }
        d_rows = {
            _strip_tag_suffix(label): fields
            for label, fields in _scored_rows(dupes["score_info"])
        }
        common = sorted(set(u_rows) & set(d_rows))
        assert common, "no labels in common between the two runs"

        mismatches = []
        for label in common:
            u, d = u_rows[label], d_rows[label]
            for col in range(1, len(u)):
                if u[col] != d[col]:
                    mismatches.append((label, col, u[col], d[col]))
                    break
        if mismatches:
            sample = mismatches[:5]
            pytest.fail(
                f"{len(mismatches)} introns scored differently across "
                f"unique vs dupe-heavy BG samples:\n"
                + "\n".join(
                    f"  [{lab}] col{c}: unique={u} vs dupe={d}"
                    for lab, c, u, d in sample
                )
            )


@pytest.mark.skipif(
    not TEST_GENOME.exists() or not TEST_ANNOTATION.exists(),
    reason="Chr19 test data not found",
)
class TestSequenceInputDoesNotDedup:
    """Sequence-input mode (-q) does NOT dedupe — each .iic entry is
    scored as-is, even when two entries have identical sequence content.

    Documented design: -q has no coordinate metadata to dedupe by, so
    every row in the input is treated as a distinct intron. Callers
    that want dedup must do it before passing the file in.
    """

    @pytest.fixture(scope="class")
    def dup_seqs(self, annotation_classify, tmp_path_factory):
        """Take the first 2 sequences and append two more copies under
        different name prefixes, yielding 6 rows from 2 unique sequences."""
        out_dir = tmp_path_factory.mktemp("dup_seqs")
        out = out_dir / "dups.introns.iic"
        with open(annotation_classify["seqs"]) as f:
            head = [next(f), next(f)]
        with open(out, "w") as w:
            for line in head:
                w.write(line)
            for line in head:
                w.write(re.sub(r"^[^-]*-", "DUPA-", line))
            for line in head:
                w.write(re.sub(r"^[^-]*-", "DUPB-", line))
        return out

    def test_every_input_entry_is_scored(self, intronIC_bin, dup_seqs, tmp_path):
        out_dir = tmp_path / "SeqDup"
        out_dir.mkdir()
        cmd = [
            intronIC_bin, "classify",
            "-q", str(dup_seqs),
            "-n", "SeqDup",
            "-o", str(out_dir),
            "-p", "1",
            "--in-memory",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            pytest.fail(
                "seq dup classify failed:\n"
                f"STDOUT:\n{result.stdout[-1500:]}\n"
                f"STDERR:\n{result.stderr[-1500:]}"
            )
        score_path = next(Path(out_dir).glob("*score_info*"))
        rows = _scored_rows(score_path)
        # 6 input rows, all scored — sequence mode does not dedupe
        assert len(rows) == 6, (
            f"expected 6 scored rows (no dedup in -q mode), got {len(rows)}: "
            f"{[r[0] for r in rows]}"
        )
