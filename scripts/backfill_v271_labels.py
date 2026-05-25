"""Backfill v2.7.1 unified labels into existing score_info.iic files.

For each species output directory that has a score_info.iic but lacks the
new type_id / confidence / history columns, compute the labels from the
existing first_pass_svm (or svm_score fallback for Tier F) + adjusted_score
columns, write them back into score_info.iic, and update metrics.iic.json
with the corresponding counts.

Idempotent: a species whose score_info.iic already has the new columns is
left untouched (unless --force is passed).

Usage:
    python scripts/backfill_v271_labels.py --root /path/to/v2.7_runs [--force]
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from intronIC.scoring.labeling import assign_labels, count_by_label


NEW_COLS = ("type_id", "confidence", "history")


def backfill_one(out_dir: Path, force: bool = False) -> tuple[str, dict | None]:
    """Returns (status, label_counts).

    status is one of "skipped_no_si", "skipped_already_labeled", "ok", "error".
    """
    si_files = list(out_dir.glob("*.score_info.iic"))
    if not si_files:
        return ("skipped_no_si", None)
    si = si_files[0]

    # Cheap sniff: only read the header first
    with open(si) as f:
        header = f.readline().rstrip("\n").split("\t")
    if all(c in header for c in NEW_COLS) and not force:
        return ("skipped_already_labeled", None)

    df = pd.read_csv(si, sep="\t", low_memory=False)
    if "adjusted_score" not in df.columns:
        # Can't compute labels without adjusted_score
        return ("error", None)
    df["adjusted_score"] = pd.to_numeric(df["adjusted_score"], errors="coerce").fillna(0)
    if "first_pass_svm" in df.columns:
        fp = pd.to_numeric(df["first_pass_svm"], errors="coerce")
        fp = fp.fillna(df["adjusted_score"])
    else:
        # Tier F species — first_pass is implicit (== svm_score == adjusted pre-discount)
        # Use svm_score if present, else adjusted_score
        if "svm_score" in df.columns:
            fp = pd.to_numeric(df["svm_score"], errors="coerce").fillna(df["adjusted_score"])
        else:
            fp = df["adjusted_score"]

    labels = [assign_labels(float(f), float(a)) for f, a in zip(fp, df["adjusted_score"])]
    df["type_id"] = [L.type_id for L in labels]
    df["confidence"] = [L.confidence for L in labels]
    df["history"] = [L.history for L in labels]
    df.to_csv(si, sep="\t", index=False, float_format="%.6f")

    counts = count_by_label(labels)

    # Also update metrics.iic.json if present
    metrics_files = list(out_dir.glob("*.metrics.iic.json"))
    if metrics_files:
        mfile = metrics_files[0]
        try:
            m = json.load(open(mfile))
            m.update(counts)
            with open(mfile, "w") as f:
                json.dump(m, f, indent=2, default=str)
        except Exception as e:
            print(f"  WARN: failed to update {mfile.name}: {e}", file=sys.stderr)

    return ("ok", counts)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, required=True,
                    help="Root of v2.7_runs directory to backfill")
    ap.add_argument("--force", action="store_true",
                    help="Re-backfill even species that already have new columns")
    args = ap.parse_args()

    if not args.root.is_dir():
        sys.exit(f"root not a directory: {args.root}")

    species_dirs = [d for d in args.root.iterdir() if d.is_dir()]
    print(f"Scanning {len(species_dirs)} subdirs under {args.root}")

    stats = {"skipped_no_si": 0, "skipped_already_labeled": 0, "ok": 0, "error": 0}
    agg = {k: 0 for k in (
        "u12_count","u12_strong_count","u12_borderline_count","u12_promoted_count",
        "u2_count","u2_strong_count","u2_borderline_count","u2_demoted_count")}

    for i, d in enumerate(sorted(species_dirs), 1):
        status, counts = backfill_one(d, force=args.force)
        stats[status] += 1
        if counts:
            for k, v in counts.items():
                agg[k] = agg.get(k, 0) + v
        if i % 25 == 0 or status == "ok":
            print(f"  [{i:>4}/{len(species_dirs)}] {d.name[:50]:<50} status={status}")

    print(f"\n=== summary ===")
    for k, v in stats.items():
        print(f"  {k}: {v}")
    print(f"\n=== aggregate label counts across backfilled species ===")
    for k, v in agg.items():
        print(f"  {k}: {v:,}")


if __name__ == "__main__":
    main()
