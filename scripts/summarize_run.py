#!/usr/bin/env python3
"""
Summarize intronIC classification results from one or more output directories.

Reads the .metrics.iic.json files written by intronIC classify runs and
produces a compact table of key metrics.

Usage:
    # Single run directory:
    python scripts/summarize_run.py models/test_runs_linear_3d_rebuilt/

    # Multiple directories (compare models):
    python scripts/summarize_run.py models/test_runs_linear_3d/ models/test_runs_rbf_8d/

    # Filter to specific species:
    python scripts/summarize_run.py models/test_runs/ -s homo_sapiens drosophila_melanogaster

    # TSV output for downstream analysis:
    python scripts/summarize_run.py models/test_runs/ --tsv

    # Include boundary breakdown:
    python scripts/summarize_run.py models/test_runs/ --boundaries
"""

import argparse
import json
import sys
from pathlib import Path


def load_metrics(metrics_path: Path) -> dict:
    """Load a single .metrics.iic.json file."""
    with open(metrics_path) as f:
        return json.load(f)


def extract_summary(metrics: dict, species: str, run_label: str) -> dict:
    """Extract key fields from metrics into a flat summary dict."""
    cv = metrics.get("cluster_validation", {})
    boundaries = metrics.get("u12_boundaries", {})

    return {
        "run": run_label,
        "species": species,
        "u12_total": metrics.get("u12_count", 0),
        "u12_confident": cv.get("n_confident_u12", metrics.get("high_confidence_u12", 0)),
        "u12_atac": boundaries.get("AT-AC", 0),
        "u12_gtag": boundaries.get("GT-AG", 0),
        "u12_other": sum(
            v for k, v in boundaries.items() if k not in ("GT-AG", "AT-AC")
        ),
        "scored": metrics.get("total_scored", 0),
        "total_introns": metrics.get("total_introns", 0),
        "valley_depth": cv.get("valley_depth"),
        "has_valley": cv.get("has_valley"),
        "regime": cv.get("regime", ""),
        "centroid_sigma": cv.get("centroid_sigma"),
        "threshold": metrics.get("threshold", 90.0),
        "model_path": metrics.get("model_path", ""),
    }


def find_metrics_files(run_dir: Path, species_filter: list = None) -> list:
    """Find all .metrics.iic.json files in a directory."""
    files = sorted(run_dir.glob("*.metrics.iic.json"))
    if species_filter:
        files = [
            f for f in files
            if f.name.replace(".metrics.iic.json", "") in species_filter
        ]
    return files


def format_table(summaries: list, show_boundaries: bool = False) -> str:
    """Format summaries as an aligned text table."""
    if not summaries:
        return "No results found."

    multi_run = len(set(s["run"] for s in summaries)) > 1

    # When comparing multiple runs, sort by species then run for easy comparison
    if multi_run:
        summaries = sorted(summaries, key=lambda s: (s["species"], s["run"]))

    lines = []

    # Header
    cols = []
    if multi_run:
        cols.append(("Run", 30))
    cols += [
        ("Species", 35),
        ("U12", 5),
        ("AT-AC", 6),
        ("GT-AG", 6),
    ]
    if show_boundaries:
        cols.append(("Other", 6))
    cols += [
        ("Scored", 10),
        ("Valley", 7),
        ("Regime", 8),
    ]

    header = "  ".join(f"{name:<{width}s}" for name, width in cols)
    lines.append(header)
    lines.append("-" * len(header))

    for s in summaries:
        parts = []
        if multi_run:
            parts.append(f"{s['run']:<30s}")
        parts.append(f"{s['species']:<35s}")
        parts.append(f"{s['u12_confident']:>5d}")
        parts.append(f"{s['u12_atac']:>6d}")
        parts.append(f"{s['u12_gtag']:>6d}")
        if show_boundaries:
            parts.append(f"{s['u12_other']:>6d}")
        parts.append(f"{s['scored']:>10,d}")

        vd = s["valley_depth"]
        parts.append(f"{vd:>7.3f}" if vd is not None else f"{'N/A':>7s}")
        parts.append(f"{s['regime']:>8s}")

        lines.append("  ".join(parts))

    return "\n".join(lines)


def format_tsv(summaries: list) -> str:
    """Format summaries as TSV for downstream analysis."""
    if not summaries:
        return ""

    fields = [
        "run", "species", "u12_confident", "u12_atac", "u12_gtag",
        "u12_other", "scored", "valley_depth", "has_valley", "regime",
        "centroid_sigma", "threshold",
    ]
    lines = ["\t".join(fields)]
    for s in summaries:
        vals = []
        for f in fields:
            v = s.get(f, "")
            if v is None:
                v = ""
            elif isinstance(v, float):
                v = f"{v:.4f}"
            else:
                v = str(v)
            vals.append(v)
        lines.append("\t".join(vals))
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Summarize intronIC classification results."
    )
    parser.add_argument(
        "run_dirs",
        nargs="+",
        type=Path,
        help="One or more output directories containing .metrics.iic.json files",
    )
    parser.add_argument(
        "-s", "--species",
        nargs="+",
        default=None,
        help="Filter to specific species names",
    )
    parser.add_argument(
        "--tsv",
        action="store_true",
        help="Output as TSV instead of aligned table",
    )
    parser.add_argument(
        "--boundaries",
        action="store_true",
        help="Show non-canonical boundary counts in table",
    )
    args = parser.parse_args()

    all_summaries = []
    for run_dir in args.run_dirs:
        if not run_dir.is_dir():
            print(f"Warning: {run_dir} is not a directory, skipping", file=sys.stderr)
            continue

        run_label = run_dir.name
        metrics_files = find_metrics_files(run_dir, args.species)

        if not metrics_files:
            print(f"Warning: no .metrics.iic.json files in {run_dir}", file=sys.stderr)
            continue

        for mf in metrics_files:
            species = mf.name.replace(".metrics.iic.json", "")
            metrics = load_metrics(mf)
            summary = extract_summary(metrics, species, run_label)
            all_summaries.append(summary)

    if args.tsv:
        print(format_tsv(all_summaries))
    else:
        print(format_table(all_summaries, show_boundaries=args.boundaries))


if __name__ == "__main__":
    main()
