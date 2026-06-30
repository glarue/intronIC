#!/usr/bin/env python3
"""Post-hoc repair of ``metrics.iic.json`` classification COUNTS.

Recomputes ``u12_count`` / ``u2_count`` / ``high_confidence_u12`` / the two
percentages / ``u12_boundaries`` / ``u2_boundaries`` of a ``metrics.iic.json``
directly from its finalized sibling ``meta.iic``, using the authoritative
``intronIC.file_io.writers`` helpers (``count_calls_from_meta`` +
``summarize_boundaries_from_meta``) — the same single source of truth a fresh run
uses, so a patched JSON matches a fresh run's count fields.

Use to repair files written before the count-consolidation fix (stale first-pass
``u12_count``; ``high_confidence_u12`` == total called; ``u2_count`` / ``u2_boundaries``
including omitted introns). The per-intron output (score_info / meta / bed) is the
source of truth and is NOT touched — only the JSON summary counts are corrected.

Usage:
    python scripts/patch_metrics_counts.py <metrics.iic.json | dir> [...]
"""
import glob
import json
import os
import sys
from pathlib import Path

from intronIC.file_io.writers import (
    count_calls_from_meta,
    summarize_boundaries_from_meta,
)

_SUFFIX = ".metrics.iic.json"


def _sibling_meta(metrics_path: Path):
    base = metrics_path.name[: -len(_SUFFIX)] if metrics_path.name.endswith(_SUFFIX) else metrics_path.stem
    for cand in (metrics_path.with_name(base + ".meta.iic"),
                 metrics_path.with_name(base + ".meta.iic.gz")):
        if cand.exists():
            return cand
    g = (glob.glob(str(metrics_path.parent / "*.meta.iic"))
         + glob.glob(str(metrics_path.parent / "*.meta.iic.gz")))
    return Path(g[0]) if g else None


def patch_metrics(metrics_path) -> dict:
    """Repair one metrics.iic.json in place from its sibling meta.iic. Returns the
    patched dict. Raises FileNotFoundError if no sibling meta.iic is found."""
    metrics_path = Path(metrics_path)
    meta = _sibling_meta(metrics_path)
    if meta is None:
        raise FileNotFoundError(f"no sibling meta.iic for {metrics_path}")

    m = json.load(open(metrics_path))
    thr = m.get("threshold", 90.0)
    n_u12, n_hc, n_u2 = count_calls_from_meta(meta)
    u12_b, u2_b, u12_hc_b, u12_hc_feat = summarize_boundaries_from_meta(meta, thr)
    n_scored = n_u12 + n_u2

    m["u12_count"] = n_u12
    m["u2_count"] = n_u2
    m["high_confidence_u12"] = n_hc
    m["u12_percentage"] = (n_u12 / n_scored * 100) if n_scored else 0.0
    m["high_confidence_percentage"] = (n_hc / n_scored * 100) if n_scored else 0.0
    m["u12_boundaries"] = u12_b
    m["u2_boundaries"] = u2_b
    m["high_confidence_u12_boundaries"] = u12_hc_b
    m["high_confidence_u12_by_feature"] = u12_hc_feat

    # Consistency (same guarantees the live finalizer asserts).
    assert n_hc <= n_u12, f"HC ({n_hc}) > u12_count ({n_u12}) in {metrics_path}"
    assert sum(u12_hc_b.values()) <= n_hc, f"HC boundaries > HC ({n_hc}) in {metrics_path}"
    assert sum(u12_hc_feat.values()) <= n_hc, f"HC-by-feature > HC ({n_hc}) in {metrics_path}"
    ts = m.get("total_scored")
    if ts is not None and ts != n_scored:
        print(f"  WARN {metrics_path.name}: total_scored={ts} != u12+u2={n_scored}", file=sys.stderr)

    with open(metrics_path, "w") as f:
        json.dump(m, f, indent=2)
    return m


def main(argv):
    paths = []
    for a in argv:
        if os.path.isdir(a):
            paths += glob.glob(os.path.join(a, "**", f"*{_SUFFIX}"), recursive=True)
        else:
            paths.append(a)
    if not paths:
        print(__doc__)
        return 1
    for p in sorted(paths):
        m = patch_metrics(p)
        print(f"patched {p}: u12_count={m['u12_count']} u2_count={m['u2_count']} "
              f"high_confidence_u12={m['high_confidence_u12']}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
