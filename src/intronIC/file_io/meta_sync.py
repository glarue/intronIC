"""Sync final per-intron calls from a rewritten ``score_info.iic`` back into ``meta.iic``.

The file-side post-process scoring modes (``raw_gated``, ``pmotif_adjudicated``) compute their
final ``type_id`` / ``rel_score`` AFTER ``meta.iic`` is already on disk (written from the in-memory
introns' first-pass values). This utility rewrites ``meta.iic``'s ``rel_score`` and ``type_id``
columns from the on-disk ``score_info`` DataFrame so the downstream ``metrics.iic.json`` boundary
tables (which read meta ``type_id``) reflect any post-process call flips.

Extracted from the (removed) mode-separation pipeline so it survives independently of the z stack.
"""
from __future__ import annotations

import os
import shutil
import tempfile
from pathlib import Path

import pandas as pd

__all__ = ["sync_meta_from_score_info"]


def sync_meta_from_score_info(meta_path: Path, df: "pd.DataFrame", messenger=None) -> None:
    """Propagate ``rel_score`` and ``type_id`` from a score_info DataFrame to ``meta.iic``.

    ``meta.iic`` historically copies these from the in-memory ``Intron`` object, computed before the
    post-process runs. This rewrites both columns from the on-disk ``score_info`` (post-process)
    values, keyed by intron ``name``. Rows absent from ``df`` pass through unchanged.
    """
    name_to_rel = dict(zip(df["name"].astype(str), df["rel_score"].astype(float)))
    name_to_type = dict(zip(df["name"].astype(str), df["type_id"].astype(str)))

    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".iic",
        dir=meta_path.parent, delete=False,
    )
    n_synced = 0
    try:
        with open(meta_path) as f_in:
            header_line = f_in.readline()
            tmp.write(header_line)
            hdr = header_line.rstrip("\n").split("\t")
            try:
                name_col = hdr.index("name")
                rel_col = hdr.index("rel_score")
                type_col = hdr.index("type_id")
            except ValueError:
                tmp.close()
                Path(tmp.name).unlink(missing_ok=True)
                return

            for line in f_in:
                parts = line.rstrip("\n").split("\t")
                if name_col < len(parts):
                    nm = parts[name_col]
                    if nm in name_to_rel:
                        parts[rel_col] = f"{name_to_rel[nm]:.4f}"
                        parts[type_col] = name_to_type[nm]
                        n_synced += 1
                tmp.write("\t".join(parts) + "\n")
        tmp.close()
        orig_mode = meta_path.stat().st_mode if meta_path.exists() else 0o644
        shutil.move(tmp.name, meta_path)
        os.chmod(meta_path, orig_mode & 0o777)
        if messenger is not None:
            messenger.info(
                f"[meta-sync] rewrote rel_score + type_id for "
                f"{n_synced} rows in {meta_path.name}"
            )
    except Exception:
        tmp.close()
        Path(tmp.name).unlink(missing_ok=True)
        raise
