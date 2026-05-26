#!/usr/bin/env python3
"""Backfill rel_score (and type_id in meta.iic) for species run before the
v2.7.1 rel_score-propagation fix.

The bug:
    apply_continuous_per_intron_discount() rewrote adjusted_score but never
    updated rel_score. Gate-pass species shipped score_info.iic with
    rel_score = first_pass_svm - 90 (stale) instead of
    rel_score = adjusted_score - 90 (canonical).

    meta.iic copied that stale rel_score AND used the legacy 50%-threshold
    type_id, so it disagreed with score_info.iic for every gate-pass species.

This script reads each species' score_info.iic, recomputes
    rel_score = adjusted_score - 90
and propagates rel_score + type_id to the matching meta.iic.

Streaming pure-Python I/O (no pandas), parallel across species.
Idempotent: re-running is a no-op.

Usage:
    python backfill_v271_rel_score.py <root> [--workers N] [--dry-run]
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

THRESHOLD = 90.0


def _parse_header(line: str) -> dict[str, int]:
    return {n: i for i, n in enumerate(line.rstrip("\n").split("\t"))}


def backfill_species(species_dir: Path, dry_run: bool = False) -> dict:
    """Fix one species' score_info.iic + meta.iic. Returns counts dict."""
    sp = species_dir.name
    score_path = species_dir / f"{sp}.score_info.iic"
    meta_path = species_dir / f"{sp}.meta.iic"

    result = dict(species=sp, n_score=0, n_score_changed=0,
                  n_meta=0, n_meta_changed=0, status="ok", err="")

    if not score_path.exists():
        result["status"] = "no_score_info"
        return result

    # === Pass 1: read header, build name → (rel_score, type_id) map ===
    name_to_rel: dict[str, str] = {}
    name_to_type: dict[str, str] = {}
    n_changed = 0
    try:
        with open(score_path) as f:
            header_line = f.readline()
            col = _parse_header(header_line)
            rel_col = col.get("rel_score")
            adj_col = col.get("adjusted_score")
            type_col = col.get("type_id")
            name_col = col.get("name", 0)
            if rel_col is None or adj_col is None:
                result["status"] = "missing_cols"
                return result

            # Stream the file, collect (name, new_rel, type_id)
            rows = []
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(rel_col, adj_col):
                    rows.append(parts)
                    continue
                adj_str = parts[adj_col]
                try:
                    adj_val = float(adj_str)
                except ValueError:
                    rows.append(parts)
                    continue
                new_rel_str = f"{adj_val - THRESHOLD:.4f}"
                if parts[rel_col] != new_rel_str:
                    n_changed += 1
                    parts[rel_col] = new_rel_str
                name_to_rel[parts[name_col]] = new_rel_str
                if type_col is not None:
                    name_to_type[parts[name_col]] = parts[type_col]
                rows.append(parts)
    except Exception as e:
        result["status"] = "score_read_err"
        result["err"] = str(e)
        return result

    result["n_score"] = len(rows)
    result["n_score_changed"] = n_changed

    # === Write back score_info.iic only if any rows changed ===
    if not dry_run and n_changed > 0:
        try:
            tmp_fd, tmp_name = tempfile.mkstemp(
                suffix=".iic", dir=str(species_dir)
            )
            try:
                with os.fdopen(tmp_fd, "w") as tmp:
                    tmp.write(header_line)
                    for parts in rows:
                        tmp.write("\t".join(parts) + "\n")
                orig_mode = score_path.stat().st_mode
                shutil.move(tmp_name, score_path)
                os.chmod(score_path, orig_mode & 0o777)
            except Exception:
                try:
                    os.unlink(tmp_name)
                except OSError:
                    pass
                raise
        except Exception as e:
            result["status"] = "score_write_err"
            result["err"] = str(e)
            return result

    # === Pass 2: stream meta.iic, sync rel_score + type_id from the map ===
    if not meta_path.exists():
        return result

    try:
        with open(meta_path) as f:
            meta_header = f.readline()
            mh = _parse_header(meta_header)
            m_name = mh.get("name", 0)
            m_rel = mh.get("rel_score")
            m_type = mh.get("type_id")
            if m_rel is None:
                return result

            n_meta_rows = 0
            n_meta_changed = 0
            tmp_fd, tmp_name = tempfile.mkstemp(
                suffix=".iic", dir=str(species_dir)
            )
            try:
                with os.fdopen(tmp_fd, "w") as tmp:
                    tmp.write(meta_header)
                    for line in f:
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) > m_name:
                            n_meta_rows += 1
                            nm = parts[m_name]
                            changed = False
                            if nm in name_to_rel and len(parts) > m_rel:
                                if parts[m_rel] != name_to_rel[nm]:
                                    parts[m_rel] = name_to_rel[nm]
                                    changed = True
                            if (m_type is not None and nm in name_to_type
                                    and len(parts) > m_type
                                    and parts[m_type] != name_to_type[nm]):
                                parts[m_type] = name_to_type[nm]
                                changed = True
                            if changed:
                                n_meta_changed += 1
                        tmp.write("\t".join(parts) + "\n")
                if not dry_run and n_meta_changed > 0:
                    orig_mode = meta_path.stat().st_mode
                    shutil.move(tmp_name, meta_path)
                    os.chmod(meta_path, orig_mode & 0o777)
                else:
                    os.unlink(tmp_name)
            except Exception:
                try:
                    os.unlink(tmp_name)
                except OSError:
                    pass
                raise
            result["n_meta"] = n_meta_rows
            result["n_meta_changed"] = n_meta_changed
    except Exception as e:
        result["status"] = "meta_err"
        result["err"] = str(e)

    return result


def _worker(args):
    species_dir, dry_run = args
    return backfill_species(species_dir, dry_run=dry_run)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("root", type=Path,
                   help="Directory with per-species subdirs.")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--limit", type=int, default=None)
    args = p.parse_args()

    if not args.root.is_dir():
        sys.exit(f"Not a directory: {args.root}")

    species_dirs = sorted(d for d in args.root.iterdir() if d.is_dir())
    if args.limit:
        species_dirs = species_dirs[:args.limit]

    print(f"Backfilling {len(species_dirs)} species with {args.workers} workers"
          f" (dry_run={args.dry_run})", flush=True)

    n_score_changed = 0
    n_meta_changed = 0
    n_done = 0
    n_with_changes = 0
    n_errors = 0

    work = [(d, args.dry_run) for d in species_dirs]
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(_worker, w): w[0].name for w in work}
        for fut in as_completed(futs):
            sp = futs[fut]
            try:
                r = fut.result()
            except Exception as e:
                print(f"  [CRASH] {sp}: {e}", flush=True)
                n_errors += 1
                continue
            n_done += 1
            if r["status"] != "ok":
                n_errors += 1
                print(f"  [{r['status']}] {sp}: {r['err']}", flush=True)
                continue
            sc = r["n_score_changed"]
            mc = r["n_meta_changed"]
            if sc or mc:
                n_with_changes += 1
                n_score_changed += sc
                n_meta_changed += mc
                if n_done % 50 == 0 or sc > 5000:
                    print(f"  [{n_done}/{len(species_dirs)}] {sp}: "
                          f"score {sc}/{r['n_score']}  meta {mc}",
                          flush=True)

    print()
    print(f"Done. Species processed: {n_done}, with changes: {n_with_changes}, "
          f"errors: {n_errors}")
    print(f"Total score_info rewrites: {n_score_changed:,}")
    print(f"Total meta.iic rewrites:   {n_meta_changed:,}")


if __name__ == "__main__":
    main()
