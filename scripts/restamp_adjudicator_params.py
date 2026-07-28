#!/usr/bin/env python3
"""Re-stamp a pmotif_adjudicated bundle's adjudicator params — a GATE change, not a retrain.

Per CLAUDE.md: "Gate changes are a re-stamp, not a retrain (calibrated Platt/anchors preserved)". This
script is the mechanical form of that rule. It:

  1. backs up the bundle to ``*.bak_<tag>_<date>``;
  2. ADDS any AdjudicatorParams field the code has but the bundle lacks, taking the code default;
  3. leaves EVERY already-stored value untouched — including the calibrated ones;
  4. updates ``params_version`` + ``provenance.adjudicator_params_version``, recording the previous
     value in ``provenance.supersedes_adjudicator_params_version``;
  5. verifies the ensemble is byte-for-byte behaviourally identical (same margins on a fixed probe).

**Why (2)/(3) and not ``asdict(AdjudicatorParams())``.** The stored calibration is NOT the code defaults.
As of 2026-07 the bundle carries ``platt_a=2.7958 / platt_c=-1.1778`` while the module's documented
defaults are the rounded ``2.796 / -1.178``. Overwriting with defaults would silently degrade the Platt
calibration on every genome — the exact failure this script exists to prevent. Stored values win, always.

Usage:
    python scripts/restamp_adjudicator_params.py --tag pre_lowk [--bundle PATH] [--dry-run]
"""
from __future__ import annotations

import argparse
import dataclasses
import shutil
import sys
from datetime import date
from pathlib import Path

import joblib
import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

from intronIC.scoring.species_adjudicator import (  # noqa: E402
    AdjudicatorParams, ADJUDICATOR_PARAMS_VERSION, ensemble_margin,
)

DEFAULT_BUNDLE = REPO / "src" / "intronIC" / "data" / "default_pretrained.model.pkl"


def _probe_margins(bundle, n_features: int = 6, n_rows: int = 256) -> np.ndarray:
    """Deterministic behavioural fingerprint of the ensemble (fixed seed, fixed shape)."""
    x = np.random.RandomState(20260727).normal(0.0, 2.0, (n_rows, n_features))
    return ensemble_margin(bundle["ensemble"].models, x)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle", type=Path, default=DEFAULT_BUNDLE)
    ap.add_argument("--tag", required=True,
                    help="short backup tag, e.g. 'pre_lowk' -> *.bak_pre_lowk_<YYYY-MM-DD>")
    ap.add_argument("--dry-run", action="store_true", help="report the diff, write nothing")
    args = ap.parse_args()

    if not args.bundle.exists():
        print(f"ERROR: bundle not found: {args.bundle}", file=sys.stderr)
        return 2

    bundle = joblib.load(args.bundle)
    if bundle.get("scoring_mode") != "pmotif_adjudicated":
        print(f"ERROR: not a pmotif_adjudicated bundle (scoring_mode="
              f"{bundle.get('scoring_mode')!r})", file=sys.stderr)
        return 2

    stored = dict(bundle.get("adjudicator_params") or {})
    code = dataclasses.asdict(AdjudicatorParams())
    old_version = stored.get("params_version") or (
        bundle.get("provenance", {}).get("adjudicator_params_version"))

    added = {k: code[k] for k in code if k not in stored}
    orphaned = sorted(set(stored) - set(code) - {"params_version"})
    preserved = {k: (stored[k], code[k]) for k in set(stored) & set(code)
                 if k != "params_version" and stored[k] != code[k]}

    print(f"bundle : {args.bundle}")
    print(f"version: {old_version}  ->  {ADJUDICATOR_PARAMS_VERSION}")
    print(f"\nADDED (new code fields, taking the code default):")
    for k, v in sorted(added.items()):
        print(f"   + {k:26s} {v!r}")
    print(f"\nPRESERVED (stored calibration kept; code default NOT applied):")
    for k, (sv, cv) in sorted(preserved.items()):
        print(f"   = {k:26s} {sv!r}   (code default {cv!r} ignored)")
    if orphaned:
        print(f"\nORPHANED (stored but no longer a code field; left in place, from_dict ignores them):")
        for k in orphaned:
            print(f"   ? {k:26s} {stored[k]!r}")
    if not added and old_version == ADJUDICATOR_PARAMS_VERSION:
        print("\nNothing to do — bundle already current.")
        return 0
    if args.dry_run:
        print("\n--dry-run: nothing written.")
        return 0

    before = _probe_margins(bundle)

    backup = args.bundle.with_suffix(args.bundle.suffix + f".bak_{args.tag}_{date.today()}")
    if backup.exists():
        print(f"ERROR: backup already exists, refusing to overwrite: {backup}", file=sys.stderr)
        return 2
    shutil.copy2(args.bundle, backup)
    print(f"\nbackup : {backup.name}")

    stored.update(added)
    stored["params_version"] = ADJUDICATOR_PARAMS_VERSION
    bundle["adjudicator_params"] = stored
    prov = dict(bundle.get("provenance") or {})
    if old_version:
        prov["supersedes_adjudicator_params_version"] = old_version
    prov["adjudicator_params_version"] = ADJUDICATOR_PARAMS_VERSION
    bundle["provenance"] = prov

    joblib.dump(bundle, args.bundle, compress=3)

    # ---- verify: params applied, ensemble behaviourally untouched ----
    rt = joblib.load(args.bundle)
    got = rt["adjudicator_params"]
    assert got["params_version"] == ADJUDICATOR_PARAMS_VERSION
    for k, (sv, _) in preserved.items():
        assert got[k] == sv, f"calibration drift on {k}: {got[k]!r} != {sv!r}"
    after = _probe_margins(rt)
    assert np.array_equal(before, after), "ensemble margins changed — re-stamp corrupted the model"
    resolved = AdjudicatorParams.from_dict(got)
    print(f"verify : ensemble margins identical over {len(before)} probes; "
          f"platt=({resolved.platt_a}, {resolved.platt_c}); "
          f"anchors=({resolved.loss_ceiling_z}, {resolved.bearer_floor_z}); "
          f"lowk=({resolved.lowk_gate_enabled}, {resolved.lowk_bg_fdr_threshold})")
    print("OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
