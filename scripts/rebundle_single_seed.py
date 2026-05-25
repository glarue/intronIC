"""Rebuild default_pretrained.model.pkl with only the first seed group.

Drops seeds 1042 and 2042 from both first_pass_seeds and seeds (second-pass),
keeping seed 42's 42 models for each. Reduces total models per pass from 126
to 42 with no measurable change in classification (validated by the
ensemble-size sweep — task #189).

Usage:
    python scripts/rebundle_single_seed.py [--in PATH] [--out PATH] [--keep SEED]

Defaults match the canonical v2.7 → v2.7.1 transition: read the shipped
default bundle, write to default_pretrained.v271.model.pkl alongside it
(does NOT overwrite). Production swap is a separate manual step.
"""
import argparse
import copy
import hashlib
import json
import pickle
import shutil
import sys
import time
from pathlib import Path


def trim_bundle(bundle, keep_seed=42):
    """Return a new bundle dict with only the given seed retained in both passes."""
    b = copy.deepcopy(bundle)

    # Validate structure
    for section in ("config", "seeds", "first_pass_config", "first_pass_seeds"):
        if section not in b:
            raise ValueError(f"bundle missing required section: {section}")

    for cfg_key, seeds_key, pass_label in [
        ("config", "seeds", "second-pass"),
        ("first_pass_config", "first_pass_seeds", "first-pass"),
    ]:
        cfg = b[cfg_key]
        seeds_dict = b[seeds_key]
        if keep_seed not in seeds_dict:
            raise ValueError(
                f"{pass_label}: seed {keep_seed} not present in {seeds_key}; "
                f"available: {sorted(seeds_dict.keys())}"
            )
        before = sorted(seeds_dict.keys())
        b[cfg_key] = dict(cfg)
        b[cfg_key]["seeds"] = [keep_seed]
        b[seeds_key] = {keep_seed: seeds_dict[keep_seed]}
        after = sorted(b[seeds_key].keys())
        n_models = len(b[seeds_key][keep_seed]["models"])
        print(f"  {pass_label}: seeds {before} → {after} "
              f"({len(before)*n_models} → {n_models} models)")

    # Tag the bundle so it's identifiable post-rebundle
    b["rebundle_provenance"] = {
        "rebundled_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "from_version": b.get("version"),
        "from_model_id": b.get("model_id"),
        "kept_seed": keep_seed,
        "note": ("Single-seed v2.7.1 rebundle. Empirically equivalent to the "
                 "3-seed v2.7.0 default per ensemble_size_sweep_first_pass_extended.tsv "
                 "across the 73-species IPA-gold panel."),
    }
    return b


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    DEFAULT_IN = Path(
        "/mnt/data/u12/intronIC_v2_repo/src/intronIC/data/default_pretrained.model.pkl"
    )
    ap.add_argument("--in", dest="inp", type=Path, default=DEFAULT_IN,
                    help="Input bundle path")
    ap.add_argument("--out", type=Path, default=None,
                    help="Output bundle path (default: <in>.v271.pkl alongside)")
    ap.add_argument("--keep", type=int, default=42,
                    help="Seed to keep (default 42)")
    ap.add_argument("--update-metadata", action="store_true",
                    help="Also write a sibling .metadata.json reflecting the rebundle")
    args = ap.parse_args()

    inp = args.inp.resolve()
    if not inp.is_file():
        sys.exit(f"input bundle not found: {inp}")
    out = args.out or inp.with_suffix(".v271.pkl")
    out = out.resolve()
    if out == inp:
        sys.exit("--out must differ from --in (refusing to overwrite the input)")

    in_size = inp.stat().st_size
    in_md5 = md5(inp)
    print(f"loading: {inp}")
    print(f"  size:  {in_size/1e6:.1f} MB")
    print(f"  md5:   {in_md5[:12]}")

    with open(inp, "rb") as f:
        bundle = pickle.load(f)
    new_bundle = trim_bundle(bundle, keep_seed=args.keep)

    print(f"writing: {out}")
    with open(out, "wb") as f:
        pickle.dump(new_bundle, f, protocol=pickle.HIGHEST_PROTOCOL)
    out_size = out.stat().st_size
    out_md5 = md5(out)
    print(f"  size:  {out_size/1e6:.1f} MB ({100*out_size/in_size:.0f}% of original)")
    print(f"  md5:   {out_md5[:12]}")

    if args.update_metadata:
        meta_in = inp.with_suffix(".metadata.json")
        meta_out = out.with_suffix(".metadata.json")
        if meta_in.is_file():
            metadata = json.load(open(meta_in))
            metadata.setdefault("rebundle_provenance", {})
            metadata["rebundle_provenance"].update(new_bundle["rebundle_provenance"])
            metadata["rebundle_provenance"]["new_md5"] = out_md5
            metadata["rebundle_provenance"]["new_size_bytes"] = out_size
            with open(meta_out, "w") as f:
                json.dump(metadata, f, indent=2)
            print(f"wrote metadata: {meta_out}")
        else:
            print(f"no sibling metadata.json at {meta_in}; skipped")

    print(f"\nDONE. Next: smoke-test new bundle on chr19 / HomSap before swap.")


if __name__ == "__main__":
    main()
