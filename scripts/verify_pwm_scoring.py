#!/usr/bin/env python3
"""
Manual PWM scoring verification — diagnostic/audit tool.

Loads score_info output from an intronIC run, looks up per-position
frequencies directly from the PWM JSON, computes log2 ratios, and
compares to the scorer's raw scores.

This catches indexing bugs, log-base errors, and extraction/scoring
window mismatches.

Note: Automated regression tests in tests/unit/test_scoring/test_pwm_regression.py
cover the same invariants using fixed test PWMs. This script is for ad-hoc
verification against the installed PWMs and actual intronIC output.

Usage:
    # Run intronIC test first to generate output:
    intronIC test -o /tmp/intronIC_verify_test -p 4

    # Then verify:
    python scripts/verify_pwm_scoring.py \
        --score-info /tmp/intronIC_verify_test/test_chr19.score_info.iic \
        --meta /tmp/intronIC_verify_test/test_chr19.meta.iic \
        --pwm-json src/intronIC/data/intronIC_scoring_PWMs.json \
        --min-svm 90 --n-introns 5

See docs/pwm_scoring_verification.md for pitfalls and interpretation.
"""

import argparse
import csv
import json
import math
import sys
from pathlib import Path


def json_key(pos_int: int) -> str:
    """Convert integer position to JSON key format.

    Negative → "-3", zero → "0", positive → "+1".
    """
    if pos_int > 0:
        return f"+{pos_int}"
    return str(pos_int)


def load_pwm_sets(pwm_json_path: str):
    """Load PWM JSON and return (raw_matrices, loader_pwm_sets).

    Returns both the raw JSON matrices (for documentation/spot-checks)
    and the loaded PWM objects (for scoring).
    """
    with open(pwm_json_path) as f:
        raw = json.load(f)
    raw_matrices = raw["matrix_groups"][0]["matrices"]

    from intronIC.scoring.pwm import PWMLoader

    pwm_sets = PWMLoader.load_from_file(Path(pwm_json_path))
    return raw_matrices, pwm_sets


def manual_log2_ratio(seq, u12_pwm, u2_pwm, seq_start_position):
    """Score a sequence using the loaded PWM objects and return log2(U12/U2).

    This replicates exactly what IntronScorer._calculate_log_ratio does:
    log2(u12_pwm.score_sequence(seq, ...) / u2_pwm.score_sequence(seq, ...))
    """
    u12_score = u12_pwm.score_sequence(seq, seq_start_position=seq_start_position)
    u2_score = u2_pwm.score_sequence(seq, seq_start_position=seq_start_position)
    return math.log2(u12_score / u2_score)


def per_position_breakdown(seq, raw_matrices, u12_name, u2_name, start_index):
    """Per-position JSON lookup for human-readable output.

    Uses the raw JSON matrices (not the loader) so we can verify the
    loader's conversion independently.

    Args:
        start_index: The biological position of seq[0]. For 5'SS/3'SS this
            is the scoring window start (e.g. -3 or -6). For BPS, pass the
            first biological position from the JSON keys (e.g. -9), NOT the
            JSON start_index field (which is 0 and refers to the scorer's
            internal seq_start_position, not the biological labeling).
    """
    u12_mat = raw_matrices[u12_name]["matrix"]
    u2_mat = raw_matrices[u2_name]["matrix"]
    lines = []
    total = 0.0
    missing = 0
    for i, base in enumerate(seq):
        logical = start_index + i
        key = json_key(logical)
        base_upper = base.upper()
        u12_entry = u12_mat.get(key)
        u2_entry = u2_mat.get(key)
        if u12_entry is None or u2_entry is None:
            lines.append(
                f"    pos {logical:+3d} (key={key}): {base_upper}  ** MISSING **"
            )
            missing += 1
            continue
        u12_freq = u12_entry.get(base_upper, 0.0001)
        u2_freq = u2_entry.get(base_upper, 0.0001)
        lr = math.log2(u12_freq / u2_freq)
        total += lr
        lines.append(
            f"    pos {logical:+3d}: {base_upper}  "
            f"u12={u12_freq:.6f}  u2={u2_freq:.6f}  LR={lr:+.4f}"
        )
    return total, lines, missing


def position_0_check(raw_matrices, u12_name, base, expected_base, label):
    """Check that position 0 has the expected base and high frequency."""
    mat = raw_matrices[u12_name]["matrix"]
    entry = mat.get("0")
    if entry is None:
        return f"  {label} pos 0: MISSING from JSON"
    freq = entry.get(base.upper(), 0)
    if expected_base:
        ok = base.upper() == expected_base and freq > 0.90
        return (
            f"  {label} pos 0: {base} (expected {expected_base}) "
            f"u12_freq={freq:.4f}  → {'PASS' if ok else 'FAIL'}"
        )
    else:
        ok = 0.05 < freq < 0.95
        return (
            f"  {label} pos 0: {base} (mixed)  "
            f"u12_freq={freq:.4f}  → {'PASS' if ok else 'CHECK'}"
        )


def verify(args):
    raw_matrices, pwm_sets = load_pwm_sets(args.pwm_json)

    # Load type_id from meta.iic
    type_map = {}
    with open(args.meta) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            type_map[row["name"]] = row["dnts"]

    # Load introns from score_info
    target_rows = []
    with open(args.score_info) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            svm = float(row["svm_score"])
            if svm >= args.min_svm:
                target_rows.append(row)
                if len(target_rows) >= args.n_introns:
                    break

    if not target_rows:
        print(f"No introns found with SVM >= {args.min_svm}")
        sys.exit(1)

    print(f"Verifying {len(target_rows)} introns (SVM >= {args.min_svm})")
    print(f"PWM JSON: {args.pwm_json}")
    print()

    # Scoring windows (must match IntronScorer defaults / config)
    # 5'SS: scorer uses (-3, 9) = 12 positions; extraction uses (-3, 10) = 13
    # 3'SS: scorer uses (-6, 4) = 10 positions; extraction uses (-14, 4) = 18
    # BPS: scorer uses seq_start_position=0, full 12-char sequence
    FIVE_SCORED_LEN = 12   # five_coords = (-3, 9)
    THREE_SCORED_LEN = 10  # three_coords = (-6, 4)
    FIVE_START = -3
    THREE_START = -6

    all_ok = True

    for row in target_rows:
        name = row["name"]
        type_id = type_map.get(name, "unknown")

        if type_id == "GT-AG":
            pfx = "gtag"
        elif type_id == "AT-AC":
            pfx = "atac"
        elif type_id == "GC-AG":
            pfx = "gtag"  # GC-AG uses GT-AG PWMs
        else:
            print(f"Skipping {name} (type {type_id})")
            continue

        u12_five_name = f"u12_{pfx}_five"
        u2_five_name = f"u2_{pfx}_five"
        u12_bp_name = f"u12_{pfx}_bp"
        u2_bp_name = f"u2_{pfx}_bp"
        u12_three_name = f"u12_{pfx}_three"
        u2_three_name = f"u2_{pfx}_three"

        u12_five_pwm = pwm_sets["five"].matrices[("u12", pfx)]
        u2_five_pwm = pwm_sets["five"].matrices[("u2", pfx)]
        u12_bp_pwm = pwm_sets["bp"].matrices[("u12", pfx)]
        u2_bp_pwm = pwm_sets["bp"].matrices[("u2", pfx)]
        u12_three_pwm = pwm_sets["three"].matrices[("u12", pfx)]
        u2_three_pwm = pwm_sets["three"].matrices[("u2", pfx)]

        five_full = row["5'_seq"]
        three_full = row["3'_seq"]
        bp_full = row["bp_seq"]

        five_raw_expected = float(row["5'_raw"])
        three_raw_expected = float(row["3'_raw"])
        bp_raw_expected = float(row["bp_raw"])

        print("=" * 75)
        print(f"  {name[:60]}  [{type_id}]")
        print("=" * 75)

        # ── 5'SS ──
        five_scored = five_full[:FIVE_SCORED_LEN]
        five_lr = manual_log2_ratio(five_scored, u12_five_pwm, u2_five_pwm, FIVE_START)
        diff5 = abs(five_lr - five_raw_expected)
        ok5 = diff5 < 0.001

        print(f"\n  5'SS scored: {five_scored}  ({FIVE_SCORED_LEN} of {len(five_full)} chars)")
        if args.verbose:
            _, lines, _ = per_position_breakdown(
                five_scored, raw_matrices, u12_five_name, u2_five_name, FIVE_START
            )
            for line in lines:
                print(line)
        print(
            f"  → log2 LR: {five_lr:.6f}  |  score_info: {five_raw_expected:.6f}  "
            f"|  diff: {diff5:.6f}  → {'MATCH' if ok5 else 'MISMATCH'}"
        )

        # ── BPS ──
        bp_lr = manual_log2_ratio(bp_full, u12_bp_pwm, u2_bp_pwm, 0)
        diff_bp = abs(bp_lr - bp_raw_expected)
        ok_bp = diff_bp < 0.001

        # BPS biological positions: derive from JSON keys (e.g. -9 to +2)
        bp_json_keys_sorted = sorted(
            raw_matrices[u12_bp_name]["matrix"].keys(), key=lambda k: int(k)
        )
        bp_bio_start = int(bp_json_keys_sorted[0])
        print(f"\n  BPS scored: {bp_full}  (bio positions {bp_bio_start} to {int(bp_json_keys_sorted[-1])})")
        if args.verbose:
            _, lines, _ = per_position_breakdown(
                bp_full, raw_matrices, u12_bp_name, u2_bp_name, bp_bio_start
            )
            for line in lines:
                print(line)
        print(
            f"  → log2 LR: {bp_lr:.6f}  |  score_info: {bp_raw_expected:.6f}  "
            f"|  diff: {diff_bp:.6f}  → {'MATCH' if ok_bp else 'MISMATCH'}"
        )

        # ── 3'SS ──
        three_scored = three_full[-THREE_SCORED_LEN:]
        three_lr = manual_log2_ratio(
            three_scored, u12_three_pwm, u2_three_pwm, THREE_START
        )
        diff3 = abs(three_lr - three_raw_expected)
        ok3 = diff3 < 0.001

        print(
            f"\n  3'SS scored: {three_scored}  "
            f"(last {THREE_SCORED_LEN} of {len(three_full)} chars)"
        )
        if args.verbose:
            _, lines, _ = per_position_breakdown(
                three_scored, raw_matrices, u12_three_name, u2_three_name, THREE_START
            )
            for line in lines:
                print(line)
        print(
            f"  → log2 LR: {three_lr:.6f}  |  score_info: {three_raw_expected:.6f}  "
            f"|  diff: {diff3:.6f}  → {'MATCH' if ok3 else 'MISMATCH'}"
        )

        # ── Position 0 spot checks ──
        print("\n  Position 0 checks:")
        # 5'SS pos 0: index abs(FIVE_START) into five_scored
        p0_5_idx = abs(FIVE_START)
        p0_5 = five_scored[p0_5_idx] if p0_5_idx < len(five_scored) else "?"
        expected_5 = "G" if pfx == "gtag" else "A"
        print(position_0_check(raw_matrices, u12_five_name, p0_5, expected_5, "5'SS"))

        # 3'SS pos 0: index abs(THREE_START) into three_scored
        p0_3_idx = abs(THREE_START)
        p0_3 = three_scored[p0_3_idx] if p0_3_idx < len(three_scored) else "?"
        print(position_0_check(raw_matrices, u12_three_name, p0_3, None, "3'SS"))

        # BPS pos 0: array index = (0 - internal_start), but for BPS with
        # JSON keys -9 to +2, position 0 = branch A = array index 9
        bp_json_keys = sorted(
            raw_matrices[u12_bp_name]["matrix"].keys(), key=lambda k: int(k)
        )
        bp_bio_positions = [int(k) for k in bp_json_keys]
        if 0 in bp_bio_positions:
            p0_bp_idx = bp_bio_positions.index(0)
            p0_bp = bp_full[p0_bp_idx] if p0_bp_idx < len(bp_full) else "?"
            print(position_0_check(raw_matrices, u12_bp_name, p0_bp, "A", "BPS "))

        intron_ok = ok5 and ok_bp and ok3
        if not intron_ok:
            all_ok = False
        print()

    # ── Summary ──
    print("=" * 75)
    if all_ok:
        print("ALL REGIONS MATCH — PWM SCORING VERIFIED CORRECT")
    else:
        print("DISCREPANCY FOUND — INVESTIGATE")
        sys.exit(1)
    print("=" * 75)


def main():
    parser = argparse.ArgumentParser(
        description="Verify PWM scoring by comparing manual JSON lookups to scorer output."
    )
    parser.add_argument(
        "--score-info",
        required=True,
        help="Path to .score_info.iic from an intronIC run",
    )
    parser.add_argument(
        "--meta",
        required=True,
        help="Path to .meta.iic (for type_id / dnts column)",
    )
    parser.add_argument(
        "--pwm-json",
        default="src/intronIC/data/intronIC_scoring_PWMs.json",
        help="Path to PWM JSON file (default: installed PWMs)",
    )
    parser.add_argument(
        "--min-svm",
        type=float,
        default=90.0,
        help="Minimum SVM score to select introns for verification (default: 90)",
    )
    parser.add_argument(
        "--n-introns",
        type=int,
        default=5,
        help="Number of introns to verify (default: 5)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Show per-position frequency breakdown",
    )
    args = parser.parse_args()
    verify(args)


if __name__ == "__main__":
    main()
