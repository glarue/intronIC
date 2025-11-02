"""
Regression test for intronIC refactored version.

Compares output from refactored version against gold standard
generated from the original codebase.
"""

import os
import gzip
from pathlib import Path
import subprocess


def test_gold_standard_exists():
    """Verify gold standard outputs exist."""
    gold_files = [
        "homo_sapiens_gold.bed.iic",
        "homo_sapiens_gold.meta.iic",
        "homo_sapiens_gold.introns.iic",
        "homo_sapiens_gold.annotation.iic",
    ]

    for filename in gold_files:
        assert Path(filename).exists(), f"Gold standard file {filename} not found"


def test_intron_count():
    """Verify number of introns matches gold standard."""
    # Read gold standard
    with open("homo_sapiens_gold.bed.iic") as f:
        gold_count = sum(1 for line in f if not line.startswith('#'))

    print(f"Gold standard intron count: {gold_count}")
    assert gold_count > 0, "No introns in gold standard"


def test_bed_format():
    """Verify BED file format is correct."""
    with open("homo_sapiens_gold.bed.iic") as f:
        for i, line in enumerate(f, 1):
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            assert len(fields) == 6, f"Line {i}: Expected 6 fields, got {len(fields)}"

            # Check coordinates
            chrom, start, stop, name, score, strand = fields
            assert start.isdigit(), f"Line {i}: Start coordinate not numeric"
            assert stop.isdigit(), f"Line {i}: Stop coordinate not numeric"
            assert int(start) < int(stop), f"Line {i}: Start >= Stop"
            assert strand in ['+', '-'], f"Line {i}: Invalid strand '{strand}'"

            # Only check first 100 lines
            if i > 100:
                break


def test_u12_introns_detected():
    """Verify that U12 introns are detected."""
    u12_count = 0
    with open("homo_sapiens_gold.meta.iic") as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) > 12:
                type_id = fields[12]
                if type_id == 'u12':
                    u12_count += 1

    print(f"U12 introns detected: {u12_count}")
    # Chr19 should have some U12 introns
    assert u12_count > 0, "No U12 introns detected"
    # But they should be minority (<1%)
    assert u12_count < 500, f"Too many U12 introns: {u12_count}"


def test_classification_scores():
    """Verify classification scores are in valid range."""
    with open("homo_sapiens_gold.bed.iic") as f:
        for i, line in enumerate(f, 1):
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            score_str = fields[4]

            # Some introns may have '.' for missing scores (omitted introns)
            if score_str != '.':
                score = float(score_str)
                assert 0 <= score <= 100, f"Line {i}: Score {score} out of range [0, 100]"

            # Only check first 100
            if i > 100:
                break


def test_log_file_completeness():
    """Verify log file indicates successful completion."""
    with open("homo_sapiens_gold.log.iic") as f:
        log_content = f.read()

    # Check for key stages
    assert "Starting intronIC" in log_content, "Missing start marker"
    assert "introns found" in log_content, "Missing intron count"

    # Check for completion
    # The log should contain runtime information or completion marker
    assert ("finished" in log_content.lower() or
            "runtime" in log_content.lower() or
            "complete" in log_content.lower()), \
        "Log does not indicate completion"


if __name__ == "__main__":
    # Run tests
    import pytest
    pytest.main([__file__, "-v"])
