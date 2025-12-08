"""
TEMPORARY: Memory profiling utilities for debugging memory usage.
TODO: Remove this entire file after memory analysis is complete.

Author: Memory profiling investigation
Date: 2025-11-26
"""

import psutil
import os


def print_memory_usage(label: str) -> float:
    """
    Print current memory usage with a label.

    Args:
        label: Description of the checkpoint

    Returns:
        Memory usage in MB
    """
    process = psutil.Process(os.getpid())
    mem_mb = process.memory_info().rss / 1024 / 1024
    print(f"[MEMORY CHECKPOINT] {label}: {mem_mb:.1f} MB")
    return mem_mb
