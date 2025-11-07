#!/usr/bin/env python3
"""Profile intronIC_refactored to identify performance bottlenecks."""

import cProfile
import pstats
import sys
from pathlib import Path

# Add intronIC_refactored to path
sys.path.insert(0, str(Path(__file__).parent))

# Import main after path setup
from intronIC_refactored.cli.main import main

if __name__ == '__main__':
    print("Starting profiled run...")

    # Run with profiling
    profiler = cProfile.Profile()
    profiler.enable()

    try:
        sys.exit(main())
    finally:
        profiler.disable()

        # Save detailed profile
        profiler.dump_stats('profile_output.pstats')

        # Print top 50 time-consuming functions
        print("\n" + "="*80)
        print("PROFILE RESULTS - Top 50 functions by cumulative time")
        print("="*80)
        stats = pstats.Stats(profiler)
        stats.strip_dirs()
        stats.sort_stats('cumulative')
        stats.print_stats(50)

        print("\n" + "="*80)
        print("PROFILE RESULTS - Top 50 functions by total time")
        print("="*80)
        stats.sort_stats('tottime')
        stats.print_stats(50)
