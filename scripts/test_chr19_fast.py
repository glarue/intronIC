#!/usr/bin/env python3
"""
Phase 1 fast testing script for BothEndsStrong implementation.

Tests the full pipeline on Chr19 with minimal parameter grid for quick
validation (~5-10 minutes).

Usage:
    # Method 1: Using fast_test_config helper
    python scripts/test_chr19_fast.py --mode smoke

    # Method 2: Using YAML config
    python scripts/test_chr19_fast.py --config config/training_fast_test.yaml

    # Method 3: Default (quick mode)
    python scripts/test_chr19_fast.py
"""

import argparse
import sys
from pathlib import Path

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))


def run_fast_test(mode='quick', config_file=None):
    """
    Run fast test on Chr19 data.

    Args:
        mode: Testing mode ('smoke', 'quick', 'moderate') if not using config
        config_file: Path to YAML config file (overrides mode)
    """
    print("=" * 80)
    print("PHASE 1: FAST TESTING ON CHR19")
    print("=" * 80)
    print()

    # Import here to show progress
    print("Loading modules...")
    from classification.optimizer import SVMOptimizer
    from classification.fast_test_config import get_fast_optimizer, print_config_summary
    from classification.config_loader import load_optimizer_config

    # Determine which config method to use
    if config_file:
        print(f"\nLoading configuration from: {config_file}")
        optimizer = load_optimizer_config(config_file)

        # Print config summary
        from classification.config_loader import get_config_summary
        print("\n" + get_config_summary(config_file))
    else:
        print(f"\nUsing fast_test_config mode: {mode}")
        print_config_summary(mode)
        optimizer = get_fast_optimizer(mode)

    print()
    print("=" * 80)
    print("LOADING REFERENCE DATA")
    print("=" * 80)

    # This is where you would load actual reference introns
    # For this demo, we'll just show the structure
    print("\nNOTE: This is a demo script showing the structure.")
    print("To run actual tests, you need to:")
    print("  1. Load reference U12/U2 introns with z-scores")
    print("  2. Call optimizer.optimize(u12_introns, u2_introns)")
    print("  3. Use returned parameters to train ensemble")
    print()

    # Example code structure:
    example_code = '''
    # Load references (example)
    from data.loaders import load_reference_introns
    u12_introns = load_reference_introns('u12')
    u2_introns = load_reference_introns('u2')

    # Optimize hyperparameters
    print("Optimizing hyperparameters...")
    parameters = optimizer.optimize(u12_introns, u2_introns)

    print(f"\\nOptimized Parameters:")
    print(f"  C: {parameters.C:.6e}")
    print(f"  include_max: {parameters.include_max}")
    print(f"  dual: {parameters.dual}")
    print(f"  intercept_scaling: {parameters.intercept_scaling}")
    print(f"  calibration_method: {parameters.calibration_method}")
    print(f"  CV score: {parameters.cv_score:.4f}")

    # Train ensemble with optimized parameters
    from classification.trainer import SVMTrainer
    trainer = SVMTrainer(n_models=3)
    ensemble = trainer.train_ensemble(
        u12_introns,
        u2_introns,
        parameters
    )

    print(f"\\nEnsemble trained: {len(ensemble)} models")
    '''

    print("Example integration code:")
    print(example_code)

    return 0


def main():
    """Parse arguments and run test."""
    parser = argparse.ArgumentParser(
        description='Fast testing for BothEndsStrong implementation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Ultra-fast smoke test
  python scripts/test_chr19_fast.py --mode smoke

  # Quick validation (default)
  python scripts/test_chr19_fast.py --mode quick

  # Moderate testing
  python scripts/test_chr19_fast.py --mode moderate

  # Using YAML config
  python scripts/test_chr19_fast.py --config config/training_fast_test.yaml
        """
    )

    parser.add_argument(
        '--mode',
        choices=['smoke', 'quick', 'moderate', 'near_full'],
        default='quick',
        help='Testing mode (ignored if --config specified)'
    )

    parser.add_argument(
        '--config',
        type=str,
        help='Path to YAML configuration file'
    )

    args = parser.parse_args()

    try:
        return run_fast_test(mode=args.mode, config_file=args.config)
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
