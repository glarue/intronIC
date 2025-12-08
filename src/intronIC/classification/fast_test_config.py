"""
Fast testing configurations for BothEndsStrong feature development.

Provides pre-configured optimizer settings and parameter grids for quick
validation without modifying core codebase.

Usage:
    from intronIC.classification.fast_test_config import get_fast_optimizer

    # Phase 1: Ultra-fast smoke test (2-5 minutes)
    optimizer = get_fast_optimizer(mode='smoke')

    # Phase 2: Quick validation (10-20 minutes)
    optimizer = get_fast_optimizer(mode='quick')

    # Phase 3: Moderate testing (30-60 minutes)
    optimizer = get_fast_optimizer(mode='moderate')

    # Then use as normal:
    parameters = optimizer.optimize(u12_introns, u2_introns)
"""

from typing import Dict, Any
from intronIC.classification.optimizer import SVMOptimizer


# Pre-configured parameter grids for different testing phases
PARAM_GRIDS = {
    # Phase 1: Ultra-fast smoke test
    # Tests basic functionality with minimal grid
    # Expected time: 2-5 minutes on Chr19
    # Combinations: 1 include_max × 1 dual × 1 intercept × 1 method = 1
    'smoke': {
        'estimator__augment__include_max': [False],     # Min features only
        'estimator__svc__dual': [False],                # Primal formulation (typical best)
        'estimator__svc__intercept_scaling': [1000.0],  # Typical good value
        'method': ['sigmoid']                           # Faster than isotonic
    },

    # Phase 2: Quick validation
    # Tests min-only vs min+max features
    # Expected time: 10-20 minutes on Chr19
    # Combinations: 2 include_max × 1 dual × 1 intercept × 1 method = 2
    'quick': {
        'estimator__augment__include_max': [False, True],  # Test both
        'estimator__svc__dual': [False],                   # Fixed
        'estimator__svc__intercept_scaling': [1000.0],     # Fixed
        'method': ['sigmoid']                              # Fixed
    },

    # Phase 3: Moderate testing
    # Tests min+max features with calibration methods
    # Expected time: 30-60 minutes on Chr19
    # Combinations: 2 include_max × 1 dual × 1 intercept × 2 method = 4
    'moderate': {
        'estimator__augment__include_max': [False, True],  # Test both
        'estimator__svc__dual': [False],                   # Fixed
        'estimator__svc__intercept_scaling': [1000.0],     # Fixed
        'method': ['sigmoid', 'isotonic']                  # Test both
    },

    # Phase 4: Near-full testing
    # Tests all parameter combinations
    # Expected time: 1-2 hours on Chr19
    # Combinations: 2 include_max × 2 dual × 3 intercept × 2 method = 24
    'near_full': {
        'estimator__augment__include_max': [False, True],          # Test both
        'estimator__svc__dual': [False, True],                     # Test both
        'estimator__svc__intercept_scaling': [10.0, 100.0, 1000.0],  # Full range
        'method': ['sigmoid', 'isotonic']                          # Both methods
    },
}


# Pre-configured optimizer settings for different testing phases
OPTIMIZER_CONFIGS = {
    'smoke': {
        'n_rounds': 1,              # Single round
        'n_points_initial': 7,      # Fewer C values
        'n_points_refine': 20,      # Not used with 1 round
        'cv_folds': 3,              # 3-fold CV
        'n_jobs': -1,               # Use all cores
        'verbose': True,
        'max_iter': 10000           # May see convergence warnings
    },

    'quick': {
        'n_rounds': 2,              # 2 rounds
        'n_points_initial': 9,      # Modest C grid
        'n_points_refine': 30,      # Modest refinement
        'cv_folds': 3,              # 3-fold CV
        'n_jobs': -1,
        'verbose': True,
        'max_iter': 20000
    },

    'moderate': {
        'n_rounds': 2,              # 2 rounds
        'n_points_initial': 11,     # Larger C grid
        'n_points_refine': 50,      # More refinement
        'cv_folds': 5,              # Full 5-fold CV
        'n_jobs': -1,
        'verbose': True,
        'max_iter': 50000
    },

    'near_full': {
        'n_rounds': 3,              # Full 3 rounds
        'n_points_initial': 13,     # Full C grid
        'n_points_refine': 100,     # Full refinement
        'cv_folds': 5,              # 5-fold CV
        'n_jobs': -1,
        'verbose': True,
        'max_iter': 100000
    },
}


def get_fast_optimizer(mode: str = 'quick', **kwargs) -> SVMOptimizer:
    """
    Create an SVMOptimizer configured for fast testing.

    Args:
        mode: Testing mode - 'smoke', 'quick', 'moderate', or 'near_full'
        **kwargs: Additional keyword arguments override config values

    Returns:
        Configured SVMOptimizer with reduced parameter grid

    Raises:
        ValueError: If mode is not recognized

    Examples:
        # Ultra-fast smoke test
        >>> optimizer = get_fast_optimizer('smoke')

        # Quick test with custom n_jobs
        >>> optimizer = get_fast_optimizer('quick', n_jobs=4)

        # Moderate test with verbose disabled
        >>> optimizer = get_fast_optimizer('moderate', verbose=False)
    """
    if mode not in PARAM_GRIDS:
        raise ValueError(
            f"Unknown mode '{mode}'. "
            f"Choose from: {', '.join(PARAM_GRIDS.keys())}"
        )

    # Get base configuration
    config = OPTIMIZER_CONFIGS[mode].copy()
    param_grid = PARAM_GRIDS[mode]

    # Apply any user overrides
    config.update(kwargs)

    # Create optimizer with custom param grid
    return SVMOptimizer(
        param_grid_override=param_grid,
        **config
    )


def get_param_grid(mode: str) -> Dict[str, Any]:
    """
    Get just the parameter grid for a testing mode.

    Useful if you want to create a custom optimizer configuration
    but use a pre-defined parameter grid.

    Args:
        mode: Testing mode - 'smoke', 'quick', 'moderate', or 'near_full'

    Returns:
        Parameter grid dictionary

    Example:
        >>> param_grid = get_param_grid('quick')
        >>> optimizer = SVMOptimizer(
        ...     n_rounds=1,
        ...     cv_folds=3,
        ...     param_grid_override=param_grid
        ... )
    """
    if mode not in PARAM_GRIDS:
        raise ValueError(
            f"Unknown mode '{mode}'. "
            f"Choose from: {', '.join(PARAM_GRIDS.keys())}"
        )

    return PARAM_GRIDS[mode].copy()


def print_config_summary(mode: str = None):
    """
    Print a summary of available testing configurations.

    Args:
        mode: Optional specific mode to show details for.
              If None, shows summary of all modes.

    Example:
        >>> from classification.fast_test_config import print_config_summary
        >>> print_config_summary()  # Show all modes
        >>> print_config_summary('quick')  # Show details for 'quick' mode
    """
    if mode is not None:
        if mode not in PARAM_GRIDS:
            print(f"Unknown mode '{mode}'")
            print(f"Available modes: {', '.join(PARAM_GRIDS.keys())}")
            return

        print(f"\n{'='*80}")
        print(f"CONFIGURATION: {mode.upper()}")
        print(f"{'='*80}")

        config = OPTIMIZER_CONFIGS[mode]
        param_grid = PARAM_GRIDS[mode]

        print(f"\nOptimizer Settings:")
        for key, value in config.items():
            print(f"  {key}: {value}")

        print(f"\nParameter Grid:")
        # Calculate combinations
        n_combinations = 1
        for key, values in param_grid.items():
            n_combinations *= len(values)
            print(f"  {key}: {values} ({len(values)} values)")

        print(f"\nTotal combinations per C value: {n_combinations}")
        print(f"Expected C values tested: ~{config['n_points_initial']} (round 1)")
        total_cv_fits = n_combinations * config['n_points_initial'] * config['cv_folds']
        print(f"Estimated CV fits: ~{total_cv_fits:,}")

    else:
        print(f"\n{'='*80}")
        print("AVAILABLE FAST TEST CONFIGURATIONS")
        print(f"{'='*80}\n")

        for mode_name in PARAM_GRIDS.keys():
            config = OPTIMIZER_CONFIGS[mode_name]
            param_grid = PARAM_GRIDS[mode_name]

            n_combinations = 1
            for values in param_grid.values():
                n_combinations *= len(values)

            print(f"{mode_name.upper()}:")
            print(f"  Rounds: {config['n_rounds']}, CV folds: {config['cv_folds']}")
            print(f"  Param combinations: {n_combinations}")
            print(f"  Estimated time on Chr19: ", end="")

            if mode_name == 'smoke':
                print("2-5 minutes")
            elif mode_name == 'quick':
                print("10-20 minutes")
            elif mode_name == 'moderate':
                print("30-60 minutes")
            elif mode_name == 'near_full':
                print("1-2 hours")

            print()

        print("Usage:")
        print("  from classification.fast_test_config import get_fast_optimizer")
        print("  optimizer = get_fast_optimizer('quick')")
        print()


if __name__ == "__main__":
    # When run as script, print summary
    print_config_summary()
