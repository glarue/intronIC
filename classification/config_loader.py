"""
Configuration file loader for SVM optimizer settings.

Allows specification of optimizer parameters and parameter grids via YAML
configuration files for easy testing without code modification.

Usage:
    from classification.config_loader import load_optimizer_config

    # Load from config file
    optimizer = load_optimizer_config('config/training_fast_test.yaml')

    # Then use as normal
    parameters = optimizer.optimize(u12_introns, u2_introns)
"""

from pathlib import Path
from typing import Dict, Any, Optional
import warnings

from classification.optimizer import SVMOptimizer


def load_optimizer_config(
    config_path: str,
    **override_kwargs
) -> SVMOptimizer:
    """
    Load SVMOptimizer configuration from YAML file.

    Args:
        config_path: Path to YAML configuration file
        **override_kwargs: Additional kwargs to override config values

    Returns:
        Configured SVMOptimizer instance

    Raises:
        FileNotFoundError: If config file doesn't exist
        ValueError: If config file is invalid

    Examples:
        >>> # Load fast test config
        >>> optimizer = load_optimizer_config('config/training_fast_test.yaml')

        >>> # Load with override
        >>> optimizer = load_optimizer_config(
        ...     'config/training_fast_test.yaml',
        ...     n_jobs=4  # Override n_jobs
        ... )
    """
    config_path = Path(config_path)

    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    try:
        import yaml
    except ImportError:
        # Fallback to simple key-value parsing if PyYAML not available
        config = _parse_simple_config(config_path)
    else:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)

    # Extract optimizer settings
    optimizer_config = config.get('optimizer', {})

    # Extract parameter grid
    param_grid = config.get('param_grid', None)

    # Apply overrides
    optimizer_config.update(override_kwargs)

    # Create optimizer
    return SVMOptimizer(
        param_grid_override=param_grid,
        **optimizer_config
    )


def _parse_simple_config(config_path: Path) -> Dict[str, Any]:
    """
    Simple YAML-like parser for when PyYAML is not available.

    Supports basic key: value syntax with sections. Does not support
    full YAML features like lists, nested structures, etc.

    Args:
        config_path: Path to config file

    Returns:
        Parsed configuration dictionary

    Note:
        This is a fallback. Install PyYAML for full feature support:
        pip install pyyaml
    """
    warnings.warn(
        "PyYAML not found. Using simplified config parser. "
        "Install PyYAML for full config support: pip install pyyaml",
        UserWarning
    )

    config = {'optimizer': {}, 'param_grid': {}}
    current_section = None

    with open(config_path, 'r') as f:
        for line in f:
            line = line.strip()

            # Skip comments and empty lines
            if not line or line.startswith('#'):
                continue

            # Section headers
            if line.endswith(':') and ' ' not in line:
                section_name = line[:-1]
                if section_name in ['optimizer', 'param_grid']:
                    current_section = section_name
                continue

            # Key-value pairs
            if ':' in line and current_section:
                key, value = line.split(':', 1)
                key = key.strip()
                value = value.strip()

                # Parse value
                parsed_value = _parse_value(value)
                config[current_section][key] = parsed_value

    return config


def _parse_value(value_str: str) -> Any:
    """
    Parse a value string to appropriate Python type.

    Args:
        value_str: String representation of value

    Returns:
        Parsed value (bool, int, float, list, or str)
    """
    value_str = value_str.strip()

    # Boolean
    if value_str.lower() in ('true', 'yes'):
        return True
    if value_str.lower() in ('false', 'no'):
        return False

    # Null
    if value_str.lower() in ('null', 'none'):
        return None

    # List (simplified: [1, 2, 3] or [true, false])
    if value_str.startswith('[') and value_str.endswith(']'):
        items = value_str[1:-1].split(',')
        return [_parse_value(item.strip()) for item in items if item.strip()]

    # Number
    try:
        if '.' in value_str or 'e' in value_str.lower():
            return float(value_str)
        else:
            return int(value_str)
    except ValueError:
        pass

    # String (remove quotes if present)
    if value_str.startswith(("'", '"')) and value_str.endswith(("'", '"')):
        return value_str[1:-1]

    return value_str


def get_config_summary(config_path: str) -> str:
    """
    Get a summary of a configuration file.

    Args:
        config_path: Path to YAML configuration file

    Returns:
        Human-readable summary string

    Example:
        >>> summary = get_config_summary('config/training_fast_test.yaml')
        >>> print(summary)
    """
    optimizer = load_optimizer_config(config_path)

    lines = []
    lines.append(f"Configuration: {config_path}")
    lines.append("")
    lines.append("Optimizer Settings:")
    lines.append(f"  Rounds: {optimizer.n_rounds}")
    lines.append(f"  Initial grid points: {optimizer.n_points_initial}")
    lines.append(f"  CV folds: {optimizer.cv_folds}")
    lines.append(f"  Parallel jobs: {optimizer.n_jobs}")
    lines.append(f"  Max iterations: {optimizer.max_iter}")
    lines.append("")

    if optimizer.param_grid_override:
        lines.append("Parameter Grid:")
        n_combinations = 1
        for key, values in optimizer.param_grid_override.items():
            n_combinations *= len(values)
            lines.append(f"  {key}: {values} ({len(values)} values)")

        lines.append("")
        lines.append(f"Total combinations per C value: {n_combinations}")

        # Estimate total CV fits
        total_fits = n_combinations * optimizer.n_points_initial * optimizer.cv_folds
        lines.append(f"Estimated CV fits: ~{total_fits:,}")

    return "\n".join(lines)


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        config_file = sys.argv[1]
        try:
            summary = get_config_summary(config_file)
            print(summary)
        except Exception as e:
            print(f"Error loading config: {e}")
            sys.exit(1)
    else:
        print("Usage: python -m classification.config_loader <config_file.yaml>")
        print("\nExample:")
        print("  python -m classification.config_loader config/training_fast_test.yaml")
        sys.exit(1)
