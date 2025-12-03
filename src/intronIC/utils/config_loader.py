"""
Configuration file loader for intronIC settings.

Supports unified YAML configuration with auto-loading from standard paths.

Usage:
    from intronIC.utils.config_loader import load_config, find_config

    # Auto-load from standard paths
    config = load_config()

    # Or load from explicit path
    config = load_config('config/config.yaml')

    # Find config file without loading
    path = find_config()

Legacy usage (deprecated):
    from intronIC.utils.config_loader import load_optimizer_config
    optimizer = load_optimizer_config('config/training_fast_test.yaml')
"""

from pathlib import Path
from typing import Dict, Any, Optional, Union
import warnings
import os

from intronIC.classification.optimizer import SVMOptimizer


def find_config(explicit_path: Optional[str] = None) -> Optional[Path]:
    """
    Find configuration file using priority order.

    Search order (highest to lowest priority):
    1. explicit_path (if provided)
    2. ./.intronIC.yaml (current directory)
    3. ~/.config/intronIC/config.yaml (XDG config dir)
    4. ~/.intronIC.yaml (user home)
    5. <install_dir>/config/config.yaml (built-in defaults)

    Args:
        explicit_path: Optional explicit path to config file

    Returns:
        Path to config file, or None if not found

    Examples:
        >>> # Find using auto-discovery
        >>> config_path = find_config()
        >>> if config_path:
        ...     print(f"Found: {config_path}")

        >>> # Check explicit path
        >>> config_path = find_config('my_config.yaml')
    """
    # Priority 1: Explicit path
    if explicit_path:
        path = Path(explicit_path)
        if path.exists():
            return path
        else:
            return None  # Let caller handle the error

    # Priority 2: Project directory (./.intronIC.yaml)
    project_config = Path('.intronIC.yaml')
    if project_config.exists():
        return project_config

    # Priority 3: XDG config directory (~/.config/intronIC/config.yaml)
    xdg_config_home = os.environ.get('XDG_CONFIG_HOME')
    if xdg_config_home:
        xdg_config = Path(xdg_config_home) / 'intronIC' / 'config.yaml'
    else:
        xdg_config = Path.home() / '.config' / 'intronIC' / 'config.yaml'
    if xdg_config.exists():
        return xdg_config

    # Priority 4: User home (~/.intronIC.yaml)
    home_config = Path.home() / '.intronIC.yaml'
    if home_config.exists():
        return home_config

    # Priority 5: Built-in defaults (<install_dir>/config/config.yaml)
    # Find install directory (repo root, 4 levels up from src/intronIC/classification/)
    install_dir = Path(__file__).parent.parent.parent.parent
    builtin_config = install_dir / 'config' / 'config.yaml'
    if builtin_config.exists():
        return builtin_config

    # No config found
    return None


def load_config(config_path: Optional[str] = None) -> Dict[str, Any]:
    """
    Load unified configuration file with auto-discovery.

    If config_path is not provided, searches standard locations in priority order:
    1. ./.intronIC.yaml (project directory)
    2. ~/.config/intronIC/config.yaml (XDG config dir)
    3. ~/.intronIC.yaml (user home)
    4. <install_dir>/config/config.yaml (built-in defaults)

    Args:
        config_path: Optional explicit path to config file

    Returns:
        Configuration dictionary with all sections (scoring, training, optimizer, etc.)

    Raises:
        FileNotFoundError: If explicit path provided but file doesn't exist
        ValueError: If config file is invalid YAML

    Examples:
        >>> # Auto-load from standard paths
        >>> config = load_config()
        >>> threshold = config['scoring']['threshold']

        >>> # Load from explicit path
        >>> config = load_config('config/profiles/quick.yaml')

        >>> # Load and access nested values
        >>> config = load_config()
        >>> n_models = config['training']['ensemble']['n_models']
    """
    # Find config file
    found_path = find_config(config_path)

    if found_path is None:
        if config_path:
            raise FileNotFoundError(f"Config file not found: {config_path}")
        else:
            # No config found anywhere - return empty dict (caller will use defaults)
            return {}

    # Load YAML
    try:
        import yaml
    except ImportError:
        raise ImportError(
            "PyYAML is required for config loading. "
            "Install with: pip install pyyaml"
        )

    with open(found_path, 'r') as f:
        config = yaml.safe_load(f)

    if config is None:
        config = {}

    return config


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

    # Extract C bounds (for passing to optimize() method later)
    c_bounds = config.get('c_bounds', {})

    # Extract training settings (ensemble configuration)
    training_config = config.get('training', {})

    # Apply overrides
    optimizer_config.update(override_kwargs)

    # Create optimizer
    optimizer = SVMOptimizer(
        param_grid_override=param_grid,
        **optimizer_config
    )

    # Attach c_bounds as attributes for later use
    # These are passed to optimize(), not to the constructor
    if c_bounds:
        eff_C_pos_range = c_bounds.get('eff_C_pos_range')
        if eff_C_pos_range:
            optimizer.eff_C_pos_range = tuple(eff_C_pos_range)  # Convert list to tuple
        eff_C_neg_max = c_bounds.get('eff_C_neg_max')
        if eff_C_neg_max is not None:
            optimizer.eff_C_neg_max = eff_C_neg_max

    # Attach training settings as attributes for later use
    # These control ensemble training (n_models, subsampling, etc.)
    if training_config:
        for key, value in training_config.items():
            setattr(optimizer, f'training_{key}', value)

    return optimizer


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
        print("Usage: python -m intronIC.utils.config_loader <config_file.yaml>")
        print("\nExample:")
        print("  python -m intronIC.utils.config_loader config/training_fast_test.yaml")
        sys.exit(1)
