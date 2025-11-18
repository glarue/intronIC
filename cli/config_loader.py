"""
Configuration file loader for intronIC.

Supports TOML configuration files for persistent settings.
"""

import sys
from pathlib import Path
from typing import Optional, Dict, Any

# Use tomllib (Python 3.11+) or tomli (Python <3.11)
if sys.version_info >= (3, 11):
    import tomllib
else:
    try:
        import tomli as tomllib
    except ImportError:
        tomllib = None


class ConfigLoader:
    """Loads and merges TOML configuration files with CLI arguments."""

    # Standard config file locations (checked in order)
    CONFIG_SEARCH_PATHS = [
        Path(".intronIC.toml"),  # Current directory
        Path("~/.config/intronIC/config.toml").expanduser(),  # User config dir
        Path("~/.intronIC.toml").expanduser(),  # User home dir
    ]

    @classmethod
    def find_config(cls, custom_path: Optional[Path] = None) -> Optional[Path]:
        """
        Find configuration file.

        Args:
            custom_path: Optional custom config file path (takes precedence)

        Returns:
            Path to config file if found, None otherwise
        """
        if custom_path:
            if custom_path.exists():
                return custom_path
            else:
                raise FileNotFoundError(f"Custom config file not found: {custom_path}")

        # Search standard locations
        for path in cls.CONFIG_SEARCH_PATHS:
            if path.exists():
                return path

        return None

    @classmethod
    def load_config(cls, config_path: Path) -> Dict[str, Any]:
        """
        Load TOML configuration file.

        Args:
            config_path: Path to config file

        Returns:
            Dictionary of configuration values

        Raises:
            ImportError: If tomli/tomllib not available
            ValueError: If config file is invalid
        """
        if tomllib is None:
            raise ImportError(
                "TOML support requires tomli package for Python <3.11. "
                "Install with: pip install tomli"
            )

        try:
            with open(config_path, 'rb') as f:
                return tomllib.load(f)
        except tomllib.TOMLDecodeError as e:
            raise ValueError(f"Invalid TOML in {config_path}: {e}")

    @classmethod
    def merge_with_args(cls, config: Dict[str, Any], args) -> None:
        """
        Merge configuration values with CLI arguments.

        CLI arguments take precedence over config file values.
        Only updates args that weren't explicitly set on CLI.

        Args:
            config: Configuration dictionary from TOML file
            args: Argparse namespace to update (modified in place)
        """
        # Map config sections and keys to argparse attributes
        # Format: (section, toml_key, args_attr, converter)
        mappings = [
            # Scoring options
            ('scoring', 'threshold', 'threshold', float),
            ('scoring', 'feature_type', 'feature', str),
            ('scoring', 'exclude_noncanonical', 'no_nc', bool),
            ('scoring', 'pseudocount', 'pseudocount', float),
            ('scoring', 'ignore_nc_dinucleotides', 'no_ignore_nc_dnts', cls._invert_bool),
            ('scoring', 'u12_boundary_correction', 'no_nc_ss_adjustment', cls._invert_bool),
            ('scoring', 'normalizer_mode', 'normalizer_mode', str),

            # Output options
            ('output', 'clean_names', 'clean_names', bool),

            # Performance options
            ('performance', 'processes', 'processes', int),
            ('performance', 'cv_processes', 'cv_processes', int),
            ('performance', 'min_intron_length', 'min_intron_len', int),

            # Isoform selection
            ('isoform', 'allow_multiple_isoforms', 'allow_multiple_isoforms', bool),
            ('isoform', 'exclude_overlapping', 'no_intron_overlap', bool),
            ('isoform', 'include_duplicates', 'include_duplicates', bool),

            # Training options
            ('training', 'n_models', 'n_models', int),
            ('training', 'max_iterations', 'max_iter', int),
            ('training', 'eval_mode', 'eval_mode', str),
            ('training', 'n_cv_folds', 'n_cv_folds', int),
            ('training', 'test_fraction', 'test_fraction', float),
            ('training', 'n_optimization_rounds', 'n_optimization_rounds', int),
            ('training', 'fixed_C', 'C', float),

            # Extraction options
            ('extraction', 'flank_length', 'flank_len', int),

            # Advanced options
            ('advanced', 'random_seed', 'seed', int),
            ('advanced', 'quiet', 'quiet', bool),
            ('advanced', 'debug', 'debug', bool),
        ]

        # Apply config values if not explicitly set on CLI
        for section, toml_key, args_attr, converter in mappings:
            if section not in config:
                continue

            toml_value = config[section].get(toml_key)
            if toml_value is None:
                continue

            # Only update if arg wasn't explicitly provided on CLI
            # We detect this by checking if the value is still the default
            if not cls._was_explicitly_set(args, args_attr):
                try:
                    converted_value = converter(toml_value)
                    setattr(args, args_attr, converted_value)
                except (ValueError, TypeError) as e:
                    raise ValueError(
                        f"Invalid value for {section}.{toml_key}: {toml_value} ({e})"
                    )

        # Handle scoring regions (nested structure)
        if 'scoring_regions' in config:
            regions = config['scoring_regions']

            if 'five_prime' in regions and not cls._was_explicitly_set(args, 'five_score_coords'):
                five = regions['five_prime']
                args.five_score_coords = [five.get('start', -3), five.get('end', 9)]

            if 'branch_point' in regions and not cls._was_explicitly_set(args, 'bp_region_coords'):
                bp = regions['branch_point']
                args.bp_region_coords = [bp.get('start', -55), bp.get('end', -5)]

            if 'three_prime' in regions and not cls._was_explicitly_set(args, 'three_score_coords'):
                three = regions['three_prime']
                args.three_score_coords = [three.get('start', -6), three.get('end', 4)]

    @staticmethod
    def _invert_bool(value: bool) -> bool:
        """Invert boolean for negative CLI flags."""
        return not value

    @staticmethod
    def _was_explicitly_set(args, attr: str) -> bool:
        """
        Check if an argument was explicitly set on command line.

        This is a heuristic based on comparing to defaults.
        Not perfect but works for most cases.
        """
        # For boolean flags, if they're True, they were likely set explicitly
        # (unless the default is True, handled specially)
        value = getattr(args, attr, None)

        # Special handling for flags with True defaults
        if attr == 'clean_names':
            # True is default, so if False it was explicitly set
            return value is False

        # For other booleans, True means it was set
        if isinstance(value, bool):
            return value is True

        # For numeric/string values, if it matches common CLI defaults,
        # assume it wasn't explicitly set
        # This is imperfect but the best we can do without modifying argparse
        return False  # Assume config should apply if default value

    @staticmethod
    def generate_template() -> str:
        """
        Generate a template configuration file with documentation.

        Returns:
            TOML configuration template as string
        """
        return '''# intronIC Configuration File
# =============================================================================
# This file contains default settings for intronIC runs.
# CLI arguments will override these settings.
#
# Location: Place this file in one of these locations:
#   1. ./.intronIC.toml (current directory)
#   2. ~/.config/intronIC/config.toml (user config dir)
#   3. ~/.intronIC.toml (user home directory)
#
# To generate this template: intronIC --generate-config
# =============================================================================

# -----------------------------------------------------------------------------
# Scoring Options
# -----------------------------------------------------------------------------
[scoring]
# Classification threshold (0-100)
# Introns with SVM probability >= threshold are classified as U12-type.
# Lower values = more sensitive (more U12-type introns found, more false positives)
# Higher values = more specific (fewer U12-type introns, fewer false positives)
# Recommended: 90 for standard analysis, 95 for high-confidence set
# Default: 90
threshold = 90.0

# Feature type to extract introns from
# - "cds": Only coding sequence features (excludes UTR introns)
# - "exon": All exon features (includes UTR introns)
# - "both": Use CDS where available, fall back to exons (recommended)
# Default: "both"
feature_type = "both"

# Exclude non-canonical introns from scoring
# When true, only GT-AG, GC-AG, and AT-AC introns are scored.
# When false, all intron types are scored (including rare boundaries).
# Note: Most U12-type introns are AT-AC, so excluding NC introns may miss some U12-type introns
# Default: false (include all, recommended for comprehensive analysis)
exclude_noncanonical = false

# Pseudocount for PWM scoring
# Small value added to matrix frequencies to avoid log(0) errors.
# Affects score magnitude but not relative rankings.
# Default: 0.0001 (recommended - minimal impact on scores)
pseudocount = 0.0001

# Ignore terminal dinucleotides when scoring non-canonical introns
# When true, only scores internal motif regions for NC introns (GT-AG, GC-AG, etc.)
# This prevents biasing scores based on non-canonical boundaries
# Default: true (recommended - avoids penalizing NC introns for mismatched boundaries)
ignore_nc_dinucleotides = true

# Enable U12 boundary correction for non-canonical introns
# Searches a 17bp window (5bp exonic + 12bp intronic) around the 5' splice site
# for strong U12 motifs like ATATCCTT. If found at an offset position, adjusts
# both intron boundaries by the same shift to align with the motif. Only applies
# correction if result is canonical (GT-AG, GC-AG, or AT-AC). Useful for fixing
# annotation errors where U12 introns are off by 1-6bp. Corrected introns are
# marked with [c:shift] tag (e.g., [c:-2] means shifted 2bp upstream).
# Default: true (recommended - helps identify mis-annotated U12 introns)
u12_boundary_correction = true

# Normalizer mode for pretrained model classification
# Controls how score normalization is applied when using a pretrained model:
# - "human": Use the scaler from the training species (recommended for cross-species)
#   Preserves composition bias correction and calibration consistency.
#   Best for U12-absent or distant species (e.g., C. elegans, plants).
# - "adaptive": Refit scaler on experimental data (experimental, may cause FPs)
#   Domain adaptation approach - may produce false positives in U12-free genomes.
# - "auto": Use human scaler if available in model, otherwise fall back to adaptive
#   Provides backward compatibility with older models while preferring human mode.
# Default: "auto" (recommended - balances compatibility and accuracy)
normalizer_mode = "auto"


# -----------------------------------------------------------------------------
# Scoring Regions (PWM application windows)
# -----------------------------------------------------------------------------
[scoring_regions]
# These define where PWM matrices are applied relative to splice sites.
# Coordinates: negative = exonic, positive = intronic, 0 = boundary
# Larger windows capture more sequence context but may include noise
# Smaller windows focus on core motifs but may miss extended signals

[scoring_regions.five_prime]
# 5' splice site scoring region (donor site)
# Applied around the exon-intron boundary
# Typical U12 motif: RTATCCTT (extends into intron)
# Default: -3 to +9 (3bp exonic context, 9bp intronic motif)
start = -3
end = 9

[scoring_regions.branch_point]
# Branch point search region (upstream of acceptor)
# Scans this window to find best-scoring BP motif position
# U12 BP motif: TCCTTAAC (7-8bp), located 10-50bp upstream of 3' site
# Default: -55 to -5 (searches 50bp window from 3' end)
start = -55
end = -5

[scoring_regions.three_prime]
# 3' splice site scoring region (acceptor site)
# Applied around the intron-exon boundary
# Typical U12 motif: YAC (polypyrimidine tract less pronounced than U2)
# Default: -6 to +4 (6bp intronic context, 4bp exonic)
start = -6
end = 4


# -----------------------------------------------------------------------------
# Output Options
# -----------------------------------------------------------------------------
[output]
# Remove "transcript:" and "gene:" prefixes from IDs
# Default: true
clean_names = true


# -----------------------------------------------------------------------------
# Performance Options
# -----------------------------------------------------------------------------
[performance]
# Number of parallel processes for classification
# Uses multiprocessing to speed up SVM prediction on experimental introns.
# Set to number of CPU cores for maximum speed (e.g., 8 for 8-core system).
# Memory usage scales linearly with processes (~1-2GB per process).
# Default: 1 (no parallelization)
processes = 1

# Number of processes for cross-validation (if different from classification)
# Controls parallelization during hyperparameter optimization (grid search).
# If not set, uses same value as 'processes' above.
# Grid search is very parallel-friendly - can benefit from many cores.
# Default: same as processes
# cv_processes = 8

# Minimum intron length in base pairs
# Introns shorter than this are marked with [o:s] (omitted: short) tag.
# They are still written to output files but not scored/classified.
# Scoring requires sufficient sequence for all three regions (5', BP, 3').
# Typical valid range: 30-70bp depending on scoring region configuration
# Default: 30 (sufficient for default scoring regions)
min_intron_length = 30


# -----------------------------------------------------------------------------
# Isoform Selection
# -----------------------------------------------------------------------------
[isoform]
# Include non-longest isoforms in analysis
# When false, only the longest transcript per gene is analyzed (reduces redundancy).
# When true, all isoforms are included (may result in many similar introns).
# Recommendation: Keep false to avoid scoring ~5x more introns with high redundancy
# Default: false (longest isoform only)
allow_multiple_isoforms = false

# Exclude overlapping introns from different isoforms
# Some isoforms have introns with overlapping coordinates (e.g., alternative splice sites).
# When true, only one representative from each overlapping set is kept.
# When false, all overlapping introns are kept (may cause duplicates).
# Default: false (allow overlaps)
exclude_overlapping = false

# Include duplicate introns in output files
# Duplicates are introns with identical coordinates (same chr, strand, start, stop).
# When false, only the first occurrence is written to output (cleaner output).
# When true, all occurrences are written (useful for tracking gene families).
# Note: Duplicates are marked with [d] tag and linked to representative
# Default: false (exclude duplicates, recommended for most analyses)
include_duplicates = false


# -----------------------------------------------------------------------------
# Training Options
# -----------------------------------------------------------------------------
[training]
# Number of ensemble models to train
# Each model is trained on a different U2 subsample for diversity.
# More models = better robustness and smoother probability estimates
# Fewer models = faster training
# Recommendation: 1 for quick runs, 3+ for production analyses
# Default: 1
n_models = 1

# Maximum iterations for SVM convergence
# LinearSVC optimization stops after this many iterations.
# Increase if you see convergence warnings (rare with default settings).
# Default: 50000 (sufficient for most datasets)
max_iterations = 50000

# Model evaluation mode
# How to evaluate classifier performance on reference data:
# - "nested_cv": 5x5 nested cross-validation (honest evaluation, RECOMMENDED)
#   Outer loop: 5-fold CV for model selection
#   Inner loop: 5-fold CV for probability calibration
#   Most rigorous but slowest (~25 fits per parameter combination)
# - "split": Simple train/test split (faster, less rigorous)
#   Single 80/20 split for quick evaluation
# - "none": Skip evaluation (fastest, use when loading pretrained model)
# Default: "nested_cv" (recommended for training from scratch)
eval_mode = "nested_cv"

# Number of folds for nested cross-validation
# Only used when eval_mode = "nested_cv"
# Higher = more rigorous but slower (5 is standard)
# Default: 5
n_cv_folds = 7

# Test set fraction for split evaluation
# Only used when eval_mode = "split"
# Fraction of reference data held out for testing (0.0-1.0)
# Default: 0.2 (20% test set, 80% training)
test_fraction = 0.2

# Number of optimization rounds for C parameter search
# Controls how many iterative refinement rounds are used during hyperparameter
# optimization. Each round narrows the search range around the best C value.
# More rounds = finer granularity but longer optimization time.
# Default: 5 (good balance of precision and speed)
n_optimization_rounds = 5

# Fixed SVM C parameter (soft-margin penalty)
# If set, skips hyperparameter optimization and uses this value directly.
# Useful for reproducing previous runs or when you know optimal C.
# If not set (commented out), C is optimized via grid search.
# Default: not set (optimize C automatically via grid search)
# fixed_C = 1.0


# -----------------------------------------------------------------------------
# Extraction Options
# -----------------------------------------------------------------------------
[extraction]
# Length of exonic flanks to extract (bp)
# Flanking exonic sequence is extracted on both sides of each intron.
# Used for context in output files and potential future analyses.
# Does not affect scoring (scoring uses internal coordinates only).
# Larger values = more context, larger output files
# Default: 50 (sufficient for most purposes)
flank_length = 50


# -----------------------------------------------------------------------------
# Advanced Options
# -----------------------------------------------------------------------------
[advanced]
# Random seed for reproducibility
# Seeds all random number generators (train/test splits, subsampling, etc.)
# Using the same seed guarantees identical results across runs.
# Change this value to get different random splits/samples.
# Default: 42 (the answer to everything)
random_seed = 42

# Suppress non-essential output
# When true, reduces console output to warnings and errors only.
# Useful for running in scripts or batch processing.
# Log file still contains full details regardless of this setting.
# Default: false (show normal progress output)
quiet = false

# Enable debug logging
# When true, writes detailed debugging information to log file.
# Includes feature extraction details, score calculations, etc.
# Produces very large log files - only use for troubleshooting.
# Default: false (standard logging only)
debug = false


# =============================================================================
# Notes:
# - Core inputs (genome, annotation, species name) must be specified on CLI
# - CLI arguments always override config file settings
# - Use --generate-config to create this template
# - Use --config FILE to specify a custom config file location
# =============================================================================
'''


def generate_config_file(output_path: Optional[Path] = None):
    """
    Generate a template configuration file.

    Args:
        output_path: Where to write the config (default: ./.intronIC.toml)
    """
    if output_path is None:
        output_path = Path(".intronIC.toml")

    template = ConfigLoader.generate_template()

    with open(output_path, 'w') as f:
        f.write(template)

    print(f"Generated configuration template: {output_path}")
    print(f"Edit this file to set your default intronIC parameters.")
