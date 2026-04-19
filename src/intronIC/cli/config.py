"""
Configuration management for intronIC CLI.

Handles pipeline configuration derived from command-line arguments.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


def get_default_pretrained_model_path() -> Optional[Path]:
    """Get path to default pretrained model in data directory.

    Returns:
        Path to default pretrained model if it exists, None otherwise
    """
    # Path relative to this file: cli/config.py -> data/default_pretrained.model.pkl
    default_path = (
        Path(__file__).parent.parent / "data" / "default_pretrained.model.pkl"
    )
    if default_path.exists():
        return default_path
    return None


@dataclass(frozen=True, slots=True)
class ScoringRegions:
    """Coordinates for scoring regions."""

    five_start: int
    five_end: int
    bp_start: int
    bp_end: int
    three_start: int
    three_end: int


@dataclass(frozen=True, slots=True)
class InputConfig:
    """Input file configuration."""

    genome: Optional[Path] = None
    annotation: Optional[Path] = None
    bed: Optional[Path] = None
    sequence_file: Optional[Path] = None

    @property
    def mode(self) -> str:
        """Determine input mode.

        Returns:
            One of: 'annotation', 'bed', 'sequences'
        """
        if self.sequence_file:
            return "sequences"
        elif self.bed:
            return "bed"
        elif self.annotation:
            return "annotation"
        else:
            raise ValueError("No valid input source specified")


@dataclass(frozen=True, slots=True)
class ExtractionConfig:
    """Configuration for intron extraction."""

    feature_type: str = "both"  # 'cds', 'exon', or 'both'
    min_intron_len: int = 30
    flank_len: int = (
        100  # Length of flanking sequence to extract (upstream and downstream)
    )
    allow_multiple_isoforms: bool = False
    no_intron_overlap: bool = False
    include_duplicates: bool = False
    u12_boundary_correction: bool = True  # Enable U12 correction by default
    u12_correction_require_canonical: bool = True  # Only correct if result is canonical (GT-AG, GC-AG, AT-AC)


@dataclass(frozen=True, slots=True)
class ScoringConfig:
    """Configuration for scoring and classification."""

    threshold: float = 90.0
    exclude_noncanonical: bool = False
    scoring_regions: ScoringRegions = field(
        default_factory=lambda: ScoringRegions(
            five_start=-3,
            five_end=10,
            bp_start=-55,
            bp_end=-5,
            three_start=-14,
            three_end=4,
        )
    )
    pwm_file: Optional[Path] = None
    reference_u12s: Optional[Path] = None
    reference_u2s: Optional[Path] = None
    generate_u2_bps_pwm: bool = False
    pseudocount: float = 0.0001
    ignore_nc_dnts: bool = (
        True  # Ignore terminal dinucleotides for NC introns by default
    )
    normalizer_mode: str = "auto"  # 'human', 'adaptive', or 'auto'
    species_prior: Optional[float] = None  # Expected U12 prior for target species
    load_normalizer: Optional[Path] = (
        None  # Load saved normalizer for reproducible normalization
    )
    save_normalizer: bool = False  # Save fitted normalizer for future runs


@dataclass(frozen=True, slots=True)
class TrainingConfig:
    """Configuration for SVM training."""

    fixed_C: Optional[float] = None
    n_models: int = 1
    recursive: bool = False
    recursive_subset: Optional[int] = None
    seed: int = 42
    max_iter: int = 50000
    eval_mode: str = "nested_cv"
    n_cv_folds: int = 7
    test_fraction: float = 0.2
    n_optimization_rounds: int = 5
    pretrained_model_path: Optional[Path] = None
    config_path: Optional[Path] = (
        None  # Path to unified config file (auto-loads if None)
    )
    use_fold_averaged_params: bool = (
        False  # Use fold-averaged hyperparameters (better cross-species)
    )


@dataclass(frozen=True, slots=True)
class OptimizerConfig:
    """Configuration for hyperparameter optimization.

    Loaded from YAML config 'optimizer' section. This consolidates all
    optimizer settings in one place, loaded once at startup.
    """

    # Grid search parameters
    n_rounds: int = 5
    n_points_initial: int = 13
    n_points_refine: int = 50
    cv_folds: int = 5
    max_iter: int = 60000
    random_state: int = 42
    n_jobs: int = -1
    verbose: bool = True

    # Scoring and regularization
    scoring_metric: str = "balanced_accuracy"
    penalty_options: tuple = ("l2",)
    loss_options: tuple = ("squared_hinge",)
    class_weight_multipliers: tuple = (1.0,)
    use_multiplier_tiebreaker: bool = True

    # Feature transformation
    features: Optional[tuple] = None  # None = default 4D (absdiff_bp_3)
    extra_features: Optional[tuple] = None  # Additional IntronScores fields for classifier

    # Gamma imbalance scaling
    gamma_imbalance_options: Optional[tuple] = None

    # SVM kernel configuration
    kernel: str = 'rbf'  # 'rbf' or 'linear'
    gamma: str = 'scale'  # RBF gamma: 'scale', 'auto', or float (ignored for linear)
    gamma_search: Optional[tuple] = None  # Grid search values for gamma (RBF only)

    # C parameter bounds
    eff_C_pos_range: tuple = (1e-3, 1e3)
    eff_C_neg_max: Optional[float] = None

    # Early stopping
    early_stop_on_plateau: bool = True
    plateau_rounds: int = 3

    # Parameter grid override (for custom grids)
    param_grid_override: Optional[dict] = None

    @classmethod
    def from_yaml(cls, yaml_config: dict) -> "OptimizerConfig":
        """Create OptimizerConfig from YAML config dict.

        Args:
            yaml_config: Full YAML config dictionary

        Returns:
            OptimizerConfig instance with values from YAML
        """
        opt = yaml_config.get("optimizer", {})

        # Convert lists to tuples for immutability
        def to_tuple(val, default):
            if val is None:
                return default
            if isinstance(val, (list, tuple)):
                return tuple(val)
            return (val,)

        # Extract C bounds
        c_bounds = opt.get("c_bounds", {})

        # Extract feature transform
        ft = opt.get("feature_transform", {})
        features = ft.get("features")
        if features is not None:
            features = tuple(features) if isinstance(features, list) else (features,)

        extra_features = ft.get("extra_features")
        if extra_features is not None:
            extra_features = tuple(extra_features) if isinstance(extra_features, list) else (extra_features,)

        return cls(
            n_rounds=opt.get("n_rounds", 5),
            n_points_initial=opt.get("n_points_initial", 13),
            n_points_refine=opt.get("n_points_refine", 50),
            cv_folds=opt.get("cv_folds", 5),
            max_iter=opt.get("max_iter", 60000),
            random_state=opt.get("random_state", 42),
            n_jobs=opt.get("n_jobs", -1),
            verbose=opt.get("verbose", True),
            scoring_metric=opt.get("scoring_metric", "balanced_accuracy"),
            penalty_options=to_tuple(opt.get("penalty_options"), ("l2",)),
            loss_options=to_tuple(opt.get("loss_options"), ("squared_hinge",)),
            class_weight_multipliers=to_tuple(
                opt.get("class_weight_multipliers"), (1.0,)
            ),
            use_multiplier_tiebreaker=opt.get("use_multiplier_tiebreaker", True),
            features=features,
            extra_features=extra_features,
            gamma_imbalance_options=to_tuple(opt.get("gamma_imbalance_options"), None)
            if opt.get("gamma_imbalance_options")
            else None,
            kernel=opt.get("kernel", "rbf"),
            gamma=str(opt.get("gamma", "scale")),
            gamma_search=to_tuple(opt.get("gamma_search"), None)
            if opt.get("gamma_search")
            else None,
            eff_C_pos_range=tuple(c_bounds.get("eff_C_pos_range", [1e-3, 1e3])),
            eff_C_neg_max=c_bounds.get("eff_C_neg_max"),
            early_stop_on_plateau=opt.get("early_stop_on_plateau", True),
            plateau_rounds=opt.get("plateau_rounds", 3),
            param_grid_override=yaml_config.get("param_grid"),
        )


@dataclass(frozen=True, slots=True)
class EnsembleConfig:
    """Configuration for ensemble training.

    Loaded from YAML config 'training.ensemble' section.
    """

    n_models: int = 1
    subsample_u2: bool = True
    subsample_ratio: float = 0.85
    max_iter: int = 50000
    random_state: int = 42

    @classmethod
    def from_yaml(cls, yaml_config: dict) -> "EnsembleConfig":
        """Create EnsembleConfig from YAML config dict.

        Args:
            yaml_config: Full YAML config dictionary

        Returns:
            EnsembleConfig instance with values from YAML
        """
        training = yaml_config.get("training", {})
        ensemble = training.get("ensemble", {})

        return cls(
            n_models=ensemble.get("n_models", 1),
            subsample_u2=ensemble.get("subsample_u2", True),
            subsample_ratio=ensemble.get("subsample_ratio", 0.85),
            max_iter=ensemble.get("max_iter", 50000),
            random_state=ensemble.get("random_state", 42),
        )


@dataclass(frozen=True, slots=True)
class SpeciesBackgroundConfig:
    """Configuration for species-specific U2 background correction.

    Controls whether and how species-level nucleotide frequencies are
    blended into U2 PWMs to correct for composition bias.  Defaults are
    conservative (enabled, no trimming, no iteration).

    Loaded from YAML config 'scoring.species_background' section.
    """

    enabled: bool = True
    n0: int = 1000  # Bayesian shrinkage: w = n / (n + n0)
    trim_percentile: float = 0  # Fraction of top-scoring introns to exclude (0 = none)
    pseudocount_per_base: float = 1.0  # Dirichlet pseudocount per base per position
    n_iterations: int = 0  # Iterative re-trimming passes (0 = single pass, no trim)
    min_introns: int = 500  # Minimum introns required; below this, skip correction

    @classmethod
    def from_yaml(cls, yaml_config: dict) -> "SpeciesBackgroundConfig":
        """Create SpeciesBackgroundConfig from YAML config dict.

        Args:
            yaml_config: Full YAML config dictionary

        Returns:
            SpeciesBackgroundConfig instance with values from YAML
        """
        scoring = yaml_config.get("scoring", {})
        bg = scoring.get("species_background", {})

        return cls(
            enabled=bg.get("enabled", True),
            n0=bg.get("n0", 1000),
            trim_percentile=bg.get("trim_percentile", 0),
            pseudocount_per_base=bg.get("pseudocount_per_base", 1.0),
            n_iterations=bg.get("n_iterations", 0),
            min_introns=bg.get("min_introns", 500),
        )


@dataclass(frozen=True, slots=True)
class PerformanceConfig:
    """Configuration for performance settings."""

    processes: int = 1
    cv_processes: int = 1
    streaming: bool = True  # Streaming mode (default): ~85% memory savings, faster with parallelization


@dataclass(frozen=True, slots=True)
class OutputConfig:
    """Configuration for output."""

    output_dir: Path
    species_name: str
    clean_names: bool = True
    quiet: bool = False
    debug: bool = False
    uninformative_naming: bool = False  # Use simple naming (species-i_ID)
    no_abbreviate: bool = False  # Use full species name in output (not abbreviated)
    abbreviate_filenames: bool = False  # Use abbreviated species name in filenames
    no_headers: bool = False  # Omit column headers from output files

    @property
    def base_filename(self) -> str:
        """Generate base filename for outputs.

        Returns:
            Base filename (e.g., 'homo_sapiens' or 'HomSap' if abbreviated)
        """
        if self.abbreviate_filenames:
            # Abbreviate for filenames (3+3 format)
            from intronIC.file_io.writers import generate_species_abbreviation

            return generate_species_abbreviation(self.species_name)
        return self.species_name

    def get_output_path(self, suffix: str) -> Path:
        """Get output file path with given suffix.

        Args:
            suffix: File suffix (e.g., '.meta.iic')

        Returns:
            Full output path
        """
        return self.output_dir / f"{self.base_filename}{suffix}"


@dataclass(frozen=True)
class IntronICConfig:
    """Complete intronIC pipeline configuration.

    This is the single source of truth for all configuration. It is created once
    at startup from YAML config + CLI args, with CLI taking precedence.

    Attributes:
        input: Input file configuration
        extraction: Intron extraction settings
        scoring: PWM scoring settings
        training: Basic training settings (n_models, eval_mode, etc.)
        performance: Parallelization and memory settings
        output: Output file settings
        optimizer: Hyperparameter optimization settings (from YAML 'optimizer' section)
        ensemble: Ensemble training settings (from YAML 'training.ensemble' section)
        species_background: Species-specific U2 background correction settings
            (from YAML 'scoring.species_background' section)
        yaml_config: Raw YAML config dict (for any custom/advanced settings)
        config_path: Path to the loaded config file (for logging)
    """

    input: InputConfig
    extraction: ExtractionConfig
    scoring: ScoringConfig
    training: TrainingConfig
    performance: PerformanceConfig
    output: OutputConfig
    optimizer: OptimizerConfig = field(default_factory=OptimizerConfig)
    ensemble: EnsembleConfig = field(default_factory=EnsembleConfig)
    species_background: SpeciesBackgroundConfig = field(
        default_factory=SpeciesBackgroundConfig
    )
    yaml_config: dict = field(default_factory=dict)
    config_path: Optional[Path] = None

    @classmethod
    def from_yaml_and_args(cls, args, yaml_config: dict = None) -> "IntronICConfig":
        """Create configuration from YAML config and CLI arguments.

        CLI arguments take precedence over YAML config values.
        This is the primary entry point - creates a unified config from both sources.

        The config is loaded ONCE here and stored. All optimizer/ensemble settings
        are extracted from the YAML and stored in typed dataclasses. The raw YAML
        is also kept for any advanced/custom settings.

        Args:
            args: Parsed argument namespace from CLI
            yaml_config: Optional YAML config dict (if None, will auto-load)

        Returns:
            IntronICConfig instance with all settings populated
        """
        from intronIC.utils.config_loader import find_config, load_config

        # Load YAML config if not provided
        config_path_arg = getattr(args, "config_path", None)
        if yaml_config is None:
            yaml_config = load_config(config_path_arg) or {}

        # Find the actual config path for logging
        found_config_path = find_config(
            str(config_path_arg) if config_path_arg else None
        )

        # Merge YAML into args (CLI args take precedence)
        # Only set arg from YAML if it wasn't explicitly set on CLI
        merged_args = cls._merge_yaml_into_args(args, yaml_config)

        # Create optimizer, ensemble, and species background configs from YAML
        optimizer_config = OptimizerConfig.from_yaml(yaml_config)
        ensemble_config = EnsembleConfig.from_yaml(yaml_config)
        species_background_config = SpeciesBackgroundConfig.from_yaml(yaml_config)

        # Now create config from merged args using existing logic
        # Pass the extra configs to from_args
        return cls.from_args(
            merged_args,
            optimizer_config=optimizer_config,
            ensemble_config=ensemble_config,
            species_background_config=species_background_config,
            yaml_config=yaml_config,
            config_path=found_config_path,
        )

    @staticmethod
    def _merge_yaml_into_args(args, yaml_config: dict):
        """Merge YAML config into argparse namespace.

        Only updates args that have their default values (weren't explicitly set on CLI).

        Args:
            args: Argparse namespace with _arg_defaults attribute from parser
            yaml_config: YAML config dictionary

        Returns:
            Updated args namespace
        """
        # Map YAML paths to arg attributes and their types
        # Format: (yaml_path, arg_attr, converter)
        # Default values are extracted from args._arg_defaults (single source of truth)
        mappings = [
            # Scoring
            ("scoring.threshold", "threshold", float),
            ("scoring.feature_type", "feature", str),
            ("scoring.exclude_noncanonical", "no_nc", bool),
            ("scoring.pseudocount", "pseudocount", float),
            ("scoring.normalizer_mode", "normalizer_mode", str),
            # Extraction
            ("extraction.flank_length", "flank_len", int),
            ("extraction.min_intron_length", "min_intron_len", int),
            # Performance
            ("performance.processes", "processes", int),
            ("performance.cv_processes", "cv_processes", int),
            # Training
            ("training.n_models", "n_models", int),
            ("training.max_iterations", "max_iter", int),
            ("training.eval_mode", "eval_mode", str),
            ("training.n_cv_folds", "n_cv_folds", int),
            ("training.fixed_C", "C", float),
            # Advanced
            ("advanced.random_seed", "seed", int),
        ]

        # Get argparse defaults from the namespace (single source of truth)
        arg_defaults = getattr(args, "_arg_defaults", {})

        for yaml_path, arg_attr, converter in mappings:
            # Get value from YAML (handle nested paths)
            yaml_value = yaml_config
            for key in yaml_path.split("."):
                if isinstance(yaml_value, dict):
                    yaml_value = yaml_value.get(key)
                else:
                    yaml_value = None
                    break

            if yaml_value is None:
                continue

            # Get the argparse default for this argument (single source of truth)
            default = arg_defaults.get(arg_attr)

            # Only apply if arg still has default value (CLI didn't override)
            # For None defaults, check if attribute is None
            # For other defaults, check if value equals default
            current_value = getattr(args, arg_attr, default)
            if default is None:
                should_apply = current_value is None
            else:
                should_apply = current_value == default

            if should_apply:
                try:
                    setattr(args, arg_attr, converter(yaml_value))
                except (ValueError, TypeError):
                    pass  # Silently skip invalid values

        # Handle scoring region overrides (list-valued, can't use scalar mapping)
        regions = yaml_config.get("scoring", {}).get("regions", {})
        if regions:
            region_mappings = [
                ("five_prime", "five_score_coords"),
                ("branch_point", "bp_region_coords"),
                ("three_prime", "three_score_coords"),
            ]
            for yaml_key, arg_attr in region_mappings:
                region = regions.get(yaml_key)
                if region and isinstance(region, dict):
                    start = region.get("start")
                    end = region.get("end")
                    if start is not None and end is not None:
                        default = arg_defaults.get(arg_attr)
                        current = getattr(args, arg_attr, default)
                        if current == default:
                            setattr(args, arg_attr, [int(start), int(end)])

        return args

    @classmethod
    def from_args(
        cls,
        args,
        optimizer_config: OptimizerConfig = None,
        ensemble_config: EnsembleConfig = None,
        species_background_config: SpeciesBackgroundConfig = None,
        yaml_config: dict = None,
        config_path: Path = None,
    ) -> "IntronICConfig":
        """Create configuration from parsed arguments.

        Args:
            args: Parsed argument namespace
            optimizer_config: Pre-built OptimizerConfig (from YAML)
            ensemble_config: Pre-built EnsembleConfig (from YAML)
            species_background_config: Pre-built SpeciesBackgroundConfig (from YAML)
            yaml_config: Raw YAML config dict
            config_path: Path to loaded config file

        Returns:
            IntronICConfig instance
        """
        # Use defaults if not provided
        if optimizer_config is None:
            optimizer_config = OptimizerConfig()
        if ensemble_config is None:
            ensemble_config = EnsembleConfig()
        if species_background_config is None:
            species_background_config = SpeciesBackgroundConfig()
        if yaml_config is None:
            yaml_config = {}

        # Input configuration
        input_config = InputConfig(
            genome=args.genome,
            annotation=args.annotation,
            bed=args.bed,
            sequence_file=args.sequence_file,
        )

        # Extraction configuration
        # Get require_canonical from YAML (not a CLI arg)
        _require_canonical = True  # default
        if yaml_config:
            _scoring = yaml_config.get("scoring", {})
            _require_canonical = _scoring.get(
                "u12_correction_require_canonical", True
            )

        extraction_config = ExtractionConfig(
            feature_type=args.feature,
            min_intron_len=args.min_intron_len,
            flank_len=args.flank_len,
            allow_multiple_isoforms=args.allow_multiple_isoforms,
            no_intron_overlap=args.no_intron_overlap,
            include_duplicates=args.include_duplicates,
            u12_boundary_correction=not args.no_nc_ss_adjustment,
            u12_correction_require_canonical=_require_canonical,
        )

        # Scoring regions
        scoring_regions = ScoringRegions(
            five_start=args.five_score_coords[0],
            five_end=args.five_score_coords[1],
            bp_start=args.bp_region_coords[0],
            bp_end=args.bp_region_coords[1],
            three_start=args.three_score_coords[0],
            three_end=args.three_score_coords[1],
        )

        # Scoring configuration
        scoring_config = ScoringConfig(
            threshold=args.threshold,
            exclude_noncanonical=args.no_nc,
            scoring_regions=scoring_regions,
            pwm_file=getattr(args, "pwms", None),
            reference_u12s=getattr(args, "reference_u12s", None),
            reference_u2s=getattr(args, "reference_u2s", None),
            generate_u2_bps_pwm=getattr(
                args, "generate_u2_bps_pwm", False
            ),  # Not currently implemented
            pseudocount=args.pseudocount,
            ignore_nc_dnts=not args.no_ignore_nc_dnts,
            normalizer_mode=args.normalizer_mode,
            species_prior=getattr(args, "species_prior", None),
            load_normalizer=getattr(args, "load_normalizer", None),
            save_normalizer=getattr(args, "save_normalizer", False),
        )

        # Training configuration
        recursive_subset = None
        recursive_arg = getattr(args, "recursive", None)  # Not currently implemented
        if recursive_arg and isinstance(recursive_arg, str):
            try:
                recursive_subset = int(recursive_arg)
            except ValueError:
                pass

        # Determine pretrained model path
        # Priority:
        # 1. If train subcommand OR --train flag: pretrained_model_path = None (force training)
        # 2. If --model <path>: use that specific model
        # 3. If --pretrained_model <path>: use that specific model (backward compat)
        # 4. Default: use default pretrained model

        # Check for training mode:
        # - args.command == 'train': train subcommand (no genome needed)
        # - args.train: classify --train flag (train on-the-fly during classification)
        is_training = getattr(args, "command", None) == "train" or getattr(
            args, "train", False
        )

        if is_training:
            # User explicitly wants to train a new model
            pretrained_model_path = None
        elif args.model:
            # User specified a custom pretrained model (new flag)
            pretrained_model_path = args.model
        elif args.pretrained_model:
            # User specified a custom pretrained model (deprecated flag)
            pretrained_model_path = args.pretrained_model
        else:
            # Default: use pretrained model
            pretrained_model_path = get_default_pretrained_model_path()
            if pretrained_model_path is None:
                raise FileNotFoundError(
                    "Default pretrained model not found at data/default_pretrained.model.pkl. "
                    "Use --train to train a new model instead."
                )

        training_config = TrainingConfig(
            fixed_C=getattr(args, "C", None),
            n_models=getattr(args, "n_models", 1),
            recursive=bool(recursive_arg),
            recursive_subset=recursive_subset,
            seed=args.seed,
            max_iter=getattr(args, "max_iter", 50000),
            eval_mode=getattr(args, "eval_mode", "nested_cv"),
            n_cv_folds=getattr(args, "n_cv_folds", 5),
            test_fraction=getattr(args, "test_fraction", 0.2),
            n_optimization_rounds=getattr(args, "n_optimization_rounds", 5),
            pretrained_model_path=pretrained_model_path,
            config_path=getattr(args, "config_path", None),
            use_fold_averaged_params=getattr(args, "use_fold_averaged_params", False),
        )

        # Performance configuration
        streaming = getattr(args, "streaming", False)

        # Handle processes: CLI args override YAML, default to 1 if neither specified
        processes = args.processes if args.processes is not None else yaml_config.get("performance", {}).get("processes", 1)
        cv_processes = args.cv_processes if args.cv_processes is not None else yaml_config.get("performance", {}).get("cv_processes")

        performance_config = PerformanceConfig(
            processes=processes,
            cv_processes=cv_processes if cv_processes is not None else processes,
            streaming=streaming,
        )

        # Output configuration
        output_config = OutputConfig(
            output_dir=args.output_dir,
            species_name=args.species_name,
            clean_names=args.clean_names,
            quiet=args.quiet,
            debug=args.debug,
            uninformative_naming=args.uninformative_naming,
            no_abbreviate=args.no_abbreviate,
            abbreviate_filenames=args.abbreviate_filenames,
            no_headers=getattr(args, "no_headers", False),
        )

        return cls(
            input=input_config,
            extraction=extraction_config,
            scoring=scoring_config,
            training=training_config,
            performance=performance_config,
            output=output_config,
            optimizer=optimizer_config,
            ensemble=ensemble_config,
            species_background=species_background_config,
            yaml_config=yaml_config,
            config_path=config_path,
        )
