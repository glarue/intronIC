"""
Configuration management for intronIC CLI.

Handles pipeline configuration derived from command-line arguments.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


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
            return 'sequences'
        elif self.bed:
            return 'bed'
        elif self.annotation:
            return 'annotation'
        else:
            raise ValueError("No valid input source specified")


@dataclass(frozen=True, slots=True)
class ExtractionConfig:
    """Configuration for intron extraction."""
    feature_type: str = 'both'  # 'cds', 'exon', or 'both'
    min_intron_len: int = 30
    flank_len: int = 50
    allow_multiple_isoforms: bool = False
    no_intron_overlap: bool = False


@dataclass(frozen=True, slots=True)
class ScoringConfig:
    """Configuration for scoring and classification."""
    sequences_only: bool = False
    threshold: float = 90.0
    exclude_noncanonical: bool = False
    scoring_regions: ScoringRegions = field(default_factory=lambda: ScoringRegions(
        five_start=-3, five_end=9,
        bp_start=-55, bp_end=-5,
        three_start=-6, three_end=4
    ))
    pwm_file: Optional[Path] = None
    reference_u12s: Optional[Path] = None
    reference_u2s: Optional[Path] = None
    generate_u2_bps_pwm: bool = False


@dataclass(frozen=True, slots=True)
class TrainingConfig:
    """Configuration for SVM training."""
    fixed_C: Optional[float] = None
    n_models: int = 1
    recursive: bool = False
    recursive_subset: Optional[int] = None
    seed: int = 42


@dataclass(frozen=True, slots=True)
class PerformanceConfig:
    """Configuration for performance settings."""
    processes: int = 1
    cv_processes: int = 1


@dataclass(frozen=True, slots=True)
class OutputConfig:
    """Configuration for output."""
    output_dir: Path
    species_name: str
    quiet: bool = False
    debug: bool = False

    @property
    def base_filename(self) -> str:
        """Generate base filename for outputs.

        Returns:
            Base filename (e.g., 'homo_sapiens')
        """
        return self.species_name

    def get_output_path(self, suffix: str) -> Path:
        """Get output file path with given suffix.

        Args:
            suffix: File suffix (e.g., '.meta.iic')

        Returns:
            Full output path
        """
        return self.output_dir / f"{self.base_filename}{suffix}"


@dataclass(frozen=True, slots=True)
class IntronICConfig:
    """Complete intronIC pipeline configuration."""
    input: InputConfig
    extraction: ExtractionConfig
    scoring: ScoringConfig
    training: TrainingConfig
    performance: PerformanceConfig
    output: OutputConfig

    @classmethod
    def from_args(cls, args) -> 'IntronICConfig':
        """Create configuration from parsed arguments.

        Args:
            args: Parsed argument namespace

        Returns:
            IntronICConfig instance
        """
        # Input configuration
        input_config = InputConfig(
            genome=args.genome,
            annotation=args.annotation,
            bed=args.bed,
            sequence_file=args.sequence_file
        )

        # Extraction configuration
        extraction_config = ExtractionConfig(
            feature_type=args.feature,
            min_intron_len=args.min_intron_len,
            flank_len=args.flank_len,
            allow_multiple_isoforms=args.allow_multiple_isoforms,
            no_intron_overlap=args.no_intron_overlap
        )

        # Scoring regions
        scoring_regions = ScoringRegions(
            five_start=args.five_score_coords[0],
            five_end=args.five_score_coords[1],
            bp_start=args.bp_region_coords[0],
            bp_end=args.bp_region_coords[1],
            three_start=args.three_score_coords[0],
            three_end=args.three_score_coords[1]
        )

        # Scoring configuration
        scoring_config = ScoringConfig(
            sequences_only=args.sequences_only,
            threshold=args.threshold,
            exclude_noncanonical=args.no_nc,
            scoring_regions=scoring_regions,
            pwm_file=args.pwms,
            reference_u12s=args.reference_u12s,
            reference_u2s=args.reference_u2s,
            generate_u2_bps_pwm=args.generate_u2_bps_pwm
        )

        # Training configuration
        recursive_subset = None
        if args.recursive and isinstance(args.recursive, str):
            try:
                recursive_subset = int(args.recursive)
            except ValueError:
                pass

        training_config = TrainingConfig(
            fixed_C=args.C,
            n_models=args.n_models,
            recursive=bool(args.recursive),
            recursive_subset=recursive_subset,
            seed=args.seed
        )

        # Performance configuration
        performance_config = PerformanceConfig(
            processes=args.processes,
            cv_processes=args.cv_processes
        )

        # Output configuration
        output_config = OutputConfig(
            output_dir=args.output_dir,
            species_name=args.species_name,
            quiet=args.quiet,
            debug=args.debug
        )

        return cls(
            input=input_config,
            extraction=extraction_config,
            scoring=scoring_config,
            training=training_config,
            performance=performance_config,
            output=output_config
        )
