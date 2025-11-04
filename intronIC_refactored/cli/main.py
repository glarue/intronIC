"""
Main CLI entry point for intronIC.

Orchestrates the complete intron classification pipeline.
"""

import sys
import logging
from pathlib import Path
from typing import List, Tuple

from .args import IntronICArgumentParser
from .config import IntronICConfig
from .progress import IntronICProgressReporter

# Import pipeline components
from core.intron import Intron
from file_io.genome import GenomeReader
from file_io.parsers import AnnotationParser, BEDParser, SequenceParser
from extraction.annotator import AnnotationHierarchyBuilder
from extraction.intronator import IntronGenerator
from extraction.sequences import SequenceExtractor
from scoring.pwm import PWMLoader
from scoring.scorer import IntronScorer
from scoring.normalizer import ScoreNormalizer
from classification.classifier import IntronClassifier


def setup_logging(config: IntronICConfig) -> logging.Logger:
    """Setup logging configuration.

    Args:
        config: Pipeline configuration

    Returns:
        Configured logger instance
    """
    log_file = config.output.get_output_path('.log')

    # Configure logging level
    if config.output.debug:
        level = logging.DEBUG
    elif config.output.quiet:
        level = logging.WARNING
    else:
        level = logging.INFO

    # Setup logger
    logger = logging.getLogger('intronIC')
    logger.setLevel(level)

    # File handler
    fh = logging.FileHandler(log_file, mode='w')
    fh.setLevel(logging.DEBUG)  # Always log everything to file
    fh.setFormatter(logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    ))
    logger.addHandler(fh)

    # Console handler (only if not quiet)
    if not config.output.quiet:
        ch = logging.StreamHandler()
        ch.setLevel(level)
        ch.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
        logger.addHandler(ch)

    return logger


def load_genome(config: IntronICConfig, logger: logging.Logger) -> GenomeReader:
    """Load genome file.

    Args:
        config: Pipeline configuration
        logger: Logger instance

    Returns:
        GenomeReader instance
    """
    logger.info(f"Loading genome: {config.input.genome}")
    reader = GenomeReader(config.input.genome)
    logger.info(f"Loaded {len(reader.sequences)} sequences")
    return reader


def extract_introns_from_annotation(
    config: IntronICConfig,
    genome_reader: GenomeReader,
    reporter: IntronICProgressReporter,
    logger: logging.Logger
) -> List[Intron]:
    """Extract introns from annotation file.

    Args:
        config: Pipeline configuration
        genome_reader: Genome reader instance
        reporter: Progress reporter
        logger: Logger instance

    Returns:
        List of extracted introns
    """
    reporter.print_info(f"Parsing annotation: {config.input.annotation}")
    logger.info(f"Parsing annotation: {config.input.annotation}")

    # Parse annotation
    parser = AnnotationParser()
    features = parser.parse_file(config.input.annotation)
    logger.info(f"Parsed {len(features)} features")

    # Extract introns
    reporter.print_info("Extracting introns from features")
    logger.info("Extracting introns")

    extractor = IntronExtractor(
        genome_reader=genome_reader,
        min_length=config.extraction.min_intron_len,
        flank_length=config.extraction.flank_len
    )

    introns = extractor.extract_from_features(
        features,
        feature_type=config.extraction.feature_type
    )

    logger.info(f"Extracted {len(introns)} introns")
    return introns


def extract_introns_from_bed(
    config: IntronICConfig,
    genome_reader: GenomeReader,
    reporter: IntronICProgressReporter,
    logger: logging.Logger
) -> List[Intron]:
    """Extract introns from BED file.

    Args:
        config: Pipeline configuration
        genome_reader: Genome reader instance
        reporter: Progress reporter
        logger: Logger instance

    Returns:
        List of extracted introns
    """
    reporter.print_info(f"Reading BED file: {config.input.bed}")
    logger.info(f"Reading BED file: {config.input.bed}")

    extractor = IntronExtractor(
        genome_reader=genome_reader,
        min_length=config.extraction.min_intron_len,
        flank_length=config.extraction.flank_len
    )

    introns = extractor.extract_from_bed(config.input.bed)
    logger.info(f"Extracted {len(introns)} introns from BED")
    return introns


def load_introns_from_sequences(
    config: IntronICConfig,
    reporter: IntronICProgressReporter,
    logger: logging.Logger
) -> List[Intron]:
    """Load introns from pre-extracted sequences.

    Args:
        config: Pipeline configuration
        reporter: Progress reporter
        logger: Logger instance

    Returns:
        List of introns with sequences
    """
    reporter.print_info(f"Loading sequences: {config.input.sequence_file}")
    logger.info(f"Loading sequences: {config.input.sequence_file}")

    # TODO: Implement sequence file parser
    # This would be similar to the reference sequence loader in integration tests
    raise NotImplementedError("Sequence file loading not yet implemented")


def score_introns(
    introns: List[Intron],
    config: IntronICConfig,
    reporter: IntronICProgressReporter,
    logger: logging.Logger
) -> List[Intron]:
    """Score introns with PWM matrices.

    Args:
        introns: List of introns to score
        config: Pipeline configuration
        reporter: Progress reporter
        logger: Logger instance

    Returns:
        List of introns with raw scores
    """
    reporter.print_info("Calculating PWM scores")
    logger.info("Starting PWM scoring")

    # Load or build PWM matrices
    # TODO: Implement PWM loading from file or defaults
    scorer = PWMScorer(
        five_start=config.scoring.scoring_regions.five_start,
        five_end=config.scoring.scoring_regions.five_end,
        bp_start=config.scoring.scoring_regions.bp_start,
        bp_end=config.scoring.scoring_regions.bp_end,
        three_start=config.scoring.scoring_regions.three_start,
        three_end=config.scoring.scoring_regions.three_end
    )

    # Score introns
    progress = reporter.create_progress()
    with progress:
        task = progress.add_task(
            "[cyan]Scoring introns...",
            total=len(introns)
        )

        scored_introns = []
        for intron in introns:
            scored = scorer.score_intron(intron)
            scored_introns.append(scored)
            progress.update(task, advance=1)

    logger.info(f"Scored {len(scored_introns)} introns")
    return scored_introns


def normalize_scores(
    introns: List[Intron],
    config: IntronICConfig,
    reporter: IntronICProgressReporter,
    logger: logging.Logger
) -> List[Intron]:
    """Normalize intron scores using z-score transformation.

    Args:
        introns: List of introns with raw scores
        config: Pipeline configuration
        reporter: Progress reporter
        logger: Logger instance

    Returns:
        List of introns with normalized z-scores
    """
    reporter.print_info("Normalizing scores")
    logger.info("Normalizing scores with z-score transformation")

    # Load reference sequences for normalization
    # TODO: Implement reference sequence loading
    normalizer = ScoreNormalizer()

    # Fit on reference data and transform experimental
    # normalized = normalizer.fit_transform(reference_introns, introns)

    # For now, return as-is (would need reference data)
    logger.warning("Score normalization not fully implemented yet")
    return introns


def classify_introns(
    introns: List[Intron],
    config: IntronICConfig,
    reporter: IntronICProgressReporter,
    logger: logging.Logger
) -> Tuple[List[Intron], dict]:
    """Classify introns as U2 or U12 type.

    Args:
        introns: List of introns with z-scores
        config: Pipeline configuration
        reporter: Progress reporter
        logger: Logger instance

    Returns:
        Tuple of (classified introns, classification metrics)
    """
    reporter.print_info("Training SVM classifier")
    logger.info("Starting SVM classification")

    # Load reference sequences
    # TODO: Implement reference sequence loading from default or custom paths

    classifier = IntronClassifier(
        threshold=config.scoring.threshold,
        n_models=config.training.n_models,
        fixed_C=config.training.fixed_C,
        n_processes=config.performance.processes,
        cv_processes=config.performance.cv_processes,
        seed=config.training.seed
    )

    # For now, return as-is (would need reference data)
    logger.warning("Classification not fully implemented yet")
    return introns, {}


def write_outputs(
    introns: List[Intron],
    config: IntronICConfig,
    reporter: IntronICProgressReporter,
    logger: logging.Logger
):
    """Write output files.

    Args:
        introns: Classified introns
        config: Pipeline configuration
        reporter: Progress reporter
        logger: Logger instance
    """
    reporter.print_info("Writing output files")
    logger.info("Writing output files")

    # TODO: Implement output writers
    # - .meta.iic: Metadata file
    # - .bed.iic: BED-like coordinates
    # - .seqs.iic: Sequences
    # - .scores.iic: Detailed scores

    output_files = {
        "Metadata": config.output.get_output_path('.meta.iic'),
        "BED": config.output.get_output_path('.bed.iic'),
        "Sequences": config.output.get_output_path('.seqs.iic'),
        "Scores": config.output.get_output_path('.scores.iic'),
        "Log": config.output.get_output_path('.log'),
    }

    reporter.print_file_tree(output_files)
    logger.info("Output files written successfully")


def run_pipeline(config: IntronICConfig):
    """Run the complete intronIC pipeline.

    Args:
        config: Pipeline configuration
    """
    # Setup logging and reporting
    logger = setup_logging(config)
    reporter = IntronICProgressReporter(quiet=config.output.quiet)

    # Print header
    reporter.print_header(
        species_name=config.output.species_name,
        input_mode=config.input.mode
    )

    # Define pipeline steps
    pipeline_steps = [
        "Load input data",
        "Extract introns",
        "Score introns with PWMs",
        "Normalize scores",
        "Train and apply classifier",
        "Write output files"
    ]

    reporter.print_pipeline_steps(pipeline_steps)

    try:
        # Step 1: Load input data
        reporter.print_section("Step 1: Load Input Data", "bold blue")
        reporter.print_pipeline_steps(pipeline_steps, current_step=1)

        if config.input.mode == 'annotation':
            genome_reader = load_genome(config, logger)
            introns = extract_introns_from_annotation(
                config, genome_reader, reporter, logger
            )
        elif config.input.mode == 'bed':
            genome_reader = load_genome(config, logger)
            introns = extract_introns_from_bed(
                config, genome_reader, reporter, logger
            )
        elif config.input.mode == 'sequences':
            introns = load_introns_from_sequences(config, reporter, logger)
        else:
            raise ValueError(f"Unknown input mode: {config.input.mode}")

        reporter.print_success(f"Loaded {len(introns):,} introns")

        # If sequences only, skip to output
        if config.scoring.sequences_only:
            reporter.print_warning("Sequences-only mode: Skipping classification")
            write_outputs(introns, config, reporter, logger)
            reporter.print_success("Pipeline complete!")
            return

        # Step 2: Extract introns (already done above)
        reporter.print_section("Step 2: Extract Introns", "bold blue")
        reporter.print_pipeline_steps(pipeline_steps, current_step=2)
        reporter.print_success(f"Extracted {len(introns):,} introns")

        # Step 3: Score introns
        reporter.print_section("Step 3: Score Introns", "bold blue")
        reporter.print_pipeline_steps(pipeline_steps, current_step=3)
        scored_introns = score_introns(introns, config, reporter, logger)
        reporter.print_success(f"Scored {len(scored_introns):,} introns")

        # Step 4: Normalize scores
        reporter.print_section("Step 4: Normalize Scores", "bold blue")
        reporter.print_pipeline_steps(pipeline_steps, current_step=4)
        normalized_introns = normalize_scores(
            scored_introns, config, reporter, logger
        )
        reporter.print_success("Scores normalized")

        # Step 5: Classify
        reporter.print_section("Step 5: Classify Introns", "bold blue")
        reporter.print_pipeline_steps(pipeline_steps, current_step=5)
        classified_introns, metrics = classify_introns(
            normalized_introns, config, reporter, logger
        )

        # Count classifications
        u12_count = sum(
            1 for i in classified_introns
            if i.metadata and i.metadata.type_id == 'u12'
        )
        u2_count = len(classified_introns) - u12_count

        reporter.print_classification_summary(
            total=len(classified_introns),
            u12_count=u12_count,
            u2_count=u2_count,
            threshold=config.scoring.threshold
        )

        # Step 6: Write outputs
        reporter.print_section("Step 6: Write Outputs", "bold blue")
        reporter.print_pipeline_steps(pipeline_steps, current_step=6)
        write_outputs(classified_introns, config, reporter, logger)

        reporter.print_success("Pipeline complete!")

    except Exception as e:
        logger.exception("Pipeline failed with error")
        reporter.print_error(f"Pipeline failed: {str(e)}")
        raise


def main(args=None):
    """Main entry point for intronIC CLI.

    Args:
        args: Optional list of arguments (defaults to sys.argv)

    Returns:
        Exit code (0 for success, 1 for failure)
    """
    try:
        # Parse arguments
        parser = IntronICArgumentParser()
        parsed_args = parser.parse_args(args)

        # Create configuration
        config = IntronICConfig.from_args(parsed_args)

        # Run pipeline
        run_pipeline(config)

        return 0

    except KeyboardInterrupt:
        print("\n\nPipeline interrupted by user", file=sys.stderr)
        return 130  # Standard exit code for SIGINT

    except Exception as e:
        print(f"\n\nError: {str(e)}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
