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
from core.intron import Intron, IntronSequences
from utils.coordinates import GenomicCoordinate
from file_io.genome import GenomeReader
from file_io.parsers import AnnotationParser, BEDParser, SequenceParser
from file_io.writers import BEDWriter, MetaWriter, SequenceWriter, ScoreWriter, MappingWriter
from extraction.annotator import AnnotationHierarchyBuilder
from extraction.intronator import IntronGenerator
from extraction.sequences import SequenceExtractor
from scoring.pwm import PWMLoader
from scoring.scorer import IntronScorer
from scoring.normalizer import ScoreNormalizer
from classification.classifier import IntronClassifier
import gzip


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


def load_reference_sequences(filepath: Path, max_count: int = None, min_length: int = 55, logger: logging.Logger = None) -> List[Intron]:
    """
    Load reference intron sequences from .iic.gz file.

    Format (tab-delimited):
    - intron_id
    - score
    - upstream_flank
    - intron_seq
    - downstream_flank
    - sources
    - source_count

    Args:
        filepath: Path to .iic.gz file
        max_count: Maximum number to load (None = all)
        min_length: Minimum intron length required (default: 55bp for scoring regions)
        logger: Optional logger instance

    Returns:
        List of Intron objects with sequences
    """
    introns = []
    skipped_short = 0

    open_fn = gzip.open if str(filepath).endswith('.gz') else open
    with open_fn(filepath, 'rt') as f:
        for line_num, line in enumerate(f, 1):
            # Skip comments
            if line.startswith('#'):
                continue

            # Parse line
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue

            intron_id = fields[0]
            upstream_flank = fields[2]
            intron_seq = fields[3]
            downstream_flank = fields[4]

            # Create dummy coordinates (not needed for scoring)
            coord = GenomicCoordinate(
                chromosome="ref",
                start=1,
                stop=len(intron_seq),
                strand="+",
                system="1-based"
            )

            # Extract terminal dinucleotides
            five_dnt = intron_seq[:2] if len(intron_seq) >= 2 else None
            three_dnt = intron_seq[-2:] if len(intron_seq) >= 2 else None

            # Create IntronSequences
            sequences = IntronSequences(
                seq=intron_seq,
                upstream_flank=upstream_flank,
                downstream_flank=downstream_flank,
                five_prime_dnt=five_dnt,
                three_prime_dnt=three_dnt
            )

            # Skip if too short for scoring regions
            if len(intron_seq) < min_length:
                skipped_short += 1
                continue

            # Create Intron
            intron = Intron(
                intron_id=intron_id,
                coordinates=coord,
                sequences=sequences
            )
            introns.append(intron)

            # Check max count
            if max_count and len(introns) >= max_count:
                break

    if logger:
        logger.info(f"Loaded {len(introns)} reference sequences from {filepath.name}")
        if skipped_short > 0:
            logger.info(f"Skipped {skipped_short} reference introns shorter than {min_length}bp")

    return introns


def load_genome(config: IntronICConfig, logger: logging.Logger) -> GenomeReader:
    """Load genome file.

    Args:
        config: Pipeline configuration
        logger: Logger instance

    Returns:
        GenomeReader instance
    """
    logger.info(f"Loading genome: {config.input.genome}")
    # Use cached mode for faster repeated access
    reader = GenomeReader(config.input.genome, cached=True)
    if reader.cache:
        logger.info(f"Loaded {len(reader.cache)} sequences into memory")
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

    # Step 1: Build annotation hierarchy
    feature_types = ['exon'] if config.extraction.feature_type == 'exon' else ['cds']
    builder = AnnotationHierarchyBuilder(child_features=feature_types)
    genes = builder.build_from_file(config.input.annotation)
    logger.info(f"Built hierarchy with {len(genes)} genes")

    # Step 2: Generate introns from exon pairs
    reporter.print_info("Generating introns from exon pairs")
    generator = IntronGenerator()
    introns_iter = generator.generate_from_genes(genes, builder.feature_index)
    introns_list = list(introns_iter)  # Materialize iterator
    logger.info(f"Generated {len(introns_list)} introns")

    # Step 3: Extract sequences for introns
    reporter.print_info("Extracting intron sequences from genome")
    sequence_extractor = SequenceExtractor(
        genome_file=str(config.input.genome),
        use_cache=True
    )
    introns_with_seq = sequence_extractor.extract_sequences(
        introns_list,
        flank_size=config.extraction.flank_len
    )
    # Materialize generator to list
    introns_all = list(introns_with_seq)
    logger.info(f"Extracted sequences for {len(introns_all)} introns")

    # Calculate actual minimum length needed for scoring regions
    # 5' region: relative to start, needs positions up to five_end
    # 3' region: relative to end, needs positions from three_start (negative)
    # BP region: relative to end, needs positions from bp_start (negative)
    min_for_five = abs(config.scoring.scoring_regions.five_end)
    min_for_three = abs(config.scoring.scoring_regions.three_start)
    min_for_bp = abs(config.scoring.scoring_regions.bp_start)

    # Intron must be long enough for all regions
    calculated_min = max(min_for_five + min_for_three, min_for_bp)
    actual_min_length = max(config.extraction.min_intron_len, calculated_min)

    # Filter by minimum length
    introns = [
        i for i in introns_all
        if i.sequences and len(i.sequences.seq) >= actual_min_length
    ]
    filtered_count = len(introns_all) - len(introns)
    if filtered_count > 0:
        logger.info(f"Filtered out {filtered_count} introns shorter than {actual_min_length}bp "
                   f"(required for scoring regions; user min was {config.extraction.min_intron_len}bp)")

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

    # Parse BED file
    parser = BEDParser()
    introns_no_seq = list(parser.parse_file(config.input.bed))
    logger.info(f"Parsed {len(introns_no_seq)} introns from BED")

    # Extract sequences
    reporter.print_info("Extracting sequences from genome")
    sequence_extractor = SequenceExtractor(
        genome_file=str(config.input.genome),
        use_cache=True
    )
    introns_with_seq = sequence_extractor.extract_sequences(
        introns_no_seq,
        flank_size=config.extraction.flank_len
    )
    # Materialize generator to list
    introns_all = list(introns_with_seq)
    logger.info(f"Extracted sequences for {len(introns_all)} introns")

    # Calculate actual minimum length needed for scoring regions
    # 5' region: relative to start, needs positions up to five_end
    # 3' region: relative to end, needs positions from three_start (negative)
    # BP region: relative to end, needs positions from bp_start (negative)
    min_for_five = abs(config.scoring.scoring_regions.five_end)
    min_for_three = abs(config.scoring.scoring_regions.three_start)
    min_for_bp = abs(config.scoring.scoring_regions.bp_start)

    # Intron must be long enough for all regions
    calculated_min = max(min_for_five + min_for_three, min_for_bp)
    actual_min_length = max(config.extraction.min_intron_len, calculated_min)

    # Filter by minimum length
    introns = [
        i for i in introns_all
        if i.sequences and len(i.sequences.seq) >= actual_min_length
    ]
    filtered_count = len(introns_all) - len(introns)
    if filtered_count > 0:
        logger.info(f"Filtered out {filtered_count} introns shorter than {actual_min_length}bp "
                   f"(required for scoring regions; user min was {config.extraction.min_intron_len}bp)")

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

    # Parse sequence file
    parser = SequenceParser()
    introns = list(parser.parse_file(config.input.sequence_file))
    logger.info(f"Loaded {len(introns)} introns from sequence file")

    return introns


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
    reporter.print_info("Loading PWM matrices")
    logger.info("Loading PWM matrices")

    # Load PWM matrices from data directory
    pwm_file = Path(__file__).parent.parent / "intronIC" / "data" / "scoring_matrices.fasta.iic"
    if not pwm_file.exists():
        raise FileNotFoundError(f"PWM file not found: {pwm_file}")

    pwm_sets = PWMLoader.load_from_file(pwm_file)
    logger.info(f"Loaded PWM matrices for {len(pwm_sets)} regions")

    # Load U2 BP matrix from separate file (fallback/conserved matrix)
    u2_bp_file = Path(__file__).parent.parent / "intronIC" / "data" / "u2.conserved_empirical_bp_pwm.iic"
    if u2_bp_file.exists():
        from dataclasses import replace
        u2_bp_matrices = PWMLoader.load_from_file(u2_bp_file)
        if 'bp' in u2_bp_matrices and u2_bp_matrices['bp'].u2_canonical:
            # PWMSet is frozen, so use replace() to create updated copy
            pwm_sets['bp'] = replace(
                pwm_sets['bp'],
                u2_canonical=u2_bp_matrices['bp'].u2_canonical
            )
            logger.info("Loaded conserved U2 BP matrix")
        else:
            logger.warning("U2 BP matrix file found but couldn't extract U2 canonical PWM")
    else:
        logger.warning(f"U2 BP matrix file not found: {u2_bp_file}")

    # Create scorer
    reporter.print_info("Calculating PWM scores")
    logger.info("Starting PWM scoring")

    scorer = IntronScorer(
        pwm_sets=pwm_sets,
        five_coords=(config.scoring.scoring_regions.five_start, config.scoring.scoring_regions.five_end),
        bp_coords=(config.scoring.scoring_regions.bp_start, config.scoring.scoring_regions.bp_end),
        three_coords=(config.scoring.scoring_regions.three_start, config.scoring.scoring_regions.three_end)
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
) -> Tuple[List[Intron], List[Intron], List[Intron]]:
    """Normalize intron scores using z-score transformation.

    Args:
        introns: List of introns with raw scores
        config: Pipeline configuration
        reporter: Progress reporter
        logger: Logger instance

    Returns:
        Tuple of (normalized experimental introns, u12 reference, u2 reference)
    """
    reporter.print_info("Loading reference sequences")
    logger.info("Loading reference sequences for normalization")

    # Load reference data
    data_dir = Path(__file__).parent.parent / "intronIC" / "data"
    u12_file = data_dir / "u12_reference.introns.iic.gz"
    u2_file = data_dir / "u2_reference.introns.iic.gz"

    if not u12_file.exists() or not u2_file.exists():
        raise FileNotFoundError(f"Reference data not found in {data_dir}")

    u12_reference = load_reference_sequences(u12_file, logger=logger)
    u2_reference = load_reference_sequences(u2_file, logger=logger)

    logger.info(f"Loaded {len(u12_reference)} U12 and {len(u2_reference)} U2 reference introns")

    # Score reference introns
    reporter.print_info("Scoring reference sequences")
    logger.info("Scoring reference sequences with PWMs")

    # Load PWM matrices
    pwm_file = data_dir / "scoring_matrices.fasta.iic"
    pwm_sets = PWMLoader.load_from_file(pwm_file)

    # Load U2 BP matrix from separate file (fallback/conserved matrix)
    from dataclasses import replace
    u2_bp_file = data_dir / "u2.conserved_empirical_bp_pwm.iic"
    if u2_bp_file.exists():
        u2_bp_matrices = PWMLoader.load_from_file(u2_bp_file)
        if 'bp' in u2_bp_matrices and u2_bp_matrices['bp'].u2_canonical:
            pwm_sets['bp'] = replace(
                pwm_sets['bp'],
                u2_canonical=u2_bp_matrices['bp'].u2_canonical
            )
            logger.info("Loaded conserved U2 BP matrix for reference scoring")

    scorer = IntronScorer(
        pwm_sets=pwm_sets,
        five_coords=(config.scoring.scoring_regions.five_start, config.scoring.scoring_regions.five_end),
        bp_coords=(config.scoring.scoring_regions.bp_start, config.scoring.scoring_regions.bp_end),
        three_coords=(config.scoring.scoring_regions.three_start, config.scoring.scoring_regions.three_end)
    )

    # Score reference introns
    u12_scored = [scorer.score_intron(i) for i in u12_reference]
    u2_scored = [scorer.score_intron(i) for i in u2_reference]
    reference_introns = u12_scored + u2_scored

    # Normalize scores
    reporter.print_info("Normalizing scores with z-score transformation")
    logger.info("Fitting normalizer on reference data")

    normalizer = ScoreNormalizer()
    normalizer.fit(reference_introns, dataset_type='reference')

    # Transform experimental introns
    normalized_introns = [normalizer.transform(i, dataset_type='experimental') for i in introns]

    logger.info(f"Normalized {len(normalized_introns)} experimental introns")
    return normalized_introns, u12_scored, u2_scored


def classify_introns(
    introns: List[Intron],
    u12_reference: List[Intron],
    u2_reference: List[Intron],
    config: IntronICConfig,
    reporter: IntronICProgressReporter,
    logger: logging.Logger
) -> Tuple[List[Intron], dict]:
    """Classify introns as U2 or U12 type.

    Args:
        introns: List of introns with z-scores
        u12_reference: U12 reference introns (scored and normalized)
        u2_reference: U2 reference introns (scored and normalized)
        config: Pipeline configuration
        reporter: Progress reporter
        logger: Logger instance

    Returns:
        Tuple of (classified introns, classification metrics)
    """
    reporter.print_info("Training SVM classifier")
    logger.info("Starting SVM classification")

    # Create classifier
    classifier = IntronClassifier(
        threshold=config.scoring.threshold,
        n_models=config.training.n_models,
        fixed_C=config.training.fixed_C,
        n_processes=config.performance.processes,
        cv_processes=config.performance.cv_processes,
        seed=config.training.seed
    )

    # Train on reference data
    reporter.print_info("Training ensemble SVM models")
    logger.info(f"Training on {len(u12_reference)} U12 and {len(u2_reference)} U2 reference introns")

    training_results = classifier.train(
        u12_introns=u12_reference,
        u2_introns=u2_reference
    )

    logger.info(f"Training complete. F1: {training_results['f1_score']:.3f}, "
                f"PR-AUC: {training_results['pr_auc']:.3f}")

    # Predict on experimental introns
    reporter.print_info("Classifying experimental introns")
    logger.info(f"Classifying {len(introns)} experimental introns")

    classified_introns = classifier.predict(introns)

    metrics = {
        'training_f1': training_results['f1_score'],
        'training_pr_auc': training_results['pr_auc']
    }

    logger.info("Classification complete")
    return classified_introns, metrics


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

    output_dir = config.output.output_dir
    base_name = config.output.base_filename

    # Write BED file
    bed_path = output_dir / f"{base_name}.bed.iic"
    logger.info(f"Writing BED file: {bed_path}")
    bed_writer = BEDWriter(bed_path)
    with bed_writer:
        for intron in introns:
            bed_writer.write_intron(intron)
    logger.info(f"Wrote {len(introns)} introns to BED file")

    # Write metadata file
    meta_path = output_dir / f"{base_name}.meta.iic"
    logger.info(f"Writing metadata file: {meta_path}")
    meta_writer = MetaWriter(meta_path)
    with meta_writer:
        for intron in introns:
            meta_writer.write_intron(intron)
    logger.info(f"Wrote metadata for {len(introns)} introns")

    # Write sequences file
    seq_path = output_dir / f"{base_name}.introns.iic"
    logger.info(f"Writing sequences file: {seq_path}")
    seq_writer = SequenceWriter(seq_path)
    with seq_writer:
        for intron in introns:
            seq_writer.write_intron(intron)
    logger.info(f"Wrote sequences for {len(introns)} introns")

    # Write score info file
    score_path = output_dir / f"{base_name}.score_info.iic"
    logger.info(f"Writing score info file: {score_path}")
    score_writer = ScoreWriter(score_path)
    with score_writer:
        for intron in introns:
            score_writer.write_intron(intron)
    logger.info(f"Wrote score info for {len(introns)} introns")

    output_files = {
        "Metadata": meta_path,
        "BED": bed_path,
        "Sequences": seq_path,
        "Scores": score_path,
        "Log": config.output.get_output_path('.log'),
    }

    reporter.print_file_tree(output_files)
    logger.info("All output files written successfully")


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
        normalized_introns, u12_reference, u2_reference = normalize_scores(
            scored_introns, config, reporter, logger
        )
        reporter.print_success("Scores normalized")

        # Step 5: Classify
        reporter.print_section("Step 5: Classify Introns", "bold blue")
        reporter.print_pipeline_steps(pipeline_steps, current_step=5)
        classified_introns, metrics = classify_introns(
            normalized_introns, u12_reference, u2_reference, config, reporter, logger
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
