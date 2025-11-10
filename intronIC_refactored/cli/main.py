"""
Main CLI entry point for intronIC.

Orchestrates the complete intron classification pipeline.
"""

import sys
import logging
import json
import time
import joblib
from pathlib import Path
from typing import List, Tuple

from .args import IntronICArgumentParser
from .config import IntronICConfig, ScoringRegions
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
from extraction.filters import IntronFilter
from scoring.pwm import PWMLoader
from scoring.scorer import IntronScorer
from scoring.normalizer import ScoreNormalizer
from classification.classifier import IntronClassifier
from visualization.plots import plot_classification_results
import gzip


def calculate_minimum_intron_length(
    scoring_regions: ScoringRegions,
    bp_matrix_length: int
) -> int:
    """Calculate minimum intron length needed for scoring.

    Port from: intronIC.py:4600-4607

    The minimum length is determined by how much of the intron is consumed
    by the scoring regions. The key insight: we need the actual BP matrix
    length (e.g., 12bp), not the search window size (e.g., 50bp).

    Args:
        scoring_regions: Scoring region coordinates
        bp_matrix_length: Length of branch point PWM matrix

    Returns:
        Minimum intron length in bp
    """
    # Calculate intronic positions in 5' region (positions >= 0)
    # Port from: intronIC.py:4601
    five_range = range(scoring_regions.five_start, scoring_regions.five_end)
    intronic_five = len([e for e in five_range if e >= 0])

    # Calculate intronic positions in 3' region (positions < 0)
    # Port from: intronIC.py:4602
    three_range = range(scoring_regions.three_start, scoring_regions.three_end)
    intronic_three_from_score = len([e for e in three_range if e < 0])

    # BP region needs: BP_MATRIX_LENGTH + distance from 3' end
    # Port from: intronIC.py:4598, 4603-4604
    bp_margin = abs(scoring_regions.bp_end)
    intronic_three_from_bp = bp_matrix_length + bp_margin

    # Use whichever is larger
    intronic_three = max(intronic_three_from_bp, intronic_three_from_score)

    # Total minimum = positions needed at 5' end + positions needed at 3' end
    # Port from: intronIC.py:4607
    minimum_length = intronic_five + intronic_three

    return minimum_length


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


def load_reference_sequences(filepath: Path, max_count: int = None, logger: logging.Logger = None) -> List[Intron]:
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
        logger: Optional logger instance

    Returns:
        List of Intron objects with sequences

    Note:
        No length filtering is applied here. Reference introns are filtered
        later using omit_check() after scoring regions are extracted, matching
        the original intronIC behavior.
    """
    introns = []

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

            # Create Intron (no length filtering - done later via omit_check)
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
    builder = AnnotationHierarchyBuilder(
        child_features=feature_types,
        clean_names=config.output.clean_names
    )
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

    # Step 3b: Apply U12 boundary correction to non-canonical introns (if enabled)
    # Port from: intronIC.py:2692 (u12_nc_ss_adjustment and u12_correction)
    # CRITICAL: This happens AFTER initial sequence extraction but BEFORE scoring
    # Corrected introns need sequence re-extraction with new coordinates
    if config.extraction.u12_boundary_correction:
        from extraction.boundary_correction import correct_intron_if_needed

        reporter.print_info("Checking non-canonical introns for U12 boundary corrections")
        corrected_count = 0
        corrected_introns = []

        for intron in introns_all:
            # Check and apply correction if needed
            corrected_intron, was_corrected = correct_intron_if_needed(
                intron,
                correction_enabled=True,
                use_strict_motif=True
            )

            # If corrected, re-extract sequences with new coordinates
            if was_corrected:
                # Re-extract sequences using corrected coordinates
                corrected_with_seq = sequence_extractor.extract_sequences(
                    [corrected_intron],
                    flank_size=config.extraction.flank_len
                )
                corrected_intron = list(corrected_with_seq)[0]
                corrected_count += 1
                logger.debug(
                    f"Corrected {intron.intron_id}: "
                    f"shift={corrected_intron.metadata.corrected}bp, "
                    f"new coords={corrected_intron.coordinates}"
                )

            corrected_introns.append(corrected_intron)

        introns_all = corrected_introns
        if corrected_count > 0:
            logger.info(f"Applied U12 boundary corrections to {corrected_count} non-canonical introns")
    else:
        logger.info("U12 boundary correction disabled (--no_nc_ss_adjustment)")

    # Load PWM matrices to get BP matrix length for minimum calculation
    # Port from: intronIC.py:4591-4592
    pwm_file = Path(__file__).parent.parent / "intronIC" / "data" / "scoring_matrices.fasta.iic"
    if not pwm_file.exists():
        raise FileNotFoundError(f"PWM file not found: {pwm_file}")

    pwm_sets = PWMLoader.load_from_file(pwm_file, pseudocount=config.scoring.pseudocount)
    # Get any U12 BP matrix to determine length (all should be same length)
    bp_matrix_length = next(iter(pwm_sets['bp'].matrices.values())).length
    logger.debug(f"BP matrix length: {bp_matrix_length}bp")

    # Calculate actual minimum length needed for scoring regions
    # Port from: intronIC.py:4600-4617
    calculated_min = calculate_minimum_intron_length(
        config.scoring.scoring_regions,
        bp_matrix_length
    )
    actual_min_length = max(config.extraction.min_intron_len, calculated_min)

    logger.info(f"Minimum intron length: {actual_min_length}bp "
               f"(user: {config.extraction.min_intron_len}bp, "
               f"scoring regions: {calculated_min}bp)")

    # Keep ALL introns at extraction time (matching original behavior)
    # Let IntronFilter decide what to omit during the filtering phase
    # This ensures short introns are marked as omitted and written to output files
    # NOTE: The original keeps all introns through extraction, then marks short ones
    # as omitted during filtering (see intronIC.py lines 4772-4786)
    introns = introns_all

    if config.scoring.sequences_only:
        logger.info("Sequences-only mode: all introns will be output regardless of length")
    else:
        logger.info(f"Extracted {len(introns):,} introns (length filtering during scoring filter phase)")

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

    # Apply U12 boundary correction (if enabled)
    if config.extraction.u12_boundary_correction:
        from extraction.boundary_correction import correct_intron_if_needed

        reporter.print_info("Checking non-canonical introns for U12 boundary corrections")
        corrected_count = 0
        corrected_introns = []

        for intron in introns_all:
            corrected_intron, was_corrected = correct_intron_if_needed(
                intron,
                correction_enabled=True,
                use_strict_motif=True
            )

            if was_corrected:
                corrected_with_seq = sequence_extractor.extract_sequences(
                    [corrected_intron],
                    flank_size=config.extraction.flank_len
                )
                corrected_intron = list(corrected_with_seq)[0]
                corrected_count += 1

            corrected_introns.append(corrected_intron)

        introns_all = corrected_introns
        if corrected_count > 0:
            logger.info(f"Applied U12 boundary corrections to {corrected_count} non-canonical introns")

    # Load PWM matrices to get BP matrix length for minimum calculation
    # Port from: intronIC.py:4591-4592
    pwm_file = Path(__file__).parent.parent / "intronIC" / "data" / "scoring_matrices.fasta.iic"
    if not pwm_file.exists():
        raise FileNotFoundError(f"PWM file not found: {pwm_file}")

    pwm_sets = PWMLoader.load_from_file(pwm_file, pseudocount=config.scoring.pseudocount)
    # Get any U12 BP matrix to determine length (all should be same length)
    bp_matrix_length = next(iter(pwm_sets['bp'].matrices.values())).length
    logger.debug(f"BP matrix length: {bp_matrix_length}bp")

    # Calculate actual minimum length needed for scoring regions
    # Port from: intronIC.py:4600-4617
    calculated_min = calculate_minimum_intron_length(
        config.scoring.scoring_regions,
        bp_matrix_length
    )
    actual_min_length = max(config.extraction.min_intron_len, calculated_min)

    logger.info(f"Minimum intron length: {actual_min_length}bp "
               f"(user: {config.extraction.min_intron_len}bp, "
               f"scoring regions: {calculated_min}bp)")

    # Keep ALL introns at extraction time (matching original behavior)
    # Let IntronFilter decide what to omit during the filtering phase
    # This ensures short introns are marked as omitted and written to output files
    # NOTE: The original keeps all introns through extraction, then marks short ones
    # as omitted during filtering (see intronIC.py lines 4772-4786)
    introns = introns_all

    if config.scoring.sequences_only:
        logger.info("Sequences-only mode: all introns will be output regardless of length")
    else:
        logger.info(f"Extracted {len(introns):,} introns (length filtering during scoring filter phase)")

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

    pwm_sets = PWMLoader.load_from_file(pwm_file, pseudocount=config.scoring.pseudocount)
    logger.info(f"Loaded PWM matrices for {len(pwm_sets)} regions")

    # Load U2 BP matrix from separate file (fallback/conserved matrix)
    u2_bp_file = Path(__file__).parent.parent / "intronIC" / "data" / "u2.conserved_empirical_bp_pwm.iic"
    if u2_bp_file.exists():
        from scoring.pwm import PWMSet
        u2_bp_matrices = PWMLoader.load_from_file(u2_bp_file, pseudocount=config.scoring.pseudocount)
        if 'bp' in u2_bp_matrices and u2_bp_matrices['bp'].u2_gtag:
            # PWMSet is frozen, so create new PWMSet with updated U2 GTAG matrix
            updated_matrices = dict(pwm_sets['bp'].matrices)
            updated_matrices[('u2', 'gtag')] = u2_bp_matrices['bp'].u2_gtag
            pwm_sets['bp'] = PWMSet(matrices=updated_matrices)
            logger.info("Loaded conserved U2 BP matrix")
        else:
            logger.warning("U2 BP matrix file found but couldn't extract U2 GTAG PWM")
    else:
        logger.warning(f"U2 BP matrix file not found: {u2_bp_file}")

    # Create scorer
    reporter.print_info("Calculating PWM scores")
    logger.info("Starting PWM scoring")

    scorer = IntronScorer(
        pwm_sets=pwm_sets,
        five_coords=(config.scoring.scoring_regions.five_start, config.scoring.scoring_regions.five_end),
        bp_coords=(config.scoring.scoring_regions.bp_start, config.scoring.scoring_regions.bp_end),
        three_coords=(config.scoring.scoring_regions.three_start, config.scoring.scoring_regions.three_end),
        ignore_nc_dnts=config.scoring.ignore_nc_dnts
    )

    # Score introns with error handling
    progress = reporter.create_progress()
    scored_introns = []
    failed_count = 0

    with progress:
        task = progress.add_task(
            "[cyan]Scoring introns...",
            total=len(introns)
        )

        for intron in introns:
            try:
                scored = scorer.score_intron(intron)
                scored_introns.append(scored)
            except Exception as e:
                # Log but continue - don't let one bad intron crash the pipeline
                logger.warning(
                    f"Failed to score intron {intron.intron_id}: {str(e)}. Skipping."
                )
                failed_count += 1

            progress.update(task, advance=1)

    logger.info(f"Scored {len(scored_introns)} introns successfully")
    if failed_count > 0:
        logger.warning(f"Failed to score {failed_count} introns (see warnings above)")
        reporter.print_warning(f"Skipped {failed_count} introns due to scoring errors")

    return scored_introns


def normalize_scores(
    introns: List[Intron],
    config: IntronICConfig,
    reporter: IntronICProgressReporter,
    logger: logging.Logger
) -> Tuple[List[Intron], List[Intron], List[Intron], 'ScoreNormalizer']:
    """Normalize intron scores using z-score transformation.

    Args:
        introns: List of introns with raw scores
        config: Pipeline configuration
        reporter: Progress reporter
        logger: Logger instance

    Returns:
        Tuple of (normalized experimental introns, u12 reference, u2 reference, normalizer)
    """
    reporter.print_info("Loading reference sequences")
    logger.info("Loading reference sequences for normalization")

    # Load reference data - use custom paths if provided, otherwise use defaults
    data_dir = Path(__file__).parent.parent / "intronIC" / "data"

    if config.scoring.reference_u12s:
        u12_file = config.scoring.reference_u12s
    else:
        u12_file = data_dir / "u12_reference.introns.iic.gz"

    if config.scoring.reference_u2s:
        u2_file = config.scoring.reference_u2s
    else:
        u2_file = data_dir / "u2_reference.introns.iic.gz"

    if not u12_file.exists() or not u2_file.exists():
        raise FileNotFoundError(
            f"Reference data not found. U12: {u12_file}, U2: {u2_file}"
        )

    u12_reference = load_reference_sequences(u12_file, logger=logger)
    u2_reference = load_reference_sequences(u2_file, logger=logger)

    logger.info(f"Loaded {len(u12_reference)} U12 and {len(u2_reference)} U2 reference introns")

    # Score reference introns
    reporter.print_info("Scoring reference sequences")
    logger.info("Scoring reference sequences with PWMs")

    # Load PWM matrices
    pwm_file = data_dir / "scoring_matrices.fasta.iic"
    pwm_sets = PWMLoader.load_from_file(pwm_file, pseudocount=config.scoring.pseudocount)

    # Load U2 BP matrix from separate file (fallback/conserved matrix)
    from scoring.pwm import PWMSet
    u2_bp_file = data_dir / "u2.conserved_empirical_bp_pwm.iic"
    if u2_bp_file.exists():
        u2_bp_matrices = PWMLoader.load_from_file(u2_bp_file, pseudocount=config.scoring.pseudocount)
        if 'bp' in u2_bp_matrices and u2_bp_matrices['bp'].u2_gtag:
            # PWMSet is frozen, so create new PWMSet with updated U2 GTAG matrix
            updated_matrices = dict(pwm_sets['bp'].matrices)
            updated_matrices[('u2', 'gtag')] = u2_bp_matrices['bp'].u2_gtag
            pwm_sets['bp'] = PWMSet(matrices=updated_matrices)
            logger.info("Loaded conserved U2 BP matrix for reference scoring")

    scorer = IntronScorer(
        pwm_sets=pwm_sets,
        five_coords=(config.scoring.scoring_regions.five_start, config.scoring.scoring_regions.five_end),
        bp_coords=(config.scoring.scoring_regions.bp_start, config.scoring.scoring_regions.bp_end),
        three_coords=(config.scoring.scoring_regions.three_start, config.scoring.scoring_regions.three_end),
        ignore_nc_dnts=config.scoring.ignore_nc_dnts
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

    # Transform reference introns (needed for classification)
    # Note: transform() takes an iterable and returns an iterator, so we materialize with list()
    u12_normalized = list(normalizer.transform(u12_scored, dataset_type='reference'))
    u2_normalized = list(normalizer.transform(u2_scored, dataset_type='reference'))

    # Transform experimental introns
    normalized_introns = list(normalizer.transform(introns, dataset_type='experimental'))

    logger.info(f"Normalized {len(normalized_introns)} experimental introns")
    logger.info(f"Normalized {len(u12_normalized)} U12 and {len(u2_normalized)} U2 reference introns")
    return normalized_introns, u12_normalized, u2_normalized, normalizer


def classify_with_pretrained_model(
    introns: List[Intron],
    model_path: Path,
    config: IntronICConfig,
    reporter: IntronICProgressReporter,
    logger: logging.Logger
) -> Tuple[List[Intron], dict]:
    """Classify introns using a pretrained model with cross-species domain adaptation.

    This function enables applying a trained SVM to new species without curated references.
    It implements unsupervised domain adaptation by fitting the normalizer on unlabeled
    experimental data, which is statistically valid because:
    - 99.5% of introns are U2-type → robust estimators learn species-specific U2 baseline
    - No label leakage (normalization uses only marginal feature distribution)
    - Corrects covariate shift from species-specific sequence characteristics

    Args:
        introns: Experimental introns (scored, not normalized)
        model_path: Path to pretrained model file (.model.pkl)
        config: Pipeline configuration
        reporter: Progress reporter
        logger: Logger instance

    Returns:
        Tuple of (classified introns, classification metrics)

    Note:
        The saved normalizer from training is NOT used. Instead, we fit a new normalizer
        on the experimental data to correct for species-specific score distributions.
    """
    reporter.print_info(f"Loading pretrained model from {model_path}")
    logger.info(f"Loading pretrained model from {model_path}")

    # Load model bundle
    if not model_path.exists():
        raise FileNotFoundError(f"Pretrained model not found: {model_path}")

    model_data = joblib.load(model_path)

    # Handle both old format (SVMEnsemble directly) and new format (dict bundle)
    if isinstance(model_data, dict):
        # New format: {'ensemble': ..., 'threshold': ...}
        ensemble = model_data['ensemble']
        saved_threshold = model_data.get('threshold', config.scoring.threshold)
        logger.info("Loaded model bundle (dict format)")
    else:
        # Old format: SVMEnsemble directly (backward compatibility)
        ensemble = model_data
        saved_threshold = config.scoring.threshold
        logger.info("Loaded model ensemble (legacy format - backward compatibility)")

    logger.info(f"Loaded ensemble with {len(ensemble.models)} models")
    logger.info(f"Using threshold: {config.scoring.threshold}")

    # Fit normalizer on experimental data (cross-species domain adaptation)
    # This is statistically valid for pretrained models:
    # - Corrects covariate shift from species-specific score distributions
    # - No label leakage (labels not used, only marginal feature distribution)
    # - RobustScaler minimally affected by rare U12s (~0.5% of data)
    reporter.print_info("Fitting normalizer on experimental data (domain adaptation)")
    logger.info(f"Fitting normalizer on {len(introns)} experimental introns")

    normalizer = ScoreNormalizer()
    normalizer.fit(introns, dataset_type='unlabeled')

    normalized_introns = list(normalizer.transform(introns, dataset_type='experimental'))
    logger.info(f"Normalized {len(normalized_introns)} experimental introns")

    # Classify using loaded ensemble
    reporter.print_info("Classifying with pretrained model")
    logger.info("Classifying introns with loaded ensemble")

    from classification.predictor import SVMPredictor
    predictor = SVMPredictor(
        threshold=config.scoring.threshold,  # Use config threshold (can override saved)
        n_jobs=config.performance.processes
    )

    classified_introns = list(predictor.predict(ensemble, normalized_introns))
    logger.info(f"Classified {len(classified_introns)} introns")

    # Create metrics (limited since we skipped training)
    metrics = {
        'optimized_C': 'N/A (pretrained)',
        'cv_score': 'N/A (pretrained)',
        'n_models': len(ensemble.models),
        'pretrained': True,
        'model_path': str(model_path)
    }

    return classified_introns, metrics


def classify_introns(
    introns: List[Intron],
    u12_reference: List[Intron],
    u2_reference: List[Intron],
    normalizer: 'ScoreNormalizer',
    config: IntronICConfig,
    reporter: IntronICProgressReporter,
    logger: logging.Logger
) -> Tuple[List[Intron], dict]:
    """Classify introns as U2 or U12 type.

    Args:
        introns: List of introns with z-scores
        u12_reference: U12 reference introns (scored and normalized)
        u2_reference: U2 reference introns (scored and normalized)
        normalizer: Fitted normalizer used for score transformation
        config: Pipeline configuration
        reporter: Progress reporter
        logger: Logger instance

    Returns:
        Tuple of (classified introns, classification metrics)
    """
    reporter.print_info("Training SVM classifier")
    logger.info("Starting SVM classification")

    # Create classifier with correct parameter names
    # IntronClassifier API uses:
    # - classification_threshold (not threshold)
    # - n_ensemble_models (not n_models)
    # - fixed_c (not fixed_C)
    # - random_state (not seed)
    # - cv_processes (for cross-validation parallelization)
    # - classification_processes (for prediction parallelization)
    classifier = IntronClassifier(
        classification_threshold=config.scoring.threshold,
        n_ensemble_models=config.training.n_models,
        fixed_c=config.training.fixed_C,
        optimize_c=(config.training.fixed_C is None),
        random_state=config.training.seed,
        cv_processes=config.performance.cv_processes,
        classification_processes=config.performance.processes,
        max_iter=config.training.max_iter,
        eval_mode=config.training.eval_mode,
        n_cv_folds=config.training.n_cv_folds,
        test_fraction=config.training.test_fraction
    )

    # Run complete classification pipeline (optimize + train + classify)
    logger.info(f"Running classification on {len(introns)} experimental introns")
    logger.info(f"Reference data: {len(u12_reference)} U12, {len(u2_reference)} U2")

    result = classifier.classify(
        u12_reference=u12_reference,
        u2_reference=u2_reference,
        experimental=introns
    )

    # Extract metrics from result
    metrics = {
        'optimized_C': result.parameters.C,
        'cv_score': result.parameters.cv_score,
        'n_models': len(result.ensemble.models)
    }

    # Add evaluation metrics if nested CV or split eval was performed
    if result.eval_result is not None:
        # Check if it's NestedCVResult (has mean_f1) or SplitEvalResult (has f1_score)
        if hasattr(result.eval_result, 'mean_f1'):
            # Nested CV result
            metrics['mean_f1'] = result.eval_result.mean_f1
            metrics['std_f1'] = result.eval_result.std_f1
            metrics['mean_pr_auc'] = result.eval_result.mean_pr_auc
            metrics['std_pr_auc'] = result.eval_result.std_pr_auc
            metrics['n_cv_folds'] = result.eval_result.n_folds
        elif hasattr(result.eval_result, 'f1_score'):
            # Split evaluation result
            metrics['f1'] = result.eval_result.f1_score
            metrics['pr_auc'] = result.eval_result.pr_auc

    logger.info(f"Classification complete")
    logger.info(f"  Optimized C: {metrics['optimized_C']:.6e}")
    logger.info(f"  Models trained: {metrics['n_models']}")

    # Log evaluation metrics if available
    if 'mean_f1' in metrics:
        logger.info(f"  Nested CV results ({metrics['n_cv_folds']} folds):")
        logger.info(f"    Mean F1: {metrics['mean_f1']:.4f} ± {metrics['std_f1']:.4f}")
        logger.info(f"    Mean PR-AUC: {metrics['mean_pr_auc']:.4f} ± {metrics['std_pr_auc']:.4f}")
    elif 'f1' in metrics:
        logger.info(f"  Test set evaluation:")
        logger.info(f"    F1: {metrics['f1']:.4f}")
        logger.info(f"    PR-AUC: {metrics['pr_auc']:.4f}")

    # Save trained model (ensemble only - normalizer fitted per-species during inference)
    model_path = config.output.get_output_path('.model.pkl')
    logger.info(f"Saving trained model to {model_path}")
    model_bundle = {
        'ensemble': result.ensemble,
        'threshold': config.scoring.threshold
    }
    joblib.dump(model_bundle, model_path, compress=3)
    logger.info(f"Model saved successfully")

    return list(result.classified_introns), metrics


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

    # Filter duplicates if not including them
    # Port from: intronIC.py uses -d/--include_duplicates flag
    if not config.extraction.include_duplicates:
        original_count = len(introns)
        introns = [i for i in introns if not (i.metadata and i.metadata.duplicate)]
        filtered_count = original_count - len(introns)
        if filtered_count > 0:
            logger.info(f"Filtered out {filtered_count} duplicate introns (use -d to include)")

    output_dir = config.output.output_dir
    base_name = config.output.base_filename
    species_name = config.output.species_name

    # Write BED file
    bed_path = output_dir / f"{base_name}.bed.iic"
    logger.info(f"Writing BED file: {bed_path}")
    bed_writer = BEDWriter(bed_path)
    with bed_writer:
        for intron in introns:
            bed_writer.write_intron(intron, species_name=species_name)
    logger.info(f"Wrote {len(introns)} introns to BED file")

    # Write metadata file
    meta_path = output_dir / f"{base_name}.meta.iic"
    logger.info(f"Writing metadata file: {meta_path}")
    meta_writer = MetaWriter(meta_path)
    with meta_writer:
        for intron in introns:
            meta_writer.write_intron(intron, species_name=species_name)
    logger.info(f"Wrote metadata for {len(introns)} introns")

    # Write sequences file
    seq_path = output_dir / f"{base_name}.introns.iic"
    logger.info(f"Writing sequences file: {seq_path}")
    seq_writer = SequenceWriter(seq_path)
    with seq_writer:
        for intron in introns:
            seq_writer.write_intron(intron, species_name=species_name)
    logger.info(f"Wrote sequences for {len(introns)} introns")

    # Write score info file
    score_path = output_dir / f"{base_name}.score_info.iic"
    logger.info(f"Writing score info file: {score_path}")
    score_writer = ScoreWriter(score_path)
    with score_writer:
        for intron in introns:
            score_writer.write_intron(intron, species_name=species_name)
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
    # Track start time for runtime reporting
    start_time = time.time()

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

        # Filter introns before scoring (duplicates, short introns, longest isoform)
        # This matches original intronIC behavior where filtering happens BEFORE scoring
        # to avoid scoring 5x more introns than necessary (which causes O(n²) slowdown)
        reporter.print_info("Filtering introns before scoring")
        logger.info("Applying pre-scoring filters (duplicates, omissions, longest isoform)")

        # Create filter with scoring-appropriate settings:
        # - longest_only=True: Only score longest isoform per gene (filters ~8k introns)
        # - include_duplicates=False: Don't score duplicates (filters ~38k introns)
        # - min_length: Filter short introns
        # - allow_noncanonical: Based on exclude_noncanonical flag
        # - allow_overlap: Based on no_intron_overlap flag
        intron_filter = IntronFilter(
            min_length=config.extraction.min_intron_len,
            bp_matrix_length=7,  # Default from original
            scoring_regions=['five', 'three'],  # Check these for ambiguous bases
            allow_noncanonical=not config.scoring.exclude_noncanonical,
            allow_overlap=not config.extraction.no_intron_overlap,
            longest_only=True,  # ALWAYS filter to longest isoform for scoring
            include_duplicates=False  # ALWAYS filter duplicates for scoring
        )

        filtered_introns = intron_filter.filter_introns(introns)

        # Report filtering statistics
        stats = intron_filter.stats
        logger.info(
            f"Filtering results: {stats.kept_introns:,}/{stats.total_introns:,} introns "
            f"kept for scoring"
        )
        logger.info(
            f"Omitted: {stats.omitted_short} short, {stats.omitted_ambiguous} ambiguous, "
            f"{stats.omitted_noncanonical} non-canonical, {stats.omitted_isoform} non-longest isoform, "
            f"{stats.omitted_overlap} overlapping"
        )
        logger.info(f"Duplicates marked: {stats.duplicates}")

        reporter.print_success(
            f"Filtered to {len(filtered_introns):,} introns for scoring "
            f"(removed {len(introns) - len(filtered_introns):,})"
        )

        # Important: Use filtered_introns for scoring, but keep original introns list
        # for potential output (user may want duplicates via -d flag)
        introns_for_scoring = filtered_introns

        # Step 3: Score introns
        reporter.print_section("Step 3: Score Introns", "bold blue")
        reporter.print_pipeline_steps(pipeline_steps, current_step=3)
        scored_introns = score_introns(introns_for_scoring, config, reporter, logger)
        reporter.print_success(f"Scored {len(scored_introns):,} introns")

        # Check if using pretrained model
        if config.training.pretrained_model_path:
            # Skip normalization/training - use pretrained model
            reporter.print_section("Step 4-5: Classify with Pretrained Model", "bold blue")
            reporter.print_pipeline_steps(pipeline_steps, current_step=4)
            classified_introns, metrics = classify_with_pretrained_model(
                scored_introns, config.training.pretrained_model_path, config, reporter, logger
            )
        else:
            # Normal flow: normalize + train + classify
            # Step 4: Normalize scores
            reporter.print_section("Step 4: Normalize Scores", "bold blue")
            reporter.print_pipeline_steps(pipeline_steps, current_step=4)
            normalized_introns, u12_reference, u2_reference, normalizer = normalize_scores(
                scored_introns, config, reporter, logger
            )
            reporter.print_success("Scores normalized")

            # Step 5: Classify
            reporter.print_section("Step 5: Classify Introns", "bold blue")
            reporter.print_pipeline_steps(pipeline_steps, current_step=5)
            classified_introns, metrics = classify_introns(
                normalized_introns, u12_reference, u2_reference, normalizer, config, reporter, logger
            )

        # Count classifications
        u12_count = sum(
            1 for i in classified_introns
            if i.metadata and i.metadata.type_id == 'u12'
        )
        u2_count = len(classified_introns) - u12_count

        # Count AT-AC introns (characteristic U12 boundaries)
        atac_count = sum(
            1 for i in classified_introns
            if (i.metadata and i.metadata.type_id == 'u12' and
                i.sequences and i.sequences.terminal_dinucleotides == 'AT-AC')
        )

        reporter.print_classification_summary(
            total=len(classified_introns),
            u12_count=u12_count,
            u2_count=u2_count,
            threshold=config.scoring.threshold
        )

        # Log classification summary to log file
        logger.info(f"Classification results:")
        logger.info(f"  {atac_count} putative AT-AC U12-type introns found")
        logger.info(f"  {u12_count} putative U12-type introns found with scores > {config.scoring.threshold}%")
        logger.info(f"  {u2_count} introns classified as U2-type")

        # Collect and log non-canonical boundary statistics
        from collections import Counter
        nc_types = Counter()
        for intron in classified_introns:
            if (intron.metadata and intron.metadata.noncanonical and
                intron.sequences and intron.sequences.terminal_dinucleotides):
                nc_types[intron.sequences.terminal_dinucleotides] += 1

        if nc_types:
            logger.info("Most common non-canonical splice sites:")
            total_nc = sum(nc_types.values())
            for dnts, count in nc_types.most_common(10):  # Top 10
                percentage = (count / total_nc) * 100
                logger.info(f"  * {dnts} ({count}/{total_nc}, {percentage:.2f}%)")

        # Save classification metrics to JSON file
        if metrics:
            metrics_path = config.output.get_output_path('.metrics.json')
            logger.info(f"Saving classification metrics to {metrics_path}")
            with open(metrics_path, 'w') as f:
                json.dump(metrics, f, indent=2)

        # Generate visualization plots
        logger.info("Generating visualization plots")
        try:
            plot_classification_results(
                introns=classified_introns,
                output_dir=config.output.output_dir,
                species_name=config.output.base_filename,
                threshold=config.scoring.threshold,
                fig_dpi=300
            )
            logger.info("Successfully generated classification plots")
        except Exception as plot_error:
            logger.warning(f"Failed to generate plots: {plot_error}")
            # Continue even if plotting fails

        # Step 6: Write outputs
        reporter.print_section("Step 6: Write Outputs", "bold blue")
        reporter.print_pipeline_steps(pipeline_steps, current_step=6)
        write_outputs(classified_introns, config, reporter, logger)

        # Calculate and log total runtime
        elapsed_seconds = time.time() - start_time
        hours, remainder = divmod(int(elapsed_seconds), 3600)
        minutes, seconds = divmod(remainder, 60)

        if hours > 0:
            runtime_str = f"{hours}h {minutes}m {seconds}s"
        elif minutes > 0:
            runtime_str = f"{minutes}m {seconds}s"
        else:
            runtime_str = f"{seconds}s"

        logger.info(f"Run finished in {runtime_str}")

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
