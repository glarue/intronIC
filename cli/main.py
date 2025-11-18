"""
Main CLI entry point for intronIC.

Orchestrates the complete intron classification pipeline.
"""

import sys
import logging
import json
import time
import gc
import joblib
import numpy as np
from pathlib import Path
from typing import List, Tuple, Optional
from multiprocessing import Pool
from itertools import repeat
from dataclasses import replace

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
from smart_open import open as smart_open
from utils.metadata import generate_training_metadata, generate_pretrained_metadata, write_metadata


# ============================================================================
# PHASE 3: Parallel Scoring Worker Function
# ============================================================================
# This function must be at module level (not nested) to be picklable for multiprocessing


def _score_intron_worker(
    intron: 'Intron',
    pwm_file: Path,
    u2_bp_file: Path,
    five_coords: Tuple[int, int],
    bp_coords: Tuple[int, int],
    three_coords: Tuple[int, int],
    ignore_nc_dnts: bool,
    pseudocount: float
) -> Tuple[Optional['Intron'], Optional[str]]:
    """
    Worker function for parallel intron scoring.

    Each worker loads PWMs (OS caches file reads for efficiency) and creates
    its own scorer instance to avoid pickling issues.

    Args:
        intron: Intron to score
        pwm_file: Path to PWM matrices file
        u2_bp_file: Path to U2 branch point matrix file
        five_coords: 5' splice site scoring region coordinates
        bp_coords: Branch point search region coordinates
        three_coords: 3' splice site scoring region coordinates
        ignore_nc_dnts: Whether to ignore non-canonical dinucleotides
        pseudocount: Pseudocount for PWM scoring

    Returns:
        (scored_intron, error_message) tuple
        If scoring succeeds: (intron, None)
        If scoring fails: (None, error_message)
    """
    try:
        # Import here to avoid issues with multiprocessing pickling
        from scoring.scorer import IntronScorer
        from scoring.pwm import PWMLoader, PWMSet
        from pathlib import Path

        # Load PWMs (OS file system cache makes this efficient across workers)
        pwm_sets = PWMLoader.load_from_file(pwm_file, pseudocount=pseudocount)

        # Load U2 BP matrix if available
        if u2_bp_file.exists():
            u2_bp_matrices = PWMLoader.load_from_file(u2_bp_file, pseudocount=pseudocount)
            if 'bp' in u2_bp_matrices and u2_bp_matrices['bp'].u2_gtag:
                # Update BP PWM set with U2 GTAG matrix
                updated_matrices = dict(pwm_sets['bp'].matrices)
                updated_matrices[('u2', 'gtag')] = u2_bp_matrices['bp'].u2_gtag
                pwm_sets['bp'] = PWMSet(matrices=updated_matrices)

        # Create scorer for this worker
        scorer = IntronScorer(
            pwm_sets=pwm_sets,
            five_coords=five_coords,
            bp_coords=bp_coords,
            three_coords=three_coords,
            ignore_nc_dnts=ignore_nc_dnts
        )

        # Score intron
        scored = scorer.score_intron(intron)
        return (scored, None)

    except Exception as e:
        # Return error instead of crashing worker
        # Include intron ID in error message for debugging
        intron_id = intron.intron_id if hasattr(intron, 'intron_id') else 'unknown'
        return (None, f"{intron_id}: {str(e)}")


def log_data_block(logger: logging.Logger, header: str, lines: List[str], use_separator: bool = True):
    """
    Log a block of data with a header as a single multi-line message.

    This ensures only the header gets a timestamp, while the data lines
    appear without individual timestamps for cleaner log formatting.

    Args:
        logger: Logger instance
        header: Header text (will get timestamp)
        lines: List of data lines to display
        use_separator: Whether to include separator lines (default: True)

    Example:
        log_data_block(
            logger,
            "Top 20 splice site boundaries (U12-type introns)",
            ["   1. GT-AG    11,684 (97.02%)", "   2. GC-AG       158 ( 1.31%)"]
        )

        Results in:
        [2025-11-13 00:18:48] INFO     Top 20 splice site boundaries (U12-type introns)
        ------------------------------------...
           1. GT-AG    11,684 (97.02%)
           2. GC-AG       158 ( 1.31%)
    """
    # Build the complete message starting with header
    parts = [header]

    if use_separator:
        parts.append("-" * 100)

    parts.extend(lines)

    # Log as single multi-line message (only first line gets timestamp)
    logger.info('\n'.join(parts))


def clear_large_sequences_for_classification(introns: List[Intron]) -> List[Intron]:
    """
    Clear large sequence fields before classification to reduce memory.

    After scoring completes, the large sequence fields (seq, upstream_flank,
    downstream_flank, bp_region_seq) are no longer needed for classification.
    Only the small scored sequences (five_seq, three_seq, bp_seq, bp_seq_u2)
    and terminal dinucleotides are needed.

    This function creates new intron objects with large sequences cleared,
    reducing memory usage by ~10-12 GB for 1M introns (from ~14 GB to ~2 GB).

    Clears:
    - seq (full intron sequence, ~500 bytes avg)
    - upstream_flank (exonic context, ~200 bytes avg)
    - downstream_flank (exonic context, ~200 bytes avg)
    - bp_region_seq (branch point search region, ~50 bytes avg)

    Keeps:
    - five_seq, three_seq, bp_seq, bp_seq_u2 (scored sequences, ~40 bytes total)
    - five_prime_dnt, three_prime_dnt (terminal dinucleotides, ~4 bytes total)
    - bp_relative_coords (branch point position, ~16 bytes)
    - All scores and metadata

    Args:
        introns: List of scored introns with full sequences

    Returns:
        New list of introns with large sequences cleared (functional style)
    """
    from dataclasses import replace

    cleared = []
    for intron in introns:
        # Create new intron with cleared sequences
        # Uses functional style - returns new object, original unchanged
        new_sequences = replace(
            intron.sequences,
            seq=None,
            upstream_flank=None,
            downstream_flank=None,
            bp_region_seq=None
        )
        cleared_intron = replace(intron, sequences=new_sequences)
        cleared.append(cleared_intron)

    return cleared


def format_count_with_percentage(count: int, total: int) -> str:
    """Format a count with percentage of total.

    Args:
        count: The count to format
        total: The total to calculate percentage from

    Returns:
        Formatted string like "1,234 (5.67%)" or "0 (0.00%)"

    Examples:
        >>> format_count_with_percentage(12074, 58933)
        '12,074 (20.49%)'
        >>> format_count_with_percentage(31, 12074)
        '31 (0.26%)'
    """
    if total == 0:
        percentage = 0.0
    else:
        percentage = (count / total) * 100
    return f"{count:,} ({percentage:.2f}%)"


def merge_scored_and_omitted_introns(
    scored_introns: List[Intron],
    all_introns: List[Intron],
    messenger: 'UnifiedMessenger'
) -> List[Intron]:
    """
    Merge scored introns with unique omitted introns for complete meta output.

    Matches original intronIC behavior where .meta.iic includes both:
    - Scored introns (with classification results)
    - Omitted introns (with [o:X] tags but no scores)

    Duplicates are excluded from output (they're already filtered).

    Args:
        scored_introns: Introns that went through scoring/classification (have scores)
        all_introns: All extracted introns (including omitted)
        messenger: Unified messenger instance

    Returns:
        Combined list: scored introns + unique omitted introns (no duplicates)
    """
    # Create set of scored intron IDs for fast lookup
    scored_ids = {id(intron) for intron in scored_introns}

    # Find omitted introns that aren't duplicates
    # These should have metadata.omitted set and NOT be duplicates
    omitted_introns = [
        intron for intron in all_introns
        if id(intron) not in scored_ids
        and intron.metadata
        and intron.metadata.omitted
        and not intron.metadata.duplicate
    ]

    total_output = len(scored_introns) + len(omitted_introns)
    total_all = len(all_introns)

    messenger.log_only(
        f"Merging output: {format_count_with_percentage(len(scored_introns), total_all)} scored + "
        f"{format_count_with_percentage(len(omitted_introns), total_all)} omitted = "
        f"{format_count_with_percentage(total_output, total_all)} total introns for output files"
    )

    # Return scored + omitted (duplicates already excluded)
    return scored_introns + omitted_introns


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


def setup_logging(config: IntronICConfig) -> tuple[logging.Logger, 'Console']:
    """Setup logging configuration with ANSI color support.

    Args:
        config: Pipeline configuration

    Returns:
        Tuple of (configured logger instance, Rich console for log file)
    """
    from rich.console import Console
    from rich.logging import RichHandler

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

    # Clear any existing handlers
    logger.handlers.clear()

    # Create Rich console for log file with ANSI color support
    # force_terminal=True ensures ANSI codes are written even to files
    log_console = Console(
        file=open(log_file, 'w', encoding='utf-8'),
        force_terminal=True,
        width=120,
        legacy_windows=False
    )

    # File handler using Rich - preserves colors and formatting
    # NOTE: We only add a file handler. Console output is handled by UnifiedMessenger
    # to avoid duplicate messages. The logger is used ONLY for the log file.
    file_handler = RichHandler(
        console=log_console,
        show_time=True,
        show_level=True,
        show_path=False,
        markup=True,  # Enable Rich markup interpretation for proper ANSI formatting
        rich_tracebacks=True,
        tracebacks_show_locals=config.output.debug,
        level=logging.DEBUG  # Always log everything to file
    )
    logger.addHandler(file_handler)

    # NO console handler - UnifiedMessenger handles console output to avoid duplication

    return logger, log_console


def load_reference_sequences(filepath: Path, max_count: int = None, messenger: 'UnifiedMessenger' = None) -> List[Intron]:
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
        messenger: Optional messenger instance

    Returns:
        List of Intron objects with sequences

    Note:
        No length filtering is applied here. Reference introns are filtered
        later using omit_check() after scoring regions are extracted, matching
        the original intronIC behavior.
    """
    introns = []

    with smart_open(filepath, 'rt') as f:
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

    if messenger:
        messenger.log_only(f"Loaded {len(introns)} reference sequences from {filepath.name}")

    return introns


def load_genome(config: IntronICConfig, messenger: 'UnifiedMessenger') -> GenomeReader:
    """Load genome file.

    Args:
        config: Pipeline configuration
        messenger: Unified messenger for output

    Returns:
        GenomeReader instance
    """
    messenger.info(f"Loading genome: {config.input.genome}")
    # Use cached mode for faster repeated access
    reader = GenomeReader(config.input.genome, cached=True)
    if reader.cache:
        messenger.info(f"Loaded {len(reader.cache)} sequences into memory")
    return reader


def extract_introns_from_annotation(
    config: IntronICConfig,
    genome_reader: GenomeReader,
    messenger: 'UnifiedMessenger',
    reporter: IntronICProgressReporter
) -> List[Intron]:
    """Extract introns from annotation file.

    Args:
        config: Pipeline configuration
        genome_reader: Genome reader instance
        messenger: Unified messenger for output
        reporter: Progress reporter (for progress bars only)

    Returns:
        List of extracted introns
    """
    messenger.info(f"Parsing annotation: {config.input.annotation}")

    # Step 1: Build annotation hierarchy
    # CRITICAL: Always extract BOTH cds and exon features (order matters: CDS first!)
    # This matches original intronIC behavior (v1.5.1 line 1859-1865):
    # - Extract from both feature types to get consistent family_size
    # - Filter afterward based on user's --feature preference
    # - CDS must come first for proper prioritization in IntronGenerator
    builder = AnnotationHierarchyBuilder(
        child_features=['cds', 'exon'],  # Always both, order matters!
        clean_names=config.output.clean_names
    )
    genes = builder.build_from_file(config.input.annotation)
    messenger.log_only(f"Built hierarchy with {len(genes)} genes")

    # Step 2: Generate introns from exon pairs
    messenger.info("Generating introns from exon pairs")
    generator = IntronGenerator()
    introns_iter = generator.generate_from_genes(genes, builder.feature_index)
    introns_all = list(introns_iter)  # Materialize iterator
    messenger.log_only(f"Generated {len(introns_all)} introns from both CDS and exon features")

    # Step 2b: Filter introns based on --feature flag
    # Port from: intronIC.py:1865 - filter by defined_by after extraction
    # This ensures family_size is consistent regardless of feature type preference
    if config.extraction.feature_type == 'cds':
        # User wants CDS-defined introns only
        introns_list = [i for i in introns_all if i.metadata.defined_by == 'cds']
        messenger.log_only(f"Filtered to {len(introns_list)} CDS-defined introns (--feature cds)")
    elif config.extraction.feature_type == 'exon':
        # User wants exon-defined introns only
        introns_list = [i for i in introns_all if i.metadata.defined_by == 'exon']
        messenger.log_only(f"Filtered to {len(introns_list)} exon-defined introns (--feature exon)")
    else:
        # Default 'both': keep all introns (CDS-defined + exon-only)
        introns_list = introns_all
        cds_count = sum(1 for i in introns_all if i.metadata.defined_by == 'cds')
        exon_count = sum(1 for i in introns_all if i.metadata.defined_by == 'exon')
        messenger.log_only(
            f"Keeping all introns: {cds_count} CDS-defined + {exon_count} exon-only (--feature both)"
        )

    # Step 3: Extract sequences for introns
    messenger.info("Extracting intron sequences from genome")
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
    messenger.log_only(f"Extracted sequences for {len(introns_all)} introns")

    # Free memory from coordinate-only list (no longer needed)
    del introns_list
    gc.collect()

    # Step 3b: Apply U12 boundary correction to non-canonical introns (if enabled)
    # Port from: intronIC.py:2692 (u12_nc_ss_adjustment and u12_correction)
    # CRITICAL: This happens AFTER initial sequence extraction but BEFORE scoring
    # Corrected introns need sequence re-extraction with new coordinates
    if config.extraction.u12_boundary_correction:
        from extraction.boundary_correction import correct_intron_if_needed

        messenger.info("Checking non-canonical introns for U12-type boundary corrections")
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
                messenger.log_only(
                    f"Corrected {intron.intron_id}: "
                    f"shift={corrected_intron.metadata.correction_distance}bp, "
                    f"new coords={corrected_intron.coordinates}",
                    level="debug"
                )

            corrected_introns.append(corrected_intron)

        introns_all = corrected_introns
        if corrected_count > 0:
            total_introns = len(introns_all)
            messenger.log_only(
                f"Applied U12-type boundary corrections to "
                f"{format_count_with_percentage(corrected_count, total_introns)} non-canonical introns"
            )
    else:
        messenger.log_only("U12-type boundary correction disabled (--no_nc_ss_adjustment)")

    # Load PWM matrices to get BP matrix length for minimum calculation
    # Port from: intronIC.py:4591-4592
    pwm_file = Path(__file__).parent.parent / "data" / "scoring_matrices.fasta.iic"
    if not pwm_file.exists():
        raise FileNotFoundError(f"PWM file not found: {pwm_file}")

    pwm_sets = PWMLoader.load_from_file(pwm_file, pseudocount=config.scoring.pseudocount)
    # Get any U12 BP matrix to determine length (all should be same length)
    bp_matrix_length = next(iter(pwm_sets['bp'].matrices.values())).length
    messenger.log_only(f"BP matrix length: {bp_matrix_length}bp", level="debug")

    # Calculate actual minimum length needed for scoring regions
    # Port from: intronIC.py:4600-4617
    calculated_min = calculate_minimum_intron_length(
        config.scoring.scoring_regions,
        bp_matrix_length
    )
    actual_min_length = max(config.extraction.min_intron_len, calculated_min)

    messenger.log_only(
        f"Minimum intron length: {actual_min_length}bp "
        f"(user: {config.extraction.min_intron_len}bp, "
        f"scoring regions: {calculated_min}bp)"
    )

    # Keep ALL introns at extraction time (matching original behavior)
    # Let IntronFilter decide what to omit during the filtering phase
    # This ensures short introns are marked as omitted and written to output files
    # NOTE: The original keeps all introns through extraction, then marks short ones
    # as omitted during filtering (see intronIC.py lines 4772-4786)
    introns = introns_all

    if config.scoring.sequences_only:
        messenger.log_only("Sequences-only mode: all introns will be output regardless of length")
    else:
        messenger.log_only(f"Extracted {len(introns):,} introns (length filtering during scoring filter phase)")

    return introns


def extract_introns_from_bed(
    config: IntronICConfig,
    genome_reader: GenomeReader,
    messenger: 'UnifiedMessenger',
    reporter: IntronICProgressReporter
) -> List[Intron]:
    """Extract introns from BED file.

    Args:
        config: Pipeline configuration
        genome_reader: Genome reader instance
        messenger: Unified messenger for output
        reporter: Progress reporter (for progress bars only)

    Returns:
        List of extracted introns
    """
    messenger.info(f"Reading BED file: {config.input.bed}")

    # Parse BED file
    parser = BEDParser()
    introns_no_seq = list(parser.parse_file(config.input.bed))
    messenger.log_only(f"Parsed {len(introns_no_seq)} introns from BED")

    # Extract sequences
    messenger.info("Extracting sequences from genome")
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
    messenger.log_only(f"Extracted sequences for {len(introns_all)} introns")

    # Free memory from coordinate-only list (no longer needed)
    del introns_no_seq
    gc.collect()

    # Apply U12 boundary correction (if enabled)
    if config.extraction.u12_boundary_correction:
        from extraction.boundary_correction import correct_intron_if_needed

        messenger.info("Checking non-canonical introns for U12-type boundary corrections")
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
            total_introns = len(introns_all)
            messenger.log_only(
                f"Applied U12-type boundary corrections to "
                f"{format_count_with_percentage(corrected_count, total_introns)} non-canonical introns"
            )

    # Load PWM matrices to get BP matrix length for minimum calculation
    # Port from: intronIC.py:4591-4592
    pwm_file = Path(__file__).parent.parent / "data" / "scoring_matrices.fasta.iic"
    if not pwm_file.exists():
        raise FileNotFoundError(f"PWM file not found: {pwm_file}")

    pwm_sets = PWMLoader.load_from_file(pwm_file, pseudocount=config.scoring.pseudocount)
    # Get any U12 BP matrix to determine length (all should be same length)
    bp_matrix_length = next(iter(pwm_sets['bp'].matrices.values())).length
    messenger.log_only(f"BP matrix length: {bp_matrix_length}bp", level="debug")

    # Calculate actual minimum length needed for scoring regions
    # Port from: intronIC.py:4600-4617
    calculated_min = calculate_minimum_intron_length(
        config.scoring.scoring_regions,
        bp_matrix_length
    )
    actual_min_length = max(config.extraction.min_intron_len, calculated_min)

    messenger.log_only(
        f"Minimum intron length: {actual_min_length}bp "
        f"(user: {config.extraction.min_intron_len}bp, "
        f"scoring regions: {calculated_min}bp)"
    )

    # Keep ALL introns at extraction time (matching original behavior)
    # Let IntronFilter decide what to omit during the filtering phase
    # This ensures short introns are marked as omitted and written to output files
    # NOTE: The original keeps all introns through extraction, then marks short ones
    # as omitted during filtering (see intronIC.py lines 4772-4786)
    introns = introns_all

    if config.scoring.sequences_only:
        messenger.log_only("Sequences-only mode: all introns will be output regardless of length")
    else:
        messenger.log_only(f"Extracted {len(introns):,} introns (length filtering during scoring filter phase)")

    return introns


def load_introns_from_sequences(
    config: IntronICConfig,
    messenger: 'UnifiedMessenger'
) -> List[Intron]:
    """Load introns from pre-extracted sequences.

    Args:
        config: Pipeline configuration
        messenger: Unified messenger for output

    Returns:
        List of introns with sequences
    """
    messenger.info(f"Loading sequences: {config.input.sequence_file}")

    # Parse sequence file
    parser = SequenceParser()
    introns = list(parser.parse_file(config.input.sequence_file))
    messenger.log_only(f"Loaded {len(introns)} introns from sequence file")

    return introns


def score_introns(
    introns: List[Intron],
    config: IntronICConfig,
    messenger: 'UnifiedMessenger',
    reporter: IntronICProgressReporter
) -> List[Intron]:
    """Score introns with PWM matrices (parallel or sequential).

    PHASE 3: Restored parallel scoring from original intronIC v1.5.1
    Uses multiprocessing.Pool when config.performance.processes > 1

    Args:
        introns: List of introns to score
        config: Pipeline configuration
        messenger: Unified messenger for output
        reporter: Progress reporter (for progress bars only)

    Returns:
        List of introns with raw scores
    """
    messenger.info("Loading PWM matrices")

    # Determine PWM file paths
    pwm_file = Path(__file__).parent.parent / "data" / "scoring_matrices.fasta.iic"
    if not pwm_file.exists():
        raise FileNotFoundError(f"PWM file not found: {pwm_file}")

    u2_bp_file = Path(__file__).parent.parent / "data" / "u2.conserved_empirical_bp_pwm.iic"

    # Extract scoring configuration
    five_coords = (
        config.scoring.scoring_regions.five_start,
        config.scoring.scoring_regions.five_end
    )
    bp_coords = (
        config.scoring.scoring_regions.bp_start,
        config.scoring.scoring_regions.bp_end
    )
    three_coords = (
        config.scoring.scoring_regions.three_start,
        config.scoring.scoring_regions.three_end
    )
    ignore_nc_dnts = config.scoring.ignore_nc_dnts
    pseudocount = config.scoring.pseudocount
    n_workers = config.performance.processes

    # ========================================================================
    # PARALLEL SCORING (n_workers > 1)
    # ========================================================================
    if n_workers > 1:
        messenger.info(f"Calculating PWM scores (parallel, {n_workers} workers)")

        progress = reporter.create_progress()
        scored_introns = []
        failed_count = 0

        with progress:
            task = progress.add_task(
                "[cyan]Scoring introns...",
                total=len(introns)
            )

            with Pool(processes=n_workers) as pool:
                try:
                    # Use imap_unordered with chunking to avoid memory explosion
                    # This streams results back instead of collecting all at once
                    # Chunksize: balance between overhead and memory usage
                    chunksize = max(1, min(1000, len(introns) // (n_workers * 4)))

                    results_iter = pool.starmap(
                        _score_intron_worker,
                        zip(
                            introns,
                            repeat(pwm_file),
                            repeat(u2_bp_file),
                            repeat(five_coords),
                            repeat(bp_coords),
                            repeat(three_coords),
                            repeat(ignore_nc_dnts),
                            repeat(pseudocount)
                        ),
                        chunksize=chunksize
                    )

                    # Process results as they arrive
                    for scored_intron, error in results_iter:
                        if error is not None:
                            messenger.warning(f"Failed to score intron: {error}")
                            failed_count += 1
                        else:
                            scored_introns.append(scored_intron)

                        progress.update(task, advance=1)

                except KeyboardInterrupt:
                    messenger.warning("User interrupt - terminating workers")
                    pool.terminate()
                    raise
                finally:
                    pool.close()
                    pool.join()

    # ========================================================================
    # SEQUENTIAL SCORING (n_workers == 1)
    # ========================================================================
    else:
        messenger.info("Calculating PWM scores (sequential)")

        # Load PWMs once for sequential mode
        pwm_sets = PWMLoader.load_from_file(pwm_file, pseudocount=pseudocount)

        # Load U2 BP matrix if available
        if u2_bp_file.exists():
            from scoring.pwm import PWMSet
            u2_bp_matrices = PWMLoader.load_from_file(u2_bp_file, pseudocount=pseudocount)
            if 'bp' in u2_bp_matrices and u2_bp_matrices['bp'].u2_gtag:
                # PWMSet is frozen, so create new PWMSet with updated U2 GTAG matrix
                updated_matrices = dict(pwm_sets['bp'].matrices)
                updated_matrices[('u2', 'gtag')] = u2_bp_matrices['bp'].u2_gtag
                pwm_sets['bp'] = PWMSet(matrices=updated_matrices)
                messenger.log_only("Loaded conserved U2-type BP matrix")

        # Create scorer
        scorer = IntronScorer(
            pwm_sets=pwm_sets,
            five_coords=five_coords,
            bp_coords=bp_coords,
            three_coords=three_coords,
            ignore_nc_dnts=ignore_nc_dnts
        )

        # Score introns sequentially with error handling
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
                    messenger.warning(
                        f"Failed to score intron {intron.intron_id}: {str(e)}. Skipping."
                    )
                    failed_count += 1

                progress.update(task, advance=1)

    # ========================================================================
    # COMPLETION (both modes)
    # ========================================================================
    total_attempted = len(introns)
    messenger.log_only(
        f"Scored {format_count_with_percentage(len(scored_introns), total_attempted)} introns successfully"
    )
    if failed_count > 0:
        messenger.warning(
            f"Failed to score {format_count_with_percentage(failed_count, total_attempted)} introns "
            f"(see warnings above)"
        )

    return scored_introns


def normalize_scores(
    introns: List[Intron],
    config: IntronICConfig,
    messenger: 'UnifiedMessenger',
    reporter: IntronICProgressReporter
) -> Tuple[List[Intron], List[Intron], List[Intron], 'ScoreNormalizer']:
    """Normalize intron scores using z-score transformation.

    Args:
        introns: List of introns with raw scores
        config: Pipeline configuration
        messenger: Unified messenger for console and log output
        reporter: Progress reporter

    Returns:
        Tuple of (normalized experimental introns, u12 reference, u2 reference, normalizer)
    """
    messenger.info("Loading reference sequences")

    # Load reference data - use custom paths if provided, otherwise use defaults
    data_dir = Path(__file__).parent.parent / "data"

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

    u12_reference = load_reference_sequences(u12_file, messenger=messenger)
    u2_reference = load_reference_sequences(u2_file, messenger=messenger)

    messenger.log_only(f"Loaded {len(u12_reference)} U12-type and {len(u2_reference)} U2-type reference introns")

    # Score reference introns
    messenger.info("Scoring reference sequences")

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
            messenger.log_only("Loaded conserved U2-type BP matrix for reference scoring")

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
    messenger.info("Normalizing scores with z-score transformation")
    messenger.log_only("Fitting normalizer on reference data")

    normalizer = ScoreNormalizer()
    normalizer.fit(reference_introns, dataset_type='reference')

    # Transform reference introns (needed for classification)
    # Note: transform() takes an iterable and returns an iterator, so we materialize with list()
    u12_normalized = list(normalizer.transform(u12_scored, dataset_type='reference'))
    u2_normalized = list(normalizer.transform(u2_scored, dataset_type='reference'))

    # Transform experimental introns
    normalized_introns = list(normalizer.transform(introns, dataset_type='experimental'))

    messenger.log_only(f"Normalized {len(normalized_introns)} experimental introns")
    messenger.log_only(f"Normalized {len(u12_normalized)} U12-type and {len(u2_normalized)} U2-type reference introns")
    return normalized_introns, u12_normalized, u2_normalized, normalizer


def classify_with_pretrained_model(
    introns: List[Intron],
    model_path: Path,
    config: IntronICConfig,
    messenger: 'UnifiedMessenger',
    reporter: IntronICProgressReporter
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
        messenger: Unified messenger for console and log output
        reporter: Progress reporter

    Returns:
        Tuple of (classified introns, classification metrics)

    Note:
        The saved normalizer from training is NOT used. Instead, we fit a new normalizer
        on the experimental data to correct for species-specific score distributions.
    """
    messenger.info(f"Loading pretrained model from {model_path}")

    # Load model bundle
    if not model_path.exists():
        raise FileNotFoundError(f"Pretrained model not found: {model_path}")

    model_data = joblib.load(model_path)

    # Handle both old format (SVMEnsemble directly) and new format (dict bundle)
    if isinstance(model_data, dict):
        # New format: {'ensemble': ..., 'normalizer': ..., 'threshold': ...}
        ensemble = model_data['ensemble']
        saved_normalizer = model_data.get('normalizer', None)
        saved_threshold = model_data.get('threshold', config.scoring.threshold)
        messenger.log_only("Loaded model bundle (dict format)")
    else:
        # Old format: SVMEnsemble directly (backward compatibility)
        ensemble = model_data
        saved_normalizer = None
        saved_threshold = config.scoring.threshold
        messenger.log_only("Loaded model ensemble (legacy format - backward compatibility)")

    messenger.log_only(f"Loaded ensemble with {len(ensemble.models)} models")
    messenger.log_only(f"Using threshold: {config.scoring.threshold}")

    # Determine normalizer mode
    normalizer_mode = config.scoring.normalizer_mode
    if normalizer_mode == 'auto':
        # Use human scaler if available, otherwise adaptive
        if saved_normalizer is not None:
            normalizer_mode = 'human'
            messenger.log_only("Auto mode: Using saved human scaler (recommended)")
        else:
            normalizer_mode = 'adaptive'
            messenger.log_only("Auto mode: Falling back to adaptive (no saved scaler)")

    # Apply normalizer strategy
    if normalizer_mode == 'human':
        if saved_normalizer is None:
            raise ValueError(
                "Normalizer mode 'human' requested but model has no saved scaler. "
                "Retrain model or use '--normalizer_mode adaptive'."
            )
        messenger.info("Using human-trained normalizer (scaler from training species)")
        messenger.log_only("This preserves composition bias correction across species")
        normalizer = saved_normalizer
    else:  # adaptive
        messenger.info("Fitting normalizer on experimental data (domain adaptation)")
        messenger.log_only("Note: May cause FPs in U12-absent species (e.g., C. elegans)")
        messenger.log_only(f"Fitting normalizer on {len(introns)} experimental introns")
        normalizer = ScoreNormalizer()
        normalizer.fit(introns, dataset_type='unlabeled')

    normalized_introns = list(normalizer.transform(introns, dataset_type='experimental'))
    messenger.log_only(f"Normalized {len(normalized_introns)} experimental introns")

    # Classify using loaded ensemble
    messenger.info("Classifying with pretrained model")

    from classification.predictor import SVMPredictor
    predictor = SVMPredictor(
        threshold=config.scoring.threshold,  # Use config threshold (can override saved)
        n_jobs=config.performance.processes
    )

    classified_introns = list(predictor.predict(ensemble, normalized_introns))
    messenger.log_only(f"Classified {len(classified_introns)} introns")

    # Create metrics (limited since we skipped training)
    metrics = {
        'optimized_C': 'N/A (pretrained)',
        'cv_score': 'N/A (pretrained)',
        'n_models': len(ensemble.models),
        'pretrained': True,
        'model_path': str(model_path)
    }

    # Generate metadata for pretrained model usage
    run_metadata_path = config.output.get_output_path('.run_metadata.json')
    messenger.log_only(f"Recording pretrained model usage")
    run_metadata = generate_pretrained_metadata(
        model_path=model_path,
        threshold=config.scoring.threshold
    )
    write_metadata(run_metadata, run_metadata_path)
    messenger.log_only(f"Run metadata saved to {run_metadata_path}")

    return classified_introns, metrics


def _write_training_log(
    log_path: Path,
    classification_result,
    u12_reference: List[Intron],
    u2_reference: List[Intron],
    config: IntronICConfig
):
    """
    Write detailed training log with evaluation metrics and fold results.

    Args:
        log_path: Path to training log file
        classification_result: ClassificationResult from IntronClassifier
        u12_reference: U12 reference introns
        u2_reference: U2 reference introns
        config: Pipeline configuration
    """
    from utils.logging_utils import TrainingLogger

    with TrainingLogger(str(log_path)) as tlog:
        # Overview section
        tlog.section("TRAINING CONFIGURATION")
        tlog.metric("Species", config.output.species_name)
        tlog.metric("Classification Threshold", f"{config.scoring.threshold}%")
        tlog.metric("Random Seed", config.training.seed)
        tlog.metric("Max Iterations", config.training.max_iter)
        tlog.blank()

        tlog.subsection("Reference Data")
        tlog.metric("U12 reference introns", f"{len(u12_reference):,}")
        tlog.metric("U2 reference introns", f"{len(u2_reference):,}")
        tlog.metric("Total reference", f"{len(u12_reference) + len(u2_reference):,}")
        tlog.blank()

        tlog.subsection("Model Configuration")
        tlog.metric("Ensemble models", config.training.n_models)
        tlog.metric("Evaluation mode", config.training.eval_mode)

        if config.training.fixed_C:
            tlog.metric("C parameter", f"{config.training.fixed_C:.6e} (fixed)")
        else:
            tlog.metric("C parameter", "Optimized via grid search")
            tlog.metric("Optimization rounds", config.training.n_optimization_rounds)

        # Hyperparameter optimization section
        tlog.section("HYPERPARAMETER OPTIMIZATION")
        tlog.metric("Optimized C", f"{classification_result.parameters.C:.6e}")
        tlog.metric("CV score (balanced accuracy)", f"{classification_result.parameters.cv_score:.4f}")
        tlog.metric("Calibration method", classification_result.parameters.calibration_method)
        tlog.metric("Dual formulation", classification_result.parameters.dual)
        tlog.metric("Intercept scaling", classification_result.parameters.intercept_scaling)
        tlog.blank()

        # Evaluation results section
        eval_result = classification_result.eval_result
        if eval_result is not None:
            if hasattr(eval_result, 'mean_f1'):
                # Nested CV results
                tlog.section("NESTED CROSS-VALIDATION RESULTS")
                tlog.metric("Number of folds", eval_result.n_folds)
                tlog.blank()

                tlog.subsection("Aggregate Performance")
                tlog.metric("Mean F1 Score", f"{eval_result.mean_f1:.4f} ± {eval_result.std_f1:.4f}")
                tlog.metric("Mean PR-AUC", f"{eval_result.mean_pr_auc:.4f} ± {eval_result.std_pr_auc:.4f}")
                tlog.blank()

                tlog.subsection("Per-Fold Results")
                headers = ["Fold", "F1 Score", "PR-AUC", "Train U12", "Train U2", "Test U12", "Test U2"]
                rows = []
                for fold in eval_result.fold_results:
                    rows.append([
                        f"{fold.fold_idx + 1}/{eval_result.n_folds}",
                        f"{fold.f1_score:.4f}",
                        f"{fold.pr_auc:.4f}",
                        f"{fold.n_u12_train:,}",
                        f"{fold.n_u2_train:,}",
                        f"{fold.n_u12_test:,}",
                        f"{fold.n_u2_test:,}"
                    ])
                tlog.table(headers, rows)

            elif hasattr(eval_result, 'test_f1'):
                # Split evaluation results
                tlog.section("TRAIN/TEST SPLIT EVALUATION")
                tlog.blank()

                tlog.subsection("Data Split")
                tlog.metric("Training set", f"{eval_result.n_u2_train + eval_result.n_u12_train:,} introns "
                                           f"({eval_result.n_u2_train:,} U2, {eval_result.n_u12_train:,} U12)")
                if hasattr(eval_result, 'n_u2_val'):
                    tlog.metric("Validation set", f"{eval_result.n_u2_val + eval_result.n_u12_val:,} introns "
                                                  f"({eval_result.n_u2_val:,} U2, {eval_result.n_u12_val:,} U12)")
                tlog.metric("Test set", f"{eval_result.n_u2_test + eval_result.n_u12_test:,} introns "
                                        f"({eval_result.n_u2_test:,} U2, {eval_result.n_u12_test:,} U12)")
                tlog.blank()

                tlog.subsection("Test Set Performance (Honest Evaluation)")
                tlog.metric("F1 Score", f"{eval_result.test_f1:.4f}")
                tlog.metric("PR-AUC", f"{eval_result.test_pr_auc:.4f}")

        # Ensemble summary
        tlog.section("TRAINED ENSEMBLE")
        tlog.metric("Number of models", len(classification_result.ensemble.models))
        tlog.blank()

        tlog.subsection("Model Details")
        for i, model in enumerate(classification_result.ensemble.models, 1):
            tlog.info(f"Model {i}/{len(classification_result.ensemble.models)}:", indent=1)
            tlog.metric(f"  Training samples", f"{model.train_size:,}")
            tlog.metric(f"  U12 samples", f"{model.u12_count:,}")
            tlog.metric(f"  U2 samples", f"{model.u2_count:,}")
            tlog.metric(f"  C parameter", f"{model.parameters.C:.6e}")
            tlog.blank()


def classify_introns(
    introns: List[Intron],
    u12_reference: List[Intron],
    u2_reference: List[Intron],
    normalizer: 'ScoreNormalizer',
    config: IntronICConfig,
    messenger: 'UnifiedMessenger',
    reporter: IntronICProgressReporter
) -> Tuple[List[Intron], dict]:
    """Classify introns as U2 or U12 type.

    Args:
        introns: List of introns with z-scores
        u12_reference: U12 reference introns (scored and normalized)
        u2_reference: U2 reference introns (scored and normalized)
        normalizer: Fitted normalizer used for score transformation
        config: Pipeline configuration
        messenger: Unified messenger for console and log output
        reporter: Progress reporter

    Returns:
        Tuple of (classified introns, classification metrics)
    """
    messenger.info("Training SVM classifier")

    # Create training log if we're actually training (not using pretrained model or fixed C)
    training_log_path = None
    if config.training.pretrained_model_path is None and config.training.eval_mode != 'none':
        training_log_path = config.output.get_output_path('.training.log')
        messenger.log_only(f"Detailed training log will be written to: {training_log_path}")

    # Load optimizer configuration if specified
    param_grid_override = None
    optimizer_n_rounds = config.training.n_optimization_rounds
    optimizer_cv_folds = config.training.n_cv_folds
    optimizer_cv_processes = config.performance.cv_processes
    optimizer_max_iter = config.training.max_iter
    optimizer_n_points_initial = 13  # Default value
    eff_C_pos_range = (1e-3, 1e3)  # Default C bounds
    eff_C_neg_max = None  # Default C bounds

    # Initialize YAML training settings to None (may be overridden if YAML config provided)
    yaml_n_models = None
    yaml_subsample_u2 = None
    yaml_subsample_ratio = None
    yaml_training_max_iter = None
    yaml_training_random_state = None

    # Auto-load default config if it exists and user didn't specify --optimizer-config
    # Only auto-load when actually training (not using a pretrained model)
    if not config.training.optimizer_config_path and config.training.pretrained_model_path is None:
        from pathlib import Path
        default_config = Path(__file__).parent.parent / "config" / "training_default.yaml"
        if default_config.exists():
            # Update config to use default
            config = replace(
                config,
                training=replace(config.training, optimizer_config_path=default_config)
            )
            messenger.log_only(f"Auto-loading default configuration: {default_config}")

    if config.training.optimizer_config_path:
        messenger.log_only(f"Loading optimizer configuration from: {config.training.optimizer_config_path}")
        try:
            from classification.config_loader import load_optimizer_config
            optimizer_from_yaml = load_optimizer_config(config.training.optimizer_config_path)

            # Extract ALL settings from YAML optimizer (not just param_grid!)
            param_grid_override = optimizer_from_yaml.param_grid_override
            optimizer_n_rounds = optimizer_from_yaml.n_rounds
            optimizer_cv_folds = optimizer_from_yaml.cv_folds
            optimizer_cv_processes = optimizer_from_yaml.n_jobs  # n_jobs → cv_processes
            optimizer_max_iter = optimizer_from_yaml.max_iter
            optimizer_n_points_initial = optimizer_from_yaml.n_points_initial

            # Extract C bounds if specified in config
            eff_C_pos_range = getattr(optimizer_from_yaml, 'eff_C_pos_range', (1e-3, 1e3))
            eff_C_neg_max = getattr(optimizer_from_yaml, 'eff_C_neg_max', None)

            # Extract training settings if present in YAML (override CLI defaults)
            yaml_n_models = getattr(optimizer_from_yaml, 'training_n_models', None)
            yaml_subsample_u2 = getattr(optimizer_from_yaml, 'training_subsample_u2', None)
            yaml_subsample_ratio = getattr(optimizer_from_yaml, 'training_subsample_ratio', None)
            yaml_training_max_iter = getattr(optimizer_from_yaml, 'training_max_iter', None)
            yaml_training_random_state = getattr(optimizer_from_yaml, 'training_random_state', None)

            if yaml_n_models is not None:
                config = replace(
                    config,
                    training=replace(config.training, n_models=yaml_n_models)
                )
                messenger.log_only(f"  Ensemble models: {yaml_n_models} (from YAML)")
            if yaml_subsample_u2 is not None:
                messenger.log_only(f"  Subsample U2: {yaml_subsample_u2} (from YAML)")
            if yaml_subsample_ratio is not None:
                messenger.log_only(f"  Subsample ratio: {yaml_subsample_ratio} (from YAML)")
            if yaml_training_max_iter is not None:
                messenger.log_only(f"  Training max_iter: {yaml_training_max_iter} (from YAML)")
            if yaml_training_random_state is not None:
                messenger.log_only(f"  Training random_state: {yaml_training_random_state} (from YAML)")

            messenger.log_only(f"Loaded custom optimizer configuration:")
            messenger.log_only(f"  Optimization rounds: {optimizer_n_rounds}")
            messenger.log_only(f"  Initial grid points: {optimizer_n_points_initial}")
            messenger.log_only(f"  CV folds: {optimizer_cv_folds}")
            messenger.log_only(f"  Parallel jobs: {optimizer_cv_processes}")
            messenger.log_only(f"  Max iterations: {optimizer_max_iter}")
            messenger.log_only(f"  Parameter grid: {len(param_grid_override)} hyperparameter sets")
            if hasattr(optimizer_from_yaml, 'eff_C_pos_range'):
                messenger.log_only(f"  C bounds: eff_C_pos_range={eff_C_pos_range}, eff_C_neg_max={eff_C_neg_max}")
        except Exception as e:
            messenger.warning(f"Failed to load optimizer config: {e}")
            messenger.warning("Continuing with default optimizer settings...")

    # Create classifier with correct parameter names
    # IntronClassifier API uses:
    # - classification_threshold (not threshold)
    # - n_ensemble_models (not n_models)
    # - fixed_c (not fixed_C)
    # - random_state (not seed)
    # - cv_processes (for cross-validation parallelization)
    # - classification_processes (for prediction parallelization)
    # - param_grid_override (optional custom parameter grid)
    # - n_points_initial (initial grid points for optimization)
    # - subsample_u2, subsample_ratio (ensemble diversity settings)
    # - max_iter (for LinearSVC convergence during ensemble training)

    # Use YAML training settings if available, otherwise use CLI defaults
    ensemble_subsample_u2 = yaml_subsample_u2 if yaml_subsample_u2 is not None else True
    ensemble_subsample_ratio = yaml_subsample_ratio if yaml_subsample_ratio is not None else 0.8
    ensemble_max_iter = yaml_training_max_iter if yaml_training_max_iter is not None else optimizer_max_iter
    ensemble_random_state = yaml_training_random_state if yaml_training_random_state is not None else config.training.seed

    # Get use_fold_averaged_params: CLI arg takes precedence over YAML
    # Priority: CLI arg (if explicitly set) > YAML value > default (False)
    yaml_use_fold_averaged = getattr(optimizer_from_yaml, 'training_use_fold_averaged_params', None) if optimizer_from_yaml else None
    if config.training.use_fold_averaged_params is not None:
        # CLI arg was explicitly provided
        use_fold_averaged_params = config.training.use_fold_averaged_params
    elif yaml_use_fold_averaged is not None:
        # Use YAML value as fallback
        use_fold_averaged_params = yaml_use_fold_averaged
    else:
        # Default
        use_fold_averaged_params = False

    classifier = IntronClassifier(
        n_optimization_rounds=optimizer_n_rounds,
        classification_threshold=config.scoring.threshold,
        n_ensemble_models=config.training.n_models,
        subsample_u2=ensemble_subsample_u2,
        subsample_ratio=ensemble_subsample_ratio,
        fixed_c=config.training.fixed_C,
        optimize_c=(config.training.fixed_C is None),
        random_state=ensemble_random_state,
        cv_processes=optimizer_cv_processes,
        classification_processes=config.performance.processes,
        max_iter=ensemble_max_iter,
        eval_mode=config.training.eval_mode,
        n_cv_folds=optimizer_cv_folds,
        test_fraction=config.training.test_fraction,
        param_grid_override=param_grid_override,
        n_points_initial=optimizer_n_points_initial,
        eff_C_pos_range=eff_C_pos_range,
        eff_C_neg_max=eff_C_neg_max,
        use_fold_averaged_params=use_fold_averaged_params
    )

    # Run complete classification pipeline (optimize + train + classify)
    messenger.log_only(f"Running classification on {len(introns)} experimental introns")
    messenger.log_only(f"Reference data: {len(u12_reference)} U12-type, {len(u2_reference)} U2-type")

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

    messenger.log_only(f"Classification complete")
    messenger.log_only(f"  Optimized C: {metrics['optimized_C']:.6e}")
    messenger.log_only(f"  Models trained: {metrics['n_models']}")

    # Log evaluation metrics if available
    if 'mean_f1' in metrics:
        messenger.log_only(f"  Nested CV results ({metrics['n_cv_folds']} folds):")
        messenger.log_only(f"    Mean F1: {metrics['mean_f1']:.4f} ± {metrics['std_f1']:.4f}")
        messenger.log_only(f"    Mean PR-AUC: {metrics['mean_pr_auc']:.4f} ± {metrics['std_pr_auc']:.4f}")
    elif 'f1' in metrics:
        messenger.log_only(f"  Test set evaluation:")
        messenger.log_only(f"    F1: {metrics['f1']:.4f}")
        messenger.log_only(f"    PR-AUC: {metrics['pr_auc']:.4f}")

    # Write detailed training log if evaluation was performed
    if training_log_path and result.eval_result is not None:
        try:
            from utils.logging_utils import TrainingLogger
            _write_training_log(
                training_log_path,
                result,
                u12_reference,
                u2_reference,
                config
            )
            messenger.log_only(f"Detailed training log written to: {training_log_path}")
        except Exception as e:
            messenger.warning(f"Failed to write training log: {e}")

    # Generate training reference plots if evaluation was performed
    if result.eval_result is not None:
        messenger.log_only("Generating training reference plots")
        try:
            from visualization.plots import plot_training_results

            # Extract normalized z-scores from reference introns
            u2_scores = np.array([
                [i.scores.five_z_score, i.scores.bp_z_score]
                for i in u2_reference
                if i.scores and i.scores.five_z_score is not None and i.scores.bp_z_score is not None
            ])
            u12_scores = np.array([
                [i.scores.five_z_score, i.scores.bp_z_score]
                for i in u12_reference
                if i.scores and i.scores.five_z_score is not None and i.scores.bp_z_score is not None
            ])

            # Get PR curves and AUC based on evaluation type
            if hasattr(result.eval_result, 'mean_f1'):
                # Nested CV - multiple curves
                pr_curves = result.eval_result.pr_curves
                pr_auc = result.eval_result.mean_pr_auc
            else:
                # Split eval - single curve
                pr_curves = [(result.eval_result.precision, result.eval_result.recall)]
                pr_auc = result.eval_result.test_pr_auc

            plot_training_results(
                u2_scores=u2_scores,
                u12_scores=u12_scores,
                pr_curves=pr_curves,
                pr_auc=pr_auc,
                output_dir=config.output.output_dir,
                species_name=config.output.base_filename,
                fig_dpi=300
            )
            messenger.log_only("Successfully generated training reference plots")
        except Exception as plot_error:
            messenger.warning(f"Failed to generate training plots: {plot_error}")
            # Continue even if plotting fails

    # Save trained model with human-trained normalizer for cross-species classification
    model_path = config.output.get_output_path('.model.pkl')
    messenger.log_only(f"Saving trained model to {model_path}")
    model_bundle = {
        'ensemble': result.ensemble,
        'normalizer': normalizer,  # Save human-trained scaler for cross-species use
        'threshold': config.scoring.threshold
    }
    joblib.dump(model_bundle, model_path, compress=3)
    messenger.log_only(f"Model saved successfully with human-trained normalizer")
    messenger.log_only(f"Use '--normalizer_mode human' for cross-species classification")

    # Generate and save training metadata
    metadata_path = model_path.with_suffix('.metadata.json')
    messenger.log_only(f"Generating training metadata")
    metadata = generate_training_metadata(
        model_name=model_path.stem,
        u12_reference_path=config.scoring.reference_u12s,
        u2_reference_path=config.scoring.reference_u2s,
        u12_introns=u12_reference,
        u2_introns=u2_reference,
        optimized_C=result.parameters.C,
        calibration_method=result.parameters.calibration_method,
        cv_score=result.parameters.cv_score,
        n_models=len(result.ensemble.models),
        threshold=config.scoring.threshold,
        eval_result=result.eval_result,
        max_iter=config.training.max_iter,
        kernel='linear',
        seed=config.training.seed
    )
    write_metadata(metadata, metadata_path)
    messenger.log_only(f"Training metadata saved to {metadata_path}")

    return list(result.classified_introns), metrics


def write_outputs(
    introns: List[Intron],
    config: IntronICConfig,
    messenger: 'UnifiedMessenger',
    reporter: IntronICProgressReporter,
    scored_only: Optional[List[Intron]] = None,
    skip_sequences: bool = False
):
    """Write output files.

    Port from: intronIC.py:4820-4912 (filter_introns_write_files) and 5232-5267 (main)

    Args:
        introns: All introns for .bed.iic, .meta.iic, .introns.iic (scored + omitted)
        config: Pipeline configuration
        messenger: Unified messenger for console and log output
        reporter: Progress reporter
        scored_only: Introns for .score_info.iic (scored only, no omitted). If None, uses introns.
        skip_sequences: If True, skip writing .introns.iic (already written earlier)

    Notes:
        Original intronIC writes different intron sets to different files:
        - .introns.iic/.seqs.iic: ALL introns except duplicates (unless -d)
        - .bed.iic, .meta.iic: Scored + omitted non-duplicates
        - .score_info.iic: ONLY scored introns (no omitted)

        This is achieved by:
        1. Two-phase writing for .bed/.meta (omitted in filter function, scored in main)
        2. Single write for .seqs (in filter function)
        3. Single write for .score_info (in main, only finalized_introns)
    """
    messenger.info("Writing output files")

    # Filter duplicates if not including them
    # Port from: intronIC.py:4806-4807
    # Note: For normal mode, duplicates are already excluded by merge_scored_and_omitted_introns()
    # This is only needed for sequences_only mode where we don't call that merge function
    if not config.extraction.include_duplicates:
        original_count = len(introns)
        introns = [i for i in introns if not (i.metadata and i.metadata.duplicate)]
        filtered_count = original_count - len(introns)
        if filtered_count > 0:
            messenger.log_only(f"Filtered out {filtered_count} duplicate introns (use -d to include)")

    output_dir = config.output.output_dir
    base_name = config.output.base_filename
    species_name = config.output.species_name
    simple_name = config.output.uninformative_naming
    no_abbreviate = config.output.no_abbreviate

    # Write BED file (scored + omitted non-duplicates)
    # Port from: intronIC.py:4823-4835 (omitted) + 5232-5237 (scored)
    bed_path = output_dir / f"{base_name}.bed.iic"
    messenger.log_only(f"Writing BED file: {bed_path}")
    bed_writer = BEDWriter(bed_path)
    with bed_writer:
        for intron in introns:
            bed_writer.write_intron(intron, species_name=species_name, simple_name=simple_name, no_abbreviate=no_abbreviate)
    messenger.log_only(f"Wrote {len(introns)} introns to BED file")

    # Write metadata file (scored + omitted non-duplicates)
    # Port from: intronIC.py:4823-4835 (omitted) + 5240, 5262-5264 (scored)
    meta_path = output_dir / f"{base_name}.meta.iic"
    messenger.log_only(f"Writing metadata file: {meta_path}")
    meta_writer = MetaWriter(meta_path)
    with meta_writer:
        for intron in introns:
            meta_writer.write_intron(intron, species_name=species_name, simple_name=simple_name, no_abbreviate=no_abbreviate)
    messenger.log_only(f"Wrote metadata for {len(introns)} introns")

    # Write sequences file (all introns except duplicates, unless -d)
    # Port from: intronIC.py:4842-4845
    seq_path = output_dir / f"{base_name}.introns.iic"
    if not skip_sequences:
        messenger.log_only(f"Writing sequences file: {seq_path}")
        seq_writer = SequenceWriter(seq_path)
        with seq_writer:
            for intron in introns:
                seq_writer.write_intron(intron, species_name=species_name, simple_name=simple_name, no_abbreviate=no_abbreviate)
        messenger.log_only(f"Wrote sequences for {len(introns)} introns")
    else:
        messenger.log_only(f"Sequences already written during scoring phase: {seq_path}")

    # Write score info file (ONLY scored introns, no omitted)
    # Port from: intronIC.py:5239-5261
    # CRITICAL: Only write scored introns, not omitted ones
    score_introns = scored_only if scored_only is not None else introns
    score_path = output_dir / f"{base_name}.score_info.iic"
    messenger.log_only(f"Writing score info file: {score_path}")
    score_writer = ScoreWriter(score_path)
    with score_writer:
        for intron in score_introns:
            score_writer.write_intron(intron, species_name=species_name, simple_name=simple_name, no_abbreviate=no_abbreviate)
    messenger.log_only(f"Wrote score info for {len(score_introns)} introns")

    output_files = {
        "Metadata": meta_path,
        "BED": bed_path,
        "Sequences": seq_path,
        "Scores": score_path,
        "Log": config.output.get_output_path('.log'),
    }

    reporter.print_file_tree(output_files)
    messenger.success("All output files written successfully")


def main_train(config: IntronICConfig):
    """Train a model on reference data only (no genome/annotation needed).

    This is a pure training workflow:
    1. Load reference U12/U2 sequences
    2. Score reference sequences with PWM matrices
    3. Normalize scores
    4. Train ensemble of models with cross-validation
    5. Save model and metadata to disk

    Args:
        config: Training configuration

    No genome/annotation required!
    """
    # Track start time
    start_time = time.time()

    # Setup logging
    logger, log_console = setup_logging(config)
    reporter = IntronICProgressReporter(quiet=config.output.quiet)

    from .messenger import UnifiedMessenger
    messenger = UnifiedMessenger(
        console=reporter.console,
        log_console=log_console,
        logger=logger,
        quiet=config.output.quiet
    )

    # Print header
    reporter.console.print(
        "\n[bold cyan]═══════════════════════════════════════════════════════[/bold cyan]"
    )
    reporter.console.print(
        "[bold cyan]                  intronIC TRAINING MODE                  [/bold cyan]"
    )
    reporter.console.print(
        "[bold cyan]═══════════════════════════════════════════════════════[/bold cyan]\n"
    )
    messenger.info(f"Training model: {config.output.species_name}")
    messenger.info("No genome/annotation needed - using reference sequences only")

    pipeline_steps = [
        "Load reference data",
        "Score reference sequences",
        "Normalize scores",
        "Train classifier"
    ]
    messenger.console_only("")
    reporter.print_pipeline_steps(pipeline_steps)

    try:
        # Steps 1-3: Load, score, and normalize reference data
        messenger.step(1, "Load Reference Data", pipeline_steps)
        messenger.step(2, "Score Reference Sequences", pipeline_steps)
        messenger.step(3, "Normalize Scores", pipeline_steps)

        # Load and process reference sequences
        # This duplicates some logic from normalize_scores() but avoids issues with empty experimental introns
        data_dir = Path(__file__).parent.parent / "data"

        u12_file = config.scoring.reference_u12s or (data_dir / "u12_reference.introns.iic.gz")
        u2_file = config.scoring.reference_u2s or (data_dir / "u2_reference.introns.iic.gz")

        if not u12_file.exists() or not u2_file.exists():
            raise FileNotFoundError(f"Reference data not found. U12: {u12_file}, U2: {u2_file}")

        messenger.log_only("Loading reference sequences")
        u12_reference = load_reference_sequences(u12_file, messenger=messenger)
        u2_reference = load_reference_sequences(u2_file, messenger=messenger)
        messenger.log_only(f"Loaded {len(u12_reference)} U12 and {len(u2_reference)} U2 reference introns")

        # Score reference introns
        messenger.log_only("Scoring reference sequences")
        all_reference = u12_reference + u2_reference
        scored_reference = score_introns(all_reference, config, messenger, reporter)

        # Split back into U12 and U2
        u12_scored = scored_reference[:len(u12_reference)]
        u2_scored = scored_reference[len(u12_reference):]

        # Normalize scores
        messenger.log_only("Normalizing scores with z-score transformation")
        from scoring.normalizer import ScoreNormalizer
        normalizer = ScoreNormalizer()
        normalizer.fit(scored_reference, dataset_type='reference')

        u12_ref_norm = list(normalizer.transform(u12_scored, dataset_type='reference'))
        u2_ref_norm = list(normalizer.transform(u2_scored, dataset_type='reference'))

        messenger.success(
            f"Loaded, scored, and normalized {len(u12_ref_norm)} U12 and {len(u2_ref_norm)} U2 reference introns"
        )

        # Step 4: Train classifier (model is saved internally by classify_introns)
        messenger.step(4, "Train Classifier", pipeline_steps)
        # Pass empty list for experimental introns - we're only training on references
        # classify_introns() will train the model and save it to disk
        classified_introns, metrics = classify_introns(
            introns=[],  # No experimental introns in train mode
            u12_reference=u12_ref_norm,
            u2_reference=u2_ref_norm,
            normalizer=normalizer,
            config=config,
            messenger=messenger,
            reporter=reporter
        )
        messenger.success("Model trained and saved")

        # Print final summary
        elapsed = time.time() - start_time
        hours, remainder = divmod(int(elapsed), 3600)
        minutes, seconds = divmod(remainder, 60)

        messenger.console_only("")
        if hours > 0:
            time_str = f"{hours}h {minutes}m {seconds}s"
        elif minutes > 0:
            time_str = f"{minutes}m {seconds}s"
        else:
            time_str = f"{seconds}s"

        messenger.success(f"Training complete! (Runtime: {time_str})")

    except Exception as e:
        messenger.error(f"Training failed: {str(e)}")
        raise


def main_classify(config: IntronICConfig):
    """Run the complete intronIC classification pipeline.

    This is the standard workflow:
    1. Load input data (genome/annotation/bed/sequences)
    2. Extract introns
    3. Filter introns
    4. Score introns
    5. Load or train model
    6. Classify introns
    7. Write outputs

    Args:
        config: Pipeline configuration
    """
    # Track start time for runtime reporting
    start_time = time.time()

    # Setup logging and reporting with ANSI color support in log files
    logger, log_console = setup_logging(config)
    reporter = IntronICProgressReporter(quiet=config.output.quiet)

    # Create unified messenger for synchronized console + log output
    # Both destinations get Rich formatting with ANSI colors
    from .messenger import UnifiedMessenger
    messenger = UnifiedMessenger(
        console=reporter.console,
        log_console=log_console,
        logger=logger,
        quiet=config.output.quiet
    )

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

    # Show initial pipeline overview (console only, not in log)
    messenger.console_only("")  # Blank line
    reporter.print_pipeline_steps(pipeline_steps)

    try:
        # Step 1: Load input data
        messenger.step(1, "Load Input Data", pipeline_steps)

        if config.input.mode == 'annotation':
            genome_reader = load_genome(config, messenger)
            introns = extract_introns_from_annotation(
                config, genome_reader, messenger, reporter
            )
        elif config.input.mode == 'bed':
            genome_reader = load_genome(config, messenger)
            introns = extract_introns_from_bed(
                config, genome_reader, messenger, reporter
            )
        elif config.input.mode == 'sequences':
            introns = load_introns_from_sequences(config, messenger)
        else:
            raise ValueError(f"Unknown input mode: {config.input.mode}")

        messenger.success(f"Loaded {len(introns):,} introns")

        # If sequences only, skip to output
        if config.scoring.sequences_only:
            messenger.warning("Sequences-only mode: Skipping classification")
            # In sequences-only mode, no introns are scored, so scored_only=None
            # This means .score_info.iic will be empty or contain all introns (but they have no scores)
            # Original behavior: writes all to .bed/.meta/.seqs, exits before .score_info is created
            write_outputs(introns, config, messenger, reporter, scored_only=[])
            messenger.success("Pipeline complete!")
            return

        # Step 2: Extract introns (already done above)
        messenger.step(2, "Extract Introns", pipeline_steps)
        messenger.success(f"Extracted {len(introns):,} introns")

        # Filter introns before scoring (duplicates, short introns, longest isoform)
        # This matches original intronIC behavior where filtering happens BEFORE scoring
        # to avoid scoring 5x more introns than necessary (which causes O(n²) slowdown)
        messenger.info("Filtering introns before scoring")

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
        total = stats.total_introns

        messenger.log_only(
            f"Filtering results: {format_count_with_percentage(stats.kept_introns, total)} "
            f"introns kept for scoring"
        )
        messenger.log_only(
            f"Omitted: {format_count_with_percentage(stats.omitted_short, total)} short, "
            f"{format_count_with_percentage(stats.omitted_ambiguous, total)} ambiguous, "
            f"{format_count_with_percentage(stats.omitted_noncanonical, total)} non-canonical, "
            f"{format_count_with_percentage(stats.omitted_isoform, total)} non-longest isoform, "
            f"{format_count_with_percentage(stats.omitted_overlap, total)} overlapping"
        )
        messenger.log_only(f"Duplicates marked: {format_count_with_percentage(stats.duplicates, total)}")

        messenger.success(
            f"Filtered to {len(filtered_introns):,} introns for scoring "
            f"(removed {len(introns) - len(filtered_introns):,})"
        )

        # Important: Use filtered_introns for scoring, but keep original introns list
        # for potential output (user may want duplicates via -d flag)
        introns_for_scoring = filtered_introns

        # Step 3: Score introns
        messenger.step(3, "Score Introns", pipeline_steps)
        scored_introns = score_introns(introns_for_scoring, config, messenger, reporter)
        messenger.success(f"Scored {len(scored_introns):,} introns")

        # PHASE 2 OPTIMIZATION: Write sequences to file now (while in memory)
        # Then clear large sequences before classification to save ~10-12 GB
        messenger.info("Writing intron sequences to file")
        seq_output_path = config.output.output_dir / f"{config.output.base_filename}.introns.iic"
        seq_writer = SequenceWriter(seq_output_path)

        # Write ALL introns (scored + omitted), matching write_outputs() behavior
        # Duplicates will be filtered later in write_outputs() if not -d flag
        introns_written = 0
        with seq_writer:
            for intron in introns:  # ALL introns, not just scored_introns
                # Filter duplicates if not -d flag (matches write_outputs line 1536-1544)
                if not config.extraction.include_duplicates:
                    if intron.metadata and intron.metadata.duplicate:
                        continue
                seq_writer.write_intron(
                    intron,
                    species_name=config.output.species_name,
                    simple_name=config.output.uninformative_naming,
                    no_abbreviate=config.output.no_abbreviate,
                    include_score=False  # Scores not final yet (not classified)
                )
                introns_written += 1
        messenger.log_only(f"Wrote sequences for {introns_written} introns to {seq_output_path.name}")

        # Clear large sequences to reduce memory before classification
        messenger.log_only("Clearing large sequences to reduce memory (keeping scored regions)")
        scored_introns = clear_large_sequences_for_classification(scored_introns)
        gc.collect()
        messenger.log_only("Memory freed from full sequences (~10-12 GB for large genomes)")

        # Check if using pretrained model
        if config.training.pretrained_model_path:
            # Skip normalization/training - use pretrained model
            messenger.step(4, "Classify with Pretrained Model", pipeline_steps)
            classified_introns, metrics = classify_with_pretrained_model(
                scored_introns, config.training.pretrained_model_path, config, messenger, reporter
            )
        else:
            # Normal flow: normalize + train + classify
            # Step 4: Normalize scores
            messenger.step(4, "Normalize Scores", pipeline_steps)
            normalized_introns, u12_reference, u2_reference, normalizer = normalize_scores(
                scored_introns, config, messenger, reporter
            )
            messenger.success("Scores normalized")

            # Step 5: Classify
            messenger.step(5, "Classify Introns", pipeline_steps)
            classified_introns, metrics = classify_introns(
                normalized_introns, u12_reference, u2_reference, normalizer, config, messenger, reporter
            )

        # Count classifications based on threshold (for reporting "high confidence" U12s)
        # Note: type_id is based on raw classifier (>50%), but reporting uses threshold
        u12_count = sum(
            1 for i in classified_introns
            if i.scores and i.scores.svm_score >= config.scoring.threshold
        )
        u2_count = len(classified_introns) - u12_count

        # Count AT-AC introns among high-confidence U12s (score >= threshold)
        atac_count = sum(
            1 for i in classified_introns
            if (i.scores and i.scores.svm_score >= config.scoring.threshold and
                i.sequences and i.sequences.terminal_dinucleotides == 'AT-AC')
        )

        reporter.print_classification_summary(
            total=len(classified_introns),
            u12_count=u12_count,
            u2_count=u2_count,
            threshold=config.scoring.threshold
        )

        # Log classification summary to log file
        total_classified = len(classified_introns)
        messenger.log_only(f"Classification results:")
        messenger.log_only(f"  {format_count_with_percentage(atac_count, total_classified)} putative AT-AC U12-type introns found")
        messenger.log_only(f"  {format_count_with_percentage(u12_count, total_classified)} putative U12-type introns found with scores > {config.scoring.threshold}%")
        messenger.log_only(f"  {format_count_with_percentage(u2_count, total_classified)} introns classified as U2-type")

        # Collect and log splice site boundary statistics (separate by U12/U2)
        from collections import Counter
        boundaries_u12 = Counter()
        boundaries_u2 = Counter()

        for intron in classified_introns:
            # Count all boundaries (canonical and non-canonical) by type
            # Use threshold to match the main counts, not raw classifier type_id
            if intron.metadata and intron.sequences and intron.sequences.terminal_dinucleotides:
                dnts = intron.sequences.terminal_dinucleotides
                if intron.scores and intron.scores.svm_score >= config.scoring.threshold:
                    boundaries_u12[dnts] += 1
                else:
                    boundaries_u2[dnts] += 1

        # Log U12 boundary statistics
        if boundaries_u12:
            total_u12 = sum(boundaries_u12.values())
            # Sort by count (descending), then alphabetically by dinucleotide
            sorted_boundaries = sorted(boundaries_u12.items(), key=lambda x: (-x[1], x[0]))[:20]
            lines = [
                f"  {i:2d}. {dnts:8s} {count:6,} ({(count / total_u12) * 100:5.2f}%)"
                for i, (dnts, count) in enumerate(sorted_boundaries, 1)
            ]
            log_data_block(
                logger,
                "Top 20 splice site boundaries (U12-type introns)",
                lines
            )

        # Log U2 boundary statistics
        if boundaries_u2:
            total_u2 = sum(boundaries_u2.values())
            # Sort by count (descending), then alphabetically by dinucleotide
            sorted_boundaries = sorted(boundaries_u2.items(), key=lambda x: (-x[1], x[0]))[:20]
            lines = [
                f"  {i:2d}. {dnts:8s} {count:6,} ({(count / total_u2) * 100:5.2f}%)"
                for i, (dnts, count) in enumerate(sorted_boundaries, 1)
            ]
            log_data_block(
                logger,
                "Top 20 splice site boundaries (U2-type introns)",
                lines
            )

        # Save classification metrics to JSON file
        if metrics:
            metrics_path = config.output.get_output_path('.metrics.json')
            messenger.log_only(f"Saving classification metrics to {metrics_path}")
            with open(metrics_path, 'w') as f:
                json.dump(metrics, f, indent=2)

        # Generate visualization plots
        messenger.log_only("Generating visualization plots")
        try:
            plot_classification_results(
                introns=classified_introns,
                output_dir=config.output.output_dir,
                species_name=config.output.base_filename,
                threshold=config.scoring.threshold,
                fig_dpi=300
            )
            messenger.log_only("Successfully generated classification plots")
        except Exception as plot_error:
            messenger.warning(f"Failed to generate plots: {plot_error}")
            # Continue even if plotting fails

        # Step 6: Write outputs
        messenger.step(6, "Write Outputs", pipeline_steps)

        # Merge classified introns with omitted introns for complete meta output
        # This matches original intronIC behavior where .meta.iic includes all introns
        # (scored + omitted), not just the ones that went through classification
        all_introns_for_output = merge_scored_and_omitted_introns(
            classified_introns, introns, messenger
        )

        # Write outputs with different intron sets for different files:
        # - all_introns_for_output: for .bed.iic, .meta.iic, .introns.iic (scored + omitted)
        # - classified_introns: for .score_info.iic (scored only, no omitted)
        # Port from: intronIC.py writes finalized_introns to .score_info (line 5239-5261)
        write_outputs(
            all_introns_for_output,
            config,
            messenger,
            reporter,
            scored_only=classified_introns,
            skip_sequences=True  # PHASE 2: Sequences already written after scoring
        )

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

        messenger.success(f"Pipeline complete! (Runtime: {runtime_str})")

    except Exception as e:
        logger.exception("Pipeline failed with error")
        messenger.error(f"Pipeline failed: {str(e)}")
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

        # Handle --generate-config
        if getattr(parsed_args, 'generate_config', False):
            from .config_loader import generate_config_file
            generate_config_file()
            return 0

        # Load and merge configuration file (if exists)
        from .config_loader import ConfigLoader

        config_path = ConfigLoader.find_config(
            custom_path=getattr(parsed_args, 'config', None)
        )

        if config_path:
            try:
                config_data = ConfigLoader.load_config(config_path)
                ConfigLoader.merge_with_args(config_data, parsed_args)
                print(f"Loaded configuration from: {config_path}")
            except ImportError as e:
                print(f"Warning: Could not load config file: {e}", file=sys.stderr)
                print("Install tomli for config file support: pip install tomli", file=sys.stderr)
            except Exception as e:
                print(f"Warning: Error loading config file {config_path}: {e}", file=sys.stderr)
                print("Continuing with CLI arguments only...", file=sys.stderr)

        # Create configuration
        config = IntronICConfig.from_args(parsed_args)

        # Route to appropriate entry point based on command
        command = getattr(parsed_args, 'command', 'classify')

        if command == 'train':
            # Train mode: No genome/annotation needed
            main_train(config)
        elif command == 'classify':
            # Classify mode: Standard pipeline
            main_classify(config)
        else:
            raise ValueError(f"Unknown command: {command}")

        return 0

    except KeyboardInterrupt:
        print("\n\nPipeline interrupted by user", file=sys.stderr)
        return 130  # Standard exit code for SIGINT

    except Exception as e:
        print(f"\n\nError: {str(e)}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
