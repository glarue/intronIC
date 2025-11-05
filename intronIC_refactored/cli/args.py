"""
Argument parser for intronIC CLI.

Handles command-line argument parsing and validation.
"""

import argparse
from pathlib import Path
from typing import Optional


class IntronICArgumentParser:
    """Argument parser for intronIC command-line interface."""

    def __init__(self):
        self.parser = self._create_parser()

    def parse_args(self, args: Optional[list] = None):
        """Parse command-line arguments.

        Args:
            args: Optional list of arguments (defaults to sys.argv)

        Returns:
            Namespace object with parsed arguments
        """
        parsed = self.parser.parse_args(args)
        self._validate_args(parsed)
        return parsed

    def _create_parser(self) -> argparse.ArgumentParser:
        """Create the argument parser with all options."""
        parser = argparse.ArgumentParser(
            description="intronIC: Intron classification and extraction tool",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="""
Examples:
  # Standard classification
  intronIC -g genome.fa -a annotation.gff3 -n species_name

  # Extract sequences only (no classification)
  intronIC -g genome.fa -a annotation.gff3 -n species_name -s

  # Parallel processing with 8 cores
  intronIC -g genome.fa -a annotation.gff3 -n species_name -p 8

  # Stricter threshold (95%)
  intronIC -g genome.fa -a annotation.gff3 -n species_name -t 95
"""
        )

        # Version
        parser.add_argument(
            '--version',
            action='version',
            version='intronIC 1.5.1'
        )

        # Required arguments
        required = parser.add_argument_group('required arguments')
        required.add_argument(
            '-n', '--species_name',
            required=True,
            help='Species name (e.g., homo_sapiens)'
        )

        # Input selection
        input_group = parser.add_argument_group(
            'input selection',
            'Choose one mode: (1) -g + -a for annotation-based extraction, '
            '(2) -g + -b for BED-based extraction, or (3) -q for pre-extracted sequences'
        )
        input_group.add_argument(
            '-g', '--genome',
            type=Path,
            help='Path to genome FASTA file (required with -a or -b)'
        )
        input_group.add_argument(
            '-a', '--annotation',
            type=Path,
            help='Path to GFF3/GTF annotation file (requires -g)'
        )
        input_group.add_argument(
            '-b', '--bed',
            type=Path,
            help='Path to BED file with intron coordinates (requires -g)'
        )
        input_group.add_argument(
            '-q', '--sequence_file',
            type=Path,
            help='Path to pre-extracted intron sequences (.iic format, standalone)'
        )

        # Output options
        output_group = parser.add_argument_group('output options')
        output_group.add_argument(
            '-o', '--output_dir',
            type=Path,
            default=Path.cwd(),
            help='Output directory (default: current directory)'
        )

        # Scoring options
        scoring_group = parser.add_argument_group('scoring options')
        scoring_group.add_argument(
            '-s', '--sequences_only',
            action='store_true',
            help='Extract sequences only, skip classification'
        )
        scoring_group.add_argument(
            '-t', '--threshold',
            type=float,
            default=90.0,
            help='U12 probability threshold (0-100, default: 90)'
        )
        scoring_group.add_argument(
            '-f', '--feature',
            choices=['cds', 'exon', 'both'],
            default='both',
            help='Feature type to extract introns from (default: both)'
        )
        scoring_group.add_argument(
            '--no_nc',
            action='store_true',
            help='Exclude non-canonical introns from scoring'
        )
        scoring_group.add_argument(
            '--pseudocount',
            type=float,
            default=0.0001,
            help='Pseudocount value for PWM scoring to avoid division by zero (default: 0.0001)'
        )
        scoring_group.add_argument(
            '--no_ignore_nc_dnts',
            action='store_true',
            help='Include terminal dinucleotides when scoring non-canonical introns (default: ignore them)'
        )
        scoring_group.add_argument(
            '--no_nc_ss_adjustment',
            action='store_true',
            help='Disable U12 boundary correction for non-canonical introns (default: enabled)'
        )

        # Performance options
        perf_group = parser.add_argument_group('performance options')
        perf_group.add_argument(
            '-p', '--processes',
            type=int,
            default=1,
            help='Number of parallel processes (default: 1)'
        )
        perf_group.add_argument(
            '--min_intron_len',
            type=int,
            default=30,
            help='Minimum intron length in bp (default: 30)'
        )

        # Isoform selection
        isoform_group = parser.add_argument_group('isoform selection')
        isoform_group.add_argument(
            '-i', '--allow_multiple_isoforms',
            action='store_true',
            help='Include non-longest isoforms'
        )
        isoform_group.add_argument(
            '-v', '--no_intron_overlap',
            action='store_true',
            help='Exclude overlapping introns'
        )
        isoform_group.add_argument(
            '-d', '--include_duplicates',
            action='store_true',
            help='Include introns with duplicate coordinates in output (default: exclude)'
        )

        # Training options
        training_group = parser.add_argument_group('training options')
        training_group.add_argument(
            '--pwms',
            type=Path,
            help='Custom PWM matrix file'
        )
        training_group.add_argument(
            '--reference_u12s',
            type=Path,
            help='Custom U12 reference sequences'
        )
        training_group.add_argument(
            '--reference_u2s',
            type=Path,
            help='Custom U2 reference sequences'
        )
        training_group.add_argument(
            '--generate_u2_bps_pwm',
            action='store_true',
            help='Generate U2 branch point PWM from data'
        )
        training_group.add_argument(
            '--recursive',
            nargs='?',
            const=True,
            help='Perform recursive training (optional: U2 subset size)'
        )
        training_group.add_argument(
            '-C',
            type=float,
            help='Fixed SVM C parameter (skips optimization)'
        )
        training_group.add_argument(
            '--n_models',
            type=int,
            default=1,
            help='Number of ensemble models to train (default: 1)'
        )

        # Scoring region coordinates
        coords_group = parser.add_argument_group('scoring region coordinates')
        coords_group.add_argument(
            '--five_score_coords',
            nargs=2,
            type=int,
            default=[-3, 9],
            metavar=('START', 'END'),
            help="5' splice site scoring region (default: -3 9)"
        )
        coords_group.add_argument(
            '--bp_region_coords',
            nargs=2,
            type=int,
            default=[-55, -5],
            metavar=('START', 'END'),
            help='Branch point region coordinates (default: -55 -5)'
        )
        coords_group.add_argument(
            '--three_score_coords',
            nargs=2,
            type=int,
            default=[-6, 4],
            metavar=('START', 'END'),
            help="3' splice site scoring region (default: -6 4)"
        )

        # Advanced options
        advanced_group = parser.add_argument_group('advanced options')
        advanced_group.add_argument(
            '--flank_len',
            type=int,
            default=50,
            help='Length of exonic flanks to extract (default: 50)'
        )
        advanced_group.add_argument(
            '--cv_processes',
            type=int,
            help='Number of processes for cross-validation (default: same as -p)'
        )
        advanced_group.add_argument(
            '--seed',
            type=int,
            default=42,
            help='Random seed for reproducibility (default: 42)'
        )
        advanced_group.add_argument(
            '--quiet',
            action='store_true',
            help='Suppress non-essential output'
        )
        advanced_group.add_argument(
            '--debug',
            action='store_true',
            help='Enable debug logging'
        )

        return parser

    def _validate_args(self, args):
        """Validate parsed arguments.

        Args:
            args: Parsed argument namespace

        Raises:
            argparse.ArgumentTypeError: If validation fails
        """
        # Determine input mode and validate
        has_annotation = args.annotation is not None
        has_bed = args.bed is not None
        has_genome = args.genome is not None
        has_sequences = args.sequence_file is not None

        # Check for genome without annotation/bed first (more specific error)
        if has_genome and not has_annotation and not has_bed:
            self.parser.error(
                "Genome (-g) requires annotation (-a) or BED file (-b)"
            )

        # Count input modes
        input_modes = sum([has_annotation, has_bed, has_sequences])

        if input_modes == 0:
            self.parser.error(
                "No input specified. Choose one mode:\n"
                "  Annotation mode: -g GENOME -a ANNOTATION\n"
                "  BED mode: -g GENOME -b BED\n"
                "  Sequences mode: -q SEQUENCE_FILE"
            )

        if input_modes > 1:
            self.parser.error(
                "Multiple input modes specified. Choose only one:\n"
                "  -a (annotation), -b (bed), or -q (sequences)"
            )

        # Validate annotation/bed mode requires genome
        if (has_annotation or has_bed) and not has_genome:
            mode = "annotation (-a)" if has_annotation else "BED (-b)"
            self.parser.error(f"{mode} mode requires genome file (-g)")

        # Threshold validation
        if not 0 <= args.threshold <= 100:
            self.parser.error("Threshold must be between 0 and 100")

        # Process count validation
        if args.processes < 1:
            self.parser.error("Process count must be >= 1")

        # Set cv_processes to processes if not specified
        if args.cv_processes is None:
            args.cv_processes = args.processes

        # Validate file paths exist
        for attr_name in ['genome', 'annotation', 'bed', 'sequence_file',
                          'pwms', 'reference_u12s', 'reference_u2s']:
            filepath = getattr(args, attr_name)
            if filepath and not filepath.exists():
                self.parser.error(f"File not found: {filepath}")

        # Create output directory if it doesn't exist
        args.output_dir.mkdir(parents=True, exist_ok=True)
