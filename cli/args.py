"""
Argument parser for intronIC CLI with subcommand support.

Supports two modes:
1. train - Train a model on reference data (no genome/annotation needed)
2. classify - Extract and classify introns (requires genome/annotation)

Backward compatible: old CLI without subcommands defaults to 'classify' mode.
"""

import argparse
from pathlib import Path
from typing import Optional


class IntronICArgumentParser:
    """Argument parser for intronIC command-line interface with subcommands."""

    def __init__(self):
        self.parser = self._create_parser()

    def parse_args(self, args: Optional[list] = None):
        """Parse command-line arguments.

        Args:
            args: Optional list of arguments (defaults to sys.argv)

        Returns:
            Namespace object with parsed arguments
        """
        # Helpful error for common mistake: --config after subcommand
        # --config is a global argument and must come BEFORE the subcommand
        import sys
        args_to_check = args if args is not None else sys.argv[1:]
        self._check_config_position(args_to_check)

        parsed = self.parser.parse_args(args)

        # Backward compatibility: if no subcommand, default to classify
        if not hasattr(parsed, 'command') or parsed.command is None:
            parsed.command = 'classify'
            # Map old --pretrained_model to --model for backward compat
            if hasattr(parsed, 'pretrained_model') and parsed.pretrained_model:
                parsed.model = parsed.pretrained_model

        self._validate_args(parsed)
        return parsed

    def _check_config_position(self, args: list):
        """Check if --config appears after subcommand and provide helpful error.

        Args:
            args: Command-line argument list

        Raises:
            SystemExit: If --config appears after subcommand
        """
        if '--config' not in args:
            return

        # Find positions
        config_idx = args.index('--config')
        subcommand_idx = None

        for idx, arg in enumerate(args):
            if arg in ('train', 'classify'):
                subcommand_idx = idx
                break

        # If --config appears after subcommand, show helpful error
        if subcommand_idx is not None and config_idx > subcommand_idx:
            print("\n❌ Error: --config must come BEFORE the subcommand\n")
            print("Incorrect usage:")
            print("  intronIC train --config config/config.yaml -n model  ❌\n")
            print("Correct usage:")
            print("  intronIC --config config/config.yaml train -n model  ✅")
            print("  intronIC --config config/profiles/quick.yaml train -n test  ✅\n")
            print("Note: --config is a global argument and must appear before 'train' or 'classify'\n")
            import sys
            sys.exit(2)

    def _create_parser(self) -> argparse.ArgumentParser:
        """Create the main parser with subcommands."""

        # Get version
        try:
            from importlib.metadata import version
            __version__ = version("intronIC")
        except:
            __version__ = "2.0.0"

        # Main parser
        parser = argparse.ArgumentParser(
            prog='intronIC',
            description="intronIC: Intron classification and extraction tool",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="""
Examples:
  # Train a model on reference data (no genome needed!)
  intronIC train -n homo_sapiens

  # Classify introns with pretrained model
  intronIC classify -g genome.fa -a annotation.gff -n species --model species.model.pkl

  # Classify and train on-the-fly
  intronIC classify -g genome.fa -a annotation.gff -n species --train

  # Backward compatible (no subcommand = classify)
  intronIC -g genome.fa -a annotation.gff -n species --model species.model.pkl
"""
        )

        # Global options (apply to all subcommands)
        parser.add_argument('--version', action='version', version=f'intronIC {__version__}')
        parser.add_argument('--quiet', action='store_true', help='Suppress non-essential output')
        parser.add_argument('--debug', action='store_true', help='Enable debug logging')
        parser.add_argument(
            '--config',
            type=Path,
            dest='config_path',
            help='Path to configuration file (auto-loads from standard paths if not specified)'
        )
        parser.add_argument('--generate-config', action='store_true',
                           help='Generate configuration file template and exit')

        # Create subparsers
        subparsers = parser.add_subparsers(
            dest='command',
            help='Command to run (default: classify if not specified)'
        )

        # ===================================================================
        # TRAIN SUBCOMMAND
        # ===================================================================
        train_parser = subparsers.add_parser(
            'train',
            help='Train a classifier on reference data only (no genome/annotation needed)',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="""
Train mode examples:
  # Basic training with built-in reference data
  intronIC train -n homo_sapiens

  # Parallel training (faster)
  intronIC train -n homo_sapiens -p 8

  # Custom config (recommended for reproducibility)
  # Note: --config is a global arg and must come BEFORE the subcommand
  intronIC --config config/config.yaml train -n homo_sapiens -p 12

  # Use testing profile
  intronIC --config config/profiles/quick.yaml train -n test_model

  # Custom reference sequences
  intronIC train -n species --reference_u12s custom_u12.iic --reference_u2s custom_u2.iic

  # Custom training parameters
  intronIC train -n species --n_models 8 --n_cv_folds 10 -p 12

  # Fast training (skip optimization)
  intronIC train -n species -C 0.1 --eval_mode none
"""
        )
        self._add_train_arguments(train_parser)

        # ===================================================================
        # CLASSIFY SUBCOMMAND
        # ===================================================================
        classify_parser = subparsers.add_parser(
            'classify',
            help='Extract and classify introns from genome/annotation',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="""
Classify mode examples:
  # With pretrained model
  intronIC classify -g genome.fa -a annotation.gff -n species --model species.model.pkl

  # Train on-the-fly
  intronIC classify -g genome.fa -a annotation.gff -n species --train

  # Extract sequences only (no classification)
  intronIC classify -g genome.fa -a annotation.gff -n species -s

  # From BED file
  intronIC classify -g genome.fa -b introns.bed -n species --model species.model.pkl
"""
        )
        self._add_classify_arguments(classify_parser)

        # ===================================================================
        # BACKWARD COMPATIBILITY
        # ===================================================================
        # Add all classify arguments to main parser for old CLI compatibility
        # This allows: intronIC -g genome.fa -a annotation.gff -n species
        self._add_classify_arguments(parser, for_backward_compat=True)

        return parser

    def _add_common_arguments(self, parser: argparse.ArgumentParser):
        """Add arguments common to both train and classify."""

        # Species name (required for both)
        parser.add_argument(
            '-n', '--species_name', '--species-name',
            required=True,
            help='Species name for output files (e.g., homo_sapiens)'
        )

        # Output directory
        parser.add_argument(
            '-o', '--output_dir', '--output-dir',
            type=Path,
            default=Path.cwd(),
            help='Output directory (default: current directory)'
        )

    def _add_train_arguments(self, parser: argparse.ArgumentParser):
        """Add arguments specific to train subcommand."""

        # Common arguments
        self._add_common_arguments(parser)

        # === Reference Data ===
        reference = parser.add_argument_group('reference data')
        reference.add_argument(
            '--reference_u12s', '--reference-u12s',
            type=Path,
            help='Custom U12 reference sequences (.iic format). Default: use built-in reference.'
        )
        reference.add_argument(
            '--reference_u2s', '--reference-u2s',
            type=Path,
            help='Custom U2 reference sequences (.iic format). Default: use built-in reference.'
        )

        # === Training Parameters ===
        training = parser.add_argument_group('training parameters')
        training.add_argument(
            '-C',
            type=float,
            help='Fixed SVM C parameter (skips hyperparameter optimization)'
        )
        training.add_argument(
            '--n_models', '--n-models',
            type=int,
            default=1,
            help='Number of ensemble models to train (default: 1, use YAML config for production)'
        )
        training.add_argument(
            '--eval_mode', '--eval-mode',
            choices=['nested_cv', 'split', 'none'],
            default='nested_cv',
            help='Evaluation mode: nested_cv (default), split, or none'
        )
        training.add_argument(
            '--n_cv_folds', '--n-cv-folds',
            type=int,
            default=5,
            help='Number of cross-validation folds (default: 5)'
        )
        training.add_argument(
            '--test_fraction', '--test-fraction',
            type=float,
            default=0.2,
            help='Test set fraction for split mode (default: 0.2)'
        )
        training.add_argument(
            '--n_optimization_rounds', '--n-optimization-rounds',
            type=int,
            default=5,
            help='Hyperparameter optimization rounds (default: 5)'
        )
        training.add_argument(
            '--max_iter', '--max-iter',
            type=int,
            default=50000,
            help='Maximum SVM iterations (default: 50000)'
        )
        training.add_argument(
            '--cv_processes', '--cv-processes',
            type=int,
            default=1,
            help='Processes for cross-validation (default: 1)'
        )
        training.add_argument(
            '-p', '--processes',
            type=int,
            default=1,
            help='Number of parallel processes for scoring reference sequences (default: 1)'
        )
        training.add_argument(
            '--use-fold-averaged-params',
            action='store_true',
            default=None,
            help='Use fold-averaged hyperparameters from nested CV instead of re-optimizing on full dataset (better cross-species generalization)'
        )

        # === Scoring Parameters ===
        # (Needed to score reference sequences during training)
        scoring = parser.add_argument_group('scoring parameters (for reference sequences)')
        scoring.add_argument(
            '--five_score_coords', '--five-score-coords',
            nargs=2,
            type=int,
            default=[-3, 9],
            metavar=('START', 'END'),
            help="5' splice site region (default: -3 9)"
        )
        scoring.add_argument(
            '--bp_region_coords', '--bp-region-coords',
            nargs=2,
            type=int,
            default=[-55, -5],
            metavar=('START', 'END'),
            help='Branch point region (default: -55 -5)'
        )
        scoring.add_argument(
            '--three_score_coords', '--three-score-coords',
            nargs=2,
            type=int,
            default=[-6, 4],
            metavar=('START', 'END'),
            help="3' splice site region (default: -6 4)"
        )
        scoring.add_argument(
            '--pseudocount',
            type=float,
            default=0.0001,
            help='PWM pseudocount (default: 0.0001)'
        )

        # === Advanced ===
        advanced = parser.add_argument_group('advanced options')
        advanced.add_argument(
            '--seed',
            type=int,
            default=42,
            help='Random seed for reproducibility (default: 42)'
        )

    def _add_classify_arguments(self, parser: argparse.ArgumentParser, for_backward_compat=False):
        """Add arguments specific to classify subcommand.

        Args:
            parser: Parser to add arguments to
            for_backward_compat: If True, makes species_name not required (for main parser)
        """

        # Common arguments (species name required only if not backward compat)
        if not for_backward_compat:
            self._add_common_arguments(parser)
        else:
            # For backward compat, don't require species_name yet (validated later)
            parser.add_argument(
                '-n', '--species_name', '--species-name',
                help='Species name for output files (e.g., homo_sapiens)'
            )
            parser.add_argument(
                '-o', '--output_dir', '--output-dir',
                type=Path,
                default=Path.cwd(),
                help='Output directory (default: current directory)'
            )

        # === Input Selection ===
        input_group = parser.add_argument_group(
            'input selection',
            'Choose one mode: (1) -g + -a for annotation, (2) -g + -b for BED, or (3) -q for sequences'
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
            '-q', '--sequence_file', '--sequence-file',
            type=Path,
            help='Path to pre-extracted intron sequences (.iic format)'
        )

        # === Model Source ===
        model_group = parser.add_argument_group('model source (choose one)')
        model_exclusive = model_group.add_mutually_exclusive_group()
        model_exclusive.add_argument(
            '--model',
            type=Path,
            help='Path to pretrained model (.model.pkl)'
        )
        model_exclusive.add_argument(
            '--train',
            action='store_true',
            help='Train model on-the-fly (slower, includes full training and evaluation)'
        )

        # Normalizer mode (for pretrained model classification)
        parser.add_argument(
            '--normalizer_mode', '--normalizer-mode',
            choices=['human', 'adaptive', 'auto'],
            default='auto',
            help='''Normalizer mode for pretrained model classification (default: auto):
              human: Use scaler from training species (recommended for U12-absent genomes)
              adaptive: Refit scaler on experimental data (experimental, may cause FPs in U12-free species)
              auto: Use human if available in model, otherwise adaptive'''
        )

        # Backward compatibility: old --pretrained_model flag
        if for_backward_compat:
            model_group.add_argument(
                '--pretrained_model', '--pretrained-model',
                type=Path,
                help='(Deprecated: use --model) Path to pretrained model'
            )

        # === Extraction Parameters ===
        extraction = parser.add_argument_group('extraction parameters')
        extraction.add_argument(
            '-f', '--feature',
            choices=['cds', 'exon', 'both'],
            default='both',
            help='Feature type to extract from (default: both)'
        )
        extraction.add_argument(
            '--min_intron_len', '--min-intron-len',
            type=int,
            default=30,
            help='Minimum intron length (default: 30)'
        )
        extraction.add_argument(
            '-i', '--allow_multiple_isoforms', '--allow-multiple-isoforms',
            action='store_true',
            help='Include non-longest isoforms'
        )
        extraction.add_argument(
            '-v', '--no_intron_overlap', '--no-intron-overlap',
            action='store_true',
            help='Exclude overlapping introns'
        )
        extraction.add_argument(
            '-d', '--include_duplicates', '--include-duplicates',
            action='store_true',
            help='Include duplicate coordinate introns'
        )
        extraction.add_argument(
            '--flank_len', '--flank-len',
            type=int,
            default=50,
            help='Exonic flank length (default: 50)'
        )
        extraction.add_argument(
            '--no_nc_ss_adjustment', '--no-nc-ss-adjustment',
            action='store_true',
            help='Disable U12 boundary correction'
        )

        # === Scoring Parameters ===
        scoring = parser.add_argument_group('scoring parameters')
        scoring.add_argument(
            '-s', '--sequences_only', '--sequences-only',
            action='store_true',
            help='Extract sequences only, skip classification'
        )
        scoring.add_argument(
            '-t', '--threshold',
            type=float,
            default=90.0,
            help='U12 probability threshold 0-100 (default: 90)'
        )
        scoring.add_argument(
            '--no_nc', '--no-nc',
            action='store_true',
            help='Exclude non-canonical introns from scoring'
        )
        scoring.add_argument(
            '--pseudocount',
            type=float,
            default=0.0001,
            help='PWM pseudocount (default: 0.0001)'
        )
        scoring.add_argument(
            '--no_ignore_nc_dnts', '--no-ignore-nc-dnts',
            action='store_true',
            help='Include terminal dinucleotides in non-canonical scoring'
        )
        scoring.add_argument(
            '--five_score_coords', '--five-score-coords',
            nargs=2,
            type=int,
            default=[-3, 9],
            metavar=('START', 'END'),
            help="5' splice site region (default: -3 9)"
        )
        scoring.add_argument(
            '--bp_region_coords', '--bp-region-coords',
            nargs=2,
            type=int,
            default=[-55, -5],
            metavar=('START', 'END'),
            help='Branch point region (default: -55 -5)'
        )
        scoring.add_argument(
            '--three_score_coords', '--three-score-coords',
            nargs=2,
            type=int,
            default=[-6, 4],
            metavar=('START', 'END'),
            help="3' splice site region (default: -6 4)"
        )

        # === Training Parameters (only if --train flag used) ===
        training = parser.add_argument_group('training parameters (only with --train)')
        training.add_argument(
            '--reference_u12s', '--reference-u12s',
            type=Path,
            help='Custom U12 reference sequences'
        )
        training.add_argument(
            '--reference_u2s', '--reference-u2s',
            type=Path,
            help='Custom U2 reference sequences'
        )
        training.add_argument(
            '--pwms',
            type=Path,
            help='Custom PWM matrix file'
        )
        training.add_argument(
            '--generate_u2_bps_pwm', '--generate-u2-bps-pwm',
            action='store_true',
            help='Generate U2 branch point PWM from data'
        )
        training.add_argument(
            '--recursive',
            nargs='?',
            const=True,
            help='Perform recursive training'
        )
        training.add_argument(
            '-C',
            type=float,
            help='Fixed SVM C parameter (skips optimization)'
        )
        training.add_argument(
            '--n_models', '--n-models',
            type=int,
            default=1,
            help='Number of ensemble models (default: 1)'
        )
        training.add_argument(
            '--max_iter', '--max-iter',
            type=int,
            default=50000,
            help='Maximum SVM iterations (default: 50000)'
        )
        training.add_argument(
            '--eval_mode', '--eval-mode',
            choices=['nested_cv', 'split', 'none'],
            default='nested_cv',
            help='Evaluation mode (default: nested_cv)'
        )
        training.add_argument(
            '--n_cv_folds', '--n-cv-folds',
            type=int,
            default=5,
            help='Cross-validation folds (default: 5)'
        )
        training.add_argument(
            '--test_fraction', '--test-fraction',
            type=float,
            default=0.2,
            help='Test fraction for split mode (default: 0.2)'
        )
        training.add_argument(
            '--n_optimization_rounds', '--n-optimization-rounds',
            type=int,
            default=5,
            help='Optimization rounds (default: 5)'
        )
        training.add_argument(
            '--use-fold-averaged-params',
            action='store_true',
            default=None,
            help='Use fold-averaged hyperparameters from nested CV instead of re-optimizing on full dataset (better cross-species generalization)'
        )

        # === Performance ===
        perf = parser.add_argument_group('performance options')
        perf.add_argument(
            '-p', '--processes',
            type=int,
            default=1,
            help='Parallel processes for scoring (default: 1)'
        )
        perf.add_argument(
            '--cv_processes', '--cv-processes',
            type=int,
            help='Processes for cross-validation (default: same as -p)'
        )

        # === Output Options ===
        output = parser.add_argument_group('output options')
        output.add_argument(
            '--clean_names', '--clean-names',
            action='store_true',
            default=True,
            help='Remove "transcript:" and "gene:" prefixes (default: True)'
        )
        output.add_argument(
            '--no_clean_names', '--no-clean-names',
            dest='clean_names',
            action='store_false',
            help='Keep ID prefixes'
        )
        output.add_argument(
            '-u', '--uninformative_naming', '--uninformative-naming',
            action='store_true',
            help='Use simple naming scheme'
        )
        output.add_argument(
            '--no_abbreviate', '--no-abbreviate', '--na',
            action='store_true',
            help='Use full species name in outputs'
        )
        output.add_argument(
            '--abbreviate_filenames', '--abbreviate-filenames', '--afn',
            action='store_true',
            help='Abbreviate species name in filenames'
        )

        # === Advanced ===
        advanced = parser.add_argument_group('advanced options')
        advanced.add_argument(
            '--seed',
            type=int,
            default=42,
            help='Random seed (default: 42)'
        )

    def _validate_args(self, args):
        """Validate parsed arguments based on command.

        Args:
            args: Parsed argument namespace

        Raises:
            argparse.ArgumentTypeError: If validation fails
        """
        # Skip validation if generating config
        if getattr(args, 'generate_config', False):
            return

        # Get command (train or classify)
        command = getattr(args, 'command', 'classify')

        # ===============================================================
        # TRAIN MODE VALIDATION
        # ===============================================================
        if command == 'train':
            # Species name required
            if not args.species_name:
                self.parser.error("train: -n/--species_name is required")

            # Validate CV parameters
            if args.n_cv_folds < 2:
                self.parser.error("train: n_cv_folds must be >= 2")

            if not 0 < args.test_fraction < 1:
                self.parser.error("train: test_fraction must be between 0 and 1")

            # Validate custom reference files exist
            if args.reference_u12s and not args.reference_u12s.exists():
                self.parser.error(f"train: reference U12 file not found: {args.reference_u12s}")
            if args.reference_u2s and not args.reference_u2s.exists():
                self.parser.error(f"train: reference U2 file not found: {args.reference_u2s}")

        # ===============================================================
        # CLASSIFY MODE VALIDATION
        # ===============================================================
        elif command == 'classify':
            # Species name required
            if not args.species_name:
                self.parser.error("classify: -n/--species_name is required")

            # Input validation
            has_annotation = args.annotation is not None
            has_bed = args.bed is not None
            has_genome = args.genome is not None
            has_sequences = args.sequence_file is not None

            # Check for genome without annotation/bed
            if has_genome and not has_annotation and not has_bed:
                self.parser.error(
                    "classify: genome (-g) requires annotation (-a) or BED file (-b)"
                )

            # Count input modes
            input_modes = sum([has_annotation, has_bed, has_sequences])

            if input_modes == 0:
                self.parser.error(
                    "classify: no input specified. Choose one:\n"
                    "  Annotation: -g GENOME -a ANNOTATION\n"
                    "  BED: -g GENOME -b BED\n"
                    "  Sequences: -q SEQUENCE_FILE"
                )

            if input_modes > 1:
                self.parser.error(
                    "classify: multiple inputs specified. Choose only one: -a, -b, or -q"
                )

            # Validate annotation/bed requires genome
            if (has_annotation or has_bed) and not has_genome:
                mode = "annotation (-a)" if has_annotation else "BED (-b)"
                self.parser.error(f"classify: {mode} requires genome file (-g)")

            # Model source validation (skip if sequences_only mode)
            if not args.sequences_only:
                has_model = hasattr(args, 'model') and args.model is not None
                has_train = hasattr(args, 'train') and args.train

                if not has_model and not has_train:
                    self.parser.error(
                        "classify: must specify model source:\n"
                        "  --model PATH  (use pretrained model)\n"
                        "  --train       (train on-the-fly)"
                    )

            # Validate file paths exist
            file_attrs = ['genome', 'annotation', 'bed', 'sequence_file', 'model',
                         'pwms', 'reference_u12s', 'reference_u2s', 'optimizer_config']
            for attr_name in file_attrs:
                if hasattr(args, attr_name):
                    filepath = getattr(args, attr_name)
                    if filepath and not filepath.exists():
                        self.parser.error(f"classify: file not found: {filepath}")

            # Threshold validation
            if not 0 <= args.threshold <= 100:
                self.parser.error("classify: threshold must be between 0 and 100")

            # Training parameters validation (if --train used)
            if hasattr(args, 'train') and args.train:
                if args.n_cv_folds < 2:
                    self.parser.error("classify: n_cv_folds must be >= 2")
                if not 0 < args.test_fraction < 1:
                    self.parser.error("classify: test_fraction must be between 0 and 1")

            # Set cv_processes to processes if not specified
            if args.cv_processes is None:
                args.cv_processes = args.processes

            # Process count validation
            if args.processes < 1:
                self.parser.error("classify: processes must be >= 1")

        # Create output directory
        args.output_dir.mkdir(parents=True, exist_ok=True)
