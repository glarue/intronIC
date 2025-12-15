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

        # Store argparse defaults in the namespace for config merge logic
        # This provides a single source of truth for default values
        parsed._arg_defaults = {}
        for action in self.parser._actions:
            if action.dest != 'help' and hasattr(action, 'default'):
                parsed._arg_defaults[action.dest] = action.default

        # Backward compatibility: if no subcommand, default to classify
        if not hasattr(parsed, "command") or parsed.command is None:
            parsed.command = "classify"
            # Map old --pretrained_model to --model for backward compat
            if hasattr(parsed, "pretrained_model") and parsed.pretrained_model:
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
        if "--config" not in args:
            return

        # Find positions
        config_idx = args.index("--config")
        subcommand_idx = None

        for idx, arg in enumerate(args):
            if arg in ("train", "classify", "extract", "test"):
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
            print(
                "Note: --config is a global argument and must appear before 'train' or 'classify'\n"
            )
            import sys

            sys.exit(2)

    def _create_parser(self) -> argparse.ArgumentParser:
        """Create the main parser with subcommands."""

        # Get version
        try:
            from importlib.metadata import version

            __version__ = version("intronIC")
        except:
            __version__ = "dev"

        # Main parser
        parser = argparse.ArgumentParser(
            prog="intronIC",
            description="intronIC: Intron classification and extraction tool",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="""
Examples:
  # Train a model on reference data (no genome needed!)
  intronIC train -n homo_sapiens

  # Classify introns with pretrained model (uses streaming mode by default)
  intronIC classify -g genome.fa -a annotation.gff -n species --model species.model.pkl

  # Extract introns without classification
  intronIC extract -g genome.fa -a annotation.gff -n species

  # Use in-memory mode for small genomes (higher memory, potentially faster)
  intronIC classify -g genome.fa -a annotation.gff -n species --in-memory

  # Backward compatible (no subcommand = classify)
  intronIC -g genome.fa -a annotation.gff -n species --model species.model.pkl
""",
        )

        # Global options (apply to all subcommands)
        parser.add_argument(
            "--version", action="version", version=f"intronIC {__version__}"
        )
        parser.add_argument(
            "--quiet", action="store_true", help="Suppress non-essential output"
        )
        parser.add_argument("--debug", action="store_true", help="Enable debug logging")
        parser.add_argument(
            "--config",
            type=Path,
            dest="config_path",
            help="Path to configuration file (auto-loads from standard paths if not specified)",
        )
        parser.add_argument(
            "--generate-config",
            action="store_true",
            help="Generate configuration file template and exit",
        )

        # Create subparsers
        subparsers = parser.add_subparsers(
            dest="command", help="Command to run (default: classify if not specified)"
        )

        # ===================================================================
        # TRAIN SUBCOMMAND
        # ===================================================================
        train_parser = subparsers.add_parser(
            "train",
            help="Train a classifier on reference data only (no genome/annotation needed)",
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
""",
        )
        self._add_train_arguments(train_parser)

        # ===================================================================
        # CLASSIFY SUBCOMMAND
        # ===================================================================
        classify_parser = subparsers.add_parser(
            "classify",
            help="Extract and classify introns from genome/annotation",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="""
Classify mode examples:
  # With pretrained model
  intronIC classify -g genome.fa -a annotation.gff -n species --model species.model.pkl

  # From BED file
  intronIC classify -g genome.fa -b introns.bed -n species --model species.model.pkl

  # Use default model if available
  intronIC classify -g genome.fa -a annotation.gff -n species
""",
        )
        self._add_classify_arguments(classify_parser)

        # ===================================================================
        # EXTRACT SUBCOMMAND
        # ===================================================================
        extract_parser = subparsers.add_parser(
            "extract",
            help="Extract intron sequences without classification",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="""
Extract mode examples:
  # Extract introns from annotation (uses streaming mode by default)
  intronIC extract -g genome.fa -a annotation.gff -n species

  # Extract with in-memory mode (higher memory usage)
  intronIC extract -g genome.fa -a annotation.gff -n species --in-memory

  # Extract from BED file
  intronIC extract -g genome.fa -b introns.bed -n species

  # Extract with custom flanking length
  intronIC extract -g genome.fa -a annotation.gff -n species --flank-len 20

Note: This command extracts intron sequences but does not perform classification.
      Use 'intronIC classify' with a pretrained model to classify extracted introns.
""",
        )
        self._add_extract_arguments(extract_parser)

        # ===================================================================
        # TEST SUBCOMMAND
        # ===================================================================
        test_parser = subparsers.add_parser(
            "test",
            help="Run a quick installation test with bundled test data",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="""
Test mode examples:
  # Show test data location only
  intronIC test --show-only

  # Run quick classification test (default)
  intronIC test

  # Run test with custom number of processes
  intronIC test -p 4

Note: This command uses the bundled test data (Homo sapiens Chr19) to verify
      your intronIC installation is working correctly (runtime: ~1 min with -p 4)
""",
        )
        self._add_test_arguments(test_parser)

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
            "-n",
            "--species-name",
            "--species_name",
            required=True,
            help="Species name for output files (e.g., homo_sapiens)",
        )

        # Output directory
        parser.add_argument(
            "-o",
            "--output-dir",
            "--output_dir",
            type=Path,
            default=Path.cwd(),
            help="Output directory (default: current directory)",
        )

    def _add_extract_arguments(self, parser: argparse.ArgumentParser):
        """Add arguments specific to extract subcommand.

        Extract mode extracts intron sequences without classification.
        It takes similar arguments to classify but without model/training/scoring params.
        """
        # Common arguments (species name and output dir)
        self._add_common_arguments(parser)

        # === Input Selection ===
        input_group = parser.add_argument_group(
            "input selection",
            "Choose one mode: (1) -g + -a for annotation, or (2) -g + -b for BED",
        )
        input_group.add_argument(
            "-g",
            "--genome",
            type=Path,
            required=True,
            help="Path to genome FASTA file (required)",
        )
        input_group.add_argument(
            "-a",
            "--annotation",
            type=Path,
            help="Path to GFF3/GTF annotation file (requires -g)",
        )
        input_group.add_argument(
            "-b",
            "--bed",
            type=Path,
            help="Path to BED file with intron coordinates (requires -g)",
        )

        # === Extraction Parameters ===
        extraction = parser.add_argument_group("extraction parameters")
        extraction.add_argument(
            "-f",
            "--feature",
            choices=["cds", "exon", "both"],
            default="both",
            help="Feature type to extract from (default: both)",
        )
        extraction.add_argument(
            "--min-intron-len",
            "--min_intron_len",
            type=int,
            default=30,
            help="Minimum intron length (default: 30)",
        )
        extraction.add_argument(
            "-i",
            "--allow-multiple-isoforms",
            "--allow_multiple_isoforms",
            action="store_true",
            help="Include non-longest isoforms",
        )
        extraction.add_argument(
            "-v",
            "--no-intron-overlap",
            "--no_intron_overlap",
            action="store_true",
            help="Exclude overlapping introns",
        )
        extraction.add_argument(
            "-d",
            "--include-duplicates",
            "--include_duplicates",
            action="store_true",
            help="Include duplicate coordinate introns",
        )
        extraction.add_argument(
            "--flank-len",
            "--flank_len",
            type=int,
            default=100,
            help="Exonic flank length (default: 100)",
        )

        # === Performance Parameters ===
        performance = parser.add_argument_group("performance parameters")
        performance.add_argument(
            "-p",
            "--processes",
            type=int,
            default=None,
            help="Number of parallel processes (default: 1)",
        )

        # Memory mode (streaming is default)
        memory_mode = performance.add_mutually_exclusive_group()
        memory_mode.add_argument(
            "--in-memory",
            "--no-streaming",
            action="store_false",
            dest="streaming",
            help="Use in-memory mode (loads full genome, higher memory usage)",
        )
        memory_mode.add_argument(
            "--streaming",
            action="store_true",
            dest="streaming",
            default=True,
            help="Use streaming mode (default: ~85%% memory savings, faster with -p)",
        )

        # === Output Parameters ===
        output = parser.add_argument_group("output parameters")
        output.add_argument(
            "--clean-names",
            action="store_true",
            help="Remove version numbers and special characters from gene/transcript names",
        )
        output.add_argument(
            "--no-clean-names",
            action="store_true",
            help="Keep original names (overrides config file)",
        )
        output.add_argument(
            "-u",
            "--no-abbreviate",
            "--no_abbreviate",
            action="store_true",
            help="Use full gene/transcript names (no abbreviation)",
        )
        output.add_argument(
            "--abbreviate-filenames",
            "--abbreviate_filenames",
            action="store_true",
            help="Use species abbreviation in filenames",
        )
        output.add_argument(
            "--no-headers",
            "--no_headers",
            action="store_true",
            help="Omit headers from output files",
        )
        output.add_argument(
            "--seed",
            type=int,
            default=42,
            help="Random seed for reproducibility (default: 42)",
        )

    def _add_test_arguments(self, parser: argparse.ArgumentParser):
        """Add arguments specific to test subcommand."""

        parser.add_argument(
            "--show-only",
            "--show_only",
            action="store_true",
            help="Show test data location without running test",
        )
        parser.add_argument(
            "-p",
            "--processes",
            type=int,
            default=1,
            help="Number of parallel processes (default: 1)",
        )
        parser.add_argument(
            "-o",
            "--output-dir",
            "--output_dir",
            type=Path,
            default=None,
            help="Output directory (default: temporary directory)",
        )

    def _add_train_arguments(self, parser: argparse.ArgumentParser):
        """Add arguments specific to train subcommand."""

        # Common arguments
        self._add_common_arguments(parser)

        # === Reference Data ===
        reference = parser.add_argument_group("reference data")
        reference.add_argument(
            "--reference-u12s",
            "--reference_u12s",
            type=Path,
            help="Custom U12 reference sequences (.iic format). Default: use built-in reference.",
        )
        reference.add_argument(
            "--reference-u2s",
            "--reference_u2s",
            type=Path,
            help="Custom U2 reference sequences (.iic format). Default: use built-in reference.",
        )

        # === Training Parameters ===
        training = parser.add_argument_group("training parameters")
        training.add_argument(
            "-C",
            type=float,
            help="Fixed SVM C parameter (skips hyperparameter optimization)",
        )
        training.add_argument(
            "--n-models",
            "--n_models",
            type=int,
            default=1,
            help="Number of ensemble models to train (default: 1, use YAML config for production)",
        )
        training.add_argument(
            "--eval-mode",
            "--eval_mode",
            choices=["nested_cv", "split", "none"],
            default="nested_cv",
            help="Evaluation mode: nested_cv (default), split, or none",
        )
        training.add_argument(
            "--n-cv-folds",
            "--n_cv_folds",
            type=int,
            default=5,
            help="Number of cross-validation folds (default: 5)",
        )
        training.add_argument(
            "--test-fraction",
            "--test_fraction",
            type=float,
            default=0.2,
            help="Test set fraction for split mode (default: 0.2)",
        )
        training.add_argument(
            "--n-optimization-rounds",
            "--n_optimization_rounds",
            type=int,
            default=5,
            help="Hyperparameter optimization rounds (default: 5)",
        )
        training.add_argument(
            "--max-iter",
            "--max_iter",
            type=int,
            default=50000,
            help="Maximum SVM iterations (default: 50000)",
        )
        training.add_argument(
            "--cv-processes",
            "--cv_processes",
            type=int,
            default=1,
            help="Processes for cross-validation (default: 1)",
        )
        training.add_argument(
            "-p",
            "--processes",
            type=int,
            default=None,
            help="Number of parallel processes for scoring reference sequences (default: 1)",
        )
        training.add_argument(
            "--use-fold-averaged-params",
            action="store_true",
            default=None,
            help="Use fold-averaged hyperparameters from nested CV instead of re-optimizing on full dataset (better cross-species generalization)",
        )

        # === Scoring Parameters ===
        # (Needed to score reference sequences during training)
        scoring = parser.add_argument_group(
            "scoring parameters (for reference sequences)"
        )
        scoring.add_argument(
            "--five-score-coords",
            "--five_score_coords",
            nargs=2,
            type=int,
            default=[-3, 9],
            metavar=("START", "END"),
            help="5' splice site region (default: -3 9)",
        )
        scoring.add_argument(
            "--bp-region-coords",
            "--bp_region_coords",
            nargs=2,
            type=int,
            default=[-55, -5],
            metavar=("START", "END"),
            help="Branch point region (default: -55 -5)",
        )
        scoring.add_argument(
            "--three-score-coords",
            "--three_score_coords",
            nargs=2,
            type=int,
            default=[-6, 4],
            metavar=("START", "END"),
            help="3' splice site region (default: -6 4)",
        )
        scoring.add_argument(
            "--pseudocount",
            type=float,
            default=0.0001,
            help="PWM pseudocount (default: 0.0001)",
        )
        scoring.add_argument(
            "--pwms", type=Path, help="Custom PWM matrix file (.iic or .yaml format)"
        )

        # === Advanced ===
        advanced = parser.add_argument_group("advanced options")
        advanced.add_argument(
            "--seed",
            type=int,
            default=42,
            help="Random seed for reproducibility (default: 42)",
        )

    def _add_classify_arguments(
        self, parser: argparse.ArgumentParser, for_backward_compat=False
    ):
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
                "-n",
                "--species-name",
                "--species_name",
                help="Species name for output files (e.g., homo_sapiens)",
            )
            parser.add_argument(
                "-o",
                "--output-dir",
                "--output_dir",
                type=Path,
                default=Path.cwd(),
                help="Output directory (default: current directory)",
            )

        # === Input Selection ===
        input_group = parser.add_argument_group(
            "input selection",
            "Choose one mode: (1) -g + -a for annotation, (2) -g + -b for BED, or (3) -q for sequences",
        )
        input_group.add_argument(
            "-g",
            "--genome",
            type=Path,
            help="Path to genome FASTA file (required with -a or -b)",
        )
        input_group.add_argument(
            "-a",
            "--annotation",
            type=Path,
            help="Path to GFF3/GTF annotation file (requires -g)",
        )
        input_group.add_argument(
            "-b",
            "--bed",
            type=Path,
            help="Path to BED file with intron coordinates (requires -g)",
        )
        input_group.add_argument(
            "-q",
            "--sequence-file",
            "--sequence_file",
            type=Path,
            help="Path to pre-extracted intron sequences (.iic format)",
        )

        # === Model Source ===
        model_group = parser.add_argument_group("model source")
        model_group.add_argument(
            "--model", type=Path, help="Path to pretrained model (.model.pkl)"
        )

        # Normalizer mode (for pretrained model classification)
        parser.add_argument(
            "--normalizer-mode",
            "--normalizer_mode",
            choices=["human", "adaptive", "auto"],
            default="auto",
            help="""Normalizer mode for pretrained model classification (default: auto):
              human: Use scaler from training species (recommended for U12-absent genomes)
              adaptive: Refit scaler on experimental data (experimental, may cause FPs in U12-free species)
              auto: Use human if available in model, otherwise adaptive""",
        )

        # Species-specific prior adjustment (for U12-absent species)
        parser.add_argument(
            "--species-prior",
            type=float,
            default=None,
            metavar="PRIOR",
            help="""Expected U12 prior for target species (0 to 1). Adjusts classification
              probabilities via Bayes rule to account for different U12 base rates.
              By default (None), no adjustment is applied - classifier uses raw probabilities
              from the model trained on human data without making assumptions about target species.
              Recommended values when adjustment is desired:
                - 0.005: Human-like species (similar to training data)
                - 1e-6: U12-absent species (C. elegans, many fungi)
                - 1e-4: U12-poor species
              Lower values reduce false positives in U12-free lineages.""",
        )

        # Load saved normalizer (for reproducible adaptive normalization)
        parser.add_argument(
            "--load-normalizer",
            type=Path,
            default=None,
            metavar="PATH",
            help="""Load a saved normalizer from a previous run (for reproducible normalization).
              When using adaptive mode on a species, the first run on the full genome should
              fit and save a normalizer. Subsequent runs on subsets can reuse this normalizer
              for consistency. Only applies when using --normalizer-mode adaptive.
              The normalizer is automatically saved as <output_prefix>.normalizer.pkl""",
        )

        # Save fitted normalizer (for future reproducibility)
        parser.add_argument(
            "--save-normalizer",
            action="store_true",
            default=False,
            help="""Save the fitted normalizer for future runs (adaptive mode only).
              Use this on your first full-genome run for a species to establish a reference
              normalization. Future runs can use --load-normalizer to reuse this normalization.
              Saved to <output_prefix>.normalizer.pkl""",
        )

        # Backward compatibility: old --pretrained_model flag
        if for_backward_compat:
            model_group.add_argument(
                "--pretrained-model",
                "--pretrained_model",
                type=Path,
                help="(Deprecated: use --model) Path to pretrained model",
            )

        # === Extraction Parameters ===
        extraction = parser.add_argument_group("extraction parameters")
        extraction.add_argument(
            "-f",
            "--feature",
            choices=["cds", "exon", "both"],
            default="both",
            help="Feature type to extract from (default: both)",
        )
        extraction.add_argument(
            "--min-intron-len",
            "--min_intron_len",
            type=int,
            default=30,
            help="Minimum intron length (default: 30)",
        )
        extraction.add_argument(
            "-i",
            "--allow-multiple-isoforms",
            "--allow_multiple_isoforms",
            action="store_true",
            help="Include non-longest isoforms",
        )
        extraction.add_argument(
            "-v",
            "--no-intron-overlap",
            "--no_intron_overlap",
            action="store_true",
            help="Exclude overlapping introns",
        )
        extraction.add_argument(
            "-d",
            "--include-duplicates",
            "--include_duplicates",
            action="store_true",
            help="Include duplicate coordinate introns",
        )
        extraction.add_argument(
            "--flank-len",
            "--flank_len",
            type=int,
            default=100,
            help="Exonic flank length (default: 100)",
        )
        extraction.add_argument(
            "--no-nc-ss-adjustment",
            "--no_nc_ss_adjustment",
            action="store_true",
            help="Disable U12 boundary correction",
        )

        # === Scoring Parameters ===
        scoring = parser.add_argument_group("scoring parameters")
        scoring.add_argument(
            "-t",
            "--threshold",
            type=float,
            default=90.0,
            help="U12 probability threshold 0-100 (default: 90)",
        )
        scoring.add_argument(
            "--no-nc",
            "--no_nc",
            action="store_true",
            help="Exclude non-canonical introns from scoring",
        )
        scoring.add_argument(
            "--pseudocount",
            type=float,
            default=0.0001,
            help="PWM pseudocount (default: 0.0001)",
        )
        scoring.add_argument(
            "--no-ignore-nc-dnts",
            "--no_ignore_nc_dnts",
            action="store_true",
            help="Include terminal dinucleotides in non-canonical scoring",
        )
        scoring.add_argument(
            "--five-score-coords",
            "--five_score_coords",
            nargs=2,
            type=int,
            default=[-3, 9],
            metavar=("START", "END"),
            help="5' splice site region (default: -3 9)",
        )
        scoring.add_argument(
            "--bp-region-coords",
            "--bp_region_coords",
            nargs=2,
            type=int,
            default=[-55, -5],
            metavar=("START", "END"),
            help="Branch point region (default: -55 -5)",
        )
        scoring.add_argument(
            "--three-score-coords",
            "--three_score_coords",
            nargs=2,
            type=int,
            default=[-6, 4],
            metavar=("START", "END"),
            help="3' splice site region (default: -6 4)",
        )

        # === Performance ===
        perf = parser.add_argument_group("performance options")
        perf.add_argument(
            "-p",
            "--processes",
            type=int,
            default=None,
            help="Parallel processes for scoring (default: 1)",
        )
        perf.add_argument(
            "--cv-processes",
            "--cv_processes",
            type=int,
            help="Processes for cross-validation (default: same as -p)",
        )

        # Memory mode (streaming is default)
        memory_mode = perf.add_mutually_exclusive_group()
        memory_mode.add_argument(
            "--in-memory",
            "--no-streaming",
            action="store_false",
            dest="streaming",
            help="Use in-memory mode: load full genome into memory. Higher memory usage "
            "but may be slightly faster for very small genomes.",
        )
        memory_mode.add_argument(
            "--streaming",
            action="store_true",
            dest="streaming",
            default=True,
            help="Use streaming mode (default): stores sequences in temp storage, "
            "keeps only scoring motifs in memory. ~85%% memory savings "
            "(e.g., 11 GB → 2 GB for human), faster with parallel processing.",
        )

        # === Output Options ===
        output = parser.add_argument_group("output options")
        output.add_argument(
            "--clean-names",
            "--clean_names",
            action="store_true",
            default=True,
            help='Remove "transcript:" and "gene:" prefixes (default: True)',
        )
        output.add_argument(
            "--no-clean-names",
            "--no_clean_names",
            dest="clean_names",
            action="store_false",
            help="Keep ID prefixes",
        )
        output.add_argument(
            "-u",
            "--uninformative-naming",
            "--uninformative_naming",
            action="store_true",
            help="Use simple naming scheme",
        )
        output.add_argument(
            "--no-abbreviate",
            "--no_abbreviate",
            "--na",
            action="store_true",
            help="Use full species name in outputs",
        )
        output.add_argument(
            "--abbreviate-filenames",
            "--abbreviate_filenames",
            "--afn",
            action="store_true",
            help="Abbreviate species name in filenames",
        )
        output.add_argument(
            "--no-headers",
            "--no_headers",
            action="store_true",
            help="Omit column headers from output files (default: include headers)",
        )

        # === Advanced ===
        advanced = parser.add_argument_group("advanced options")
        advanced.add_argument(
            "--seed", type=int, default=42, help="Random seed (default: 42)"
        )

    def _require_file_exists(self, path: Optional[Path], command: str, label: str):
        """Validate file exists, error if not."""
        if path and not path.exists():
            self.parser.error(f"{command}: {label} not found: {path}")

    def _validate_files(self, args, command: str, file_specs: list):
        """
        Validate multiple file attributes exist.

        Args:
            args: Parsed arguments
            command: Command name for error messages
            file_specs: List of (attr_name, label) tuples
        """
        for attr, label in file_specs:
            path = getattr(args, attr, None)
            self._require_file_exists(path, command, label)

    def _validate_args(self, args):
        """Validate parsed arguments based on command.

        Args:
            args: Parsed argument namespace

        Raises:
            argparse.ArgumentTypeError: If validation fails
        """
        # Skip validation if generating config
        if getattr(args, "generate_config", False):
            return

        # Get command (train or classify)
        command = getattr(args, "command", "classify")

        # ===============================================================
        # TRAIN MODE VALIDATION
        # ===============================================================
        if command == "train":
            # Species name required
            if not args.species_name:
                self.parser.error("train: -n/--species_name is required")

            # Validate CV parameters
            if args.n_cv_folds < 2:
                self.parser.error("train: n_cv_folds must be >= 2")

            if not 0 < args.test_fraction < 1:
                self.parser.error("train: test_fraction must be between 0 and 1")

            # Validate custom reference files exist
            self._validate_files(
                args,
                "train",
                [
                    ("reference_u12s", "U12 reference file"),
                    ("reference_u2s", "U2 reference file"),
                ],
            )

        # ===============================================================
        # EXTRACT MODE VALIDATION
        # ===============================================================
        elif command == "extract":
            # Species name required
            if not args.species_name:
                self.parser.error("extract: -n/--species_name is required")

            # Genome required
            if not args.genome:
                self.parser.error("extract: genome file (-g) is required")

            # Input validation
            has_annotation = args.annotation is not None
            has_bed = args.bed is not None

            # Must have either annotation or bed
            if not has_annotation and not has_bed:
                self.parser.error(
                    "extract: must specify input:\n"
                    "  Annotation: -a ANNOTATION\n"
                    "  BED: -b BED"
                )

            # Can't have both
            if has_annotation and has_bed:
                self.parser.error(
                    "extract: specify only one: -a (annotation) or -b (BED)"
                )

            # Validate file paths exist
            self._validate_files(
                args,
                "extract",
                [
                    ("genome", "genome file"),
                    ("annotation", "annotation file"),
                    ("bed", "BED file"),
                ],
            )

        # ===============================================================
        # CLASSIFY MODE VALIDATION
        # ===============================================================
        elif command == "classify":
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

            # Model source validation
            has_model = hasattr(args, "model") and args.model is not None

            if not has_model:
                # Try to use default pretrained model if available
                from intronIC.cli.config import get_default_pretrained_model_path

                default_model = get_default_pretrained_model_path()

                if default_model:
                    # Auto-populate with default model
                    args.model = default_model
                    has_model = True
                else:
                    # No default model available - require explicit specification
                    self.parser.error(
                        "classify: must specify model:\n"
                        "  --model PATH  (path to pretrained model)\n"
                        "\n"
                        "To train a new model, use: intronIC train -n species_name\n"
                        "Note: Default pretrained model not found at data/default_pretrained.model.pkl"
                    )

            # Validate file paths exist
            self._validate_files(
                args,
                "classify",
                [
                    ("genome", "genome file"),
                    ("annotation", "annotation file"),
                    ("bed", "BED file"),
                    ("sequence_file", "sequence file"),
                    ("model", "model file"),
                ],
            )

            # Threshold validation
            if not 0 <= args.threshold <= 100:
                self.parser.error("classify: threshold must be between 0 and 100")

            # Process count validation
            if args.processes is not None and args.processes < 1:
                self.parser.error("classify: processes must be >= 1")

        # Create output directory (skip for test command which may have None output_dir)
        if hasattr(args, 'output_dir') and args.output_dir is not None:
            args.output_dir.mkdir(parents=True, exist_ok=True)
