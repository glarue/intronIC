# Next Session Guide: M2.1 - Command-Line Interface (CLI)

**Current Status:** M1.1-M1.5 COMPLETE ✅ (Phase 1 Done!)
**Last Session:** 2025-11-04 (Classification System - M1.5 COMPLETE!)
**Branch:** `refactor/phase-0-foundation`
**Ready to Start:** M2.1 - CLI Implementation
**Tests Passing:** 420/420 ✅ (3 skipped)

---

## 🎯 Session 2025-11-04 Summary

### What We Completed: M1.5 Classification System ✅

**Production Code Added:** ~1,535 lines
- `classification/optimizer.py` (356 lines) - Geometric grid search
- `classification/trainer.py` (278 lines) - Ensemble training with U2 subsampling
- `classification/predictor.py` (181 lines) - F1-weighted classification
- `classification/classifier.py` (360 lines) - High-level orchestrator

**Test Code Added:** ~2,746 lines (81 new tests, all passing)
- `test_optimizer.py` (520 lines, 19 tests)
- `test_trainer.py` (540 lines, 19 tests)
- `test_predictor.py` (547 lines, 20 tests)
- `test_classifier.py` (559 lines, 18 tests)
- `test_classification_pipeline.py` (580 lines, 5 integration tests)

**Critical Achievements:**
- ✅ **100% Classification Accuracy** on synthetic test data
- ✅ **Issue #1 FIXED:** No re-normalization after classification (data leakage prevented)
- ✅ **Real Data Validation:** Tested with 387 U12 + 20,690 U2 reference introns
- ✅ **Reproducible Results:** Deterministic with random seeds
- ✅ **Production Ready:** Comprehensive error handling, batch processing, logging

---

## 🚀 What to Build: M2.1 - Command-Line Interface (CLI)

### Overview

M2.1 implements a **production-ready command-line interface** that orchestrates the complete intronIC pipeline from genome annotation to classified intron output.

**Goal:** Create a user-facing CLI that:
1. Matches original intronIC's interface (for backwards compatibility)
2. Uses **rich library** for beautiful, colored terminal output 🎨
3. Provides clear progress reporting
4. Validates inputs and provides helpful error messages
5. Supports all input modes (annotation + genome, BED file, sequences)

**Key Benefits:**
- Makes the refactored codebase actually usable
- Validates end-to-end integration
- Enables gold standard comparison with original intronIC
- Professional user experience with rich formatting

**Time Estimate:** 10-14 hours

---

## 📦 New Dependencies

### Rich Library for Beautiful CLI Output

We're adding the [`rich`](https://rich.readthedocs.io/) library for enhanced terminal output:

```bash
# Already added to pyproject.toml
# Install with:
pixi install

# Or if using pip:
pip install rich>=10.0
```

**Rich Features We'll Use:**
- **Console:** Colored text, markup, styled output
- **Progress:** Beautiful progress bars with multiple tasks
- **Table:** Display results in formatted tables
- **Panel:** Highlight important information
- **Tree:** Show file structure and hierarchy
- **Syntax:** Code highlighting for error messages
- **Logging:** Rich-formatted log messages

---

## 📋 Implementation Steps

### Step 1: Project Setup (10 minutes)

```bash
cd /home/glarue/code/intronIC/intronIC_refactored

# Create CLI module
mkdir -p cli tests/integration/test_cli

# Create module files
touch cli/__init__.py
touch cli/main.py
touch cli/args.py
touch cli/config.py
touch cli/progress.py
touch cli/validation.py

# Create test files
touch tests/integration/test_cli/__init__.py
touch tests/integration/test_cli/test_args.py
touch tests/integration/test_cli/test_main.py

# Install rich
pixi install
```

---

### Step 2: Argument Parser (2-3 hours)

**Port from:** `intronIC.py:5908-6093` (argparse setup in `main()`)

**File:** `cli/args.py`

**What to Build:**

```python
"""
Argument parsing for intronIC CLI.

Maintains backwards compatibility with original intronIC interface
while adding new features.
"""

import argparse
from pathlib import Path
from typing import Optional
from rich.console import Console

console = Console()


def create_parser() -> argparse.ArgumentParser:
    """
    Create argument parser matching original intronIC interface.

    Categories:
    - Required: species name
    - Input: genome, annotation, BED, sequences
    - Filtering: length, non-canonical, isoforms
    - Scoring: PWMs, coordinates
    - Classification: threshold, fixed C, ensemble size
    - Output: species name, output directory
    - Performance: parallel processes
    """
    parser = argparse.ArgumentParser(
        prog='intronIC',
        description='[bold blue]intronIC[/bold blue] - Classify introns as U2-type or U12-type',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic: annotation + genome
  intronIC -g genome.fa -a annotation.gff3 -n species_name

  # With custom threshold
  intronIC -g genome.fa -a annotation.gff3 -n species_name -t 85

  # Parallel processing
  intronIC -g genome.fa -a annotation.gff3 -n species_name -p 8

  # Extract sequences only (no classification)
  intronIC -g genome.fa -a annotation.gff3 -n species_name -s

For more information: https://github.com/glarue/intronIC
        """
    )

    # Required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-n', '--species_name',
        required=True,
        help='Species name (used for output file naming)'
    )

    # Input mode (one required)
    input_group = parser.add_argument_group('input (one required)')
    input_mode = input_group.add_mutually_exclusive_group(required=True)

    input_mode.add_argument(
        '-a', '--annotation',
        type=Path,
        help='Annotation file (GFF3/GTF format, gzip supported)'
    )

    input_mode.add_argument(
        '-b', '--bed',
        type=Path,
        help='BED file with intron coordinates'
    )

    input_mode.add_argument(
        '-q', '--sequence_file',
        type=Path,
        help='Pre-extracted intron sequences (.iic format)'
    )

    # Genome (required for annotation/BED modes)
    input_group.add_argument(
        '-g', '--genome',
        type=Path,
        help='Genome FASTA file (required with -a or -b)'
    )

    # Filtering options
    filtering = parser.add_argument_group('filtering')
    filtering.add_argument(
        '--min_intron_len',
        type=int,
        default=30,
        help='Minimum intron length (default: 30)'
    )

    filtering.add_argument(
        '-i', '--allow_multiple_isoforms',
        action='store_true',
        help='Include introns from all isoforms (not just longest)'
    )

    filtering.add_argument(
        '-v', '--no_intron_overlap',
        action='store_true',
        help='Exclude overlapping introns from different isoforms'
    )

    filtering.add_argument(
        '--no_nc',
        action='store_true',
        help='Exclude non-canonical introns from scoring'
    )

    filtering.add_argument(
        '-f', '--feature',
        choices=['cds', 'exon'],
        help='Feature type to extract introns from (default: both)'
    )

    # Scoring options
    scoring = parser.add_argument_group('scoring')
    scoring.add_argument(
        '--pwms',
        type=Path,
        help='Custom PWM file (default: use built-in matrices)'
    )

    scoring.add_argument(
        '--five_score_coords',
        nargs=2,
        type=int,
        metavar=('START', 'STOP'),
        default=(-3, 9),
        help="5' splice site scoring coordinates (default: -3 9)"
    )

    scoring.add_argument(
        '--bp_coords',
        nargs=2,
        type=int,
        metavar=('START', 'STOP'),
        default=(-55, -5),
        help='Branch point search region (default: -55 -5)'
    )

    scoring.add_argument(
        '--three_score_coords',
        nargs=2,
        type=int,
        metavar=('START', 'STOP'),
        default=(-6, 4),
        help="3' splice site scoring coordinates (default: -6 4)"
    )

    scoring.add_argument(
        '--ignore_nc_dnts',
        action='store_true',
        help='Ignore non-canonical dinucleotides when scoring'
    )

    # Classification options
    classification = parser.add_argument_group('classification')
    classification.add_argument(
        '-s', '--sequences_only',
        action='store_true',
        help='Extract sequences only (skip scoring and classification)'
    )

    classification.add_argument(
        '-t', '--threshold',
        type=float,
        default=90.0,
        help='U12 classification threshold 0-100 (default: 90)'
    )

    classification.add_argument(
        '-C',
        type=float,
        help='Fixed SVM C parameter (skips optimization)'
    )

    classification.add_argument(
        '--n_models',
        type=int,
        default=3,
        help='Number of ensemble models (default: 3)'
    )

    classification.add_argument(
        '--n_optimization_rounds',
        type=int,
        default=5,
        help='Hyperparameter optimization rounds (default: 5)'
    )

    classification.add_argument(
        '--reference_u12s',
        type=Path,
        help='Custom U12 reference sequences'
    )

    classification.add_argument(
        '--reference_u2s',
        type=Path,
        help='Custom U2 reference sequences'
    )

    # Output options
    output = parser.add_argument_group('output')
    output.add_argument(
        '-o', '--output_dir',
        type=Path,
        default=Path('.'),
        help='Output directory (default: current directory)'
    )

    output.add_argument(
        '--no_plots',
        action='store_true',
        help='Skip generating plots'
    )

    # Performance options
    performance = parser.add_argument_group('performance')
    performance.add_argument(
        '-p', '--processes',
        type=int,
        default=1,
        help='Number of parallel processes (default: 1)'
    )

    # Logging
    logging_group = parser.add_argument_group('logging')
    logging_group.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose output'
    )

    logging_group.add_argument(
        '--quiet',
        action='store_true',
        help='Minimal output (errors only)'
    )

    logging_group.add_argument(
        '--log_file',
        type=Path,
        help='Write log to file'
    )

    return parser


def validate_args(args: argparse.Namespace) -> None:
    """
    Validate argument combinations and file existence.

    Raises:
        ValueError: If arguments are invalid or inconsistent
    """
    # Check genome is provided for annotation/BED modes
    if (args.annotation or args.bed) and not args.genome:
        raise ValueError(
            "Genome file (-g/--genome) is required when using "
            "annotation (-a) or BED file (-b) input"
        )

    # Check input files exist
    if args.annotation and not args.annotation.exists():
        raise FileNotFoundError(f"Annotation file not found: {args.annotation}")

    if args.bed and not args.bed.exists():
        raise FileNotFoundError(f"BED file not found: {args.bed}")

    if args.sequence_file and not args.sequence_file.exists():
        raise FileNotFoundError(f"Sequence file not found: {args.sequence_file}")

    if args.genome and not args.genome.exists():
        raise FileNotFoundError(f"Genome file not found: {args.genome}")

    if args.pwms and not args.pwms.exists():
        raise FileNotFoundError(f"PWM file not found: {args.pwms}")

    # Validate threshold
    if not 0 <= args.threshold <= 100:
        raise ValueError(f"Threshold must be 0-100, got {args.threshold}")

    # Validate coordinates
    if args.five_score_coords[0] >= args.five_score_coords[1]:
        raise ValueError("Invalid five_score_coords: start must be < stop")

    if args.bp_coords[0] >= args.bp_coords[1]:
        raise ValueError("Invalid bp_coords: start must be < stop")

    if args.three_score_coords[0] >= args.three_score_coords[1]:
        raise ValueError("Invalid three_score_coords: start must be < stop")

    # Create output directory if needed
    if not args.output_dir.exists():
        args.output_dir.mkdir(parents=True, exist_ok=True)

    # Validate conflicting options
    if args.quiet and args.verbose:
        raise ValueError("Cannot use both --quiet and --verbose")
```

**Tests:** `tests/integration/test_cli/test_args.py`

---

### Step 3: Rich Progress Reporter (2 hours)

**File:** `cli/progress.py`

**What to Build:**

```python
"""
Rich-formatted progress reporting for intronIC pipeline.

Uses rich library for beautiful terminal output with:
- Colored status messages
- Progress bars for long operations
- Formatted results tables
- Highlighted warnings and errors
"""

from typing import Optional, Dict, Any
from pathlib import Path
from rich.console import Console
from rich.progress import (
    Progress,
    SpinnerColumn,
    TextColumn,
    BarColumn,
    TaskProgressColumn,
    TimeRemainingColumn,
    TimeElapsedColumn,
)
from rich.table import Table
from rich.panel import Panel
from rich.syntax import Syntax
from rich.tree import Tree
from rich import box


class IntronICProgress:
    """Progress reporter with rich formatting."""

    def __init__(self, verbose: bool = False, quiet: bool = False):
        self.console = Console()
        self.verbose = verbose
        self.quiet = quiet

        # Create progress bar
        self.progress = Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
            console=self.console,
            transient=False  # Keep progress bars after completion
        )

    def print_header(self):
        """Print intronIC header."""
        if self.quiet:
            return

        header = """
[bold blue]╔═══════════════════════════════════════════════════╗[/bold blue]
[bold blue]║[/bold blue]   [bold white]intronIC v1.5.1 (Refactored)[/bold white]                 [bold blue]║[/bold blue]
[bold blue]║[/bold blue]   U2/U12 Intron Classification Pipeline       [bold blue]║[/bold blue]
[bold blue]╚═══════════════════════════════════════════════════╝[/bold blue]
        """
        self.console.print(header)

    def print_config(self, config: Dict[str, Any]):
        """Print configuration summary."""
        if self.quiet:
            return

        table = Table(title="Configuration", box=box.ROUNDED, show_header=False)
        table.add_column("Parameter", style="cyan", no_wrap=True)
        table.add_column("Value", style="white")

        for key, value in config.items():
            table.add_row(key, str(value))

        self.console.print(table)
        self.console.print()

    def start_stage(self, stage_name: str, description: str = ""):
        """Start a new pipeline stage."""
        if self.quiet:
            return

        self.console.print()
        self.console.rule(f"[bold cyan]{stage_name}[/bold cyan]")
        if description:
            self.console.print(f"[dim]{description}[/dim]")
        self.console.print()

    def success(self, message: str):
        """Print success message."""
        if not self.quiet:
            self.console.print(f"[bold green]✓[/bold green] {message}")

    def info(self, message: str):
        """Print info message."""
        if self.verbose:
            self.console.print(f"[blue]ℹ[/blue] {message}")

    def warning(self, message: str):
        """Print warning message."""
        if not self.quiet:
            self.console.print(f"[yellow]⚠[/yellow] {message}")

    def error(self, message: str):
        """Print error message."""
        self.console.print(f"[bold red]✗[/bold red] {message}")

    def print_stats(self, title: str, stats: Dict[str, Any]):
        """Print statistics table."""
        if self.quiet:
            return

        table = Table(title=title, box=box.SIMPLE)
        table.add_column("Metric", style="cyan")
        table.add_column("Count", justify="right", style="white")
        table.add_column("Percentage", justify="right", style="dim")

        total = stats.get('total', 0)
        for key, value in stats.items():
            if key == 'total':
                continue
            percentage = f"{100 * value / total:.1f}%" if total > 0 else "N/A"
            table.add_row(key.replace('_', ' ').title(), str(value), percentage)

        if 'total' in stats:
            table.add_section()
            table.add_row("[bold]Total[/bold]", f"[bold]{stats['total']}[/bold]", "100.0%")

        self.console.print(table)
        self.console.print()

    def print_results(self, results: Dict[str, Any]):
        """Print classification results."""
        if self.quiet:
            return

        panel = Panel.fit(
            f"""[bold green]Classification Complete![/bold green]

[cyan]U12-type introns:[/cyan] [white]{results.get('u12_count', 0)}[/white] ([dim]{results.get('u12_percentage', 0):.2f}%[/dim])
[cyan]U2-type introns:[/cyan]  [white]{results.get('u2_count', 0)}[/white] ([dim]{results.get('u2_percentage', 0):.2f}%[/dim])
[cyan]Total classified:[/cyan] [white]{results.get('total', 0)}[/white]

[dim]Threshold:[/dim] [yellow]{results.get('threshold', 90.0)}%[/yellow]
[dim]Mean SVM score:[/dim] [yellow]{results.get('mean_svm_score', 0):.2f}[/yellow]
            """,
            title="Results",
            border_style="green"
        )
        self.console.print(panel)

    def print_files(self, output_dir: Path, species_name: str):
        """Print output files in tree structure."""
        if self.quiet:
            return

        tree = Tree(
            f"[bold cyan]Output Files[/bold cyan] ({output_dir})",
            guide_style="dim"
        )

        files = [
            (f"{species_name}.meta.iic", "Metadata with all intron information"),
            (f"{species_name}.bed.iic", "BED format for genome browsers"),
            (f"{species_name}.seqs.iic", "Intron sequences with flanks"),
            (f"{species_name}.scores.iic", "Detailed scoring information"),
            (f"{species_name}.dupe_map.iic", "Duplicate intron mappings"),
            (f"{species_name}.overlap_map.iic", "Overlapping intron mappings"),
            (f"{species_name}.log", "Execution log"),
        ]

        for filename, description in files:
            filepath = output_dir / filename
            if filepath.exists():
                size = filepath.stat().st_size
                size_str = self._format_size(size)
                tree.add(f"[green]{filename}[/green] [dim]({size_str})[/dim] - {description}")

        self.console.print(tree)
        self.console.print()

    def _format_size(self, size_bytes: int) -> str:
        """Format file size in human-readable format."""
        for unit in ['B', 'KB', 'MB', 'GB']:
            if size_bytes < 1024:
                return f"{size_bytes:.1f} {unit}"
            size_bytes /= 1024
        return f"{size_bytes:.1f} TB"

    def create_progress_bar(self, description: str, total: Optional[int] = None):
        """Create a task in the progress bar."""
        return self.progress.add_task(description, total=total)

    def update_progress(self, task_id, advance: int = 1, **kwargs):
        """Update progress bar."""
        self.progress.update(task_id, advance=advance, **kwargs)

    def __enter__(self):
        """Enter context (start progress bars)."""
        if not self.quiet:
            self.progress.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Exit context (stop progress bars)."""
        if not self.quiet:
            self.progress.stop()
```

---

### Step 4: Main CLI Entry Point (4-5 hours)

**File:** `cli/main.py`

**What to Build:**

```python
"""
Main CLI entry point for intronIC.

Orchestrates the complete pipeline:
1. Parse arguments
2. Load inputs (genome, annotation, etc.)
3. Extract introns
4. Filter and deduplicate
5. Score with PWMs
6. Normalize z-scores
7. Classify with SVM ensemble
8. Write outputs
9. Generate plots (optional)
"""

import sys
from pathlib import Path
from typing import Optional
import logging

from rich.console import Console
from rich.logging import RichHandler

from cli.args import create_parser, validate_args
from cli.progress import IntronICProgress
from cli.config import Config

# Import pipeline components
from io_operations.parsers import GenomeReader, AnnotationParser, BEDParser
from extraction.extractor import IntronExtractor
from extraction.filter import IntronFilter
from scoring.scorer import IntronScorer
from scoring.normalizer import ScoreNormalizer
from classification.classifier import IntronClassifier


def setup_logging(args) -> logging.Logger:
    """Setup logging with rich formatting."""
    # Create logger
    logger = logging.getLogger("intronIC")
    logger.setLevel(logging.DEBUG if args.verbose else logging.INFO)

    # Console handler with rich formatting
    if not args.quiet:
        console_handler = RichHandler(
            rich_tracebacks=True,
            markup=True,
            show_time=True,
            show_path=args.verbose
        )
        console_handler.setLevel(logging.DEBUG if args.verbose else logging.INFO)
        logger.addHandler(console_handler)

    # File handler if requested
    if args.log_file:
        file_handler = logging.FileHandler(args.log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(
            logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        )
        logger.addHandler(file_handler)

    return logger


def main():
    """Main CLI entry point."""
    console = Console()

    try:
        # Parse arguments
        parser = create_parser()
        args = parser.parse_args()

        # Validate arguments
        try:
            validate_args(args)
        except (ValueError, FileNotFoundError) as e:
            console.print(f"[bold red]Error:[/bold red] {e}")
            sys.exit(1)

        # Setup logging
        logger = setup_logging(args)

        # Create progress reporter
        with IntronICProgress(verbose=args.verbose, quiet=args.quiet) as progress:
            progress.print_header()

            # Build configuration
            config = Config.from_args(args)
            progress.print_config(config.to_display_dict())

            # Run pipeline
            result = run_pipeline(args, config, progress, logger)

            # Print results
            progress.print_results(result)
            progress.print_files(args.output_dir, args.species_name)

            progress.success("Pipeline complete!")

        return 0

    except KeyboardInterrupt:
        console.print("\n[yellow]Interrupted by user[/yellow]")
        return 130

    except Exception as e:
        console.print(f"\n[bold red]Fatal error:[/bold red] {e}")
        if args.verbose if 'args' in locals() else False:
            console.print_exception()
        return 1


def run_pipeline(args, config, progress, logger):
    """
    Execute the complete intronIC pipeline.

    Returns:
        Dict with results and statistics
    """
    # Stage 1: Load inputs
    progress.start_stage("Stage 1: Loading Inputs", "Reading genome and annotation files")

    # ... implementation continues

    # See full implementation in Step 5


if __name__ == '__main__':
    sys.exit(main())
```

---

### Step 5: Pipeline Orchestration (3-4 hours)

Complete the `run_pipeline()` function that ties all components together:

1. **Load inputs** (genome, annotation, PWMs, references)
2. **Extract introns** (from annotation/BED/sequences)
3. **Filter and deduplicate**
4. **Score with PWMs**
5. **Normalize z-scores**
6. **Train classifier** (or load pre-trained)
7. **Classify introns**
8. **Write outputs**
9. **Generate plots** (optional)

Each stage should:
- Use progress reporter for visual feedback
- Log important information
- Handle errors gracefully
- Report statistics

---

### Step 6: Configuration Management (1-2 hours)

**File:** `cli/config.py`

Create a `Config` dataclass that holds all pipeline parameters:

```python
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Tuple, Dict, Any

@dataclass
class Config:
    """Pipeline configuration."""

    # Input
    species_name: str
    genome: Optional[Path] = None
    annotation: Optional[Path] = None
    bed: Optional[Path] = None
    sequence_file: Optional[Path] = None

    # Filtering
    min_intron_len: int = 30
    allow_multiple_isoforms: bool = False
    no_intron_overlap: bool = False
    no_nc: bool = False
    feature: Optional[str] = None

    # Scoring
    pwms: Optional[Path] = None
    five_score_coords: Tuple[int, int] = (-3, 9)
    bp_coords: Tuple[int, int] = (-55, -5)
    three_score_coords: Tuple[int, int] = (-6, 4)
    ignore_nc_dnts: bool = False

    # Classification
    sequences_only: bool = False
    threshold: float = 90.0
    fixed_c: Optional[float] = None
    n_models: int = 3
    n_optimization_rounds: int = 5
    reference_u12s: Optional[Path] = None
    reference_u2s: Optional[Path] = None

    # Output
    output_dir: Path = Path('.')
    no_plots: bool = False

    # Performance
    processes: int = 1

    @classmethod
    def from_args(cls, args):
        """Create Config from argparse namespace."""
        return cls(
            species_name=args.species_name,
            genome=args.genome,
            annotation=args.annotation,
            bed=args.bed,
            sequence_file=args.sequence_file,
            min_intron_len=args.min_intron_len,
            allow_multiple_isoforms=args.allow_multiple_isoforms,
            no_intron_overlap=args.no_intron_overlap,
            no_nc=args.no_nc,
            feature=args.feature,
            pwms=args.pwms,
            five_score_coords=tuple(args.five_score_coords),
            bp_coords=tuple(args.bp_coords),
            three_score_coords=tuple(args.three_score_coords),
            ignore_nc_dnts=args.ignore_nc_dnts,
            sequences_only=args.sequences_only,
            threshold=args.threshold,
            fixed_c=args.C,
            n_models=args.n_models,
            n_optimization_rounds=args.n_optimization_rounds,
            reference_u12s=args.reference_u12s,
            reference_u2s=args.reference_u2s,
            output_dir=args.output_dir,
            no_plots=args.no_plots,
            processes=args.processes,
        )

    def to_display_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for display."""
        # Format for rich table display
        return {
            "Species": self.species_name,
            "Input Mode": self._get_input_mode(),
            "Min Intron Length": f"{self.min_intron_len} bp",
            "Classification Threshold": f"{self.threshold}%",
            "Ensemble Models": self.n_models,
            "Parallel Processes": self.processes,
            "Output Directory": str(self.output_dir),
        }

    def _get_input_mode(self) -> str:
        """Determine input mode."""
        if self.annotation:
            return f"Annotation ({self.annotation.name})"
        elif self.bed:
            return f"BED ({self.bed.name})"
        elif self.sequence_file:
            return f"Sequences ({self.sequence_file.name})"
        return "Unknown"
```

---

### Step 7: Integration Tests (2-3 hours)

**File:** `tests/integration/test_cli/test_main.py`

Write integration tests that:
1. Test CLI with chr19 test data
2. Verify all output files are created
3. Check output file formats
4. Validate classification results
5. Test error handling

---

### Step 8: Update Entry Point (15 minutes)

Update `pyproject.toml` to use new CLI:

```toml
[project.scripts]
intronIC = "cli.main:main"
```

---

## 📊 Rich Output Examples

Here's what the user will see:

```
╔═══════════════════════════════════════════════════╗
║   intronIC v1.5.1 (Refactored)                    ║
║   U2/U12 Intron Classification Pipeline           ║
╚═══════════════════════════════════════════════════╝

┌─ Configuration ─────────────────────────────────┐
│ Species                  │ homo_sapiens          │
│ Input Mode               │ Annotation (chr19....)│
│ Min Intron Length        │ 30 bp                 │
│ Classification Threshold │ 90.0%                 │
│ Ensemble Models          │ 3                     │
│ Parallel Processes       │ 1                     │
│ Output Directory         │ .                     │
└─────────────────────────────────────────────────┘

━━━━━━━━━━━━━━━━━━━━━ Stage 1: Loading Inputs ━━━━━━━━━━━━━━━━━━━━━
Reading genome and annotation files

⠹ Loading genome...      ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:02
✓ Loaded 1 chromosome, 58.6 MB

⠹ Parsing annotation...  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 0:00:03
✓ Parsed 15,234 genes, 45,678 transcripts

━━━━━━━━━━━━━━━━━━ Stage 2: Extracting Introns ━━━━━━━━━━━━━━━━━━━

...

┌─ Classification Complete! ──────────────────────┐
│                                                  │
│ U12-type introns: 42 (0.32%)                    │
│ U2-type introns:  13,158 (99.68%)               │
│ Total classified: 13,200                        │
│                                                  │
│ Threshold: 90.0%                                 │
│ Mean SVM score: 15.3                             │
└──────────────────────────────────────────────────┘

Output Files (.)
├── homo_sapiens.meta.iic (1.2 MB) - Metadata...
├── homo_sapiens.bed.iic (512 KB) - BED format...
├── homo_sapiens.seqs.iic (8.3 MB) - Sequences...
└── homo_sapiens.log (45 KB) - Execution log

✓ Pipeline complete!
```

---

## ✅ Success Criteria

### M2.1 Complete When:
- [ ] CLI accepts all original intronIC arguments
- [ ] Rich formatting displays correctly
- [ ] Progress bars show for long operations
- [ ] All pipeline stages execute successfully
- [ ] Output files match expected formats
- [ ] Error messages are clear and helpful
- [ ] Integration tests pass
- [ ] Can run on chr19 test data
- [ ] Performance is comparable to original

---

## 🎨 Rich Library Features to Use

### Console Styling
```python
from rich.console import Console
console = Console()

console.print("[bold blue]intronIC[/bold blue]")
console.print("[green]✓[/green] Success!")
console.print("[yellow]⚠[/yellow] Warning")
console.print("[red]✗[/red] Error")
```

### Progress Bars
```python
from rich.progress import track

for item in track(items, description="Processing"):
    process(item)
```

### Tables
```python
from rich.table import Table

table = Table(title="Results")
table.add_column("Type", style="cyan")
table.add_column("Count", justify="right")
table.add_row("U12", "42")
table.add_row("U2", "13,158")
console.print(table)
```

### Panels
```python
from rich.panel import Panel

console.print(Panel("Important message", title="Note"))
```

---

## 📚 Reference Documentation

- **Rich Library:** https://rich.readthedocs.io/
- **Original intronIC CLI:** `intronIC.py:5908-6093`
- **Argparse Guide:** https://docs.python.org/3/library/argparse.html

---

## 🚀 Next Steps After M2.1

Once the CLI is complete, we can:

1. **M2.2 - Gold Standard Validation**
   - Run both versions on chr19
   - Compare outputs systematically
   - Document any differences

2. **M2.3 - Performance Optimization**
   - Profile the code
   - Optimize bottlenecks
   - Implement parallel processing

3. **M2.4 - Documentation**
   - User guide
   - API documentation
   - Migration guide

---

## 💡 Implementation Tips

1. **Start with basic functionality** - Get a simple pipeline working first
2. **Add rich formatting gradually** - Start with plain output, enhance iteratively
3. **Test frequently** - Run on chr19 data after each stage
4. **Match original behavior** - Compare outputs with original intronIC
5. **Handle errors gracefully** - Rich formatting should work even when things fail

---

## ⏰ Estimated Timeline

- **Day 1 (4 hours):** Argument parsing + validation + config
- **Day 2 (4 hours):** Rich progress reporter + basic pipeline structure
- **Day 3 (4 hours):** Pipeline orchestration (stages 1-5)
- **Day 4 (3 hours):** Pipeline completion (stages 6-9) + testing

**Total: 12-14 hours**

---

**Ready to build a beautiful CLI! 🎨✨**
