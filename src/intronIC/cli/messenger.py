"""
Unified messaging system for intronIC CLI.

Provides synchronized output to both console (with rich formatting) and log file
(with preserved ANSI colors) to eliminate redundant message calls.
"""

import logging
from datetime import datetime
from typing import List, Optional

from rich import box
from rich.console import Console
from rich.progress import Progress
from rich.table import Table
from rich.text import Text

from intronIC.cli.colors import PALETTE


class UnifiedMessenger:
    """
    Unified messaging that sends to both console and log with Rich formatting.

    Eliminates redundancy by providing single method calls that output to both
    destinations with preserved ANSI colors and formatting in the log file.
    """

    def __init__(
        self,
        console: Console,
        log_console: Console,
        logger: logging.Logger,
        quiet: bool = False,
        show_timestamps: bool = True,
    ):
        """
        Initialize unified messenger.

        Args:
            console: Rich console for terminal output
            log_console: Rich console for log file output (with ANSI)
            logger: Python logger (for log levels, but we use consoles for output)
            quiet: If True, suppress console output (but not logs)
            show_timestamps: If True, show timestamps on key operations
        """
        self.console = console
        self.log_console = log_console
        self.logger = logger
        self.quiet = quiet
        self.show_timestamps = show_timestamps
        self.last_timestamp = None  # Track last timestamp to avoid spam

    def _get_timestamp(self, force: bool = False) -> str:
        """
        Get formatted timestamp if enough time has passed.

        Args:
            force: If True, always return timestamp

        Returns:
            Formatted timestamp string or empty string
        """
        if not self.show_timestamps and not force:
            return ""

        now = datetime.now()

        # Show timestamp if forced or if >1 second since last timestamp
        if (
            force
            or self.last_timestamp is None
            or (now - self.last_timestamp).total_seconds() >= 1.0
        ):
            self.last_timestamp = now
            return f"[{PALETTE.timestamp}][{now.strftime('%H:%M:%S')}][/{PALETTE.timestamp}] "

        return ""

    def info(self, message: str, with_timestamp: bool = False):
        """
        Send info message to both console (with ℹ) and log (with preserved formatting).

        Args:
            message: Message text
            with_timestamp: If True, force timestamp on this message
        """
        timestamp = self._get_timestamp(force=with_timestamp)
        # Apply color to the entire message text for visibility in logs
        formatted = f"{timestamp}[{PALETTE.info}]ℹ {message}[/{PALETTE.info}]"

        # Console: with formatting
        if not self.quiet:
            self.console.print(formatted)

        # Log: same formatting with ANSI preserved, written directly to log_console
        self.log_console.print(formatted)

    def success(self, message: str, with_timestamp: bool = False):
        """
        Send success message to both console (with ✓) and log (with preserved formatting).

        Args:
            message: Message text
            with_timestamp: If True, force timestamp on this message
        """
        timestamp = self._get_timestamp(force=with_timestamp)
        # Apply color to the entire message text for visibility in logs
        formatted = f"{timestamp}[{PALETTE.success}]✓ {message}[/{PALETTE.success}]"

        # Console: with formatting
        if not self.quiet:
            self.console.print(formatted)

        # Log: same formatting with ANSI preserved, written directly to log_console
        self.log_console.print(formatted)

    def warning(self, message: str, with_timestamp: bool = False):
        """
        Send warning message to both console (with ⚠) and log (with preserved formatting).

        Args:
            message: Message text
            with_timestamp: If True, force timestamp on this message
        """
        timestamp = self._get_timestamp(force=with_timestamp)
        # Apply color to the entire message text for visibility in logs
        formatted = f"{timestamp}[{PALETTE.warning}]⚠ {message}[/{PALETTE.warning}]"

        # Console: with formatting
        if not self.quiet:
            self.console.print(formatted)

        # Log: same formatting with ANSI preserved, written directly to log_console
        self.log_console.print(formatted)

    def error(self, message: str, with_timestamp: bool = False):
        """
        Send error message to both console (with ✗) and log (with preserved formatting).

        Args:
            message: Message text
            with_timestamp: If True, force timestamp on this message
        """
        timestamp = self._get_timestamp(force=with_timestamp)
        # Apply color and bold to the entire message text for visibility in logs
        formatted = (
            f"{timestamp}[bold {PALETTE.error}]✗ {message}[/bold {PALETTE.error}]"
        )

        # Console: with formatting, bold
        if not self.quiet:
            self.console.print(formatted)

        # Log: same formatting with ANSI preserved, written directly to log_console
        self.log_console.print(formatted)

    def step(
        self,
        step_num: int,
        step_name: str,
        all_steps: List[str],
        show_tree: bool = True,
    ):
        """
        Display pipeline step header with optional tree view.

        Args:
            step_num: Step number (1-based)
            step_name: Name of current step
            all_steps: List of all pipeline steps
            show_tree: Whether to show the pipeline tree (default: True)
        """
        # Force timestamp for step transitions
        timestamp = self._get_timestamp(force=True)

        # Console: Show section divider + optional tree
        if not self.quiet:
            # Section divider
            self.console.print()
            self.console.rule(
                f"[{PALETTE.header}]Step {step_num}: {step_name}", style=PALETTE.header
            )

            # Optional: Show pipeline tree with current step highlighted
            if show_tree:
                from rich.tree import Tree

                tree = Tree(
                    f"🔄 [{PALETTE.highlight}]Pipeline Steps", guide_style="dim"
                )

                for i, step in enumerate(all_steps, 1):
                    if i < step_num:
                        # Completed
                        tree.add(f"[{PALETTE.step_complete}]✓ {step}")
                    elif i == step_num:
                        # Current
                        tree.add(f"[{PALETTE.step_current}]▶ {step}")
                    else:
                        # Pending
                        tree.add(f"[{PALETTE.step_pending}]○ {step}")

                self.console.print(tree)

        # Log: Simple section marker with color
        sep_line = "=" * 80
        step_text = f"STEP {step_num}: {step_name.upper()}"
        self.log_console.print("")  # Blank line
        self.log_console.print(f"[{PALETTE.header}]{sep_line}[/{PALETTE.header}]")
        self.log_console.print(
            f"[bold {PALETTE.header}]{step_text}[/bold {PALETTE.header}]"
        )
        self.log_console.print(f"[{PALETTE.header}]{sep_line}[/{PALETTE.header}]")
        self.log_console.print(
            f"{timestamp}[{PALETTE.info}]Starting step {step_num} of {len(all_steps)}[/{PALETTE.info}]"
        )

    def log_only(self, message: str, level: str = "info"):
        """
        Send message only to log file (not console).

        Useful for detailed debugging info that shouldn't clutter console.

        Args:
            message: Message text
            level: Log level (info, warning, error, debug) - currently ignored, all written as plain text
        """
        # Write directly to log_console without Rich formatting
        self.log_console.print(message)

    def console_only(self, message: str, style: str = ""):
        """
        Send message only to console (not log).

        Useful for interactive elements that don't belong in logs.

        Args:
            message: Message text
            style: Rich style string
        """
        if not self.quiet:
            if style:
                self.console.print(message, style=style)
            else:
                self.console.print(message)

    def print_startup_banner(
        self,
        species_name: str,
        input_mode: str,
        output_dir: str,
        threshold: float,
        command: Optional[str] = None,
        working_dir: Optional[str] = None,
        model_path: Optional[str] = None,
        genome_path: Optional[str] = None,
        annotation_path: Optional[str] = None,
        bed_path: Optional[str] = None,
        sequences_path: Optional[str] = None,
    ):
        """
        Print startup banner with key metadata to both console and log.

        Args:
            species_name: Name of species being analyzed
            input_mode: Input mode (annotation, bed, sequences)
            output_dir: Output directory path
            threshold: Classification threshold
            command: Full command line (optional)
            working_dir: Working directory (optional)
            model_path: Path to model file (optional)
            genome_path: Path to genome file (optional)
            annotation_path: Path to annotation file (optional)
            bed_path: Path to BED file (optional)
            sequences_path: Path to sequences file (optional)
        """
        from datetime import datetime

        from rich import box
        from rich.panel import Panel
        from rich.text import Text

        # Get version
        try:
            from importlib.metadata import version

            ver = version("intronIC")
        except Exception:
            ver = "2.0.0"

        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Console: Rich panel with metadata
        if not self.quiet:
            header = Text()
            header.append("intronIC ", style=f"bold {PALETTE.highlight}")
            header.append(f"v{ver}\n\n", style="dim")
            header.append("Run name: ", style=PALETTE.table_value)
            header.append(f"{species_name}\n", style=PALETTE.mustard)
            header.append("Input: ", style=PALETTE.table_value)
            header.append(f"{input_mode}\n", style=PALETTE.success)
            header.append("Output: ", style=PALETTE.table_value)
            header.append(f"{output_dir}\n", style=PALETTE.path)
            header.append("Threshold: ", style=PALETTE.table_value)
            header.append(f"{threshold}%\n", style=PALETTE.highlight)

            if model_path:
                # Show just the filename for brevity
                import os

                model_name = os.path.basename(model_path)
                header.append("Model: ", style=PALETTE.table_value)
                header.append(f"{model_name}\n", style=PALETTE.path)

            header.append("Started: ", style=PALETTE.table_value)
            header.append(f"{timestamp}", style=PALETTE.timestamp)

            panel = Panel(
                header, border_style=PALETTE.highlight, box=box.DOUBLE, padding=(1, 2)
            )
            self.console.print(panel)

        # Log: Clean text format with full paths and color
        sep_line = "=" * 80
        self.log_console.print(f"[{PALETTE.header}]{sep_line}[/{PALETTE.header}]")
        self.log_console.print(
            f"[bold {PALETTE.highlight}]intronIC v{ver}[/bold {PALETTE.highlight}]"
        )
        self.log_console.print(
            f"[{PALETTE.timestamp}]Started: {timestamp}[/{PALETTE.timestamp}]"
        )
        self.log_console.print(f"[{PALETTE.header}]{sep_line}[/{PALETTE.header}]")
        self.log_console.print("")

        # Combined Command and Configuration section
        self.log_console.print(
            f"[bold {PALETTE.header}]Command and Configuration:[/bold {PALETTE.header}]"
        )

        # Command and working directory
        if command:
            self.log_console.print(
                f"  [{PALETTE.table_header}]Command:[/{PALETTE.table_header}] {command}"
            )
        if working_dir:
            self.log_console.print(
                f"  [{PALETTE.table_header}]Working directory:[/{PALETTE.table_header}] [{PALETTE.path}]{working_dir}[/{PALETTE.path}]"
            )

        # Configuration details
        self.log_console.print(
            f"  [{PALETTE.table_header}]Run name:[/{PALETTE.table_header}] {species_name}"
        )
        self.log_console.print(
            f"  [{PALETTE.table_header}]Input mode:[/{PALETTE.table_header}] {input_mode}"
        )
        self.log_console.print(
            f"  [{PALETTE.table_header}]Classification threshold:[/{PALETTE.table_header}] [{PALETTE.highlight}]{threshold}%[/{PALETTE.highlight}]"
        )
        self.log_console.print(
            f"  [{PALETTE.table_header}]Output directory:[/{PALETTE.table_header}] [{PALETTE.path}]{output_dir}[/{PALETTE.path}]"
        )

        # Log input files with full paths
        if genome_path:
            self.log_console.print(
                f"  [{PALETTE.table_header}]Genome:[/{PALETTE.table_header}] [{PALETTE.path}]{genome_path}[/{PALETTE.path}]"
            )
        if annotation_path:
            self.log_console.print(
                f"  [{PALETTE.table_header}]Annotation:[/{PALETTE.table_header}] [{PALETTE.path}]{annotation_path}[/{PALETTE.path}]"
            )
        if bed_path:
            self.log_console.print(
                f"  [{PALETTE.table_header}]BED file:[/{PALETTE.table_header}] [{PALETTE.path}]{bed_path}[/{PALETTE.path}]"
            )
        if sequences_path:
            self.log_console.print(
                f"  [{PALETTE.table_header}]Sequences:[/{PALETTE.table_header}] [{PALETTE.path}]{sequences_path}[/{PALETTE.path}]"
            )
        if model_path:
            self.log_console.print(
                f"  [{PALETTE.table_header}]Model:[/{PALETTE.table_header}] [{PALETTE.path}]{model_path}[/{PALETTE.path}]"
            )

        self.log_console.print("")

    def print_filtering_summary(
        self,
        total: int,
        omitted_short: int = 0,
        omitted_ambiguous: int = 0,
        omitted_noncanonical: int = 0,
        omitted_isoform: int = 0,
        omitted_overlap: int = 0,
        duplicates: int = 0,
        duplicates_only: int = 0,
        total_removed: int = 0,
        kept: int = 0,
    ):
        """
        Print filtering summary table to both console and log.

        Args:
            total: Total number of introns before filtering
            omitted_short: Number omitted for being too short
            omitted_ambiguous: Number omitted for ambiguous bases
            omitted_noncanonical: Number omitted for non-canonical splice sites
            omitted_isoform: Number omitted for being non-longest isoform
            omitted_overlap: Number omitted for overlapping
            duplicates: Total number of duplicates
            duplicates_only: Number of duplicates not also omitted
            total_removed: Total removed (omitted + duplicates_only)
            kept: Number retained for scoring
        """
        from rich import box
        from rich.table import Table

        total_omitted = (
            omitted_short
            + omitted_ambiguous
            + omitted_noncanonical
            + omitted_isoform
            + omitted_overlap
        )

        # Console: Rich table
        if not self.quiet:
            table = Table(
                title="Intron Filtering Summary",
                box=box.ROUNDED,
                title_style=f"bold {PALETTE.highlight}",
            )
            table.add_column(
                "Category", style=PALETTE.table_header, no_wrap=False, width=30
            )
            table.add_column(
                "Count", style=PALETTE.table_value, justify="right", width=12
            )
            table.add_column(
                "Percent", style=PALETTE.table_value, justify="right", width=10
            )

            # Show omitted section if any were omitted
            if total_omitted > 0:
                table.add_row(
                    "Omitted",
                    f"{total_omitted:,}",
                    f"{(total_omitted / total * 100) if total > 0 else 0:.2f}%",
                )
                # Show breakdown items that have non-zero counts
                if omitted_short > 0:
                    table.add_row(
                        "  - Too short",
                        f"{omitted_short:,}",
                        f"{(omitted_short / total * 100) if total > 0 else 0:.2f}%",
                    )
                if omitted_ambiguous > 0:
                    table.add_row(
                        "  - Ambiguous bases",
                        f"{omitted_ambiguous:,}",
                        f"{(omitted_ambiguous / total * 100) if total > 0 else 0:.2f}%",
                    )
                if omitted_noncanonical > 0:
                    table.add_row(
                        "  - Non-canonical",
                        f"{omitted_noncanonical:,}",
                        f"{(omitted_noncanonical / total * 100) if total > 0 else 0:.2f}%",
                    )
                if omitted_isoform > 0:
                    table.add_row(
                        "  - Non-longest isoform",
                        f"{omitted_isoform:,}",
                        f"{(omitted_isoform / total * 100) if total > 0 else 0:.2f}%",
                    )
                if omitted_overlap > 0:
                    table.add_row(
                        "  - Overlapping",
                        f"{omitted_overlap:,}",
                        f"{(omitted_overlap / total * 100) if total > 0 else 0:.2f}%",
                    )

            # Show duplicates row if there are duplicates that weren't also omitted
            if duplicates_only > 0:
                table.add_row(
                    "Duplicates only",
                    f"{duplicates_only:,}",
                    f"{(duplicates_only / total * 100) if total > 0 else 0:.2f}%",
                )

            # Separator and totals
            table.add_section()
            table.add_row(
                "Total filtered out",
                f"{total_removed:,}",
                f"{(total_removed / total * 100) if total > 0 else 0:.2f}%",
                style=f"bold {PALETTE.warning}",
            )
            table.add_row(
                "Retained for scoring",
                f"{kept:,}",
                f"{(kept / total * 100) if total > 0 else 0:.2f}%",
                style=f"bold {PALETTE.success}",
            )

            self.console.print()
            self.console.print(table)
            self.console.print()

        # Log: ASCII table with color
        self.log_console.print("")
        self.log_console.print(
            f"[bold {PALETTE.highlight}]Intron Filtering Summary:[/bold {PALETTE.highlight}]"
        )
        self.log_console.print(
            f"[{PALETTE.table_header}]┌──────────────────────────────┬────────────┬───────────┐[/{PALETTE.table_header}]"
        )
        self.log_console.print(
            f"[{PALETTE.table_header}]│ Category                     │ Count      │ Percent   │[/{PALETTE.table_header}]"
        )
        self.log_console.print(
            f"[{PALETTE.table_header}]├──────────────────────────────┼────────────┼───────────┤[/{PALETTE.table_header}]"
        )

        # Show omitted section if any were omitted
        if total_omitted > 0:
            self.log_console.print(
                f"│ Omitted                      │ {total_omitted:>10,} │ {(total_omitted / total * 100) if total > 0 else 0:>8.2f}% │"
            )
            # Show breakdown items that have non-zero counts
            if omitted_short > 0:
                self.log_console.print(
                    f"│   - Too short                │ {omitted_short:>10,} │ {(omitted_short / total * 100) if total > 0 else 0:>8.2f}% │"
                )
            if omitted_ambiguous > 0:
                self.log_console.print(
                    f"│   - Ambiguous bases          │ {omitted_ambiguous:>10,} │ {(omitted_ambiguous / total * 100) if total > 0 else 0:>8.2f}% │"
                )
            if omitted_noncanonical > 0:
                self.log_console.print(
                    f"│   - Non-canonical            │ {omitted_noncanonical:>10,} │ {(omitted_noncanonical / total * 100) if total > 0 else 0:>8.2f}% │"
                )
            if omitted_isoform > 0:
                self.log_console.print(
                    f"│   - Non-longest isoform      │ {omitted_isoform:>10,} │ {(omitted_isoform / total * 100) if total > 0 else 0:>8.2f}% │"
                )
            if omitted_overlap > 0:
                self.log_console.print(
                    f"│   - Overlapping              │ {omitted_overlap:>10,} │ {(omitted_overlap / total * 100) if total > 0 else 0:>8.2f}% │"
                )

        # Show duplicates row if there are duplicates that weren't also omitted
        if duplicates_only > 0:
            self.log_console.print(
                f"│ Duplicates only              │ {duplicates_only:>10,} │ {(duplicates_only / total * 100) if total > 0 else 0:>8.2f}% │"
            )

        self.log_console.print(
            f"[{PALETTE.table_header}]├──────────────────────────────┼────────────┼───────────┤[/{PALETTE.table_header}]"
        )
        self.log_console.print(
            f"│ [{PALETTE.warning}]Total filtered out[/{PALETTE.warning}]           │ [{PALETTE.warning}]{total_removed:>10,}[/{PALETTE.warning}] │ [{PALETTE.warning}]{(total_removed / total * 100) if total > 0 else 0:>8.2f}%[/{PALETTE.warning}] │"
        )
        self.log_console.print(
            f"│ [{PALETTE.success}]Retained for scoring[/{PALETTE.success}]         │ [{PALETTE.success}]{kept:>10,}[/{PALETTE.success}] │ [{PALETTE.success}]{(kept / total * 100) if total > 0 else 0:>8.2f}%[/{PALETTE.success}] │"
        )
        self.log_console.print(
            f"[{PALETTE.table_header}]└──────────────────────────────┴────────────┴───────────┘[/{PALETTE.table_header}]"
        )
        self.log_console.print("")

    def print_classification_results(
        self,
        total: int,
        u12_count: int,
        u2_count: int,
        atac_count: int,
        threshold: float,
    ):
        """
        Print classification results table to both console and log.

        Args:
            total: Total number of introns classified
            u12_count: Number of U12-type introns
            u2_count: Number of U2-type introns
            atac_count: Number of AT-AC U12-type introns
            threshold: Classification threshold used
        """
        u12_pct = (u12_count / total * 100) if total > 0 else 0
        u2_pct = (u2_count / total * 100) if total > 0 else 0
        atac_pct = (atac_count / total * 100) if total > 0 else 0

        # Console: Rich table
        if not self.quiet:
            table = Table(
                title=f"Classification Results (threshold: {threshold}%)",
                box=box.DOUBLE_EDGE,
                title_style=f"bold {PALETTE.highlight}",
            )
            table.add_column("Type", style=PALETTE.table_header, no_wrap=True)
            table.add_column("Count", style=PALETTE.table_value, justify="right")
            table.add_column("Percentage", style=PALETTE.table_value, justify="right")

            table.add_row(
                f"[{PALETTE.u12_highlight}]U12-type (total)",
                f"[{PALETTE.u12_highlight}]{u12_count:,}",
                f"[{PALETTE.u12_highlight}]{u12_pct:.2f}%",
            )
            table.add_row(
                f"[{PALETTE.u12_highlight}]U12-type (AT-AC)",
                f"[{PALETTE.u12_highlight}]{atac_count:,}",
                f"[{PALETTE.u12_highlight}]{atac_pct:.2f}%",
            )
            table.add_row("U2-type", f"{u2_count:,}", f"{u2_pct:.2f}%")
            table.add_row("[bold]Total", f"[bold]{total:,}", "[bold]100.00%")

            self.console.print(table)

        # Log: ASCII table with color
        self.log_console.print("")
        self.log_console.print(
            f"[bold {PALETTE.highlight}]Classification Results (threshold: {threshold}%):[/bold {PALETTE.highlight}]"
        )
        self.log_console.print(
            f"[{PALETTE.table_header}]┌──────────────────────┬───────────┬────────────┐[/{PALETTE.table_header}]"
        )
        self.log_console.print(
            f"[{PALETTE.table_header}]│ Type                 │ Count     │ Percentage │[/{PALETTE.table_header}]"
        )
        self.log_console.print(
            f"[{PALETTE.table_header}]├──────────────────────┼───────────┼────────────┤[/{PALETTE.table_header}]"
        )
        self.log_console.print(
            f"│ [{PALETTE.u12_highlight}]U12-type (total)[/{PALETTE.u12_highlight}]     │ [{PALETTE.u12_highlight}]{u12_count:>9,}[/{PALETTE.u12_highlight}] │ [{PALETTE.u12_highlight}]{u12_pct:>9.2f}%[/{PALETTE.u12_highlight}] │"
        )
        self.log_console.print(
            f"│ [{PALETTE.u12_highlight}]U12-type (AT-AC)[/{PALETTE.u12_highlight}]     │ [{PALETTE.u12_highlight}]{atac_count:>9,}[/{PALETTE.u12_highlight}] │ [{PALETTE.u12_highlight}]{atac_pct:>9.2f}%[/{PALETTE.u12_highlight}] │"
        )
        self.log_console.print(
            f"│ U2-type              │ {u2_count:>9,} │ {u2_pct:>9.2f}% │"
        )
        self.log_console.print(
            f"[{PALETTE.table_header}]├──────────────────────┼───────────┼────────────┤[/{PALETTE.table_header}]"
        )
        self.log_console.print(
            f"│ [bold]Total[/bold]                │ [bold]{total:>9,}[/bold] │ [bold]{100.0:>9.2f}%[/bold] │"
        )
        self.log_console.print(
            f"[{PALETTE.table_header}]└──────────────────────┴───────────┴────────────┘[/{PALETTE.table_header}]"
        )
        self.log_console.print("")

    def print_dinucleotide_boundaries(
        self, intron_type: str, boundaries: list[tuple[str, int]], top_n: int = 20
    ):
        """
        Print dinucleotide boundary table to both console and log.

        Args:
            intron_type: Type of introns (e.g., "U12-type", "U2-type")
            boundaries: List of (dinucleotide, count) tuples, sorted by count
            top_n: Number of top boundaries to show
        """
        if not boundaries:
            return

        total = sum(count for _, count in boundaries)
        boundaries_to_show = boundaries[:top_n]

        # Console: Rich table
        if not self.quiet:
            table = Table(
                title=f"Top {top_n} Splice Site Boundaries ({intron_type} introns)",
                box=box.ROUNDED,
                title_style=f"bold {PALETTE.highlight}",
            )
            table.add_column("Rank", style=PALETTE.table_header, justify="right")
            table.add_column("Dinucleotide", style=PALETTE.table_value, justify="right")
            table.add_column("Count", style=PALETTE.table_value, justify="right")
            table.add_column("Percent", style=PALETTE.table_value, justify="right")

            for i, (dnts, count) in enumerate(boundaries_to_show, 1):
                pct = (count / total * 100) if total > 0 else 0
                table.add_row(str(i), dnts, f"{count:,}", f"{pct:.2f}%")

            self.console.print()
            self.console.print(table)

        # Log: ASCII table with color
        self.log_console.print("")
        self.log_console.print(
            f"[bold {PALETTE.highlight}]Top {top_n} Splice Site Boundaries ({intron_type} introns):[/bold {PALETTE.highlight}]"
        )
        self.log_console.print(
            f"[{PALETTE.table_header}]┌──────┬──────────────┬──────────┬───────────┐[/{PALETTE.table_header}]"
        )
        self.log_console.print(
            f"[{PALETTE.table_header}]│ Rank │ Dinucleotide │ Count    │ Percent   │[/{PALETTE.table_header}]"
        )
        self.log_console.print(
            f"[{PALETTE.table_header}]├──────┼──────────────┼──────────┼───────────┤[/{PALETTE.table_header}]"
        )

        for i, (dnts, count) in enumerate(boundaries_to_show, 1):
            pct = (count / total * 100) if total > 0 else 0
            self.log_console.print(
                f"│ {i:>4} │ {dnts:>12} │ {count:>8,} │ {pct:>8.2f}% │"
            )

        self.log_console.print(
            f"[{PALETTE.table_header}]└──────┴──────────────┴──────────┴───────────┘[/{PALETTE.table_header}]"
        )
        self.log_console.print("")

    def print_file_tree(self, output_files: dict[str, str]):
        """
        Print output files tree to both console and log.

        Args:
            output_files: Dictionary mapping file type to file path
        """
        # Console: Rich tree
        if not self.quiet:
            from rich.tree import Tree

            tree = Tree(f"📁 [bold {PALETTE.highlight}]Output Files", guide_style="dim")

            for file_type, filepath in output_files.items():
                tree.add(
                    f"[{PALETTE.success}]{file_type}:[/{PALETTE.success}] [{PALETTE.path}]{filepath}"
                )

            self.console.print(tree)

        # Log: ASCII tree with box-drawing characters
        self.log_console.print("")
        self.log_console.print(
            f"[bold {PALETTE.highlight}]Output Files:[/bold {PALETTE.highlight}]"
        )

        # Convert dict to list for easier handling of first/last items
        items = list(output_files.items())
        for i, (file_type, filepath) in enumerate(items):
            is_last = i == len(items) - 1
            prefix = "└──" if is_last else "├──"
            self.log_console.print(
                f"  {prefix} [{PALETTE.success}]{file_type}:[/{PALETTE.success}] [{PALETTE.path}]{filepath}[/{PALETTE.path}]"
            )

        self.log_console.print("")

    def create_progress(self) -> Progress:
        """
        Create and return a rich Progress instance.

        Returns:
            Configured Progress object for progress bars
        """
        from rich.progress import (
            BarColumn,
            Progress,
            SpinnerColumn,
            TaskProgressColumn,
            TextColumn,
            TimeElapsedColumn,
            TimeRemainingColumn,
        )

        return Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
            console=self.console,
            transient=False,
        )

    def print_training_config(
        self,
        species: str,
        threshold: float,
        seed: int,
        max_iter: int,
        u12_count: int,
        u2_count: int,
        n_models: int,
        eval_mode: str,
        fixed_c: float = None,
        n_optimization_rounds: int = None,
    ):
        """Print training configuration summary to log only."""
        # Log: Training configuration section
        self.log_console.print("")
        self.log_console.print(
            f"[bold {PALETTE.header}]Training Configuration:[/bold {PALETTE.header}]"
        )
        self.log_console.print(f"  Species: {species}")
        self.log_console.print(f"  Classification threshold: {threshold}%")
        self.log_console.print(f"  Random seed: {seed}")
        self.log_console.print(f"  Max iterations: {max_iter:,}")
        self.log_console.print(
            f"  Reference data: {u12_count:,} U12-type, {u2_count:,} U2-type"
        )
        self.log_console.print(f"  Ensemble models: {n_models}")
        self.log_console.print(f"  Evaluation mode: {eval_mode}")
        if fixed_c:
            self.log_console.print(f"  C parameter: {fixed_c:.6e} (fixed)")
        else:
            self.log_console.print(f"  C parameter: Optimized via grid search")
            if n_optimization_rounds:
                self.log_console.print(
                    f"  Optimization rounds: {n_optimization_rounds}"
                )
        self.log_console.print("")

    def print_hyperparameter_results(
        self,
        optimized_c: float,
        cv_score: float,
        calibration_method: str,
        dual: bool,
        intercept_scaling: float,
    ):
        """Print hyperparameter optimization results to log only."""
        self.log_console.print("")
        self.log_console.print(
            f"[bold {PALETTE.header}]Hyperparameter Optimization:[/bold {PALETTE.header}]"
        )
        self.log_console.print(f"  Optimized C: {optimized_c:.6e}")
        self.log_console.print(f"  CV score (balanced accuracy): {cv_score:.4f}")
        self.log_console.print(f"  Calibration method: {calibration_method}")
        self.log_console.print(f"  Dual formulation: {dual}")
        self.log_console.print(f"  Intercept scaling: {intercept_scaling}")
        self.log_console.print("")

    def print_nested_cv_results(
        self,
        n_folds: int,
        mean_f1: float,
        std_f1: float,
        mean_pr_auc: float,
        std_pr_auc: float,
        fold_results: list,
    ):
        """Print nested cross-validation results to log only.

        Args:
            n_folds: Number of CV folds
            mean_f1: Mean F1 score across folds
            std_f1: Standard deviation of F1 score
            mean_pr_auc: Mean PR-AUC across folds
            std_pr_auc: Standard deviation of PR-AUC
            fold_results: List of fold result objects with attributes:
                fold_idx, f1_score, pr_auc, n_u12_train, n_u2_train, n_u12_test, n_u2_test
        """
        self.log_console.print("")
        self.log_console.print(
            f"[bold {PALETTE.header}]Nested Cross-Validation Results ({n_folds} folds):[/bold {PALETTE.header}]"
        )
        self.log_console.print(f"  Mean F1 score: {mean_f1:.4f} ± {std_f1:.4f}")
        self.log_console.print(f"  Mean PR-AUC: {mean_pr_auc:.4f} ± {std_pr_auc:.4f}")
        self.log_console.print("")
        self.log_console.print("  Per-fold results:")

        # Group identical results
        from collections import defaultdict

        grouped = defaultdict(list)
        for fold in fold_results:
            # Create a key from the metrics (rounded to avoid floating point issues)
            key = (
                round(fold.f1_score, 4),
                round(fold.pr_auc, 4),
                fold.n_u12_train,
                fold.n_u2_train,
                fold.n_u12_test,
                fold.n_u2_test,
            )
            grouped[key].append(fold.fold_idx + 1)

        # Helper function to format fold indices into ranges
        def format_fold_indices(indices):
            """Format list of fold indices into compact range notation."""
            if len(indices) == 1:
                return str(indices[0])

            sorted_indices = sorted(indices)
            ranges = []
            start = sorted_indices[0]
            end = sorted_indices[0]

            for i in range(1, len(sorted_indices)):
                if sorted_indices[i] == end + 1:
                    end = sorted_indices[i]
                else:
                    # End of a range
                    if start == end:
                        ranges.append(str(start))
                    else:
                        ranges.append(f"{start}-{end}")
                    start = sorted_indices[i]
                    end = sorted_indices[i]

            # Add the last range
            if start == end:
                ranges.append(str(start))
            else:
                ranges.append(f"{start}-{end}")

            return ", ".join(ranges)

        # Create ASCII table with box-drawing characters (wider fold column)
        self.log_console.print(
            f"  ┌───────────────┬──────────┬────────┬───────────┬──────────┬──────────┬─────────┐"
        )
        self.log_console.print(
            f"  │ Fold(s)       │ F1 Score │ PR-AUC │ Train U12 │ Train U2 │ Test U12 │ Test U2 │"
        )
        self.log_console.print(
            f"  ├───────────────┼──────────┼────────┼───────────┼──────────┼──────────┼─────────┤"
        )

        # Sort by first fold index in each group
        for key in sorted(grouped.keys(), key=lambda k: min(grouped[k])):
            f1, pr_auc, n_u12_train, n_u2_train, n_u12_test, n_u2_test = key
            fold_indices = grouped[key]

            # Format fold label
            fold_str = format_fold_indices(fold_indices)
            if len(fold_indices) > 1:
                fold_label = f"{fold_str} (x{len(fold_indices)})"
            else:
                fold_label = f"{fold_str}/{n_folds}"

            self.log_console.print(
                f"  │ {fold_label:13} │ {f1:8.4f} │ {pr_auc:6.4f} │ "
                f"{n_u12_train:>9,} │ {n_u2_train:>8,} │ {n_u12_test:>8,} │ {n_u2_test:>7,} │"
            )

        self.log_console.print(
            f"  └───────────────┴──────────┴────────┴───────────┴──────────┴──────────┴─────────┘"
        )
        self.log_console.print("")

    def print_split_eval_results(
        self,
        n_u12_train: int,
        n_u2_train: int,
        n_u12_test: int,
        n_u2_test: int,
        test_f1: float,
        test_pr_auc: float,
        n_u12_val: int = None,
        n_u2_val: int = None,
    ):
        """Print train/test split evaluation results to log only."""
        self.log_console.print("")
        self.log_console.print(
            f"[bold {PALETTE.header}]Train/Test Split Evaluation:[/bold {PALETTE.header}]"
        )
        self.log_console.print(
            f"  Training set: {n_u12_train + n_u2_train:,} introns "
            f"({n_u2_train:,} U2, {n_u12_train:,} U12)"
        )
        if n_u12_val is not None and n_u2_val is not None:
            self.log_console.print(
                f"  Validation set: {n_u12_val + n_u2_val:,} introns "
                f"({n_u2_val:,} U2, {n_u12_val:,} U12)"
            )
        self.log_console.print(
            f"  Test set: {n_u12_test + n_u2_test:,} introns "
            f"({n_u2_test:,} U2, {n_u12_test:,} U12)"
        )
        self.log_console.print("")
        self.log_console.print(f"  Test set performance:")
        self.log_console.print(f"    F1 score: {test_f1:.4f}")
        self.log_console.print(f"    PR-AUC: {test_pr_auc:.4f}")
        self.log_console.print("")

    def print_ensemble_summary(self, n_models: int, models: list):
        """Print trained ensemble summary to log only.

        Args:
            n_models: Number of models in ensemble
            models: List of model objects with attributes:
                train_size, u12_count, u2_count, parameters.C
        """
        self.log_console.print("")
        self.log_console.print(
            f"[bold {PALETTE.header}]Trained Ensemble ({n_models} models):[/bold {PALETTE.header}]"
        )

        for i, model in enumerate(models, 1):
            self.log_console.print(f"  Model {i}/{n_models}:")
            self.log_console.print(
                f"    Training samples: {model.train_size:,} ({model.u12_count:,} U12, {model.u2_count:,} U2)"
            )
            self.log_console.print(f"    C parameter: {model.parameters.C:.6e}")

        self.log_console.print("")

        self.log_console.print("")
