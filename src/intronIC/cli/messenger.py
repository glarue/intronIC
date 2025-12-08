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
    destinations. Uses dual color scheme: hex colors for console (truecolor),
    named ANSI colors for log files (better compatibility).
    """

    # Message type configurations: emoji, palette color key, bold flag
    _MESSAGE_TYPES = {
        "info": ("ℹ", "info", False),
        "success": ("✓", "success", False),
        "warning": ("⚠", "warning", False),
        "error": ("✗", "error", True),
    }

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

    def _get_timestamp(self, force: bool = False, for_log: bool = False) -> str:
        """
        Get formatted timestamp if enough time has passed.

        Args:
            force: If True, always return timestamp
            for_log: If True, use log-appropriate color

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
            ts_color = (
                PALETTE.timestamp["log"] if for_log else PALETTE.timestamp["console"]
            )
            return f"[{ts_color}][{now.strftime('%H:%M:%S')}][/{ts_color}] "

        return ""

    def _styled(
        self, text: str, color_key: str, for_log: bool = False, bold: bool = False
    ) -> str:
        """
        Wrap text with Rich markup using appropriate color scheme.

        Args:
            text: Text to wrap
            color_key: PALETTE attribute name (e.g., "header", "success")
            for_log: If True, use log colors; otherwise console colors
            bold: If True, add bold modifier

        Returns:
            Rich-formatted string
        """
        scheme = "log" if for_log else "console"
        color = getattr(PALETTE, color_key)[scheme]
        tag = f"bold {color}" if bold else color
        return f"[{tag}]{text}[/{tag}]"

    def _emit(self, msg_type: str, message: str, with_timestamp: bool = False):
        """
        Unified message emitter for both console and log.

        Args:
            msg_type: Key into _MESSAGE_TYPES (info, success, warning, error)
            message: Message text
            with_timestamp: If True, force timestamp on this message
        """
        emoji, color_key, bold = self._MESSAGE_TYPES[msg_type]

        # Console output
        if not self.quiet:
            ts = self._get_timestamp(force=with_timestamp, for_log=False)
            styled = self._styled(
                f"{emoji} {message}", color_key, for_log=False, bold=bold
            )
            self.console.print(f"{ts}{styled}")

        # Log output
        ts = self._get_timestamp(force=with_timestamp, for_log=True)
        styled = self._styled(f"{emoji} {message}", color_key, for_log=True, bold=bold)
        self.log_console.print(f"{ts}{styled}")

    def info(self, message: str, with_timestamp: bool = False):
        """Send info message (ℹ) to both console and log."""
        self._emit("info", message, with_timestamp)

    def success(self, message: str, with_timestamp: bool = False):
        """Send success message (✓) to both console and log."""
        self._emit("success", message, with_timestamp)

    def warning(self, message: str, with_timestamp: bool = False):
        """Send warning message (⚠) to both console and log."""
        self._emit("warning", message, with_timestamp)

    def error(self, message: str, with_timestamp: bool = False):
        """Send error message (✗) to both console and log."""
        self._emit("error", message, with_timestamp)

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
        # Console: Show section divider + optional tree
        if not self.quiet:
            header_c = PALETTE.header["console"]

            # Section divider
            self.console.print()
            self.console.rule(
                f"[{header_c}]Step {step_num}: {step_name}", style=header_c
            )

            # Optional: Show pipeline tree with current step highlighted
            if show_tree:
                from rich.tree import Tree

                tree = Tree(
                    f"🔄 [{PALETTE.highlight['console']}]Pipeline Steps",
                    guide_style="dim",
                )

                for i, step in enumerate(all_steps, 1):
                    if i < step_num:
                        tree.add(f"[{PALETTE.step_complete['console']}]✓ {step}")
                    elif i == step_num:
                        tree.add(f"[{PALETTE.step_current['console']}]▶ {step}")
                    else:
                        tree.add(f"[{PALETTE.step_pending['console']}]○ {step}")

                self.console.print(tree)

        # Log: Simple section marker with color
        sep_line = "=" * 80
        step_text = f"STEP {step_num}: {step_name.upper()}"
        timestamp = self._get_timestamp(force=True, for_log=True)
        self.log_console.print("")  # Blank line
        self.log_console.print(self._styled(sep_line, "header", for_log=True))
        self.log_console.print(
            self._styled(step_text, "header", for_log=True, bold=True)
        )
        self.log_console.print(self._styled(sep_line, "header", for_log=True))
        self.log_console.print(
            f"{timestamp}{self._styled(f'Starting step {step_num} of {len(all_steps)}', 'info', for_log=True)}"
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
            header.append("intronIC ", style=f"bold {PALETTE.highlight['console']}")
            header.append(f"v{ver}\n\n", style="dim")
            header.append("Run name: ", style=PALETTE.table_value["console"])
            header.append(f"{species_name}\n", style=PALETTE.warning["console"])
            header.append("Input: ", style=PALETTE.table_value["console"])
            header.append(f"{input_mode}\n", style=PALETTE.success["console"])
            header.append("Output: ", style=PALETTE.table_value["console"])
            header.append(f"{output_dir}\n", style=PALETTE.path["console"])
            header.append("Threshold: ", style=PALETTE.table_value["console"])
            header.append(f"{threshold}%\n", style=PALETTE.highlight["console"])

            if model_path:
                # Show just the filename for brevity
                import os

                model_name = os.path.basename(model_path)
                header.append("Model: ", style=PALETTE.table_value["console"])
                header.append(f"{model_name}\n", style=PALETTE.path["console"])

            header.append("Started: ", style=PALETTE.table_value["console"])
            header.append(f"{timestamp}", style=PALETTE.timestamp["console"])

            panel = Panel(
                header,
                border_style=PALETTE.highlight["console"],
                box=box.DOUBLE,
                padding=(1, 2),
            )
            self.console.print(panel)

        # Log: Clean text format with full paths and color
        sep_line = "=" * 80
        self.log_console.print(self._styled(sep_line, "header", for_log=True))
        self.log_console.print(
            self._styled(f"intronIC v{ver}", "highlight", for_log=True, bold=True)
        )
        self.log_console.print(
            self._styled(f"Started: {timestamp}", "timestamp", for_log=True)
        )
        self.log_console.print(self._styled(sep_line, "header", for_log=True))
        self.log_console.print("")

        # Combined Command and Configuration section
        self.log_console.print(
            self._styled(
                "Command and Configuration:", "header", for_log=True, bold=True
            )
        )

        # Helper for labeled log lines
        def log_item(label: str, value: str, value_style: str = None):
            styled_label = self._styled(f"{label}:", "table_header", for_log=True)
            if value_style:
                styled_value = self._styled(value, value_style, for_log=True)
            else:
                styled_value = value
            self.log_console.print(f"  {styled_label} {styled_value}")

        # Command and working directory
        if command:
            log_item("Command", command)
        if working_dir:
            log_item("Working directory", working_dir, "path")

        # Configuration details
        log_item("Run name", species_name)
        log_item("Input mode", input_mode)
        log_item("Classification threshold", f"{threshold}%", "highlight")
        log_item("Output directory", output_dir, "path")

        # Log input files with full paths
        if genome_path:
            log_item("Genome", genome_path, "path")
        if annotation_path:
            log_item("Annotation", annotation_path, "path")
        if bed_path:
            log_item("BED file", bed_path, "path")
        if sequences_path:
            log_item("Sequences", sequences_path, "path")
        if model_path:
            log_item("Model", model_path, "path")

        self.log_console.print("")

    def print_filtering_summary(
        self,
        total: int,
        short: int = 0,
        ambiguous: int = 0,
        noncanonical: int = 0,
        isoform: int = 0,
        overlap: int = 0,
        duplicates: int = 0,
        kept: int = 0,
        # User options that affect what gets excluded
        include_duplicates: bool = False,
        include_isoforms: bool = False,
        exclude_noncanonical: bool = False,
        exclude_overlap: bool = False,
    ):
        """
        Print filtering summary table with Included and Excluded columns.

        Each category count appears in either Included or Excluded column
        based on user options, making it clear what action was taken.

        Args:
            total: Total number of introns before filtering
            short: Number too short for scoring
            ambiguous: Number with ambiguous bases in scoring regions
            noncanonical: Number with non-canonical splice sites
            isoform: Number from alternative isoforms
            overlap: Number with overlapping coordinates
            duplicates: Number with duplicate coordinates
            kept: Number retained for scoring
            include_duplicates: User specified -d/--include-duplicates flag
            include_isoforms: User specified -i/--allow-multiple-isoforms flag
            exclude_noncanonical: User specified --no-nc flag
            exclude_overlap: User specified -v/--no-intron-overlap flag
        """
        from rich import box
        from rich.table import Table

        # Helper to format count in appropriate column
        def fmt_inc_exc(count: int, is_excluded: bool) -> tuple:
            """Return (included_str, excluded_str) for a count."""
            if count == 0:
                return ("0", "0")
            if is_excluded:
                return ("0", f"{count:,}")
            else:
                return (f"{count:,}", "0")

        # Determine which column each category goes in
        # Always excluded: short, ambiguous
        # Conditionally based on flags:
        dup_inc, dup_exc = fmt_inc_exc(duplicates, not include_duplicates)
        short_inc, short_exc = fmt_inc_exc(short, True)  # always excluded
        ambig_inc, ambig_exc = fmt_inc_exc(ambiguous, True)  # always excluded
        nc_inc, nc_exc = fmt_inc_exc(noncanonical, exclude_noncanonical)
        overlap_inc, overlap_exc = fmt_inc_exc(overlap, exclude_overlap)
        iso_inc, iso_exc = fmt_inc_exc(isoform, not include_isoforms)

        # Calculate totals
        total_excluded = (
            (short if True else 0)
            + (ambiguous if True else 0)
            + (noncanonical if exclude_noncanonical else 0)
            + (overlap if exclude_overlap else 0)
            + (isoform if not include_isoforms else 0)
            + (duplicates if not include_duplicates else 0)
        )

        # Console: Rich table
        if not self.quiet:
            table = Table(
                title="Intron Filtering Summary",
                box=box.ROUNDED,
                title_style=f"bold {PALETTE.highlight['console']}",
            )
            table.add_column(
                "Category",
                style=PALETTE.table_header["console"],
                no_wrap=False,
                width=28,
            )
            table.add_column(
                "Included",
                style=PALETTE.table_value["console"],
                justify="right",
                width=10,
            )
            table.add_column(
                "Excluded",
                style=PALETTE.table_value["console"],
                justify="right",
                width=10,
            )

            # Always show all categories to indicate which checks were performed
            table.add_row("  Duplicates", dup_inc, dup_exc)
            table.add_row("  Too short", short_inc, short_exc)
            table.add_row("  Ambiguous bases", ambig_inc, ambig_exc)
            table.add_row("  Non-canonical", nc_inc, nc_exc)
            table.add_row("  Overlapping", overlap_inc, overlap_exc)
            table.add_row("  Alternative isoform", iso_inc, iso_exc)

            # Separator and totals
            table.add_section()
            table.add_row(
                "Total excluded",
                "",
                f"{total_excluded:,}",
                style=f"bold {PALETTE.warning['console']}",
            )
            table.add_row(
                "Retained for scoring",
                "",
                f"{kept:,}",
                style=f"bold {PALETTE.success['console']}",
            )

            self.console.print()
            self.console.print(table)
            self.console.print()

        # Log: ASCII table with color
        def pad(s: str, width: int = 10) -> str:
            return s.rjust(width)

        border = self._styled
        tbl = "table_header"

        self.log_console.print("")
        self.log_console.print(
            self._styled(
                "Intron Filtering Summary:", "highlight", for_log=True, bold=True
            )
        )
        self.log_console.print(
            border(
                "┌────────────────────────────┬────────────┬────────────┐",
                tbl,
                for_log=True,
            )
        )
        self.log_console.print(
            border(
                "│ Category                   │ Included   │ Excluded   │",
                tbl,
                for_log=True,
            )
        )
        self.log_console.print(
            border(
                "├────────────────────────────┼────────────┼────────────┤",
                tbl,
                for_log=True,
            )
        )

        # Data rows (styled consistently with table borders)
        self.log_console.print(
            border(
                f"│   Duplicates               │ {pad(dup_inc)} │ {pad(dup_exc)} │",
                tbl,
                for_log=True,
            )
        )
        self.log_console.print(
            border(
                f"│   Too short                │ {pad(short_inc)} │ {pad(short_exc)} │",
                tbl,
                for_log=True,
            )
        )
        self.log_console.print(
            border(
                f"│   Ambiguous bases          │ {pad(ambig_inc)} │ {pad(ambig_exc)} │",
                tbl,
                for_log=True,
            )
        )
        self.log_console.print(
            border(
                f"│   Non-canonical            │ {pad(nc_inc)} │ {pad(nc_exc)} │",
                tbl,
                for_log=True,
            )
        )
        self.log_console.print(
            border(
                f"│   Overlapping              │ {pad(overlap_inc)} │ {pad(overlap_exc)} │",
                tbl,
                for_log=True,
            )
        )
        self.log_console.print(
            border(
                f"│   Alternative isoform      │ {pad(iso_inc)} │ {pad(iso_exc)} │",
                tbl,
                for_log=True,
            )
        )

        self.log_console.print(
            border(
                "├────────────────────────────┼────────────┼────────────┤",
                tbl,
                for_log=True,
            )
        )
        # Summary rows with highlighting (borders styled, values highlighted)
        border_char = self._styled("│", tbl, for_log=True)
        warn = self._styled("Total excluded", "warning", for_log=True)
        warn_val = self._styled(f"{total_excluded:>10,}", "warning", for_log=True)
        self.log_console.print(
            f"{border_char} {warn}             {border_char}            {border_char} {warn_val} {border_char}"
        )
        succ = self._styled("Retained for scoring", "success", for_log=True)
        succ_val = self._styled(f"{kept:>10,}", "success", for_log=True)
        self.log_console.print(
            f"{border_char} {succ}       {border_char}            {border_char} {succ_val} {border_char}"
        )
        self.log_console.print(
            border(
                "└────────────────────────────┴────────────┴────────────┘",
                tbl,
                for_log=True,
            )
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
                title_style=f"bold {PALETTE.highlight['console']}",
            )
            table.add_column(
                "Type", style=PALETTE.table_header["console"], no_wrap=True
            )
            table.add_column(
                "Count", style=PALETTE.table_value["console"], justify="right"
            )
            table.add_column(
                "Percentage", style=PALETTE.table_value["console"], justify="right"
            )

            table.add_row(
                f"[{PALETTE.u12_highlight['log']}]U12-type (total)",
                f"[{PALETTE.u12_highlight['log']}]{u12_count:,}",
                f"[{PALETTE.u12_highlight['log']}]{u12_pct:.2f}%",
            )
            table.add_row(
                f"[{PALETTE.u12_highlight['log']}]U12-type (AT-AC)",
                f"[{PALETTE.u12_highlight['log']}]{atac_count:,}",
                f"[{PALETTE.u12_highlight['log']}]{atac_pct:.2f}%",
            )
            table.add_row("U2-type", f"{u2_count:,}", f"{u2_pct:.2f}%")
            table.add_row("[bold]Total", f"[bold]{total:,}", "[bold]100.00%")

            self.console.print(table)

        # Log: ASCII table with color
        border = self._styled
        tbl = "table_header"
        u12 = lambda t: self._styled(t, "u12_highlight", for_log=True)

        self.log_console.print("")
        self.log_console.print(
            self._styled(
                f"Classification Results (threshold: {threshold}%):",
                "highlight",
                for_log=True,
                bold=True,
            )
        )
        self.log_console.print(
            border(
                "┌──────────────────────┬───────────┬────────────┐", tbl, for_log=True
            )
        )
        self.log_console.print(
            border(
                "│ Type                 │ Count     │ Percentage │", tbl, for_log=True
            )
        )
        self.log_console.print(
            border(
                "├──────────────────────┼───────────┼────────────┤", tbl, for_log=True
            )
        )
        self.log_console.print(
            f"│ {u12('U12-type (total)')}     │ {u12(f'{u12_count:>9,}')} │ {u12(f'{u12_pct:>9.2f}%')} │"
        )
        self.log_console.print(
            f"│ {u12('U12-type (AT-AC)')}     │ {u12(f'{atac_count:>9,}')} │ {u12(f'{atac_pct:>9.2f}%')} │"
        )
        self.log_console.print(
            f"│ U2-type              │ {u2_count:>9,} │ {u2_pct:>9.2f}% │"
        )
        self.log_console.print(
            border(
                "├──────────────────────┼───────────┼────────────┤", tbl, for_log=True
            )
        )
        self.log_console.print(
            f"│ [bold]Total[/bold]                │ [bold]{total:>9,}[/bold] │ [bold]{100.0:>9.2f}%[/bold] │"
        )
        self.log_console.print(
            border(
                "└──────────────────────┴───────────┴────────────┘", tbl, for_log=True
            )
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
                title_style=f"bold {PALETTE.highlight['console']}",
            )
            table.add_column(
                "Rank", style=PALETTE.table_header["console"], justify="right"
            )
            table.add_column(
                "Dinucleotide", style=PALETTE.table_value["console"], justify="right"
            )
            table.add_column(
                "Count", style=PALETTE.table_value["console"], justify="right"
            )
            table.add_column(
                "Percent", style=PALETTE.table_value["console"], justify="right"
            )

            for i, (dnts, count) in enumerate(boundaries_to_show, 1):
                pct = (count / total * 100) if total > 0 else 0
                table.add_row(str(i), dnts, f"{count:,}", f"{pct:.2f}%")

            self.console.print()
            self.console.print(table)

        # Log: ASCII table with color
        border = self._styled
        tbl = "table_header"

        self.log_console.print("")
        self.log_console.print(
            self._styled(
                f"Top {top_n} Splice Site Boundaries ({intron_type} introns):",
                "highlight",
                for_log=True,
                bold=True,
            )
        )
        self.log_console.print(
            border("┌──────┬──────────────┬──────────┬───────────┐", tbl, for_log=True)
        )
        self.log_console.print(
            border("│ Rank │ Dinucleotide │ Count    │ Percent   │", tbl, for_log=True)
        )
        self.log_console.print(
            border("├──────┼──────────────┼──────────┼───────────┤", tbl, for_log=True)
        )

        for i, (dnts, count) in enumerate(boundaries_to_show, 1):
            pct = (count / total * 100) if total > 0 else 0
            self.log_console.print(
                f"│ {i:>4} │ {dnts:>12} │ {count:>8,} │ {pct:>8.2f}% │"
            )

        self.log_console.print(
            border("└──────┴──────────────┴──────────┴───────────┘", tbl, for_log=True)
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

            tree = Tree(
                f"📁 [bold {PALETTE.highlight['console']}]Output Files",
                guide_style="dim",
            )

            for file_type, filepath in output_files.items():
                tree.add(
                    f"[{PALETTE.success['console']}]{file_type}:[/{PALETTE.success['console']}] [{PALETTE.path['console']}]{filepath}"
                )

            self.console.print(tree)

        # Log: ASCII tree with box-drawing characters
        self.log_console.print("")
        self.log_console.print(
            self._styled("Output Files:", "highlight", for_log=True, bold=True)
        )

        # Convert dict to list for easier handling of first/last items
        items = list(output_files.items())
        for i, (file_type, filepath) in enumerate(items):
            is_last = i == len(items) - 1
            prefix = "└──" if is_last else "├──"
            label = self._styled(f"{file_type}:", "success", for_log=True)
            path = self._styled(filepath, "path", for_log=True)
            self.log_console.print(f"  {prefix} {label} {path}")

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
        self.log_console.print("")
        self.log_console.print(
            self._styled("Training Configuration:", "header", for_log=True, bold=True)
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
            self._styled(
                "Hyperparameter Optimization:", "header", for_log=True, bold=True
            )
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
            self._styled(
                f"Nested Cross-Validation Results ({n_folds} folds):",
                "header",
                for_log=True,
                bold=True,
            )
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
            self._styled(
                "Train/Test Split Evaluation:", "header", for_log=True, bold=True
            )
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
            self._styled(
                f"Trained Ensemble ({n_models} models):",
                "header",
                for_log=True,
                bold=True,
            )
        )

        for i, model in enumerate(models, 1):
            self.log_console.print(f"  Model {i}/{n_models}:")
            self.log_console.print(
                f"    Training samples: {model.train_size:,} ({model.u12_count:,} U12, {model.u2_count:,} U2)"
            )
            self.log_console.print(f"    C parameter: {model.parameters.C:.6e}")

        self.log_console.print("")
