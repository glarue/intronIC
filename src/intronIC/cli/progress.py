"""
Rich progress reporter for intronIC CLI.

Provides colored console output and progress tracking using the rich library.
"""

from typing import Optional, Dict, Any
from rich.console import Console

from intronIC.cli.colors import PALETTE

# Import version from package
try:
    from importlib.metadata import version
    __version__ = version("intronIC")
except Exception:
    __version__ = "dev"  # Fallback for development
from rich.progress import (
    Progress,
    SpinnerColumn,
    TextColumn,
    BarColumn,
    TaskProgressColumn,
    TimeRemainingColumn,
    TimeElapsedColumn,
)
from rich.panel import Panel
from rich.table import Table
from rich.tree import Tree
from rich.text import Text
from rich import box


class IntronICProgressReporter:
    """Rich-based progress reporter for intronIC pipeline."""

    def __init__(self, quiet: bool = False):
        """Initialize progress reporter.

        Args:
            quiet: If True, suppress non-essential output
        """
        # Disable automatic highlighting to avoid false matches (e.g. "PUT" in "OUTPUT")
        self.console = Console(highlight=False)
        self.quiet = quiet
        self._progress: Optional[Progress] = None

    def print_header(self, species_name: str, input_mode: str):
        """Print styled header banner.

        Args:
            species_name: Name of species being analyzed
            input_mode: Input mode (annotation, bed, sequences)
        """
        if self.quiet:
            return

        header = Text()
        header.append("intronIC ", style=f"bold {PALETTE.highlight['console']}")
        header.append(f"v{__version__}", style="dim")
        header.append(" - Intron Classification Pipeline\n", style=f"bold {PALETTE.highlight['console']}")
        header.append(f"Species: ", style=PALETTE.table_value["console"])
        header.append(f"{species_name}\n", style=PALETTE.mustard_hex)
        header.append(f"Input: ", style=PALETTE.table_value["console"])
        header.append(f"{input_mode}", style=PALETTE.success["console"])

        panel = Panel(
            header,
            border_style=PALETTE.highlight["console"],
            box=box.DOUBLE,
            padding=(1, 2)
        )
        self.console.print(panel)

    def print_section(self, title: str, style: str = None):
        """Print section header.

        Args:
            title: Section title
            style: Rich style string (defaults to PALETTE.header)
        """
        if self.quiet:
            return

        if style is None:
            style = f"bold {PALETTE.header['console']}"

        self.console.print()
        self.console.rule(f"[{style}]{title}", style=style)

    def print_info(self, message: str, prefix: str = "ℹ"):
        """Print info message.

        Args:
            message: Message to print
            prefix: Prefix icon/text
        """
        if self.quiet:
            return

        color = PALETTE.info["console"]
        self.console.print(f"[{color}]{prefix}[/{color}] {message}")

    def print_success(self, message: str, prefix: str = "✓"):
        """Print success message.

        Args:
            message: Message to print
            prefix: Prefix icon/text
        """
        color = PALETTE.success["console"]
        self.console.print(f"[{color}]{prefix}[/{color}] {message}")

    def print_warning(self, message: str, prefix: str = "⚠"):
        """Print warning message.

        Args:
            message: Message to print
            prefix: Prefix icon/text
        """
        color = PALETTE.warning["console"]
        self.console.print(f"[{color}]{prefix}[/{color}] {message}")

    def print_error(self, message: str, prefix: str = "✗"):
        """Print error message.

        Args:
            message: Message to print
            prefix: Prefix icon/text
        """
        color = PALETTE.error["console"]
        self.console.print(f"[{color}]{prefix}[/{color}] {message}", style="bold")

    def print_stats_table(self, stats: Dict[str, Any], title: str = "Statistics"):
        """Print statistics as a formatted table.

        Args:
            stats: Dictionary of statistics
            title: Table title
        """
        if self.quiet:
            return

        table = Table(title=title, box=box.ROUNDED, title_style=f"bold {PALETTE.highlight['console']}")
        table.add_column("Metric", style=PALETTE.table_header["console"], no_wrap=True)
        table.add_column("Value", style=PALETTE.table_value["console"], justify="right")

        for key, value in stats.items():
            # Format values nicely
            if isinstance(value, float):
                if value < 1:
                    formatted = f"{value:.4f}"
                else:
                    formatted = f"{value:,.2f}"
            elif isinstance(value, int):
                formatted = f"{value:,}"
            else:
                formatted = str(value)

            table.add_row(key, formatted)

        self.console.print(table)

    def print_classification_summary(
        self,
        total: int,
        u12_count: int,
        u2_count: int,
        threshold: float
    ):
        """Print classification summary.

        Args:
            total: Total introns classified
            u12_count: Number of U12 introns
            u2_count: Number of U2 introns
            threshold: Classification threshold used
        """
        if self.quiet:
            return

        u12_pct = (u12_count / total * 100) if total > 0 else 0
        u2_pct = (u2_count / total * 100) if total > 0 else 0

        table = Table(
            title=f"Classification Results (threshold: {threshold}%)",
            box=box.DOUBLE_EDGE,
            title_style=f"bold {PALETTE.highlight['console']}"
        )
        table.add_column("Type", style=PALETTE.table_header["console"], no_wrap=True)
        table.add_column("Count", style=PALETTE.table_value["console"], justify="right")
        table.add_column("Percentage", style=PALETTE.table_value["console"], justify="right")

        table.add_row("[bold]U12-type", f"[bold]{u12_count:,}", f"[bold]{u12_pct:.2f}%")
        table.add_row("U2-type", f"{u2_count:,}", f"{u2_pct:.2f}%")
        table.add_row("Total", f"{total:,}", "100.00%")

        self.console.print(table)

    def print_file_tree(self, output_files: Dict[str, str]):
        """Print output files as a tree.

        Args:
            output_files: Dictionary mapping file type to file path
        """
        if self.quiet:
            return

        tree = Tree(f"📁 [bold {PALETTE.highlight['console']}]Output Files", guide_style="dim")

        for file_type, filepath in output_files.items():
            success_color = PALETTE.success["console"]
            path_color = PALETTE.path["console"]
            tree.add(f"[{success_color}]{file_type}:[/{success_color}] [{path_color}]{filepath}")

        self.console.print(tree)

    def create_progress(self) -> Progress:
        """Create and return a rich Progress instance.

        Returns:
            Configured Progress object
        """
        self._progress = Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
            console=self.console,
            transient=False,
        )
        return self._progress

    def print_pipeline_steps(self, steps: list[str], current_step: Optional[int] = None):
        """Print pipeline steps with current step highlighted.

        Args:
            steps: List of step descriptions
            current_step: Index of current step (0-based), or None for overview
        """
        if self.quiet:
            return

        tree = Tree(f"🔄 [bold {PALETTE.highlight['console']}]Pipeline Steps", guide_style="dim")

        for i, step in enumerate(steps, 1):
            if current_step is not None:
                if i < current_step:
                    # Completed
                    tree.add(f"[{PALETTE.step_complete['console']}]✓ {step}")
                elif i == current_step:
                    # Current
                    tree.add(f"[bold {PALETTE.step_current['console']}]▶ {step}")
                else:
                    # Pending
                    tree.add(f"[{PALETTE.step_pending['console']}]○ {step}")
            else:
                # Overview mode
                tree.add(f"[{PALETTE.table_value['console']}]{i}. {step}")

        self.console.print(tree)

    def print_model_metrics(self, metrics: Dict[str, float], model_index: int):
        """Print model training metrics.

        Args:
            metrics: Dictionary of metric names to values
            model_index: Model number (1-based)
        """
        if self.quiet:
            return

        table = Table(
            title=f"Model {model_index} Metrics",
            box=box.SIMPLE,
            title_style=f"bold {PALETTE.header['console']}"
        )
        table.add_column("Metric", style=PALETTE.table_header["console"])
        table.add_column("Value", style=PALETTE.table_value["console"], justify="right")

        for metric, value in metrics.items():
            if isinstance(value, float):
                formatted = f"{value:.4f}"
            else:
                formatted = str(value)
            table.add_row(metric, formatted)

        self.console.print(table)
