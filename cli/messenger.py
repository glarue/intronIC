"""
Unified messaging system for intronIC CLI.

Provides synchronized output to both console (with rich formatting) and log file
(with preserved ANSI colors) to eliminate redundant message calls.
"""

import logging
from typing import Optional, List
from rich.console import Console
from rich.progress import Progress
from rich.text import Text


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
        quiet: bool = False
    ):
        """
        Initialize unified messenger.

        Args:
            console: Rich console for terminal output
            log_console: Rich console for log file output (with ANSI)
            logger: Python logger (for log levels, but we use consoles for output)
            quiet: If True, suppress console output (but not logs)
        """
        self.console = console
        self.log_console = log_console
        self.logger = logger
        self.quiet = quiet

    def info(self, message: str):
        """
        Send info message to both console (with ℹ) and log (with preserved formatting).

        Args:
            message: Message text
        """
        formatted = f"[blue]ℹ[/blue] {message}"

        # Console: with formatting
        if not self.quiet:
            self.console.print(formatted)

        # Log: same formatting with ANSI preserved + timestamp via logger
        self.logger.info(message)

    def success(self, message: str):
        """
        Send success message to both console (with ✓) and log (with preserved formatting).

        Args:
            message: Message text
        """
        formatted = f"[green]✓[/green] {message}"

        # Console: with formatting
        if not self.quiet:
            self.console.print(formatted)

        # Log: same formatting with ANSI preserved
        self.logger.info(message)

    def warning(self, message: str):
        """
        Send warning message to both console (with ⚠) and log (with preserved formatting).

        Args:
            message: Message text
        """
        formatted = f"[yellow]⚠[/yellow] {message}"

        # Console: with formatting
        if not self.quiet:
            self.console.print(formatted)

        # Log: same formatting with ANSI preserved
        self.logger.warning(message)

    def error(self, message: str):
        """
        Send error message to both console (with ✗) and log (with preserved formatting).

        Args:
            message: Message text
        """
        formatted = f"[red]✗[/red] {message}"

        # Console: with formatting, bold
        if not self.quiet:
            self.console.print(formatted, style="bold")

        # Log: same formatting with ANSI preserved
        self.logger.error(message)

    def step(
        self,
        step_num: int,
        step_name: str,
        all_steps: List[str],
        show_tree: bool = True
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
            # Section divider
            self.console.print()
            self.console.rule(f"[bold blue]Step {step_num}: {step_name}", style="bold blue")

            # Optional: Show pipeline tree with current step highlighted
            if show_tree:
                from rich.tree import Tree
                tree = Tree("🔄 [bold cyan]Pipeline Steps", guide_style="dim")

                for i, step in enumerate(all_steps, 1):
                    if i < step_num:
                        # Completed
                        tree.add(f"[dim green]✓ {step}[/dim green]")
                    elif i == step_num:
                        # Current
                        tree.add(f"[bold yellow]▶ {step}[/bold yellow]")
                    else:
                        # Pending
                        tree.add(f"[dim white]○ {step}[/dim white]")

                self.console.print(tree)

        # Log: Section marker with formatting preserved
        step_header = f"[bold blue]{'=' * 40}[/bold blue]"
        step_text = f"[bold cyan]STEP {step_num}: {step_name.upper()}[/bold cyan]"
        self.log_console.print(f"\n{step_header}")
        self.log_console.print(step_text)
        self.log_console.print(step_header)

    def log_only(self, message: str, level: str = "info"):
        """
        Send message only to log file (not console).

        Useful for detailed debugging info that shouldn't clutter console.

        Args:
            message: Message text
            level: Log level (info, warning, error, debug)
        """
        log_func = getattr(self.logger, level.lower(), self.logger.info)
        log_func(message)

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

    def create_progress(self) -> Progress:
        """
        Create and return a rich Progress instance.

        Returns:
            Configured Progress object for progress bars
        """
        from rich.progress import (
            Progress,
            SpinnerColumn,
            TextColumn,
            BarColumn,
            TaskProgressColumn,
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
