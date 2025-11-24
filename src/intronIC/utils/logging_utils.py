"""
Enhanced logging utilities for intronIC.

Provides formatted logging with visual appeal (without requiring colors in files).
"""

import logging
from typing import Optional, TextIO
from datetime import datetime


class EnhancedFormatter(logging.Formatter):
    """
    Custom log formatter with simple visual enhancements.

    Makes plain text logs more readable with:
    - Section separators
    - Clean formatting
    """

    def __init__(self, width: int = 100):
        """
        Initialize formatter.

        Args:
            width: Total line width for formatting (default: 100)
        """
        super().__init__()
        self.width = width

    def format(self, record: logging.LogRecord) -> str:
        """
        Format a log record with visual enhancements.

        Args:
            record: Log record to format

        Returns:
            Formatted log message
        """
        timestamp = datetime.fromtimestamp(record.created).strftime('%Y-%m-%d %H:%M:%S')
        level = record.levelname.ljust(8)
        message = record.getMessage()

        # Check if this is a section header (starts with "==")
        if message.startswith('=='):
            return self._format_section(message)

        # Check if this is a subsection header (starts with "--")
        if message.startswith('--'):
            return self._format_subsection(message)

        # Regular message
        return f"[{timestamp}] {level} {message}"

    def _format_section(self, message: str) -> str:
        """Format a major section header."""
        # Extract the actual title (remove "==" markers)
        title = message.strip('= \n')

        # Simple format with separators
        sep = "=" * self.width
        return f"\n{sep}\n{title}\n{sep}"

    def _format_subsection(self, message: str) -> str:
        """Format a subsection header."""
        # Extract the actual title (remove "--" markers)
        title = message.strip('- \n')

        # Simple format with lighter separator
        sep = "-" * self.width
        return f"\n{title}\n{sep}"


class TrainingLogger:
    """
    Specialized logger for detailed training information.

    Captures all training details to a separate log file including:
    - Hyperparameter optimization progress
    - CV fold results
    - Model performance metrics
    - Detailed timing information
    """

    def __init__(self, filepath: str, width: int = 100):
        """
        Initialize training logger.

        Args:
            filepath: Path to training log file
            width: Line width for formatting (default: 100)
        """
        self.filepath = filepath
        self.width = width
        self.file: Optional[TextIO] = None
        self._indent_level = 0

    def __enter__(self):
        """Open log file."""
        self.file = open(self.filepath, 'w', encoding='utf-8')
        self._write_header()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Close log file."""
        if self.file:
            self._write_footer()
            self.file.close()
            self.file = None

    def _write_header(self):
        """Write file header."""
        BOX_H = '═'
        BOX_TL = '╔'
        BOX_TR = '╗'

        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        title = " intronIC Model Training Log "

        top = f"{BOX_TL}{BOX_H * (self.width - 2)}{BOX_TR}"
        mid = f"║{title.center(self.width - 2)}║"
        bot = f"║{f'Generated: {timestamp}'.center(self.width - 2)}║"
        sep = f"╚{BOX_H * (self.width - 2)}╝"

        self.file.write(f"{top}\n{mid}\n{bot}\n{sep}\n\n")

    def _write_footer(self):
        """Write file footer."""
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        self.file.write(f"\n{'─' * self.width}\n")
        self.file.write(f"Log completed: {timestamp}\n")

    def section(self, title: str):
        """Write a major section header."""
        BOX_H = '═'
        separator = f"\n{BOX_H * self.width}\n"
        self.file.write(f"{separator}{title.upper()}\n{separator}\n")
        self._indent_level = 0

    def subsection(self, title: str):
        """Write a subsection header."""
        separator = f"\n{'─' * self.width}\n"
        self.file.write(f"{separator}{title}\n{'─' * 40}\n")
        self._indent_level = 1

    def info(self, message: str, indent: Optional[int] = None):
        """Write an info message."""
        if indent is None:
            indent = self._indent_level
        prefix = '  ' * indent

        # Handle multi-line messages
        for line in message.split('\n'):
            self.file.write(f"{prefix}{line}\n")

    def table(self, headers: list, rows: list, indent: Optional[int] = None):
        """
        Write a formatted table.

        Args:
            headers: Column headers
            rows: List of row data (each row is a list)
            indent: Indentation level
        """
        if indent is None:
            indent = self._indent_level
        prefix = '  ' * indent

        # Calculate column widths
        col_widths = [len(str(h)) for h in headers]
        for row in rows:
            for i, cell in enumerate(row):
                col_widths[i] = max(col_widths[i], len(str(cell)))

        # Format header
        header_line = prefix + ' │ '.join(
            str(h).ljust(w) for h, w in zip(headers, col_widths)
        )
        separator = prefix + '─' * (len(header_line) - len(prefix))

        self.file.write(f"{header_line}\n{separator}\n")

        # Format rows
        for row in rows:
            row_line = prefix + ' │ '.join(
                str(cell).ljust(w) for cell, w in zip(row, col_widths)
            )
            self.file.write(f"{row_line}\n")

        self.file.write('\n')

    def metric(self, name: str, value: any, indent: Optional[int] = None):
        """Write a metric line (aligned key-value pair)."""
        if indent is None:
            indent = self._indent_level
        prefix = '  ' * indent

        # Align at column 40
        name_part = f"{name}:"
        self.file.write(f"{prefix}{name_part.ljust(40)} {value}\n")

    def blank(self):
        """Write a blank line."""
        self.file.write('\n')

    def flush(self):
        """Flush buffer to disk."""
        if self.file:
            self.file.flush()
