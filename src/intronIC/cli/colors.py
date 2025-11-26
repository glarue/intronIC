"""
Color palette for intronIC CLI output.

Provides a unified, consistent color scheme across console and log outputs.
Palette: Yale Blue, Baltic Blue, Sky Aqua, Mustard, Deep Saffron
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class ColorPalette:
    """
    Unified color scheme for intronIC output.

    Colors from palette:
    - Yale Blue (#094074): Dark blue, authoritative
    - Baltic Blue (#3c6997): Medium blue, readable
    - Sky Aqua (#5adbff): Bright cyan, attention-grabbing
    - Mustard (#ffdd4a): Warm yellow, noticeable
    - Deep Saffron (#fe9000): Warm orange, positive
    """

    # Primary colors (Rich color names or hex)
    yale_blue: str = "#094074"
    baltic_blue: str = "#3c6997"
    sky_aqua: str = "#5adbff"
    mustard: str = "#ffdd4a"
    deep_saffron: str = "#fe9000"

    # Semantic mappings for message types
    @property
    def info(self) -> str:
        """Info messages (ℹ)."""
        return self.baltic_blue

    @property
    def success(self) -> str:
        """Success messages (✓)."""
        return self.deep_saffron

    @property
    def warning(self) -> str:
        """Warning messages (⚠)."""
        return self.mustard

    @property
    def error(self) -> str:
        """Error messages (✗)."""
        return "red"  # Keep red for errors (universal standard)

    @property
    def header(self) -> str:
        """Headers and titles."""
        return self.yale_blue

    @property
    def highlight(self) -> str:
        """Numbers, important values."""
        return self.sky_aqua

    @property
    def path(self) -> str:
        """File paths."""
        return self.baltic_blue

    @property
    def step_current(self) -> str:
        """Current step in progress."""
        return self.mustard

    @property
    def step_complete(self) -> str:
        """Completed steps."""
        return f"dim {self.deep_saffron}"

    @property
    def step_pending(self) -> str:
        """Pending steps."""
        return "dim white"

    @property
    def table_header(self) -> str:
        """Table column headers."""
        return self.sky_aqua

    @property
    def table_value(self) -> str:
        """Table values."""
        return "white"

    @property
    def u12_highlight(self) -> str:
        """U12 introns (rare, important)."""
        return f"bold {self.mustard}"

    @property
    def timestamp(self) -> str:
        """Timestamps."""
        return "dim white"


# Global palette instance
PALETTE = ColorPalette()
