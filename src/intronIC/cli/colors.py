"""
Color palette for intronIC CLI output.

Provides dual color schemes optimized for different output destinations:
- Console: Rich truecolor hex codes for beautiful rendering
- Log files: Named ANSI colors that downgrade gracefully

This ensures the best visual experience in both contexts.
"""

from dataclasses import dataclass
from typing import Literal


@dataclass(frozen=True)
class ColorPalette:
    """
    Dual color scheme for intronIC output.

    Each color property returns a dict with 'console' (hex) and 'log' (ANSI) keys.
    This provides a clean, unified interface for accessing context-appropriate colors.

    Console colors (hex, truecolor):
    - Yale Blue #094074: Dark blue, authoritative
    - Baltic Blue #3c6997: Medium blue, readable
    - Sky Aqua #5adbff: Bright cyan, attention-grabbing
    - Mustard #ffdd4a: Warm yellow, noticeable
    - Deep Saffron #fe9000: Warm orange, positive

    Log colors (named ANSI):
    - blue (34), bright_blue (94): Blues
    - bright_cyan (96): Highlights
    - yellow (33), bright_yellow (93): Warnings and success
    """

    # Primary hex colors (console)
    yale_blue_hex: str = "#094074"
    baltic_blue_hex: str = "#3c6997"
    sky_aqua_hex: str = "#5adbff"
    mustard_hex: str = "#ffdd4a"
    deep_saffron_hex: str = "#fe9000"

    # Primary ANSI colors (log)
    yale_blue_ansi: str = "blue"
    baltic_blue_ansi: str = "bright_blue"
    sky_aqua_ansi: str = "bright_cyan"
    mustard_ansi: str = "yellow"
    deep_saffron_ansi: str = "bright_yellow"

    def _color_dict(self, console: str, log: str) -> dict[Literal["console", "log"], str]:
        """Helper to create color dictionary."""
        return {"console": console, "log": log}

    # Semantic color mappings
    @property
    def info(self) -> dict[Literal["console", "log"], str]:
        """Info messages (ℹ)."""
        return self._color_dict(self.baltic_blue_hex, self.baltic_blue_ansi)

    @property
    def success(self) -> dict[Literal["console", "log"], str]:
        """Success messages (✓)."""
        return self._color_dict(self.deep_saffron_hex, self.deep_saffron_ansi)

    @property
    def warning(self) -> dict[Literal["console", "log"], str]:
        """Warning messages (⚠)."""
        return self._color_dict(self.mustard_hex, self.mustard_ansi)

    @property
    def error(self) -> dict[Literal["console", "log"], str]:
        """Error messages (✗)."""
        return self._color_dict("red", "red")

    @property
    def header(self) -> dict[Literal["console", "log"], str]:
        """Headers and titles."""
        return self._color_dict(self.yale_blue_hex, self.yale_blue_ansi)

    @property
    def highlight(self) -> dict[Literal["console", "log"], str]:
        """Numbers, important values."""
        return self._color_dict(self.sky_aqua_hex, self.sky_aqua_ansi)

    @property
    def path(self) -> dict[Literal["console", "log"], str]:
        """File paths."""
        return self._color_dict(self.baltic_blue_hex, self.baltic_blue_ansi)

    @property
    def step_current(self) -> dict[Literal["console", "log"], str]:
        """Current step in progress."""
        return self._color_dict(self.mustard_hex, self.mustard_ansi)

    @property
    def step_complete(self) -> dict[Literal["console", "log"], str]:
        """Completed steps."""
        return self._color_dict(
            f"dim {self.deep_saffron_hex}",
            f"dim {self.deep_saffron_ansi}"
        )

    @property
    def step_pending(self) -> dict[Literal["console", "log"], str]:
        """Pending steps."""
        return self._color_dict("dim white", "dim white")

    @property
    def table_header(self) -> dict[Literal["console", "log"], str]:
        """Table column headers."""
        return self._color_dict(self.sky_aqua_hex, self.sky_aqua_ansi)

    @property
    def table_value(self) -> dict[Literal["console", "log"], str]:
        """Table values."""
        return self._color_dict("white", "white")

    @property
    def u12_highlight(self) -> dict[Literal["console", "log"], str]:
        """U12 introns (rare, important)."""
        return self._color_dict(
            f"bold {self.mustard_hex}",
            f"bold {self.mustard_ansi}"
        )

    @property
    def timestamp(self) -> dict[Literal["console", "log"], str]:
        """Timestamps."""
        return self._color_dict("dim white", "dim white")


# Global palette instance
PALETTE = ColorPalette()
