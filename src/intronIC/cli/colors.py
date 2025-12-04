"""
Color palette for intronIC CLI output.

Provides dual color schemes optimized for different output destinations:
- Console: Rich truecolor hex codes for beautiful rendering
- Log files: Named ANSI colors that downgrade gracefully

This ensures the best visual experience in both contexts.
"""

from dataclasses import dataclass
from typing import ClassVar, Literal


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

    # Semantic color mappings: name -> (hex_attr, ansi_attr, prefix)
    # prefix is optional style modifier (e.g., "dim", "bold")
    _SEMANTIC_COLORS: ClassVar[dict] = {
        "info": ("baltic_blue_hex", "baltic_blue_ansi", ""),
        "success": ("deep_saffron_hex", "deep_saffron_ansi", ""),
        "warning": ("mustard_hex", "mustard_ansi", ""),
        "error": ("red", "red", ""),  # literals, not attrs
        "header": ("yale_blue_hex", "yale_blue_ansi", ""),
        "highlight": ("sky_aqua_hex", "sky_aqua_ansi", ""),
        "path": ("baltic_blue_hex", "baltic_blue_ansi", ""),
        "step_current": ("mustard_hex", "mustard_ansi", ""),
        "step_complete": ("deep_saffron_hex", "deep_saffron_ansi", "dim"),
        "step_pending": ("white", "white", "dim"),
        "table_header": ("sky_aqua_hex", "sky_aqua_ansi", ""),
        "table_value": ("white", "white", ""),
        "u12_highlight": ("mustard_hex", "mustard_ansi", "bold"),
        "timestamp": ("white", "white", "dim"),
    }

    def __getattr__(self, name: str) -> dict[Literal["console", "log"], str]:
        """
        Get semantic color dict by name.

        Returns dict with 'console' (hex) and 'log' (ANSI) keys.
        """
        if name.startswith("_"):
            raise AttributeError(f"'{type(self).__name__}' has no attribute '{name}'")

        if name in self._SEMANTIC_COLORS:
            hex_key, ansi_key, prefix = self._SEMANTIC_COLORS[name]

            # Get actual color values (attrs or literals)
            hex_val = getattr(self, hex_key, hex_key)
            ansi_val = getattr(self, ansi_key, ansi_key)

            # Apply prefix if present
            if prefix:
                return {"console": f"{prefix} {hex_val}", "log": f"{prefix} {ansi_val}"}
            return {"console": hex_val, "log": ansi_val}

        raise AttributeError(f"'{type(self).__name__}' has no attribute '{name}'")


# Global palette instance
PALETTE = ColorPalette()
