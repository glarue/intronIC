"""Tests for visualization helpers (pure logic only — no figure rendering)."""

import numpy as np

from intronIC.visualization.plots import _legend_loc_by_density

XLIM = (0.0, 10.0)
YLIM = (0.0, 10.0)


class TestLegendLocByDensity:
    """Density-weighted legend placement: emptiest allowed corner, never upper-right,
    tie-break upper left > lower left > lower right."""

    def test_never_upper_right_even_when_emptiest(self):
        # Fill the three allowed corners; leave the upper-right region empty. It must NOT be chosen.
        pts = np.array([
            [0.2, 0.2], [0.3, 0.3],   # lower left
            [0.2, 9.8], [0.3, 9.7],   # upper left
            [9.8, 0.2], [9.7, 0.3],   # lower right
        ])
        assert _legend_loc_by_density(pts, XLIM, YLIM) != "upper right"

    def test_picks_emptiest_allowed_corner(self):
        # Pack upper-left and lower-left; leave lower-right empty -> lower right wins on count.
        pts = np.array([[0.5, 9.5]] * 20 + [[0.5, 0.5]] * 20)
        assert _legend_loc_by_density(pts, XLIM, YLIM) == "lower right"

    def test_tiebreak_prefers_upper_left(self):
        # All three allowed corners equally empty (no points in any corner box) -> upper left wins.
        pts = np.array([[5.0, 5.0]] * 50)  # dead-centre, outside every corner box
        assert _legend_loc_by_density(pts, XLIM, YLIM) == "upper left"

    def test_tiebreak_second_is_lower_left(self):
        # Put points only in the upper-left corner; lower-left and lower-right tie at 0 -> lower left wins.
        pts = np.array([[0.5, 9.5]] * 10)
        assert _legend_loc_by_density(pts, XLIM, YLIM) == "lower left"

    def test_empty_points_default_first_priority(self):
        assert _legend_loc_by_density(np.empty((0, 2)), XLIM, YLIM) == "upper left"

    def test_returns_valid_loc_string(self):
        rng = np.random.RandomState(0)
        pts = rng.normal(5, 2, size=(500, 2))
        assert _legend_loc_by_density(pts, XLIM, YLIM) in {"upper left", "lower left", "lower right"}
