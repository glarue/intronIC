"""Tests for visualization helpers (pure logic + the raw-motif scatter gating)."""

import numpy as np

from intronIC.visualization.plots import (
    _legend_loc_by_density,
    plot_raw_motif_scatter_from_file,
    scatter_plot_from_arrays,
)

XLIM = (0.0, 10.0)
YLIM = (0.0, 10.0)

# score_info.iic columns the raw-motif scatter needs (order irrelevant — the function resolves by name)
_RAW_COLS = ["name", "5'_raw", "bp_raw", "P_motif", "motif_category", "type_id"]


def _write_score_info(path, rows, cols=_RAW_COLS):
    """Write a minimal tab-delimited score_info.iic from a list of dict rows."""
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "NA")) for c in cols) + "\n")


def _synthetic_rows(category, n=40, n_strong=6):
    """n introns in a U2 blob; the last n_strong are motif-strong (high P_motif). type_id honours
    the species call: all u2 when suppressed (NOT_DETECTED), else the strong ones are u12."""
    rng = np.random.RandomState(0)
    rows = []
    for i in range(n):
        strong = i >= n - n_strong
        pm = 0.98 if strong else 0.001
        if category == "NOT_DETECTED":
            tid = "u2"  # all suppressed
        else:
            tid = "u12" if strong else "u2"
        rows.append({
            "name": f"i{i}",
            "5'_raw": round(float(rng.normal(3 if strong else 0, 0.4)), 3),
            "bp_raw": round(float(rng.normal(2 if strong else 0, 0.4)), 3),
            "P_motif": pm,
            "motif_category": category,
            "type_id": tid,
        })
    return rows


class TestRawMotifScatterGating:
    """plot_raw_motif_scatter_from_file emits ONLY for NOT_DETECTED genomes and only when the
    adjudicator columns exist."""

    def test_emits_for_not_detected(self, tmp_path):
        sf = tmp_path / "sp.score_info.iic"
        _write_score_info(sf, _synthetic_rows("NOT_DETECTED"))
        out = plot_raw_motif_scatter_from_file(sf, tmp_path, "sp", threshold=90.0, fig_dpi=60)
        assert out is not None and out.exists()
        assert out.name == "sp.plot.scatter_raw.iic.png"

    def test_noop_for_detected(self, tmp_path):
        sf = tmp_path / "sp.score_info.iic"
        _write_score_info(sf, _synthetic_rows("DETECTED"))
        out = plot_raw_motif_scatter_from_file(sf, tmp_path, "sp", threshold=90.0, fig_dpi=60)
        assert out is None
        assert not (tmp_path / "sp.plot.scatter_raw.iic.png").exists()

    def test_noop_for_inconclusive(self, tmp_path):
        sf = tmp_path / "sp.score_info.iic"
        _write_score_info(sf, _synthetic_rows("INCONCLUSIVE"))
        assert plot_raw_motif_scatter_from_file(sf, tmp_path, "sp", threshold=90.0, fig_dpi=60) is None

    def test_noop_when_pmotif_column_absent(self, tmp_path):
        # raw_gated / legacy schema without the adjudicator columns -> nothing to draw
        cols = ["name", "5'_raw", "bp_raw", "type_id"]
        rows = [{"name": "i0", "5'_raw": 0.1, "bp_raw": 0.2, "type_id": "u2"}]
        sf = tmp_path / "sp.score_info.iic"
        _write_score_info(sf, rows, cols=cols)
        assert plot_raw_motif_scatter_from_file(sf, tmp_path, "sp", threshold=90.0, fig_dpi=60) is None


class TestDegenerateConfidenceBand:
    """A fully-suppressed genome has zero score spread; the legend must not render the nonsensical
    'X < U12-type <= X' bin (the C. elegans NOT_DETECTED regression)."""

    def test_zero_spread_scatter_renders(self, tmp_path):
        # All-equal tier scores -> score_stdev == 0 -> med_val == high_val (the degenerate case).
        rng = np.random.RandomState(1)
        xy = rng.normal(0, 1, size=(60, 2))
        svm = [0.0] * 60  # every point below threshold, zero spread
        out = tmp_path / "deg.png"
        scatter_plot_from_arrays(
            score_vector=xy, svm_scores=svm, species_name="deg", output_path=out,
            xlab="x", ylab="y", threshold=90.0, type_ids=["u2"] * 60, fig_dpi=60,
        )
        assert out.exists()  # collapsed-band guard: renders cleanly rather than a degenerate label


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
