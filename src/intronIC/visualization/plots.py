"""
Plotting functions for intronIC visualization.

Recreates the plots from original intronIC:
- Density hexplot of intron scores
- Scatter plot with U12 classifications
- Score histogram
- Training reference plots
- Precision-Recall AUC curves
"""

import matplotlib
import numpy as np

matplotlib.use("Agg")  # Allow to run without X display server
from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from intronIC.core.intron import Intron


def _motif_standardizer_from_bundle(model_path) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """``(mean, scale)`` (length-3 arrays for [5'_raw, bp_raw, 3'_raw]) from the bundle's fitted
    in-pipeline GLOBAL ``StandardScaler`` — the very scaler the classifier applies internally.

    Plotting in this standardized space puts the three motif axes on a common scale (the space the RBF
    actually sees), un-stretching the disparate-range raw axes, while staying a faithful *per-axis affine*
    transform (so relative positions within each axis are preserved). The transformer's first three output
    columns are the base motif features passed through unchanged and ``StandardScaler`` is per-feature, so
    ``(x - mean[i]) / scale[i]`` standardizes each axis exactly. Returns ``None`` (caller falls back to raw
    axes) if the bundle/scaler can't be read.
    """
    if not model_path:
        return None
    try:
        import pickle

        with open(model_path, "rb") as fh:
            bundle = pickle.load(fh)
        ensemble = bundle["ensemble"] if isinstance(bundle, dict) else getattr(bundle, "ensemble", None)
        est = ensemble.models[0].model.calibrated_classifiers_[0].estimator
        scaler = est.named_steps.get("scale") if hasattr(est, "named_steps") else None
        if scaler is None:
            return None
        mean = np.asarray(scaler.mean_, dtype=float)
        scale = np.asarray(scaler.scale_, dtype=float)
        if mean.size < 3 or scale.size < 3 or not np.all(scale[:3] > 0):
            return None
        return mean[:3], scale[:3]
    except Exception:
        return None


#: Axis labels for the two plotting spaces (raw motif log-odds vs the classifier's standardized space).
_RAW_AXIS_LABELS = ("5'SS raw score", "BPS raw score", "3'SS raw score")
_STD_AXIS_LABELS = ("5'SS score (standardized)", "BPS score (standardized)", "3'SS score (standardized)")


def plot_classification_results_from_file(
    score_file: Path,
    output_dir: Path,
    species_name: str,
    threshold: float,
    fig_dpi: int = 300,
    model_path=None,
):
    """
    Generate classification plots by reading scores from output file.

    This is used in streaming mode where introns are not kept in memory.
    Reads the .score_info.iic file to extract z-scores and SVM scores.

    Args:
        score_file: Path to .score_info.iic file
        output_dir: Directory to save plots
        species_name: Species name for plot titles
        threshold: U12 classification threshold
        fig_dpi: Figure DPI for output images
    """
    # Read scores from file. Post z-stack removal the classifier feature space is the
    # background-corrected RAW motif log-odds (5'/BP/3'_raw), so we plot in that space.
    five_raw_scores = []
    bp_raw_scores = []
    three_raw_scores = []
    svm_scores = []
    type_id_list = []

    with open(score_file, "r") as f:
        header = f.readline().strip().split("\t")

        # Find column indices
        try:
            five_raw_idx = header.index("5'_raw")
            bp_raw_idx = header.index("bp_raw")
            svm_idx = header.index("svm_score")
        except ValueError as e:
            raise ValueError(f"Missing required column in score file: {e}")

        three_raw_idx = header.index("3'_raw") if "3'_raw" in header else None
        rel_idx = header.index("rel_score") if "rel_score" in header else None
        # v2.7+: prefer adjusted_score for the visual calling tier; falls
        # back to svm_score on rows where adjusted_score is "NA" or column
        # is absent. svm_score remains the raw classifier output.
        adj_idx = header.index("adjusted_score") if "adjusted_score" in header else None

        for line in f:
            fields = line.strip().split("\t")
            if len(fields) <= max(five_raw_idx, bp_raw_idx, svm_idx):
                continue

            # Parse raw motif scores (may be "NA")
            five_raw = fields[five_raw_idx]
            bp_raw = fields[bp_raw_idx]
            svm = fields[svm_idx]
            three_raw = fields[three_raw_idx] if three_raw_idx is not None else "NA"

            if five_raw != "NA" and bp_raw != "NA":
                try:
                    five_raw_scores.append(float(five_raw))
                    bp_raw_scores.append(float(bp_raw))
                    three_raw_scores.append(float(three_raw) if three_raw != "NA" else 0.0)
                except ValueError:
                    continue

            # Pick the score used for color tiers + counts: prefer
            # adjusted_score (v2.7+ continuous-discount output); fall back
            # to svm_score (raw classifier output).
            score_for_tier = svm
            if adj_idx is not None and adj_idx < len(fields):
                adj_val = fields[adj_idx]
                if adj_val not in ("NA", "", "null"):
                    score_for_tier = adj_val

            if score_for_tier != "NA":
                try:
                    svm_scores.append(float(score_for_tier))
                except ValueError:
                    continue

    if not five_raw_scores:
        # No valid scores to plot
        return

    score_vector = np.array(list(zip(five_raw_scores, bp_raw_scores)))
    score_vector_3d = np.array(list(zip(five_raw_scores, bp_raw_scores, three_raw_scores)))

    # Plot in the standardized space the classifier sees (the bundle's global StandardScaler), which
    # un-stretches the disparate-range raw motif axes. Faithful per-axis affine transform; falls back to
    # raw axes if the scaler can't be read.
    std = _motif_standardizer_from_bundle(model_path)
    if std is not None:
        mean, scale = std
        score_vector = (score_vector - mean[:2]) / scale[:2]
        score_vector_3d = (score_vector_3d - mean[:3]) / scale[:3]
        xlab, ylab, zlab = _STD_AXIS_LABELS
    else:
        xlab, ylab, zlab = _RAW_AXIS_LABELS

    # 1. Density hexplot
    hexplot_path = output_dir / f"{species_name}.plot.hex.iic.png"
    density_hexplot(
        score_vector,
        species_name=species_name,
        output_path=hexplot_path,
        xlab=xlab,
        ylab=ylab,
        fig_dpi=fig_dpi,
    )

    # 2. Scatter plot with U12 classification
    # Create minimal data structure for scatter plot
    scatter_path = output_dir / f"{species_name}.plot.scatter.iic.png"
    scatter_plot_from_arrays(
        score_vector=score_vector,
        svm_scores=svm_scores,
        species_name=species_name,
        output_path=scatter_path,
        xlab=xlab,
        ylab=ylab,
        threshold=threshold,
        fig_dpi=fig_dpi,
    )

    # 3. 3D scatter plot (5', BP, 3' motif axes)
    if three_raw_idx is not None:
        scatter_3d_path = output_dir / f"{species_name}.plot.scatter3d.iic.png"
        scatter_3d(
            score_vector_3d=score_vector_3d,
            svm_scores=svm_scores,
            species_name=species_name,
            output_path=scatter_3d_path,
            threshold=threshold,
            xlab=xlab,
            ylab=ylab,
            zlab=zlab,
            fig_dpi=fig_dpi,
        )

    # 4. Score histogram
    if svm_scores:
        hist_path = output_dir / f"{species_name}.plot.score_histogram.iic.png"
        histogram(
            svm_scores,
            threshold=threshold,
            species_name=species_name,
            output_path=hist_path,
            fig_dpi=fig_dpi,
        )


def plot_classification_results(
    introns: List[Intron],
    output_dir: Path,
    species_name: str,
    threshold: float,
    fig_dpi: int = 300,
    model_path=None,
):
    """
    Generate all classification result plots.

    Creates three plots matching original intronIC output:
    1. Density hexplot of all intron scores
    2. Scatter plot with U12s colored by confidence
    3. Histogram of SVM scores

    Args:
        introns: List of classified introns
        output_dir: Directory to save plots
        species_name: Species name for plot titles
        threshold: U12 classification threshold
        fig_dpi: Figure DPI for output images
        model_path: Optional path to the bundle whose global StandardScaler defines the
            standardized plotting space (falls back to raw motif axes if absent).
    """
    # Extract score vectors (5' raw motif, BP raw motif)
    score_vector = []
    for intron in introns:
        if (
            intron.scores
            and intron.scores.five_raw_score is not None
            and intron.scores.bp_raw_score is not None
        ):
            score_vector.append([intron.scores.five_raw_score, intron.scores.bp_raw_score])

    score_vector = np.array(score_vector)

    # Plot in the classifier's standardized space (bundle's global StandardScaler) so the disparate-range
    # raw motif axes are on a common scale; faithful per-axis affine transform, raw fallback if absent.
    std = _motif_standardizer_from_bundle(model_path)
    if std is not None and score_vector.size:
        mean, scale = std
        score_vector = (score_vector - mean[:2]) / scale[:2]
        xlab, ylab, zlab = _STD_AXIS_LABELS
    else:
        xlab, ylab, zlab = _RAW_AXIS_LABELS

    # 1. Density hexplot
    hexplot_path = output_dir / f"{species_name}.plot.hex.iic.png"
    density_hexplot(
        score_vector,
        species_name=species_name,
        output_path=hexplot_path,
        xlab=xlab,
        ylab=ylab,
        fig_dpi=fig_dpi,
    )

    # 2. Scatter plot with U12 classification
    scatter_path = output_dir / f"{species_name}.plot.scatter.iic.png"
    scatter_plot(
        introns,
        score_vector,
        species_name=species_name,
        output_path=scatter_path,
        xlab=xlab,
        ylab=ylab,
        threshold=threshold,
        fig_dpi=fig_dpi,
    )

    # 3. 3D scatter plot (5'z, BPz, 3'z)
    score_vector_3d = []
    for intron in introns:
        if (
            intron.scores
            and intron.scores.five_raw_score is not None
            and intron.scores.bp_raw_score is not None
            and intron.scores.three_raw_score is not None
        ):
            score_vector_3d.append([
                intron.scores.five_raw_score,
                intron.scores.bp_raw_score,
                intron.scores.three_raw_score,
            ])

    if score_vector_3d:
        score_vector_3d = np.array(score_vector_3d)
        if std is not None:
            score_vector_3d = (score_vector_3d - std[0][:3]) / std[1][:3]
        # v2.7+: prefer adjusted_score (continuous-discount-adjusted call
        # score) for color tiers + counts; fall back to svm_score.
        def _tier_score(intron):
            if intron.scores is None:
                return None
            adj = getattr(intron.scores, "adjusted_score", None)
            if adj is not None:
                return adj
            return intron.scores.svm_score
        svm_scores_for_3d = [
            _tier_score(i)
            for i in introns
            if i.scores and i.scores.svm_score is not None
            and i.scores.five_raw_score is not None
            and i.scores.bp_raw_score is not None
            and i.scores.three_raw_score is not None
        ]
        type_ids_for_3d = [
            i.metadata.type_id if i.metadata else None
            for i in introns
            if i.scores and i.scores.svm_score is not None
            and i.scores.five_raw_score is not None
            and i.scores.bp_raw_score is not None
            and i.scores.three_raw_score is not None
        ]
        scatter_3d_path = output_dir / f"{species_name}.plot.scatter3d.iic.png"
        scatter_3d(
            score_vector_3d=score_vector_3d,
            svm_scores=svm_scores_for_3d,
            species_name=species_name,
            output_path=scatter_3d_path,
            threshold=threshold,
            type_ids=type_ids_for_3d,
            xlab=xlab,
            ylab=ylab,
            zlab=zlab,
            fig_dpi=fig_dpi,
        )

    # 4. Score histogram (use adjusted_score if available, else svm_score)
    svm_scores = [
        (getattr(i.scores, "adjusted_score", None)
         if getattr(i.scores, "adjusted_score", None) is not None
         else i.scores.svm_score)
        for i in introns
        if i.scores and i.scores.svm_score is not None
    ]
    hist_path = output_dir / f"{species_name}.plot.score_histogram.iic.png"
    histogram(
        svm_scores,
        threshold=threshold,
        species_name=species_name,
        output_path=hist_path,
        fig_dpi=fig_dpi,
    )


def density_hexplot(
    scores: np.ndarray,
    species_name: str,
    output_path: Path,
    xlab: Optional[str] = None,
    ylab: Optional[str] = None,
    fsize: int = 14,
    fig_dpi: int = 300,
):
    """
    Create a density hexbin plot of intron scores.

    Args:
        scores: Nx2 array of (x, y) scores
        species_name: Species name from -n argument (for title)
        output_path: Full path where plot should be saved
        xlab: X-axis label
        ylab: Y-axis label
        fsize: Font size
        fig_dpi: Figure DPI
    """
    plt.figure(figsize=(8, 8))
    ax = plt.gca()

    hx = ax.hexbin(*scores.T, mincnt=1, cmap="inferno", bins="log", linewidths=0)

    # Note: Don't set aspect='equal' - let matplotlib auto-scale based on data ranges
    # Original v1.5.1 used auto aspect for hexplot

    # Clean title: species_name + description + count
    plot_title = f"{species_name} - Motif Score Density (n={len(scores)})"

    if xlab:
        plt.xlabel(xlab, fontsize=fsize)
    if ylab:
        plt.ylabel(ylab, fontsize=fsize)
    plt.title(plot_title, fontsize=fsize)

    # Add colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(hx, cax=cax)
    cb.set_label("Bin density (log10(n))")

    # Save figure
    plt.savefig(output_path, dpi=fig_dpi, bbox_inches="tight")
    plt.close()


def scatter_plot(
    introns: List[Intron],
    scores: np.ndarray,
    species_name: str,
    output_path: Path,
    xlab: str,
    ylab: str,
    threshold: float,
    fsize: int = 14,
    fig_dpi: int = 300,
):
    """
    Create a scatter plot with U12s colored by confidence level and marginal distributions.

    Wrapper around scatter_plot_from_arrays that extracts scores from Intron objects.

    Args:
        introns: List of classified introns
        scores: Nx2 array of (x, y) scores
        species_name: Species name from -n argument (for title)
        output_path: Full path where plot should be saved
        xlab: X-axis label
        ylab: Y-axis label
        threshold: U12 classification threshold
        fsize: Font size
        fig_dpi: Figure DPI
    """
    # Extract SVM scores and type classifications from introns
    svm_scores = []
    type_ids = []
    for intron in introns:
        if intron.scores and intron.scores.svm_score is not None:
            svm_scores.append(intron.scores.svm_score)
            type_id = intron.metadata.type_id if intron.metadata else None
            type_ids.append(type_id)

    scatter_plot_from_arrays(
        score_vector=scores,
        svm_scores=svm_scores,
        species_name=species_name,
        output_path=output_path,
        xlab=xlab,
        ylab=ylab,
        threshold=threshold,
        type_ids=type_ids,
        fsize=fsize,
        fig_dpi=fig_dpi,
    )


def scatter_plot_from_arrays(
    score_vector: np.ndarray,
    svm_scores: List[float],
    species_name: str,
    output_path: Path,
    xlab: str,
    ylab: str,
    threshold: float,
    type_ids: Optional[List[Optional[str]]] = None,
    fsize: int = 14,
    fig_dpi: int = 300,
):
    """
    Create a scatter plot with U12s colored by confidence level and marginal distributions.

    Args:
        score_vector: Nx2 array of (5' z-score, BP z-score)
        svm_scores: List of SVM scores (same length as score_vector)
        species_name: Species name for title
        output_path: Full path where plot should be saved
        xlab: X-axis label
        ylab: Y-axis label
        threshold: U12 classification threshold
        type_ids: Optional list of type classifications ('u2', 'u12', or None).
            If provided, uses these for U2/U12 classification instead of score threshold.
        fsize: Font size
        fig_dpi: Figure DPI
    """
    # Create figure with GridSpec for marginal distributions
    # Use 2x2 grid with symmetric ratios so the main plot region is square
    fig = plt.figure(figsize=(8, 8))
    gs = gridspec.GridSpec(
        2,
        2,
        figure=fig,
        width_ratios=[4, 1],  # Main plot : right marginal
        height_ratios=[1, 4],  # Top marginal : main plot
        hspace=0.02,
        wspace=0.02,
        top=0.92,  # Leave space for suptitle
        bottom=0.08,
        left=0.10,
        right=0.95,
    )

    ax_main = fig.add_subplot(gs[1, 0])
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax_main)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main)

    # Calculate threshold boundaries for U12 confidence levels
    score_stdev = np.std(svm_scores) if svm_scores else 10.0
    high_val = threshold
    med_val = threshold - score_stdev

    # Assign each point to a confidence tier; plotted in z-order below (U2 cloud first,
    # putative U12s last) so low/med/high-confidence U12s are never buried under the
    # dense U2 markers.
    tier_idx = {"u2": [], "low": [], "med": [], "high": []}
    for i, score in enumerate(svm_scores):
        if i >= len(score_vector):
            break
        # Determine if U2 or U12 based on type_ids if provided, else use score
        if type_ids is not None and i < len(type_ids):
            is_u2 = type_ids[i] == "u2"
        else:
            is_u2 = score < 50  # Raw classifier threshold

        if is_u2:
            tier_idx["u2"].append(i)
        elif score > high_val:
            tier_idx["high"].append(i)
        elif med_val < score <= high_val:
            tier_idx["med"].append(i)
        else:
            tier_idx["low"].append(i)

    u2_count = len(tier_idx["u2"])
    u12_low = len(tier_idx["low"])
    u12_med = len(tier_idx["med"])
    u12_high = len(tier_idx["high"])
    plot_scores = score_vector

    # Create legend
    legend_colors = ["xkcd:medium grey", "xkcd:red", "xkcd:orange", "xkcd:green"]
    legend_labels = [
        "U2",
        f"U12<={int(med_val)}",
        f"{int(med_val)}<U12<={int(high_val)}",
        f"U12>{int(high_val)}",
    ]
    legend_counts = [u2_count, u12_low, u12_med, u12_high]

    legend_patches = []
    for label, count, color in zip(legend_labels, legend_counts, legend_colors):
        label_with_count = f"{label} ({count})"
        patch = mpatches.Patch(color=color, label=label_with_count)
        legend_patches.append(patch)

    # Compute symmetric square limits UP FRONT: the U2 hexbin needs an explicit extent to align with
    # the U12 scatter, and the markers + marginals all share these limits.
    x_data = plot_scores[:, 0]
    y_data = plot_scores[:, 1]
    max_range = max(x_data.max() - x_data.min(), y_data.max() - y_data.min())
    x_center = (x_data.max() + x_data.min()) / 2
    y_center = (y_data.max() + y_data.min()) / 2
    margin = max_range * 0.05  # 5% margin
    xlim = (x_center - max_range / 2 - margin, x_center + max_range / 2 + margin)
    ylim = (y_center - max_range / 2 - margin, y_center + max_range / 2 + margin)

    # U2 set: grayscale hexbin DENSITY (was a flat grey scatter) — shows where the U2 bulk concentrates
    # without burying the U12s, and reads cleanly even where the cloud is millions-dense. log color scale
    # since the U2 core is orders of magnitude denser than its fringe.
    u2_idx = tier_idx["u2"]
    if u2_idx:
        u2_pts = plot_scores[u2_idx]
        ax_main.hexbin(
            u2_pts[:, 0], u2_pts[:, 1],
            gridsize=60, cmap="Greys", bins="log", mincnt=1,
            extent=(xlim[0], xlim[1], ylim[0], ylim[1]),
            linewidths=0, zorder=1,
        )

    # U12 tiers as a colored scatter ON TOP of the U2 density (low -> med -> high z-order), so putative
    # U12s of every confidence are drawn above the U2 background.
    for name, color, z, alpha in (
        ("low", "xkcd:red", 3, 0.85),
        ("med", "xkcd:orange", 4, 0.9),
        ("high", "xkcd:green", 5, 0.9),
    ):
        idx = tier_idx[name]
        if idx:
            pts = plot_scores[idx]
            ax_main.scatter(
                pts[:, 0], pts[:, 1], s=20, c=color, alpha=alpha,
                zorder=z, edgecolors="none", rasterized=True,
            )

    ax_main.legend(handles=legend_patches, fontsize=fsize - 2)
    ax_main.set_xlabel(xlab, fontsize=fsize)
    ax_main.set_ylabel(ylab, fontsize=fsize)

    ax_main.set_xlim(xlim)
    ax_main.set_ylim(ylim)

    # Set equal aspect with adjustable='box' to make plot physically square
    ax_main.set_aspect("equal", adjustable="box")

    # Plot marginal distributions with explicit range to match symmetric limits
    ax_top.hist(
        plot_scores[:, 0],
        bins=50,
        range=xlim,
        color="steelblue",
        alpha=0.7,
        edgecolor="none",
    )
    ax_top.set_xlim(xlim)  # Explicitly set to match main plot
    ax_top.set_ylabel("Count", fontsize=fsize - 2)
    ax_top.tick_params(labelbottom=False)
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)

    ax_right.hist(
        plot_scores[:, 1],
        bins=50,
        range=ylim,
        orientation="horizontal",
        color="steelblue",
        alpha=0.7,
        edgecolor="none",
    )
    ax_right.set_ylim(ylim)  # Explicitly set to match main plot
    ax_right.set_xlabel("Count", fontsize=fsize - 2)
    ax_right.tick_params(labelleft=False, labelrotation=45)
    ax_right.spines["top"].set_visible(False)
    ax_right.spines["right"].set_visible(False)

    fig.suptitle(
        f"{species_name} - U12 Classification Results",
        fontsize=fsize + 2,
        y=0.98,
        weight="bold",
    )

    plt.savefig(output_path, dpi=fig_dpi, bbox_inches="tight")
    plt.close()


def scatter_3d(
    score_vector_3d: np.ndarray,
    svm_scores: List[float],
    species_name: str,
    output_path: Path,
    threshold: float,
    type_ids: Optional[List[Optional[str]]] = None,
    xlab: str = "5'SS raw score",
    ylab: str = "BPS raw score",
    zlab: str = "3'SS raw score",
    fsize: int = 12,
    fig_dpi: int = 300,
):
    """
    Create a 3D scatter plot showing 5'z, BPz, and 3'z simultaneously.

    Displays two viewing angles side by side so all three axes are visible.
    U12 calls are colored by confidence level; U2 introns are shown as a
    translucent cloud for context.

    Args:
        score_vector_3d: Nx3 array of (5'z, BPz, 3'z)
        svm_scores: List of SVM scores (same length as score_vector_3d)
        species_name: Species name for title
        output_path: Full path where plot should be saved
        threshold: U12 classification threshold
        type_ids: Optional list of type classifications ('u2', 'u12', or None)
        fsize: Font size
        fig_dpi: Figure DPI
    """
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    n_points = min(len(score_vector_3d), len(svm_scores))
    scores_3d = score_vector_3d[:n_points]
    svms = svm_scores[:n_points]

    # Classify points
    score_stdev = np.std(svms) if svms else 10.0
    high_val = threshold
    med_val = threshold - score_stdev

    u2_idx, u12_low_idx, u12_med_idx, u12_high_idx = [], [], [], []

    for i, score in enumerate(svms):
        if type_ids is not None and i < len(type_ids):
            is_u2 = type_ids[i] == "u2"
        else:
            is_u2 = score < 50

        if is_u2:
            u2_idx.append(i)
        elif score > high_val:
            u12_high_idx.append(i)
        elif med_val < score <= high_val:
            u12_med_idx.append(i)
        else:
            u12_low_idx.append(i)

    # Subsample U2 for performance (max 20K points)
    if len(u2_idx) > 20000:
        rng = np.random.RandomState(42)
        u2_idx = list(rng.choice(u2_idx, 20000, replace=False))

    fig = plt.figure(figsize=(16, 7))

    # Two viewing angles
    views = [(25, -50), (25, 40)]
    view_labels = ["View 1", "View 2"]

    # Pre-compute data extents (so origin-axis lines span the visible range)
    if len(scores_3d) > 0:
        xmin, ymin, zmin = scores_3d.min(axis=0)
        xmax, ymax, zmax = scores_3d.max(axis=0)
        # Pad slightly so axis lines extend to the edges of the bounding box
        pad = 0.05
        ranges = [
            (xmin - pad * (xmax - xmin), xmax + pad * (xmax - xmin)),
            (ymin - pad * (ymax - ymin), ymax + pad * (ymax - ymin)),
            (zmin - pad * (zmax - zmin), zmax + pad * (zmax - zmin)),
        ]
    else:
        ranges = [(-1, 1)] * 3

    for panel, (elev, azim) in enumerate(views):
        ax = fig.add_subplot(1, 2, panel + 1, projection="3d")

        # Origin axis lines (z=0 planes): three darker lines through the
        # 3D origin, spanning the data extent on each axis. Draw FIRST so
        # data points overlay them.
        origin_color = "0.30"   # dark grey
        origin_lw = 1.5
        origin_alpha = 0.85
        ax.plot([ranges[0][0], ranges[0][1]], [0, 0], [0, 0],
                color=origin_color, linewidth=origin_lw,
                alpha=origin_alpha, zorder=0)
        ax.plot([0, 0], [ranges[1][0], ranges[1][1]], [0, 0],
                color=origin_color, linewidth=origin_lw,
                alpha=origin_alpha, zorder=0)
        ax.plot([0, 0], [0, 0], [ranges[2][0], ranges[2][1]],
                color=origin_color, linewidth=origin_lw,
                alpha=origin_alpha, zorder=0)
        # Origin marker — small black dot at (0,0,0) to fix the eye
        ax.scatter([0], [0], [0], s=30, c="black",
                   marker="o", zorder=1, alpha=0.9, edgecolors="white",
                   linewidth=0.5)

        # Plot U2 as translucent grey
        if u2_idx:
            u2_pts = scores_3d[u2_idx]
            ax.scatter(
                u2_pts[:, 0], u2_pts[:, 1], u2_pts[:, 2],
                s=1, c="xkcd:medium grey", alpha=0.03, rasterized=True,
            )

        # Plot U12 by confidence tier (low → high so high draws on top)
        for idx_list, color, label, size in [
            (u12_low_idx, "xkcd:red", f"U12<={int(med_val)}", 15),
            (u12_med_idx, "xkcd:orange", f"{int(med_val)}<U12<={int(high_val)}", 20),
            (u12_high_idx, "xkcd:green", f"U12>{int(high_val)}", 25),
        ]:
            if idx_list:
                pts = scores_3d[idx_list]
                ax.scatter(
                    pts[:, 0], pts[:, 1], pts[:, 2],
                    s=size, c=color, alpha=0.6, edgecolors="none",
                    label=f"{label} ({len(idx_list)})",
                )

        ax.set_xlabel(xlab, fontsize=fsize, labelpad=8)
        ax.set_ylabel(ylab, fontsize=fsize, labelpad=8)
        ax.set_zlabel(zlab, fontsize=fsize, labelpad=8)
        ax.tick_params(labelsize=fsize - 2)
        ax.view_init(elev=elev, azim=azim)

        if panel == 0:
            ax.legend(
                fontsize=fsize - 2, loc="upper left",
                markerscale=2, framealpha=0.8,
            )

    # Note: color tiers + counts derive from whatever score the caller
    # passed in `svm_scores`. Production callers (v2.7+) pass adjusted_score
    # (continuous-discount-adjusted call score). Point positions are always
    # 3D z-scores so the geometry is invariant across runs.
    _space = "standardized" if "standardized" in xlab else "raw"
    fig.suptitle(
        f"{species_name} - U12 Classification (3D {_space} scores)",
        fontsize=fsize + 2, weight="bold", y=0.97,
    )

    plt.savefig(output_path, dpi=fig_dpi, bbox_inches="tight")
    plt.close()


def histogram(
    data_list: List[float],
    threshold: float,
    species_name: str,
    output_path: Path,
    grid: bool = True,
    bins: int = 100,
    log: bool = True,
    fig_dpi: int = 300,
):
    """
    Create a histogram of SVM scores with threshold line.

    Args:
        data_list: List of SVM scores
        threshold: U12 classification threshold
        species_name: Species name from -n argument (for title)
        output_path: Full path where plot should be saved
        grid: Show grid lines
        bins: Number of histogram bins
        log: Use log scale for y-axis
        fig_dpi: Figure DPI
    """
    plt.figure(figsize=(10, 6))

    if log:
        plt.yscale("log")

    plt.hist(data_list, bins=bins)

    if grid:
        plt.grid(True, which="both", ls="--", alpha=0.7)

    # Clean title: species_name + description
    plt.title(f"{species_name} - U12 Score Distribution", fontsize=14)

    plt.xlabel("U12 score", fontsize=14)
    plt.ylabel("Number of introns", fontsize=14)

    # Add threshold line
    plt.axvline(
        threshold, color="orange", linestyle="--", label=f"U12 threshold: {threshold}"
    )

    plt.legend()
    plt.tight_layout()

    # Save figure
    plt.savefig(output_path, dpi=fig_dpi)
    plt.close()


def plot_training_results(
    u2_scores: np.ndarray,
    u12_scores: np.ndarray,
    pr_curves: List[Tuple[np.ndarray, np.ndarray]],
    pr_auc: float,
    output_dir: Path,
    species_name: str,
    fig_dpi: int = 300,
):
    """
    Generate training reference plots.

    Creates three plots:
    1. Reference scatter plot (U2 vs U12 training data)
    2. Reference hexplot (density of reference data)
    3. Precision-Recall AUC curve with aggregated statistics

    For nested CV (multiple curves):
    - Individual fold curves shown in light gray
    - Mean PR curve computed via interpolation (bold blue line)
    - Confidence bands (±1 std dev) as shaded region

    For split evaluation (single curve):
    - Single PR curve shown in bold blue

    Args:
        u2_scores: Nx2 array of U2 reference scores (5' z, BP z)
        u12_scores: Nx2 array of U12 reference scores (5' z, BP z)
        pr_curves: List of (precision, recall) tuples. Can be single curve or multiple (from CV folds)
        pr_auc: Average Precision-Recall AUC score
        output_dir: Directory to save plots
        species_name: Species name for plot titles
        fig_dpi: Figure DPI
    """
    # 1. Reference scatter plot
    ref_scatter(
        u2_scores,
        u12_scores,
        species_name=species_name,
        output_dir=output_dir,
        fig_dpi=fig_dpi,
    )

    # 2. Reference hexplot
    ref_hex_path = output_dir / f"{species_name}.ref_hex.iic.png"
    combined_scores = np.concatenate((u2_scores, u12_scores))
    density_hexplot(
        combined_scores,
        species_name=species_name,
        output_path=ref_hex_path,
        xlab="5'SS raw score",
        ylab="BPS raw score",
        fig_dpi=fig_dpi,
    )

    # 3. Precision-Recall curve with aggregated statistics
    plt.figure(figsize=(8, 8))

    # If multiple curves (nested CV), show individual folds + mean + confidence bands
    if len(pr_curves) > 1:
        # Common recall grid for interpolation (0 to 1)
        recall_grid = np.linspace(0, 1, 100)

        # Interpolate each fold to common grid
        interp_precisions = []
        for precision, recall in pr_curves:
            # Ensure monotonic decreasing recall (required for interpolation)
            # PR curves typically go from high recall (1.0) to low recall (0.0)
            if recall[0] < recall[-1]:
                # Already increasing, reverse it
                recall = recall[::-1]
                precision = precision[::-1]

            # Interpolate to common grid
            # Fill with edge values for out-of-bounds recall
            interp_p = np.interp(recall_grid[::-1], recall, precision)[::-1]
            interp_precisions.append(interp_p)

        # Compute mean and std
        prec_array = np.array(interp_precisions)
        mean_prec = prec_array.mean(axis=0)
        std_prec = prec_array.std(axis=0)

        # Plot individual fold curves in light gray
        for precision, recall in pr_curves:
            plt.plot(recall, precision, color="lightgray", alpha=0.5, linewidth=1)

        # Plot confidence band (±1 std)
        plt.fill_between(
            recall_grid,
            np.clip(mean_prec - std_prec, 0, 1),
            np.clip(mean_prec + std_prec, 0, 1),
            alpha=0.2,
            color="steelblue",
            label=f"±1 std (n={len(pr_curves)} folds)",
        )

        # Plot mean curve in bold
        plt.plot(
            recall_grid,
            mean_prec,
            color="steelblue",
            linewidth=2.5,
            label=f"Mean PR curve (AUC={pr_auc:.3f})",
        )
    else:
        # Single curve (split eval) - just plot it
        precision, recall = pr_curves[0]
        plt.plot(
            recall,
            precision,
            color="steelblue",
            linewidth=2.5,
            label=f"PR curve (AUC={pr_auc:.3f})",
        )

    plt.xlabel("Recall", fontsize=14)
    plt.ylabel("Precision", fontsize=14)
    plt.title(f"{species_name} - Precision-Recall Curve", fontsize=14)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.legend(loc="best", fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    auc_path = output_dir / f"{species_name}.AUC.iic.png"
    plt.savefig(auc_path, dpi=fig_dpi)
    plt.close()


def ref_scatter(
    u2_vector: np.ndarray,
    u12_vector: np.ndarray,
    species_name: str,
    output_dir: Path,
    fsize: int = 14,
    fig_dpi: int = 300,
):
    """
    Create scatter plot of reference training data.

    Shows U2 and U12 reference introns in 2D score space.

    Args:
        u2_vector: Nx2 array of U2 scores (5' z, BP z)
        u12_vector: Nx2 array of U12 scores (5' z, BP z)
        species_name: Species name for title
        output_dir: Directory to save plot
        fsize: Font size
        fig_dpi: Figure DPI
    """
    plt.figure(figsize=(8, 8))

    plt.scatter(
        *u2_vector[:, :2].T,
        c="xkcd:medium grey",
        alpha=0.5,
        s=42,
        label=f"U2 (n={len(u2_vector)})",
        rasterized=True,
    )

    plt.scatter(
        *u12_vector[:, :2].T,
        c="xkcd:green",
        alpha=0.5,
        s=42,
        label=f"U12 (n={len(u12_vector)})",
        rasterized=True,
    )

    plt.xlabel("5'SS raw score", fontsize=fsize)
    plt.ylabel("BPS raw score", fontsize=fsize)
    plt.title(f"{species_name} - Training Reference Data", fontsize=fsize)

    # Set equal aspect ratio to match original intronIC
    plt.gca().set_aspect("equal", adjustable="box")
    plt.tight_layout()

    output_path = output_dir / f"{species_name}.plot.training_scatter.iic.png"
    plt.savefig(output_path, format="png", dpi=fig_dpi)
    plt.close()


def plot_decision_surface(
    ensemble,
    u12_z: np.ndarray,
    u2_z: np.ndarray,
    output_dir: Path,
    species_name: str,
    hard_neg_threshold: Tuple[float, float] = (2.0, 1.0),
    fig_dpi: int = 150,
):
    """
    Generate 4-panel decision surface visualization for a trained ensemble.

    Panels:
    1. Decision boundary with training data density (hexbin U2, scatter U12/hard neg)
    2. Decision function heatmap with margin contours
    3. Decision function vs 5'z at fixed BPz values
    4. Per-feature contribution breakdown for representative intron profiles

    The plot projects the 6D decision surface onto the 5'z × BPz plane,
    fixing 3'z at the U12 median. Shows the piecewise-linear boundary
    structure created by min/sqrt composite features.

    Args:
        ensemble: Trained SVMEnsemble with models
        u12_z: Nx3 array of U12 z-scores (5'z, BPz, 3'z)
        u2_z: Mx3 array of U2 z-scores (5'z, BPz, 3'z)
        output_dir: Directory for output file
        species_name: Species name for filename
        hard_neg_threshold: (5'z_min, BPz_max) for identifying hard negatives
        fig_dpi: Output resolution
    """
    from matplotlib.colors import TwoSlopeNorm

    hard_mask = (u2_z[:, 0] > hard_neg_threshold[0]) & (u2_z[:, 1] < hard_neg_threshold[1])
    s3_u12_med = float(np.median(u12_z[:, 2]))

    # Compute median values for any extra feature columns (beyond the 3 z-scores)
    n_extra = u12_z.shape[1] - 3
    if n_extra > 0:
        extra_medians = np.median(u12_z[:, 3:], axis=0)
    else:
        extra_medians = np.array([])

    def ensemble_decision(X_3d):
        """Mean decision function across ensemble models."""
        decisions = []
        for m in ensemble.models:
            model = m.model
            if hasattr(model, "calibrated_classifiers_"):
                pipeline = model.calibrated_classifiers_[0].estimator
            else:
                pipeline = model
            decisions.append(pipeline.decision_function(X_3d))
        return np.mean(decisions, axis=0)

    # Grid for surface
    s5_range = np.linspace(-2, 4.5, 200)
    bp_range = np.linspace(-2, 3.5, 200)
    s5_grid, bp_grid = np.meshgrid(s5_range, bp_range)

    X_grid = np.column_stack([
        s5_grid.ravel(), bp_grid.ravel(),
        np.full(s5_grid.size, s3_u12_med)
    ])
    # Pad with median extra-feature values if model uses them
    if n_extra > 0:
        extra_cols = np.tile(extra_medians, (X_grid.shape[0], 1))
        X_grid = np.column_stack([X_grid, extra_cols])
    dec = ensemble_decision(X_grid).reshape(s5_grid.shape)

    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    # Panel 1: Boundary + data density
    ax = axes[0, 0]
    ax.hexbin(u2_z[~hard_mask, 0], u2_z[~hard_mask, 1],
              gridsize=50, cmap="Blues", alpha=0.6, mincnt=1)
    ax.scatter(u2_z[hard_mask, 0], u2_z[hard_mask, 1],
               c="darkorange", s=25, alpha=0.9, marker="x", linewidths=1.5,
               zorder=5, label=f"Hard neg U2 (n={hard_mask.sum()})")
    ax.scatter(u12_z[:, 0], u12_z[:, 1],
               c="red", s=12, alpha=0.5, edgecolors="darkred", linewidths=0.3,
               zorder=6, label=f"U12 (n={len(u12_z)})")
    ax.contour(s5_range, bp_range, dec, levels=[0], colors="black", linewidths=2.5)
    ax.plot([-2, 4.5], [-2, 4.5], "k--", alpha=0.3, linewidth=1)
    ax.set_xlabel("5'z score", fontsize=12)
    ax.set_ylabel("BPz score", fontsize=12)
    ax.set_title("Decision boundary with data density", fontsize=12)
    ax.legend(loc="upper left", fontsize=9)
    ax.set_xlim(-2, 4.5)
    ax.set_ylim(-2, 3.5)

    # Panel 2: Decision function heatmap
    ax = axes[0, 1]
    vmax = 3
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    im = ax.pcolormesh(s5_range, bp_range, np.clip(dec, -vmax, vmax),
                       cmap="RdBu_r", norm=norm, shading="auto")
    ax.contour(s5_range, bp_range, dec, levels=[0], colors="black", linewidths=2.5)
    ax.contour(s5_range, bp_range, dec, levels=[-1, 1],
               colors="gray", linewidths=1, linestyles="dashed")
    ax.plot([-2, 4.5], [-2, 4.5], "k--", alpha=0.3, linewidth=1)
    ax.set_xlabel("5'z score", fontsize=12)
    ax.set_ylabel("BPz score", fontsize=12)
    ax.set_title(f"Decision function heatmap\n(3'z={s3_u12_med:.2f}, dashed=margin \u00b11)", fontsize=12)
    plt.colorbar(im, ax=ax, label="Decision function", shrink=0.8)
    ax.set_xlim(-2, 4.5)
    ax.set_ylim(-2, 3.5)

    # Panel 3: Decision function slices
    ax = axes[1, 0]
    for bp_val, color in [(0.0, "#2166ac"), (1.0, "#4dac26"), (2.0, "#d62728")]:
        X_line = np.column_stack([
            s5_range, np.full(len(s5_range), bp_val),
            np.full(len(s5_range), s3_u12_med)
        ])
        if n_extra > 0:
            X_line = np.column_stack([X_line, np.tile(extra_medians, (len(s5_range), 1))])
        ax.plot(s5_range, ensemble_decision(X_line),
                color=color, linewidth=2, label=f"BPz={bp_val:.0f}")
    ax.axhline(0, color="black", linewidth=1, alpha=0.5)
    ax.set_xlabel("5'z score", fontsize=12)
    ax.set_ylabel("Decision function", fontsize=12)
    ax.set_title("Decision function vs 5'z at fixed BPz", fontsize=12)
    ax.legend(fontsize=10)
    ax.set_xlim(-1, 4.5)

    # Panel 4: Per-feature contributions (linear) or decision function histogram (RBF)
    ax = axes[1, 1]
    # Get pipeline components from first model
    m0 = ensemble.models[0].model
    if hasattr(m0, "calibrated_classifiers_"):
        pipe = m0.calibrated_classifiers_[0].estimator
    else:
        pipe = m0
    transformer = pipe.named_steps["transform"]
    scaler = pipe.named_steps.get("scale")
    svc = pipe.named_steps["svc"]

    if hasattr(svc, "coef_"):
        # Linear kernel: show per-feature contribution breakdown
        feature_names = list(transformer.get_feature_names_out())
        w = svc.coef_[0]

        # Build profiles with extra feature medians appended
        _em = list(extra_medians) if n_extra > 0 else []
        profiles = {
            f"Real U12\n(3.0, 2.5, {s3_u12_med:.1f})": [3.0, 2.5, s3_u12_med] + _em,
            "Hard neg\n(3.0, 0.0, -0.5)": [3.0, 0.0, -0.5] + _em,
            "Marginal\n(2.0, 0.5, 0.0)": [2.0, 0.5, 0.0] + _em,
            "Bulk U2\n(-0.5, -0.5, -0.5)": [-0.5, -0.5, -0.5] + _em,
        }

        x_pos = np.arange(len(profiles))
        bar_width = 0.12
        feat_colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628"]

        for f_idx, (fname, fcolor) in enumerate(zip(feature_names, feat_colors)):
            vals = []
            for raw_3d in profiles.values():
                X_3d = np.array([raw_3d])
                X_t = transformer.transform(X_3d)
                X_s = scaler.transform(X_t) if scaler else X_t
                vals.append(w[f_idx] * X_s[0, f_idx])
            ax.bar(x_pos + f_idx * bar_width, vals, bar_width,
                   label=fname, color=fcolor, alpha=0.8)

        for i, raw_3d in enumerate(profiles.values()):
            total = pipe.decision_function(np.array([raw_3d]))[0]
            ax.plot(i + 2.5 * bar_width, total, "kD", markersize=8, zorder=10)

        ax.axhline(0, color="black", linewidth=1, alpha=0.5)
        ax.set_xticks(x_pos + 2.5 * bar_width)
        ax.set_xticklabels(profiles.keys(), fontsize=9)
        ax.set_ylabel("Weighted contribution (scaled)", fontsize=12)
        ax.set_title("Per-feature contributions (\u25c6 = total decision)", fontsize=12)
        ax.legend(fontsize=8, ncol=2, loc="upper right")
    else:
        # RBF or other non-linear kernel: show decision function histogram
        # Compute decision function values for U12 and U2 reference data
        # Use all columns (z-scores + any extra features)
        u12_3d = u12_z
        u2_3d = u2_z

        u12_dec = ensemble_decision(u12_3d)
        u2_dec = ensemble_decision(u2_3d)

        # Plot overlapping histograms
        bins = np.linspace(min(u2_dec.min(), u12_dec.min()),
                          max(u2_dec.max(), u12_dec.max()), 60)
        ax.hist(u2_dec, bins=bins, alpha=0.6, color="steelblue",
                label=f"U2 (n={len(u2_dec)})", density=True)
        ax.hist(u12_dec, bins=bins, alpha=0.6, color="red",
                label=f"U12 (n={len(u12_dec)})", density=True)
        ax.axvline(0, color="black", linewidth=1.5, linestyle="--", alpha=0.7,
                   label="Decision boundary")
        ax.set_xlabel("Decision function value", fontsize=12)
        ax.set_ylabel("Density", fontsize=12)
        kernel_name = getattr(svc, "kernel", "non-linear")
        n_sv = sum(svc.n_support_) if hasattr(svc, "n_support_") else "?"
        ax.set_title(f"Decision function distribution ({kernel_name}, {n_sv} SVs)", fontsize=12)
        ax.legend(fontsize=10)

    plt.tight_layout()
    output_path = output_dir / f"{species_name}.plot.decision_surface.iic.png"
    plt.savefig(output_path, format="png", dpi=fig_dpi)
    plt.close()
    return output_path
