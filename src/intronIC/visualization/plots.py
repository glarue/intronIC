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


def plot_classification_results_from_file(
    score_file: Path,
    output_dir: Path,
    species_name: str,
    threshold: float,
    fig_dpi: int = 300,
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
    # Read scores from file
    five_z_scores = []
    bp_z_scores = []
    svm_scores = []

    with open(score_file, "r") as f:
        header = f.readline().strip().split("\t")

        # Find column indices
        try:
            five_z_idx = header.index("5'_z")
            bp_z_idx = header.index("bp_z")
            svm_idx = header.index("svm_score")
        except ValueError as e:
            raise ValueError(f"Missing required column in score file: {e}")

        for line in f:
            fields = line.strip().split("\t")
            if len(fields) <= max(five_z_idx, bp_z_idx, svm_idx):
                continue

            # Parse z-scores (may be "NA")
            five_z = fields[five_z_idx]
            bp_z = fields[bp_z_idx]
            svm = fields[svm_idx]

            if five_z != "NA" and bp_z != "NA":
                try:
                    five_z_scores.append(float(five_z))
                    bp_z_scores.append(float(bp_z))
                except ValueError:
                    continue

            if svm != "NA":
                try:
                    svm_scores.append(float(svm))
                except ValueError:
                    continue

    if not five_z_scores:
        # No valid scores to plot
        return

    score_vector = np.array(list(zip(five_z_scores, bp_z_scores)))

    # 1. Density hexplot
    hexplot_path = output_dir / f"{species_name}.plot.hex.iic.png"
    density_hexplot(
        score_vector,
        species_name=species_name,
        output_path=hexplot_path,
        xlab="5' z-score",
        ylab="BPS z-score",
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
        xlab="5' z-score",
        ylab="BPS z-score",
        threshold=threshold,
        fig_dpi=fig_dpi,
    )

    # 3. Score histogram
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
    """
    # Extract score vectors (5' z-score, BP z-score)
    score_vector = []
    for intron in introns:
        if (
            intron.scores
            and intron.scores.five_z_score is not None
            and intron.scores.bp_z_score is not None
        ):
            score_vector.append([intron.scores.five_z_score, intron.scores.bp_z_score])

    score_vector = np.array(score_vector)

    # 1. Density hexplot
    hexplot_path = output_dir / f"{species_name}.plot.hex.iic.png"
    density_hexplot(
        score_vector,
        species_name=species_name,
        output_path=hexplot_path,
        xlab="5' z-score",
        ylab="BPS z-score",
        fig_dpi=fig_dpi,
    )

    # 2. Scatter plot with U12 classification
    scatter_path = output_dir / f"{species_name}.plot.scatter.iic.png"
    scatter_plot(
        introns,
        score_vector,
        species_name=species_name,
        output_path=scatter_path,
        xlab="5' z-score",
        ylab="BPS z-score",
        threshold=threshold,
        fig_dpi=fig_dpi,
    )

    # 3. Score histogram
    svm_scores = [
        i.scores.svm_score
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

    # Assign colors based on classification and confidence
    cluster_colors = []
    u2_count, u12_low, u12_med, u12_high = 0, 0, 0, 0

    for i, score in enumerate(svm_scores):
        if i >= len(score_vector):
            break

        # Determine if U2 or U12 based on type_ids if provided, else use score
        if type_ids is not None and i < len(type_ids):
            is_u2 = type_ids[i] == "u2"
        else:
            is_u2 = score < 50  # Raw classifier threshold

        if is_u2:
            u2_count += 1
            color = "xkcd:medium grey"
        elif score > high_val:
            u12_high += 1
            color = "xkcd:green"
        elif med_val < score <= high_val:
            u12_med += 1
            color = "xkcd:orange"
        else:
            u12_low += 1
            color = "xkcd:red"

        cluster_colors.append(color)

    # Trim score_vector to match colors if needed
    n_points = len(cluster_colors)
    plot_scores = score_vector[:n_points]

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

    # Plot main scatter
    ax_main.scatter(
        *plot_scores[:, :2].T, s=20, c=cluster_colors, alpha=0.5, rasterized=True
    )

    ax_main.legend(handles=legend_patches, fontsize=fsize - 2)
    ax_main.set_xlabel(xlab, fontsize=fsize)
    ax_main.set_ylabel(ylab, fontsize=fsize)

    # First calculate symmetric limits BEFORE plotting anything
    x_data = plot_scores[:, 0]
    y_data = plot_scores[:, 1]
    x_range = x_data.max() - x_data.min()
    y_range = y_data.max() - y_data.min()
    max_range = max(x_range, y_range)

    # Center on data and extend by max range
    x_center = (x_data.max() + x_data.min()) / 2
    y_center = (y_data.max() + y_data.min()) / 2
    margin = max_range * 0.05  # 5% margin

    # Calculate the symmetric limits
    xlim = (x_center - max_range / 2 - margin, x_center + max_range / 2 + margin)
    ylim = (y_center - max_range / 2 - margin, y_center + max_range / 2 + margin)

    # Set limits on main plot first
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
        xlab="5' z-score",
        ylab="BPS z-score",
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

    plt.xlabel("5' z-score", fontsize=fsize)
    plt.ylabel("BPS z-score", fontsize=fsize)
    plt.title(f"{species_name} - Training Reference Data", fontsize=fsize)

    # Set equal aspect ratio to match original intronIC
    plt.gca().set_aspect("equal", adjustable="box")
    plt.tight_layout()

    output_path = output_dir / f"{species_name}.plot.training_scatter.iic.png"
    plt.savefig(output_path, format="png", dpi=fig_dpi)
    plt.close()
