"""
Plotting functions for intronIC visualization.

Recreates the plots from original intronIC:
- Density hexplot of intron scores
- Scatter plot with U12 classifications
- Score histogram
- Training reference plots
- Precision-Recall AUC curves
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Allow to run without X display server
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pathlib import Path
from typing import List, Optional, Tuple

from core.intron import Intron


def plot_classification_results(
    introns: List[Intron],
    output_dir: Path,
    species_name: str,
    threshold: float,
    fig_dpi: int = 300
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
        if intron.scores and intron.scores.five_z_score is not None and intron.scores.bp_z_score is not None:
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
        fig_dpi=fig_dpi
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
        fig_dpi=fig_dpi
    )

    # 3. Score histogram
    svm_scores = [i.scores.svm_score for i in introns if i.scores and i.scores.svm_score is not None]
    hist_path = output_dir / f"{species_name}.plot.score_histogram.iic.png"
    histogram(
        svm_scores,
        threshold=threshold,
        species_name=species_name,
        output_path=hist_path,
        fig_dpi=fig_dpi
    )


def density_hexplot(
    scores: np.ndarray,
    species_name: str,
    output_path: Path,
    xlab: Optional[str] = None,
    ylab: Optional[str] = None,
    fsize: int = 14,
    fig_dpi: int = 300
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

    hx = ax.hexbin(
        *scores.T,
        mincnt=1,
        cmap='inferno',
        bins='log',
        linewidths=0
    )

    # Note: Don't set aspect='equal' - let matplotlib auto-scale based on data ranges
    # Original v1.5.1 used auto aspect for hexplot

    # Clean title: species_name + description + count
    plot_title = f'{species_name} - Motif Score Density (n={len(scores)})'

    if xlab:
        plt.xlabel(xlab, fontsize=fsize)
    if ylab:
        plt.ylabel(ylab, fontsize=fsize)
    plt.title(plot_title, fontsize=fsize)

    # Add colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(hx, cax=cax)
    cb.set_label('Bin density (log10(n))')

    # Save figure
    plt.savefig(output_path, dpi=fig_dpi, bbox_inches='tight')
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
    fig_dpi: int = 300
):
    """
    Create a scatter plot with U12s colored by confidence level and marginal density distributions.

    U12 confidence levels:
    - High (green): score > threshold
    - Medium (orange): (threshold - stdev) < score <= threshold
    - Low (red): score <= (threshold - stdev)
    - U2 (grey): classified as U2

    Marginal density distributions are shown on the top and right axes.

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
    # Create figure with GridSpec for marginal distributions
    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(
        3, 3,
        figure=fig,
        width_ratios=[1, 4, 0.8],  # Increased right margin width for count labels
        height_ratios=[1, 4, 0.1],
        hspace=0.02,
        wspace=0.02,
        top=0.95,  # Leave space for suptitle
        bottom=0.05,
        left=0.05,
        right=0.95
    )

    # Main scatter plot (center)
    ax_main = fig.add_subplot(gs[1, 1])
    # Top marginal (x distribution)
    ax_top = fig.add_subplot(gs[0, 1], sharex=ax_main)
    # Right marginal (y distribution)
    ax_right = fig.add_subplot(gs[1, 2], sharey=ax_main)

    # Calculate threshold boundaries
    svm_scores = [i.scores.svm_score for i in introns if i.scores and i.scores.svm_score is not None]
    score_stdev = np.std(svm_scores)
    high_val = threshold
    med_val = threshold - score_stdev

    # Assign colors based on classification and confidence
    cluster_colors = []
    u2_count, u12_low, u12_med, u12_high = 0, 0, 0, 0

    for intron in introns:
        if not intron.metadata or not intron.scores:
            continue

        itype = intron.metadata.type_id
        score = intron.scores.svm_score

        if itype == 'u2':
            u2_count += 1
            color = 'xkcd:medium grey'
        elif score > high_val:
            u12_high += 1
            color = 'xkcd:green'
        elif med_val < score <= high_val:
            u12_med += 1
            color = 'xkcd:orange'
        else:  # score <= med_val
            u12_low += 1
            color = 'xkcd:red'

        cluster_colors.append(color)

    # Create legend
    legend_colors = ['xkcd:medium grey', 'xkcd:red', 'xkcd:orange', 'xkcd:green']
    legend_labels = [
        'U2',
        f'U12<={int(med_val)}',
        f'{int(med_val)}<U12<={int(high_val)}',
        f'U12>{int(high_val)}'
    ]
    legend_counts = [u2_count, u12_low, u12_med, u12_high]

    legend_patches = []
    for label, count, color in zip(legend_labels, legend_counts, legend_colors):
        label_with_count = f'{label} ({count})'
        patch = mpatches.Patch(color=color, label=label_with_count)
        legend_patches.append(patch)

    # Plot main scatter
    ax_main.scatter(
        *scores[:, :2].T,
        s=20,
        c=cluster_colors,
        alpha=0.5,
        rasterized=True
    )

    ax_main.legend(handles=legend_patches, fontsize=fsize-2)
    ax_main.set_xlabel(xlab, fontsize=fsize)
    ax_main.set_ylabel(ylab, fontsize=fsize)

    # Note: Don't set aspect='equal' here - z-scores can have different ranges
    # (BP variance is often 3× higher than 5' variance in training data)
    # Original v1.5.1 used auto aspect for classification plots

    # Plot marginal distributions
    # Top: X distribution (5' z-score)
    ax_top.hist(scores[:, 0], bins=50, color='steelblue', alpha=0.7, edgecolor='none')
    ax_top.set_ylabel('Count', fontsize=fsize-2)
    ax_top.tick_params(labelbottom=False)
    ax_top.spines['top'].set_visible(False)
    ax_top.spines['right'].set_visible(False)

    # Right: Y distribution (BP z-score)
    ax_right.hist(scores[:, 1], bins=50, orientation='horizontal', color='steelblue', alpha=0.7, edgecolor='none')
    ax_right.set_xlabel('Count', fontsize=fsize-2)
    ax_right.tick_params(labelleft=False, labelrotation=45)  # Rotate labels to prevent overlap
    ax_right.spines['top'].set_visible(False)
    ax_right.spines['right'].set_visible(False)

    # Add figure-level title above marginal distributions
    # Clean title: species_name + description
    fig.suptitle(f'{species_name} - U12 Classification Results', fontsize=fsize+2, y=0.98, weight='bold')

    # Save figure
    plt.savefig(output_path, dpi=fig_dpi, bbox_inches='tight')
    plt.close()


def histogram(
    data_list: List[float],
    threshold: float,
    species_name: str,
    output_path: Path,
    grid: bool = True,
    bins: int = 100,
    log: bool = True,
    fig_dpi: int = 300
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
        plt.yscale('log')

    plt.hist(data_list, bins=bins)

    if grid:
        plt.grid(True, which="both", ls="--", alpha=0.7)

    # Clean title: species_name + description
    plt.title(f'{species_name} - U12 Score Distribution', fontsize=14)

    plt.xlabel('U12 score', fontsize=14)
    plt.ylabel('Number of introns', fontsize=14)

    # Add threshold line
    plt.axvline(
        threshold,
        color='orange',
        linestyle='--',
        label=f'U12 threshold: {threshold}'
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
    fig_dpi: int = 300
):
    """
    Generate training reference plots.

    Creates three plots:
    1. Reference scatter plot (U2 vs U12 training data)
    2. Reference hexplot (density of reference data)
    3. Precision-Recall AUC curve

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
        fig_dpi=fig_dpi
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
        fig_dpi=fig_dpi
    )

    # 3. Precision-Recall curve
    # Plot all curves (matches original intronIC behavior for multiple CV folds)
    plt.figure(figsize=(8, 8))
    for precision, recall in pr_curves:
        plt.plot(recall, precision)
    plt.xlabel('Recall', fontsize=14)
    plt.ylabel('Precision', fontsize=14)
    plt.title(f'{species_name} - Precision-Recall AUC: {pr_auc:.3f}', fontsize=14)
    plt.tight_layout()

    auc_path = output_dir / f'{species_name}.AUC.iic.png'
    plt.savefig(auc_path, dpi=fig_dpi)
    plt.close()


def ref_scatter(
    u2_vector: np.ndarray,
    u12_vector: np.ndarray,
    species_name: str,
    output_dir: Path,
    fsize: int = 14,
    fig_dpi: int = 300
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
        c='xkcd:medium grey',
        alpha=0.5,
        s=42,
        label=f'U2 (n={len(u2_vector)})',
        rasterized=True
    )

    plt.scatter(
        *u12_vector[:, :2].T,
        c='xkcd:green',
        alpha=0.5,
        s=42,
        label=f'U12 (n={len(u12_vector)})',
        rasterized=True
    )

    plt.xlabel("5' z-score", fontsize=fsize)
    plt.ylabel('BPS z-score', fontsize=fsize)
    plt.title(f'{species_name} - Training Reference Data', fontsize=fsize)

    # Set equal aspect ratio to match original intronIC
    plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout()

    output_path = output_dir / f'{species_name}.plot.training_scatter.iic.png'
    plt.savefig(output_path, format='png', dpi=fig_dpi)
    plt.close()
