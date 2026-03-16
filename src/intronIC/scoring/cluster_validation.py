"""
Species-level U12 cluster validation using multi-bandwidth density valley detection.

Assesses whether the SVM's confident U12 calls form a distinct cluster
separated from the U2 distribution by a density valley. Projects all scored
introns onto the U2→U12 discriminating axis (2D: 5'z, BPz), estimates
1D density at multiple KDE bandwidths, and checks for a valley between
the clusters that persists under smoothing.

Two-regime prior framework:
  - Valley detected (median depth > 0.3): U12 cluster confirmed, trust SVM
  - No valley (median depth ≤ 0.3): No distinct U12 cluster, discount prior
    to suppress borderline false positives

Feature space: 2D (5'z, BPz) — the two dimensions where U12 and U2 populations
are most distinct, corresponding to the two U12-specific spliceosomal recognition
events (U11 snRNA ↔ 5'SS, U12 snRNA ↔ BPS). The 3' z-score adds within-U12
variance (AT-AC vs GT-AG subtype variation) that degrades the valley signal.

Validated on 10 species: all confirmed U12 species show median depth ≥ 0.645;
all confirmed U12-absent species show median depth ≤ 0.001.
"""

import numpy as np
from scipy.stats import gaussian_kde
from typing import Optional


def compute_valley_depth(
    u12_points: np.ndarray,
    u2_points: np.ndarray,
    bandwidth_multipliers: tuple = (0.5, 1.0, 1.5, 2.0, 3.0, 5.0),
    n_eval: int = 300,
    max_u2_sample: int = 10000,
    random_seed: int = 42,
) -> dict:
    """
    Multi-bandwidth density valley detection between U12 and U2 clusters.

    Projects all points onto the U2→U12 centroid axis, then estimates 1D
    density at multiple bandwidths. A real valley persists across bandwidths;
    tail artifacts vanish under smoothing.

    Args:
        u12_points: (n_u12, 2) array of [5'z, BPz] for confident U12 calls
        u2_points: (n_u2, 2) array of [5'z, BPz] for U2 introns
        bandwidth_multipliers: Multipliers of Silverman bandwidth to test
        n_eval: Number of points for density evaluation
        max_u2_sample: Max U2 introns to subsample for KDE
        random_seed: Random seed for subsampling

    Returns:
        Dictionary with:
        - median_depth: Median valley depth across bandwidths [0, 1]
        - per_bandwidth_depths: List of depths at each bandwidth
        - has_valley: True if median_depth > 0.3
        - gap_width: Distance between U2 99th pctl and U12 25th pctl on projection axis
        - centroid_sigma: U12 centroid distance from U2 in U2 std units
    """
    if u12_points is None or len(u12_points) < 3:
        return {
            'median_depth': float('nan'),
            'per_bandwidth_depths': [],
            'has_valley': None,
            'gap_width': float('nan'),
            'centroid_sigma': float('nan'),
            'reason': 'insufficient U12 points',
        }

    # Direction vector from U2 centroid to U12 centroid (2D)
    u2_centroid = u2_points.mean(axis=0)
    u12_centroid = u12_points.mean(axis=0)
    direction = u12_centroid - u2_centroid
    norm = np.linalg.norm(direction)
    if norm < 1e-10:
        return {
            'median_depth': 0.0,
            'per_bandwidth_depths': [],
            'has_valley': False,
            'gap_width': 0.0,
            'centroid_sigma': 0.0,
            'reason': 'U12 and U2 centroids coincide',
        }
    direction = direction / norm

    # Project all points onto discriminating axis
    all_points = np.vstack([u2_points, u12_points])
    projections = (all_points - u2_centroid) @ direction
    u2_proj = projections[:len(u2_points)]
    u12_proj = projections[len(u2_points):]

    # Define gap region
    u2_p95 = np.percentile(u2_proj, 95)
    u2_p99 = np.percentile(u2_proj, 99)
    u12_p25 = np.percentile(u12_proj, 25)
    u12_median = np.median(u12_proj)

    # Centroid distance in U2 std units
    u2_std = np.std(u2_proj)
    centroid_sigma = (np.mean(u12_proj) - np.mean(u2_proj)) / u2_std if u2_std > 0 else 0

    gap_width = u12_p25 - u2_p99

    # Check for overlap
    if u12_p25 < u2_p95:
        return {
            'median_depth': 0.0,
            'per_bandwidth_depths': [],
            'has_valley': False,
            'gap_width': gap_width,
            'centroid_sigma': centroid_sigma,
            'reason': 'U12 cluster overlaps with U2 bulk',
        }

    # Subsample U2 for KDE
    np.random.seed(random_seed)
    if len(u2_proj) > max_u2_sample:
        u2_sample = u2_proj[np.random.choice(len(u2_proj), max_u2_sample, replace=False)]
    else:
        u2_sample = u2_proj
    all_sample = np.concatenate([u2_sample, u12_proj])

    # Silverman bandwidth (data-driven)
    silverman_bw = 1.06 * np.std(all_sample) * len(all_sample) ** (-1 / 5)

    # Evaluate density at each bandwidth
    x_eval = np.linspace(u2_p95, u12_median, n_eval)

    depths = []
    for mult in bandwidth_multipliers:
        bw = silverman_bw * mult
        try:
            kde = gaussian_kde(all_sample, bw_method=bw / np.std(all_sample))
            density = kde(x_eval)
        except Exception:
            depths.append(0.0)
            continue

        # Find gap boundaries in evaluation grid
        gap_start = np.searchsorted(x_eval, u2_p99)
        gap_end = np.searchsorted(x_eval, u12_p25)

        if gap_end <= gap_start + 3:
            depths.append(0.0)
            continue

        # Valley depth: proportional drop relative to min of both sides
        density_u2_side = density[gap_start]
        density_u12_side = density[gap_end] if gap_end < len(density) else density[-1]
        gap_min = density[gap_start:gap_end].min()

        reference = min(density_u2_side, density_u12_side)
        if reference > 0:
            depth = max(0.0, 1.0 - gap_min / reference)
        else:
            depth = 0.0

        depths.append(depth)

    median_depth = float(np.median(depths)) if depths else 0.0

    return {
        'median_depth': median_depth,
        'per_bandwidth_depths': depths,
        'has_valley': median_depth > 0.3,
        'gap_width': gap_width,
        'centroid_sigma': centroid_sigma,
    }


def compute_valley_adjusted_prior(
    valley_result: dict,
    empirical_prior: float,
    k: float = 5.0,
) -> float:
    """
    Map valley depth to an adjusted prior for Bayesian probability correction.

    Two-regime framework:
      Valley detected (median depth > 0.3): return 0.5 (no adjustment, trust SVM)
      No valley (median depth ≤ 0.3): return discounted prior via sigmoid

    Args:
        valley_result: Output from compute_valley_depth()
        empirical_prior: Observed U12 fraction in this species
        k: Sigmoid steepness for negative-regime discounting

    Returns:
        Adjusted prior for use with adjust_probabilities_for_prior()
    """
    median_depth = valley_result.get('median_depth', 0.0)

    if np.isnan(median_depth):
        # Too few U12 calls — conservative discount
        return min(empirical_prior, 0.01)

    if valley_result.get('has_valley', False):
        # Valley detected — trust the SVM
        return 0.5

    # No valley — discount the prior
    # Use median_depth as a continuous confidence measure
    scale = 1.0 / (1.0 + np.exp(-k * (median_depth - 0.3)))
    adjusted = empirical_prior * scale

    return max(adjusted, 1e-10)


def validate_u12_cluster(
    five_z_scores: np.ndarray,
    bp_z_scores: np.ndarray,
    svm_scores: np.ndarray,
    type_ids: np.ndarray,
    confidence_threshold: float = 90.0,
) -> dict:
    """
    Full cluster validation pipeline: compute valley depth and adjusted prior.

    Main entry point for species-level U12 cluster validation. Extracts
    confident U12 calls, computes the multi-bandwidth density valley depth,
    and returns an adjusted prior suitable for Bayesian probability correction.

    Args:
        five_z_scores: Array of 5' z-scores for all scored introns
        bp_z_scores: Array of BP z-scores for all scored introns
        svm_scores: Array of SVM scores (0-100) for all scored introns
        type_ids: Array of type assignments ('u12' or 'u2')
        confidence_threshold: SVM score threshold for confident calls (default: 90.0)

    Returns:
        Dictionary with valley depth, regime, adjusted prior, and diagnostics
    """
    # Separate confident U12 and all U2
    u12_mask = (type_ids == 'u12') & (svm_scores >= confidence_threshold)
    u2_mask = type_ids == 'u2'

    n_u12_conf = int(u12_mask.sum())
    n_u12_total = int((type_ids == 'u12').sum())
    n_u2 = int(u2_mask.sum())
    n_total = n_u12_total + n_u2

    empirical_prior = n_u12_total / n_total if n_total > 0 else 0.0

    if n_u12_conf < 3 or n_u2 < 10:
        return {
            'valley_depth': float('nan'),
            'has_valley': None,
            'n_confident_u12': n_u12_conf,
            'n_u12_total': n_u12_total,
            'n_u2': n_u2,
            'empirical_prior': empirical_prior,
            'adjusted_prior': min(empirical_prior, 0.01),
            'regime': 'insufficient',
            'centroid_sigma': float('nan'),
        }

    # Extract 2D points (5'z, BPz)
    u12_points = np.column_stack([five_z_scores[u12_mask], bp_z_scores[u12_mask]])
    u2_points = np.column_stack([five_z_scores[u2_mask], bp_z_scores[u2_mask]])

    # Compute valley depth
    valley_result = compute_valley_depth(u12_points, u2_points)

    # Compute adjusted prior
    adjusted_prior = compute_valley_adjusted_prior(valley_result, empirical_prior)

    regime = 'valley' if valley_result.get('has_valley', False) else 'no_valley'

    return {
        'valley_depth': valley_result['median_depth'],
        'has_valley': valley_result.get('has_valley', False),
        'per_bandwidth_depths': valley_result.get('per_bandwidth_depths', []),
        'n_confident_u12': n_u12_conf,
        'n_u12_total': n_u12_total,
        'n_u2': n_u2,
        'empirical_prior': empirical_prior,
        'adjusted_prior': adjusted_prior,
        'regime': regime,
        'centroid_sigma': valley_result.get('centroid_sigma', float('nan')),
        'gap_width': valley_result.get('gap_width', float('nan')),
    }
