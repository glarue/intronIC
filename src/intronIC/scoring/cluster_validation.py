"""
Species-level U12 cluster validation using kNN-median silhouette score.

Computes whether the SVM's confident U12 calls form a coherent cluster
separate from the U2 distribution. This provides a species-level confidence
measure that can be used to adjust the classification prior.

Two-regime prior framework:
  - Positive silhouette (> 0): U12 cluster detected, trust SVM as-is
  - Negative silhouette (≤ 0): No confident U12 cluster, discount prior
    to suppress borderline false positives

The kNN-median silhouette is a robust variant of the standard silhouette
coefficient that uses k-nearest-neighbor distances (robust to outliers
with degraded motifs) and median aggregation (robust to individual
outlier introns).

Validated on 30 species spanning vertebrates, insects, plants, fungi,
protists, and known U12-absent lineages. All confirmed U12 species
score positive; all confirmed U12-absent species score negative.
"""

import numpy as np
from scipy.spatial import cKDTree
from typing import Tuple, Optional


def compute_knn_median_silhouette(
    u12_points: np.ndarray,
    u2_points: np.ndarray,
    k: int = 5,
) -> float:
    """
    Compute the kNN-median silhouette score for U12 calls.

    For each U12 point:
      a_i = distance to k-th nearest U12 neighbor (intra-cluster)
      b_i = distance to k-th nearest U2 neighbor (inter-cluster)
      s_i = (b_i - a_i) / max(a_i, b_i)

    Returns the median of per-point silhouette values.

    Args:
        u12_points: (n_u12, 2) array of [5'z, BPz] for confident U12 calls
        u2_points: (n_u2, 2) array of [5'z, BPz] for U2 introns
        k: Number of nearest neighbors (default: 5)

    Returns:
        Median silhouette score in [-1, 1].
        Positive = U12 calls form a distinct cluster.
        Negative = U12 calls are dispersed in the U2 tail.
        NaN if fewer than 3 U12 points.
    """
    if u12_points is None or len(u12_points) < 3:
        return float('nan')

    k_use = min(k, len(u12_points) - 1)
    if k_use < 1:
        return float('nan')

    # Inter-cluster: distance from each U12 to k-th nearest U2
    u2_tree = cKDTree(u2_points)
    inter_dists, _ = u2_tree.query(u12_points, k=k_use)
    b_vals = inter_dists[:, -1]  # k-th nearest

    # Intra-cluster: distance from each U12 to k-th nearest other U12
    u12_tree = cKDTree(u12_points)
    intra_dists, _ = u12_tree.query(u12_points, k=min(k_use + 1, len(u12_points)))
    a_vals = intra_dists[:, -1]  # k-th nearest (query includes self, so -1 is k-th neighbor)

    # Per-point silhouette
    per_point = (b_vals - a_vals) / np.maximum(a_vals, b_vals)

    return float(np.median(per_point))


def compute_silhouette_adjusted_prior(
    silhouette: float,
    empirical_prior: float,
    k: float = 5.0,
) -> float:
    """
    Map silhouette score to an adjusted prior for Bayesian probability correction.

    Two-regime framework:
      silhouette > 0:  return 0.5 (no adjustment, trust SVM calibration)
      silhouette ≤ 0:  return empirical_prior × sigmoid(k × silhouette)

    For positive silhouette, the SVM's own calibration is appropriate —
    there IS a real U12 cluster and the decision boundary is meaningful.
    A prior of 0.5 means no Bayesian adjustment (the SVM operates as trained).

    For negative silhouette, the SVM is likely calling U2 tail introns.
    The sigmoid-scaled empirical prior requires very high SVM confidence
    (>99.9%) for a call to survive, effectively suppressing borderline FPs.

    Args:
        silhouette: kNN-median silhouette score (from compute_knn_median_silhouette)
        empirical_prior: Observed U12 fraction in this species (n_u12 / n_total)
        k: Sigmoid steepness parameter (default: 5.0).
           Higher k = sharper transition at silhouette=0.

    Returns:
        Adjusted prior for use with adjust_probabilities_for_prior().
        Range: [~0, 0.5]
    """
    if np.isnan(silhouette):
        # Too few U12 calls to compute silhouette — use empirical prior
        # with a conservative discount
        return min(empirical_prior, 0.01)

    if silhouette > 0:
        # Positive silhouette: U12 cluster detected, trust SVM
        return 0.5

    # Negative silhouette: discount the prior via sigmoid
    scale = 1.0 / (1.0 + np.exp(-k * silhouette))
    adjusted = empirical_prior * scale

    # Floor to prevent exactly zero (which breaks Bayes' rule)
    return max(adjusted, 1e-10)


def validate_u12_cluster(
    five_z_scores: np.ndarray,
    bp_z_scores: np.ndarray,
    svm_scores: np.ndarray,
    type_ids: np.ndarray,
    confidence_threshold: float = 90.0,
    k: int = 5,
    sigmoid_k: float = 5.0,
) -> dict:
    """
    Full cluster validation pipeline: compute silhouette and adjusted prior.

    This is the main entry point for species-level U12 cluster validation.
    It extracts confident U12 calls, computes the kNN-median silhouette,
    and returns an adjusted prior suitable for Bayesian probability correction.

    Args:
        five_z_scores: Array of 5' z-scores for all scored introns
        bp_z_scores: Array of BP z-scores for all scored introns
        svm_scores: Array of SVM scores (0-100) for all scored introns
        type_ids: Array of type assignments ('u12' or 'u2') for all scored introns
        confidence_threshold: SVM score threshold for "confident" calls (default: 90.0)
        k: kNN parameter for silhouette (default: 5)
        sigmoid_k: Steepness for silhouette-to-prior mapping (default: 5.0)

    Returns:
        Dictionary with:
        - silhouette: kNN-median silhouette score
        - n_confident_u12: Number of confident U12 calls
        - n_u2: Number of U2 calls
        - empirical_prior: Observed U12 fraction
        - adjusted_prior: Prior for Bayesian adjustment
        - regime: 'positive' or 'negative' or 'insufficient'
        - u12_5z_mean: Mean 5' z-score of confident U12 calls
        - u12_bpz_mean: Mean BP z-score of confident U12 calls
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
            'silhouette': float('nan'),
            'n_confident_u12': n_u12_conf,
            'n_u12_total': n_u12_total,
            'n_u2': n_u2,
            'empirical_prior': empirical_prior,
            'adjusted_prior': min(empirical_prior, 0.01),
            'regime': 'insufficient',
            'u12_5z_mean': float('nan'),
            'u12_bpz_mean': float('nan'),
        }

    # Extract 2D points
    u12_points = np.column_stack([five_z_scores[u12_mask], bp_z_scores[u12_mask]])
    u2_points = np.column_stack([five_z_scores[u2_mask], bp_z_scores[u2_mask]])

    # Compute silhouette
    silhouette = compute_knn_median_silhouette(u12_points, u2_points, k=k)

    # Compute adjusted prior
    adjusted_prior = compute_silhouette_adjusted_prior(
        silhouette, empirical_prior, k=sigmoid_k
    )

    regime = 'positive' if silhouette > 0 else 'negative'

    return {
        'silhouette': silhouette,
        'n_confident_u12': n_u12_conf,
        'n_u12_total': n_u12_total,
        'n_u2': n_u2,
        'empirical_prior': empirical_prior,
        'adjusted_prior': adjusted_prior,
        'regime': regime,
        'u12_5z_mean': float(u12_points[:, 0].mean()),
        'u12_bpz_mean': float(u12_points[:, 1].mean()),
    }
