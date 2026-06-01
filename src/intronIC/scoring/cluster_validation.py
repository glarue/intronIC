"""
Species-level U12 cluster validation using multi-bandwidth density valley detection.

Assesses whether the SVM's confident U12 calls form a distinct cluster
separated from the U2 distribution by a density valley. Projects all scored
introns onto the U2→U12 discriminating axis (3D: 5'z, BPz, 3'z), estimates
1D density along that axis at multiple KDE bandwidths, and checks for a
valley between the clusters that persists under smoothing.

Separation statistic (v2.7 gate + π_species decision input)
-----------------------------------------------------------
The decision statistic is `gap_fraction` = gap_width / Δmean on the Fisher
axis (dimensionless, clade-invariant), with a deterministic bootstrap SD as
its confidence. It replaced the KDE `median_depth` ("valley depth"), which
was seed-fragile and spuriously failed real low-N U12 bearers. `median_depth`
is still computed and surfaced as a DIAGNOSTIC, but no longer drives routing
or score adjustment.

Prior framework (species-level π_species, confidence-shrunk):
  - Clear separation (high gap_fraction): π_species → π_train, trust the SVM
  - Poor separation, CONFIDENT (low gap_fraction, tight bootstrap): suppress
  - Poor separation, UNCERTAIN (wide bootstrap, typically low-N): shrink back
    toward no-suppression — don't over-kill a possibly-real low-N bearer

Feature space and projection direction
--------------------------------------
The discriminating direction is computed via Fisher's linear discriminant:

    w = Σ⁻¹ · (μ_U12 − μ_U2)

where Σ is the diagonal-shrunk pooled within-class covariance. This is the
optimal 1D projection under Gaussian assumptions: instead of just pointing
from one centroid to the other (the naive choice), it down-weights features
with high within-class variance — which is what makes adding 3'z safe here.

Including 3'z with the *naive* centroid direction degrades the valley signal,
because 3'z carries substantial within-U12 variance (driven by AT-AC vs GT-AG
subtype differences). Fisher's reweighting absorbs that variance instead of
inheriting it into the projection. Empirically, switching to 3D Fisher gains
+2–10% valley depth on real-U12 species (largest where the cluster is
diffuse — Arabidopsis, vertebrates) without inflating depth on any
U12-absent species.

Validated on 16 species: all confirmed U12 species show median depth ≥ 0.57;
all confirmed U12-absent species show median depth ≤ 0.005, well below the
0.3 valley threshold.
"""

import csv
import math
import tempfile
import shutil
import numpy as np
from dataclasses import replace
from pathlib import Path
from scipy.stats import gaussian_kde
from typing import Optional, List, TYPE_CHECKING

if TYPE_CHECKING:
    from intronIC.core.intron import Intron


# gap_fraction → π_species mapping (the species-level prior). midpoint/width bracket the separation
# band where real-U12 and U12-absent species' UCLs begin to interleave. The UCL quantile sets how
# conservative the low-N guard is: π_species keys on the q-th bootstrap percentile of gap_fraction (an
# OPTIMISTIC estimate), so higher q ⇒ more benefit-of-the-doubt ⇒ we only suppress when CONFIDENT the
# separation is poor.
#
# STEP-6 CALIBRATED (2026-05-31; conservation_corpus/scripts/pi_species_refine.py) on the 14-species
# panel + green-algae anchors, net-quality objective. Key empirical finding: real and absent species'
# UCLs INTERLEAVE at the edges (clean chlorophyte-absents Coccomyxa/Chlorella 0.27-0.29 < all reals,
# but TetThe-absent 0.46 sits inside the vertebrate real range 0.45-0.51, and real Klebs 0.37 < real
# AmbTri/OrySat 0.40-0.44). So π_species CANNOT be aggressive — a high midpoint over-kills moderate-UCL
# reals (esp. diverged plants) to chase TetThe. midpoint=0.35 is the knee: captures the safely-available
# FP reduction (chlorophytes well below it) without reaching into the real range; below it FP is flat
# (the per-intron discount, not π_species, is the FP lever there); 0.42 cut FP 21→12 but jumped hard-TP
# loss 12→17 (all diverged-plant AmbTri). width=0.10 (sharper) strictly dominated 0.14 (lower FP, ≥ TP
# retention per-species, 0 new hard loss). q: UCLs barely move over q 0.7-0.9 (bootstrap converged) → 0.8.
# The motif-calibrated absent residual (TetThe etc.) is irreducible by any species-internal prior →
# offline-snRNA clean-negative case. See [[tp-loss-vs-fp-tradeoff-preference]] / calibration_plan §6b.
DEFAULT_GAP_FRACTION_MIDPOINT = 0.35
DEFAULT_GAP_FRACTION_WIDTH = 0.10
DEFAULT_GAP_FRACTION_UCL_QUANTILE = 0.8


def _fishers_discriminant_direction(
    u12_points: np.ndarray,
    u2_points: np.ndarray,
    shrinkage: float = 0.05,
) -> np.ndarray:
    """Fisher's linear discriminant direction Σ⁻¹(μ_U12 − μ_U2), normalized.

    Σ is the pooled within-class covariance with diagonal shrinkage:

        Σ = (1 − α) · pooled_cov + α · diag(pooled_cov)

    Shrinkage stabilizes when the U12 cluster is small or features are
    highly correlated. A small ridge term is also added to guarantee
    invertibility under all input conditions.

    Returns a unit vector with the same dimensionality as the input
    points; if the centroids coincide (degenerate case) returns a zero
    vector.
    """
    n12, n2 = len(u12_points), len(u2_points)
    mu12 = u12_points.mean(axis=0)
    mu2 = u2_points.mean(axis=0)

    if n12 > 1:
        cov12 = np.cov(u12_points, rowvar=False)
    else:
        cov12 = np.zeros((u12_points.shape[1], u12_points.shape[1]))
    cov2 = np.cov(u2_points, rowvar=False)
    pooled = ((n12 - 1) * cov12 + (n2 - 1) * cov2) / max(n12 + n2 - 2, 1)

    diag = np.diag(np.diag(pooled))
    Sigma = (1 - shrinkage) * pooled + shrinkage * diag
    Sigma = Sigma + 1e-6 * np.eye(Sigma.shape[0])

    diff = mu12 - mu2
    direction = np.linalg.solve(Sigma, diff)
    norm = np.linalg.norm(direction)
    if norm < 1e-10:
        return np.zeros_like(diff)
    return direction / norm


def compute_valley_depth(
    u12_points: np.ndarray,
    u2_points: np.ndarray,
    bandwidth_multipliers: tuple = (0.5, 1.0, 1.5, 2.0, 3.0, 5.0),
    n_eval: int = 300,
    max_u2_sample: int = 10000,
    random_seed: int = 42,
    ucl_quantile: float = DEFAULT_GAP_FRACTION_UCL_QUANTILE,
) -> dict:
    """
    Multi-bandwidth density valley detection between U12 and U2 clusters.

    Projects all points onto Fisher's linear discriminant direction
    Σ⁻¹(μ_U12 − μ_U2), then estimates 1D density at multiple bandwidths.
    A real valley persists across bandwidths; tail artifacts vanish under
    smoothing.

    Args:
        u12_points: (n_u12, D) array of features for confident U12 calls.
            Production input is D=3 with columns (5'z, BPz, 3'z).
        u2_points: (n_u2, D) array of the same features for U2 introns.
        bandwidth_multipliers: Multipliers of Silverman bandwidth to test
        n_eval: Number of points for density evaluation
        max_u2_sample: Max U2 introns to subsample for KDE
        random_seed: Random seed for subsampling

    Returns:
        Dictionary with:
        - gap_fraction: gap_width / Δmean on the Fisher axis — the dimensionless,
          clade-invariant separation statistic that drives the gate and π_species.
          Finite > 0 = measurable separation (graded); finite ≤ 0 = overlap; -inf =
          structural no-separation (coincident/inverted centroids); NaN = cannot
          assess (too few confident U12 calls).
        - gap_fraction_ucl: deterministic bootstrap UPPER confidence limit (q-th
          percentile) of gap_fraction — the OPTIMISTIC estimate π_species keys on.
          Wide bootstrap (low-N/unstable) ⇒ UCL well above the point estimate ⇒ less
          suppression (the over-kill guard); -inf when most resamples show no separation.
        - median_depth: Median KDE valley depth across bandwidths [0, 1] — now a
          DIAGNOSTIC only (no longer a decision input; superseded by gap_fraction)
        - per_bandwidth_depths: List of depths at each bandwidth (diagnostic)
        - has_valley: True if median_depth > 0.3 (diagnostic)
        - gap_width: Distance between U2 99th pctl and U12 25th pctl on projection axis
        - centroid_sigma: U12 centroid distance from U2 in U2 std units
        - u2_std: Std of the U2 projection (scale of the Fisher axis)
    """
    if u12_points is None or len(u12_points) < 3:
        return {
            'median_depth': float('nan'),
            'per_bandwidth_depths': [],
            'has_valley': None,
            'gap_width': float('nan'),
            'centroid_sigma': float('nan'),
            # gap_fraction = NaN ⇒ "cannot assess separation" ⇒ downstream π_species
            # leaves the score unchanged (conservative no-op), matching the historical
            # NaN-valley_depth skip for species with too few confident U12 calls.
            'gap_fraction': float('nan'),
            'gap_fraction_ucl': float('nan'),
            'u2_std': float('nan'),
            'reason': 'insufficient U12 points',
        }

    # Determinism: canonically sort both point sets at entry. np.mean/np.cov/np.percentile
    # accumulate in array order, so float non-associativity makes every downstream reduction
    # (Fisher direction → projections → gap_fraction/centroid_sigma) depend on the caller's row
    # order at the ULP level — which differs between the streaming path (imap_unordered over contigs)
    # and the in-memory path. Lexsorting the rows fixes the accumulation order, so the result is
    # bit-identical regardless of input order (the values are unchanged; only their order is
    # canonicalized — every statistic here is mathematically permutation-invariant).
    u12_points = u12_points[np.lexsort(u12_points.T)]
    u2_points = u2_points[np.lexsort(u2_points.T)]

    # Fisher's linear discriminant: Σ⁻¹(μ_U12 − μ_U2), shrunk and normalized.
    # Down-weights features with large within-class variance (e.g., 3'z, where
    # AT-AC / GT-AG subtype differences inflate the U12 spread) so they don't
    # contaminate the projection direction.
    u2_centroid = u2_points.mean(axis=0)
    u12_centroid = u12_points.mean(axis=0)
    direction = _fishers_discriminant_direction(u12_points, u2_points)
    if np.linalg.norm(direction) < 1e-10:
        return {
            'median_depth': 0.0,
            'per_bandwidth_depths': [],
            'has_valley': False,
            'gap_width': 0.0,
            'centroid_sigma': 0.0,
            # Centroids coincide ⇒ no discriminating direction ⇒ structurally no U12
            # cluster. gap_fraction = -inf flags confident suppression in π_species
            # (independent of bootstrap spread).
            'gap_fraction': float('-inf'),
            'gap_fraction_ucl': float('nan'),
            'u2_std': 0.0,
            'reason': 'U12 and U2 centroids coincide',
        }

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

    # gap_fraction = gap_width / Δmean on the Fisher axis — the dimensionless, clade-invariant
    # separation statistic that replaces median_depth as the gate / π_species decision input.
    # RNG-free; bit-identical across input orders thanks to the lexsort at entry (without it, float
    # non-associativity of the mean/percentile reductions leaks the caller's row order into the last
    # ULP, breaking streaming-vs-in-memory bit-identity).
    delta_mu = float(np.mean(u12_proj) - np.mean(u2_proj))
    gap_fraction = float(gap_width / delta_mu) if delta_mu > 0 else float('-inf')

    # gap_fraction confidence — deterministic bootstrap UPPER confidence limit (UCL). π_species keys
    # on THIS (not the point estimate): the q-th percentile (q = ucl_quantile) of the bootstrap
    # gap_fraction distribution — an OPTIMISTIC estimate. When the bootstrap is tight (well-sampled)
    # the UCL ≈ the point estimate; when it is wide (low-N / unstable) the UCL sits well above the
    # point estimate ⇒ less suppression ⇒ we suppress only when CONFIDENT the separation is poor (the
    # low-N over-kill guard, and it rescues a real species whose point estimate dips ≤ 0 by sampling
    # noise). Resamples with Δmean ≤ 0 (centroids inverted ⇒ no separation) form a point mass at the
    # BOTTOM of the distribution; if they exceed the (1−q) tail the UCL is -inf (confident
    # no-separation). The UCL is bounded/percentile-based, so unlike σ/μ it stays well-behaved when
    # gap_fraction crosses zero or Δmean → 0. LOCAL np.random.default_rng(random_seed) over SORTED
    # points + B=500 ⇒ bit-identical across runs and streaming-vs-in-memory. Computed before the
    # gap-based early returns so the overlap path carries it too.
    gap_fraction_ucl = float('nan')
    if len(u12_proj) >= 3:
        _rng = np.random.default_rng(random_seed)
        _u12s = np.sort(u12_proj)
        _n, _B = len(_u12s), 500
        _bs = _u12s[_rng.integers(0, _n, size=(_B, _n))]
        _dmu = _bs.mean(axis=1) - np.mean(u2_proj)
        _pos = _dmu > 0
        _n_pos = int(_pos.sum())
        _frac_bad = 1.0 - _n_pos / _B            # mass of "no separation" (Δmean ≤ 0) resamples
        if _n_pos == 0 or _frac_bad >= ucl_quantile:
            # the q-th percentile falls inside the no-separation point mass
            gap_fraction_ucl = float('-inf')
        else:
            _ratio = (np.percentile(_bs[_pos], 25, axis=1) - u2_p99) / _dmu[_pos]
            # remap q into the positive-Δmean sub-distribution (point mass excluded from the bottom)
            _q_adj = (ucl_quantile - _frac_bad) / (_n_pos / _B)
            gap_fraction_ucl = float(np.percentile(_ratio, 100.0 * _q_adj))

    # Check for overlap
    if u12_p25 < u2_p95:
        return {
            'median_depth': 0.0,
            'per_bandwidth_depths': [],
            'has_valley': False,
            'gap_width': gap_width,
            'centroid_sigma': centroid_sigma,
            # Overlap ⇒ gap_fraction is finite-negative (u12_p25 < u2_p95 < u2_p99); the
            # graded π_species suppresses it iff the bootstrap UCL is also low (confident).
            'gap_fraction': gap_fraction,
            'gap_fraction_ucl': gap_fraction_ucl,
            'u2_std': float(u2_std),
            'reason': 'U12 cluster overlaps with U2 bulk',
        }

    # Subsample U2 for KDE.
    #
    # Sort the projection before sampling so that the random indices produce
    # the same subsample regardless of the upstream iteration order. Without
    # this, the streaming path (parallel imap_unordered over contigs) and the
    # in-memory path (sequential per-contig) feed introns in different orders,
    # so np.random.choice picks different points and the KDE valley depth
    # diverges by a few percent — which then leaks into adjusted_score and
    # rel_score for thousands of introns in score_info.iic.
    np.random.seed(random_seed)
    u2_proj_sorted = np.sort(u2_proj)
    if len(u2_proj_sorted) > max_u2_sample:
        u2_sample = u2_proj_sorted[
            np.random.choice(len(u2_proj_sorted), max_u2_sample, replace=False)
        ]
    else:
        u2_sample = u2_proj_sorted
    all_sample = np.concatenate([u2_sample, np.sort(u12_proj)])

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
        'gap_fraction': gap_fraction,
        'gap_fraction_ucl': gap_fraction_ucl,
        'u2_std': float(u2_std),
    }


def _confidence_shrunk_pi_species(
    gap_fraction_ucl: Optional[float],
    *,
    gap_fraction_midpoint: float = DEFAULT_GAP_FRACTION_MIDPOINT,
    gap_fraction_width: float = DEFAULT_GAP_FRACTION_WIDTH,
    prior_floor: float = 0.001,
    pi_train: float = 0.5,
):
    """Map the bootstrap-UCL of gap_fraction to a species prior π_species.

    Returns ``(pi_species, applies)``. ``applies=False`` means "cannot assess the
    separation" — the caller should leave the score unchanged (conservative no-op),
    matching the historical NaN-valley_depth skip.

    π_species keys on the OPTIMISTIC bootstrap upper-confidence-limit of gap_fraction
    (see compute_valley_depth), so the low-N over-kill guard is already baked into the
    input:
      * UCL None/NaN → (pi_train, False): too few confident U12 calls to assess.
      * UCL == -inf  → (prior_floor, True): even the optimistic estimate shows no
        separation (most resamples inverted/coincident) ⇒ confident suppression.
      * UCL finite   → graded:
            w = sigmoid((4.394/width)·(UCL − midpoint))   # optimistic separation → trust
            π_species = prior_floor + w·(π_train − prior_floor)
        A wide bootstrap lifts the UCL above the point estimate ⇒ w→1 ⇒ no suppression,
        so we suppress only when CONFIDENT (even optimistically) the separation is poor.
        The UCL is percentile-based, hence well-behaved where gap_fraction crosses zero
        / Δmean → 0 (unlike σ/μ).
    """
    if gap_fraction_ucl is None:
        return pi_train, False
    ucl = float(gap_fraction_ucl)
    if math.isnan(ucl):
        return pi_train, False
    if math.isinf(ucl):
        # UCL = -inf: even the optimistic estimate shows no separation → suppress hard.
        return prior_floor, True

    # Overflow-safe logistic: UCL can be large-negative (deep overlap), which would
    # overflow math.exp; clamp to the saturated tails.
    z = (4.394 / gap_fraction_width) * (ucl - gap_fraction_midpoint)
    if z <= -700.0:
        w = 0.0
    elif z >= 700.0:
        w = 1.0
    else:
        w = 1.0 / (1.0 + math.exp(-z))

    pi_species = prior_floor + w * (pi_train - prior_floor)
    return pi_species, True


def compute_adjusted_score(
    svm_score: float,
    gap_fraction_ucl: Optional[float],
    ensemble_sigma: float,
    gap_fraction_midpoint: float = DEFAULT_GAP_FRACTION_MIDPOINT,
    gap_fraction_width: float = DEFAULT_GAP_FRACTION_WIDTH,
    prior_floor: float = 0.001,
    k_sigma: float = 3.0,
    pi_train: float = 0.5,
) -> float:
    """Compute the gap_fraction-prior + σ adjusted confidence score for one intron.

    Combines three signals in log-odds space:
        logit(p_adj) = logit(p_svm) + log(π_species / π_train) − k_σ · σ

    where π_species is the species prior derived from the confidence-shrunk
    `gap_fraction` upper-confidence-limit (see `_confidence_shrunk_pi_species`).
    When the separation cannot be assessed (UCL None/NaN) the prior correction is
    0 (the SVM score passes through; the caller normally skips these rows entirely).

    Args:
        svm_score: Raw SVM ensemble mean probability (0-100 scale)
        gap_fraction_ucl: Bootstrap upper-confidence-limit of the species separation
            statistic (Fisher-axis gap/Δmean); the optimistic, low-N-shrunk estimate
        ensemble_sigma: Std of per-model probabilities (0-100 scale)
        gap_fraction_midpoint: Separation at which suppression is half-on
        gap_fraction_width: Width of the gap_fraction→trust transition
        prior_floor: Minimum species-level prior
        k_sigma: Ensemble disagreement penalty coefficient (0 disables)
        pi_train: Training prior (0.5 for balanced)

    Returns:
        Adjusted probability (0-100 scale)
    """
    eps = 1e-9
    p = max(eps, min(1 - eps, svm_score / 100.0))

    # Term 1: instance-level evidence
    logit_svm = math.log(p / (1 - p))

    # Term 2: confidence-shrunk species-level prior correction
    pi_species, applies = _confidence_shrunk_pi_species(
        gap_fraction_ucl,
        gap_fraction_midpoint=gap_fraction_midpoint,
        gap_fraction_width=gap_fraction_width,
        prior_floor=prior_floor, pi_train=pi_train,
    )
    prior_correction = math.log(pi_species / pi_train) if applies else 0.0

    # Term 3: epistemic uncertainty penalty (σ on 0-1 scale)
    sigma_penalty = k_sigma * (ensemble_sigma / 100.0)

    logit_adj = logit_svm + prior_correction - sigma_penalty
    return 100.0 / (1.0 + math.exp(-logit_adj))


def compute_adjusted_scores_batch(
    svm_scores: np.ndarray,
    gap_fraction_ucl: Optional[float],
    ensemble_sigmas: np.ndarray,
    gap_fraction_midpoint: float = DEFAULT_GAP_FRACTION_MIDPOINT,
    gap_fraction_width: float = DEFAULT_GAP_FRACTION_WIDTH,
    prior_floor: float = 0.001,
    k_sigma: float = 3.0,
    pi_train: float = 0.5,
) -> np.ndarray:
    """Vectorized score adjustment for a batch of introns.

    Same formula as compute_adjusted_score but operates on arrays. The species
    prior correction is a single scalar (gap_fraction_ucl is per-species), so it
    is computed once and broadcast.

    Args:
        svm_scores: Raw SVM probabilities (0-100 scale), shape (n,)
        gap_fraction_ucl: Bootstrap UCL of the species separation statistic (scalar)
        ensemble_sigmas: Per-intron ensemble σ (0-100 scale), shape (n,)
        Other args: see compute_adjusted_score

    Returns:
        Adjusted probabilities (0-100 scale), shape (n,)
    """
    eps = 1e-9
    p = np.clip(svm_scores / 100.0, eps, 1 - eps)

    logit_svm = np.log(p / (1 - p))

    pi_species, applies = _confidence_shrunk_pi_species(
        gap_fraction_ucl,
        gap_fraction_midpoint=gap_fraction_midpoint,
        gap_fraction_width=gap_fraction_width,
        prior_floor=prior_floor, pi_train=pi_train,
    )
    prior_correction = math.log(pi_species / pi_train) if applies else 0.0

    sigma_penalty = k_sigma * (ensemble_sigmas / 100.0)

    logit_adj = logit_svm + prior_correction - sigma_penalty
    return 100.0 / (1.0 + np.exp(-logit_adj))


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
    three_z_scores: Optional[np.ndarray] = None,
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
        three_z_scores: Array of 3' z-scores for all scored introns. When
            provided, valley detection runs in 3D (5'z, BPz, 3'z) using
            Fisher's discriminant. When omitted, falls back to 2D (5'z, BPz)
            for callers that don't have 3'z available.

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
            # NaN gap_fraction ⇒ π_species cannot assess ⇒ no-op (conservative).
            'gap_fraction': float('nan'),
            'gap_fraction_ucl': float('nan'),
        }

    # Stack feature columns. 3D (5'z, BPz, 3'z) is the production default;
    # 2D fallback for callers that lack 3'z (kept for backwards compat).
    if three_z_scores is not None:
        u12_points = np.column_stack([
            five_z_scores[u12_mask], bp_z_scores[u12_mask],
            three_z_scores[u12_mask],
        ])
        u2_points = np.column_stack([
            five_z_scores[u2_mask], bp_z_scores[u2_mask],
            three_z_scores[u2_mask],
        ])
    else:
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
        'gap_fraction': valley_result.get('gap_fraction'),
        'gap_fraction_ucl': valley_result.get('gap_fraction_ucl'),
    }


def apply_valley_prior_adjustment(
    introns: List["Intron"],
    validation_result: dict,
    threshold: float = 90.0,
) -> List["Intron"]:
    """
    Apply valley-based prior adjustment to classified introns.

    If the valley detection found no distinct U12-type cluster (no_valley or
    insufficient regime), adjusts SVM probabilities downward using the
    computed adjusted prior. Introns whose adjusted probability drops below
    50% are reclassified as U2-type.

    If a valley was detected (valley regime), returns introns unchanged.

    Args:
        introns: List of classified Intron objects
        validation_result: Output from validate_u12_cluster()
        threshold: Original classification threshold (default: 90.0)

    Returns:
        List of Intron objects with potentially adjusted scores and type_ids
    """
    from intronIC.scoring.prior_adjustment import adjust_probabilities_for_prior
    from intronIC.core.intron import IntronScores, IntronMetadata

    adjusted_prior = validation_result.get('adjusted_prior', 0.5)

    # No adjustment needed if valley detected (prior = 0.5)
    if adjusted_prior >= 0.5:
        return introns

    training_prior = 0.5  # SVM's implicit balanced prior

    result = []
    for intron in introns:
        # Only adjust scored introns with U12-type calls
        if (intron.scores is None or
                intron.scores.svm_score is None or
                intron.metadata is None or
                intron.metadata.type_id != 'u12'):
            result.append(intron)
            continue

        # Apply Bayesian prior adjustment
        original_prob = intron.scores.svm_score / 100.0
        adjusted_prob = adjust_probabilities_for_prior(
            original_prob, training_prior, adjusted_prior
        )
        adjusted_svm = adjusted_prob * 100.0
        adjusted_rel = adjusted_svm - threshold

        # Reclassify if adjusted probability < 50%
        if adjusted_prob < 0.5:
            new_type_id = 'u2'
        else:
            new_type_id = 'u12'

        # Compute log-odds for adjusted probability
        if 0 < adjusted_prob < 1:
            adjusted_dd = math.log(adjusted_prob / (1 - adjusted_prob))
        elif adjusted_prob >= 1:
            adjusted_dd = 23.03  # cap
        else:
            adjusted_dd = -23.03

        new_scores = replace(
            intron.scores,
            svm_score=adjusted_svm,
            relative_score=adjusted_rel,
            decision_distance=adjusted_dd,
        )
        new_metadata = replace(
            intron.metadata,
            type_id=new_type_id,
        )
        result.append(replace(intron, scores=new_scores, metadata=new_metadata))

    return result


def rewrite_outputs_with_prior_adjustment(
    meta_path: Path,
    score_info_path: Path,
    validation_result: dict,
    threshold: float = 90.0,
) -> int:
    """
    Rewrite meta.iic and score_info.iic files with valley-adjusted scores.

    For the sequential streaming path where introns are written per-contig
    and not buffered. Reads the output files, applies prior adjustment to
    U12-type calls, and rewrites with adjusted scores and type assignments.

    Only modifies files if the valley detection found no distinct cluster
    (adjusted_prior < 0.5). Returns the number of introns reclassified.

    Args:
        meta_path: Path to the .meta.iic output file
        score_info_path: Path to the .score_info.iic output file
        validation_result: Output from validate_u12_cluster()
        threshold: Classification threshold (default: 90.0)

    Returns:
        Number of introns reclassified from U12-type to U2-type
    """
    from intronIC.scoring.prior_adjustment import adjust_probabilities_for_prior

    adjusted_prior = validation_result.get('adjusted_prior', 0.5)
    if adjusted_prior >= 0.5:
        return 0  # No adjustment needed

    training_prior = 0.5
    n_reclassified = 0

    # Rewrite meta.iic: adjust type_id and rel_score columns
    if meta_path.exists():
        tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.iic',
                                           dir=meta_path.parent, delete=False)
        try:
            with open(meta_path) as f_in:
                header_line = None
                type_col = None
                rel_col = None

                for line in f_in:
                    if line.startswith('#'):
                        tmp.write(line)
                        continue

                    parts = line.rstrip('\n').split('\t')

                    # Parse header
                    if header_line is None:
                        header_line = parts
                        try:
                            type_col = parts.index('type_id')
                            rel_col = parts.index('rel_score')
                        except ValueError:
                            pass
                        tmp.write(line)
                        continue

                    # Adjust U12 calls
                    if type_col is not None and parts[type_col] == 'u12':
                        try:
                            rel_score = float(parts[rel_col])
                            svm_score = rel_score + threshold
                            original_prob = svm_score / 100.0
                            adjusted_prob = adjust_probabilities_for_prior(
                                original_prob, training_prior, adjusted_prior
                            )
                            if adjusted_prob < 0.5:
                                parts[type_col] = 'u2'
                                adjusted_svm = adjusted_prob * 100.0
                                parts[rel_col] = f"{adjusted_svm - threshold:.4f}"
                                n_reclassified += 1
                        except (ValueError, IndexError):
                            pass

                    tmp.write('\t'.join(parts) + '\n')

            tmp.close()
            shutil.move(tmp.name, meta_path)
        except Exception:
            tmp.close()
            Path(tmp.name).unlink(missing_ok=True)
            raise

    # Rewrite score_info.iic: adjust svm_score, rel_score, decision_dist
    if score_info_path.exists():
        tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.iic',
                                           dir=score_info_path.parent, delete=False)
        try:
            with open(score_info_path) as f_in:
                header_line = None
                svm_col = rel_col = dd_col = None

                for line in f_in:
                    if line.startswith('#'):
                        tmp.write(line)
                        continue

                    parts = line.rstrip('\n').split('\t')

                    if header_line is None:
                        header_line = parts
                        try:
                            svm_col = parts.index('svm_score')
                            rel_col = parts.index('rel_score')
                            dd_col = parts.index('decision_dist')
                        except ValueError:
                            pass
                        tmp.write(line)
                        continue

                    # Adjust scores for introns that were reclassified
                    if svm_col is not None:
                        try:
                            svm_score = float(parts[svm_col])
                            original_prob = svm_score / 100.0
                            if original_prob >= threshold / 100.0:  # was U12
                                adjusted_prob = adjust_probabilities_for_prior(
                                    original_prob, training_prior, adjusted_prior
                                )
                                adjusted_svm = adjusted_prob * 100.0
                                parts[svm_col] = f"{adjusted_svm:.4f}"
                                if rel_col is not None:
                                    parts[rel_col] = f"{adjusted_svm - threshold:.4f}"
                                if dd_col is not None:
                                    if 0 < adjusted_prob < 1:
                                        dd = math.log(adjusted_prob / (1 - adjusted_prob))
                                    else:
                                        dd = 23.03 if adjusted_prob >= 1 else -23.03
                                    parts[dd_col] = f"{dd:.4f}"
                        except (ValueError, IndexError):
                            pass

                    tmp.write('\t'.join(parts) + '\n')

            tmp.close()
            shutil.move(tmp.name, score_info_path)
        except Exception:
            tmp.close()
            Path(tmp.name).unlink(missing_ok=True)
            raise

    return n_reclassified
