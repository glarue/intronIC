"""
Prior-aware probability adjustment for cross-species U12 classification.

This module implements Bayesian prior adjustment to account for species-specific
U12 base rates. This is particularly important for U12-absent species where the
training prior (human ≈ 0.5%) significantly overestimates the true prior.

Theory:
    If a classifier is trained on data with prior π_train and applied to data with
    prior π_target, the posterior probabilities need to be adjusted via Bayes' rule.

    Given:
    - p_train: Probability from model trained on human data (π_train ≈ 0.005)
    - π_target: True prior for target species

    Compute:
    - r = odds_target / odds_train = [π_target/(1-π_target)] / [π_train/(1-π_train)]
    - p_adjusted = (r * p_train) / (r * p_train + (1 - p_train))

Use cases:
    1. U12-absent species (C. elegans, many fungi): π_target ≈ 1e-6
       → Drastically reduces false positive rate
    2. U12-rich species (some vertebrates): π_target ≈ 0.01
       → Slightly increases sensitivity
    3. Uncertain species: Use training prior (π_target = π_train)
       → No adjustment

Design principles:
    1. Modular: Separate from model and normalization
    2. Conservative: Only adjust if user explicitly requests
    3. Transparent: Log effect of adjustment
    4. Numerically stable: Prevent division by zero

Reference:
    Expert guidance on cross-species adaptation (2025)
    "Add guardrails so adaptation cannot quietly destroy false-positive control"
"""

import numpy as np
from typing import Union


def adjust_probabilities_for_prior(
    probabilities: Union[np.ndarray, float],
    training_prior: float,
    target_prior: float,
    epsilon: float = 1e-9
) -> Union[np.ndarray, float]:
    """
    Adjust classification probabilities for different class priors.

    Uses Bayes' rule to shift posterior probabilities when the target
    species has a different U12 base rate than the training species.

    Args:
        probabilities: Probability or probabilities from model [0, 1]
                      Can be scalar or array
        training_prior: U12 prior in training data (e.g., 0.005 for human)
        target_prior: Expected U12 prior in target species
                     Examples:
                     - 1e-6: U12-absent (C. elegans, many fungi)
                     - 0.005: Human-like
                     - 0.01: U12-rich vertebrates
        epsilon: Small value to prevent division by zero (default: 1e-9)

    Returns:
        Adjusted probabilities (same shape as input)

    Algorithm:
        1. Compute odds ratio: r = [π_t/(1-π_t)] / [π_h/(1-π_h)]
        2. Apply Bayes' rule: p' = (r*p) / (r*p + (1-p))
        3. Clamp to [0, 1] for numerical stability

    Example:
        >>> # Human model gives p=0.90 (90% confident U12)
        >>> # But target is C. elegans (U12-absent)
        >>> p_adj = adjust_probabilities_for_prior(
        ...     probabilities=0.90,
        ...     training_prior=0.005,
        ...     target_prior=1e-6
        ... )
        >>> # Result: p_adj ≈ 0.002 (dramatically reduced)
        >>> # This prevents false positives in U12-free species

    Example:
        >>> # Batch adjustment
        >>> probs = np.array([0.95, 0.80, 0.50, 0.20])
        >>> probs_adj = adjust_probabilities_for_prior(
        ...     probs, training_prior=0.005, target_prior=1e-6
        ... )
        >>> # All probabilities reduced proportionally

    Mathematical Derivation:
        Let:
        - H: hypothesis "intron is U12"
        - E: evidence (features)
        - π_h: P(H) in training data
        - π_t: P(H) in target data

        Training classifier gives:
            P(H|E, π_h) = p_train

        Want:
            P(H|E, π_t) = ?

        By Bayes' rule:
            P(H|E) ∝ P(E|H) * P(H)

        So:
            P(E|H) = P(H|E, π_h) * P(E) / π_h

        Therefore:
            P(H|E, π_t) = P(E|H) * π_t / P(E)
                        = [P(H|E, π_h) / π_h] * π_t * [P(E)/P(E)]

        After normalization:
            P(H|E, π_t) = (r * p_train) / (r * p_train + (1 - p_train))

        where r = π_t * (1-π_h) / [π_h * (1-π_t)]
    """
    # Input validation
    if not (0.0 < training_prior < 1.0):
        raise ValueError(
            f"training_prior must be in (0, 1), got {training_prior}"
        )
    if not (0.0 < target_prior < 1.0):
        raise ValueError(
            f"target_prior must be in (0, 1), got {target_prior}"
        )

    # Handle scalar input
    is_scalar = np.isscalar(probabilities)
    probs_array = np.atleast_1d(probabilities)

    # Validate probability range
    if np.any((probs_array < 0) | (probs_array > 1)):
        raise ValueError(
            "probabilities must be in [0, 1], got values outside this range"
        )

    # Compute odds ratio
    # r = [π_t/(1-π_t)] / [π_h/(1-π_h)]
    odds_train = training_prior / max(1.0 - training_prior, epsilon)
    odds_target = target_prior / max(1.0 - target_prior, epsilon)
    r = odds_target / max(odds_train, epsilon)

    # Apply Bayes' rule
    # p_adj = (r * p) / (r * p + (1 - p))
    numerator = r * probs_array
    denominator = numerator + (1.0 - probs_array)

    # Prevent division by zero
    # This can happen if p=1.0 (certain U12) and r is very small
    p_adjusted = numerator / np.maximum(denominator, epsilon)

    # Clamp to [0, 1] for numerical stability
    p_adjusted = np.clip(p_adjusted, 0.0, 1.0)

    # Return scalar if input was scalar
    if is_scalar:
        return float(p_adjusted[0])
    return p_adjusted


def estimate_prior_from_counts(
    n_u12: int,
    n_u2: int,
    pseudocount: float = 0.5
) -> float:
    """
    Estimate U12 prior from reference counts with Laplace smoothing.

    Args:
        n_u12: Number of U12 introns
        n_u2: Number of U2 introns
        pseudocount: Pseudocount for Laplace smoothing (default: 0.5)

    Returns:
        Estimated prior: (n_u12 + α) / (n_u12 + n_u2 + 2α)

    Example:
        >>> # Human: ~500 U12, ~100,000 U2
        >>> prior = estimate_prior_from_counts(500, 100000)
        >>> # Result: ≈ 0.005 (0.5%)
    """
    if n_u12 < 0 or n_u2 < 0:
        raise ValueError(
            f"Counts must be non-negative, got n_u12={n_u12}, n_u2={n_u2}"
        )

    if n_u12 + n_u2 == 0:
        raise ValueError("Cannot estimate prior from zero counts")

    # Laplace smoothing to avoid zero probabilities
    return (n_u12 + pseudocount) / (n_u12 + n_u2 + 2 * pseudocount)


def compute_prior_adjustment_diagnostics(
    probabilities: np.ndarray,
    adjusted_probabilities: np.ndarray,
    training_prior: float,
    target_prior: float,
    threshold: float = 0.9
) -> dict:
    """
    Compute diagnostics for prior adjustment.

    Useful for logging and understanding the effect of adjustment.

    Args:
        probabilities: Original probabilities
        adjusted_probabilities: After prior adjustment
        training_prior: Training U12 prior
        target_prior: Target U12 prior
        threshold: Classification threshold (default: 0.9)

    Returns:
        Dictionary with diagnostics:
        - odds_ratio: r = odds_target / odds_train
        - n_u12_before: Number predicted U12 before adjustment
        - n_u12_after: Number predicted U12 after adjustment
        - mean_prob_before: Mean probability before
        - mean_prob_after: Mean probability after
        - median_prob_before: Median probability before
        - median_prob_after: Median probability after

    Example:
        >>> probs = np.random.beta(2, 8, size=10000)  # Skewed toward 0
        >>> probs_adj = adjust_probabilities_for_prior(probs, 0.005, 1e-6)
        >>> diag = compute_prior_adjustment_diagnostics(
        ...     probs, probs_adj, 0.005, 1e-6, threshold=0.9
        ... )
        >>> print(f"U12 predictions: {diag['n_u12_before']} → {diag['n_u12_after']}")
    """
    # Compute odds ratio
    epsilon = 1e-9
    odds_train = training_prior / max(1.0 - training_prior, epsilon)
    odds_target = target_prior / max(1.0 - target_prior, epsilon)
    odds_ratio = odds_target / max(odds_train, epsilon)

    # Count classifications
    n_u12_before = int(np.sum(probabilities >= threshold))
    n_u12_after = int(np.sum(adjusted_probabilities >= threshold))

    # Summary statistics
    mean_before = float(np.mean(probabilities))
    mean_after = float(np.mean(adjusted_probabilities))
    median_before = float(np.median(probabilities))
    median_after = float(np.median(adjusted_probabilities))

    return {
        'odds_ratio': odds_ratio,
        'n_u12_before': n_u12_before,
        'n_u12_after': n_u12_after,
        'frac_u12_before': n_u12_before / len(probabilities),
        'frac_u12_after': n_u12_after / len(probabilities),
        'mean_prob_before': mean_before,
        'mean_prob_after': mean_after,
        'median_prob_before': median_before,
        'median_prob_after': median_after
    }
