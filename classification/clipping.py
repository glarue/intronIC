"""
Symmetric clipping transformer for z-space features.

Clips extreme values at specified quantiles to prevent outliers from
dominating the linear SVM decision boundary.

Architectural context:
In the redesigned pipeline, this operates AFTER zero-anchored scaling,
clipping z-scores to prevent cross-species composition artifacts from
overwhelming balance features.

Expert guidance: "Rare huge LLRs (e.g., 5' outliers) can dominate a linear
margin and calibration. Aggressive clipping (and optionally a saturating
transform) so a single insane PWM score can't swamp balance features."

Redesign: SCALER_ARCHITECTURE_REVIEW.md (Expert feedback 2025)
Original: FP_REDUCTION_IMPLEMENTATION.md (Tier 2, Step 5)
"""

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin


class SymmetricClipper(BaseEstimator, TransformerMixin):
    """
    Symmetric clipping of z-scores to prevent extreme outliers.

    In the redesigned pipeline, this operates in z-space AFTER zero-anchored
    scaling. Clips extreme z-scores to prevent cross-species composition
    artifacts from overwhelming balance features.

    Pipeline position:
        ZeroAnchoredRobustScaler → SymmetricClipper → SaturatingTransform → ...

    Expert guidance: "Aggressive clipping (and optionally a saturating
    transform) so a single insane PWM score can't swamp balance features."

    **DOMAIN-ADAPTIVE MODE (2025 Cross-Species Fix):**

    When `domain_adaptive=True` (RECOMMENDED for cross-species deployment):
    - During training: Learns quantile q (dimensionless hyperparameter)
    - At inference: Recomputes clip bounds from **target species** z-scores
    - Fixes catastrophic C. elegans failure where human bounds ±1.48σ
      were applied to C. elegans z∈[-63, +16], destroying all signal

    Key insight: "q transfers, not numeric thresholds"
    - Human: q=0.95 → caps from human z-score distribution
    - C. elegans: q=0.95 → caps from C. elegans z-score distribution

    When `domain_adaptive=False` (legacy behavior):
    - Learns numeric clip bounds from training data
    - Applies those fixed bounds at inference (fails cross-species)

    The clipping is:
    - Symmetric: Uses quantile of |z| to get cap, clips to [-cap, +cap]
    - Per-feature: Each z-score (5'SS, BPS, 3'SS) gets independent threshold
    - Aggressive: Default Q_0.975 clips ~2.5% of extreme values (vs old 0.1%)
    - Zero-preserving: Does not affect the semantic zero point (z=0 → z=0)
    - Hyperparameter: quantile tuned via grid search for optimal FP reduction

    Attributes:
        quantile: Quantile for symmetric clipping (default: 0.975)
                 Recommended range: [0.95, 0.975, 0.99] (grid search)
                 Lower = more aggressive clipping
        domain_adaptive: Whether to recompute caps per-species (default: True)
        caps_: Learned clipping thresholds per feature (set during fit)
               Only used when domain_adaptive=False

    Example (Domain-Adaptive Mode - RECOMMENDED):
        >>> clipper = SymmetricClipper(quantile=0.95, domain_adaptive=True)
        >>> clipper.fit(z_human)  # Store quantile=0.95
        >>> # C. elegans with extreme composition bias
        >>> z_celegans = np.array([[-0.15, 0.30, 11.45]])
        >>> z_clipped = clipper.transform(z_celegans)
        >>> # Caps recomputed from C. elegans data using q=0.95
        >>> # Preserves signal instead of destroying it

    Example (Legacy Mode - Cross-Species Fails):
        >>> clipper = SymmetricClipper(quantile=0.95, domain_adaptive=False)
        >>> clipper.fit(z_human)  # Learn numeric caps [2.0, 1.9, 1.5]
        >>> z_celegans = np.array([[-0.15, 0.30, 11.45]])
        >>> z_clipped = clipper.transform(z_celegans)
        >>> # Applies human caps to C. elegans → destroys signal

    Redesign: SCALER_CENTERING_FIX_SUMMARY.md (Cross-species fix, 2025)
    """

    def __init__(self, quantile: float = 0.975, domain_adaptive: bool = True):
        """
        Initialize SymmetricClipper.

        Args:
            quantile: Quantile for clipping (0.975 = clip at 97.5th percentile)
                     Lower values = more aggressive clipping
                     Grid search recommendations: [0.95, 0.975, 0.99]
                     - 0.95: Very aggressive (clips ~5% of extreme values)
                     - 0.975: Moderate (clips ~2.5%, default)
                     - 0.99: Conservative (clips ~1%)
            domain_adaptive: Whether to recompute caps from input data (default: True)
                     True: Recompute caps at transform time (cross-species safe)
                     False: Use fixed caps from fit (legacy, cross-species fails)

        Raises:
            ValueError: If quantile not in (0, 1)
        """
        if not 0 < quantile < 1:
            raise ValueError(f"quantile must be in (0, 1), got {quantile}")
        self.quantile = quantile
        self.domain_adaptive = domain_adaptive
        self.caps_ = None

    def fit(self, X, y=None):
        """
        Learn clipping thresholds from z-score training data.

        Computes symmetric caps at the specified quantile of absolute values
        for each feature independently. This preserves the zero-anchoring of
        z-scores (z=0 means "U12≈U2").

        Args:
            X: Training z-scores of shape (n_samples, n_features)
               For intronIC: (n_introns, 3) for [z5, zBP, z3]
               Typically from ZeroAnchoredRobustScaler output
            y: Target values (ignored, for sklearn compatibility)

        Returns:
            self (fitted clipper with caps_ attribute set)

        Example:
            >>> # z-scores from scaler
            >>> z_train = np.array([[0.5, 1.2], [1.5, 1.8], [0.3, 5.5]])
            >>> clipper = SymmetricClipper(quantile=0.95)
            >>> clipper.fit(z_train)
            >>> clipper.caps_  # [1.5, 5.5] (95th percentile of |z|)
        """
        X = np.asarray(X)

        # Compute symmetric caps: quantile of |z| per feature
        # Shape: (n_features,)
        # Uses absolute values to ensure symmetric clipping
        self.caps_ = np.quantile(np.abs(X), self.quantile, axis=0)

        return self

    def transform(self, X):
        """
        Apply clipping thresholds to z-scores.

        **Domain-Adaptive Mode (domain_adaptive=True):**
        Recomputes caps from input data X using the learned quantile.
        This allows the same quantile hyperparameter to adapt to different
        species' z-score distributions.

        **Legacy Mode (domain_adaptive=False):**
        Uses fixed caps learned during fit() on training data.

        Args:
            X: Z-scores to clip of shape (n_samples, n_features)
               Same format as fit() input

        Returns:
            X_clipped: Clipped z-scores of same shape
                      Values outside [-cap_i, +cap_i] are clamped to bounds

        Raises:
            ValueError: If transform called before fit (legacy mode only)

        Example (Domain-Adaptive):
            >>> clipper = SymmetricClipper(quantile=0.95, domain_adaptive=True)
            >>> clipper.fit(z_human)  # Stores quantile=0.95
            >>> # C. elegans with different z-score distribution
            >>> z_celegans = np.array([[0.2, 11.45]])
            >>> z_clipped = clipper.transform(z_celegans)
            >>> # Caps computed from C. elegans using q=0.95

        Example (Legacy):
            >>> clipper = SymmetricClipper(quantile=0.95, domain_adaptive=False)
            >>> clipper.fit(z_human)  # Learn numeric caps
            >>> z_celegans = np.array([[0.2, 11.45]])
            >>> z_clipped = clipper.transform(z_celegans)
            >>> # Uses human caps (fails cross-species)
        """
        X = np.asarray(X)

        if self.domain_adaptive:
            # Domain-adaptive mode: Recompute caps from input data
            # This allows quantile to transfer across species
            caps = np.quantile(np.abs(X), self.quantile, axis=0)
        else:
            # Legacy mode: Use fixed caps from fit()
            if self.caps_ is None:
                raise ValueError("SymmetricClipper must be fit before transform")
            caps = self.caps_

        # Clip each feature independently
        X_clipped = np.copy(X)
        for i, cap in enumerate(caps):
            X_clipped[:, i] = np.clip(X_clipped[:, i], -cap, cap)

        return X_clipped

    def get_feature_names_out(self, input_features=None):
        """
        Get output feature names (same as input).

        Clipping is a transform that doesn't change feature names,
        only their value ranges.

        Args:
            input_features: Input feature names (optional)
                          If None, generates generic names

        Returns:
            Array of output feature names (same as input)

        Example:
            >>> clipper = SymmetricClipper()
            >>> clipper.caps_ = np.array([3.5, 3.2, 3.8])
            >>> clipper.get_feature_names_out(['z5', 'zBP', 'z3'])
            array(['z5', 'zBP', 'z3'], dtype='<U3')
        """
        if input_features is None:
            # If not provided, use generic names
            n_features = len(self.caps_) if self.caps_ is not None else 0
            return np.array([f'x{i}' for i in range(n_features)])
        return np.asarray(input_features)


# Backward compatibility alias (for old code that imports OutlierClipper)
OutlierClipper = SymmetricClipper
