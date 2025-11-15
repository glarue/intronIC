"""
Custom sklearn transformers for feature augmentation.

This module provides transformers for augmenting base features with
composite features that help the linear SVM better separate U12 from U2.
"""

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin


class BothEndsStrongTransformer(BaseEstimator, TransformerMixin):
    """
    Augment 3D features (5'SS, BPS, 3'SS z-scores) with composite features
    to capture "both ends strong" patterns and penalize imbalance.

    This transformer builds composite features using min operations and negated
    absolute differences to help the linear SVM identify and reject false positives
    that show "one-end-strong" patterns.

    Port from: Expert recommendations for reducing false positives
               See: FP_REDUCTION_IMPLEMENTATION.md

    Background:
        U12-type introns have strong correlation between 5'SS and BPS motifs.
        False positives often show "one-end-strong" patterns: one very high
        score compensating for weak/negative scores elsewhere.

        By explicitly computing min features and imbalance penalties, we make
        it easy for the linear model to:
        1. Require all signals to be strong (via min_all)
        2. Penalize imbalance between signals (via neg_absdiff_*)
        3. Optionally include "at least one strong" signals (via max_*)

    Mathematical insight:
        min(a, b) = 0.5 * ((a + b) - |a - b|)
        max(a, b) = 0.5 * ((a + b) + |a - b|)
        neg_absdiff(a, b) = -|a - b|  ← Penalty for imbalance

        Expert guidance: "model will mostly weight min" and "even with min_*,
        the model may not penalize one-end-strong enough" → add neg_absdiff.

    Attributes:
        include_max: Whether to include max features (default: False)
                    If True, adds max_5_bp and max_5_3 features
                    Expert: "you can drop max_* entirely"

    Input features (3D):
        - s5:  5' splice site z-score (LLR, zero-anchored)
        - sBP: Branch point z-score (LLR, zero-anchored)
        - s3:  3' splice site z-score (LLR, zero-anchored)

    Output features (9D default, 11D with include_max=True):
        1. s5 (original, passed through)
        2. sBP (original, passed through)
        3. s3 (original, passed through)
        4. min_5_bp = min(s5, sBP)              ← Both 5' and BP must be strong
        5. min_5_3 = min(s5, s3)                ← Both 5' and 3' must be strong
        6. min_all = min(s5, sBP, s3)           ← ALL THREE must be strong (3-way AND)
        7. neg_absdiff_5_bp = -|s5 - sBP|       ← Penalty for 5'/BP imbalance
        8. neg_absdiff_5_3 = -|s5 - s3|         ← Penalty for 5'/3' imbalance
        9. neg_absdiff_bp_3 = -|sBP - s3|       ← Penalty for BP/3' imbalance
        [10. max_5_bp = max(s5, sBP)]           ← At least one of 5'/BP strong (optional)
        [11. max_5_3 = max(s5, s3)]             ← At least one of 5'/3' strong (optional)

    Example (balanced U12 intron):
        >>> transformer = BothEndsStrongTransformer(include_max=False)
        >>> X = np.array([[2.0, 1.8, 1.5]])  # All strong and balanced
        >>> X_aug = transformer.transform(X)
        >>> # min_all = 1.5 (all three strong)
        >>> # neg_absdiff_5_bp = -0.2 (small penalty, signals consistent)

    Example (one-end-strong FP):
        >>> X_bad = np.array([[5.0, -1.0, -0.5]])  # Only 5'SS strong
        >>> X_bad_aug = transformer.transform(X_bad)
        >>> # min_all = -1.0 (low, not all strong)
        >>> # neg_absdiff_5_bp = -6.0 (large penalty, huge imbalance)
        >>> # These features help SVM reject this as FP
    """

    def __init__(
        self,
        include_max: bool = False,
        gamma_5_bp=None,  # Deprecated - kept for backward compatibility
        gamma_5_3=None    # Deprecated - kept for backward compatibility
    ):
        """
        Initialize BothEndsStrongTransformer.

        Args:
            include_max: Whether to include max features (default: False)
                        Expert guidance: "model will mostly weight min"
            gamma_5_bp: DEPRECATED - Ignored, kept for backward compatibility with old models
            gamma_5_3: DEPRECATED - Ignored, kept for backward compatibility with old models
        """
        self.include_max = include_max

        # Backward compatibility: Store gamma parameters as attributes even though they're unused
        # This allows old pickled models to load without errors
        # Old models used gamma-weighted sum+absdiff, new models use min/max
        self.gamma_5_bp = gamma_5_bp
        self.gamma_5_3 = gamma_5_3

        if gamma_5_bp is not None or gamma_5_3 is not None:
            import warnings
            warnings.warn(
                "gamma_5_bp and gamma_5_3 are deprecated. "
                "The transformer now uses min/max features instead of gamma-weighted features. "
                "Old models will work but should be retrained with the new approach.",
                DeprecationWarning,
                stacklevel=2
            )

    def fit(self, X, y=None):
        """
        Fit transformer (no-op, this transformer is stateless).

        Args:
            X: Training data (ignored)
            y: Target values (ignored)

        Returns:
            self
        """
        return self

    def transform(self, X):
        """
        Augment 3D features to 9D (or 11D with include_max=True).

        Args:
            X: Array of shape (n_samples, 3) with [s5, sBP, s3]

        Returns:
            Array of shape (n_samples, 9) or (n_samples, 11) with augmented features
            Features: [s5, sBP, s3, min_5_bp, min_5_3, min_all,
                      neg_absdiff_5_bp, neg_absdiff_5_3, neg_absdiff_bp_3,
                      [max_5_bp, max_5_3]]

        Raises:
            ValueError: If input doesn't have exactly 3 features
        """
        # Ensure input is numpy array
        X = np.asarray(X)

        # Validate input shape
        if X.shape[1] != 3:
            raise ValueError(
                f"BothEndsStrongTransformer expects 3 input features "
                f"[s5, sBP, s3], got {X.shape[1]}"
            )

        # Extract base features
        s5 = X[:, 0]
        sBP = X[:, 1]
        s3 = X[:, 2]

        # Build pairwise min/max features
        # For 5'SS-BPS (strongest correlation)
        # min(a, b) = 0.5 * ((a + b) - |a - b|)
        # max(a, b) = 0.5 * ((a + b) + |a - b|)
        sum_5_bp = s5 + sBP
        absdiff_5_bp = np.abs(s5 - sBP)
        min_5_bp = 0.5 * (sum_5_bp - absdiff_5_bp)

        # For 5'SS-3'SS (secondary correlation)
        sum_5_3 = s5 + s3
        absdiff_5_3 = np.abs(s5 - s3)
        min_5_3 = 0.5 * (sum_5_3 - absdiff_5_3)

        # 3-way AND: ALL THREE must be strong (NEW)
        # min_all = min(s5, sBP, s3)
        min_all = np.minimum(np.minimum(s5, sBP), s3)

        # Imbalance penalties (NEW)
        # Negative absdiff = penalty for inconsistency
        # Model puts positive weight on these → rewards balanced signals
        neg_absdiff_5_bp = -absdiff_5_bp
        neg_absdiff_5_3 = -absdiff_5_3

        # For BPS-3'SS imbalance
        absdiff_bp_3 = np.abs(sBP - s3)
        neg_absdiff_bp_3 = -absdiff_bp_3

        # Stack base + min + penalty features
        features = [
            s5[:, np.newaxis],
            sBP[:, np.newaxis],
            s3[:, np.newaxis],
            min_5_bp[:, np.newaxis],
            min_5_3[:, np.newaxis],
            min_all[:, np.newaxis],           # NEW
            neg_absdiff_5_bp[:, np.newaxis],  # NEW
            neg_absdiff_5_3[:, np.newaxis],   # NEW
            neg_absdiff_bp_3[:, np.newaxis]   # NEW
        ]

        # Optional: max features (expert: "you can drop max_* entirely")
        if self.include_max:
            max_5_bp = 0.5 * (sum_5_bp + absdiff_5_bp)
            max_5_3 = 0.5 * (sum_5_3 + absdiff_5_3)
            features.append(max_5_bp[:, np.newaxis])
            features.append(max_5_3[:, np.newaxis])

        return np.hstack(features)

    def get_feature_names_out(self, input_features=None):
        """
        Get output feature names for display.

        Args:
            input_features: Input feature names (ignored, we know it's [s5, sBP, s3])

        Returns:
            Array of output feature names (9 or 11 strings)
        """
        names = [
            's5',
            'sBP',
            's3',
            'min_5_bp',
            'min_5_3',
            'min_all',           # NEW
            'neg_absdiff_5_bp',  # NEW
            'neg_absdiff_5_3',   # NEW
            'neg_absdiff_bp_3'   # NEW
        ]

        if self.include_max:
            names.append('max_5_bp')
            names.append('max_5_3')

        return np.array(names)
