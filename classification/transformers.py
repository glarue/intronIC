"""
Custom sklearn transformers for feature augmentation.

This module provides transformers for augmenting base features with
composite features that help the linear SVM better separate U12 from U2.
"""

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin


class BothEndsStrongTransformer(BaseEstimator, TransformerMixin):
    """
    Augment 3D features (5'SS, BPS, 3'SS z-scores) with pairwise
    min/max features to capture "both ends strong" patterns.

    This transformer builds composite features using min and max operations
    to directly express "both must be strong" (min) and "at least one is strong" (max).

    Port from: Expert recommendations for reducing false positives
               See: BOTHENDS_IMPLEMENTATION_PLAN.md

    Background:
        U12-type introns have strong correlation between 5'SS and BPS motifs.
        False positives often show "one-end-strong" patterns: one very high
        score compensating for weak/negative scores elsewhere.

        By explicitly computing min/max features, we make it easy for the
        linear model to express "both 5' AND BP must be strong" by simply
        weighting min_5_bp positively.

    Mathematical insight:
        min(a, b) = 0.5 * ((a + b) - |a - b|)
        max(a, b) = 0.5 * ((a + b) + |a - b|)

        The linear model will primarily weight min features (expert guidance),
        since min directly captures "both must be strong for U12".

    Attributes:
        include_max: Whether to include max features (default: False)
                    If True, adds max_5_bp and max_5_3 features
                    Expert notes: "model will mostly weight min"

    Input features (3D):
        - s5:  5' splice site z-score (LLR, zero-anchored)
        - sBP: Branch point z-score (LLR, zero-anchored)
        - s3:  3' splice site z-score (LLR, zero-anchored)

    Output features (7D default, 9D with include_max=True):
        1. s5 (original, passed through)
        2. sBP (original, passed through)
        3. s3 (original, passed through)
        4. min_5_bp = 0.5 * ((s5 + sBP) - |s5 - sBP|)  ← Both must be strong
        5. min_5_3 = 0.5 * ((s5 + s3) - |s5 - s3|)     ← Both must be strong
        [6. max_5_bp = 0.5 * ((s5 + sBP) + |s5 - sBP|)]  ← At least one strong (optional)
        [7. max_5_3 = 0.5 * ((s5 + s3) + |s5 - s3|)]     ← At least one strong (optional)

    Example:
        >>> transformer = BothEndsStrongTransformer(include_max=False)
        >>> X = np.array([[2.0, 1.5, -0.5]])  # s5=2, sBP=1.5, s3=-0.5
        >>> X_aug = transformer.transform(X)
        >>> # X_aug = [2.0, 1.5, -0.5, 1.75, -1.25]
        >>> #          [s5,  sBP, s3,  min_5_bp, min_5_3]
        >>> # min_5_bp = 0.5 * ((2.0 + 1.5) - |2.0 - 1.5|) = 0.5 * 3.25 = 1.625
        >>> # min_5_3 = 0.5 * ((2.0 + (-0.5)) - |2.0 - (-0.5)|) = 0.5 * (-1.0) = -0.5

        One-end-strong example (should have low min):
        >>> X_bad = np.array([[-1.0, -0.5, 10.0]])  # Only 3'SS strong
        >>> X_bad_aug = transformer.transform(X_bad)
        >>> # min_5_3 = 0.5 * ((-1.0 + 10.0) - |-1.0 - 10.0|) = 0.5 * (9.0 - 11.0) = -1.0
        >>> # Low min indicates one-end-strong pattern!
    """

    def __init__(
        self,
        include_max: bool = False
    ):
        """
        Initialize BothEndsStrongTransformer.

        Args:
            include_max: Whether to include max features (default: False)
                        Expert guidance: "model will mostly weight min"
        """
        self.include_max = include_max

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
        Augment 3D features to 5D (or 7D with include_max=True).

        Args:
            X: Array of shape (n_samples, 3) with [s5, sBP, s3]

        Returns:
            Array of shape (n_samples, 5) or (n_samples, 7) with augmented features

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

        # Stack base + min features
        features = [
            s5[:, np.newaxis],
            sBP[:, np.newaxis],
            s3[:, np.newaxis],
            min_5_bp[:, np.newaxis],
            min_5_3[:, np.newaxis]
        ]

        # Optional: max features (expert notes: "model will mostly weight min")
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
            Array of output feature names
        """
        names = [
            's5',
            'sBP',
            's3',
            'min_5_bp',
            'min_5_3'
        ]

        if self.include_max:
            names.append('max_5_bp')
            names.append('max_5_3')

        return np.array(names)
