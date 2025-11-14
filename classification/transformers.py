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
    both-ends-strong features.

    This transformer builds composite features that reward joint strength
    (sum) and penalize imbalance (abs difference) for correlated regions.

    The γ (gamma) parameters control how strongly to weight the imbalance
    penalties. Higher γ makes the model more sensitive to one-end-strong
    patterns by changing the effective L2 regularization for that feature.

    Port from: Expert recommendations for reducing false positives
               See: BOTHENDS_IMPLEMENTATION_PLAN.md

    Background:
        U12-type introns have strong correlation between 5'SS and BPS motifs.
        False positives often show "one-end-strong" patterns: one very high
        score compensating for weak/negative scores elsewhere.

        By explicitly adding asymmetry penalties as features (rather than
        relying on the linear model to learn this implicitly), we make the
        classification problem more linearly separable.

    How γ works:
        γ is NOT a loss weight - it's a feature-scaling hyperparameter.

        - Larger γ makes the imbalance feature larger
        - L2 regularization makes it "cheaper" to use larger features
        - Model learns to heavily penalize imbalance
        - Tuned via GridSearchCV to find optimal balance

    Attributes:
        gamma_5_bp: Weight for 5'SS-BPS imbalance penalty (default: 1.0)
                   Suggested range for tuning: [1, 2, 4, 8]
                   Higher values emphasize 5'SS-BPS correlation
        gamma_5_3: Weight for 5'SS-3'SS imbalance penalty (default: 1.0)
                  Suggested range for tuning: [1, 2, 4]
        include_min: Whether to include min(s5, sBP, s3) feature (default: False)
        include_total_imbalance: Whether to include total imbalance (default: False)

    Input features (3D):
        - s5:  5' splice site z-score (LLR, zero-anchored)
        - sBP: Branch point z-score (LLR, zero-anchored)
        - s3:  3' splice site z-score (LLR, zero-anchored)

    Output features (7D by default, or 9D with optional features):
        1. s5 (original, passed through)
        2. sBP (original, passed through)
        3. s3 (original, passed through)
        4. sum_5_bp = s5 + sBP
        5. absdiff_5_bp = |s5 - sBP| × gamma_5_bp  ← Scaled penalty
        6. sum_5_3 = s5 + s3
        7. absdiff_5_3 = |s5 - s3| × gamma_5_3    ← Scaled penalty
        [8. min_all = min(s5, sBP, s3)]            ← Optional
        [9. imbalance_all = |s5-sBP| + |s5-s3| + |sBP-s3|]  ← Optional

    Mathematical insight:
        A linear model on (sum, absdiff) can emulate an AND gate:
            min(a, b) = ((a + b) - |a - b|) / 2

        So features (sum_5_bp, absdiff_5_bp) let the model learn:
            "Both 5'SS AND BPS must be strong for U12"

    Example:
        >>> transformer = BothEndsStrongTransformer(gamma_5_bp=4, gamma_5_3=2)
        >>> X = np.array([[2.0, 1.5, -0.5]])  # s5=2, sBP=1.5, s3=-0.5
        >>> X_aug = transformer.transform(X)
        >>> # X_aug = [2.0, 1.5, -0.5, 3.5, 2.0, 1.5, 5.0]
        >>> #          [s5,  sBP, s3,  sum, |diff|×4, sum, |diff|×2]

        One-end-strong example (should be penalized):
        >>> X_bad = np.array([[-1.0, -0.5, 10.0]])  # Only 3'SS strong
        >>> X_bad_aug = transformer.transform(X_bad)
        >>> # absdiff_5_3 = |-1.0 - 10.0| × 2 = 22.0  ← Large penalty!
    """

    def __init__(
        self,
        gamma_5_bp: float = 1.0,
        gamma_5_3: float = 1.0,
        include_min: bool = False,
        include_total_imbalance: bool = False
    ):
        """
        Initialize BothEndsStrongTransformer.

        Args:
            gamma_5_bp: Weight for 5'SS-BPS imbalance penalty
            gamma_5_3: Weight for 5'SS-3'SS imbalance penalty
            include_min: Whether to add min(s5, sBP, s3) feature
            include_total_imbalance: Whether to add total imbalance feature
        """
        self.gamma_5_bp = gamma_5_bp
        self.gamma_5_3 = gamma_5_3
        self.include_min = include_min
        self.include_total_imbalance = include_total_imbalance

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
        Augment 3D features to 7+D.

        Args:
            X: Array of shape (n_samples, 3) with [s5, sBP, s3]

        Returns:
            Array of shape (n_samples, 7+) with augmented features

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

        # Build pairwise both-ends-strong features
        # For 5'SS-BPS (strongest correlation)
        sum_5_bp = s5 + sBP
        absdiff_5_bp = np.abs(s5 - sBP) * self.gamma_5_bp

        # For 5'SS-3'SS (secondary correlation)
        sum_5_3 = s5 + s3
        absdiff_5_3 = np.abs(s5 - s3) * self.gamma_5_3

        # Stack base + pairwise features
        features = [
            s5[:, np.newaxis],
            sBP[:, np.newaxis],
            s3[:, np.newaxis],
            sum_5_bp[:, np.newaxis],
            absdiff_5_bp[:, np.newaxis],
            sum_5_3[:, np.newaxis],
            absdiff_5_3[:, np.newaxis]
        ]

        # Optional: min of all three (all-three-strong feature)
        if self.include_min:
            min_all = np.minimum(np.minimum(s5, sBP), s3)
            features.append(min_all[:, np.newaxis])

        # Optional: total imbalance (sum of all pairwise differences)
        if self.include_total_imbalance:
            imbalance_all = (
                np.abs(s5 - sBP) +
                np.abs(s5 - s3) +
                np.abs(sBP - s3)
            )
            features.append(imbalance_all[:, np.newaxis])

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
            'sum_5_bp',
            'absdiff_5_bp',
            'sum_5_3',
            'absdiff_5_3'
        ]

        if self.include_min:
            names.append('min_all')

        if self.include_total_imbalance:
            names.append('imbalance_all')

        return np.array(names)
