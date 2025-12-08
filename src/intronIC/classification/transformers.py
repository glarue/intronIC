"""
Custom sklearn transformers for feature augmentation.

This module provides transformers for augmenting base features with
composite features that help the linear SVM better separate U12 from U2.
"""

from typing import Optional

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
        include_pairwise_mins: Whether to include pairwise min features (default: False)
                    If False, only min_all is included (recommended - non-redundant)
                    If True, includes min_5_bp, min_5_3, AND min_all (backward compatible)

    Input features (3D):
        - s5:  5' splice site z-score (LLR, zero-anchored)
        - sBP: Branch point z-score (LLR, zero-anchored)
        - s3:  3' splice site z-score (LLR, zero-anchored)

    Output features (varies by configuration):
        Base configuration (7D, recommended: include_pairwise_mins=False, include_max=False):
        1. s5 (original, passed through)
        2. sBP (original, passed through)
        3. s3 (original, passed through)
        4. min_all = min(s5, sBP, s3)           ← ALL THREE must be strong (3-way AND)
        5. neg_absdiff_5_bp = -|s5 - sBP|       ← Penalty for 5'/BP imbalance
        6. neg_absdiff_5_3 = -|s5 - s3|         ← Penalty for 5'/3' imbalance
        7. neg_absdiff_bp_3 = -|sBP - s3|       ← Penalty for BP/3' imbalance

        With include_pairwise_mins=True (9D, adds redundant features):
        4. min_5_bp = min(s5, sBP)              ← Both 5' and BP must be strong (redundant)
        5. min_5_3 = min(s5, s3)                ← Both 5' and 3' must be strong (redundant)
        6. min_all = min(s5, sBP, s3)           ← ALL THREE must be strong
        7-9. [neg_absdiff features as above]

        With include_max=True (adds 2 features):
        - max_5_bp = max(s5, sBP)               ← At least one of 5'/BP strong (optional)
        - max_5_3 = max(s5, s3)                 ← At least one of 5'/3' strong (optional)

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
        include_pairwise_mins: bool = False,
        features: Optional[list] = ['absdiff_bp_3'],  # Default: minimal 4D set from L1 analysis
        gamma_imbalance: float = 1.0,
        gamma_5_bp=None,  # Deprecated - kept for backward compatibility
        gamma_5_3=None    # Deprecated - kept for backward compatibility
    ):
        """
        Initialize BothEndsStrongTransformer.

        Args:
            include_max: Whether to include max features (default: False)
                        DEPRECATED: Use features parameter instead for fine-grained control
            include_pairwise_mins: Whether to include pairwise min features (default: False)
                        DEPRECATED: Use features parameter instead for fine-grained control
            features: List of composite feature names to include
                     DEFAULT: ['absdiff_bp_3'] - minimal 4D set identified by L1 regularization
                     Set to None for legacy 7D behavior (backward compatibility with old models)
                     Available features: 'min_5_bp', 'min_5_3', 'min_all',
                                        'absdiff_5_bp', 'absdiff_5_3', 'absdiff_bp_3',
                                        'max_5_bp', 'max_5_3'
                     Deprecated (backward compat): 'neg_absdiff_5_bp', 'neg_absdiff_5_3', 'neg_absdiff_bp_3'
                     Example: ['absdiff_bp_3'] for minimal 4D feature space (DEFAULT)
                     Example: ['min_all', 'absdiff_5_bp', 'absdiff_5_3', 'absdiff_bp_3'] for 7D
            gamma_imbalance: Scaling factor for imbalance features (default: 1.0)
                           Multiplies all absdiff_* features to increase penalty for imbalance.
                           Higher values (e.g., 2.0, 4.0) make imbalance more costly.
                           Use with grid search to find optimal scaling.
            gamma_5_bp: DEPRECATED - Ignored, kept for backward compatibility with old models
            gamma_5_3: DEPRECATED - Ignored, kept for backward compatibility with old models
        """
        # Store original flags for backward compatibility with pickled models
        self.include_max = include_max
        self.include_pairwise_mins = include_pairwise_mins
        self.gamma_imbalance = gamma_imbalance

        # If features is explicitly None, use the default
        if features is None:
            features = ['absdiff_bp_3']

        # Validate feature names
        valid_features = {
            'min_5_bp', 'min_5_3', 'min_all',
            'absdiff_5_bp', 'absdiff_5_3', 'absdiff_bp_3',
            'max_5_bp', 'max_5_3'
        }
        invalid = set(features) - valid_features
        if invalid:
            raise ValueError(
                f"Invalid feature names: {invalid}. "
                f"Valid features: {valid_features}"
            )

        self.features = features

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
        Augment 3D features to 7D, 9D, or 11D depending on configuration.

        Args:
            X: Array of shape (n_samples, 3) with [s5, sBP, s3]

        Returns:
            Array of shape (n_samples, N) where N depends on configuration:
            - 7D: include_pairwise_mins=False, include_max=False (recommended)
                  [s5, sBP, s3, min_all, neg_absdiff_5_bp, neg_absdiff_5_3, neg_absdiff_bp_3]
            - 9D: include_pairwise_mins=True, include_max=False
                  [s5, sBP, s3, min_5_bp, min_5_3, min_all, neg_absdiff_5_bp, neg_absdiff_5_3, neg_absdiff_bp_3]
            - 9D: include_pairwise_mins=False, include_max=True
                  [s5, sBP, s3, min_all, neg_absdiff_5_bp, neg_absdiff_5_3, neg_absdiff_bp_3, max_5_bp, max_5_3]
            - 11D: include_pairwise_mins=True, include_max=True
                   [s5, sBP, s3, min_5_bp, min_5_3, min_all, neg_absdiff_5_bp, neg_absdiff_5_3, neg_absdiff_bp_3, max_5_bp, max_5_3]

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

        # Imbalance penalties (CORRECTED 2025-01-20)
        # Absolute difference = penalty for inconsistency
        # Model should learn negative weights → penalizes imbalanced signals
        # (Previously used neg_absdiff which caused wrong sign on coefficients)
        absdiff_5_bp_feat = absdiff_5_bp
        absdiff_5_3_feat = absdiff_5_3
        absdiff_bp_3 = np.abs(sBP - s3)

        # Stack base features (always included)
        features = [
            s5[:, np.newaxis],
            sBP[:, np.newaxis],
            s3[:, np.newaxis]
        ]

        # Use explicit feature list
        # Compute max features if needed
        if 'max_5_bp' in self.features or 'max_5_3' in self.features:
            max_5_bp = 0.5 * (sum_5_bp + absdiff_5_bp)
            max_5_3 = 0.5 * (sum_5_3 + absdiff_5_3)

        # Add composite features in deterministic order
        # This ensures consistent feature ordering regardless of input list order
        # Apply gamma scaling to imbalance features (absdiff_*)
        composite_feature_map = {
            'min_5_bp': min_5_bp[:, np.newaxis],
            'min_5_3': min_5_3[:, np.newaxis],
            'min_all': min_all[:, np.newaxis],
            'absdiff_5_bp': self.gamma_imbalance * absdiff_5_bp_feat[:, np.newaxis],
            'absdiff_5_3': self.gamma_imbalance * absdiff_5_3_feat[:, np.newaxis],
            'absdiff_bp_3': self.gamma_imbalance * absdiff_bp_3[:, np.newaxis],
        }
        # Add max features if computed
        if 'max_5_bp' in self.features or 'max_5_3' in self.features:
            composite_feature_map['max_5_bp'] = max_5_bp[:, np.newaxis]
            composite_feature_map['max_5_3'] = max_5_3[:, np.newaxis]

        # Add features in canonical order (not self.features order)
        for feature_name in ['min_5_bp', 'min_5_3', 'min_all',
                             'absdiff_5_bp', 'absdiff_5_3', 'absdiff_bp_3',
                             'max_5_bp', 'max_5_3']:
            if feature_name in self.features:
                features.append(composite_feature_map[feature_name])

        return np.hstack(features)

    def get_feature_names_out(self, input_features=None):
        """
        Get output feature names for display.

        Args:
            input_features: Input feature names (ignored, we know it's [s5, sBP, s3])

        Returns:
            Array of output feature names (varies depending on config)
        """
        names = [
            's5',
            'sBP',
            's3'
        ]

        # Handle old pickled models that don't have 'features' attribute
        if hasattr(self, 'features') and self.features is not None:
            # Use explicit feature list in canonical order
            # Support both new (absdiff_*) and deprecated (neg_absdiff_*) names
            for feature_name in ['min_5_bp', 'min_5_3', 'min_all',
                                 'absdiff_5_bp', 'absdiff_5_3', 'absdiff_bp_3',
                                 'neg_absdiff_5_bp', 'neg_absdiff_5_3', 'neg_absdiff_bp_3',
                                 'max_5_bp', 'max_5_3']:
                if feature_name in self.features:
                    names.append(feature_name)
        else:
            # Legacy flag-based behavior
            # Optional: pairwise min features
            if self.include_pairwise_mins:
                names.extend(['min_5_bp', 'min_5_3'])

            # Always include: 3-way min and imbalance penalties
            names.extend([
                'min_all',
                'neg_absdiff_5_bp',
                'neg_absdiff_5_3',
                'neg_absdiff_bp_3'
            ])

            # Optional: max features
            if self.include_max:
                names.extend(['max_5_bp', 'max_5_3'])

        return np.array(names)
