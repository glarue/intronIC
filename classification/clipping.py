"""
Outlier clipping transformer for robust feature scaling.

Clips extreme values at high quantiles to prevent rare outliers from
dominating the linear SVM decision boundary and calibration.

Expert guidance: "Rare huge LLRs (e.g., 5' outliers) can dominate a linear
margin and calibration."

Port from: FP_REDUCTION_IMPLEMENTATION.md (Tier 2, Step 5)
"""

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin


class OutlierClipper(BaseEstimator, TransformerMixin):
    """
    Clip extreme outlier values at specified quantiles.

    This transformer learns symmetric clipping thresholds from training data
    and applies them to all data. Helps prevent rare huge LLR scores from
    dominating the linear model.

    Expert guidance: "Rare huge LLRs (e.g., 5' outliers) can dominate a
    linear margin and calibration. Winsorize each base LLR at a high quantile
    learned on human training (e.g., Q_0.999)."

    The clipping is:
    - Symmetric: Uses quantile of |X| to get cap, then clips to [-cap, +cap]
    - Per-feature: Each feature (5'SS, BPS, 3'SS) gets independent threshold
    - Conservative: Default Q_0.999 only clips extreme 0.1% outliers
    - Zero-preserving: Does not affect the semantic zero point

    Attributes:
        quantile: Quantile for symmetric clipping (default: 0.999)
                 Values clipped to [-cap, +cap] where cap = quantile of |X|
        caps_: Learned clipping thresholds per feature (set during fit)
               Shape: (n_features,)

    Example:
        >>> clipper = OutlierClipper(quantile=0.999)
        >>> # Fit on reference introns only (prevent data leakage)
        >>> clipper.fit(X_reference)
        >>> # Apply same caps to reference and experimental
        >>> X_reference_clipped = clipper.transform(X_reference)
        >>> X_experimental_clipped = clipper.transform(X_experimental)

    Example (manual inspection):
        >>> X = np.array([
        ...     [1.0, 2.0, 1.5],
        ...     [1.2, 2.1, 1.4],
        ...     [1.1, 2.0, 1.6],
        ...     [10.0, 2.2, 1.5]  # Outlier in first feature
        ... ])
        >>> clipper = OutlierClipper(quantile=0.95)
        >>> clipper.fit(X)
        >>> X_clipped = clipper.transform(X)
        >>> # caps_ ≈ [10.0, 2.2, 1.6]  (95th percentile of abs values)
        >>> # Outlier (10.0) gets clipped to cap, others pass through

    Port from: intronIC expert recommendations for FP reduction
    """

    def __init__(self, quantile: float = 0.999):
        """
        Initialize OutlierClipper.

        Args:
            quantile: Quantile for clipping (0.999 = clip at 99.9th percentile)
                     Higher values = less aggressive clipping
                     Expert recommendation: 0.999 (clips only extreme 0.1%)

        Raises:
            ValueError: If quantile not in (0, 1)
        """
        if not 0 < quantile < 1:
            raise ValueError(f"quantile must be in (0, 1), got {quantile}")
        self.quantile = quantile
        self.caps_ = None

    def fit(self, X, y=None):
        """
        Learn clipping thresholds from training data.

        Computes symmetric caps at the specified quantile of absolute values
        for each feature independently. This preserves the zero-anchoring of
        LLR scores (s=0 means "U12≈U2").

        Args:
            X: Training data of shape (n_samples, n_features)
               For intronIC: typically (n_introns, 3) for [s5, sBP, s3]
            y: Target values (ignored, for sklearn compatibility)

        Returns:
            self

        Example:
            >>> X = np.array([[1.0, 2.0], [1.5, 3.0], [10.0, 2.5]])
            >>> clipper = OutlierClipper(quantile=0.95)
            >>> clipper.fit(X)
            >>> clipper.caps_  # [10.0, 3.0] (95th percentile of |X|)
        """
        X = np.asarray(X)

        # Compute symmetric caps: quantile of |X| per feature
        # Shape: (n_features,)
        # Uses absolute values to ensure symmetric clipping
        self.caps_ = np.quantile(np.abs(X), self.quantile, axis=0)

        return self

    def transform(self, X):
        """
        Apply learned clipping thresholds to data.

        Clips each feature to [-cap, +cap] where cap was learned during fit.
        This is applied element-wise per feature, preserving the zero point.

        Args:
            X: Data to clip of shape (n_samples, n_features)

        Returns:
            X_clipped: Clipped data of same shape
                      Values outside [-cap_i, +cap_i] are clamped to bounds

        Raises:
            ValueError: If transform called before fit

        Example:
            >>> X = np.array([[1.0, 2.0], [1.5, 3.0], [10.0, 2.5]])
            >>> clipper = OutlierClipper(quantile=0.90)
            >>> clipper.fit(X[[0, 1]])  # Fit on subset
            >>> X_clipped = clipper.transform(X)
            >>> # Outlier (10.0) gets clipped to learned cap
        """
        if self.caps_ is None:
            raise ValueError("OutlierClipper must be fit before transform")

        X = np.asarray(X)

        # Clip each feature independently
        X_clipped = np.copy(X)
        for i, cap in enumerate(self.caps_):
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
            >>> clipper = OutlierClipper()
            >>> clipper.caps_ = np.array([10.0, 5.0, 3.0])
            >>> clipper.get_feature_names_out(['s5', 'sBP', 's3'])
            array(['s5', 'sBP', 's3'], dtype='<U3')
        """
        if input_features is None:
            # If not provided, use generic names
            n_features = len(self.caps_) if self.caps_ is not None else 0
            return np.array([f'x{i}' for i in range(n_features)])
        return np.asarray(input_features)
