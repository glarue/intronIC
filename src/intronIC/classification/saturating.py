"""
Saturating transformation for extreme z-scores.

Applies a log-based compression to prevent extreme outliers from
dominating the linear SVM, while preserving ordering and sign.

Expert guidance: "Optional saturating transform (if needed): After clipping
(or instead of very large Z_max), you can use f(z) = sign(z) * log(1 + |z|).
This preserves ordering but compresses huge values—so 11.45 vs 4 is much less
different in the model space."

Redesign: SCALER_ARCHITECTURE_REVIEW.md (Expert feedback 2025)
"""

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin


class SaturatingTransform(BaseEstimator, TransformerMixin):
    """
    Optional saturating transformation for z-scores.

    Compresses extreme values using f(z) = sign(z) * log(1 + |z|) to prevent
    outliers from dominating the linear model. This is applied AFTER symmetric
    clipping as an additional layer of outlier protection.

    Pipeline position:
        ... → SymmetricClipper → SaturatingTransform → BothEndsStrong → SVC

    Expert guidance: "Preserves ordering but compresses huge values—so 11.45
    vs 4 is much less different in the model space."

    Use Case:
        When SymmetricClipper alone isn't sufficient (e.g., clipping at 0.99
        allows z=5 to pass through), saturation further compresses extreme
        values so they can't overwhelm balance features.

    Transformation:
        f(z) = sign(z) * log(1 + |z|)

    Properties:
        - Preserves sign: negative z → negative f(z)
        - Preserves ordering: if z1 < z2, then f(z1) < f(z2)
        - Preserves zero: f(0) = 0
        - Compresses extremes:
            z=0 → 0
            z=1 → 0.69
            z=2 → 1.10
            z=4 → 1.61  (clipped typical max)
            z=11.45 → 2.52  (C. elegans extreme artifact)
        - Difference compression:
            Before: 11.45 - 4 = 7.45 (huge gap)
            After:  2.52 - 1.61 = 0.91 (moderate gap)

    Hyperparameter:
        enabled: Boolean, controlled via grid search
            - False: Identity transform (pass through unchanged)
            - True: Apply saturation

    Attributes:
        enabled: Whether to apply transformation (default: True)

    Example (Typical use):
        >>> # After clipping at quantile=0.975
        >>> z_clipped = np.array([[0.5, 1.8, -0.3],   # Normal intron
        ...                       [2.0, 1.5,  3.8]])  # At clip boundary
        >>> saturate = SaturatingTransform(enabled=True)
        >>> z_sat = saturate.fit_transform(z_clipped)
        >>> z_sat  # [[0.41, 1.03, -0.26], [1.10, 0.92, 1.57]]
        >>> # Extreme values compressed, typical values mostly preserved

    Example (Cross-species artifact):
        >>> # C. elegans with extreme value (after clipping to cap=4)
        >>> z_clipped = np.array([[-0.15, 0.30, 4.0]])  # Clipped from 11.45
        >>> saturate = SaturatingTransform(enabled=True)
        >>> z_sat = saturate.fit_transform(z_clipped)
        >>> z_sat  # [[-0.14, 0.26, 1.61]]
        >>> # 4.0 → 1.61 (compressed), preventing it from overwhelming penalties

    Example (Disabled - identity):
        >>> z = np.array([[1.5, -2.0, 0.5]])
        >>> saturate = SaturatingTransform(enabled=False)
        >>> z_out = saturate.fit_transform(z)
        >>> np.array_equal(z, z_out)  # True (pass through)

    Redesign: SCALER_ARCHITECTURE_REVIEW.md (Expert feedback, 2025)
    """

    def __init__(self, enabled: bool = True):
        """
        Initialize SaturatingTransform.

        Args:
            enabled: Whether to apply saturation (default: True)
                    - True: Apply f(z) = sign(z) * log(1 + |z|)
                    - False: Identity transform (pass through)
                    Grid search recommendation: [False, True]
        """
        self.enabled = enabled

    def fit(self, X, y=None):
        """
        Fit the transformer (no-op, stateless transform).

        This transformer has no parameters to learn. The fit method is
        provided only for sklearn Pipeline compatibility.

        Args:
            X: Training data (ignored)
            y: Target values (ignored)

        Returns:
            self
        """
        # Stateless transform, nothing to fit
        return self

    def transform(self, X):
        """
        Apply saturating transformation (or pass through if disabled).

        Args:
            X: Z-scores to transform of shape (n_samples, n_features)

        Returns:
            X_transformed: Saturated z-scores (or unchanged if disabled)

        Algorithm:
            if enabled:
                f(z) = sign(z) * log(1 + |z|)
            else:
                f(z) = z  (identity)

        Example:
            >>> z = np.array([[0, 1, -2, 4, 11.45]])
            >>> saturate = SaturatingTransform(enabled=True)
            >>> saturate.fit_transform(z).round(2)
            array([[ 0.  ,  0.69, -1.1 ,  1.61,  2.52]])
        """
        X = np.asarray(X)

        if not self.enabled:
            # Identity transform (pass through)
            return X.copy()

        # Apply saturating transform: f(z) = sign(z) * log(1 + |z|)
        # Using np.sign and np.log1p (log1p(x) = log(1+x), more numerically stable)
        return np.sign(X) * np.log1p(np.abs(X))

    def inverse_transform(self, X):
        """
        Inverse saturating transformation (for interpretability).

        This can be used to map saturated values back to z-score space
        for interpretation, though it's rarely needed in practice.

        Args:
            X: Saturated values

        Returns:
            Z: Original z-scores (approximately)

        Algorithm:
            if enabled:
                z = sign(X) * (exp(|X|) - 1)
            else:
                z = X  (identity)

        Example:
            >>> saturate = SaturatingTransform(enabled=True)
            >>> z_original = np.array([[4.0, -2.0, 11.45]])
            >>> z_sat = saturate.fit_transform(z_original)
            >>> z_recovered = saturate.inverse_transform(z_sat)
            >>> np.allclose(z_original, z_recovered)  # True
        """
        X = np.asarray(X)

        if not self.enabled:
            # Identity transform
            return X.copy()

        # Inverse: z = sign(X) * (exp(|X|) - 1)
        # Using np.expm1 (expm1(x) = exp(x)-1, more numerically stable)
        return np.sign(X) * np.expm1(np.abs(X))

    def get_feature_names_out(self, input_features=None):
        """
        Get output feature names (same as input).

        Saturation doesn't change feature names, only their scale.

        Args:
            input_features: Input feature names (optional)

        Returns:
            Array of output feature names (same as input)

        Example:
            >>> saturate = SaturatingTransform()
            >>> saturate.get_feature_names_out(['z5', 'zBP', 'z3'])
            array(['z5', 'zBP', 'z3'], dtype='<U3')
        """
        if input_features is None:
            # Generic names (won't know n_features without fitting)
            return np.array([])
        return np.asarray(input_features)
