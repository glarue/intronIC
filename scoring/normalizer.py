"""
Score normalization with Zero-Anchored Robust (ZAR) scaling.

This module implements ZAR normalization for intron PWM log-odds ratio scores.
ZAR preserves the semantic zero of LLRs while correcting for dispersion shifts.

CRITICAL: This implementation:
1. Preserves semantic zero (zero = "U12 and U2 equally plausible")
2. Prevents data leakage through explicit dataset_type API
3. Enables cross-species deployment via unlabeled domain adaptation

Design Principles:
1. Zero-anchored: Scale by spread around zero, never center
2. Robust: Use median of absolute values with winsorization
3. Explicit dataset_type parameter forces conscious decision
4. Statistics frozen after fitting (immutable)

Expert Guidance:
For LLR features s = log(P(seq|PWM_U12)) - log(P(seq|PWM_U2)):
- Zero means "U12 ≈ U2" (preserve this!)
- Centering destroys semantic information
- Scale by median(|s|) after winsorizing at 99.5th percentile
"""

from typing import Literal, Optional, Iterable, Iterator
from dataclasses import replace
import numpy as np

from core.intron import Intron


# Type alias for dataset classification
DatasetType = Literal["reference", "experimental", "unlabeled"]


class ZeroAnchoredRobustScaler:
    """
    Zero-Anchored Robust (ZAR) scaler for log-odds ratio features.

    Scales features by robust spread around zero WITHOUT centering,
    preserving the semantic zero point of log-likelihood ratios.

    Standard approach (RobustScaler):
        z = (s - median) / IQR  ❌ Destroys semantic zero

    ZAR approach:
        z = s / median(|s|)     ✅ Preserves zero point

    The transformation:
    1. Winsorize |s| at high quantile (default: 99.5%) to reduce outlier impact
    2. Compute scale = median(winsorized |s|)
    3. Transform: z = s / scale

    Properties:
    - Zero stays zero (s=0 → z=0)
    - Sign preserved (s>0 → z>0, s<0 → z<0)
    - Robust to outliers (winsorization + median)
    - Minimal contamination from rare U12s (~0.5% of data)

    Attributes:
        scales_: Per-feature scale factors (shape: n_features)
        winsor_quantile_: Quantile for winsorization (default: 0.995)

    Example:
        >>> scaler = ZeroAnchoredRobustScaler()
        >>> X_scaled = scaler.fit_transform(X_train)
        >>> X_test_scaled = scaler.transform(X_test)
    """

    def __init__(self, winsor_quantile: float = 0.995):
        """
        Initialize ZAR scaler.

        Args:
            winsor_quantile: Quantile for winsorizing |s| (default: 0.995)
                            Clamps extreme values before computing median
        """
        self.winsor_quantile_ = winsor_quantile
        self.scales_: Optional[np.ndarray] = None

    def fit(self, X: np.ndarray) -> "ZeroAnchoredRobustScaler":
        """
        Fit ZAR scaler on training data.

        Args:
            X: Feature matrix of shape (n_samples, n_features)
               For intron scoring: (n_introns, 3) for [five, bp, three] LLRs

        Returns:
            self (for method chaining)

        Algorithm:
            For each feature d:
            1. Take absolute values: |s_d|
            2. Winsorize at quantile q: clip(|s_d|, 0, Q_q)
            3. Compute scale: σ_d = median(winsorized)
            4. Store for transform: s → s/σ_d
        """
        X = np.asarray(X)

        # Compute per-feature scales
        scales = []
        for feature_idx in range(X.shape[1]):
            feature_values = X[:, feature_idx]

            # Take absolute values
            abs_values = np.abs(feature_values)

            # Winsorize: clip at quantile to reduce outlier impact
            clip_value = np.quantile(abs_values, self.winsor_quantile_)
            winsorized = np.clip(abs_values, 0, clip_value)

            # Compute robust spread: median of winsorized absolute values
            scale = np.median(winsorized)

            # Handle edge case: all values near zero
            if scale == 0.0 or not np.isfinite(scale):
                scale = 1.0  # Fallback to no scaling

            scales.append(scale)

        self.scales_ = np.array(scales)
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        """
        Transform features using fitted scales.

        Args:
            X: Feature matrix of shape (n_samples, n_features)

        Returns:
            Scaled features: z = s / scale (zero preserved)

        Raises:
            RuntimeError: If fit() hasn't been called yet
        """
        if self.scales_ is None:
            raise RuntimeError("Must call fit() before transform()")

        X = np.asarray(X)

        # Divide by scales (element-wise per feature)
        # Broadcasting: (n_samples, n_features) / (n_features,)
        return X / self.scales_

    def fit_transform(self, X: np.ndarray) -> np.ndarray:
        """
        Fit and transform in one step.

        Args:
            X: Feature matrix

        Returns:
            Scaled features
        """
        return self.fit(X).transform(X)


class ScoreNormalizer:
    """
    Normalize PWM scores to z-scores with ML integrity guarantees.

    This class enforces ML best practices by preventing data leakage.
    The API is designed to make Issue #1 (post-classification re-normalization)
    impossible to replicate.

    Example Usage:
        >>> normalizer = ScoreNormalizer()
        >>> normalizer.fit(reference_introns, dataset_type="reference")  #  OK
        >>> normalized_refs = normalizer.transform(reference_introns, "reference")
        >>> normalized_exp = normalizer.transform(experimental_introns, "experimental")
        >>>
        >>> # This will raise ValueError:
        >>> normalizer.fit(experimental_introns, dataset_type="experimental")  # L Error!

    Attributes:
        _scaler: sklearn RobustScaler (None until fitted)
        _fitted_on: Which dataset type was used for fitting
        _is_fitted: Whether fit() has been called
    """

    def __init__(self):
        """Initialize an unfitted normalizer."""
        self._scaler: Optional[ZeroAnchoredRobustScaler] = None
        self._fitted_on: Optional[DatasetType] = None
        self._is_fitted: bool = False

    def fit(
        self,
        introns: Iterable[Intron],
        dataset_type: DatasetType = "reference"
    ) -> "ScoreNormalizer":
        """
        Fit scaler on intron scores.

        This method enforces ML best practices with three dataset types:
        - 'reference': Labeled training data (U12/U2 references) - standard case
        - 'unlabeled': Unlabeled experimental data for domain adaptation (cross-species)
        - 'experimental': FORBIDDEN - prevents post-classification re-normalization

        Args:
            introns: Introns with raw scores populated.
            dataset_type:
                - "reference": Fit on labeled training data (default)
                - "unlabeled": Fit on unlabeled experimental data for cross-species
                               domain adaptation. Valid because:
                               * No label leakage (labels not used)
                               * 99.5% U2 majority → robust estimators learn U2 baseline
                               * Corrects covariate shift from species differences
                - "experimental": FORBIDDEN - raises ValueError

        Returns:
            self (for method chaining)

        Raises:
            ValueError: If dataset_type is "experimental" (prevents Issue #1)
            ValueError: If introns list is empty
            ValueError: If introns have missing raw scores

        Cross-Species Use Case:
            When applying pretrained models to new species without curated references,
            fitting on unlabeled experimental data is statistically valid:
            - It's unsupervised domain adaptation (covariate-shift correction)
            - RobustScaler (median/IQR) is minimally affected by rare U12s (~0.5%)
            - Captures species-specific baseline for PWM score distributions
            - Allows pretrained SVM to detect U12s in correct coordinate system

        Port from: intronIC.py:3727-3731 (scale_scores - fitting part)
        DO NOT PORT: intronIC.py:5247-5251 (bad re-normalization)
        """
        # CRITICAL: Prevent Issue #1 by rejecting experimental data
        # (experimental = labeled data after classification)
        if dataset_type == "experimental":
            raise ValueError(
                "Cannot fit normalizer on experimental data! "
                "This would cause data leakage and invalidate ML evaluation. "
                "Use dataset_type='reference' for training/reference data, or "
                "dataset_type='unlabeled' for cross-species domain adaptation."
            )

        # Convert to list to check length and iterate multiple times
        intron_list = list(introns)

        # Validate input
        if len(intron_list) == 0:
            raise ValueError("Cannot fit on empty intron list")

        # Extract raw scores into matrix
        # Port from: intronIC.py:5696-5699 (get_score_vector)
        score_matrix = self._extract_score_matrix(intron_list)

        # Fit Zero-Anchored Robust (ZAR) scaler
        # Scales by median(|s|) after winsorization, preserving semantic zero
        self._scaler = ZeroAnchoredRobustScaler().fit(score_matrix)
        self._fitted_on = dataset_type
        self._is_fitted = True

        return self

    def transform(
        self,
        introns: Iterable[Intron],
        dataset_type: DatasetType
    ) -> Iterator[Intron]:
        """
        Transform intron scores to z-scores.

        Can be called on either reference or experimental data,
        but scaler must have been fit on reference data first.

        Args:
            introns: Introns with raw scores populated
            dataset_type: "reference" or "experimental" (for tracking only)

        Yields:
            Introns with z-scores populated (five_z_score, bp_z_score, three_z_score)

        Raises:
            RuntimeError: If fit() has not been called yet

        Port from: intronIC.py:3655-3665 (scale_scores - transformation part)
        """
        # Check if fitted
        if not self._is_fitted:
            raise RuntimeError(
                "Must call fit() before transform(). "
                "The normalizer needs to be fitted on reference data first."
            )

        # Extract raw scores into matrix
        # Port from: intronIC.py:5696-5699 (get_score_vector)
        intron_list = list(introns)
        score_matrix = self._extract_score_matrix(intron_list)

        # Transform using fitted scaler
        # Port from: intronIC.py:3662
        z_score_matrix = self._scaler.transform(score_matrix)

        # Yield introns with updated z-scores
        # Port from: intronIC.py:3663 (set_attributes) but using dataclasses.replace
        for intron, z_scores in zip(intron_list, z_score_matrix):
            yield self._update_intron_with_zscores(intron, z_scores)

    def fit_transform(
        self,
        introns: Iterable[Intron],
        dataset_type: DatasetType = "reference"
    ) -> Iterator[Intron]:
        """
        Convenience method: fit and transform in one step.

        Args:
            introns: Introns with raw scores populated
            dataset_type: Must be "reference" for fitting

        Yields:
            Introns with z-scores populated

        Raises:
            ValueError: If dataset_type is "experimental"
        """
        # Convert to list once (needed for both fit and transform)
        intron_list = list(introns)

        # Fit on the data
        self.fit(intron_list, dataset_type=dataset_type)

        # Transform and yield
        yield from self.transform(intron_list, dataset_type=dataset_type)

    # ========================================================================
    # Helper Methods
    # ========================================================================

    def _extract_score_matrix(self, introns: list[Intron]) -> np.ndarray:
        """
        Extract raw scores from introns into a numpy matrix.

        Port from: intronIC.py:5696-5699 (get_score_vector)

        Args:
            introns: List of introns with raw scores populated

        Returns:
            Numpy array of shape (n_introns, 3) with columns:
            [five_raw_score, bp_raw_score, three_raw_score]

        Raises:
            ValueError: If any intron has missing scores
            TypeError: If scores are None
        """
        score_matrix = []

        for intron in introns:
            # Check if scores exist
            if intron.scores is None:
                raise ValueError(
                    f"Intron {intron.intron_id} has no scores object. "
                    "Raw scores must be populated before normalization."
                )

            # Extract the three raw scores
            five_score = intron.scores.five_raw_score
            bp_score = intron.scores.bp_raw_score
            three_score = intron.scores.three_raw_score

            # Check for None values
            if five_score is None or bp_score is None or three_score is None:
                raise ValueError(
                    f"Intron {intron.intron_id} has missing raw scores: "
                    f"five={five_score}, bp={bp_score}, three={three_score}. "
                    "All raw scores must be populated before normalization."
                )

            # Add to matrix
            score_matrix.append([five_score, bp_score, three_score])

        # Convert to numpy array
        return np.asarray(score_matrix)

    def _update_intron_with_zscores(
        self,
        intron: Intron,
        z_scores: np.ndarray
    ) -> Intron:
        """
        Create a new intron with updated z-scores.

        Port from: intronIC.py:3644-3652 (set_attributes) but using immutable approach

        Args:
            intron: Original intron
            z_scores: Array of [five_z_score, bp_z_score, three_z_score]

        Returns:
            New Intron with z-scores populated
        """
        # Extract z-scores
        five_z, bp_z, three_z = z_scores

        # Create new scores object with z-scores populated
        # Port from: intronIC.py:3648-3649
        updated_scores = replace(
            intron.scores,
            five_z_score=float(five_z),
            bp_z_score=float(bp_z),
            three_z_score=float(three_z)
        )

        # Return new intron with updated scores
        return replace(intron, scores=updated_scores)
