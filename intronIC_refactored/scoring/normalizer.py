"""
Score normalization with ML integrity guarantees.

This module implements z-score normalization for intron PWM scores while
preventing data leakage through API design.

CRITICAL: This implementation fixes Issue #1 from SCORING_ANALYSIS.md by
making it impossible to accidentally fit the normalizer on experimental data.

Design Principles:
1. Explicit dataset_type parameter forces conscious decision
2. Fitting on experimental data raises ValueError (fail-fast)
3. Statistics are frozen after fitting (immutable)
4. Clear separation between reference and experimental data
"""

from typing import Literal, Optional, Iterable, Iterator
from dataclasses import replace
import numpy as np
from sklearn.preprocessing import RobustScaler

from core.intron import Intron


# Type alias for dataset classification
DatasetType = Literal["reference", "experimental"]


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
        self._scaler: Optional[RobustScaler] = None
        self._fitted_on: Optional[DatasetType] = None
        self._is_fitted: bool = False

    def fit(
        self,
        introns: Iterable[Intron],
        dataset_type: DatasetType = "reference"
    ) -> "ScoreNormalizer":
        """
        Fit scaler on intron scores.

        CRITICAL: This method enforces ML best practices by only allowing
        fitting on reference data. Attempting to fit on experimental data
        will raise ValueError.

        Args:
            introns: Introns with raw scores populated. Should be REFERENCE
                    introns only for proper ML practice.
            dataset_type: Must be "reference" to prevent accidental misuse.
                         If "experimental", raises ValueError.

        Returns:
            self (for method chaining)

        Raises:
            ValueError: If dataset_type is "experimental" (prevents Issue #1)
            ValueError: If introns list is empty
            ValueError: If introns have missing raw scores

        Port from: intronIC.py:3727-3731 (scale_scores - fitting part)
        DO NOT PORT: intronIC.py:5247-5251 (bad re-normalization)
        """
        # CRITICAL: Prevent Issue #1 by rejecting experimental data
        if dataset_type == "experimental":
            raise ValueError(
                "Cannot fit normalizer on experimental data! "
                "This would cause data leakage and invalidate ML evaluation. "
                "Use dataset_type='reference' for training/reference data only."
            )

        # Convert to list to check length and iterate multiple times
        intron_list = list(introns)

        # Validate input
        if len(intron_list) == 0:
            raise ValueError("Cannot fit on empty intron list")

        # Extract raw scores into matrix
        # Port from: intronIC.py:5696-5699 (get_score_vector)
        score_matrix = self._extract_score_matrix(intron_list)

        # Fit RobustScaler (median/IQR normalization)
        # Uses middle 50% for scaling, giving better differentiation for tail values (U12s)
        self._scaler = RobustScaler().fit(score_matrix)
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
