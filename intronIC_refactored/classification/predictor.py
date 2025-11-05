"""
SVM ensemble prediction with F1-weighted averaging.

This module implements the classification algorithm from intronIC.py:5651-5900,
which applies trained ensemble models to classify introns as U2 or U12 type.

Key features:
- F1-weighted averaging of ensemble predictions for robust classification
- Decision boundary distance calculation for confidence estimation
- Type assignment based on probability threshold (default: 90%)
- No re-normalization of z-scores (prevents data leakage)

Port from: intronIC.py:5651-5687 (average_svm_score_info), 5816-5900 (parallel_svm_score)
"""

from dataclasses import replace
from typing import Sequence
import numpy as np

from core.intron import Intron, IntronScores
from classification.trainer import SVMEnsemble


class SVMPredictor:
    """
    Apply trained SVM ensemble to classify introns.

    Uses F1-weighted averaging across models for robust predictions.
    Does NOT re-normalize z-scores (fixes Issue #1 - data leakage).

    Port from: intronIC.py:5651-5900
    """

    def __init__(self, threshold: float = 90.0):
        """
        Initialize predictor.

        Args:
            threshold: U12 probability threshold (0-100, default: 90)
        """
        if not 0 <= threshold <= 100:
            raise ValueError(f"Threshold must be between 0 and 100, got {threshold}")
        self.threshold = threshold

    def predict(
        self,
        ensemble: SVMEnsemble,
        introns: Sequence[Intron]
    ) -> Sequence[Intron]:
        """
        Classify introns using trained ensemble.

        Args:
            ensemble: Trained ensemble of SVM models
            introns: Introns to classify (must have z-scores)

        Returns:
            Introns with updated classification scores

        Raises:
            ValueError: If introns lack z-scores
        """
        if not ensemble.models:
            raise ValueError("Ensemble has no models")

        # Prepare feature matrix
        X = self._prepare_features(introns)

        # Get predictions from each model
        probas = []
        for model in ensemble.models:
            # predict_proba returns [P(U2), P(U12)]
            # We want P(U12) which is column 1
            proba = model.model.predict_proba(X)[:, 1]
            probas.append(proba)

        probas = np.array(probas)  # Shape: (n_models, n_introns)

        # F1-weighted averaging (Port from: intronIC.py:5671-5676)
        f1_scores = np.array([m.f1_score for m in ensemble.models])
        weights = f1_scores / f1_scores.sum()

        # Weighted average across models
        avg_probas = np.dot(weights, probas)  # Shape: (n_introns,)

        # Convert to 0-100 scale (Port from: intronIC.py:5678)
        svm_scores = avg_probas * 100.0

        # Calculate decision function scores for relative_score
        # Use first model for consistency (all use same C parameter)
        decision_scores = ensemble.models[0].model.decision_function(X)

        # Update introns with classification results
        classified_introns = []
        for i, intron in enumerate(introns):
            svm_score = float(svm_scores[i])
            relative_score = float(decision_scores[i])
            type_id = 'u12' if svm_score >= self.threshold else 'u2'

            # Update scores (create new IntronScores with added fields)
            new_scores = replace(
                intron.scores,
                svm_score=svm_score,
                relative_score=relative_score
            )

            # Update metadata with type_id
            # IntronMetadata is mutable, but since Intron is frozen,
            # we need to create a new metadata object
            from core.intron import IntronMetadata
            if intron.metadata is None:
                new_metadata = IntronMetadata(type_id=type_id)
            else:
                new_metadata = replace(
                    intron.metadata,
                    type_id=type_id
                )

            # Update intron with new scores and metadata
            new_intron = replace(
                intron,
                scores=new_scores,
                metadata=new_metadata
            )

            classified_introns.append(new_intron)

        return classified_introns

    def _prepare_features(self, introns: Sequence[Intron]) -> np.ndarray:
        """
        Extract feature matrix from introns.

        Features: [five_z_score, bp_z_score, three_z_score]

        CRITICAL: Does NOT re-normalize z-scores.
        Uses z-scores computed from reference data only (prevents data leakage).

        Args:
            introns: Introns with z-scores

        Returns:
            Feature matrix (n_introns, 3)

        Raises:
            ValueError: If any intron lacks z-scores
        """
        features = []
        for intron in introns:
            if intron.scores is None:
                raise ValueError(f"Intron {intron.intron_id} has no scores")
            if (intron.scores.five_z_score is None or
                intron.scores.bp_z_score is None or
                intron.scores.three_z_score is None):
                raise ValueError(f"Intron {intron.intron_id} missing z-scores")

            features.append([
                intron.scores.five_z_score,
                intron.scores.bp_z_score,
                intron.scores.three_z_score
            ])

        return np.array(features)

    def predict_batch(
        self,
        ensemble: SVMEnsemble,
        introns: Sequence[Intron],
        batch_size: int = 10000
    ) -> Sequence[Intron]:
        """
        Classify introns in batches for memory efficiency.

        Useful for very large datasets where loading all features
        into memory at once would be problematic.

        Args:
            ensemble: Trained ensemble
            introns: Introns to classify
            batch_size: Number of introns per batch (default: 10000)

        Returns:
            Classified introns
        """
        classified = []

        for i in range(0, len(introns), batch_size):
            batch = introns[i:i + batch_size]
            classified_batch = self.predict(ensemble, batch)
            classified.extend(classified_batch)

        return classified
