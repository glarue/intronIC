"""
SVM ensemble prediction with F1-weighted averaging (CORRECTED ARCHITECTURE).

This module implements the classification algorithm from intronIC.py:5651-5900,
which applies trained ensemble models to classify introns as U2 or U12 type.

Key features:
- F1-weighted averaging of ensemble predictions for robust classification
- Decision boundary distance calculation for confidence estimation
- Type assignment based on probability threshold (default: 90%)
- Parallel processing support for large datasets

CORRECTED ARCHITECTURE (Expert guidance, 2025):
- Introns must have Z-SCORES (five_z_score, bp_z_score, three_z_score)
  populated by ScoreNormalizer BEFORE prediction
- Pipeline contains NO scaler (scaling done externally by ScoreNormalizer)
- Single scaling step prevents double-scaling issues
- Domain adaptation via ScoreNormalizer refitting (per-species)

Port from: intronIC.py:5651-5687 (average_svm_score_info), 5816-5900 (parallel_svm_score)
Redesign: SCALER_ARCHITECTURE_REVIEW.md
"""

import os
from dataclasses import replace
from multiprocessing import Pool, cpu_count
from typing import Sequence

import numpy as np

from intronIC.classification.trainer import SVMEnsemble
from intronIC.core.intron import Intron, IntronMetadata, IntronScores


def _predict_chunk_worker(
    ensemble: SVMEnsemble, introns: Sequence[Intron], threshold: float
) -> Sequence[Intron]:
    """
    Worker function for parallel classification (NEW ARCHITECTURE).

    This function is defined at module level so it can be pickled for
    multiprocessing. It processes a chunk of introns using the same logic
    as SVMPredictor._predict_chunk().

    Args:
        ensemble: Trained SVM ensemble
        introns: Chunk of introns to classify
        threshold: U12 probability threshold

    Returns:
        Classified introns with scores and type_id

    Redesign: SCALER_ARCHITECTURE_REVIEW.md (Expert feedback, 2025)
    """
    if not introns:
        return []

    # CORRECTED: Extract z-scores from introns (already scaled by ScoreNormalizer)
    # See: Expert workflow doc, SCALER_ARCHITECTURE_REVIEW.md
    #
    # Data flow:
    # 1. ScoreNormalizer (in main.py) fits RobustScaler on reference LLRs
    # 2. ScoreNormalizer transforms experimental LLRs → z-scores (five_z_score, bp_z_score, three_z_score)
    # 3. Predictor extracts z-scores and passes to model
    # 4. Pipeline has NO scaler (would cause double-scaling)
    features = []
    for intron in introns:
        if intron.scores is None:
            raise ValueError(f"Intron {intron.intron_id} has no scores")
        if (
            intron.scores.five_z_score is None
            or intron.scores.bp_z_score is None
            or intron.scores.three_z_score is None
        ):
            raise ValueError(f"Intron {intron.intron_id} missing z-scores")

        features.append(
            [
                intron.scores.five_z_score,
                intron.scores.bp_z_score,
                intron.scores.three_z_score,
            ]
        )

    X_z = np.array(features)

    # Get predictions from each model (pass z-scores - NO scaling in pipeline)
    probas = []
    for model in ensemble.models:
        proba = model.model.predict_proba(X_z)[:, 1]
        probas.append(proba)

    probas = np.array(probas)

    # F1-weighted averaging (if f1_scores available), else equal weights
    # In nested CV, models don't have f1_scores, so use equal weights
    try:
        f1_scores = np.array([m.f1_score for m in ensemble.models])
        weights = f1_scores / f1_scores.sum()
    except AttributeError:
        # Models don't have f1_scores (e.g., from nested CV)
        # Use equal weights
        weights = np.ones(len(ensemble.models)) / len(ensemble.models)

    avg_probas = np.dot(weights, probas)

    # Convert to 0-100 scale
    svm_scores = avg_probas * 100.0

    # Calculate relative scores as distance from threshold
    # This ensures U12s (score >= threshold) have relative_score >= 0
    # and U2s (score < threshold) have relative_score < 0
    # Matches original intronIC: relative_score = svm_score - threshold
    relative_scores = svm_scores - threshold

    # Also calculate log-odds for decision_distance (alternative confidence metric)
    # log(p / (1-p)) gives negative for U2, positive for U12, zero at p=0.5
    epsilon = 1e-10  # Avoid log(0)
    clipped_probas = np.clip(avg_probas, epsilon, 1 - epsilon)
    log_odds = np.log(clipped_probas / (1 - clipped_probas))

    # Get include_max parameter from first model (all models use same pipeline parameters)
    include_max = ensemble.models[0].parameters.include_max

    # Update introns with classification results
    classified_introns = []
    for i, intron in enumerate(
        introns
    ):  # introns already have z-scores from ScoreNormalizer
        svm_score = float(svm_scores[i])
        relative_score = float(relative_scores[i])
        decision_distance = float(log_odds[i])

        # type_id based on decision_distance (log-odds) where > 0 indicates U12
        # This matches original intronIC behavior where threshold only affects
        # which introns are considered "high confidence" U12s for reporting/filtering
        # decision_distance > 0 is equivalent to probability > 50% (the raw classifier decision)
        type_id = "u12" if decision_distance > 0 else "u2"

        # Compute BothEndsStrong augmented features for output
        # These are computed from z-scores (what the model sees after scaling)
        five_z = intron.scores.five_z_score
        bp_z = intron.scores.bp_z_score
        three_z = intron.scores.three_z_score

        # min(a, b) = 0.5 * ((a + b) - |a - b|)
        # max(a, b) = 0.5 * ((a + b) + |a - b|)
        sum_5_bp = five_z + bp_z
        absdiff_5_bp = abs(five_z - bp_z)
        min_5_bp = 0.5 * (sum_5_bp - absdiff_5_bp)

        sum_5_3 = five_z + three_z
        absdiff_5_3 = abs(five_z - three_z)
        min_5_3 = 0.5 * (sum_5_3 - absdiff_5_3)

        max_5_bp = None
        max_5_3 = None
        if include_max:
            max_5_bp = 0.5 * (sum_5_bp + absdiff_5_bp)
            max_5_3 = 0.5 * (sum_5_3 + absdiff_5_3)

        # Update scores
        new_scores = replace(
            intron.scores,
            svm_score=svm_score,
            relative_score=relative_score,
            decision_distance=decision_distance,
            min_5_bp=min_5_bp,
            min_5_3=min_5_3,
            max_5_bp=max_5_bp,
            max_5_3=max_5_3,
        )

        # Update metadata with type_id
        if intron.metadata is None:
            new_metadata = IntronMetadata(type_id=type_id)
        else:
            new_metadata = replace(intron.metadata, type_id=type_id)

        # Update intron
        new_intron = replace(intron, scores=new_scores, metadata=new_metadata)

        classified_introns.append(new_intron)

    return classified_introns


class SVMPredictor:
    """
    Apply trained SVM ensemble to classify introns (NEW ARCHITECTURE).

    Uses F1-weighted averaging across models for robust predictions.

    NEW: Operates on RAW PWM scores; pipeline handles scaling internally.
    This prevents double-scaling and enables cross-species deployment.

    Port from: intronIC.py:5651-5900
    Redesign: SCALER_ARCHITECTURE_REVIEW.md (Expert feedback, 2025)
    """

    def __init__(self, threshold: float = 90.0, n_jobs: int = 1):
        """
        Initialize predictor.

        Args:
            threshold: U12 probability threshold (0-100, default: 90)
            n_jobs: Number of parallel processes for classification (default: 1)
        """
        if not 0 <= threshold <= 100:
            raise ValueError(f"Threshold must be between 0 and 100, got {threshold}")
        self.threshold = threshold
        self.n_jobs = n_jobs

    def predict(
        self, ensemble: SVMEnsemble, introns: Sequence[Intron]
    ) -> Sequence[Intron]:
        """
        Classify introns using trained ensemble (NEW ARCHITECTURE).

        Uses parallel processing if n_jobs > 1. For efficiency, introns are
        split into chunks (one per worker) and each chunk is processed using
        vectorized sklearn operations.

        NEW: Scaling is handled externally by ScoreNormalizer (RobustScaler).
        Introns must have z-scores pre-computed, NOT raw PWM scores.
        This architecture prevents double-scaling issues.

        Args:
            ensemble: Trained ensemble of SVM models
            introns: Introns to classify (must have raw PWM scores)

        Returns:
            Introns with updated classification scores

        Raises:
            ValueError: If introns lack raw scores

        Redesign: SCALER_ARCHITECTURE_REVIEW.md (Expert feedback, 2025)
        """
        if not ensemble.models:
            raise ValueError("Ensemble has no models")

        # Sequential processing for n_jobs=1 or small datasets
        if self.n_jobs == 1 or len(introns) < 100:
            return self._predict_chunk(ensemble, introns)

        # Prevent thread oversubscription when using parallelization
        # Pool workers may use multi-threaded BLAS (OpenBLAS, MKL, etc.), causing
        # n_jobs × BLAS_threads competing threads and severe slowdown.
        # Solution: Force single-threaded BLAS in each worker.
        os.environ["OPENBLAS_NUM_THREADS"] = "1"
        os.environ["MKL_NUM_THREADS"] = "1"
        os.environ["OMP_NUM_THREADS"] = "1"
        os.environ["NUMEXPR_NUM_THREADS"] = "1"

        # Parallel processing: split introns into chunks
        # Handle n_jobs=-1 (use all cores) like scikit-learn
        n_jobs_actual = cpu_count() if self.n_jobs == -1 else self.n_jobs
        n_workers = min(n_jobs_actual, cpu_count(), len(introns))
        chunk_size = max(1, len(introns) // n_workers)

        # Create chunks
        chunks = []
        for i in range(0, len(introns), chunk_size):
            chunk = introns[i : i + chunk_size]
            chunks.append(chunk)

        # Process chunks in parallel
        # Use Pool.starmap to pass ensemble and chunk to worker function
        with Pool(processes=n_workers) as pool:
            try:
                # Each worker processes its chunk independently
                chunk_results = pool.starmap(
                    _predict_chunk_worker,
                    [(ensemble, chunk, self.threshold) for chunk in chunks],
                )
            except KeyboardInterrupt:
                pool.terminate()
                raise
            finally:
                pool.close()
                pool.join()

        # Flatten results from all chunks
        classified_introns = []
        for chunk_result in chunk_results:
            classified_introns.extend(chunk_result)

        return classified_introns

    def _predict_chunk(
        self, ensemble: SVMEnsemble, introns: Sequence[Intron]
    ) -> Sequence[Intron]:
        """
        Classify a chunk of introns (internal helper for predict()).

        This method performs the actual classification using vectorized
        sklearn operations for efficiency.

        Args:
            ensemble: Trained ensemble of SVM models
            introns: Introns to classify

        Returns:
            Classified introns with scores and type_id
        """
        if not introns:
            return []

        # CORRECTED: Extract z-scores from introns (already scaled by ScoreNormalizer)
        # See: Expert workflow doc, SCALER_ARCHITECTURE_REVIEW.md
        features = []
        for intron in introns:
            if intron.scores is None:
                raise ValueError(f"Intron {intron.intron_id} has no scores")
            if (
                intron.scores.five_z_score is None
                or intron.scores.bp_z_score is None
                or intron.scores.three_z_score is None
            ):
                raise ValueError(f"Intron {intron.intron_id} missing z-scores")

            features.append(
                [
                    intron.scores.five_z_score,
                    intron.scores.bp_z_score,
                    intron.scores.three_z_score,
                ]
            )

        X_z = np.array(features)

        # Get predictions from each model (pass z-scores - NO scaling in pipeline)
        probas = []
        for model in ensemble.models:
            # predict_proba returns [P(U2), P(U12)]
            # We want P(U12) which is column 1
            proba = model.model.predict_proba(X_z)[:, 1]
            probas.append(proba)

        probas = np.array(probas)  # Shape: (n_models, n_introns)

        # F1-weighted averaging (Port from: intronIC.py:5671-5676)
        # Handle models that may not have f1_score (e.g., loaded from old format)
        f1_scores = np.array([getattr(m, "f1_score", 1.0) for m in ensemble.models])
        weights = f1_scores / f1_scores.sum()

        # Weighted average across models
        avg_probas = np.dot(weights, probas)  # Shape: (n_introns,)

        # Convert to 0-100 scale (Port from: intronIC.py:5678)
        svm_scores = avg_probas * 100.0

        # Calculate relative scores as distance from threshold
        # This ensures U12s (score >= threshold) have relative_score >= 0
        # and U2s (score < threshold) have relative_score < 0
        # Matches original intronIC: relative_score = svm_score - threshold
        relative_scores = svm_scores - self.threshold

        # Also calculate log-odds for decision_distance (alternative confidence metric)
        # log(p / (1-p)) gives negative for U2, positive for U12, zero at p=0.5
        epsilon = 1e-10  # Avoid log(0)
        clipped_probas = np.clip(avg_probas, epsilon, 1 - epsilon)
        log_odds = np.log(clipped_probas / (1 - clipped_probas))

        # Get include_max parameter from first model (all models use same pipeline parameters)
        include_max = ensemble.models[0].parameters.include_max

        # Update introns with classification results
        classified_introns = []
        for i, intron in enumerate(
            introns
        ):  # introns already have z-scores from ScoreNormalizer
            svm_score = float(svm_scores[i])
            relative_score = float(relative_scores[i])
            decision_distance = float(log_odds[i])

            # type_id based on decision_distance (log-odds) where > 0 indicates U12
            # This matches original intronIC behavior where threshold only affects
            # which introns are considered "high confidence" U12s for reporting/filtering
            # decision_distance > 0 is equivalent to probability > 50% (the raw classifier decision)
            type_id = "u12" if decision_distance > 0 else "u2"

            # Compute BothEndsStrong augmented features for output
            # These are computed from z-scores (what the model sees after scaling)
            five_z = intron.scores.five_z_score
            bp_z = intron.scores.bp_z_score
            three_z = intron.scores.three_z_score

            # min(a, b) = 0.5 * ((a + b) - |a - b|)
            # max(a, b) = 0.5 * ((a + b) + |a - b|)
            sum_5_bp = five_z + bp_z
            absdiff_5_bp = abs(five_z - bp_z)
            min_5_bp = 0.5 * (sum_5_bp - absdiff_5_bp)

            sum_5_3 = five_z + three_z
            absdiff_5_3 = abs(five_z - three_z)
            min_5_3 = 0.5 * (sum_5_3 - absdiff_5_3)

            max_5_bp = None
            max_5_3 = None
            if include_max:
                max_5_bp = 0.5 * (sum_5_bp + absdiff_5_bp)
                max_5_3 = 0.5 * (sum_5_3 + absdiff_5_3)

            # Update scores (create new IntronScores with added fields)
            new_scores = replace(
                intron.scores,
                svm_score=svm_score,
                relative_score=relative_score,
                decision_distance=decision_distance,
                min_5_bp=min_5_bp,
                min_5_3=min_5_3,
                max_5_bp=max_5_bp,
                max_5_3=max_5_3,
            )

            # Update metadata with type_id
            # IntronMetadata is mutable, but since Intron is frozen,
            # we need to create a new metadata object
            if intron.metadata is None:
                new_metadata = IntronMetadata(type_id=type_id)
            else:
                new_metadata = replace(intron.metadata, type_id=type_id)

            # Update intron with new scores and metadata
            new_intron = replace(intron, scores=new_scores, metadata=new_metadata)

            classified_introns.append(new_intron)

        return classified_introns

    def _prepare_features(self, introns: Sequence[Intron]) -> np.ndarray:
        """
        Extract feature matrix from introns (NEW ARCHITECTURE).

        Features: [five_raw_score, bp_raw_score, three_raw_score]

        CRITICAL: Extracts RAW LLR scores, NOT z-scores.
        Scaling is handled externally by ScoreNormalizer (RobustScaler).
        The pipeline receives pre-scaled z-scores after normalization.
        This prevents double-scaling issues and allows cross-species deployment.

        Args:
            introns: Introns with raw PWM scores

        Returns:
            Feature matrix (n_introns, 3) of raw LLR scores

        Raises:
            ValueError: If any intron lacks raw scores

        Redesign: SCALER_ARCHITECTURE_REVIEW.md (Expert feedback, 2025)
        """
        features = []
        for intron in introns:
            if intron.scores is None:
                raise ValueError(f"Intron {intron.intron_id} has no scores")
            if (
                intron.scores.five_raw_score is None
                or intron.scores.bp_raw_score is None
                or intron.scores.three_raw_score is None
            ):
                raise ValueError(f"Intron {intron.intron_id} missing raw scores")

            features.append(
                [
                    intron.scores.five_raw_score,
                    intron.scores.bp_raw_score,
                    intron.scores.three_raw_score,
                ]
            )

        return np.array(features)

    def predict_batch(
        self, ensemble: SVMEnsemble, introns: Sequence[Intron], batch_size: int = 10000
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
            batch = introns[i : i + batch_size]
            classified_batch = self.predict(ensemble, batch)
            classified.extend(classified_batch)

        return classified


def classify_introns_streaming(
    introns: "Iterator[Intron]",
    ensemble: SVMEnsemble,
    threshold: float = 90.0,
) -> "Iterator[Intron]":
    """
    Classify introns one at a time for memory-efficient streaming.

    This function yields classified introns one at a time, allowing
    immediate processing (writing to files, freeing memory) without
    accumulating all introns in memory.

    Used in true streaming mode with pre-trained models for large
    genomes where peak memory is a concern.

    Args:
        introns: Iterator of introns with z-scores populated
                (from score_and_normalize_introns)
        ensemble: Pre-trained SVM ensemble from model bundle
        threshold: U12 probability threshold (default: 90.0)

    Yields:
        Introns with classification scores (svm_score, relative_score,
        decision_distance, type_id) populated

    Example:
        >>> # Load pre-trained model
        >>> ensemble = model_bundle["ensemble"]
        >>> normalizer = model_bundle["normalizer"]
        >>> scaler = normalizer.get_frozen_scaler()
        >>>
        >>> # Stream through introns: score → normalize → classify → write
        >>> scored = score_and_normalize_introns(introns, scorer, scaler)
        >>> for intron in classify_introns_streaming(scored, ensemble):
        ...     write_intron_to_outputs(intron)
        ...     # Memory freed after each intron

    Note:
        For efficiency with batches, use classify_introns_batch instead.
        This function is optimized for memory, not speed.
    """
    from typing import Iterator

    # Get include_max parameter from first model
    include_max = ensemble.models[0].parameters.include_max

    # Pre-compute weights (same for all introns)
    try:
        f1_scores = np.array([m.f1_score for m in ensemble.models])
        weights = f1_scores / f1_scores.sum()
    except AttributeError:
        weights = np.ones(len(ensemble.models)) / len(ensemble.models)

    for intron in introns:
        # Validate intron has z-scores
        if intron.scores is None:
            raise ValueError(f"Intron {intron.intron_id} has no scores")
        if (
            intron.scores.five_z_score is None
            or intron.scores.bp_z_score is None
            or intron.scores.three_z_score is None
        ):
            raise ValueError(f"Intron {intron.intron_id} missing z-scores")

        # Build feature vector
        X_z = np.array(
            [
                [
                    intron.scores.five_z_score,
                    intron.scores.bp_z_score,
                    intron.scores.three_z_score,
                ]
            ]
        )

        # Get predictions from each model
        probas = []
        for model in ensemble.models:
            proba = model.model.predict_proba(X_z)[:, 1]
            probas.append(proba[0])

        probas = np.array(probas)

        # Weighted average
        avg_proba = float(np.dot(weights, probas))

        # Convert to 0-100 scale
        svm_score = avg_proba * 100.0

        # Relative score as distance from threshold
        relative_score = svm_score - threshold

        # Log-odds for decision_distance
        epsilon = 1e-10
        clipped_proba = max(epsilon, min(1 - epsilon, avg_proba))
        decision_distance = float(np.log(clipped_proba / (1 - clipped_proba)))

        # Type assignment based on decision_distance
        type_id = "u12" if decision_distance > 0 else "u2"

        # Compute augmented features
        five_z = intron.scores.five_z_score
        bp_z = intron.scores.bp_z_score
        three_z = intron.scores.three_z_score

        sum_5_bp = five_z + bp_z
        absdiff_5_bp = abs(five_z - bp_z)
        min_5_bp = 0.5 * (sum_5_bp - absdiff_5_bp)

        sum_5_3 = five_z + three_z
        absdiff_5_3 = abs(five_z - three_z)
        min_5_3 = 0.5 * (sum_5_3 - absdiff_5_3)

        max_5_bp = None
        max_5_3 = None
        if include_max:
            max_5_bp = 0.5 * (sum_5_bp + absdiff_5_bp)
            max_5_3 = 0.5 * (sum_5_3 + absdiff_5_3)

        # Update scores
        new_scores = replace(
            intron.scores,
            svm_score=svm_score,
            relative_score=relative_score,
            decision_distance=decision_distance,
            min_5_bp=min_5_bp,
            min_5_3=min_5_3,
            max_5_bp=max_5_bp,
            max_5_3=max_5_3,
        )

        # Update metadata with type_id
        if intron.metadata is None:
            new_metadata = IntronMetadata(type_id=type_id)
        else:
            new_metadata = replace(intron.metadata, type_id=type_id)

        yield replace(intron, scores=new_scores, metadata=new_metadata)


def classify_introns_batch(
    introns: list["Intron"],
    ensemble: SVMEnsemble,
    threshold: float = 90.0,
) -> list["Intron"]:
    """
    Classify a batch of introns efficiently (vectorized).

    This is more efficient than classify_introns_streaming when
    processing a batch (e.g., all introns from one chromosome).
    Uses vectorized numpy operations for better performance.

    Args:
        introns: List of introns with z-scores populated
        ensemble: Pre-trained SVM ensemble from model bundle
        threshold: U12 probability threshold (default: 90.0)

    Returns:
        List of introns with classification scores populated

    Example:
        >>> # Process one chromosome at a time
        >>> chromosome_introns = score_and_normalize_batch(introns, scorer, scaler)
        >>> classified = classify_introns_batch(chromosome_introns, ensemble)
        >>> write_batch_to_outputs(classified)
    """
    if not introns:
        return []

    # Use the existing _predict_chunk_worker which is already vectorized
    return list(_predict_chunk_worker(ensemble, introns, threshold))
