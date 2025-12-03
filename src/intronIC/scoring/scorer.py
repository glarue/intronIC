"""
IntronScorer: Main scoring pipeline orchestrator.

This module orchestrates the full scoring process:
1. PWM scoring for 5' splice site (with U12 and U2 matrices)
2. Branch point detection and scoring
3. PWM scoring for 3' splice site
4. Log-ratio calculation (log2(U12/U2) for each region)
5. Population of intron score objects

Port from: intronIC.py:3115-3400 (get_raw_scores, assign_raw_score, multi_matrix_score)

Design:
- Single class orchestrates all scoring
- Clear separation of concerns (5', BP, 3' scoring)
- Returns new introns with scores populated (functional style)
- No multiprocessing (keep simple for now)
"""

import math
from dataclasses import replace
from typing import Dict, Iterator, Optional, Set, Tuple

from intronIC.core.intron import Intron, IntronScores, IntronSequences
from intronIC.scoring.branch_point import BranchPointMatch, BranchPointScorer
from intronIC.scoring.pwm import PWMSet


class IntronScorer:
    """
    Main scoring pipeline: PWM scoring + branch point detection + log ratios.

    This class orchestrates the full scoring process while maintaining
    clear separation between U12 and U2 scoring.

    Port from: intronIC.py:3040-3112 (assign_raw_score), 3115-3144 (get_raw_scores)

    Attributes:
        pwm_sets: Dictionary of PWMSets for each region (five, bp, three)
        five_coords: Tuple of (start, stop) for 5' scoring region
        bp_coords: Tuple of (start, stop) for branch point search region
        three_coords: Tuple of (start, stop) for 3' scoring region
        ignore_nc_dnts: Whether to ignore non-canonical dinucleotides in scoring
        bp_scorer: BranchPointScorer instance
    """

    def __init__(
        self,
        pwm_sets: Dict[str, PWMSet],
        five_coords: Tuple[int, int] = (-3, 9),
        bp_coords: Tuple[int, int] = (-55, -5),
        three_coords: Tuple[int, int] = (-6, 4),
        ignore_nc_dnts: bool = True,
    ):
        """
        Initialize scoring pipeline.

        Args:
            pwm_sets: Dictionary with keys 'five', 'bp', 'three'
                     Values are PWMSet objects with u2/u12 canonical/noncanonical PWMs
            five_coords: Coordinates for 5' scoring region relative to intron start
                        Default: (-3, 9) matches original intronIC
            bp_coords: Coordinates for branch point search relative to intron 3' end
                      Default: (-55, -5) matches original intronIC
            three_coords: Coordinates for 3' scoring region relative to intron stop
                         Default: (-6, 4) matches original intronIC
            ignore_nc_dnts: If True, ignore dinucleotide positions when scoring
                           non-canonical introns (default: True)
        """
        self.pwm_sets = pwm_sets
        self.five_coords = five_coords
        self.bp_coords = bp_coords
        self.three_coords = three_coords
        self.ignore_nc_dnts = ignore_nc_dnts

        # Branch point scorer will be created with appropriate PWMs per intron
        # No longer storing a single scorer instance

    def score_introns(self, introns: Iterator[Intron]) -> Iterator[Intron]:
        """
        Score multiple introns (generator).

        Port from: intronIC.py:3115-3144 (get_raw_scores)

        Args:
            introns: Iterable of introns to score

        Yields:
            Introns with raw scores populated

        Example:
            >>> scorer = IntronScorer(pwm_sets)
            >>> scored = list(scorer.score_introns(introns))
        """
        for intron in introns:
            yield self.score_intron(intron)

    def score_intron(self, intron: Intron) -> Intron:
        """
        Score a single intron through the full pipeline.

        Port from: intronIC.py:3040-3112 (assign_raw_score)

        This is the main scoring method that:
        1. Scores 5' splice site with U12 and U2 PWMs
        2. Finds and scores branch point
        3. Scores 3' splice site with U12 and U2 PWMs
        4. Calculates log ratios for each region
        5. Returns intron with scores populated

        Supports two modes:
        - Default mode: Uses full sequences (intron.sequences)
        - Streaming mode: Uses pre-extracted motifs (intron.motifs)

        Args:
            intron: Intron to score

        Returns:
            New intron with scores populated

        Raises:
            ValueError: If intron lacks required sequences or motifs
        """
        # Validate intron has sequences OR motifs (streaming mode)
        has_sequences = (
            intron.sequences is not None and intron.sequences.seq is not None
        )
        has_motifs = intron.motifs is not None

        if not has_sequences and not has_motifs:
            raise ValueError(
                f"Intron {intron.intron_id} has no sequence data. "
                "Cannot score without sequences or motifs."
            )

        # Determine which positions to ignore based on canonical status
        ignore_five, ignore_three = self._get_ignore_positions(intron)

        # Score 5' splice site
        # Port from: intronIC.py:3055-3072 (calls multi_matrix_score)
        five_u12_score, five_u2_score = self._score_five_site(
            intron, ignore_positions=ignore_five
        )

        # Score branch point
        # Port from: intronIC.py:3078-3084 (extracts BP info)
        bp_u12_match, bp_u2_score = self._score_branch_point(intron)

        # Score 3' splice site
        three_u12_score, three_u2_score = self._score_three_site(
            intron, ignore_positions=ignore_three
        )

        # Calculate log ratios
        # Port from: intronIC.py:3086-3103 (log_ratio calls)
        five_raw_score = self._calculate_log_ratio(five_u12_score, five_u2_score)

        # Handle None match from short introns - use pseudocount for both U12 and U2
        # Port from: intronIC.py:2944
        if bp_u12_match is None:
            # When match is None, bp_u2_score already contains pseudocount
            # Use same pseudocount for U12 score for consistency
            bp_raw_score = self._calculate_log_ratio(bp_u2_score, bp_u2_score)
        else:
            bp_raw_score = self._calculate_log_ratio(bp_u12_match.score, bp_u2_score)

        three_raw_score = self._calculate_log_ratio(three_u12_score, three_u2_score)

        # Update intron with scores
        # Port from: intronIC.py:3099-3111
        updated_scores = replace(
            intron.scores if intron.scores else IntronScores(),
            five_raw_score=five_raw_score,
            bp_raw_score=bp_raw_score,
            three_raw_score=three_raw_score,
        )

        # Update sequences with BP information if match was found
        # Port from: intronIC.py:3080-3084, 3105-3106
        updated_sequences = intron.sequences
        if bp_u12_match is not None and intron.sequences is not None:
            updated_sequences = replace(
                intron.sequences,
                bp_seq=bp_u12_match.sequence,
                bp_seq_u2=bp_u12_match.sequence_u2,
                bp_relative_coords=(
                    bp_u12_match.start_in_region,
                    bp_u12_match.stop_in_region,
                ),
                bp_region_seq=bp_u12_match.search_region,  # Store the actual searched sequence
            )

        return replace(intron, scores=updated_scores, sequences=updated_sequences)

    def _get_ignore_positions(
        self, intron: Intron
    ) -> Tuple[Optional[Set[int]], Optional[Set[int]]]:
        """
        Determine which positions to ignore for non-canonical introns.

        Port from: intronIC.py:2885-2898 (ignore_nc_dnts logic)

        For non-canonical introns (GC-AG, etc.), we ignore the dinucleotide
        positions when scoring to focus on the surrounding motif.

        Args:
            intron: Intron to check

        Returns:
            Tuple of (ignore_five, ignore_three) sets of positions to ignore
        """
        if not self.ignore_nc_dnts:
            return None, None

        # Check if intron is non-canonical
        is_noncanonical = intron.metadata and intron.metadata.noncanonical is True

        if not is_noncanonical:
            return None, None

        # For non-canonical, ignore the dinucleotide positions
        # 5' dinucleotide is at positions 0 and 1 (relative to intron start)
        # 3' dinucleotide is at positions -2 and -1 (relative to intron end)
        ignore_five = {0, 1}
        ignore_three = {-2, -1}

        return ignore_five, ignore_three

    def _get_terminal_dinucleotides(self, intron: Intron) -> str:
        """
        Get terminal dinucleotides for PWM selection.

        Supports both default mode (from sequences) and streaming mode (from motifs).

        Args:
            intron: Intron to get dinucleotides from

        Returns:
            Lowercase dinucleotide string (e.g., 'gtag', 'gcag', 'atac')
        """
        # Streaming mode: use motifs.terminal_dnts (format: 'GT-AG')
        if intron.motifs and intron.motifs.terminal_dnts:
            return intron.motifs.terminal_dnts.replace("-", "").lower()

        # Default mode: get from sequences
        if (
            intron.sequences
            and intron.sequences.five_prime_dnt
            and intron.sequences.three_prime_dnt
        ):
            return (
                intron.sequences.five_prime_dnt + intron.sequences.three_prime_dnt
            ).lower()

        # Fallback to GT-AG
        return "gtag"

    def _score_five_site(
        self, intron: Intron, ignore_positions: Optional[Set[int]] = None
    ) -> Tuple[float, float]:
        """
        Score 5' splice site with U12 and U2 PWMs.

        Selects appropriate matrices based on intron dinucleotides.
        Tries all available PWM versions and selects the highest-scoring ones.

        Port from: intronIC.py:2862-2944 (multi_matrix_score for 'five' region)

        Args:
            intron: Intron to score
            ignore_positions: Positions to ignore in scoring

        Returns:
            Tuple of (best_u12_score, best_u2_score)
        """
        # Extract 5' region
        five_region = self._extract_five_region(intron)

        # Get dinucleotides for PWM selection (handles both default and streaming modes)
        dnts = self._get_terminal_dinucleotides(intron)

        # Get all available PWM versions
        # Port from: intronIC.py:2915-2942 (multi_matrix_score tries all versions)
        u12_pwms = self.pwm_sets["five"].get_all_versions("u12", dnts)
        u2_pwms = self.pwm_sets["five"].get_all_versions("u2", dnts)

        # Pass the starting position of the 5' region (first coordinate)
        seq_start_pos = self.five_coords[0]  # e.g., -3

        # Try all U12 PWM versions and keep the best score
        best_u12_score = 0.0
        for u12_pwm in u12_pwms:
            score = u12_pwm.score_sequence(
                five_region,
                seq_start_position=seq_start_pos,
                ignore_positions=ignore_positions,
            )
            if score > best_u12_score:
                best_u12_score = score

        # Try all U2 PWM versions and keep the best score
        best_u2_score = 0.0
        for u2_pwm in u2_pwms:
            score = u2_pwm.score_sequence(
                five_region,
                seq_start_position=seq_start_pos,
                ignore_positions=ignore_positions,
            )
            if score > best_u2_score:
                best_u2_score = score

        return best_u12_score, best_u2_score

    def _score_three_site(
        self, intron: Intron, ignore_positions: Optional[Set[int]] = None
    ) -> Tuple[float, float]:
        """
        Score 3' splice site with U12 and U2 PWMs.

        Selects appropriate matrices based on intron dinucleotides.
        Tries all available PWM versions and selects the highest-scoring ones.

        Port from: intronIC.py:2862-2944 (multi_matrix_score for 'three' region)

        Args:
            intron: Intron to score
            ignore_positions: Positions to ignore in scoring

        Returns:
            Tuple of (best_u12_score, best_u2_score)
        """
        # Extract 3' region
        three_region = self._extract_three_region(intron)

        # Get dinucleotides for PWM selection (handles both default and streaming modes)
        dnts = self._get_terminal_dinucleotides(intron)

        # Get all available PWM versions
        # Port from: intronIC.py:2915-2942 (multi_matrix_score tries all versions)
        u12_pwms = self.pwm_sets["three"].get_all_versions("u12", dnts)
        u2_pwms = self.pwm_sets["three"].get_all_versions("u2", dnts)

        # Pass the starting position of the 3' region (first coordinate)
        seq_start_pos = self.three_coords[0]  # e.g., -6

        # Try all U12 PWM versions and keep the best score
        best_u12_score = 0.0
        for u12_pwm in u12_pwms:
            score = u12_pwm.score_sequence(
                three_region,
                seq_start_position=seq_start_pos,
                ignore_positions=ignore_positions,
            )
            if score > best_u12_score:
                best_u12_score = score

        # Try all U2 PWM versions and keep the best score
        best_u2_score = 0.0
        for u2_pwm in u2_pwms:
            score = u2_pwm.score_sequence(
                three_region,
                seq_start_position=seq_start_pos,
                ignore_positions=ignore_positions,
            )
            if score > best_u2_score:
                best_u2_score = score

        return best_u12_score, best_u2_score

    def _score_branch_point(
        self, intron: Intron
    ) -> Tuple[BranchPointMatch | None, float]:
        """
        Find and score branch point with both U12 and U2 PWMs.

        Selects appropriate matrices based on intron dinucleotides.
        Tries all available PWM versions (e.g., vA9, vA10) and selects
        the highest-scoring U12 version.

        Port from: intronIC.py:2920-2944, 3078-3084

        The U12 PWM is used to find the best position, then both U12 and U2
        PWMs score that sequence for the log ratio.

        For introns where the search window is too small (e.g., short introns),
        returns None for the match and uses a pseudocount for scoring.

        Args:
            intron: Intron to score

        Returns:
            Tuple of (u12_match, u2_score). u12_match will be None if the
            search window is too small for the PWM.
        """
        # Get dinucleotides for PWM selection (handles both default and streaming modes)
        dnts = self._get_terminal_dinucleotides(intron)

        # Get all available PWM versions for U12 and U2
        # Port from: intronIC.py:2915-2942 (multi_matrix_score tries all versions)
        u12_pwms = self.pwm_sets["bp"].get_all_versions("u12", dnts)
        u2_pwms = self.pwm_sets["bp"].get_all_versions("u2", dnts)

        # Try all U12 PWM versions and select the highest-scoring one
        # Port from: intronIC.py:2934-2942 (keeps best score per matrix_category)
        best_match = None
        best_u12_score = 0.0

        for u12_pwm in u12_pwms:
            # U2 typically has only one version, but use all just in case
            for u2_pwm in u2_pwms:
                # Create a branch point scorer with both U12 and U2 PWMs
                bp_scorer = BranchPointScorer(u12_pwm=u12_pwm, u2_pwm=u2_pwm)

                # Find best matches for both U12 and U2 PWMs
                match = bp_scorer.find_best_match(
                    intron,
                    search_window=self.bp_coords,
                    five_coords=self.five_coords,
                    three_coords=self.three_coords,
                )

                # If match is None (window too small), use pseudocount
                if match is None:
                    if best_match is None:
                        # Only set pseudocount if we haven't found any match yet
                        pseudocount_score = u2_pwm.pseudocount * u2_pwm.length
                        return None, pseudocount_score
                    continue

                # Compare raw product scores (higher is better)
                # Port from: intronIC.py:2937 (if target.get('score', 0) > s)
                if match.score > best_u12_score:
                    best_u12_score = match.score
                    best_match = match

        # If no match found (all windows too small), use pseudocount
        if best_match is None:
            # Use the first U2 PWM for pseudocount calculation
            pseudocount_score = u2_pwms[0].pseudocount * u2_pwms[0].length
            return None, pseudocount_score

        # Return the best match found across all versions
        return best_match, best_match.score_u2

    def _extract_five_region(self, intron: Intron) -> str:
        """
        Extract 5' splice site region for scoring.

        Supports two modes:
        - Streaming mode: Uses pre-extracted motifs (intron.motifs.five_region)
        - Default mode: Extracts from full sequences (intron.sequences)

        The 5' coordinates are relative to the intron start. Negative values
        require upstream flanking sequence.

        Args:
            intron: Intron object

        Returns:
            Sequence for 5' scoring region

        Example:
            For coords (-3, 5) on intron starting with GTAAGT...:
            - Positions: -3,-2,-1,0,1,2,3,4 (8 positions, stop is exclusive)
            - Need 3bp upstream flank + 5bp intron sequence
            - Returns: upstream[-3:] + intron.seq[:5]
        """
        # Streaming mode: use pre-extracted motif
        if intron.motifs is not None and intron.motifs.five_region:
            return intron.motifs.five_region

        # Default mode: extract from full sequences
        if intron.sequences is None or intron.sequences.seq is None:
            raise ValueError(f"Intron {intron.intron_id} has no sequence data")

        start, stop = self.five_coords

        # Calculate how much we need from upstream flank vs intron
        if start < 0:
            # Need upstream flank
            upstream = intron.sequences.upstream_flank or ""
            upstream_needed = abs(start)
            upstream_part = (
                upstream[-upstream_needed:]
                if len(upstream) >= upstream_needed
                else "N" * upstream_needed
            )

            # Intron part (from position 0 to stop)
            if stop <= 0:
                # All from upstream
                intron_part = ""
            else:
                intron_part = intron.sequences.seq[:stop]

            return upstream_part + intron_part
        else:
            # All from intron sequence
            return intron.sequences.seq[start:stop]

    def _extract_three_region(self, intron: Intron) -> str:
        """
        Extract 3' splice site region for scoring.

        Supports two modes:
        - Streaming mode: Uses pre-extracted motifs (intron.motifs.three_region)
        - Default mode: Extracts from full sequences (intron.sequences)

        The 3' coordinates are relative to the intron end. Negative values
        are positions before the end, positive values require downstream flank.

        Args:
            intron: Intron object

        Returns:
            Sequence for 3' scoring region

        Example:
            For coords (-6, 2) on intron ending with ...TTTCAG:
            - Positions: -6,-5,-4,-3,-2,-1,0,1 (8 positions, stop is exclusive)
            - Need 6bp before end + 2bp downstream
            - Returns: intron.seq[-6:] + downstream[:2]
        """
        # Streaming mode: use pre-extracted motif
        if intron.motifs is not None and intron.motifs.three_region:
            return intron.motifs.three_region

        # Default mode: extract from full sequences
        if intron.sequences is None or intron.sequences.seq is None:
            raise ValueError(f"Intron {intron.intron_id} has no sequence data")

        start, stop = self.three_coords
        intron_length = len(intron.sequences.seq)

        # Start is negative (positions before end)
        # Stop might be positive (positions after end)

        # Calculate intron indices
        if start < 0:
            intron_start = intron_length + start  # e.g., 100 + (-6) = 94
        else:
            intron_start = intron_length  # All from downstream

        if stop <= 0:
            intron_stop = intron_length + stop  # e.g., 100 + (-2) = 98
        else:
            intron_stop = intron_length  # Take to end of intron

        # Extract intron part
        if intron_start < intron_stop:
            intron_part = intron.sequences.seq[intron_start:intron_stop]
        else:
            intron_part = ""

        # If stop > 0, need downstream flank
        if stop > 0:
            downstream = intron.sequences.downstream_flank or ""
            downstream_part = (
                downstream[:stop] if len(downstream) >= stop else "N" * stop
            )
            return intron_part + downstream_part
        else:
            # All from intron
            return intron_part

    @staticmethod
    def _calculate_log_ratio(u12_score: float, u2_score: float) -> float:
        """
        Calculate log2 ratio of U12 to U2 scores.

        Port from: intronIC.py:2975-2980 (log_ratio)

        This is the core scoring metric: how much better does the sequence
        score with U12 matrices vs U2 matrices?

        Args:
            u12_score: Score from U12 PWM
            u2_score: Score from U2 PWM

        Returns:
            log2(u12_score / u2_score)

        Example:
            >>> _calculate_log_ratio(0.8, 0.2)
            2.0  # log2(4) = 2.0
        """
        return math.log2(u12_score / u2_score)


def score_and_normalize_introns(
    introns: Iterator[Intron],
    scorer: IntronScorer,
    scaler: "RobustScaler",
) -> Iterator[Intron]:
    """Score introns and normalize using pre-fitted scaler.

    This is the core streaming function for memory-efficient classification.
    It yields introns one at a time with both raw and z-scores populated,
    without accumulating results in memory.

    Used in true streaming mode with pre-trained models where the scaler
    is already fitted on reference data.

    Args:
        introns: Iterator of introns to score (with sequences or motifs populated)
        scorer: IntronScorer configured with PWM matrices
        scaler: Pre-fitted sklearn RobustScaler from the model bundle
                (obtained via ScoreNormalizer.get_frozen_scaler())

    Yields:
        Introns with both raw_score and z_score fields populated

    Example:
        >>> # Extract scaler from pre-trained model
        >>> normalizer = model_bundle['normalizer']
        >>> scaler = normalizer.get_frozen_scaler()
        >>>
        >>> # Score and normalize in single pass
        >>> for intron in score_and_normalize_introns(introns, scorer, scaler):
        ...     # Intron has both raw and z-scores
        ...     print(intron.scores.five_z_score)
        ...     # Memory freed after processing each intron

    Note:
        The scaler expects features in order: [five, bp, three]
        This matches the training order in ScoreNormalizer.
    """
    import numpy as np

    for intron in introns:
        # Score with PWMs (populates raw scores)
        scored = scorer.score_intron(intron)

        # Extract raw scores as numpy array for scaler
        # Order: [five, bp, three] - matches ScoreNormalizer training
        raw_scores = np.array(
            [
                [
                    scored.scores.five_raw_score,
                    scored.scores.bp_raw_score,
                    scored.scores.three_raw_score,
                ]
            ]
        )

        # Transform using frozen scaler (no fitting!)
        z_scores = scaler.transform(raw_scores)

        # Update intron with z-scores
        updated_scores = replace(
            scored.scores,
            five_z_score=float(z_scores[0, 0]),
            bp_z_score=float(z_scores[0, 1]),
            three_z_score=float(z_scores[0, 2]),
        )

        yield replace(scored, scores=updated_scores)


def score_and_normalize_batch(
    introns: list[Intron],
    scorer: IntronScorer,
    scaler: "RobustScaler",
) -> list[Intron]:
    """Score and normalize a batch of introns efficiently.

    Vectorized version of score_and_normalize_introns for better performance
    when processing batches (e.g., all introns from a chromosome).

    Args:
        introns: List of introns to score
        scorer: IntronScorer configured with PWM matrices
        scaler: Pre-fitted sklearn RobustScaler from the model bundle

    Returns:
        List of introns with both raw and z-scores populated

    Example:
        >>> # Process one chromosome's worth of introns
        >>> chromosome_introns = list(generate_introns(chromosome_genes))
        >>> scored_introns = score_and_normalize_batch(chromosome_introns, scorer, scaler)
    """
    import numpy as np

    if not introns:
        return []

    # Score all introns (collect raw scores)
    scored_introns = [scorer.score_intron(intron) for intron in introns]

    # Extract all raw scores as matrix
    raw_scores = np.array(
        [
            [
                intron.scores.five_raw_score,
                intron.scores.bp_raw_score,
                intron.scores.three_raw_score,
            ]
            for intron in scored_introns
        ]
    )

    # Vectorized transform
    z_scores = scaler.transform(raw_scores)

    # Update introns with z-scores
    result = []
    for i, intron in enumerate(scored_introns):
        updated_scores = replace(
            intron.scores,
            five_z_score=float(z_scores[i, 0]),
            bp_z_score=float(z_scores[i, 1]),
            three_z_score=float(z_scores[i, 2]),
        )
        result.append(replace(intron, scores=updated_scores))

    return result
