"""
Branch point detection using sliding window PWM scoring.

This module implements the branch point search algorithm from intronIC,
which finds the optimal branch point sequence within a search region
upstream of the 3' splice site.

Port from: intronIC.py:2143-2178 (bp_score), 2097-2111 (sliding_window)

Algorithm:
1. Extract search region from intron (e.g., -55 to -5 relative to 3' end)
2. Use sliding window of PWM length to scan the region
3. Score each window position using PWM
4. Return the highest-scoring match with its position and sequence

Design:
- Immutable BranchPointMatch result object
- BranchPointScorer class encapsulates U12/U2 PWMs
- Clear separation between sequence extraction and scoring
"""

from dataclasses import dataclass
from typing import Tuple, Iterator

from scoring.pwm import PWM
from core.intron import Intron


@dataclass(frozen=True, slots=True)
class BranchPointMatch:
    """
    Result from branch point search.

    Attributes:
        sequence: The best-matching U12 branch point sequence
        score: PWM score for U12 sequence
        position: Position relative to intron 3' end (negative, e.g., -30)
        start_in_region: Start coordinate within search region (0-based)
        stop_in_region: Stop coordinate within search region (exclusive)
        sequence_u2: The best-matching U2 branch point sequence
        score_u2: PWM score for U2 sequence

    Example:
        >>> match = BranchPointMatch(
        ...     sequence="TACTAAC",
        ...     score=0.85,
        ...     position=-30,
        ...     start_in_region=20,
        ...     stop_in_region=27,
        ...     sequence_u2="CACAG",
        ...     score_u2=0.65
        ... )
        >>> match.sequence
        'TACTAAC'
        >>> match.sequence_u2
        'CACAG'
    """
    sequence: str  # U12 BP sequence
    score: float  # U12 score
    position: int  # Relative to 3' end (negative)
    start_in_region: int  # Position in search region
    stop_in_region: int  # Position in search region (exclusive)
    sequence_u2: str | None = None  # U2 BP sequence
    score_u2: float | None = None  # U2 score


class BranchPointScorer:
    """
    Find optimal branch point location using PWM scoring.

    Uses a sliding window approach to score all possible branch point
    positions within a search region, returning the best match.

    Port from: intronIC.py:2143-2178 (bp_score)

    Attributes:
        u12_pwm: U12-type branch point PWM (for scoring)
        u2_pwm: U2-type branch point PWM (for future use)
    """

    def __init__(self, u12_pwm: PWM, u2_pwm: PWM):
        """
        Initialize branch point scorer with PWMs.

        Args:
            u12_pwm: U12 branch point PWM (typically TACTAAC motif)
            u2_pwm: U2 branch point PWM (typically degenerate YNYURAY)
        """
        self.u12_pwm = u12_pwm
        self.u2_pwm = u2_pwm

    def find_best_match(
        self,
        intron: Intron,
        search_window: Tuple[int, int] = (-55, -5)
    ) -> BranchPointMatch | None:
        """
        Find the best branch point match in an intron.

        Port from: intronIC.py:2143-2178 (bp_score), 2944 (None handling)

        Args:
            intron: Intron object with sequence
            search_window: Tuple of (start, stop) positions relative to 3' end
                          Default: (-55, -5) matches original intronIC

        Returns:
            BranchPointMatch with best-scoring sequence and position, or None
            if the search window is too small for the PWM length.

        Raises:
            ValueError: If intron has no sequence

        Example:
            >>> scorer = BranchPointScorer(u12_pwm, u2_pwm)
            >>> match = scorer.find_best_match(intron, search_window=(-55, -5))
            >>> if match:
            ...     print(f"Found {match.sequence} at position {match.position}")
            Found TACTAAC at position -30
        """
        # Validate intron has sequence
        if intron.sequences is None or intron.sequences.seq is None:
            raise ValueError(
                f"Intron {intron.intron_id} has no sequence. "
                "Cannot search for branch point."
            )

        # Extract search region from intron
        search_region = self._extract_search_region(intron, search_window)

        # Check if search region is long enough for PWM
        # Port from: intronIC.py:2944 - returns None for too-short sequences
        window_size = self.u12_pwm.length
        if len(search_region) < window_size:
            # Return None instead of raising - caller will handle with pseudocount
            return None

        # Find best match in search region using U12 PWM
        # Port from: intronIC.py:2920-2944 (multi_matrix_score bp scoring)
        # Pass search_window[0] so scorer knows the genomic position of the search region
        u12_match = self._find_best_in_sequence(search_region, self.u12_pwm, search_window_start=search_window[0])

        # Calculate position relative to 3' end
        # search_window[0] is negative (e.g., -55)
        # u12_match.start_in_region is offset into region (e.g., 20)
        # position = start_of_region + offset_in_region
        position = search_window[0] + u12_match.start_in_region

        # Score the SAME sequence (U12's best match) with the U2 PWM
        # Port from: intronIC.py:3086-3095 (log_ratio using same bp_region_seq)
        # This ensures we get a proper log-odds ratio: log2(P(seq|U12) / P(seq|U2))
        # CRITICAL: BP PWMs use start_index=0, so use default seq_start_position=0
        u2_score = self.u2_pwm.score_sequence(u12_match.sequence)  # Use default seq_start_position=0

        # Also find U2's own best match for bp_seq_u2 (diagnostic/output purposes)
        # Port from: intronIC.py:3082-3084 (separate U2 BP sequence tracking)
        u2_match = self._find_best_in_sequence(search_region, self.u2_pwm, search_window_start=search_window[0])

        # Create combined match with U12 and U2 results
        # CRITICAL: score_u2 is the U2 score of the U12's best-match sequence (for log ratio)
        # sequence_u2 is U2's own best match (for informational purposes only)
        return BranchPointMatch(
            sequence=u12_match.sequence,
            score=u12_match.score,
            position=position,
            start_in_region=u12_match.start_in_region,
            stop_in_region=u12_match.stop_in_region,
            sequence_u2=u2_match.sequence,
            score_u2=u2_score  # Score of U12's sequence with U2 PWM
        )

    def _extract_search_region(
        self,
        intron: Intron,
        search_window: Tuple[int, int]
    ) -> str:
        """
        Extract search region from intron sequence.

        Port from: intronIC.py:2527-2560 (_short_bp_adjust), 2607-2610

        For short introns, the search window is automatically adjusted to fit
        within the intron boundaries. This ensures we search the maximum
        available region rather than failing for short introns.

        Args:
            intron: Intron object
            search_window: (start, stop) relative to 3' end (negative values)

        Returns:
            Substring of intron sequence for the search region

        Example:
            For 100bp intron with window (-55, -5):
            - start_pos = 100 + (-55) = 45
            - stop_pos = 100 + (-5) = 95
            - Returns intron.seq[45:95]

            For 50bp intron with window (-55, -5):
            - start_pos would be 50 + (-55) = -5 (INVALID!)
            - Clamp to 0: start_pos = 0
            - stop_pos = 50 + (-5) = 45
            - Returns intron.seq[0:45] (45bp search region)
        """
        intron_length = len(intron.sequences.seq)

        # Convert relative positions to absolute indices
        # search_window values are negative (e.g., -55, -5)
        start_pos = intron_length + search_window[0]
        stop_pos = intron_length + search_window[1]

        # Clamp start position to stay within intron boundaries
        # Port from: intronIC.py:2527-2560 (_short_bp_adjust)
        # If start_pos < 0, the window extends before the intron start
        if start_pos < 0:
            start_pos = 0

        # Extract region
        return intron.sequences.seq[start_pos:stop_pos]

    def _find_best_in_sequence(self, sequence: str, pwm: PWM, search_window_start: int) -> BranchPointMatch:
        """
        Find best-scoring subsequence using sliding window.

        Port from: intronIC.py:2143-2178 (bp_score)

        This is the core algorithm:
        1. Slide window of PWM length across sequence
        2. Score each window position
        3. Track best score, coordinates, and sequence
        4. Return best match

        Args:
            sequence: Search region sequence
            pwm: PWM to use for scoring (U12 or U2)
            search_window_start: Starting position of search region (e.g., -55)
                                Used to calculate seq_start_position for each window

        Returns:
            BranchPointMatch with best score

        Raises:
            ValueError: If sequence is shorter than PWM length
        """
        window_size = pwm.length

        # Validate sequence length
        if len(sequence) < window_size:
            raise ValueError(
                f"Sequence ({len(sequence)}bp) is too short for "
                f"PWM length ({window_size}bp)"
            )

        # Initialize tracking variables
        # Port from: intronIC.py:2158-2160
        best_score = None
        best_coords = None
        best_seq = None

        # Sliding window search
        # Port from: intronIC.py:2161-2177
        start = 0
        stop = window_size

        for sub_seq in self._sliding_window(sequence, window_size):
            # Score this window
            # Port from: intronIC.py:2164
            # CRITICAL: BP PWMs have start_index=0 and expect positions 0-11
            # The original bp_score() calls seq_score(sub_seq, matrix) with NO start_index,
            # which defaults to 0. We must do the same - seq_start_position=0 (the default).
            # DO NOT pass genomic position here!
            new_score = pwm.score_sequence(sub_seq)  # Use default seq_start_position=0
            new_coords = (start, stop)

            # Check if this is the best so far
            # Port from: intronIC.py:2166-2174
            if best_score is None or new_score > best_score:
                best_score = new_score
                best_coords = new_coords
                best_seq = sub_seq

            # Advance window
            # Port from: intronIC.py:2175-2177
            start += 1
            stop += 1

        # Return result
        # Port from: intronIC.py:2178
        return BranchPointMatch(
            sequence=best_seq,
            score=best_score,
            position=0,  # Will be updated by find_best_match
            start_in_region=best_coords[0],
            stop_in_region=best_coords[1]
        )

    @staticmethod
    def _sliding_window(sequence: str, window_size: int) -> Iterator[str]:
        """
        Generate sliding windows over sequence.

        Port from: intronIC.py:2097-2111 (sliding_window)

        This is a simple, efficient sliding window implementation that
        yields successive overlapping windows of fixed size.

        Original used itertools approach with tuples. We simplify to
        use string slicing for clarity.

        Args:
            sequence: Input sequence
            window_size: Size of sliding window

        Yields:
            Successive windows of size window_size

        Example:
            >>> list(_sliding_window("ABCDEF", 3))
            ['ABC', 'BCD', 'CDE', 'DEF']
        """
        # Yield windows from position 0 to len(seq) - window_size
        for i in range(len(sequence) - window_size + 1):
            yield sequence[i:i + window_size]
