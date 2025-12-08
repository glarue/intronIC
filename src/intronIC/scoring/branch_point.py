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

from intronIC.scoring.pwm import PWM
from intronIC.core.intron import Intron


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
        search_region: The actual sequence searched (for bp_region_seq output)

    Example:
        >>> match = BranchPointMatch(
        ...     sequence="TACTAAC",
        ...     score=0.85,
        ...     position=-30,
        ...     start_in_region=20,
        ...     stop_in_region=27,
        ...     sequence_u2="CACAG",
        ...     score_u2=0.65,
        ...     search_region="ACGT..."*12
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
    search_region: str | None = None  # The actual search region sequence


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
        search_window: Tuple[int, int] = (-55, -5),
        five_coords: Tuple[int, int] = (-3, 9),
        three_coords: Tuple[int, int] = (-6, 4)
    ) -> BranchPointMatch | None:
        """
        Find the best branch point match in an intron.

        Port from: intronIC.py:2143-2178 (bp_score), 2944 (None handling)

        Supports two modes:
        - Default mode: Extracts BP region from full sequences (intron.sequences)
        - Streaming mode: Uses pre-extracted BP region (intron.motifs.bp_region)

        Args:
            intron: Intron object with sequence or motifs
            search_window: Tuple of (start, stop) positions relative to 3' end
                          Default: (-55, -5) matches original intronIC
            five_coords: 5' splice site scoring region coordinates (start, stop)
                        Used to exclude this region from BP search for short introns
                        Default: (-3, 9) matches original intronIC
            three_coords: 3' splice site scoring region coordinates (start, stop) relative to 3' end
                         Used to exclude this region from BP search for short introns
                         Default: (-6, 4) matches original intronIC

        Returns:
            BranchPointMatch with best-scoring sequence and position, or None
            if the search window is too small for the PWM length.

        Raises:
            ValueError: If intron has no sequence or motifs

        Example:
            >>> scorer = BranchPointScorer(u12_pwm, u2_pwm)
            >>> match = scorer.find_best_match(intron, search_window=(-55, -5))
            >>> if match:
            ...     print(f"Found {match.sequence} at position {match.position}")
            Found TACTAAC at position -30
        """
        # Check for streaming mode: use pre-extracted bp_region
        if intron.motifs is not None and intron.motifs.bp_region:
            search_region = intron.motifs.bp_region
            # In streaming mode, we use the search_window start as the actual position
            # since the bp_region was extracted with those exact coordinates
            actual_start_pos = search_window[0]  # Will be adjusted below for intron_length
            # For position calculation, we need an estimated intron length
            # We can infer it from the bp_region: if extracted with (-55, -5), region is 50bp
            # and actual intron_length would make position = actual_start_pos + match.start - intron_length
            # But in streaming mode, we store search_window[0] directly and it's already relative to 3' end
            use_streaming_mode = True
        elif intron.sequences is not None and intron.sequences.seq is not None:
            # Default mode: extract from full sequences
            search_region, actual_start_pos = self._extract_search_region(intron, search_window, five_coords, three_coords)
            use_streaming_mode = False
        else:
            raise ValueError(
                f"Intron {intron.intron_id} has no sequence data. "
                "Cannot search for branch point without sequences or motifs."
            )

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
        if use_streaming_mode:
            # In streaming mode, search_window coordinates are already relative to 3' end
            # So position = search_window[0] + match.start_in_region
            # e.g., search_window=(-55, -5), match at position 10 in region = -55 + 10 = -45
            position = search_window[0] + u12_match.start_in_region
        else:
            # Default mode: calculate from absolute positions
            # actual_start_pos is the absolute position in the intron where search region starts (after clamping)
            # u12_match.start_in_region is offset into the search region
            # absolute position in intron = actual_start_pos + u12_match.start_in_region
            # position relative to 3' end = absolute_position - intron_length
            intron_length = len(intron.sequences.seq)
            absolute_position = actual_start_pos + u12_match.start_in_region
            position = absolute_position - intron_length

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
        # Store the actual search_region so it can be saved to bp_region_seq
        return BranchPointMatch(
            sequence=u12_match.sequence,
            score=u12_match.score,
            position=position,
            start_in_region=u12_match.start_in_region,
            stop_in_region=u12_match.stop_in_region,
            sequence_u2=u2_match.sequence,
            score_u2=u2_score,  # Score of U12's sequence with U2 PWM
            search_region=search_region  # The actual sequence that was searched
        )

    def _extract_search_region(
        self,
        intron: Intron,
        search_window: Tuple[int, int],
        five_coords: Tuple[int, int],
        three_coords: Tuple[int, int]
    ) -> Tuple[str, int]:
        """
        Extract search region from intron sequence.

        Port from: intronIC.py:2527-2560 (_short_bp_adjust), 2607-2610

        For short introns, the search window is automatically adjusted to fit
        within the intron boundaries AND to exclude the 5' and 3' splice site
        scoring regions. This ensures we never overlap with other scoring regions.

        Args:
            intron: Intron object
            search_window: (start, stop) relative to 3' end (negative values)
            five_coords: (start, stop) for 5' splice site region
            three_coords: (start, stop) for 3' splice site region (relative to 3' end)

        Returns:
            Tuple of (search_region_sequence, actual_start_position) where:
            - search_region_sequence: Substring of intron sequence for the search region
            - actual_start_position: Absolute position in intron where search region starts (after clamping)

        Example:
            For 100bp intron with window (-55, -5), five_coords=(-3, 9), three_coords=(-6, 4):
            - start_pos = 100 + (-55) = 45
            - stop_pos = 100 + (-5) = 95
            - five_end = 9 (end of 5' region)
            - Final: start_pos = max(45, 9) = 45, stop_pos = min(95, 100) = 95
            - Returns intron.seq[45:95]

            For 47bp intron with same coords (short intron case):
            - start_pos = 47 + (-55) = -8 (INVALID!)
            - stop_pos = 47 + (-5) = 42
            - five_end = 9
            - Clamp: start_pos = max(-8, 9) = 9, stop_pos = min(42, 47) = 42
            - Returns intron.seq[9:42] (33bp search region, excluding 5' region only)
            - Note: Position 41 overlaps with 3' region [41-46], which is allowed
        """
        intron_length = len(intron.sequences.seq)

        # Convert relative positions to absolute indices
        # search_window values are negative (e.g., -55, -5)
        start_pos = intron_length + search_window[0]
        stop_pos = intron_length + search_window[1]

        # Calculate the boundaries of the 5' and 3' scoring regions
        # Port from: intronIC.py:2527-2560 (_short_bp_adjust)
        # 5' region: [0, five_coords[1]) - ends at five_coords[1]
        five_end = five_coords[1]  # e.g., 9 for (-3, 9)

        # 3' region: [intron_length + three_coords[0], intron_length) - starts at three_coords[0] from end
        three_start = intron_length + three_coords[0]  # e.g., for 47bp intron and (-6, 4): 47 + (-6) = 41

        # Clamp start position to:
        # 1. Stay within intron boundaries (>= 0)
        # 2. Not overlap with 5' scoring region (>= five_end)
        # Port from: intronIC.py:2608-2612 (_short_bp_adjust only adjusts start, not stop)
        start_pos = max(start_pos, five_end)

        # Note: Original code does NOT clamp stop position to avoid 3' region
        # BP region is allowed to overlap with 3' scoring region
        # Only clamp to intron boundaries
        stop_pos = min(stop_pos, intron_length)

        # Extract region (may be empty string if stop_pos <= start_pos)
        search_region = intron.sequences.seq[start_pos:stop_pos]

        # Return both the search region and its actual start position in the intron
        return search_region, start_pos

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
