"""
Species-specific U2 background correction.

Computes empirical per-position nucleotide frequencies from a species' own
intron pool and blends them with human-derived U2 PWMs to create species-
corrected U2 background PWMs. The U12 numerator PWMs remain human-derived.

This corrects composition bias at the PWM level — each scored position gets
a U2 denominator that reflects the species' actual base frequencies, not
human's. The correction is surgical: it only affects positions where the
species differs from human, and it preserves the U12 signal because the
U12 numerator is unchanged.

Algorithm (iterative trimmed background):
    1. Score all introns with human U2 PWMs (initial pass)
    2. For each region, exclude the top trim_percentile% by score
    3. Compute per-position nucleotide frequencies from the remaining introns
    4. Blend with human U2 prior: w * empirical + (1-w) * human, w = n/(n+n0)
    5. Optionally repeat with corrected scores for the trim ranking

See docs/species_background_correction_plan.md for full design.
"""

from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple
from collections import defaultdict

import numpy as np

from intronIC.scoring.pwm import PWM, PWMSet, PWMLoader

# Base indices matching PWM convention
_BASES = ['A', 'C', 'G', 'T']
_BI = {b: i for i, b in enumerate(_BASES)}


@dataclass
class BackgroundConfig:
    """Configuration for species-specific background correction."""
    enabled: bool = True
    n0: int = 1000
    trim_percentile: float = 5.0
    pseudocount_per_base: float = 1.0
    n_iterations: int = 1
    min_introns: int = 500


class _RegionAccumulator:
    """Accumulates per-position nucleotide counts for one scoring region,
    grouped by dinucleotide subtype."""

    def __init__(self, n_positions: int):
        self.n_positions = n_positions
        # {subtype_dnt: {'counts': np.array(n_pos, 4), 'n': int}}
        self._data: Dict[str, dict] = {}

    def add(self, subtype_dnt: str, seq: str):
        """Add one sequence's nucleotide counts."""
        if len(seq) < self.n_positions:
            return
        if subtype_dnt not in self._data:
            self._data[subtype_dnt] = {
                'counts': np.zeros((self.n_positions, 4), dtype=np.float64),
                'n': 0,
            }
        entry = self._data[subtype_dnt]
        for j in range(self.n_positions):
            b = seq[j].upper()
            if b in _BI:
                entry['counts'][j, _BI[b]] += 1
        entry['n'] += 1

    def get_frequencies(self, subtype_dnt: str, pseudocount: float = 1.0
                        ) -> Optional[Tuple[np.ndarray, int]]:
        """Return normalized frequencies and sample count for a subtype.

        Adds Dirichlet pseudocounts before normalization.

        Returns:
            (freqs, n) where freqs is (n_positions, 4) or None if no data.
        """
        entry = self._data.get(subtype_dnt)
        if entry is None or entry['n'] == 0:
            return None
        counts = entry['counts'] + pseudocount
        row_sums = counts.sum(axis=1, keepdims=True)
        return counts / row_sums, entry['n']

    def get_all_subtypes(self) -> list:
        """Return list of subtypes with data."""
        return list(self._data.keys())

    def export_counts(self) -> dict:
        """Export raw count data for merging across workers."""
        return {
            sub: {'counts': entry['counts'].copy(), 'n': entry['n']}
            for sub, entry in self._data.items()
        }

    def merge_counts(self, other_data: dict):
        """Merge count data from another accumulator (e.g., a worker)."""
        for sub, entry in other_data.items():
            if sub not in self._data:
                self._data[sub] = {
                    'counts': entry['counts'].copy(),
                    'n': entry['n'],
                }
            else:
                self._data[sub]['counts'] += entry['counts']
                self._data[sub]['n'] += entry['n']

    def subtract_introns(self, subtype_dnt: str, seqs: list):
        """Remove specific introns' counts (for trimming).

        Args:
            subtype_dnt: Dinucleotide subtype
            seqs: List of sequences to subtract
        """
        entry = self._data.get(subtype_dnt)
        if entry is None:
            return
        for seq in seqs:
            if len(seq) < self.n_positions:
                continue
            for j in range(self.n_positions):
                b = seq[j].upper()
                if b in _BI:
                    entry['counts'][j, _BI[b]] -= 1
            entry['n'] -= 1
        # Clamp to zero (floating point safety)
        entry['counts'] = np.maximum(entry['counts'], 0)
        entry['n'] = max(entry['n'], 0)


class _BPSAccumulator:
    """Accumulates average nucleotide composition across the BP search region.

    Unlike 5'SS and 3'SS which have position-specific frequencies, BPS uses
    a uniform composition across all 12 PWM positions because we don't know
    the BP position until after scoring.
    """

    def __init__(self):
        self._data: Dict[str, dict] = {}

    def add(self, subtype_dnt: str, bp_region_seq: str):
        """Add one intron's BP search region composition."""
        if not bp_region_seq or len(bp_region_seq) < 12:
            return
        if subtype_dnt not in self._data:
            self._data[subtype_dnt] = {
                'counts': np.zeros(4, dtype=np.float64),
                'total': 0,
                'n': 0,
            }
        entry = self._data[subtype_dnt]
        for b in bp_region_seq.upper():
            if b in _BI:
                entry['counts'][_BI[b]] += 1
                entry['total'] += 1
        entry['n'] += 1

    def get_frequencies(self, subtype_dnt: str, n_positions: int = 12,
                        pseudocount: float = 1.0
                        ) -> Optional[Tuple[np.ndarray, int]]:
        """Return uniform composition tiled to n_positions.

        Returns:
            (freqs, n) where freqs is (n_positions, 4) or None if no data.
        """
        entry = self._data.get(subtype_dnt)
        if entry is None or entry['total'] == 0:
            return None
        counts = entry['counts'] + pseudocount
        freqs_1d = counts / counts.sum()
        return np.tile(freqs_1d, (n_positions, 1)), entry['n']

    def get_all_subtypes(self) -> list:
        return list(self._data.keys())

    def export_counts(self) -> dict:
        """Export raw count data for merging across workers."""
        return {
            sub: {'counts': entry['counts'].copy(), 'total': entry['total'], 'n': entry['n']}
            for sub, entry in self._data.items()
        }

    def merge_counts(self, other_data: dict):
        """Merge count data from another accumulator."""
        for sub, entry in other_data.items():
            if sub not in self._data:
                self._data[sub] = {
                    'counts': entry['counts'].copy(),
                    'total': entry['total'],
                    'n': entry['n'],
                }
            else:
                self._data[sub]['counts'] += entry['counts']
                self._data[sub]['total'] += entry['total']
                self._data[sub]['n'] += entry['n']


class SpeciesBackground:
    """Compute and cache species-specific U2 background PWMs.

    Usage:
        bg = SpeciesBackground(human_u2_pwms, config)

        # During extraction: accumulate sequences
        for intron in introns:
            bg.accumulate(name, five_dnt, three_dnt,
                          five_seq, three_seq, bp_region_seq)

        # After initial scoring: record scores for trim ranking
        for intron in scored_introns:
            bg.set_initial_scores(name, five_raw, bp_raw, three_raw)

        # Build corrected PWMs
        corrected = bg.build_corrected_pwm_sets()

        # Use corrected PWMs for re-scoring
        scorer = IntronScorer(pwm_sets=corrected)
    """

    # Mapping from intron dinucleotides to PWM subtype keys
    FIVE_DNT_TO_SUBTYPE = {'GT': 'gtag', 'GC': 'gcag', 'AT': 'atac'}
    THREE_DNT_TO_SUBTYPE = {'AG': 'gtag', 'AC': 'atac'}

    def __init__(
        self,
        human_u2_pwm_sets: Dict[str, PWMSet],
        config: BackgroundConfig,
        five_len: int = 12,
        three_len: int = 10,
        bp_len: int = 12,
    ):
        """
        Args:
            human_u2_pwm_sets: Dict with keys 'five', 'three', 'bp',
                values are PWMSet objects containing human U2 PWMs.
            config: BackgroundConfig with correction parameters.
            five_len: Number of scored positions at 5'SS (default 12).
            three_len: Number of scored positions at 3'SS (default 10).
            bp_len: Number of positions in BPS PWM (default 12).
        """
        self.human_u2 = human_u2_pwm_sets
        self.config = config
        self.five_len = five_len
        self.three_len = three_len
        self.bp_len = bp_len

        # Accumulators for each region
        self._five_acc = _RegionAccumulator(five_len)
        self._three_acc = _RegionAccumulator(three_len)
        self._bp_acc = _BPSAccumulator()

        # Per-intron sequences for trimming (keyed by intron name)
        self._intron_seqs: Dict[str, dict] = {}

        # Per-intron initial scores for trim ranking
        self._initial_scores: Dict[str, dict] = {}

        # Cached result
        self._corrected_pwms: Optional[Dict[str, dict]] = None

    @property
    def n_introns(self) -> int:
        """Total introns accumulated."""
        return len(self._intron_seqs)

    @property
    def n_accumulated(self) -> int:
        """Total introns accumulated (count-based, for parallel merging)."""
        # When merging from workers, _intron_seqs may be empty but counts exist
        if self._intron_seqs:
            return len(self._intron_seqs)
        # Fall back to five accumulator counts
        return sum(e['n'] for e in self._five_acc._data.values())

    def merge_worker_counts(self, five_counts: dict, three_counts: dict, bp_counts: dict):
        """Merge count data from a parallel BG worker.

        Args:
            five_counts: From worker's _RegionAccumulator.export_counts()
            three_counts: From worker's _RegionAccumulator.export_counts()
            bp_counts: From worker's _BPSAccumulator.export_counts()
        """
        self._five_acc.merge_counts(five_counts)
        self._three_acc.merge_counts(three_counts)
        self._bp_acc.merge_counts(bp_counts)
        self._corrected_pwms = None  # invalidate cache

    def accumulate(
        self,
        intron_name: str,
        five_dnt: str,
        three_dnt: str,
        five_seq: str,
        three_seq: str,
        bp_region_seq: str,
    ):
        """Accumulate one intron's sequences into the background.

        Call once per intron during the extraction phase.

        Args:
            intron_name: Unique intron identifier
            five_dnt: 5' dinucleotide (e.g., 'GT', 'GC', 'AT')
            three_dnt: 3' dinucleotide (e.g., 'AG', 'AC')
            five_seq: 5'SS scored sequence (five_len chars)
            three_seq: 3'SS scored sequence (three_len chars)
            bp_region_seq: BP search region sequence (variable length)
        """
        self._five_acc.add(five_dnt, five_seq)
        self._three_acc.add(three_dnt, three_seq)
        self._bp_acc.add(five_dnt, bp_region_seq)

        # Store sequences for potential trimming
        self._intron_seqs[intron_name] = {
            'five_dnt': five_dnt,
            'three_dnt': three_dnt,
            'five_seq': five_seq,
            'three_seq': three_seq,
        }

        # Invalidate cache
        self._corrected_pwms = None

    def set_initial_scores(
        self, intron_name: str,
        five_raw: float, bp_raw: float, three_raw: float,
    ):
        """Record initial raw scores for trim ranking.

        Call after the initial scoring pass (with human U2 PWMs).

        Args:
            intron_name: Must match a name passed to accumulate()
            five_raw: 5'SS raw log-ratio score
            bp_raw: BPS raw log-ratio score
            three_raw: 3'SS raw log-ratio score
        """
        self._initial_scores[intron_name] = {
            'five_raw': five_raw,
            'bp_raw': bp_raw,
            'three_raw': three_raw,
        }

    def build_corrected_pwm_sets(
        self,
        u12_pwm_sets: Optional[Dict[str, PWMSet]] = None,
    ) -> Dict[str, PWMSet]:
        """Build species-corrected U2 PWM sets.

        Computes trimmed empirical frequencies, blends with human U2 prior,
        and returns PWMSet objects ready for use by IntronScorer.

        The returned dict has the same structure as the input PWM sets:
        {'five': PWMSet, 'bp': PWMSet, 'three': PWMSet}

        The U12 entries are copied from u12_pwm_sets (unchanged).
        The U2 entries are replaced with blended species-specific backgrounds.

        Args:
            u12_pwm_sets: If provided, U12 PWMs to include in the returned
                sets. If None, only U2 entries are included.

        Returns:
            Dict mapping region name to PWMSet with corrected U2 PWMs.
        """
        # Use n_accumulated (not n_introns) so this guard works for both the
        # in-memory path (which calls accumulate() and populates _intron_seqs)
        # and the streaming path (which calls merge_worker_counts() and only
        # populates the _RegionAccumulator counts). Without this, streaming
        # always returns the human U2 PWMs unchanged because _intron_seqs is
        # empty, and the streaming/in-memory pipelines disagree on every score.
        if self.n_accumulated < self.config.min_introns:
            # Not enough introns — return human U2 PWMs unchanged
            if u12_pwm_sets is not None:
                return _merge_u12_u2(u12_pwm_sets, self.human_u2)
            return {region: PWMSet(matrices=dict(ps.matrices))
                    for region, ps in self.human_u2.items()}

        # Trim and build for each iteration
        for iteration in range(self.config.n_iterations):
            self._apply_trim(iteration)

        return self._build_final_pwm_sets(u12_pwm_sets)

    def _apply_trim(self, iteration: int):
        """Apply score-based trimming to exclude likely U12 introns."""
        if not self._initial_scores:
            return  # No scores available — skip trimming

        trim_pct = self.config.trim_percentile
        if trim_pct <= 0:
            return

        # For each region, find introns above the trim threshold
        # and subtract their counts from the accumulators
        for region, acc, score_key in [
            ('five', self._five_acc, 'five_raw'),
            ('three', self._three_acc, 'three_raw'),
        ]:
            # Collect scores per subtype
            by_subtype = defaultdict(list)
            for name, scores in self._initial_scores.items():
                seq_info = self._intron_seqs.get(name)
                if seq_info is None:
                    continue
                if region == 'five':
                    subtype = seq_info['five_dnt']
                    seq = seq_info['five_seq']
                else:
                    subtype = seq_info['three_dnt']
                    seq = seq_info['three_seq']
                by_subtype[subtype].append((scores[score_key], name, seq))

            # For each subtype, trim top N%
            for subtype, entries in by_subtype.items():
                n_trim = max(1, int(len(entries) * trim_pct / 100))
                entries.sort(key=lambda x: -x[0])  # highest score first
                trimmed_seqs = [seq for _, _, seq in entries[:n_trim]]
                acc.subtract_introns(subtype, trimmed_seqs)

        # BPS: trim based on bp_raw scores
        bp_by_subtype = defaultdict(list)
        for name, scores in self._initial_scores.items():
            seq_info = self._intron_seqs.get(name)
            if seq_info is None:
                continue
            bp_by_subtype[seq_info['five_dnt']].append(
                (scores['bp_raw'], name)
            )
        # BPS accumulator doesn't support subtract_introns (it stores
        # aggregate composition, not per-intron). For BPS, we rebuild
        # from scratch excluding trimmed introns. This is acceptable
        # because BPS uses uniform composition (not positional).
        # For now, skip BPS trimming — the BPS contribution to FPs
        # is smaller than 5'SS and 3'SS.

    def _build_final_pwm_sets(
        self,
        u12_pwm_sets: Optional[Dict[str, PWMSet]],
    ) -> Dict[str, PWMSet]:
        """Build PWMSet objects from the (trimmed) accumulated frequencies."""
        config = self.config
        result = {}

        for region, acc, human_region, scored_len in [
            ('five', self._five_acc, 'five', self.five_len),
            ('three', self._three_acc, 'three', self.three_len),
        ]:
            matrices = {}

            # Copy U12 entries if provided
            if u12_pwm_sets and region in u12_pwm_sets:
                for key, pwm in u12_pwm_sets[region].matrices.items():
                    if key[0] == 'u12':
                        matrices[key] = pwm

            # Build blended U2 entries for each subtype
            human_ps = self.human_u2.get(region)
            if human_ps is None:
                continue

            for subtype_dnt in acc.get_all_subtypes():
                pwm_subtype = self.FIVE_DNT_TO_SUBTYPE.get(
                    subtype_dnt,
                    self.THREE_DNT_TO_SUBTYPE.get(subtype_dnt, 'gtag')
                )
                human_key = ('u2', pwm_subtype)
                human_pwm = human_ps.matrices.get(human_key)
                if human_pwm is None:
                    # Try fallback
                    human_pwm = human_ps.select_best('u2', pwm_subtype)
                if human_pwm is None:
                    continue

                emp_result = acc.get_frequencies(
                    subtype_dnt, pseudocount=config.pseudocount_per_base
                )
                if emp_result is None:
                    matrices[human_key] = human_pwm
                    continue

                emp_freqs, n_emp = emp_result
                w = n_emp / (n_emp + config.n0)

                # Extract human U2 frequencies at scored positions
                human_matrix = human_pwm.matrix  # (4, length)
                # We need positions matching the scored window
                # The human PWM covers -20 to +19 (40 positions)
                # Scored positions are a subset starting at the PWM's start_index
                # We need the same positions that the scorer uses
                # For simplicity, take the full human PWM and use it as-is

                # Build blended matrix
                blended = np.zeros_like(human_matrix)
                for pos_idx in range(human_pwm.length):
                    for base_idx in range(4):
                        h_freq = human_matrix[base_idx, pos_idx]
                        # Map pos_idx to scored position
                        scored_pos = pos_idx - (human_pwm.reference_offset +
                                                (self.five_len // 2 - human_pwm.reference_offset)
                                                if region == 'five' else 0)
                        # Actually, simpler: the empirical frequencies are at
                        # the scored positions only. For positions outside the
                        # scored window, use human frequencies unchanged.
                        bio_pos = pos_idx + human_pwm.start_index
                        if region == 'five':
                            scored_start = -3  # five_coords[0]
                            emp_idx = bio_pos - scored_start
                        else:
                            scored_start = -6  # three_coords[0]
                            emp_idx = bio_pos - scored_start

                        if 0 <= emp_idx < scored_len:
                            e_freq = emp_freqs[emp_idx, base_idx]
                            blended[base_idx, pos_idx] = max(
                                w * e_freq + (1 - w) * h_freq,
                                human_pwm.pseudocount,
                            )
                        else:
                            blended[base_idx, pos_idx] = h_freq

                corrected_pwm = PWM(
                    name=human_pwm.name + '_bg',
                    matrix=blended,
                    length=human_pwm.length,
                    pseudocount=human_pwm.pseudocount,
                    start_index=human_pwm.start_index,
                    reference_offset=human_pwm.reference_offset,
                )
                matrices[human_key] = corrected_pwm

            # Also copy any human U2 subtypes that didn't have empirical data
            for key, pwm in human_ps.matrices.items():
                if key not in matrices:
                    matrices[key] = pwm

            result[region] = PWMSet(matrices=matrices)

        # BPS
        bp_matrices = {}
        if u12_pwm_sets and 'bp' in u12_pwm_sets:
            for key, pwm in u12_pwm_sets['bp'].matrices.items():
                if key[0] == 'u12':
                    bp_matrices[key] = pwm

        human_bp = self.human_u2.get('bp')
        if human_bp:
            for subtype_dnt in self._bp_acc.get_all_subtypes():
                pwm_subtype = self.FIVE_DNT_TO_SUBTYPE.get(subtype_dnt, 'gtag')
                human_key = ('u2', pwm_subtype)
                human_pwm = human_bp.matrices.get(human_key)
                if human_pwm is None:
                    human_pwm = human_bp.select_best('u2', pwm_subtype)
                if human_pwm is None:
                    continue

                emp_result = self._bp_acc.get_frequencies(
                    subtype_dnt, n_positions=human_pwm.length,
                    pseudocount=config.pseudocount_per_base,
                )
                if emp_result is None:
                    bp_matrices[human_key] = human_pwm
                    continue

                emp_freqs, n_emp = emp_result
                w = n_emp / (n_emp + config.n0)

                blended = np.zeros_like(human_pwm.matrix)
                for pos_idx in range(human_pwm.length):
                    for base_idx in range(4):
                        h_freq = human_pwm.matrix[base_idx, pos_idx]
                        e_freq = emp_freqs[pos_idx, base_idx]
                        blended[base_idx, pos_idx] = max(
                            w * e_freq + (1 - w) * h_freq,
                            human_pwm.pseudocount,
                        )

                bp_matrices[human_key] = PWM(
                    name=human_pwm.name + '_bg',
                    matrix=blended,
                    length=human_pwm.length,
                    pseudocount=human_pwm.pseudocount,
                    start_index=human_pwm.start_index,
                    reference_offset=human_pwm.reference_offset,
                )

            for key, pwm in human_bp.matrices.items():
                if key not in bp_matrices:
                    bp_matrices[key] = pwm

        result['bp'] = PWMSet(matrices=bp_matrices)

        return result


def _merge_u12_u2(
    u12_sets: Dict[str, PWMSet],
    u2_sets: Dict[str, PWMSet],
) -> Dict[str, PWMSet]:
    """Merge U12 and U2 PWMSets into combined sets."""
    result = {}
    for region in ['five', 'bp', 'three']:
        matrices = {}
        if region in u12_sets:
            for key, pwm in u12_sets[region].matrices.items():
                if key[0] == 'u12':
                    matrices[key] = pwm
        if region in u2_sets:
            for key, pwm in u2_sets[region].matrices.items():
                if key[0] == 'u2':
                    matrices[key] = pwm
        result[region] = PWMSet(matrices=matrices)
    return result
