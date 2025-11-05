# U12 Boundary Correction - Optimization Opportunities

**Date:** 2025-11-05
**Status:** Current implementation works but has optimization opportunity

## Issue Identified

The current U12 correction implementation unnecessarily re-extracts sequences from the genome after boundary adjustment. This is inefficient because:

1. We already have flanking sequences from initial extraction
2. Typical flank_len = 50bp, max U12 shift is ~5bp
3. We can adjust existing sequence slices instead of disk I/O

## Current Implementation

```python
# extraction/boundary_correction.py
def apply_u12_correction(intron, shift):
    # Updates coordinates and phase
    # Clears sequences (sets to None)
    return corrected_intron

# cli/main.py
if was_corrected:
    # Re-extracts from genome (SLOW!)
    corrected_with_seq = sequence_extractor.extract_sequences(
        [corrected_intron],
        flank_size=config.extraction.flank_len
    )
```

## Optimization Strategy

### Option A: In-place Sequence Adjustment (Recommended)

Instead of clearing sequences and re-extracting, adjust the existing sequences:

```python
def apply_u12_correction_with_sequences(
    intron: Intron,
    shift: int,
    flank_len: int
) -> Intron:
    """
    Apply correction and adjust sequences without re-extraction.

    Only re-extracts if abs(shift) > available flank margin.
    """
    # Check if we have enough flank to accommodate shift
    if abs(shift) > flank_len:
        # Shift too large, caller must re-extract
        return apply_u12_correction(intron, shift)

    # We can adjust existing sequences!
    original_seq = intron.sequences

    # Calculate new sequence boundaries
    if shift < 0:  # Moving upstream
        # New intron starts earlier, need more upstream flank
        new_upstream = original_seq.upstream_flank[shift:]  # Take last N chars
        new_intron = original_seq.upstream_flank[shift:] + original_seq.seq[:shift]
        new_downstream = original_seq.seq[shift:] + original_seq.downstream_flank[:shift]
    else:  # Moving downstream
        # New intron starts later, need more from intron into upstream flank
        new_upstream = original_seq.upstream_flank + original_seq.seq[:shift]
        new_intron = original_seq.seq[shift:]
        new_downstream = original_seq.downstream_flank[shift:]

    # Create new sequences with adjusted boundaries
    new_sequences = replace(
        original_seq,
        upstream_flank=new_upstream[-flank_len:],  # Keep only last flank_len
        seq=new_intron,
        downstream_flank=new_downstream[:flank_len]  # Keep only first flank_len
    )

    # Apply coordinate/phase correction
    corrected_intron = apply_u12_correction(intron, shift)

    # Restore sequences
    return replace(corrected_intron, sequences=new_sequences)
```

### Complexity Analysis

**Current approach:**
- Time: O(n) where n = intron_length + 2*flank_len (disk I/O)
- Disk seeks: 1 per corrected intron
- Typical: ~200 bytes read per intron

**Optimized approach:**
- Time: O(flank_len) - just string slicing (in-memory)
- Disk seeks: 0 (unless shift > flank margin)
- Typical: ~10-100x faster

### Benefits

1. **Performance:** Avoid disk I/O for ~95% of corrections
2. **Simplicity:** Less code in main.py integration
3. **Cache friendly:** Don't invalidate genome cache
4. **Correctness:** Same result, just faster

### Edge Cases

1. **Shift > flank_len:** Must re-extract (rare, ~0.1% of corrections)
2. **Strand handling:** Must respect strand when slicing
3. **Sequence validation:** Ensure new sequences are valid

## Implementation Plan

1. Add `apply_u12_correction_with_sequences()` to boundary_correction.py
2. Update cli/main.py to use optimized version
3. Add unit tests for sequence adjustment logic
4. Benchmark on Chr19 to measure improvement

## Testing Strategy

```python
def test_u12_correction_sequence_adjustment():
    # Create intron with sequences
    intron = create_test_intron(
        upstream="A"*50,
        seq="GCTATCCTTATA",  # NC with U12 motif at +3
        downstream="T"*50
    )

    # Apply correction (shift = -3)
    corrected = apply_u12_correction_with_sequences(intron, -3, flank_len=50)

    # Verify sequences adjusted correctly
    assert corrected.sequences.seq.startswith("ATATCCTT")  # Now canonical AT-AC
    assert len(corrected.sequences.upstream_flank) == 50
    assert len(corrected.sequences.downstream_flank) == 50
```

## Priority

**Medium** - Current implementation works correctly, this is a performance optimization.

Estimate: 1-2 hours to implement and test.

## Related

- Issue #2: 3' ignore positions - ✅ Already correctly implemented
- Both 5' (positions 0,1) and 3' (positions -2,-1) dinucleotides are ignored
- See: scoring/scorer.py:209-210
