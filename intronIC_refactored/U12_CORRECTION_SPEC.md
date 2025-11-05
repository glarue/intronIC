# U12 Boundary Correction Algorithm Specification

## Purpose
Detects and corrects misannotated U12-type intron boundaries by searching for strong U12-specific 5' splice site motifs in the surrounding sequence.

## When Applied
- **Only** for introns with non-canonical terminal dinucleotides
- **Only** when `--no_nc_ss_adjustment` flag is NOT set (enabled by default)
- Applied during intron sequence assignment phase

## Algorithm

### Step 1: Define Search Parameters
```python
up_n = 5        # bases of upstream flank to search
down_n = 12     # bases into intron to search
```

### Step 2: Define U12 5' Splice Site Motif
**Strict motif** (used by default):
```
[AG]TATCC([ACTG]T|T[ACTG])
```
This matches patterns like:
- ATATCCAT, ATATCCTT, GTATCCAT, etc.

**Lax motif** (commented out in original, not currently used):
```
[AG]TATC[CT]
```

### Step 3: Extract Search Region
```
search_region = upstream_flank[-5:] + intron_seq[:12]
```
- Last 5 bp of upstream exon
- First 12 bp of intron
- Total: 17 bp search window

### Step 4: Search for Motif
```python
match = strict_motif.search(search_region)
if not match:
    return False  # No correction possible
```

### Step 5: Calculate Shift
```python
match_index = match.start()  # Position where motif starts (0-17)
shift = match_index - up_n    # Shift relative to expected position (5)
```

**Examples:**
- Motif at position 5 → shift = 0 (no correction needed)
- Motif at position 3 → shift = -2 (move boundaries 2bp upstream)
- Motif at position 7 → shift = +2 (move boundaries 2bp downstream)

If shift == 0, no correction is needed, return False.

### Step 6: Apply Correction
**Update intron attributes:**
```python
intron.corrected = shift           # Store shift amount
intron.start += shift              # Adjust start coordinate
intron.stop += shift               # Adjust stop coordinate (same shift on both ends)
intron.phase = _shift_phase(intron.phase, shift)  # Rotate phase
intron.dynamic_tag.add('[c:{}]'.format(shift))    # Add correction tag
```

**Strand handling:**
- For negative strand: `shift *= -1` before applying to coordinates
- This ensures coordinates always increase left-to-right on reference

### Step 7: Phase Shift Logic
```python
def _shift_phase(phase, shift):
    phases = deque([0, 1, 2])
    try:
        index = phases.index(int(phase))
    except ValueError:  # e.g. '.' for non-coding exons
        return phase
    phases.rotate(-shift)
    return phases[index]
```

**Examples:**
- Phase 0, shift -2 → Phase 1
- Phase 1, shift -2 → Phase 2
- Phase 2, shift -2 → Phase 0

### Step 8: Re-assign Sequences
After correction, the intron's sequences must be re-extracted using the new coordinates:
```python
intron = assign_seqs(
    intron,
    region_seq,
    int_flank_size,
    five_score_coords,
    three_score_coords,
    bp_coords
)
```

## Impact on Downstream Processing

### Annotation File Correction (Optional)
The `correct_annotation()` function can modify the original GFF3/GTF file:
- Adjusts exon boundaries to match corrected intron coordinates
- Adds `;shift:[start/stop]:[shift_value]` tags to modified lines
- This ensures the annotation reflects the biological reality

### Scoring
- Corrected introns are scored with their new sequences
- The `[c:shift]` tag appears in the intron name in output files

### Statistics
- Corrected introns are tracked and reported separately
- Helps identify annotation quality issues

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--no_nc_ss_adjustment` | False | Disables U12 correction when set |
| upstream search | 5 bp | Bases of upstream flank to search |
| downstream search | 12 bp | Bases into intron to search |
| motif | `[AG]TATCC([ACTG]T\|T[ACTG])` | Strict U12 5'SS pattern |

## Implementation Locations (Original)

1. **Function definition:** `intronIC.py:2298` - `u12_correction()`
2. **Called from:** `intronIC.py:2692` - During sequence assignment
3. **Argument parsing:** `intronIC.py:1107` - `--no_nc_ss_adjustment`
4. **Config setup:** `intronIC.py:4497` - `U12_NC_SS_ADJUSTMENT = not args.no_nc_ss_adjustment`

## Refactored Implementation Plan

### Location in Refactored Code
**Module:** `extraction/boundary_correction.py` (new module)

### Key Functions
```python
def correct_u12_boundary(
    intron: Intron,
    upstream_flank: str,
    downstream_intron: str,
    strict: bool = True
) -> Optional[int]:
    """
    Search for misannotated U12 boundaries.

    Args:
        intron: Intron object to check
        upstream_flank: Upstream exonic sequence (last 5bp used)
        downstream_intron: Beginning of intron sequence (first 12bp used)
        strict: Use strict motif (default) vs lax

    Returns:
        Shift amount if correction found, None otherwise
    """

def apply_boundary_correction(
    intron: Intron,
    shift: int,
    genome: GenomeSequence
) -> Intron:
    """
    Apply boundary correction to intron.

    Args:
        intron: Original intron
        shift: Number of bases to shift (can be negative)
        genome: Genome sequence for fetching updated sequences

    Returns:
        New Intron object with corrected coordinates
    """
```

### Integration Points
1. **CLI:** Add `--no_nc_ss_adjustment` flag to `cli/args.py`
2. **Config:** Add to `cli/config.py` in `ExtractionConfig`
3. **Extraction:** Call from `extraction/extractor.py` after initial sequence assignment
4. **Testing:** Add unit tests for motif search, shift calculation, phase rotation

## Edge Cases

1. **Multiple motif matches:** First match is used
2. **Shift == 0:** No correction applied
3. **Canonical introns:** Not checked (only non-canonical)
4. **No motif found:** Intron remains unchanged
5. **Negative strand:** Shift direction reversed for coordinate adjustment
6. **Phase = '.':** Returned unchanged (non-coding)

## Example

**Input:**
```
Intron: chr1:1000-2000:+
Terminal dinucleotides: GC-AG (non-canonical 5')
Upstream flank (last 5): CGATC
Intron seq (first 12): GCTATCCTTATA
```

**Search region:**
```
CGATCGCTATCCTTATA
```

**Motif search:**
```
Pattern: [AG]TATCC([ACTG]T|T[ACTG])
Match: GTATCCTT at position 3
Shift: 3 - 5 = -2
```

**Correction:**
```
Original: chr1:1000-2000:+
Corrected: chr1:998-1998:+
Shift: -2
Tag: [c:-2]
```

**New terminal dinucleotides:** GT-AG (canonical U12)

## References
- Original implementation: `intronIC/intronIC.py:2298-2353`
- Moyer et al. (2020) NAR - describes U12 intron characteristics
- U12-specific 5' splice site: ATATCCTT consensus
