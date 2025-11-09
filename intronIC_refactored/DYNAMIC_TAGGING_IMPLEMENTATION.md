# Dynamic Tagging Implementation Plan

## Overview

Dynamic tags are state markers added to intron names based on metadata already tracked in `IntronMetadata`. We already have all the data needed - we just need to populate the `dynamic_tags` set when we set the corresponding flags.

## Current State

All required metadata is already being tracked:
- ✅ `noncanonical: bool` - Set in `extraction/sequences.py` line 260
- ✅ `duplicate: Optional[str]` - Set in `extraction/filters.py` line 298
- ✅ `longest_isoform: bool` - Set in `extraction/filters.py` line 309
- ✅ `corrected: bool` - Field exists (not yet used)
- ✅ `correction_distance: Optional[int]` - Field exists (not yet used)
- ✅ `dynamic_tags: set[str]` - Field exists, ready to populate

## Tag Mapping

| Tag | Meaning | When Added | Current Code Location |
|-----|---------|------------|----------------------|
| `[n]` | Non-canonical | When `noncanonical = True` | `extraction/sequences.py:260` |
| `[d]` | Duplicate | When `duplicate` is set | `extraction/filters.py:298` |
| `[i]` | Not longest isoform | When `longest_isoform = False` | `extraction/filters.py:309` |
| `[c:N]` | Corrected by N bases | When boundary adjusted | Not yet implemented |
| `[e]` | Edge case | For unusual situations | Not yet needed |
| Terminal DNTs | Non-standard DNTs | e.g., `GC-AG`, `AT-AC` | Could add after scoring |

## Implementation

### Step 1: Update `extraction/sequences.py` for `[n]` tag

**Location:** Line 260
**Current code:**
```python
# Update intron metadata
intron.metadata.noncanonical = not is_canonical
```

**Updated code:**
```python
# Update intron metadata
intron.metadata.noncanonical = not is_canonical

# Add dynamic tag for non-canonical introns
if intron.metadata.noncanonical:
    intron.metadata.dynamic_tags.add('[n]')
```

### Step 2: Update `extraction/filters.py` for `[d]` and `[i]` tags

**Location A:** Line 298 (duplicate tagging)
**Current code:**
```python
else:
    # Duplicate found - reference the original
    intron.metadata.duplicate = region_idx[coord_key]['unique_id']
    intron.metadata.overlap = intron.metadata.duplicate
```

**Updated code:**
```python
else:
    # Duplicate found - reference the original
    intron.metadata.duplicate = region_idx[coord_key]['unique_id']
    intron.metadata.overlap = intron.metadata.duplicate
    # Add dynamic tag for duplicates
    intron.metadata.dynamic_tags.add('[d]')
```

**Location B:** Line 309 (longest isoform tagging)
**Current code:**
```python
if grandparent and grandparent in self.longest_isoforms:
    # Check if this intron's transcript is the "longest" (first seen) for its gene
    longest_transcript = self.longest_isoforms[grandparent]
    intron.metadata.longest_isoform = (parent == longest_transcript)
else:
    # No grandparent info, assume longest
    intron.metadata.longest_isoform = True
```

**Updated code:**
```python
if grandparent and grandparent in self.longest_isoforms:
    # Check if this intron's transcript is the "longest" (first seen) for its gene
    longest_transcript = self.longest_isoforms[grandparent]
    intron.metadata.longest_isoform = (parent == longest_transcript)
    # Add dynamic tag for non-longest isoforms
    if not intron.metadata.longest_isoform:
        intron.metadata.dynamic_tags.add('[i]')
else:
    # No grandparent info, assume longest
    intron.metadata.longest_isoform = True
```

### Step 3: (Optional) Boundary Correction `[c:N]` tag

This would be added if/when we implement boundary correction for U12 introns. The code would look like:

```python
# In boundary correction function (if implemented)
if boundaries_corrected:
    shift_distance = abs(new_start - old_start)
    intron.metadata.corrected = True
    intron.metadata.correction_distance = shift_distance
    intron.metadata.dynamic_tags.add(f'[c:{shift_distance}]')
```

**Note:** The original `u12_correction()` function exists at `intronIC.py:2298` but may not be implemented in refactored version yet.

### Step 4: (Optional) Terminal Dinucleotide tags

For non-standard terminal dinucleotides (e.g., `GC-AG`), we could add them after scoring:

```python
# After scoring, when we have terminal dinucleotides
if intron.sequences and intron.sequences.terminal_dinucleotides:
    dnts = intron.sequences.terminal_dinucleotides
    # Only add if non-standard (not GT-AG)
    if dnts not in ('GT-AG', 'AT-AC'):
        intron.metadata.dynamic_tags.add(dnts)
```

## Testing

After implementation, intron names should match original format:

**Before:**
```
HomSap-gene:ENSG00000196218@transcript:ENST00000355481-intron_69(104)
```

**After (with tags):**
```
HomSap-gene:ENSG00000196218@transcript:ENST00000355481-intron_69(104);[o:i];[i]
```

## Summary

**Minimal Implementation (Steps 1-2):** 3 simple edits
- Add `[n]` tag in `sequences.py`
- Add `[d]` tag in `filters.py`
- Add `[i]` tag in `filters.py`

**Optional Enhancements (Steps 3-4):**
- Boundary correction (if needed)
- Terminal dinucleotide tags (cosmetic)

The dynamic tags will then automatically appear in output via the `format_dynamic_tags()` function we already implemented in `writers.py`.
