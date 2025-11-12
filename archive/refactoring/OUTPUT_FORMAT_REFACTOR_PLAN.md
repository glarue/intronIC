# Output Format Refactor Plan

**Date:** 2025-11-09
**Branch:** `claude/investigate-ref-performance-011CUt9YD4yNb4AzXF2TKAU9`
**Goal:** Achieve output format parity with original intronIC implementation

---

## Table of Contents

1. [Overview](#overview)
2. [Current vs Target Format](#current-vs-target-format)
3. [Missing Components](#missing-components)
4. [Implementation Plan](#implementation-plan)
5. [Files to Modify](#files-to-modify)
6. [Testing Strategy](#testing-strategy)

---

## Overview

The refactored intronIC currently produces simplified output that is missing several metadata fields present in the original implementation. This plan outlines the changes needed to achieve full parity.

### Format Discrepancy Example

**Original:**
```
HomSap-gene:ENSG00000196218@transcript:ENST00000355481-intron_69(104);[o:i];[i]
```

**Refactored (current):**
```
ENST00000397910_1(83)
```

---

## Current vs Target Format

### Intron Name Format

#### Original Format
```
{species_abbrev}-{grandparent}@{parent}-intron_{index}({family_size}){omit_tag}{dynamic_tags}
```

**Components:**
- `species_abbrev`: 3-char genus + 3-char species (e.g., "HomSap" from "Homo sapiens")
- `grandparent`: Gene ID, often with "gene:" prefix (e.g., "gene:ENSG00000196218")
- `@`: Separator between gene and transcript
- `parent`: Transcript ID, often with "transcript:" prefix (e.g., "transcript:ENST00000355481")
- `-intron_`: Literal separator
- `index`: Intron ordinal position (1-based)
- `family_size`: Total introns in transcript
- `omit_tag`: `;[o:{code}]` if omitted (code = s/a/n/i/v)
- `dynamic_tags`: `;[tag1];[tag2]...` for various states

#### Refactored Format (current)
```
[{species_name}_]{parent}_{index}({family_size})
```

**Issues:**
- Full species name instead of abbreviation
- Missing grandparent (gene ID)
- Missing @ separator
- Missing "gene:" and "transcript:" prefixes
- Using underscore instead of "-intron_"
- Missing omission tags
- Missing dynamic tags

### Metadata File (.meta.iic) Fields

#### Original Fields (tab-delimited)
```
name  rel_score  dnts  motif_schematic  bp_context  length  parent  grandparent
index  family_size  frac_pos  phase  type_id  feature
```

#### Refactored Fields (current)
```
name  rel_score  dnts  motif_schematic  bp_context  length  parent  grandparent
index  family_size  frac_pos  phase  type_id  feature  attributes
```

**Issues:**
- `motif_schematic`: Currently returns '.' instead of formatted schematic
- `bp_context`: Not properly formatted with brackets around BP sequence
- Missing proper tag representation in name field

---

## Missing Components

### 1. Species Abbreviation

**Original Implementation:**
```python
# Inferred from usage: first 3 chars of genus + first 3 chars of species
# "Homo sapiens" → "HomSap"
```

**Needed:**
- Function to generate abbreviation from full species name
- Handle edge cases (single-word names, short words)

### 2. Tagging System

#### Omission Tags

**Original Implementation** (line 629-632):
```python
if self.omitted:
    omit_tag = ';[o:{}]'.format(self.omitted)
else:
    omit_tag = ''
```

**Omission Codes** (line 687-693):
- `s`: short (length < minimum)
- `a`: ambiguous sequence (non-ACTG characters)
- `n`: noncanonical (non-standard terminal dinucleotides)
- `v`: coordinate overlap
- `i`: not in longest isoform

#### Dynamic Tags

**Original Implementation** (line 633-636):
```python
if self.dynamic_tag:
    dyn_tag = ';{}'.format(';'.join(sorted(self.dynamic_tag)))
else:
    dyn_tag = ''
```

**Dynamic Tag Types:**
- `[c:N]`: Boundary corrected by N bases (line 2347)
- `[d]`: Duplicate-related tag (lines 3589, 3639, 3879)
- `[e]`: Error/edge case (line 442)
- `[n]`: Non-canonical marker (line 4789)
- `[i]`: Isoform-related (line 4791)
- Terminal dinucleotides can also be added (line 4832)

**Refactored Model:**
Currently only tracks:
- `omitted: OmissionCode` (s/a/n/i/v)
- `corrected: bool` (no shift distance tracked)

**Missing:**
- `dynamic_tag: set` attribute to store tags
- Shift distance for correction tag

### 3. Motif Schematic

**Original Implementation** (line 725-742):
```python
def motif_string(self, exonic=3):
    """Returns a schematic of all scored motifs."""
    five_boundary = '{}|{}'.format(
        self.upstream_flank[-exonic:], self.five_display_seq)
    three_boundary = '{}|{}'.format(
        self.three_display_seq, self.downstream_flank[:exonic])
    try:
        bps_display_seq = '/'.join([self.bp_seq, self.bp_seq_u2])
    except:
        bps_display_seq = None
    schematic_bits = [five_boundary, bps_display_seq, three_boundary]
    schematic_bits = [b for b in schematic_bits if b is not None]
    schematic_string = '...'.join(schematic_bits)
    return schematic_string
```

**Format Example:**
```
AAG|GTCGGGGCTT...TACTAAC/CACAG...TTTAG|TCC
```

**Components:**
- `upstream_flank[-3:]`: Last 3 exonic bases
- `|`: Exon-intron boundary marker
- `five_display_seq`: First 10bp of intron
- `...`: Separator
- `bp_seq/bp_seq_u2`: U12 and U2 branch point sequences
- `...`: Separator
- `three_display_seq`: Intron from bp_coords[1] to end
- `|`: Intron-exon boundary marker
- `downstream_flank[:3]`: First 3 exonic bases

**Required Attributes:**
- `five_display_seq`: First 10bp of intron (line 2593-2595)
- `three_display_seq`: From bp_coords[1] to end (line 2596)
- `bp_seq_u2`: U2 branch point sequence

**Refactored Model:**
Currently missing these display sequences entirely.

### 4. Branch Point Context

**Original Implementation** (line 744-752):
```python
def bps_context(self):
    try:
        context = annotate(
            self.bp_region_seq,
            *self.bp_relative_coords) + self.three_display_seq
    except:
        context = None
    return context
```

**Helper Function** (line 3208-3211):
```python
def annotate(string, start, stop):
    """Add brackets around substring."""
    annotated = string[:start] + '[' + string[start:stop] + ']' + string[stop:]
    return annotated
```

**Format Example:**
```
TTGACAGGCAGTGATAT[TACTAAC]GACTGAGTTTAG
```

The BP sequence is wrapped in brackets within the bp_region_seq, then the three_display_seq is appended.

**Required Attributes:**
- `bp_relative_coords`: Tuple of (start, stop) positions of BP within bp_region_seq

**Refactored Model:**
Currently only stores `bp_seq` string, not its relative coordinates.

### 5. Prefix Handling

**Original Behavior:**
- Gene IDs often have "gene:" prefix
- Transcript IDs often have "transcript:" prefix
- These come from GFF3 annotation parsing

**Refactored Behavior:**
The `clean_names` feature strips these prefixes by default (enabled by default).

**Location:** `cli/args.py` line 102-113
```python
output_group.add_argument(
    '--clean_names',
    action='store_true',
    default=True,
    help='Remove "transcript:" and "gene:" prefixes from IDs (default: True)'
)
output_group.add_argument(
    '--no_clean_names',
    dest='clean_names',
    action='store_false',
    help='Keep "transcript:" and "gene:" prefixes in IDs'
)
```

**Issue:**
The original implementation appears to **preserve** these prefixes in output, but the refactored version strips them by default. We need to determine correct behavior.

---

## Implementation Plan

### Phase 1: Data Model Extensions

**File:** `intronIC_refactored/core/intron.py`

#### 1.1 Extend IntronSequences
Add missing sequence attributes:

```python
@dataclass(frozen=True, slots=True)
class IntronSequences:
    # ... existing fields ...

    # NEW FIELDS:
    five_display_seq: Optional[str] = None  # First 10bp of intron
    three_display_seq: Optional[str] = None  # From bp_coords[1] to end
    bp_seq_u2: Optional[str] = None  # U2 branch point sequence
    bp_relative_coords: Optional[tuple[int, int]] = None  # BP position within bp_region_seq
```

**Rationale:**
These fields are required for proper motif_string and bps_context generation.

#### 1.2 Extend IntronMetadata
Add dynamic tag tracking:

```python
@dataclass(slots=True)
class IntronMetadata:
    # ... existing fields ...

    # NEW FIELDS:
    dynamic_tags: set[str] = field(default_factory=set)  # [c:N], [d], [e], [n], [i], etc.
    correction_distance: Optional[int] = None  # Distance for [c:N] tag
```

**Rationale:**
The original uses a `dynamic_tag` set to track various state markers. The refactored version only has a boolean `corrected` flag, losing information about the correction distance.

### Phase 2: Sequence Processing Updates

**File:** `intronIC_refactored/sequences/extractor.py`

#### 2.1 Populate Display Sequences
Modify sequence extraction to populate new fields:

```python
def _extract_scored_sequences(...) -> IntronSequences:
    # ... existing code ...

    # NEW: Populate display sequences
    five_display_length = 10
    five_display_seq = intron_seq[:five_display_length]
    three_display_seq = intron_seq[bp_coords[1]:]  # From bp end to intron end

    return IntronSequences(
        seq=intron_seq,
        # ... other fields ...
        five_display_seq=five_display_seq,
        three_display_seq=three_display_seq,
        bp_relative_coords=bp_coords  # Store relative coordinates
    )
```

#### 2.2 Store U2 Branch Point
Modify PWM scoring to capture U2 branch point:

**File:** `intronIC_refactored/scoring/pwm.py`

```python
def score_branch_point(...) -> BranchPointResult:
    # ... existing code scores both U12 and U2 ...

    # Currently only returns best U12 BP
    # MODIFY to return both:
    return BranchPointResult(
        bp_seq=best_u12_seq,
        bp_seq_u2=best_u2_seq,  # NEW
        bp_coords=best_u12_coords,
        # ...
    )
```

### Phase 3: Tagging System Implementation

**File:** `intronIC_refactored/filtering/quality_control.py`

#### 3.1 Update Omission Tagging
Modify `omit_check()` to use consistent codes:

```python
def omit_check(intron: Intron, ...) -> Intron:
    """Check intron for omission criteria."""

    omission_reason = None

    if intron.length < min_length:
        omission_reason = 's'  # short
    elif has_ambiguous_sequence(intron):
        omission_reason = 'a'  # ambiguous
    elif not allow_noncanonical and intron.metadata.noncanonical:
        omission_reason = 'n'  # noncanonical
    elif longest_only and not intron.metadata.longest_isoform:
        omission_reason = 'i'  # not in longest isoform
    elif not allow_overlap and intron.metadata.overlap:
        omission_reason = 'v'  # coordinate overlap

    if omission_reason:
        # Update metadata with omission code
        intron.metadata.omitted = omission_reason

    return intron
```

#### 3.2 Add Dynamic Tag Helpers
Create helper functions for tag management:

```python
def add_dynamic_tag(intron: Intron, tag: str) -> Intron:
    """Add a dynamic tag to intron metadata."""
    intron.metadata.dynamic_tags.add(tag)
    return intron

def add_correction_tag(intron: Intron, distance: int) -> Intron:
    """Add boundary correction tag with distance."""
    intron.metadata.correction_distance = distance
    intron.metadata.dynamic_tags.add(f'[c:{distance}]')
    return intron

def add_noncanonical_tag(intron: Intron) -> Intron:
    """Add non-canonical marker tag."""
    if intron.metadata.noncanonical:
        intron.metadata.dynamic_tags.add('[n]')
    return intron
```

**File:** `intronIC_refactored/sequences/boundary_correction.py`

Modify boundary correction to track distance:

```python
def correct_u12_boundaries(...) -> Intron:
    # ... existing correction logic ...

    if boundaries_adjusted:
        shift_distance = abs(new_start - old_start)  # or new_end - old_end
        intron = add_correction_tag(intron, shift_distance)

    return intron
```

### Phase 4: Name Generation Updates

**File:** `intronIC_refactored/file_io/writers.py`

#### 4.1 Species Abbreviation Function

```python
def generate_species_abbreviation(species_name: str) -> str:
    """
    Generate 3+3 character abbreviation from species name.

    Args:
        species_name: Full species name (e.g., "homo_sapiens")

    Returns:
        6-character abbreviation (e.g., "HomSap")

    Examples:
        >>> generate_species_abbreviation("homo_sapiens")
        'HomSap'
        >>> generate_species_abbreviation("drosophila_melanogaster")
        'DroMel'
        >>> generate_species_abbreviation("c_elegans")
        'CEle'
    """
    # Split on underscore or space
    parts = species_name.replace('_', ' ').split()

    if len(parts) >= 2:
        # Standard binomial nomenclature
        genus = parts[0]
        species = parts[1]

        # First 3 chars of each, capitalized
        genus_abbrev = genus[:3].capitalize()
        species_abbrev = species[:3].capitalize()

        return genus_abbrev + species_abbrev
    elif len(parts) == 1:
        # Single name (unusual case)
        name = parts[0]
        if len(name) >= 6:
            return name[:6].capitalize()
        else:
            # Pad with capital first letter
            return (name[:1].upper() + name[1:]).ljust(6, 'X')
    else:
        # Empty or invalid
        return "XXXXXX"
```

#### 4.2 Tag Formatting Functions

```python
def format_omission_tag(omitted: Optional[str]) -> str:
    """
    Format omission tag.

    Args:
        omitted: Omission code (s/a/n/i/v) or None

    Returns:
        Formatted tag string like ';[o:s]' or empty string
    """
    if omitted:
        return f';[o:{omitted}]'
    return ''

def format_dynamic_tags(tags: set[str]) -> str:
    """
    Format dynamic tags.

    Args:
        tags: Set of tag strings (may or may not have brackets)

    Returns:
        Formatted tag string like ';[c:5];[d]' or empty string

    Examples:
        >>> format_dynamic_tags({'[c:5]', '[d]'})
        ';[c:5];[d]'
        >>> format_dynamic_tags(set())
        ''
    """
    if not tags:
        return ''

    # Ensure all tags have brackets
    formatted_tags = []
    for tag in sorted(tags):
        if not tag.startswith('['):
            tag = f'[{tag}]'
        formatted_tags.append(tag)

    return ';' + ';'.join(formatted_tags)
```

#### 4.3 Updated Name Generation

```python
def _generate_name(
    self,
    intron: Intron,
    species_name: Optional[str],
    simple_name: bool
) -> str:
    """
    Generate intron name matching original format.

    Format: {species_abbrev}-{grandparent}@{parent}-intron_{index}({family_size}){omit_tag}{dynamic_tags}

    Args:
        intron: Intron object
        species_name: Full species name (e.g., 'homo_sapiens')
        simple_name: Use simplified format (excludes gene info)

    Returns:
        Formatted intron name

    Examples:
        >>> _generate_name(intron, "homo_sapiens", False)
        'HomSap-gene:ENSG00000196218@transcript:ENST00000355481-intron_69(104)'
        >>> _generate_name(omitted_intron, "homo_sapiens", False)
        'HomSap-gene:ENSG00000196218@transcript:ENST00000355481-intron_69(104);[o:i]'
    """
    if simple_name:
        # Simple format: species-i_{unique_num}{tags}
        # NOTE: We don't have unique_num in refactored version yet
        # Fallback to basic format
        species_abbrev = generate_species_abbreviation(species_name) if species_name else "XXXXXX"
        omit_tag = format_omission_tag(intron.metadata.omitted) if intron.metadata else ''
        dyn_tag = format_dynamic_tags(intron.metadata.dynamic_tags) if intron.metadata else ''
        return f"{species_abbrev}-{intron.intron_id}{omit_tag}{dyn_tag}"

    # Full format
    if not intron.metadata:
        # Fallback if no metadata
        return intron.intron_id

    # Species abbreviation
    species_abbrev = generate_species_abbreviation(species_name) if species_name else "XXXXXX"

    # Gene ID (grandparent)
    grandparent = intron.metadata.grandparent if intron.metadata.grandparent else "?"

    # Transcript ID (parent)
    parent = intron.metadata.parent if intron.metadata.parent else "?"

    # Index and family size
    index = intron.metadata.index if intron.metadata.index is not None else "?"
    family_size = intron.metadata.family_size if intron.metadata.family_size is not None else "?"

    # Tags
    omit_tag = format_omission_tag(intron.metadata.omitted)
    dyn_tag = format_dynamic_tags(intron.metadata.dynamic_tags)

    # Build name
    name = f"{species_abbrev}-{grandparent}@{parent}-intron_{index}({family_size}){omit_tag}{dyn_tag}"

    return name
```

### Phase 5: Motif and Context Generation

**File:** `intronIC_refactored/file_io/writers.py`

#### 5.1 Helper Functions

```python
def annotate_sequence(sequence: str, start: int, stop: int) -> str:
    """
    Add brackets around substring.

    Args:
        sequence: Full sequence
        start: Start position (0-based)
        stop: Stop position (0-based, exclusive)

    Returns:
        Annotated sequence with [brackets]

    Examples:
        >>> annotate_sequence("ABCDEFGH", 2, 5)
        'AB[CDE]FGH'
    """
    return sequence[:start] + '[' + sequence[start:stop] + ']' + sequence[stop:]

def generate_motif_schematic(intron: Intron, exonic: int = 3) -> str:
    """
    Generate motif schematic string.

    Format: {exon_3bp}|{5'_10bp}...{bp_u12}/{bp_u2}...{3'_display}|{exon_3bp}

    Args:
        intron: Intron object with sequences
        exonic: Number of exonic bases to show (default: 3)

    Returns:
        Motif schematic string or '.' if sequences missing

    Examples:
        >>> generate_motif_schematic(intron)
        'AAG|GTCGGGGCTT...TACTAAC/CACAG...TTTAG|TCC'
    """
    if not intron.sequences:
        return '.'

    seqs = intron.sequences

    # Check for required sequences
    if not all([
        seqs.upstream_flank,
        seqs.five_display_seq,
        seqs.three_display_seq,
        seqs.downstream_flank
    ]):
        return '.'

    # Five prime boundary: exonic|intron
    five_boundary = f"{seqs.upstream_flank[-exonic:]}|{seqs.five_display_seq}"

    # Three prime boundary: intron|exonic
    three_boundary = f"{seqs.three_display_seq}|{seqs.downstream_flank[:exonic]}"

    # Branch point display: U12/U2
    bps_display = None
    if seqs.bp_seq and seqs.bp_seq_u2:
        bps_display = f"{seqs.bp_seq}/{seqs.bp_seq_u2}"
    elif seqs.bp_seq:
        bps_display = seqs.bp_seq

    # Assemble schematic
    schematic_parts = [five_boundary]
    if bps_display:
        schematic_parts.append(bps_display)
    schematic_parts.append(three_boundary)

    return '...'.join(schematic_parts)

def generate_bp_context(intron: Intron) -> str:
    """
    Generate branch point context string.

    Format: bp_region with [brackets] around bp_seq + three_display_seq

    Args:
        intron: Intron object with sequences

    Returns:
        BP context string or '.' if sequences missing

    Examples:
        >>> generate_bp_context(intron)
        'TTGACAGGCAGTGATAT[TACTAAC]GACTGAGTTTAG'
    """
    if not intron.sequences:
        return '.'

    seqs = intron.sequences

    # Check for required sequences
    if not all([
        seqs.bp_region_seq,
        seqs.bp_relative_coords,
        seqs.three_display_seq
    ]):
        return '.'

    try:
        start, stop = seqs.bp_relative_coords
        annotated_bp_region = annotate_sequence(seqs.bp_region_seq, start, stop)
        context = annotated_bp_region + seqs.three_display_seq
        return context
    except Exception:
        return '.'
```

#### 5.2 Update MetaWriter

Modify `write_intron()` method:

```python
def write_intron(
    self,
    intron: Intron,
    species_name: Optional[str] = None,
    simple_name: bool = False,
    null: str = '.'
) -> None:
    """Write a single intron's metadata."""
    # ... existing validation ...

    # Generate intron name (UPDATED)
    name = self._generate_name(intron, species_name, simple_name)

    # ... rel_score, dnts (existing) ...

    # Motif schematic (UPDATED)
    motif = generate_motif_schematic(intron)

    # Branch point context (UPDATED)
    bp_context = generate_bp_context(intron)

    # ... rest of fields (existing) ...
```

### Phase 6: Prefix Handling Decision

**Issue:** The original code appears to preserve "gene:" and "transcript:" prefixes in output, but the refactored version strips them by default with `--clean_names`.

**Options:**

1. **Change default behavior** - Set `clean_names=False` by default to match original
2. **Keep current behavior** - Document that this is an intentional improvement
3. **Make it conditional** - Only strip prefixes in internal storage, preserve in output

**Recommendation:** Option 3 - Internal storage without prefixes, but output should preserve them to match original format exactly.

**Implementation:**

```python
def _generate_name(...):
    # ... existing code ...

    # Add prefixes for output (even if internally stripped)
    grandparent_formatted = grandparent
    if grandparent and not grandparent.startswith('gene:'):
        grandparent_formatted = f"gene:{grandparent}"

    parent_formatted = parent
    if parent and not parent.startswith('transcript:'):
        parent_formatted = f"transcript:{parent}"

    name = f"{species_abbrev}-{grandparent_formatted}@{parent_formatted}-intron_{index}({family_size}){omit_tag}{dyn_tag}"
```

**Note:** This assumes prefixes were stripped during parsing. If they're already present, don't re-add.

### Phase 7: Integration Points

**Files that need updates to propagate new data:**

1. **Sequence Extractor** (`sequences/extractor.py`)
   - Populate `five_display_seq`, `three_display_seq` in IntronSequences
   - Store `bp_relative_coords` from BP scoring

2. **PWM Scorer** (`scoring/pwm.py`)
   - Return both U12 and U2 branch point sequences
   - Modify BranchPointResult dataclass

3. **Boundary Corrector** (`sequences/boundary_correction.py`)
   - Add correction tags with distance
   - Update IntronMetadata.dynamic_tags

4. **Quality Control** (`filtering/quality_control.py`)
   - Populate omission codes consistently
   - Add dynamic tags for various states

5. **Writers** (`file_io/writers.py`)
   - Implement new name generation
   - Implement motif and context generation
   - Update all output formats (META, SEQ, SCORE, BED)

---

## Files to Modify

### Primary Changes

1. **`core/intron.py`**
   - Add fields to IntronSequences and IntronMetadata
   - Estimated changes: +10 lines

2. **`sequences/extractor.py`**
   - Populate display sequences and BP coordinates
   - Estimated changes: +15 lines

3. **`scoring/pwm.py`**
   - Return U2 BP sequence alongside U12
   - Estimated changes: +10 lines

4. **`filtering/quality_control.py`**
   - Add dynamic tagging functions
   - Update omission logic
   - Estimated changes: +40 lines

5. **`sequences/boundary_correction.py`**
   - Track correction distance
   - Add correction tags
   - Estimated changes: +10 lines

6. **`file_io/writers.py`**
   - Implement species abbreviation
   - Implement tag formatting
   - Implement new name generation
   - Implement motif/context generation
   - Update MetaWriter, ScoreWriter, SequenceWriter, BedWriter
   - Estimated changes: +150 lines

### Supporting Changes

7. **`cli/args.py`**
   - Potentially adjust `--clean_names` default (if needed)
   - Estimated changes: +5 lines (documentation)

---

## Testing Strategy

### Unit Tests

1. **Species Abbreviation**
   ```python
   def test_species_abbreviation():
       assert generate_species_abbreviation("homo_sapiens") == "HomSap"
       assert generate_species_abbreviation("drosophila_melanogaster") == "DroMel"
       assert generate_species_abbreviation("c_elegans") == "CEle"
   ```

2. **Tag Formatting**
   ```python
   def test_omission_tag():
       assert format_omission_tag('s') == ';[o:s]'
       assert format_omission_tag(None) == ''

   def test_dynamic_tags():
       tags = {'[c:5]', '[d]'}
       result = format_dynamic_tags(tags)
       assert result == ';[c:5];[d]' or result == ';[d];[c:5]'  # order may vary
   ```

3. **Name Generation**
   ```python
   def test_name_generation():
       intron = create_test_intron(
           grandparent="gene:ENSG001",
           parent="transcript:ENST001",
           index=5,
           family_size=10
       )
       name = _generate_name(intron, "homo_sapiens", False)
       assert name.startswith("HomSap-")
       assert "@" in name
       assert "-intron_5(10)" in name
   ```

4. **Motif Schematic**
   ```python
   def test_motif_schematic():
       intron = create_test_intron_with_sequences()
       schematic = generate_motif_schematic(intron)
       assert '|' in schematic  # Exon-intron boundary
       assert '...' in schematic  # Separator
       assert '/' in schematic  # U12/U2 BP separator (if both present)
   ```

### Integration Tests

1. **Full Pipeline Test**
   - Run on Chr19 test data
   - Compare output format to original implementation
   - Check that all metadata fields are populated

2. **Regression Test**
   - Generate outputs from known introns
   - Compare to golden reference files
   - Ensure no information loss

### Manual Validation

1. **Visual Inspection**
   - Compare `.meta.iic` files side-by-side (original vs refactored)
   - Verify species abbreviation correctness
   - Check tag formatting

2. **Downstream Compatibility**
   - Ensure outputs can be read by downstream tools expecting original format
   - Test that R/Python parsing scripts work

---

## Implementation Order

**Recommended sequence:**

1. ✅ **Phase 1** - Data model extensions (IntronSequences, IntronMetadata)
2. ✅ **Phase 2** - Sequence processing updates (populate new fields)
3. ✅ **Phase 3** - Tagging system (omission codes, dynamic tags)
4. ✅ **Phase 4** - Name generation (species abbreviation, tag formatting)
5. ✅ **Phase 5** - Motif and context generation (schematic, bp_context)
6. ✅ **Phase 6** - Prefix handling decision (gene:/transcript:)
7. ✅ **Phase 7** - Integration and testing

**Estimated Time:**
- Implementation: 4-6 hours
- Testing: 2-3 hours
- Total: 6-9 hours

---

## Risk Assessment

### Low Risk

- Species abbreviation generation (pure function, easy to test)
- Tag formatting functions (pure functions)
- Data model extensions (backward compatible, optional fields)

### Medium Risk

- Name generation changes (central to output, needs careful testing)
- Motif schematic generation (depends on multiple sequence fields)
- BP context generation (coordinate manipulation)

### High Risk

- Prefix handling decision (affects parsing/compatibility)
- Integration across multiple modules (sequence extractor → PWM scorer → writers)
- Backward compatibility with existing pipelines

### Mitigation Strategies

1. **Extensive testing** - Unit tests for all new functions
2. **Gradual rollout** - Test on small datasets first
3. **Comparison script** - Automated diff between original and refactored outputs
4. **Feature flag** - Optional `--legacy_format` flag to enable old behavior if needed

---

## Success Criteria

Output format refactor is complete when:

1. ✅ Intron names match original format exactly
2. ✅ All metadata fields are populated correctly
3. ✅ Motif schematics display properly
4. ✅ Branch point context is formatted correctly
5. ✅ Omission and dynamic tags appear in output
6. ✅ Species abbreviation is correct
7. ✅ All tests pass
8. ✅ Manual comparison shows parity with original

---

## Notes

### Design Decisions

1. **Composition over modification** - Extending dataclasses rather than replacing them
2. **Backward compatibility** - New fields are optional, don't break existing code
3. **Clear separation** - Formatting logic in writers, not in data models
4. **Testability** - Pure functions for formatting enable easy testing

### Future Enhancements

1. **Configurable format** - Allow users to choose output format
2. **JSON output** - Modern alternative to tab-delimited format
3. **Compressed metadata** - Optional minimal format for large datasets

---

## Appendix: Original Code References

### Key Functions (intronIC.py)

- `get_name()`: Line 622-646
- `get_label()`: Line 648-669
- `omit_check()`: Line 671-723
- `motif_string()`: Line 725-742
- `bps_context()`: Line 744-752
- `add_tags()`: Line 1897-1952
- `annotate()`: Line 3208-3211
- `output_format()`: Line 4202-4282

### Omission Tag Codes

```python
omit_tags = {
    'short': 's',
    'ambiguous sequence': 'a',
    'noncanonical': 'n',
    'coordinate overlap': 'v',
    'not in longest isoform': 'i'
}
```

### Dynamic Tag Assignments

- Line 442: `[e]` - Edge case
- Line 2347: `[c:N]` - Boundary correction with distance
- Line 3589, 3639, 3879: `[d]` - Duplicate marker
- Line 4789: `[n]` - Non-canonical marker
- Line 4791: `[i]` - Isoform-related
- Line 4832: Terminal dinucleotides

---

**END OF REFACTOR PLAN**
