# Streaming Annotation Processing Implementation Plan

**Date:** 2025-11-21
**Goal:** Reduce peak memory from 12 GB → 2-3 GB during annotation processing
**Strategy:** Two-pass streaming with lightweight contig indexing

---

## Problem Statement

### Current Memory Profile
```
1. Load genome:                 ~250 MB
2. Parse full annotation:       ~10 GB   (genes, transcripts, exons, indexes)
3. Generate all introns:        +2 GB    (2.1M intron objects)
4. Peak memory:                 ~12 GB   ← Too high!
5. Free annotation hierarchy:   -10 GB
6. Continue with introns:       ~2 GB
```

**Issue:** 12 GB peak is unavoidable with current approach - we must build the entire annotation hierarchy in memory before generating introns.

### Target Memory Profile
```
1. Load genome:                 ~250 MB
2. Build lightweight index:     ~20 MB   (contig -> line numbers)
3. For each contig:
   - Load contig annotations:   ~500 MB  (one contig's hierarchy)
   - Generate contig introns:   ~100 MB
   - Free contig hierarchy:     ~0 MB
   - Accumulate introns:        ~2 GB    (grows incrementally)
4. Peak memory:                 ~2-3 GB  ← 80% reduction!
```

**Improvement:** Never load full hierarchy - process one contig at a time!

---

## Solution: Two-Pass Streaming with Pre-Indexing

### Overview

**Pass 1: Quick Index (Fast)**
- Scan annotation file once
- Record which line numbers belong to each contig
- Build lightweight index: `{contig: [line_numbers]}`
- Memory: ~10-20 MB (just integers!)
- Time: ~30 seconds for human genome

**Pass 2: Stream by Contig (Memory-Efficient)**
- For each contig:
  - Extract only that contig's lines from file
  - Build hierarchy for just that contig
  - Generate introns
  - Free contig hierarchy immediately
- Memory: ~500 MB per contig max
- Time: Similar to current (maybe 10-20% slower due to seeks)

### Key Insight

**Index stores line numbers, not data:**
```python
# Lightweight index (~20 MB)
contig_index = {
    'chr1': [120, 121, 122, ..., 45230],      # ~20K integers = 80 KB
    'chr2': [45231, 45232, ..., 89450],       # ~18K integers = 72 KB
    ...
    'chrY': [2100450, 2100451, ..., 2105000], # ~500 integers = 2 KB
}

# NOT storing full data (would be 10 GB):
# 'chr1': [<Gene>, <Gene>, ...] ❌
```

---

## Implementation Plan

### Phase 1: Create Annotation Indexer

**New file:** `file_io/annotation_index.py`

**Key components:**
1. `build_contig_index(annotation_file)` - Fast indexing pass
2. `ContigIndex` class - Store and query index
3. `extract_contig_lines(annotation_file, line_numbers)` - Extract specific lines
4. Optional: Save/load index from cache file

**Implementation:**

```python
# file_io/annotation_index.py

from typing import Dict, List, Set
from pathlib import Path
import gzip
import pickle
from collections import defaultdict
from smart_open import open as smart_open

class ContigIndex:
    """
    Lightweight index mapping contigs to annotation file line numbers.

    Memory: ~10-20 MB for human genome (vs 10+ GB for full hierarchy)
    """

    def __init__(self):
        self.contig_to_lines: Dict[str, List[int]] = defaultdict(list)
        self.total_lines: int = 0

    def add_line(self, contig: str, line_num: int) -> None:
        """Record that a line belongs to a contig."""
        self.contig_to_lines[contig].append(line_num)

    def get_contigs(self) -> List[str]:
        """Get all contigs in the index."""
        return sorted(self.contig_to_lines.keys())

    def get_line_numbers(self, contig: str) -> List[int]:
        """Get line numbers for a specific contig."""
        return self.contig_to_lines.get(contig, [])

    def memory_estimate(self) -> int:
        """Estimate memory usage in bytes."""
        # Each integer is ~8 bytes in Python
        total_entries = sum(len(lines) for lines in self.contig_to_lines.values())
        return total_entries * 8

    def save(self, cache_file: Path) -> None:
        """Save index to cache file."""
        with open(cache_file, 'wb') as f:
            pickle.dump(self, f, protocol=pickle.HIGHEST_PROTOCOL)

    @classmethod
    def load(cls, cache_file: Path) -> 'ContigIndex':
        """Load index from cache file."""
        with open(cache_file, 'rb') as f:
            return pickle.load(f)


def build_contig_index(annotation_file: str, messenger=None) -> ContigIndex:
    """
    Build lightweight index of which lines belong to which contigs.

    This is a fast scanning pass that only records line numbers,
    not full annotation data. Memory usage: ~10-20 MB.

    Args:
        annotation_file: Path to GFF3/GTF file (can be gzipped)
        messenger: Optional messenger for progress logging

    Returns:
        ContigIndex mapping contigs to line numbers

    Examples:
        >>> index = build_contig_index("annotations.gff3.gz")
        >>> index.get_contigs()
        ['chr1', 'chr2', ..., 'chrY']
        >>> chr1_lines = index.get_line_numbers('chr1')
        >>> len(chr1_lines)
        45230
    """
    index = ContigIndex()

    if messenger:
        messenger.info(f"Building annotation index: {annotation_file}")

    with smart_open(annotation_file, 'rt') as f:
        for line_num, line in enumerate(f):
            line = line.strip()

            # Skip comments and empty lines
            if not line or line.startswith('#'):
                continue

            # Extract contig from first column (GFF3/GTF format)
            fields = line.split('\t')
            if len(fields) < 8:
                continue

            contig = fields[0]
            index.add_line(contig, line_num)

            # Progress every 100K lines
            if messenger and (line_num + 1) % 100000 == 0:
                messenger.log_only(f"  Indexed {line_num + 1:,} lines", level="debug")

    index.total_lines = line_num + 1

    if messenger:
        num_contigs = len(index.get_contigs())
        memory_mb = index.memory_estimate() / (1024 * 1024)
        messenger.info(
            f"Index built: {num_contigs} contigs, "
            f"{index.total_lines:,} lines, "
            f"~{memory_mb:.1f} MB memory"
        )

    return index


def extract_contig_lines(
    annotation_file: str,
    line_numbers: List[int]
) -> List[str]:
    """
    Extract specific lines from annotation file by line number.

    This does random seeks to extract only the lines needed for one contig.
    More efficient than filtering entire file.

    Args:
        annotation_file: Path to annotation file
        line_numbers: Line numbers to extract

    Returns:
        List of line contents

    Note:
        For gzipped files, this still requires sequential reading
        (gzip doesn't support true random access), but we skip
        irrelevant lines quickly.
    """
    if not line_numbers:
        return []

    # Sort line numbers for sequential access
    line_numbers_set = set(line_numbers)
    min_line = min(line_numbers)
    max_line = max(line_numbers)

    extracted = []

    with smart_open(annotation_file, 'rt') as f:
        for line_num, line in enumerate(f):
            # Skip until we reach first target line
            if line_num < min_line:
                continue

            # Stop after last target line
            if line_num > max_line:
                break

            # Extract if this is a target line
            if line_num in line_numbers_set:
                extracted.append(line)

    return extracted


def get_index_cache_path(annotation_file: str) -> Path:
    """
    Get path for cached index file.

    Returns path like: annotations.gff3.gz.contig_index.pkl
    """
    annotation_path = Path(annotation_file)
    return annotation_path.parent / f"{annotation_path.name}.contig_index.pkl"


def load_or_build_index(
    annotation_file: str,
    force_rebuild: bool = False,
    messenger=None
) -> ContigIndex:
    """
    Load index from cache if available, otherwise build and cache it.

    Args:
        annotation_file: Path to annotation file
        force_rebuild: Force rebuilding even if cache exists
        messenger: Optional messenger for logging

    Returns:
        ContigIndex instance
    """
    cache_path = get_index_cache_path(annotation_file)
    annotation_path = Path(annotation_file)

    # Check if cache is valid (exists and newer than annotation file)
    if not force_rebuild and cache_path.exists():
        cache_mtime = cache_path.stat().st_mtime
        annot_mtime = annotation_path.stat().st_mtime

        if cache_mtime > annot_mtime:
            if messenger:
                messenger.info(f"Loading cached annotation index: {cache_path}")
            try:
                return ContigIndex.load(cache_path)
            except Exception as e:
                if messenger:
                    messenger.log_only(
                        f"Failed to load cached index: {e}, rebuilding",
                        level="warning"
                    )

    # Build new index
    index = build_contig_index(annotation_file, messenger)

    # Try to cache it
    try:
        index.save(cache_path)
        if messenger:
            messenger.log_only(f"Cached index saved: {cache_path}")
    except Exception as e:
        if messenger:
            messenger.log_only(
                f"Failed to cache index: {e} (continuing without cache)",
                level="warning"
            )

    return index
```

**Testing:**
```python
# Test with Chr19 data
index = build_contig_index("Homo_sapiens.Chr19.Ensembl_91.gff3.gz")
print(f"Contigs: {index.get_contigs()}")
print(f"Chr19 lines: {len(index.get_line_numbers('19'))}")
print(f"Memory: {index.memory_estimate() / 1024:.1f} KB")

# Extract chr19 lines
lines = extract_contig_lines(
    "Homo_sapiens.Chr19.Ensembl_91.gff3.gz",
    index.get_line_numbers('19')
)
print(f"Extracted {len(lines)} lines for chr19")
```

---

### Phase 2: Add Streaming Support to AnnotationHierarchyBuilder

**File:** `extraction/annotator.py`

**Changes needed:**

```python
class AnnotationHierarchyBuilder:
    """Build annotation hierarchy (genes, transcripts, exons)."""

    # Existing method - keep for backward compatibility
    def build_from_file(
        self,
        annotation_file: str,
        filter_contigs: Optional[Set[str]] = None  # NEW parameter
    ) -> Dict[str, Gene]:
        """
        Build hierarchy from annotation file.

        Args:
            annotation_file: Path to GFF3/GTF file
            filter_contigs: Optional set of contigs to include (skip others)

        Returns:
            Dictionary mapping gene IDs to Gene objects
        """
        # ... existing code ...

        for feature in self._parse_gff(annotation_file):
            # NEW: Skip if filtering and not in target set
            if filter_contigs and feature.chromosome not in filter_contigs:
                continue

            # ... rest of existing code ...

    # NEW method for streaming
    def build_from_lines(
        self,
        annotation_lines: List[str]
    ) -> Dict[str, Gene]:
        """
        Build hierarchy from pre-extracted lines (for streaming).

        This allows processing one contig at a time without loading
        the entire annotation file.

        Args:
            annotation_lines: List of GFF3/GTF lines (already filtered)

        Returns:
            Dictionary mapping gene IDs to Gene objects
        """
        genes = {}

        for line in annotation_lines:
            feature = self._parse_gff_line(line)
            if feature is None:
                continue

            self._add_feature_to_hierarchy(feature, genes)

        return genes

    def _parse_gff_line(self, line: str) -> Optional[Feature]:
        """Parse a single GFF3/GTF line."""
        # Extract existing line parsing logic
        # ... existing code ...
```

---

### Phase 3: Update extract_introns_from_annotation

**File:** `cli/main.py`

**Major refactoring to process by contig:**

```python
def extract_introns_from_annotation(
    config: IntronICConfig,
    genome_reader: GenomeReader,
    messenger: 'UnifiedMessenger',
    reporter: IntronICProgressReporter
) -> List[Intron]:
    """
    Extract introns contig-by-contig with streaming annotation processing.

    NEW APPROACH:
    1. Build lightweight annotation index (~20 MB)
    2. For each contig:
       - Extract that contig's annotation lines
       - Build hierarchy for just that contig
       - Generate introns
       - Free contig hierarchy immediately
    3. Pre-filter, extract sequences, etc. (existing logic)

    Memory savings: 12 GB peak → 2-3 GB peak (80% reduction!)
    """
    messenger.info(f"Loading annotation: {config.input.annotation}")

    # PHASE 1: Build lightweight index
    from file_io.annotation_index import load_or_build_index

    annotation_index = load_or_build_index(
        str(config.input.annotation),
        messenger=messenger
    )

    contigs = annotation_index.get_contigs()
    messenger.info(f"Found {len(contigs)} contigs in annotation")

    # PHASE 2: Process each contig incrementally
    messenger.info("Generating introns contig-by-contig (streaming mode)")

    all_introns = []
    builder = AnnotationHierarchyBuilder(
        child_features=['cds', 'exon'],
        clean_names=config.output.clean_names
    )
    generator = IntronGenerator()

    for contig_idx, contig in enumerate(contigs, 1):
        # Progress tracking
        if contig_idx % 10 == 0 or contig_idx == len(contigs):
            messenger.log_only(
                f"Processing contig {contig_idx}/{len(contigs)}: {contig}"
            )

        # 2a: Extract lines for this contig only
        from file_io.annotation_index import extract_contig_lines
        contig_lines = extract_contig_lines(
            str(config.input.annotation),
            annotation_index.get_line_numbers(contig)
        )

        if not contig_lines:
            continue

        # 2b: Build hierarchy for just this contig
        contig_genes = builder.build_from_lines(contig_lines)

        # 2c: Generate introns for this contig
        contig_introns_iter = generator.generate_from_genes(
            contig_genes,
            builder.feature_index
        )
        contig_introns = list(contig_introns_iter)

        # 2d: Accumulate and FREE contig hierarchy immediately
        all_introns.extend(contig_introns)
        del contig_lines
        del contig_genes
        del contig_introns
        gc.collect()

    messenger.info(f"Generated {len(all_introns):,} introns from {len(contigs)} contigs")

    # PHASE 3: Filter by feature type (existing logic)
    if config.extraction.feature_type == 'cds':
        introns_list = [i for i in all_introns if i.metadata.defined_by == 'cds']
    elif config.extraction.feature_type == 'exon':
        introns_list = [i for i in all_introns if i.metadata.defined_by == 'exon']
    else:
        introns_list = all_introns

    del all_introns
    gc.collect()

    # PHASE 4-6: Pre-filter and extract (existing logic continues)
    # ... rest of existing code (pre-filtering, sequence extraction, etc.) ...
```

---

## Memory Profile Comparison

### Before (Current)
```
Step 1: Parse full annotation       12 GB  ← Peak!
Step 2: Generate all introns        14 GB
Step 3: Free annotation hierarchy   2 GB
Step 4: Pre-filter                  2 GB
Step 5: Extract sequences           5 GB
───────────────────────────────────────────
Peak memory:                        14 GB
```

### After (Streaming)
```
Step 1: Build index                 20 MB
Step 2: For each contig:
  - Load contig annotations         500 MB
  - Generate introns                +100 MB
  - Free contig hierarchy           0 MB
  - Accumulate (grows gradually)    → 2 GB
Step 3: Pre-filter                  2 GB   ← Peak!
Step 4: Extract sequences           5 GB
───────────────────────────────────────────
Peak memory:                        5 GB  (down from 14 GB!)
```

**Improvement: ~9 GB saved during annotation phase!**

---

## Performance Considerations

### Index Building (Pass 1)
- **Time:** ~30-60 seconds for human genome (2.5M lines)
- **I/O:** Sequential read, very fast
- **Cache:** Can save index for reuse (instant load next time)

### Contig Processing (Pass 2)
- **Time:** Similar to current (maybe 10-20% slower)
- **Why slower:** More seeks vs single sequential read
- **Mitigation:** Index caching, efficient line extraction

### Overall Impact
- **First run:** +30-60 seconds (index build)
- **Subsequent runs:** ~Same speed (cached index)
- **Memory saved:** ~9 GB (65% reduction!)

**Trade-off is worth it!**

---

## Implementation Phases

### Phase 1: Core Indexing (2-3 hours)
- [ ] Create `file_io/annotation_index.py`
- [ ] Implement `ContigIndex` class
- [ ] Implement `build_contig_index()`
- [ ] Implement `extract_contig_lines()`
- [ ] Add caching (save/load)
- [ ] Unit tests

### Phase 2: Streaming Support (1-2 hours)
- [ ] Add `build_from_lines()` to `AnnotationHierarchyBuilder`
- [ ] Add `filter_contigs` parameter to `build_from_file()`
- [ ] Update `_parse_gff()` to support filtering
- [ ] Unit tests

### Phase 3: Integration (2-3 hours)
- [ ] Refactor `extract_introns_from_annotation()` to use streaming
- [ ] Add progress tracking for contig processing
- [ ] Update memory cleanup logic
- [ ] Integration tests

### Phase 4: Testing & Validation (1-2 hours)
- [ ] Test with Chr19 (verify same results)
- [ ] Test with full genome (measure memory)
- [ ] Verify index caching works
- [ ] Performance benchmarking

**Total estimate: 6-10 hours**

---

## Testing Strategy

### Unit Tests

```python
# test_annotation_index.py

def test_build_index():
    """Test index building."""
    index = build_contig_index("test.gff3")
    assert '19' in index.get_contigs()
    assert len(index.get_line_numbers('19')) > 0

def test_extract_lines():
    """Test line extraction."""
    index = build_contig_index("test.gff3")
    lines = extract_contig_lines("test.gff3", index.get_line_numbers('19'))
    assert all('19' in line for line in lines)

def test_index_caching():
    """Test save/load index."""
    index1 = build_contig_index("test.gff3")
    cache_file = Path("test.pkl")
    index1.save(cache_file)
    index2 = ContigIndex.load(cache_file)
    assert index1.get_contigs() == index2.get_contigs()
```

### Integration Tests

```python
def test_streaming_vs_original():
    """Verify streaming gives same results as original."""
    # Original approach (load everything)
    introns_original = extract_introns_from_annotation_original(...)

    # Streaming approach
    introns_streaming = extract_introns_from_annotation(...)

    # Should be identical
    assert len(introns_original) == len(introns_streaming)
    assert set(i.intron_id for i in introns_original) == \
           set(i.intron_id for i in introns_streaming)
```

### Memory Testing

```bash
# Measure peak memory with original
/usr/bin/time -v intronIC ... 2>&1 | grep "Maximum resident set"

# Measure peak memory with streaming
/usr/bin/time -v intronIC ... 2>&1 | grep "Maximum resident set"

# Should show ~9 GB reduction
```

---

## Backward Compatibility

✅ **Fully backward compatible:**
- No changes to CLI arguments
- No changes to output formats
- `AnnotationHierarchyBuilder.build_from_file()` still works (with optional filter)
- New methods are additions, not replacements

**Migration path:**
1. Implement new methods alongside old
2. Switch `extract_introns_from_annotation()` to use streaming
3. Keep old methods for other use cases (training, etc.)

---

## Future Enhancements (Optional)

### 1. Parallel Streaming
Process multiple contigs in parallel, each with their own index:

```python
with Pool(n_processes) as pool:
    results = pool.starmap(process_contig_streaming, inputs)
```

### 2. Smart Batching
Group small contigs together to reduce overhead:

```python
if len(contig_lines) < 1000:
    # Batch small contigs together
    batch.append(contig)
```

### 3. Index Optimization
Use more compact storage (numpy arrays, compressed):

```python
# Store as numpy array instead of list
contig_to_lines['chr1'] = np.array([120, 121, ...], dtype=np.uint32)
```

---

## Success Criteria

✅ **Memory:**
- Peak memory < 5 GB for human genome (down from 14 GB)
- Index memory < 50 MB

✅ **Performance:**
- First run: <2 minutes slower than original
- Cached runs: Same speed as original
- Index build: <1 minute

✅ **Correctness:**
- Produces identical results to original
- All existing tests pass
- No regression in accuracy

✅ **Usability:**
- Transparent to users (no CLI changes)
- Progress messages show streaming
- Cached index reused automatically

---

## Risk Assessment

### Low Risk
- ✅ Index building (simple, well-tested pattern)
- ✅ Line extraction (straightforward file I/O)
- ✅ Caching (standard pickle, optional feature)

### Medium Risk
- ⚠️ Integration with existing code (careful refactoring needed)
- ⚠️ Feature index handling (need to aggregate across contigs)
- ⚠️ Performance with many small contigs (batching may help)

### Mitigation
- Incremental implementation (phases)
- Comprehensive testing at each phase
- Keep old code path available during transition
- Monitor memory with `time -v` at each step

---

## Conclusion

**Streaming annotation processing with pre-indexing:**
- ✅ Reduces peak memory by ~9 GB (65% reduction!)
- ✅ Minimal performance impact (<2 minutes slower)
- ✅ Index caching makes subsequent runs fast
- ✅ Fully backward compatible
- ✅ Clean, maintainable implementation

**Combined with existing optimizations:**
- Pre-filtering: 85% reduction in extractions
- Annotation streaming: 65% reduction in peak memory
- Parallel processing: 3-8x speedup
- **Total: 30 GB → 5 GB (83% reduction!)**

This is the final major memory optimization needed to make intronIC run comfortably on standard workstations (8-16 GB RAM).
