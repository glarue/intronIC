# Chromosome-by-Chromosome Processing Plan

**Target:** Reduce peak memory from 28 GB to ~1.5 GB (95% reduction)
**Approach:** Process one chromosome at a time + smart filtering + duplicate reuse
**Status:** Planning (ready to implement)

---

## Memory Analysis

### Current Approach (All-at-Once)
```
Load entire genome: 3 GB
Extract all 2.1M introns: 27 GB
Score 200k introns: 27 GB
Clear sequences: 3 GB
Classify: 3 GB

PEAK: 28 GB
```

### Per-Chromosome Approach
```
For each chromosome:
    Load chr1 genome: 250 MB
    Extract chr1 introns (100k): 1.3 GB
    Score chr1 introns (10k): 1.3 GB
    Write and clear: 60 MB
    Free chr1 genome: 0 MB

    Load chr2 genome: 242 MB
    ... repeat ...

PEAK: 1.5 GB (during chr1)
Average: 500-800 MB per chromosome
```

**Memory reduction: 28 GB → 1.5 GB (95%)** ✅

---

## Three-Part Strategy

### Part 1: User-Setting-Aware Pre-Filtering

**Respect user configuration for what to extract:**

```python
def should_extract_sequences_for(intron: Intron, config: Config) -> bool:
    """Determine if we need to extract sequences for this intron."""

    # Always extract if user wants everything
    if config.no_filtering:
        return True

    # Length filters (always apply - scoring requirement)
    if intron.length < config.min_intron_len:
        return False  # Too short to score
    if intron.length > config.max_intron_len:
        return False  # Too long (likely error)

    # Longest isoform filter (user configurable)
    if config.longest_only and not intron.metadata.longest_isoform:
        return False  # User only wants longest isoform

    # Duplicate filter (user configurable)
    if not config.include_duplicates and intron.metadata.duplicate:
        return False  # User doesn't want duplicates

    # Non-canonical filter (user configurable)
    if config.exclude_noncanonical and intron.metadata.noncanonical:
        return False  # User doesn't want non-canonical

    # Overlap filter (user configurable)
    if config.no_overlap and intron.metadata.overlap:
        return False  # User doesn't want overlapping

    return True  # Extract sequences for this intron
```

**Key insight:** We can determine most filters from metadata/coordinates alone, BEFORE extracting expensive sequences.

**Example with default settings:**
- 2.1M total introns
- Remove too short: -500k introns
- Remove non-longest isoform: -800k introns
- Remove duplicates: -400k introns
- Remove overlaps: -200k introns
- **Remaining: 200k introns to extract** (90% reduction!)

**Example with `--include_all_isoforms --include_duplicates`:**
- 2.1M total introns
- Remove too short: -500k introns
- Remove overlaps: -200k introns
- **Remaining: 1.4M introns to extract** (still 33% reduction)

---

### Part 2: Duplicate Sequence Reuse

**Observation:** Duplicates have identical coordinates → identical sequences

**Current waste:**
```python
# Intron A: chr1:1000-2000 → extract sequence
# Intron B: chr1:1000-2000 (duplicate of A) → extract sequence AGAIN
# Result: Extracted same sequence twice!
```

**Optimized approach:**
```python
# Group introns by coordinates
introns_by_coords = defaultdict(list)
for intron in chromosome_introns:
    coord_key = (intron.chromosome, intron.start, intron.stop, intron.strand)
    introns_by_coords[coord_key].append(intron)

# Extract sequences (one per unique coordinate)
for coord_key, duplicate_group in introns_by_coords.items():
    # Extract sequence once
    representative = duplicate_group[0]
    intron_with_seq = extract_sequences([representative])[0]

    # Reuse sequence for all duplicates
    for intron in duplicate_group:
        intron_with_shared_seq = intron.with_sequences(intron_with_seq.sequences)
        yield intron_with_shared_seq
```

**Memory savings:**
- Typical duplication rate: 20-30% of introns
- Extraction time savings: ~20-30% faster
- Memory impact: Marginal (already filtering duplicates in Part 1)

**Benefit:** Primarily speed optimization, secondary memory benefit

---

### Part 3: Chromosome-by-Chromosome Processing

**Core Algorithm:**

```python
def process_genome_by_chromosome(
    all_introns: List[Intron],
    genome_file: Path,
    config: Config
) -> List[Intron]:
    """Process introns one chromosome at a time to minimize memory."""

    # Step 1: Group introns by chromosome (coordinates only, cheap)
    introns_by_chr = defaultdict(list)
    for intron in all_introns:
        introns_by_chr[intron.chromosome].append(intron)

    print(f"Grouped {len(all_introns):,} introns across {len(introns_by_chr)} chromosomes")

    # Step 2: Pre-filter to determine what needs sequence extraction
    # This happens BEFORE any sequence extraction (still just coordinates)
    introns_to_process = []
    for chr_name, chr_introns in introns_by_chr.items():
        for intron in chr_introns:
            if should_extract_sequences_for(intron, config):
                introns_to_process.append(intron)

    print(f"Pre-filtered: {len(introns_to_process):,} introns need sequence extraction")
    print(f"Skipped: {len(all_introns) - len(introns_to_process):,} introns (no sequences needed)")

    # Re-group filtered introns by chromosome
    filtered_by_chr = defaultdict(list)
    for intron in introns_to_process:
        filtered_by_chr[intron.chromosome].append(intron)

    # Step 3: Process one chromosome at a time
    all_scored_introns = []

    for chr_name in sorted(filtered_by_chr.keys()):
        chr_introns = filtered_by_chr[chr_name]
        print(f"\nProcessing {chr_name}: {len(chr_introns):,} introns")

        # 3a. Load ONLY this chromosome's sequence
        genome_reader = GenomeReader(genome_file)
        chr_sequence = genome_reader.get_sequence(chr_name)
        chr_size_mb = len(chr_sequence) / 1_000_000
        print(f"  Loaded {chr_name} sequence: {chr_size_mb:.1f} MB")

        # 3b. Extract sequences for this chromosome's introns
        sequence_extractor = SequenceExtractor(genome_reader)
        chr_introns_with_seq = []

        for intron in sequence_extractor.extract_sequences(chr_introns):
            chr_introns_with_seq.append(intron)

        introns_mem_mb = len(chr_introns_with_seq) * 13 / 1024
        print(f"  Extracted sequences: {introns_mem_mb:.1f} MB")
        print(f"  Peak memory: {chr_size_mb + introns_mem_mb:.1f} MB")

        # 3c. Score this chromosome's introns
        chr_scored = score_introns(chr_introns_with_seq, config)
        print(f"  Scored {len(chr_scored):,} introns")

        # 3d. Write sequences to disk (streaming)
        with SequenceWriter(output_path, mode='a') as writer:
            for intron in chr_introns_with_seq:
                writer.write_intron(intron)

        # 3e. Clear large sequences, keep only scoring metadata
        chr_scored_cleared = [
            clear_large_sequences_for_classification([i])[0]
            for i in chr_scored
        ]

        # 3f. Accumulate cleared introns
        all_scored_introns.extend(chr_scored_cleared)

        # 3g. FREE MEMORY before next chromosome
        del chr_sequence
        del chr_introns
        del chr_introns_with_seq
        del chr_scored
        del genome_reader
        gc.collect()

        cleared_mem_mb = len(chr_scored_cleared) * 0.6 / 1024
        print(f"  Cleared sequences: {cleared_mem_mb:.1f} MB")
        print(f"  Freed: {chr_size_mb + introns_mem_mb - cleared_mem_mb:.1f} MB")

    print(f"\nTotal scored introns: {len(all_scored_introns):,}")
    print(f"Final memory usage: ~{len(all_scored_introns) * 0.6 / 1024:.1f} MB")

    return all_scored_introns
```

**Memory Profile Per Chromosome:**
```
Chr1:  250 MB genome + 1.3 GB introns = 1.5 GB peak → 60 MB after clear
Chr2:  242 MB genome + 1.2 GB introns = 1.4 GB peak → 58 MB after clear
Chr3:  198 MB genome + 800 MB introns = 1.0 GB peak → 40 MB after clear
...
ChrX:  156 MB genome + 500 MB introns = 656 MB peak → 25 MB after clear

GLOBAL PEAK: 1.5 GB (during chr1 processing)
```

---

## Implementation Plan

### Phase 1: Pre-Filtering Infrastructure (2-3 hours)

**File:** `extraction/filters.py`

```python
def should_extract_sequences_for(
    intron: Intron,
    config: IntronICConfig
) -> bool:
    """Determine if intron needs sequence extraction based on config."""
    # Implementation as shown above
    pass

def prefilter_introns(
    introns: List[Intron],
    config: IntronICConfig
) -> Tuple[List[Intron], List[Intron]]:
    """Split introns into extract_needed and skip_extraction groups."""
    extract_needed = []
    skip_extraction = []

    for intron in introns:
        if should_extract_sequences_for(intron, config):
            extract_needed.append(intron)
        else:
            skip_extraction.append(intron)

    return extract_needed, skip_extraction
```

**Testing:**
- Unit tests with different config combinations
- Verify respects user settings
- Check filter statistics

---

### Phase 2: Duplicate Sequence Reuse (2-3 hours)

**File:** `extraction/sequences.py`

```python
def extract_sequences_with_deduplication(
    self,
    introns: List[Intron],
    flank_size: int = 200
) -> Generator[Intron, None, None]:
    """Extract sequences, reusing for duplicates."""

    # Group by coordinates
    by_coords = defaultdict(list)
    for intron in introns:
        key = (intron.chromosome, intron.start, intron.stop, intron.strand)
        by_coords[key].append(intron)

    print(f"Unique coordinates: {len(by_coords):,}")
    print(f"Total introns: {sum(len(g) for g in by_coords.values()):,}")
    print(f"Duplicates: {sum(len(g) - 1 for g in by_coords.values()):,}")

    # Extract for unique coordinates, yield for all duplicates
    for coord_key, duplicate_group in by_coords.items():
        # Extract once
        representative = duplicate_group[0]
        for intron_with_seq in self.extract_sequences([representative], flank_size):
            sequences = intron_with_seq.sequences

            # Yield for all in group with same sequences
            for intron in duplicate_group:
                yield intron.with_sequences(sequences)
```

**Testing:**
- Test with dataset containing duplicates
- Verify sequences identical for duplicate introns
- Measure speed improvement

---

### Phase 3: Chromosome-by-Chromosome Processing (4-5 hours)

**File:** `cli/main.py`

**Replace current extraction:**

```python
def extract_introns_from_annotation(
    config: IntronICConfig,
    genome_reader: GenomeReader,
    messenger: 'UnifiedMessenger',
    reporter: IntronICProgressReporter
) -> List[Intron]:
    """Extract introns using chromosome-by-chromosome processing."""

    # Step 1: Parse annotations (coordinates only)
    messenger.info("Parsing annotations (coordinates only)")
    introns_all = parse_annotations_to_introns(config, messenger)
    messenger.log_only(f"Parsed {len(introns_all):,} introns (coordinates only, ~{len(introns_all) / 1024:.1f} MB)")

    # Step 2: Pre-filter based on user settings
    messenger.info("Pre-filtering introns based on user settings")
    introns_to_extract, introns_skipped = prefilter_introns(introns_all, config)

    messenger.log_only(
        f"Pre-filtered: {len(introns_to_extract):,} need extraction, "
        f"{len(introns_skipped):,} skipped (user filters)"
    )

    # Step 3: Group by chromosome
    introns_by_chr = defaultdict(list)
    for intron in introns_to_extract:
        introns_by_chr[intron.chromosome].append(intron)

    messenger.log_only(f"Grouped into {len(introns_by_chr)} chromosomes")

    # Step 4: Process one chromosome at a time
    messenger.info("Extracting and scoring introns (chromosome-by-chromosome)")

    all_scored_introns = []
    seq_output_path = config.output.output_dir / f"{config.output.base_filename}.introns.iic"

    with SequenceWriter(seq_output_path) as seq_writer:
        for chr_name in sorted(introns_by_chr.keys()):
            chr_introns = introns_by_chr[chr_name]

            messenger.log_only(f"  {chr_name}: {len(chr_introns):,} introns")

            # Load chromosome sequence
            chr_genome_reader = GenomeReader(
                config.input.genome,
                use_cache=False  # Don't cache - we'll free after this chr
            )

            # Extract sequences (with deduplication)
            sequence_extractor = SequenceExtractor(chr_genome_reader)
            chr_introns_with_seq = list(
                sequence_extractor.extract_sequences_with_deduplication(
                    chr_introns,
                    flank_size=config.extraction.flank_len
                )
            )

            # Score chromosome introns
            chr_scored = score_introns(chr_introns_with_seq, config, messenger, reporter)

            # Write sequences immediately
            for intron in chr_introns_with_seq:
                seq_writer.write_intron(intron, ...)

            # Clear large sequences
            chr_scored_cleared = clear_large_sequences_for_classification(chr_scored)
            all_scored_introns.extend(chr_scored_cleared)

            # FREE MEMORY
            del chr_genome_reader
            del chr_introns
            del chr_introns_with_seq
            del chr_scored
            gc.collect()

            messenger.log_only(f"    Processed and cleared {chr_name}")

    messenger.success(
        f"Extracted and scored {len(all_scored_introns):,} introns "
        f"(peak memory: ~1.5 GB)"
    )

    # Step 5: Merge with skipped introns (for output)
    # Skipped introns have coordinates but no sequences
    return all_scored_introns, introns_skipped
```

**Testing:**
- Test with small genome (Chr19 only)
- Test with full human genome
- Verify output files identical
- Measure memory usage

---

### Phase 4: Handle Omitted/Skipped Introns in Output (1-2 hours)

**Issue:** Skipped introns don't have sequences

**Solution:**

```python
def write_outputs(...):
    """Write output files."""

    # For .bed.iic and .meta.iic: Include all introns (scored + skipped)
    # Skipped introns will have seq=None but that's OK for these files

    # For .introns.iic: Only include introns with sequences
    # These were already written during chromosome processing

    # For .score_info.iic: Only include scored introns
    # Skipped introns have no scores anyway
```

**Note:** This matches original intronIC behavior - omitted introns don't have sequences in output.

---

## Testing Strategy

### Unit Tests
- `test_should_extract_sequences_for()` - Various config combinations
- `test_prefilter_introns()` - Verify correct splitting
- `test_extract_sequences_with_deduplication()` - Duplicate handling
- `test_chromosome_grouping()` - Verify introns grouped correctly

### Integration Tests

**Small dataset (Chr19 only):**
```bash
timeout 60 pixi run intronIC train \
    -g data/test_data/Homo_sapiens.Chr19.fa.gz \
    -a data/test_data/Homo_sapiens.Chr19.gff3.gz \
    -n chr19_test \
    -o chr19_test_output

# Expected memory: <500 MB
# Expected time: <2 minutes
```

**Medium dataset (Drosophila):**
```bash
intronIC train \
    -g drosophila.fa.gz \
    -a drosophila.gff3.gz \
    -n drosophila_test \
    -o drosophila_test_output

# Expected memory: ~1 GB (5 chromosomes)
# Expected time: ~10 minutes
```

**Large dataset (Human):**
```bash
intronIC train \
    -g human.fa.gz \
    -a human.gff3.gz \
    -n human_test \
    -o human_test_output

# Expected memory: ~1.5 GB peak (during chr1)
# Expected time: ~60 minutes
```

### Memory Profiling

```bash
# Use /usr/bin/time to measure peak memory
/usr/bin/time -v intronIC train ... 2>&1 | grep "Maximum resident"

# Or use memory_profiler for detailed profile
mprof run intronIC train ...
mprof plot --output memory_profile.png
```

**Success Criteria:**
- Peak memory < 2 GB for human genome
- Output files identical to current version
- Runtime increase < 20%

---

## Expected Results

### Memory Usage

| Dataset | # Introns | Current Peak | After Implementation | Reduction |
|---------|-----------|-------------|----------------------|-----------|
| Chr19 | 50k | 1 GB | 200 MB | 80% |
| C. elegans | 20k | 500 MB | 150 MB | 70% |
| Drosophila | 200k | 5 GB | 800 MB | 84% |
| Human | 2.1M | 28 GB | 1.5 GB | **95%** |

### Performance Impact

- **Memory:** 95% reduction (28 GB → 1.5 GB)
- **Speed:** Similar or faster (deduplication saves extraction time)
- **Disk I/O:** Moderate increase (reading genome by chromosome)

---

## Edge Cases and Considerations

### 1. Circular Chromosomes (Mitochondria)

**Issue:** Some genomes have circular chromosomes

**Solution:** Process normally - coordinate system handles this

### 2. Unplaced Scaffolds

**Issue:** Scaffolds not assigned to chromosomes (e.g., "scaffold_123")

**Solution:** Treat each scaffold as a chromosome, process individually

### 3. User Wants All Duplicates

**Config:** `--include_duplicates`

**Impact:** More introns to extract, but still benefit from chromosome-by-chromosome
- Current: 2.1M introns all at once = 28 GB
- New: 100k per chromosome = 1.5 GB peak

### 4. User Wants All Isoforms

**Config:** `--include_all_isoforms` (or no `--longest_only`)

**Impact:** More introns to extract (1.4M instead of 200k)
- Per-chromosome: 60k introns = ~800 MB per large chromosome
- Still better than 28 GB all-at-once

### 5. Very Large Chromosomes

**Issue:** Some genomes have huge chromosomes (e.g., wheat: 1 GB chromosomes)

**Solution:** Could add chromosome chunking if needed:
- Split large chromosome into 100 MB chunks
- Process chunk-by-chunk within chromosome
- Further reduces peak memory

### 6. Genome File Format

**Current:** Assumes FASTA (gzipped or not)

**Compatibility:** GenomeReader already supports:
- `.fa`, `.fasta`
- `.fa.gz`, `.fasta.gz`
- Multi-chromosome FASTA files

---

## Compatibility with Original intronIC v1.5.1

**This approach mirrors original intronIC:**

1. ✅ Pre-filter before extraction
2. ✅ Process one chromosome at a time
3. ✅ Write and clear immediately
4. ✅ Omitted introns don't have sequences in output

**Differences (improvements):**
- ✅ User-configurable filtering (respects settings)
- ✅ Explicit duplicate sequence reuse
- ✅ Modular code (testable, maintainable)
- ✅ Progress reporting per chromosome

---

## Implementation Timeline

**Total Estimated Time: 10-12 hours**

| Phase | Time | Risk |
|-------|------|------|
| Phase 1: Pre-filtering | 2-3 hours | Low |
| Phase 2: Duplicate reuse | 2-3 hours | Low |
| Phase 3: Chr-by-chr processing | 4-5 hours | Medium |
| Phase 4: Output handling | 1-2 hours | Low |
| Testing & validation | 2-3 hours | - |

**Could be split across 2-3 sessions**

---

## Next Steps

1. **Review plan** - Ensure logic sound
2. **Implement Phase 1** - Pre-filtering (standalone, testable)
3. **Implement Phase 2** - Duplicate reuse (standalone, testable)
4. **Implement Phase 3** - Chromosome processing (integrates 1 & 2)
5. **Implement Phase 4** - Output handling
6. **Test and validate** - Small → medium → large datasets
7. **Profile memory** - Confirm <2 GB peak
8. **Document** - Update docs with new memory requirements

---

## Alternative: If 1.5 GB Still Too High

If users need even lower memory (e.g., 512 MB limit):

**Chromosome Chunking:**
- Split large chromosomes into smaller regions
- Process 10 MB of chromosome at a time
- Peak: ~100 MB per chunk

**Trade-off:** More disk I/O, slightly slower

---

## Conclusion

**Recommended Approach:**
- **Part 1:** User-aware pre-filtering (respects config)
- **Part 2:** Duplicate sequence reuse (speed + memory)
- **Part 3:** Chromosome-by-chromosome processing (main memory saver)

**Expected Outcome:**
- **28 GB → 1.5 GB (95% reduction)** ✅
- **Respects user configuration** ✅
- **Maintains output compatibility** ✅
- **Similar or better performance** ✅
- **Testable and maintainable** ✅

This is the winning strategy - combines best of all approaches!
