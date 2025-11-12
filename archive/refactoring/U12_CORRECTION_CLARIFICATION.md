# U12 Boundary Correction - Sequence Re-extraction Clarified

## User Question (2025-11-05)

> "For U12 boundary correction, we don't actually need to 're-extract' the sequence, right? 
> We should have flanking sequence we extract initially around each intron anyway, so we can 
> use that to determine the 'corrected' sequence."

## Answer: We ARE Re-extracting, But From In-Memory Cache (Efficient!)

### Original Implementation

```python
# intronIC.py:2669
for region_name, region_seq in fasta_parse(flex_open(genome)):
    # region_seq = ENTIRE CHROMOSOME loaded into memory
    
    intron = assign_seqs(intron, region_seq, ...)  # Extract from region_seq
    
    if u12_correction(intron):  # Changes intron.start and intron.stop
        intron = assign_seqs(intron, region_seq, ...)  # Re-extract from SAME region_seq
```

The key insight: `region_seq` is the **entire chromosome sequence in memory**. The second 
`assign_seqs` call doesn't go back to disk - it just re-slices the in-memory string with 
the new coordinates.

### Our Implementation

```python
# cli/main.py:261
sequence_extractor = SequenceExtractor(genome_file, use_cache=True)

# First extraction - loads chromosome into GenomeReader cache
introns_with_seq = sequence_extractor.extract_sequences(introns_list, ...)

# U12 correction changes coordinates
corrected_intron, was_corrected = correct_intron_if_needed(intron, ...)

# Re-extraction uses CACHED chromosome sequence
corrected_with_seq = sequence_extractor.extract_sequences([corrected_intron], ...)
```

**GenomeReader caching** (file_io/genome.py:82-112):
- `GenomeReader(cached=True)` loads chromosomes into memory
- First `get_sequence("chr1")` - loads from disk
- Subsequent `get_sequence("chr1")` - returns cached string
- Re-extraction is just string slicing: `region_seq[start:stop]`

### Why Not Use Existing Flanks?

**Answer:** The flanks might not be big enough!

```python
# Example: Correction shifts boundary by -3bp
Original:  [--flank(50bp)--][GCTATCCTT...intron...][--flank(50bp)--]
                            ^ Old start

Corrected: [--flank(50bp)--][ATATCCTT...intron...][--flank(50bp)--]
                         ^ New start (3bp earlier)
```

If the original flanks were extracted with `flank_len=50`, and we shift -3bp:
- We need 3 more bases on the upstream side
- Those 3 bases are NOT in `intron.sequences.upstream_flank`
- We need to go back to `region_seq` to get them

**BUT:** Since `region_seq` is cached in memory, this is just string indexing, not disk I/O!

### Performance

**Original:** 
- Entire chr loaded once: O(chr_length) disk read
- Re-extraction: O(1) - just string slicing from memory

**Ours:**
- Entire chr cached once: O(chr_length) disk read  
- Re-extraction: O(1) - GenomeReader.get_sequence() returns cached chr
- Subsequence extraction: O(1) - string slicing

**Conclusion:** Both implementations are equally efficient! The "re-extraction" is a 
misnomer - it's really just "re-slicing from cached chromosome".

### Alternative Optimization (Not Needed)

We COULD adjust existing sequences if `abs(shift) < min(upstream_flank_len, downstream_flank_len)`:

```python
if abs(shift) <= len(intron.sequences.upstream_flank):
    # Can adjust in-place without accessing genome
    if shift < 0:  # Move upstream
        new_seq = intron.sequences.upstream_flank[shift:] + intron.sequences.seq[:shift]
    else:  # Move downstream
        new_seq = intron.sequences.seq[shift:] + intron.sequences.downstream_flank[:shift]
```

**But this adds complexity for minimal benefit** since GenomeReader caching already makes 
re-extraction O(1).

### Recommendation

✅ **Keep current implementation** - it's clean, correct, and already efficient.

The documentation in U12_CORRECTION_OPTIMIZATION.md suggested a premature optimization. 
The current approach matches the original and performs well.
