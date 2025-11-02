# intronIC Refactoring - Session Notes

## Session 2025-11-02 (continued #2): Tag Format Clarification - COMPLETE ✅

### Important Tag Format Correction

**Issue:** Initial implementation used incorrect tag format for omitted introns.
- **Wrong:** `[s]`, `[a]`, `[v]` for omitted introns
- **Correct:** `[o:s]`, `[o:a]`, `[o:v]`, `[o:n]`, `[o:i]` for omitted introns

### Key Understanding (from original intronIC):
1. **Property tags and omission tags are INDEPENDENT**:
   - An intron can have both `[n]` (noncanonical property) AND `[o:s]` (omitted for being short)
   - Example: `[n][i][o:s]` = noncanonical, not longest isoform, omitted because too short

2. **Omission codes** (stored in `intron.metadata.omitted`):
   - 's' = short (below minimum length)
   - 'a' = ambiguous sequence (contains N or non-ACTG)
   - 'n' = noncanonical (non-standard splice sites, when excluded)
   - 'v' = coordinate overlap (overlapping with another intron)
   - 'i' = not in longest isoform (when excluded)

3. **Output format**: Omission codes are prefixed with `o:` to become `[o:s]`, `[o:a]`, etc.

### Changes Made:
1. **Updated TAG_TO_ATTRIBUTE dictionary** in `file_io/writers.py`:
   - Added all omission code mappings: 'o:s', 'o:a', 'o:n', 'o:v', 'o:i'
   - Maps to verbose names: omitted_short, omitted_ambiguous, etc.

2. **Verified tag generation logic**:
   - `_generate_tags()` correctly adds `[o:{code}]` format
   - Property tags (`[n]`, `[i]`, `[c]`, `[d]`) added independently
   - Both can appear together

3. **Updated tests**:
   - Added assertion for `[o:s]` tag in test_write_intron_with_tags
   - All 262 tests pass ✅

4. **Updated documentation**:
   - Fixed examples/attribute_mapping_demo.py
   - Updated docstrings to explain all omission codes

### Verification:
```bash
PYTHONPATH=. pixi run pytest  # 262 tests pass
PYTHONPATH=. pixi run python examples/attribute_mapping_demo.py  # Shows correct format
```

**Result:** Tag system now matches original intronIC exactly, with omission tags properly prefixed with `o:`.

---

## Session 2025-11-02 (continued): M1.2 COMPLETE - I/O Layer

### Completed Work

#### 1. Parser Tests (test_parsers.py - 44 tests)
- **BioGLAnnotationParser** (14 tests): GFF3/GTF parsing, phases, strands, file handling
- **BEDParser** (14 tests): BED3/BED6 formats, 0-based coordinates, extra columns
- **SequenceParser** (13 tests): .iic format, scores, empty flanks, long sequences
- **Data Structures** (3 tests): Dataclass creation and validation
- **All tests pass** with real chr19 data

#### 2. Writers Implementation (writers.py - 850 lines)
- **BEDWriter**: BED format output (.bed.iic)
  - 0-based coordinates (BED convention)
  - Intron labels with tags ([n], [i], [c], [d])
  - SVM scores, species naming
- **MetaWriter**: Comprehensive metadata (.meta.iic)
  - 14 columns: rel_score, dnts, motif, bp_context, length, parent, grandparent, index, family_size, frac_pos, phase, type_id, feature
  - Fractional position calculations
  - Null value handling
- **SequenceWriter**: Sequence output (.introns.iic)
  - Format: name, [score], upstream_flank, sequence, downstream_flank
  - Optional score column
  - Memory efficient (no buffering)
- **ScoreWriter**: Detailed scoring (.score_info.iic)
  - 14 columns: All PWM raw scores, z-scores, sequences
  - Branch point sequences (bp_seq, bp_region)
  - Decision distances
- **MappingWriter**: Mapping files
  - Duplicate mappings (.dupe_map.iic)
  - Overlap mappings (.overlap_map.iic)

**Key Features:**
- All writers support context managers
- Generator-friendly (accept iterables)
- Proper null value handling ('.')
- Format matches original intronIC exactly

#### 3. Writer Tests (test_writers.py - 700 lines, 37 tests)
- **BEDWriter** (10 tests): Basic writing, coordinates, tags, multiple introns
- **MetaWriter** (7 tests): Headers, full introns, null values, fractional positions
- **SequenceWriter** (7 tests): With/without scores, empty flanks, multiple introns
- **ScoreWriter** (6 tests): Full scores, partial scores, null values, sequences
- **MappingWriter** (5 tests): Single/multiple mappings, context managers
- **Integration** (2 tests): All writers, consistent naming

#### 4. Integration Tests (test_parser_writer_pipeline.py - 343 lines, 9 tests)
- **Annotation Parser Integration** (3 tests): Chr19 parsing, structure, parent relationships
- **Genome Reader Integration** (2 tests): Chr19 loading, subsequence extraction
- **Parser → Writer Pipeline** (3 tests): BED roundtrip, sequence roundtrip, multiple formats
- **Real Data Pipeline** (1 test): Full chr19 annotation → outputs

### Test Results
- **Total:** 262 tests (247 unit + 9 integration + 6 gold standard)
- **Status:** ✅ All passing
- **Runtime:** ~14 seconds

### Files Created/Modified
```
file_io/
├── writers.py              (850 lines) - NEW
└── __init__.py             (45 lines) - Updated to export writers

tests/unit/
└── test_writers.py         (700 lines, 37 tests) - NEW

tests/integration/
├── __init__.py             - NEW
└── test_parser_writer_pipeline.py (343 lines, 9 tests) - NEW
```

### Milestone M1.2: I/O Layer - COMPLETE ✅
**Total Lines:** 2,726 new/modified lines
- file_io/genome.py: 329 lines
- file_io/parsers.py: 459 lines
- file_io/writers.py: 850 lines
- test_parsers.py: (part of test_parsers.py)
- test_writers.py: 700 lines
- test_parser_writer_pipeline.py: 343 lines

**Capabilities:**
- ✅ Read: FASTA genomes, GFF3/GTF annotations, BED files, .iic sequences
- ✅ Write: BED, metadata, sequences, scores, mappings
- ✅ Formats match original intronIC exactly
- ✅ Memory efficient (streaming, generators)
- ✅ Fully tested with real chr19 data

---

## Session 2025-11-02 (earlier): M1.1 Complete + M1.2 Partial

### Completed Work

#### 1. Memory Optimization (M1.1 follow-up)
- **Added `slots=True`** to all 9 dataclasses (commit: `a859581`)
- **Result:** ~60-70% reduction in object overhead
- Verified: No `__dict__`, only `__slots__`
- **Critical for:** High-intron-count species (plants with 500k-1M+ introns)

#### 2. I/O Layer - Reading Components (M1.2 partial)

**Files Created:**
```
file_io/
├── __init__.py           (33 lines)
├── genome.py             (330 lines) - FASTA reading
└── parsers.py            (450 lines) - Modular parsers

tests/unit/
└── test_genome.py        (250 lines, 22 tests)
```

**file_io/genome.py:**
- `parse_fasta()`: Streaming FASTA parser with gzip support
- `GenomeReader`: Two-mode genome reader
  - **Streaming mode** (default): Memory efficient, one chromosome at a time
  - **Cached mode**: Fast random access, loads entire genome
- `extract_subsequence()`: Coordinate-based extraction with flanks
- Auto reverse complement for negative strand
- Tested with chr19 genome data ✓

**file_io/parsers.py:**
- **Modular Design** (per user request):
  - `AnnotationParser` Protocol: Abstract interface
  - `BioGLAnnotationParser`: Wraps `biogl.GxfParse`
  - **Easy to swap** implementations later
- `AnnotationLine`: Standardized parsed annotation
- `BEDParser`: BED format (0-based coordinates)
- `SequenceParser`: .iic sequence files
- Preserves exact metadata extraction from original

**tests/unit/test_genome.py:**
- 22 tests, all passing ✓
- Tests parse_fasta (5)
- Tests GenomeReader streaming (3)
- Tests GenomeReader cached (6)
- Tests subsequence extraction (5)
- Tests with real chr19 data (3)

### Key Decisions Made

1. **Renamed `io/` to `file_io/`** to avoid conflict with Python builtin `io` module

2. **Use biogl for annotations**:
   - Already a dependency
   - Permissive enough for intronIC's needs
   - Preserves exact metadata extraction
   - Wrapped in Protocol for modularity

3. **Memory strategy: Stick with in-memory approach**
   - SQLite backend discussed but deferred
   - Original approach handles even large genomes well
   - Can revisit SQLite later if needed

4. **Design principles**:
   - Modular parsers (Protocol-based, swappable)
   - Memory conscious (streaming by default)
   - Type safe (full type hints)
   - Preserves optimizations from original

### What's Next (M1.2 completion)

#### Still TODO:
1. **file_io/writers.py** (estimated 500-800 lines):
   - MetaWriter: `.meta.iic` format
   - BEDWriter: `.bed.iic` format
   - SequenceWriter: `.introns.iic` format (sequences)
   - ScoreWriter: `.score_info.iic` format
   - MappingWriter: `.dupe_map.iic`, `.overlap_map.iic`

2. **tests/unit/test_parsers.py** (estimated 30-40 tests):
   - Test BioGLAnnotationParser with real GFF3/GTF
   - Test BEDParser
   - Test SequenceParser

3. **tests/unit/test_writers.py** (estimated 25-35 tests):
   - Test all writer formats
   - Verify output matches gold standard format

4. **Integration test**:
   - Run parsers on chr19 annotation
   - Verify output structure

5. **Update MILESTONES.md**: Mark M1.2 complete

### Output Format Reference

From gold standard files:

**`.meta.iic`** (Tab-delimited metadata):
```
name    rel_score    dnts    motif_schematic    bp_context    length    parent    grandparent    index    family_size    frac_pos    phase    type_id    feature
```

**`.bed.iic`** (BED format, 0-based):
```
chrom    start    stop    name    score    strand
```

**`.introns.iic`** (Sequences with flanks):
```
name    score    upstream_flank    sequence    downstream_flank
```

**`.score_info.iic`** (Detailed scores):
```
name    rel_score    svm_score    decision_dist    5'_seq    5'_raw    5'_z    bp_seq_u12    bp_seq_u2    bp_raw    bp_z    3'_seq    3'_raw    3'_z
```

**`.dupe_map.iic`** (Duplicate mappings):
```
duplicate_name    representative_name
```

### Git Status
- All work committed (commit: `92f2d3d`)
- Branch: `refactor/phase-0-foundation`
- Tests: 166 passing (144 from Phase 0/M1.1 + 22 from M1.2)

### Performance Notes
- GenomeReader tested with chr19 (~58MB)
- Loads in ~1.4 seconds cached
- Memory efficient streaming works well
- No performance regressions from original

### Important Reminders for Next Session
1. Check original for any performance optimizations before implementing writers
2. Writers need to be generator-friendly (don't hold all data in memory)
3. Output formats must match gold standard exactly
4. Test with chr19 data for integration verification
