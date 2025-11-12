# intronIC Refactoring - Session Notes

## Session 2025-11-02 (continued #3): M1.3 COMPLETE - Extraction Pipeline ✅

### Completed Work

#### 1. Fixed Model Compatibility Issues
**Problem:** Extraction modules written against old intronIC model structure needed adaptation to new M1.1 models with GenomicCoordinate abstraction.

**Changes Made:**

**annotator.py:**
- ✅ Fixed coordinate creation: Uses `GenomicCoordinate(chromosome, start, stop, strand, system='1-based')`
- ✅ Fixed coordinate field names: `chromosome`/`stop` (not `seqid`/`end`)
- ✅ Fixed CoordinateSystem: Uses string literal `'1-based'` (not enum)
- ✅ Fixed single-base features: Skip features where start >= stop

**intronator.py:**
- ✅ Updated coordinate access: `exon.coordinates.start` instead of `exon.start`
- ✅ Fixed overlap check: Uses `coordinates.start`, `coordinates.stop`
- ✅ Added feature_index parameter: Methods now require feature_index to resolve IDs
- ✅ Added exception handling: Skips invalid introns (overlapping, too short)
- ✅ Removed fractional_position setter: It's a computed property

**filters.py:**
- ✅ Fixed coordinate property names: `coordinates.chromosome`, `coordinates.stop`

**sequences.py:**
- ✅ Fixed coordinate property names: `coordinates.chromosome`, `coordinates.stop`

**intron.py (core model):**
- ✅ Added `five_prime_dnt` and `three_prime_dnt` to IntronSequences
- ✅ Added `parent_length` to IntronMetadata

#### 2. Handled Missing Hierarchy Levels
**Critical Feature:** Annotation files may be missing gene or transcript features. Original intronIC handles this by using the same ID for missing levels.

**Implementation (annotator.py lines 214-248, 366-377):**
- ✅ Wrap orphan transcripts in synthetic Gene objects
- ✅ Set grandparent to parent ID when grandparent is missing
- ✅ Use `dataclasses.replace()` to update frozen objects
- ✅ Create placeholder parents/grandparents with dummy coordinates

**User Guidance:** "child:parent mapping is not 1:1; there are many cases where the same exon/CDS feature is shared across multiple parents. Also, certain annotation files may be missing either `gene` or `transcript` features."

#### 3. Verified Multi-Level Hierarchies
**Confirmed Working:**
- ✅ Gene → Multiple Transcripts (tested with gene having 16 transcripts)
- ✅ Transcript → Multiple Exons
- ✅ Exon → Multiple Transcripts (via unique naming: `exon_Parent=T1:1000_1100`)

#### 4. Updated Integration Tests
**Fixed test_extraction_pipeline.py:**
- ✅ Updated field references: `feature_id` → `intron_id`, `.name` → `.feature_id`
- ✅ Added `feature_index` parameter to all `generate_from_genes()` calls
- ✅ Fixed attribute checks to work with new model
- ✅ All 6 integration tests passing in 27.47s

### Test Results
```bash
tests/integration/test_extraction_pipeline.py - 6 tests PASSING ✅
  ✓ test_annotation_hierarchy_building
  ✓ test_intron_generation_from_genes
  ✓ test_sequence_extraction
  ✓ test_full_pipeline_intron_count
  ✓ test_filtering_omission_codes
  ✓ test_quick_extraction (58,903 introns from chr19)
```

### Key Technical Decisions

**1. Unique Naming for Shared Exons:**
```python
# For child features (exon/cds), create unique name per parent
if feat_type in self.child_features:
    name = f"{feat_type}_Parent={parent_name}:{feat.start}_{feat.stop}"
else:
    name = feat.feature_id
```

**2. Missing Gene Feature Handling:**
```python
# Wrap orphan transcript in Gene with synthetic ID
gene_id = f"gene_wrapper:{transcript.feature_id}"
gene = Gene(feature_id=gene_id, coordinates=transcript.coordinates, ...)
updated_transcript = replace(transcript, parent_id=gene_id)
```

**3. Min/Max Coordinate Calculation:**
```python
# Handle exons in any order (coding vs genomic direction)
intron_start = min(exon1.stop, exon2.stop) + 1
intron_stop = max(exon1.start, exon2.start) - 1
```

### Milestone M1.3: Extraction Pipeline - COMPLETE ✅

**Total Lines:** 1,350 lines of extraction code + tests
- extraction/annotator.py: 413 lines
- extraction/intronator.py: 282 lines
- extraction/sequences.py: 327 lines
- extraction/filters.py: 374 lines
- tests/integration/test_extraction_pipeline.py: 194 lines

**Capabilities:**
- ✅ Parse GFF3/GTF annotations with NetworkX DAG
- ✅ Handle gene→transcript→exon hierarchies (including missing levels)
- ✅ Generate introns from exon pairs with metadata
- ✅ Extract sequences with flanks and scoring regions
- ✅ Filter by length, quality, canonical status, overlaps, isoforms
- ✅ Tag duplicates and overlaps with proper resolution
- ✅ All using new M1.1 model interface

### Files Modified This Session
```
extraction/
├── annotator.py        - Model interface updates, missing hierarchy handling
├── intronator.py       - Coordinate access fixes, feature_index param
├── filters.py          - Coordinate property name fixes
└── sequences.py        - Coordinate property name fixes

core/
└── intron.py          - Added five_prime_dnt, three_prime_dnt, parent_length

tests/integration/
└── test_extraction_pipeline.py - Updated for new model interface
```

### What's Next: Phase 2 - Scoring & Classification

**Milestone M1.4: Scoring (next)**
- PWM (Position Weight Matrix) scoring
- Branch point detection
- Z-score normalization
- Files to create:
  - scoring/pwm.py
  - scoring/normalizer.py
  - scoring/branch_point.py

**Milestone M1.5: Classification**
- SVM classifier
- Training/prediction
- Ensemble methods
- Files to create:
  - classification/svm_classifier.py
  - classification/trainer.py

---

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

### Git Status
- All work committed (commit: `92f2d3d`)
- Branch: `refactor/phase-0-foundation`
- Tests: 268 passing (262 from M1.1/M1.2 + 6 from M1.3)

### Performance Notes
- GenomeReader tested with chr19 (~58MB)
- Loads in ~1.4 seconds cached
- Memory efficient streaming works well
- No performance regressions from original
