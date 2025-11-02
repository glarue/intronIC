# intronIC - Comprehensive Code Documentation

**Last Updated:** 2025-11-01
**Version:** 1.5.1
**Author:** Graham E. Larue
**Repository:** https://github.com/glarue/intronIC

---

## Table of Contents

1. [Project Overview](#project-overview)
2. [Architecture & Design](#architecture--design)
3. [Core Components](#core-components)
4. [Data Flow & Pipeline](#data-flow--pipeline)
5. [Key Algorithms](#key-algorithms)
6. [File Formats](#file-formats)
7. [Code Organization](#code-organization)
8. [Usage & Configuration](#usage--configuration)
9. [Performance Considerations](#performance-considerations)
10. [Future Improvement Areas](#future-improvement-areas)

---

## Project Overview

### Purpose
intronIC is a bioinformatics tool for **extracting and classifying intron sequences** from genomic data. It performs binary classification of introns as:
- **U2-type (major)**: ~99.5% of introns in most genomes
- **U12-type (minor)**: ~0.5% of introns, rare but functionally important

### Classification Method
Uses a **Support Vector Machine (SVM)** classifier trained on curated reference datasets:
- Training data: 100,000+ U2 sequences, 500+ U12 sequences
- Features: Position-weight matrix (PWM) scores for splice site regions
- Scoring regions: 5' splice site, branch point, 3' splice site

### Key Stats
- **Main module**: 6,093 lines (intronIC.py)
- **Classes**: 6 (GenomeFeature hierarchy)
- **Core functions**: 60+ major functions
- **Performance**: Human genome ~20-35 minutes (single process), ~5GB memory
- **License**: GPL v3.0

---

## Architecture & Design

### Object-Oriented Data Model

```
GenomeFeature (base class)
├── Parent (abstract parent for hierarchical features)
│   ├── Gene
│   └── Transcript
├── Exon
└── Intron (CORE CLASS - 41 optimized attributes)
```

### The Intron Class (Primary Data Structure)

The `Intron` class is the heart of the system, with **41 `__slots__` attributes** for memory optimization:

**Location & Identity:**
- `region`, `start`, `stop`, `strand` - genomic coordinates
- `parent`, `grandparent` - transcript and gene IDs
- `index` - ordinal position in transcript
- `family_size` - total introns in transcript

**Sequence Data:**
- `seq` - full intron sequence
- `five_seq`, `three_seq` - splice site sequences
- `bp_seq`, `bp_region_seq` - branch point sequences
- `upstream_flank`, `downstream_flank` - exonic context

**Scoring & Classification:**
- `five_raw_score`, `bp_raw_score`, `three_raw_score` - PWM scores
- `five_z_score`, `bp_z_score`, `three_z_score` - normalized scores
- `svm_score` - U12 probability (0-100%)
- `relative_score` - score relative to threshold
- `type_id` - 'u2' or 'u12'

**Quality & Metadata:**
- `dnts` - terminal dinucleotides (e.g., 'GT-AG')
- `noncanonical` - boolean for non-standard boundaries
- `omitted` - omission reason code ('s'=short, 'a'=ambiguous, etc.)
- `duplicate` - reference to duplicate intron
- `overlap` - overlapping coordinate record
- `longest_isoform` - boolean for longest transcript
- `corrected` - boolean for boundary-adjusted introns

### Key Design Patterns

**Factory Pattern:**
- `Intron.from_exon_pair()` - creates introns from adjacent exons

**Generator Pattern:**
- `get_sub_seqs()` yields introns one-at-a-time for memory efficiency
- Sequences cleared after output to reduce memory footprint

**Indexing Strategies:**
- Coordinate-based hashing for O(1) duplicate detection
- Region-based indexing for fast chromosome lookup
- Flat annotation dictionary for quick parent/grandparent access

---

## Core Components

### 1. Intron Extraction System

**Three Input Modes:**

#### Mode 1: Annotation + Genome
**Function:** `introns_from_annotation()` (line 4102)
**Process:**
1. Parse GFF3/GTF into hierarchical feature objects
2. Build parent-child relationships using NetworkX DAG
3. Extract introns from CDS/exon features
4. Create Intron objects from adjacent exon pairs via `intronator()`

#### Mode 2: BED File + Genome
**Function:** `introns_from_bedfile()` (line 4161)
**Process:**
1. Parse BED lines (0-based, half-open coordinates)
2. Automatically adjust to 1-based coordinates
3. Create Intron objects directly from coordinates
4. Fetch sequences from genome

#### Mode 3: Pre-extracted Sequences
**Function:** `introns_from_seqs()` (line 4077)
**Process:**
1. Load from `.iic` format file
2. Filter omitted introns
3. Skip genome fetching (already have sequences)

### 2. Filtering & Quality Control

**Function:** `filter_introns_write_files()` (line 4668)

**Omission Criteria:**
- **'s' (short)**: Length < minimum (default 30bp)
- **'a' (ambiguous)**: Contains non-ACTG characters
- **'n' (noncanonical)**: Non-standard terminal dinucleotides
- **'i' (isoform)**: Not from longest transcript
- **'v' (overlap)**: Coordinates overlap other introns

**Duplicate Detection:**
- Hash coordinates (region, strand, start, stop)
- Priority: longest isoform > transcript length > annotation order
- First occurrence becomes representative

**Non-Canonical Boundary Adjustment:**
- Function: `u12_correction()` (line 2298)
- Scans ±200bp for strong U12-like motifs
- Adjusts boundaries if better canonical sites found
- Tags corrected introns with [c:distance]

### 3. Scoring Pipeline

**Stage 1: PWM Loading**
- Default matrices in `intronIC/data/scoring_matrices.fasta.iic`
- Separate matrices for U2/U12, canonical/non-canonical
- Custom matrices supported via `--pwms` argument

**Stage 2: Raw Score Calculation**
- Function: `get_raw_scores()` (line 3115)
- Applies PWMs to three regions per intron:
  - 5' splice site: default (-3, +9) relative to donor
  - Branch point: default (-55, -5) relative to acceptor
  - 3' splice site: default (-6, +4) relative to acceptor
- Uses `seq_score()` (line 2114) for log-odds ratios
- Uses `bp_score()` (line 2143) for optimal branch point search

**Stage 3: Z-Score Normalization**
- Function: `scale_scores()` (line 3655)
- Fit StandardScaler on reference sequences only (prevents data leakage)
- Transform: z = (raw_score - mean) / stdev
- Converts features to common scale for SVM

**Stage 4: SVM Hyperparameter Optimization**
- Function: `optimize_svm()` (line 5431)
- **Method**: 5-round grid search with geometric refinement
- **Parameter**: C (soft-margin penalty), range 10^-6 to 10^6
- **Validation**: 5-fold cross-validation
- **Metric**: balanced_accuracy
- **Process**:
  1. Round 1: Coarse grid (13 points, log scale)
  2. Rounds 2-5: Refine around best value
  3. Final: Average of rank-1 parameters

**Stage 5: SVM Training**
- Function: `train_svm()` (line 5345)
- **Model**: sklearn.svm.SVC with linear kernel
- **Split**: 80/20 training/test
- **Class balancing**: Balanced class weights
- **Optional**: Multiple models via U2 subsampling
- **Metrics**: F1 score, Precision-Recall AUC

**Stage 6: Classification**
- Function: `parallel_svm_score()` (line 5816)
- Apply trained models to experimental introns
- Multiprocessing for parallelization
- Methods: `predict_proba()`, `predict()`, `decision_function()`

**Stage 7: Ensemble Aggregation**
- Function: `average_svm_score_info()` (line 5651)
- Average probabilities across models
- Weight by F1 scores
- Set `svm_score`, `relative_score`, `type_id`

### 4. Recursive Training (Optional)

**Function:** `recursive_scoring()` (line 3890)
**Trigger:** ≥5 confident U12s found (>threshold)

**Process:**
1. Extract high-confidence U12s from first pass
2. Subsample U2 set (optional, for speed)
3. Build new PWM matrices from species-specific data
4. Re-train SVM models
5. Re-score all introns
6. Compare to first-pass results

**Benefit:** Improved classification for species distant from training data

---

## Data Flow & Pipeline

```
┌─────────────────────────┐
│   INPUT SELECTION       │
│  • Annotation + Genome  │
│  • BED File             │
│  • Sequences            │
└───────────┬─────────────┘
            │
┌───────────▼─────────────┐
│  EXTRACTION PHASE       │
│  • Parse annotations    │
│  • Create Intron objects│
│  • Extract sequences    │
│  • Calculate flanks     │
└───────────┬─────────────┘
            │
┌───────────▼─────────────┐
│  FILTERING PHASE        │
│  • Length validation    │
│  • Sequence quality     │
│  • Non-canonical check  │
│  • Duplicate detection  │
│  • Longest isoform tag  │
│  • Overlap detection    │
└───────────┬─────────────┘
            │
┌───────────▼─────────────┐
│  SCORING PHASE          │
│  • Load PWM matrices    │
│  • Calculate raw scores │
│  • Normalize z-scores   │
│  • Optimize SVM params  │
│  • Train classifier     │
│  • Apply to introns     │
│  • Aggregate scores     │
└───────────┬─────────────┘
            │
┌───────────▼─────────────┐
│  OPTIONAL RECURSIVE     │
│  (if >5 confident U12s) │
│  • Build new matrices   │
│  • Re-train SVM         │
│  • Re-score introns     │
└───────────┬─────────────┘
            │
┌───────────▼─────────────┐
│  OUTPUT GENERATION      │
│  • .meta.iic            │
│  • .bed.iic             │
│  • .seqs.iic            │
│  • .scores.iic          │
│  • Mapping files        │
│  • Plots (.png)         │
│  • Log file             │
└─────────────────────────┘
```

---

## Key Algorithms

### 1. Annotation Hierarchy Parsing
**Function:** `annotation_hierarchy()` (line 1567)
- Builds directed acyclic graph (DAG) using NetworkX
- Creates parent-child relationships from flat GFF3
- Topologically sorts to identify hierarchy
- Result: gene→transcript→exon/CDS structure

### 2. Intron Generation from Exons
**Function:** `intronator()` (line 1955)
- Sort exons in coding direction
- Sliding window of size 2 over exon list
- Create Intron via `Intron.from_exon_pair()` for each pair
- Check for overlapping exons; skip if found
- **Key insight:** Intron = space between consecutive exons

### 3. Position-Weight Matrix Scoring
**Function:** `seq_score()` (line 2114)
- **Formula:** Score = ∏(frequency[base_i][position_i])
- Pseudocount added to avoid zeros (default: 0.0001)
- Can ignore specific positions (e.g., canonical dinucleotides)
- Returns log-odds ratio: U12 frequency / U2 frequency

### 4. Branch Point Motif Search
**Function:** `bp_score()` (line 2143)
- Sliding window of matrix length across region
- Score each window with `seq_score()`
- Return max score, coordinates, and sequence
- Finds optimal branch point location within search region

### 5. Duplicate Intron Resolution
**Data structure:** Coordinate-based indexing
- Hash: (region, strand, start, stop)
- First occurrence: initialize as representative
- Subsequent: mark as duplicate, link to original
- **Priority rules:**
  1. Prefer longest-isoform introns
  2. Use parent transcript length as tiebreaker
  3. Use annotation line order as final tiebreaker

### 6. Dynamic U2 Branch Point Matrix Generation
**Function:** `build_u2_bp_matrix()` (line 2775)
- Identify likely U2 introns (low U12 scores)
- Extract their branch point regions
- Find best matches to known U12 BP motif
- Build species-specific U2 BP PWM
- Replace default human PWM with species-specific version

---

## File Formats

### Input Formats

#### 1. Genome (FASTA)
```
>chr1
ACTGACTGACTG...
>chr2
...
```
- Flag: `-g` or `--genome`
- Supports gzip compression

#### 2. Annotation (GFF3/GTF)
```
chr1  source  gene  1000  2000  .  +  .  ID=gene001
chr1  source  mRNA  1000  2000  .  +  .  ID=trans001;Parent=gene001
chr1  source  exon  1000  1200  .  +  .  Parent=trans001
```
- Flag: `-a` or `--annotation`
- Supports GFF3, GTF formats
- Supports gzip compression

#### 3. BED File
```
chr1  12000  12950  intron_1  100  +
```
- Flag: `-b` or `--bed`
- 0-based start coordinate
- 6 columns: chrom, start, stop, name, score, strand

#### 4. Sequence File (.iic)
```
intron_1  ...GGATG...  GTAATACGTAG  AGGTT...
```
- Flag: `-q` or `--sequence_file`
- Tab-delimited: name, 5'_flank, intron_seq, 3'_flank

### Output Formats

#### 1. .meta.iic (Metadata)
Tab-delimited with comprehensive metadata:
- Intron name with tags
- Relative score (vs threshold)
- Terminal dinucleotides
- Motif schematic
- Branch point context
- Intron length
- Parent transcript ID
- Grandparent gene ID
- Intron index and family size
- Fractional position
- Exon phase
- Type ID (u12/u2)

#### 2. .bed.iic (BED-like)
- Chromosome
- Start (0-based for BED compatibility)
- Stop
- Intron label (name;probability)
- SVM score (0-100)
- Strand

#### 3. .seqs.iic (Sequences)
- Intron name
- 5' flanking sequence
- Intron sequence
- 3' flanking sequence
- SVM score (optional)

#### 4. .scores.iic (Detailed Scoring)
Per-intron breakdown:
- Name, relative score, SVM score
- Score distance to decision boundary
- 5' sequence, raw score, z-score
- Branch point sequences (U12 & U2 versions)
- BP raw score, z-score
- 3' sequence, raw score, z-score

#### 5. Mapping Files
- **.dupe_map.iic**: Maps duplicates to representative
- **.overlap_map.iic**: Maps overlapping introns

#### 6. Visualization (.png)
- Scatter plots (2D/3D) of training data and classifications
- Hexbin density plots of reference introns
- Precision-Recall AUC curves

---

## Code Organization

### File Structure
```
intronIC/
├── intronIC/
│   ├── __init__.py
│   ├── intronIC.py              # Main module (6,093 lines)
│   ├── pwms_from_seqs.py        # PWM generation utility
│   ├── matrices_from_seqs.py    # Matrix generation utility
│   ├── _version.py              # Version management
│   ├── data/                    # Training data
│   │   ├── scoring_matrices.fasta.iic
│   │   ├── u2_reference.introns.iic.gz (5.4 MB)
│   │   ├── u12_reference.introns.iic.gz (144 KB)
│   │   └── [additional matrix variants]
│   └── test_data/               # Example data
│       ├── Homo_sapiens.Chr19.Ensembl_91.fa.gz
│       └── Homo_sapiens.Chr19.Ensembl_91.gff3.gz
├── setup.py
├── setup.cfg
├── README.md
├── LICENSE (GPL v3.0)
└── versioneer.py
```

### Main Function Organization

**Entry Point:** `main()` (line 5038)

**Flow:**
1. Initialization (5039-5052): Args, multiprocessing setup
2. Logging setup (5055-5073): File and console loggers
3. Reference loading (5108-5119): U2/U12 training data
4. Intron collection (5121-5132): Parse inputs, create objects
5. Filtering (5124-5130): Quality control, deduplication
6. Scoring (5145-5175): PWM scoring, SVM training/application
7. Optional recursive pass (5179-5187): Second training round
8. Output generation (5189+): Write files, plots, statistics

### Critical Code Locations

| Component | Function | Line |
|-----------|----------|------|
| **Entry point** | `main()` | 5038 |
| **Annotation parsing** | `annotation_hierarchy()` | 1567 |
| **Intron creation** | `intronator()` | 1955 |
| **Intron from annotation** | `introns_from_annotation()` | 4102 |
| **Intron from BED** | `introns_from_bedfile()` | 4161 |
| **Sequence extraction** | `get_sub_seqs()` | 2645 |
| **Filtering** | `filter_introns_write_files()` | 4668 |
| **Omission check** | `omit_check()` | 671 |
| **Duplicate tagging** | `add_tags()` | 1897 |
| **NC adjustment** | `u12_correction()` | 2298 |
| **PWM scoring** | `seq_score()` | 2114 |
| **BP search** | `bp_score()` | 2143 |
| **Raw scores** | `get_raw_scores()` | 3115 |
| **Normalization** | `scale_scores()` | 3655 |
| **SVM optimization** | `optimize_svm()` | 5431 |
| **SVM training** | `train_svm()` | 5345 |
| **SVM scoring** | `parallel_svm_score()` | 5816 |
| **Score aggregation** | `average_svm_score_info()` | 5651 |
| **Recursive training** | `recursive_scoring()` | 3890 |
| **Reference loading** | `get_reference_introns()` | 1423 |
| **U2 BP matrix** | `build_u2_bp_matrix()` | 2775 |

---

## Usage & Configuration

### Basic Usage

```bash
# Standard classification
intronIC -g genome.fa.gz -a annotation.gff3.gz -n homo_sapiens

# Extract sequences only (no classification)
intronIC -g genome.fa.gz -a annotation.gff3.gz -n homo_sapiens -s

# Parallel processing (8 cores)
intronIC -g genome.fa.gz -a annotation.gff3.gz -n homo_sapiens -p 8

# Stricter threshold (95%)
intronIC -g genome.fa.gz -a annotation.gff3.gz -n homo_sapiens -t 95
```

### Advanced Options

```bash
# Include alternative isoforms
intronIC -g genome.fa.gz -a annotation.gff3.gz -n species -i

# Exclude non-canonical introns from scoring
intronIC -g genome.fa.gz -a annotation.gff3.gz -n species --no_nc

# Recursive training for distant species
intronIC -g genome.fa.gz -a annotation.gff3.gz -n species --recursive

# Custom reference sequences
intronIC -g genome.fa.gz -a annotation.gff3.gz -n species \
  --reference_u12s custom_u12s.iic --reference_u2s custom_u2s.iic

# Custom scoring region coordinates
intronIC -g genome.fa.gz -a annotation.gff3.gz -n species \
  --five_score_coords -5 10 --three_score_coords -8 5
```

### Key Arguments

**Required:**
- `-n, --species_name`: Binomial name (e.g., homo_sapiens)

**Input Selection (one required):**
- `-g, --genome`: FASTA genome file
- `-a, --annotation`: GFF3/GTF annotation file
- `-b, --bed`: BED file with coordinates
- `-q, --sequence_file`: Pre-extracted sequences

**Scoring Options:**
- `-s, --sequences_only`: Skip classification
- `-t, --threshold`: U12 probability threshold (0-100, default: 90)
- `-f, --feature`: Feature type (cds/exon, default: both)
- `--no_nc`: Exclude non-canonical introns

**Performance:**
- `-p, --processes N`: Parallel processes
- `--min_intron_len N`: Minimum length (default: 30)

**Isoform Selection:**
- `-i, --allow_multiple_isoforms`: Include non-longest
- `-v, --no_intron_overlap`: Exclude overlapping

**Training:**
- `--pwms file`: Custom PWM files
- `--generate_u2_bps_pwm`: Build U2 BP PWM from data
- `--recursive [subset]`: Second training round
- `-C value`: Fixed SVM C parameter

---

## Performance Considerations

### Memory Usage

**Per-Intron Memory:**
- Base: ~50-100 bytes per Intron object (41 slots)
- With sequences: +200-500 bytes average
- Large genomes: Millions of introns = several GB

**Optimization Strategies:**
- `__slots__` prevents dynamic attribute addition
- Generator pattern for sequence extraction
- Sequences cleared after output (line 4847)
- Region-based indexing reduces search space

**Human Genome:**
- ~100M+ introns with sequences
- Peak memory: ~5GB
- Sequences cleared during processing to stay manageable

### Runtime Performance

**Complexity:**
- Parsing: O(n) where n = annotation lines
- Extraction: O(n) where n = exons
- Scoring: O(n*m) where m = intron count
- SVM training: O(k²*log(k)) where k = reference count

**Typical Performance:**
- Human genome: 20-35 minutes (single process)
- Parallelization: 4-6x speedup with 8 processes
- Chr19 only: ~2-3 minutes

**Bottlenecks:**
1. Sequence extraction from large genomes
2. PWM scoring across all introns
3. SVM hyperparameter optimization (5 rounds)
4. Plot generation for large datasets

### Parallelization

**Supported:**
- Intron scoring: `-p N` for N processes
- Cross-validation: `--cv_processes N`

**Not Parallelized:**
- Annotation parsing (I/O bound)
- Filtering and tagging (sequential)
- SVM training (single model)
- Output writing (I/O bound)

---

## Future Improvement Areas

### Code Quality & Architecture

**Current State:**
- Written by PhD student with no formal CS background
- Single monolithic file (6,093 lines)
- Complex functions with multiple responsibilities
- Some deeply nested logic

**Potential Improvements:**

1. **Modularization**
   - Split `intronIC.py` into logical modules:
     - `models.py` - Class definitions
     - `parsers.py` - Input parsing functions
     - `extractors.py` - Intron extraction logic
     - `scorers.py` - PWM and scoring functions
     - `classifiers.py` - SVM training and prediction
     - `filters.py` - Quality control and filtering
     - `writers.py` - Output generation
   - Benefit: Easier testing, maintenance, understanding

2. **Function Simplification**
   - Break down large functions (e.g., `main()`, `filter_introns_write_files()`)
   - Extract helper functions for complex logic blocks
   - Reduce nesting depth (currently some 4-5 levels deep)
   - Benefit: Improved readability, easier debugging

3. **Type Hints**
   - Add type annotations throughout
   - Use `typing` module for complex types
   - Consider using `mypy` for type checking
   - Benefit: Better IDE support, catches bugs early

4. **Error Handling**
   - More specific exception types
   - Better error messages for users
   - Graceful degradation for non-critical failures
   - Benefit: Better user experience, easier troubleshooting

5. **Documentation**
   - Add comprehensive docstrings (many functions lack them)
   - Document algorithm complexity
   - Add inline comments for complex logic
   - Benefit: Easier onboarding, better maintenance

### Performance Optimizations

1. **Caching**
   - Cache PWM matrix loading
   - Cache genome sequence lookups
   - Memoize repeated calculations
   - Benefit: Faster repeated runs, reduced I/O

2. **Vectorization**
   - Use NumPy operations for batch scoring
   - Vectorize PWM scoring where possible
   - Benefit: 2-10x speedup for scoring phase

3. **Memory Optimization**
   - Stream processing for very large genomes
   - More aggressive sequence clearing
   - Consider using memory-mapped files for genome
   - Benefit: Handle larger datasets, reduce peak memory

4. **Parallel Parsing**
   - Parallelize annotation parsing (chunk by chromosome)
   - Parallel sequence extraction
   - Benefit: Better CPU utilization, faster overall runtime

### Testing

**Current State:**
- Test data provided (Chr19)
- No automated test suite

**Recommended:**
1. Unit tests for core functions
2. Integration tests for pipeline stages
3. Regression tests for classification accuracy
4. Performance benchmarks
5. Consider pytest framework

### Algorithm Improvements

1. **Advanced Branch Point Detection**
   - Current: Simple sliding window
   - Potential: Machine learning-based BP prediction
   - Benefit: More accurate BP identification

2. **Feature Engineering**
   - Current: 3 z-scores (five, bp, three)
   - Potential: Additional features (intron length, GC content, context)
   - Benefit: Potentially better classification

3. **Deep Learning Alternative**
   - Current: Linear SVM
   - Potential: CNN or transformer-based classifier
   - Benefit: May capture more complex patterns

4. **Active Learning**
   - Current: Fixed training set
   - Potential: Iteratively refine training data
   - Benefit: Better performance on specific species

### Usability

1. **Configuration Files**
   - Support YAML/JSON config files
   - Preset configurations for common organisms
   - Benefit: Easier repeated analysis, reproducibility

2. **Resume Capability**
   - Checkpoint system for long runs
   - Resume from interruption
   - Benefit: Robustness for large datasets

3. **Interactive Mode**
   - Simple GUI or CLI wizard for basic usage
   - Parameter tuning interface
   - Benefit: Lower barrier to entry

---

## Important Notes for Future Development

### Critical Logic - Handle with Care

These areas involve complex intron extraction and metadata tracking:

1. **Coordinate Systems**
   - Mixed 0-based (BED input) and 1-based (internal) coordinates
   - Careful adjustment in `introns_from_bedfile()` (line 4161)
   - Off-by-one errors fixed in recent versions (see git history)

2. **Strand Handling**
   - Reverse complement logic for negative strand
   - Coordinate ordering differs by strand
   - Critical in `get_sub_seqs()` and `intronator()`

3. **Non-Canonical Boundary Adjustment**
   - Complex search algorithm in `u12_correction()` (line 2298)
   - Adjusts coordinates based on motif strength
   - Must maintain metadata consistency after adjustment

4. **Duplicate Resolution**
   - Priority system must be consistent
   - Order matters for reproducibility
   - Affects downstream statistics

5. **Isoform Selection**
   - Longest isoform identification affects many downstream steps
   - Overlapping introns from different isoforms need careful handling
   - Critical for avoiding false positives

### Dependencies to Monitor

- **numpy < 2.0**: Version constraint due to API changes
- **scikit-learn >= 0.22**: Required for SVM features
- **networkx >= 2.5.1**: For graph operations
- **biogl**: Custom library, may need maintenance

### Known Issues/Limitations

1. Very large genomes (>5GB) may require HPC resources
2. Non-model organisms may benefit from recursive training
3. Annotation quality directly affects intron extraction accuracy
4. PWM matrices trained primarily on vertebrate data

---

## Conclusion

intronIC is a sophisticated, production-ready bioinformatics tool with a solid algorithmic foundation. While the codebase could benefit from modularization and additional documentation, the core logic for intron extraction, classification, and metadata tracking is comprehensive and has been refined through use in published research.

The classification approach (PWM + SVM) is well-justified for this problem domain, and the system handles edge cases (non-canonical introns, multiple isoforms, duplicates) thoughtfully. Future development should focus on code organization and maintainability while preserving the carefully implemented extraction and labeling logic.

**Key Strengths:**
- Comprehensive intron extraction with proper metadata tracking
- Robust filtering and quality control
- Well-tuned SVM classifier with hyperparameter optimization
- Flexible input/output formats
- Parallel processing support
- Extensive output for downstream analysis

**Areas for Enhancement:**
- Code modularization and documentation
- Test coverage
- Performance optimization through vectorization
- Type hints and error handling
- Configuration file support

---

**Citation:**
Moyer et al. (2020): "Comprehensive database and evolutionary dynamics of U12-type introns"
*Nucleic Acids Research*, 48(13):7066–7078
https://doi.org/10.1093/nar/gkaa464
