# intronIC Output Format Reference

**Version:** 1.5.1 (Refactored)
**Last Updated:** 2025-11-14

This document provides comprehensive specifications for all intronIC output file formats.

---

## Table of Contents

1. [Overview](#overview)
2. [Metadata Format (.meta.iic)](#metadata-format-metaiic)
3. [BED Format (.bed.iic)](#bed-format-bediic)
4. [Sequence Format (.introns.iic)](#sequence-format-intronsiic)
5. [Score Details Format (.score_info.iic)](#score-details-format-score_infoiic)
6. [Mapping Files](#mapping-files)
7. [Visualization Files](#visualization-files)
8. [Common Elements](#common-elements)

---

## Overview

All intronIC output files use:
- **Tab-delimited** text format (`.iic` extension)
- **UTF-8** encoding
- **Header row** (except BED format)
- **Dot (.)** as null/missing value placeholder
- File naming: `{species_name}.{type}.iic`

### File Purposes

| File | Purpose | Use Cases |
|------|---------|-----------|
| `.meta.iic` | Comprehensive metadata | Filtering, statistics, downstream analysis |
| `.bed.iic` | Genomic coordinates | Genome browsers, intersection analysis |
| `.introns.iic` | Sequences | Re-analysis, BLAST, motif search |
| `.score_info.iic` | Detailed scores | Score distribution, debugging, refinement |
| `.dupe_map.iic` | Duplicate mapping | Tracing alternative isoforms |
| `.overlap_map.iic` | Overlap mapping | Resolving coordinate conflicts |

---

## Metadata Format (.meta.iic)

**Filename:** `{species_name}.meta.iic`
**Purpose:** Comprehensive intron metadata for downstream analysis

### Column Specifications

| # | Column | Type | Description | Example |
|---|--------|------|-------------|---------|
| 1 | `name` | string | Full intron identifier with tags | `HomSap-gene:ENSG00000196218@transcript:ENST00000355481-intron_69(104)` |
| 2 | `rel_score` | float | Score relative to threshold (svm_score - threshold) | `5.2341` |
| 3 | `dnts` | string | Terminal dinucleotides (5'-3') | `GT-AG`, `AT-AC`, `GC-AG` |
| 4 | `motif_schematic` | string | Splice site and BP motif schematic | `AAG\|GTCGGGGCTT...TACTAAC/CACAG...TTTAG\|TCC` |
| 5 | `bp_context` | string | Branch point sequence with context | `TTGACAGGCAGTGATAT[TACTAAC]GACTGAGTTTAG` |
| 6 | `length` | int | Intron length in base pairs | `1523` |
| 7 | `parent` | string | Parent transcript ID | `transcript:ENST00000355481` |
| 8 | `grandparent` | string | Grandparent gene ID | `gene:ENSG00000196218` |
| 9 | `index` | int | Intron ordinal position (1-based) | `69` |
| 10 | `family_size` | int | Total introns in transcript | `104` |
| 11 | `frac_pos` | float | Fractional position in transcript (0-1) | `0.6635` |
| 12 | `phase` | int | Exon phase (0, 1, 2) | `0` |
| 13 | `type_id` | string | Classification: `u2`, `u12`, or `.` | `u12` |
| 14 | `feature` | string | Source feature type | `cds`, `exon` |
| 15 | `attributes` | string | Comma-separated verbose attributes | `noncanonical,not_longest_isoform` |

### Column Details

#### 1. name
Full intron identifier with embedded metadata:

**Format:** `{species}-{grandparent}@{parent}-intron_{index}({family_size}){tags}`

**Components:**
- `species`: 6-char abbreviation (3 genus + 3 species, e.g., "HomSap" from "Homo sapiens")
- `grandparent`: Gene ID (often with "gene:" prefix)
- `@`: Separator between gene and transcript
- `parent`: Transcript ID (often with "transcript:" prefix)
- `index`: Intron ordinal position (1-based)
- `family_size`: Total introns in transcript
- `tags`: Optional tags (see [Tags](#tags) section)

**Examples:**
```
HomSap-gene:ENSG00000196218@transcript:ENST00000355481-intron_69(104)
HomSap-gene:ENSG00000196218@transcript:ENST00000355481-intron_69(104);[o:i];[i]
DroMel-gene:FBgn0000123@transcript:FBtr0000456-intron_1(5);[c:3]
```

#### 2. rel_score
**Relative Score** = `svm_score - threshold`

- **Positive values**: Above threshold (classified as U12 at chosen confidence)
- **Negative values**: Below threshold (classified as U2)
- **Zero**: Exactly at threshold
- **Dot (.)**: Not scored (omitted intron)

**Purpose:** Makes filtering easier - just check if `rel_score > 0` to find U12-type introns.

**Example filtering:**
```bash
# Get all U12-type introns at chosen threshold
awk '($2!="." && $2>0)' species.meta.iic

# Get high-confidence U12s (>10 points above threshold)
awk '($2!="." && $2>10)' species.meta.iic
```

#### 3. dnts
Terminal dinucleotides showing splice site types.

**Format:** `{5'_dinucleotide}-{3'_dinucleotide}`

**Common types:**
- `GT-AG`: Canonical major spliceosome (~99% of introns)
- `GC-AG`: Non-canonical major spliceosome (~0.5%)
- `AT-AC`: Canonical minor spliceosome (~0.05% of U12 introns)
- `GT-AG`: Also used by GT-AG type U12 introns (~0.45% of U12 introns)

#### 4. motif_schematic
Visual representation of key splice motifs.

**Format:** `{exon_3bp}|{5'_10bp}...{bp_u12}/{bp_u2}...{3'_display}|{exon_3bp}`

**Components:**
- `exon_3bp`: Last 3bp of upstream exon, first 3bp of downstream exon
- `|`: Exon-intron boundary
- `5'_10bp`: First 10bp of intron (5' splice site region)
- `...`: Separator indicating skipped sequence
- `bp_u12/bp_u2`: Branch point sequences (U12 and U2 versions)
- `3'_display`: Intron sequence from BP end to 3' splice site

**Example:**
```
AAG|GTCGGGGCTT...TACTAAC/CACAG...TTTAG|TCC
│   │          │  │             │     │   │
│   │          │  │             │     │   └─ First 3bp of downstream exon
│   │          │  │             │     └───── 3' splice site region
│   │          │  │             └─────────── Alternative U2 BP
│   │          │  └───────────────────────── U12 BP consensus
│   │          └──────────────────────────── Separator
│   └─────────────────────────────────────── 5' splice site (first 10bp)
└─────────────────────────────────────────── Last 3bp of upstream exon
```

#### 5. bp_context
Branch point sequence in genomic context with brackets highlighting the BP.

**Format:** `{bp_region}[{bp_seq}]{bp_region_right}{3'_display_seq}`

**Example:**
```
TTGACAGGCAGTGATAT[TACTAAC]GACTGAGTTTAG
                  └───────┘
                  Branch point sequence
```

#### 11. frac_pos
Fractional position of intron within transcript (0 = first, 1 = last).

**Calculation:** `index / (family_size + 1)`

**Use cases:**
- Analyzing positional bias of U12 introns (often enriched in 5' regions)
- Studying intron evolution patterns

#### 15. attributes
Verbose, human-readable attributes (new in refactored version).

**Possible values:**
- `noncanonical`: Non-canonical terminal dinucleotides
- `not_longest_isoform`: Not from longest transcript isoform
- `corrected`: Splice site boundaries were adjusted
- `duplicate`: Duplicate coordinates (excluded from analysis)
- `omitted_short`: Too short (below minimum length)
- `omitted_ambiguous`: Contains ambiguous bases (N, etc.)
- `omitted_noncanonical`: Omitted due to non-canonical sites
- `omitted_overlap`: Overlapping coordinates
- `omitted_not_longest_isoform`: Omitted due to isoform filtering

**Examples:**
```
noncanonical
noncanonical,not_longest_isoform
corrected
omitted_short
.
```

### Example Entry

```
HomSap-gene:ENSG00000141510@transcript:ENST00000268261-intron_1(2)	15.2341	AT-AC	TGG|ATACCCTTCT...TCCTTAAC/CTTTG...AACAG|GGC	TGACAGGCAGTGATAT[TCCTTAAC]GACTGAGTTTAG	1523	transcript:ENST00000268261	gene:ENSG00000141510	1	2	0.3333	0	u12	cds	.
```

---

## BED Format (.bed.iic)

**Filename:** `{species_name}.bed.iic`
**Purpose:** Genomic coordinates for genome browsers and intersection tools

### Column Specifications

| # | Column | Type | Description | Example |
|---|--------|------|-------------|---------|
| 1 | `chrom` | string | Chromosome/contig name | `chr19`, `NC_000019.10` |
| 2 | `start` | int | Start coordinate (0-based, BED standard) | `12000` |
| 3 | `stop` | int | Stop coordinate (1-based, BED standard) | `13523` |
| 4 | `name` | string | Intron label with score and tags | `HomSap-ENST00000355481_69(104);95.23;[n]` |
| 5 | `score` | float | SVM probability score (0-100) | `95.23` |
| 6 | `strand` | char | Genomic strand | `+`, `-` |
| 7 | `attributes` | string | Verbose attributes (refactored addition) | `noncanonical` |

### Column Details

#### 2-3. start, stop
**Important:** BED format uses **0-based, half-open** coordinates:
- `start`: 0-based (first base is position 0)
- `stop`: 1-based (not inclusive)

**Example:**
```
Sequence:  A C G T A C G T
Position:  0 1 2 3 4 5 6 7

BED for "CGTA": start=1, stop=5
```

**Conversion to/from 1-based:**
```
bed_start = one_based_start - 1
bed_stop = one_based_stop
```

#### 4. name
Compact intron label for visualization.

**Format:** `[{species}-]{parent}_{index}({family_size});{svm_score}[;{tags}]`

**Examples:**
```
HomSap-ENST00000355481_69(104);95.23
HomSap-ENST00000355481_69(104);95.23;[n]
HomSap-ENST00000355481_69(104);95.23;[n];[i]
```

### Usage Examples

**Load into UCSC Genome Browser:**
```bash
# Sort and compress for indexing
sort -k1,1 -k2,2n species.bed.iic | bgzip > species.bed.iic.gz
tabix -p bed species.bed.iic.gz
```

**Find overlaps with bedtools:**
```bash
bedtools intersect -a species.bed.iic -b genes.bed -wa -wb
```

**Extract U12-type introns:**
```bash
awk '$5 > 90' species.bed.iic > u12_introns.bed
```

---

## Sequence Format (.introns.iic)

**Filename:** `{species_name}.introns.iic` (or `.seqs.iic`)
**Purpose:** Intron sequences with flanking regions for re-analysis

### Column Specifications

| # | Column | Type | Description | Example |
|---|--------|------|-------------|---------|
| 1 | `name` | string | Full intron identifier | `HomSap-gene:ENSG00000196218@...` |
| 2 | `score` | float | SVM score (optional, with `-s` flag) | `95.23` or `.` |
| 3 | `upstream_flank` | string | 5' exonic flanking sequence | `ATGCGTAGCTAGC...` |
| 4 | `sequence` | string | Intron sequence | `GTACGTACGT...` |
| 5 | `downstream_flank` | string | 3' exonic flanking sequence | `...GCTAGCTAGCT` |

### Column Details

#### 2. score
- **Present** when classification is performed (default)
- **Absent** when using `-s` (sequences-only mode)
- **Dot (.)** for omitted introns that weren't scored

#### 3. upstream_flank
Exonic sequence immediately upstream (5') of intron.

**Default length:** Variable (typically 20-50bp, depending on exon size)

**Purpose:**
- Provides context for 5' splice site
- Useful for codon phase analysis
- Required for some motif finding algorithms

#### 5. downstream_flank
Exonic sequence immediately downstream (3') of intron.

**Default length:** Variable (typically 20-50bp)

### Usage Examples

**Extract U12 sequences for BLAST:**
```bash
# Get U12 introns with scores > 90
awk '$2 > 90 {print ">"$1"\n"$4}' species.introns.iic > u12_seqs.fasta
```

**Re-score with different parameters:**
```bash
# First extract sequences without classification
intronIC -g genome.fa -a annotation.gff3 -n species -s

# Then classify with different threshold
intronIC -q species.introns.iic -n species -t 95
```

**Search for specific motifs:**
```bash
# Find introns with TCCTTAAC branch point
awk '$4 ~ /TCCTTAAC/ {print $1}' species.introns.iic
```

---

## Score Details Format (.score_info.iic)

**Filename:** `{species_name}.score_info.iic`
**Purpose:** Detailed scoring breakdown for each intron

### Column Specifications

| # | Column | Type | Description | Example |
|---|--------|------|-------------|---------|
| 1 | `name` | string | Full intron identifier | `HomSap-gene:ENSG00000196218@...` |
| 2 | `rel_score` | float | Relative score (svm_score - threshold) | `5.2341` |
| 3 | `svm_score` | float | SVM probability score (0-100) | `95.23` |
| 4 | `5'_seq` | string | 5' splice site sequence scored | `GTCGGGGCTT` |
| 5 | `5'_raw` | float | 5' raw PWM score (log-odds) | `0.123456` |
| 6 | `5'_z` | float | 5' normalized z-score | `2.3451` |
| 7 | `bp_seq` | string | Branch point sequence (U12) | `TACTAAC` |
| 8 | `bp_seq_u2` | string | Branch point sequence (U2) | `CACAG` |
| 9 | `bp_raw` | float | BP raw PWM score (log-odds) | `-1.234567` |
| 10 | `bp_z` | float | BP normalized z-score | `-0.8765` |
| 11 | `3'_seq` | string | 3' splice site sequence scored | `TTTAG` |
| 12 | `3'_raw` | float | 3' raw PWM score (log-odds) | `0.987654` |
| 13 | `3'_z` | float | 3' normalized z-score | `1.2345` |
| 14 | `decision_dist` | float | Distance to SVM decision boundary | `2.1234` |

### Column Details

#### 4, 7, 11. Scored sequences
The actual sequence segments that were scored by PWM matrices.

**Typical lengths:**
- `5'_seq`: 12bp (e.g., -3 to +9 relative to donor site)
- `bp_seq`: 7bp (consensus sequence found in search window)
- `3'_seq`: 10bp (e.g., -6 to +4 relative to acceptor site)

**Note:** Lengths depend on PWM matrix dimensions and can be customized.

#### 5, 9, 12. Raw scores
**Raw PWM scores** are log-odds ratios:

**Formula:** `log10(P(seq|U12) / P(seq|U2))`

**Interpretation:**
- **Positive**: Sequence favors U12 model
- **Negative**: Sequence favors U2 model
- **Magnitude**: Strength of preference

**Important:** Raw scores are NOT directly comparable between regions (5', BP, 3') because:
- Different matrix sizes
- Different information content
- Different score ranges

#### 6, 10, 13. Z-scores
**Normalized scores** (z-scores) standardize raw scores for comparison.

**Formula:** `z = (raw_score - mean) / std_dev`

**Properties:**
- Mean = 0, StdDev = 1 (for reference sequences)
- Directly comparable across regions
- Used as features for SVM classifier

**Interpretation:**
- `z > 2`: Strong U12-like signal (>95th percentile)
- `z > 3`: Very strong U12-like signal (>99.7th percentile)
- `z < -2`: Strong U2-like signal

#### 14. decision_dist
**Decision distance** = SVM decision function output (log-odds).

**Calculation:** `log(p / (1 - p))` where `p = svm_score / 100`

**Interpretation:**
- `dist = 0`: Equal probability (50%)
- `dist > 0`: Favors U12
- `dist < 0`: Favors U2
- `|dist| > 2`: High confidence

**Relationship to svm_score:**
```
svm_score = 100 / (1 + exp(-decision_dist))
```

### Usage Examples

**Plot score distributions:**
```bash
# Extract z-scores for plotting
awk 'NR>1 {print $6, $10, $13}' species.score_info.iic > z_scores.txt
```

**Find introns with specific patterns:**
```bash
# High 5' score, low BP score (unusual for U12)
awk 'NR>1 && $6 > 3 && $10 < 0' species.score_info.iic
```

**Debug classification:**
```bash
# Show all scores for specific intron
grep "intron_69" species.score_info.iic
```

---

## Mapping Files

### Duplicate Map (.dupe_map.iic)

**Filename:** `{species_name}.dupe_map.iic`
**Purpose:** Map duplicate introns to their representative

**Format:** Tab-delimited, no header
```
{representative_name}	{duplicate_name}
```

**Columns:**
1. `representative`: Name of the intron kept in analysis
2. `duplicate`: Name of the duplicate intron (excluded)

**Example:**
```
HomSap-gene:ENSG001@transcript:ENST001-intron_5(10)	HomSap-gene:ENSG001@transcript:ENST002-intron_5(12)
HomSap-gene:ENSG001@transcript:ENST001-intron_5(10)	HomSap-gene:ENSG001@transcript:ENST003-intron_3(8)
```

**Interpretation:** Multiple transcripts may produce introns at identical genomic coordinates. intronIC keeps one representative and maps others to it.

**Priority rules:**
1. Longest isoform preferred
2. Longer transcript preferred (tiebreaker)
3. First in annotation (final tiebreaker)

### Overlap Map (.overlap_map.iic)

**Filename:** `{species_name}.overlap_map.iic`
**Purpose:** Map overlapping intron coordinates

**Format:** Tab-delimited, no header
```
{representative_name}	{overlapping_name}
```

**Example:**
```
HomSap-gene:ENSG001@transcript:ENST001-intron_5(10)	HomSap-gene:ENSG002@transcript:ENST004-intron_2(7)
```

**Interpretation:** Some introns have partially overlapping coordinates (e.g., alternative splicing, overlapping genes). One is kept as representative.

---

## Visualization Files

### Scatter Plot (.plot.scatter.iic.png)

**Filename:** `{species_name}.plot.scatter.iic.png`
**Content:** 2D scatter plot of classified introns

**Axes:**
- X: 5' splice site z-score
- Y: Branch point z-score
- Color: U12 (red) vs U2 (blue)
- Size: SVM confidence

**Marginal distributions:** Histograms along each axis

**Purpose:** Visualize separation between U12 and U2 introns in feature space

### Training Scatter (.training_scatter.png)

**Filename:** `{species_name}.training_scatter.png`
**Content:** Scatter plot of training data

**Purpose:** Visualize reference sequences used for training

### Training Hexplot (.training_hexplot.png)

**Filename:** `{species_name}.training_hexplot.png`
**Content:** Hexbin density plot of reference introns

**Purpose:** Show density distribution of training data

### Precision-Recall Curve (.pr_curve.png)

**Filename:** `{species_name}.pr_curve.png`
**Content:** Precision-Recall curve from cross-validation

**Purpose:** Evaluate classifier performance during training

---

## Common Elements

### Tags

Tags appear in intron names and provide additional metadata.

#### Omission Tags

**Format:** `;[o:{code}]`

| Code | Meaning | Description |
|------|---------|-------------|
| `o:s` | Short | Length below minimum threshold |
| `o:a` | Ambiguous | Contains non-ACTG bases (N, etc.) |
| `o:n` | Noncanonical | Omitted due to non-standard dinucleotides |
| `o:v` | Overlap | Coordinates overlap with another intron |
| `o:i` | Isoform | Not in longest transcript isoform |

**Example:** `HomSap-gene:ENSG001@transcript:ENST001-intron_5(10);[o:i]`

#### Property Tags

**Format:** `;[{code}]` or `;[{code}:{value}]`

| Tag | Meaning | Description |
|-----|---------|-------------|
| `[n]` | Noncanonical | Non-canonical terminal dinucleotides |
| `[i]` | Isoform | Not from longest transcript |
| `[c:N]` | Corrected | Boundary adjusted by N bases |
| `[d]` | Duplicate | Duplicate coordinates |

**Examples:**
```
;[n]                    # Noncanonical
;[n];[i]               # Noncanonical and not longest isoform
;[c:5]                 # Boundary corrected by 5 bases
;[o:i];[i]             # Omitted due to isoform, tagged as not longest
```

### Null Values

**Dot (.)** represents missing or not applicable values.

**Common situations:**
- Intron not scored (omitted before scoring phase)
- Optional field not populated
- Calculation not applicable

**Example:**
```
# Omitted intron with no scores
HomSap-...-intron_5(10);[o:s]  .  GT-AG  .  .  25  ...  .
```

### Coordinate Systems

**Internal (1-based):** Used in most output formats
- First base of chromosome is position 1
- Ranges are inclusive [start, stop]

**BED format (0-based, half-open):**
- First base is position 0
- Stop position is exclusive [start, stop)
- Standard for genome browsers

### Species Abbreviation

**Format:** 3-char genus + 3-char species (capitalized)

**Examples:**
```
homo_sapiens          → HomSap
drosophila_melanogaster → DroMel
caenorhabditis_elegans → CaeEle
arabidopsis_thaliana  → AraTha
c_elegans             → CEle
```

---

## File Sizes and Compression

### Expected File Sizes

Approximate sizes for human genome (~1M introns):

| File | Uncompressed | Gzipped | Notes |
|------|--------------|---------|-------|
| `.meta.iic` | ~200 MB | ~40 MB | Most comprehensive |
| `.bed.iic` | ~100 MB | ~20 MB | Minimal metadata |
| `.introns.iic` | ~500 MB | ~120 MB | Full sequences |
| `.score_info.iic` | ~150 MB | ~30 MB | Numerical data |
| `.dupe_map.iic` | ~5 MB | ~1 MB | Sparse |
| `.overlap_map.iic` | ~1 MB | ~200 KB | Very sparse |

### Compression Recommendations

**For long-term storage:**
```bash
gzip species.*.iic
```

**For indexed access (BED):**
```bash
sort -k1,1 -k2,2n species.bed.iic | bgzip > species.bed.iic.gz
tabix -p bed species.bed.iic.gz
```

---

## Parsing Examples

### Python

```python
import csv

# Read metadata file
with open('species.meta.iic') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if row['rel_score'] != '.' and float(row['rel_score']) > 0:
            print(f"U12 intron: {row['name']} (score: {row['rel_score']})")

# Read score details
with open('species.score_info.iic') as f:
    reader = csv.DictReader(f, delimiter='\t')
    z_scores = [(row['5\'_z'], row['bp_z'], row['3\'_z']) for row in reader]
```

### R

```r
# Read metadata
meta <- read.delim('species.meta.iic', stringsAsFactors=FALSE)

# Filter U12 introns
u12 <- meta[!is.na(meta$rel_score) & meta$rel_score > 0, ]

# Plot score distribution
scores <- read.delim('species.score_info.iic')
plot(scores$five_z, scores$bp_z, col=ifelse(scores$svm_score > 90, 'red', 'blue'))
```

### Bash/awk

```bash
# Count U12 introns
awk '($2!="." && $2>0)' species.meta.iic | wc -l

# Extract U12 sequences
awk '$2 > 90 {print ">"$1"\n"$4}' species.introns.iic > u12.fasta

# Get high-confidence AT-AC introns
awk '($2 > 10 && $3 == "AT-AC")' species.meta.iic
```

---

## Version History

### Version 1.5.1 (Refactored)
- Added `attributes` column to `.meta.iic` and `.bed.iic`
- Improved tag formatting consistency
- Enhanced documentation

### Version 1.5.0
- Original monolithic implementation
- Established core format specifications

---

## See Also

- [README.md](../README.md) - Quick start and basic usage
- [CLAUDE.md](../CLAUDE.md) - Comprehensive code documentation
- [writers.py](../file_io/writers.py) - Implementation details
- [intronIC wiki](https://github.com/glarue/intronIC/wiki) - Additional guides

---

**Questions or issues?** Report at [GitHub Issues](https://github.com/glarue/intronIC/issues)
