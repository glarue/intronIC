# Computational Methods for Intron Classification: Field Approaches & intronIC Comparison

**Last Updated:** 2025-11-01
**Purpose:** Reference guide for U12/U2 intron classification methods used in the field and comparison to intronIC's implementation

---

## Table of Contents

1. [Overview of Computational Challenge](#overview-of-computational-challenge)
2. [Classical Methods (Pre-intronIC)](#classical-methods-pre-intronic)
3. [intronIC's Approach](#intronics-approach)
4. [Comparative Analysis](#comparative-analysis)
5. [Recommendations for Development](#recommendations-for-development)

---

## Overview of Computational Challenge

The computational distinction between major (U2-type) and minor (U12-type) introns is fundamentally based on the **highly conserved sequence motifs** of the U12-type spliceosome, which contrast sharply with the more variable U2-type signals.

### Key Biological Features

| Feature | U2-type (Major) | U12-type (Minor) |
|---------|-----------------|------------------|
| **Frequency** | ~99.5% of introns | ~0.5% of introns |
| **5' Splice Site** | KAG GTRAGT (variable) | More conserved |
| **Branch Point** | CURACU (degenerate) | TTCCTTRACYCY or UUUCCUUAACYCY (highly conserved) |
| **BP-to-3'ss Distance** | Variable | **Sharply restricted: 8-30 nt** |
| **Polypyrimidine Tract** | Present (YYYYYYYYYY) | **Absent or weak** |
| **3' Splice Site** | YYYYYYYYYYYYoYAG | Goo (lacks strong PPT) |
| **Terminal Dinucleotides** | GT-AG (~99.6%) | GT-AG (70%), AT-AC, others |

### Computational Challenge

- **Rarity:** U12 introns are <1% of total, creating severe class imbalance
- **Variability:** Even U12 introns show some sequence variation
- **Overlap:** GT-AG U12 introns use same terminal dinucleotides as U2
- **Annotation dependence:** Quality of input annotations affects extraction accuracy

---

## Classical Methods (Pre-intronIC)

### 1. The Burge/Levine/Durbin Method (U12Scan)

**Publication Context:** Early systematic U12 intron detection
**Core Strategy:** PWM-based log-odds scoring with statistical outlier detection

#### Algorithm Details

**Step 1: Model Training**
- Construct PWMs from known U12 AT-AC introns
- Primarily from published human U12 intron sets
- Define sequence properties of:
  - 5' splice site (5'ss)
  - Branch site

**Step 2: Log-Odds Scoring**
- Calculate normalized log-odds ratios:
  - **L^U12_5'ss**: Log(P(seq|U12) / P(seq|U2)) for 5' splice site
  - **L^U12_3'ss**: Log(P(seq|U12) / P(seq|U2)) for 3' splice site/BP region
- Compares probability of sequence matching U12 vs U2 consensus

**Step 3: Z-Score Normalization**
- Convert log-odds to z-scores (x, y):
  - **x = (L^U12_5'ss - μ_5') / σ_5'**
  - **y = (L^U12_BP - μ_BP) / σ_BP**
- Adjust for sample mean and standard deviation of known intron distributions
- Intron represented as vector **(x, y)** in score plane

**Step 4: Classification Criteria**

Two alternative criteria:

1. **Statistical outlier criterion:**
   - **t² = x² + y² > 20** (Mahalanobis-like distance)
   - Requires scores fall in **first quadrant** (both positive)
   - Identifies introns far from U2 cluster center

2. **Minimum threshold criterion:**
   - 5'ss and BP scores must **meet or exceed** minimum values found in reference U12 set
   - More conservative, requires both regions score above baseline

**Implementation:**
- Tools like **GeneID** modified to use positional weight matrices
- Separate matrices for U12 and U2 donor, branch, and acceptor sites

#### Strengths
- Clear statistical framework
- Geometric interpretation (distance in score space)
- Can identify outliers without hard thresholds

#### Limitations
- Relies on two-dimensional scoring (5'ss + BP only)
- Fixed t² threshold may not generalize across species
- Trained primarily on human data

---

### 2. The Sheth et al. Method (Score-Based Thresholds)

**Publication Context:** More recent approach with strict spatial constraints
**Core Strategy:** PWM scoring combined with biological constraints

#### Algorithm Details

**Step 1: Branch Point Sequence (BPS) Quality Assessment**

A BPS is deemed **"good"** if:
1. **Score > 65** for either:
   - U12-BPS-A9 PWM (Adenosine at position 9)
   - U12-BPS-A10 PWM (Adenosine at position 10)
2. **Located 8-30 nt upstream** of 3' splice site

**Step 2: Classification Rules for U12 GT-AG Introns**

An intron is classified as **U12-type GT-AG** if:

1. **Strong 5'ss criterion:**
   - U12_5'ss_score - U2_5'ss_score > **25 units**
   - Clear preference for U12 motif

2. **Moderate 5'ss + good BPS criterion:**
   - U12_5'ss_score - U2_5'ss_score > **10 units** AND
   - Possesses a "good" U12-type BPS (as defined above)

#### Key Features

**Distance Constraint (Critical):**
- BPS to 3'ss distance: **8-30 nt upstream**
- Sharply restricted compared to U2 (more variable)
- Based on biological requirement of U12 spliceosome

**Threshold Values:**
- BPS score threshold: **65**
- Strong 5'ss differential: **25**
- Moderate 5'ss differential: **10**

#### Strengths
- Incorporates biological distance constraints
- Two-tiered classification (strong vs moderate+good BP)
- Explicit handling of GT-AG U12 introns (hardest to classify)

#### Limitations
- Hard-coded thresholds may not generalize
- Requires species-specific calibration
- No explicit 3'ss scoring

---

## intronIC's Approach

### Core Methodology

**Classification Strategy:** Machine learning (SVM) with ensemble voting
**Feature Engineering:** Three-dimensional PWM-based z-scores
**Training Data:** Large curated reference sets (100k+ U2, 500+ U12)

### Detailed Implementation

#### 1. Feature Extraction via PWM Scoring

**Three Scoring Regions (vs two in classical methods):**

| Region | Default Coords | Biological Target |
|--------|----------------|-------------------|
| **5' splice site** | (-3, +9) from GT/AT | Donor site consensus |
| **Branch point** | (-55, -5) from AG | Branch point motif |
| **3' splice site** | (-6, +4) from AG | Acceptor site + PPT context |

**Code Location:** `seq_score()` (line 2114), `bp_score()` (line 2143)

**Scoring Algorithm:**
```
For each region:
  1. Apply PWM: Score = ∏(frequency[base_i][position_i])
  2. Add pseudocount to avoid zeros (default: 0.0001)
  3. For BP: sliding window to find optimal motif position
```

#### 2. Feature Normalization

**Method:** Z-score transformation via sklearn.StandardScaler
**Code Location:** `scale_scores()` (line 3655)

**Process:**
1. Fit scaler on **reference sequences only** (prevents data leakage)
2. Calculate: **z = (raw_score - μ) / σ**
3. Apply same transformation to test set
4. Results in three features per intron: **(five_z, bp_z, three_z)**

**Key Difference from U12Scan:**
- Three dimensions vs two
- Includes explicit 3'ss scoring
- Uses sklearn for standardization

#### 3. SVM Classification

**Model:** Linear SVM (sklearn.svm.SVC)
**Code Location:** `train_svm()` (line 5345), `optimize_svm()` (line 5431)

**Hyperparameter Optimization:**
- **Parameter:** C (soft-margin penalty)
- **Range:** 10^-6 to 10^6 (geometric spacing)
- **Method:** 5-round grid search with refinement
  - Round 1: Coarse grid (13 points)
  - Rounds 2-5: Geometric refinement around best value
- **Cross-validation:** 5-fold stratified
- **Metric:** balanced_accuracy (handles class imbalance)

**Training Process:**
1. Split: 80/20 training/test
2. Class weights: Balanced (compensates for U12 rarity)
3. Optional: Multiple models via U2 subsampling (ensemble)
4. Track metrics: F1 score, Precision-Recall AUC

#### 4. Classification & Thresholding

**Default Threshold:** 90% U12 probability
**Code Location:** Line 964 (argument default)

**Classification Process:**
1. SVM outputs: `predict_proba()` → P(U12) ∈ [0, 1]
2. Convert to percentage: **svm_score = P(U12) * 100**
3. Apply threshold: **type_id = 'u12' if svm_score > 90 else 'u2'**
4. Calculate relative score: **relative_score = svm_score - threshold**

**Ensemble Aggregation** (if multiple models):
- Average probabilities across models
- Weight by F1 scores (better models weighted higher)
- Compute mean distance to decision boundary

#### 5. Advanced Features

**Recursive Training** (`recursive_scoring()`, line 3890):
- Trigger: ≥5 confident U12s found
- Build species-specific PWMs from confident U12s
- Re-train SVM on species-specific data
- Re-score all introns
- Benefit: Better for distant species

**Dynamic U2 BP Matrix Generation** (`build_u2_bp_matrix()`, line 2775):
- Identify likely U2 introns (low U12 scores)
- Extract branch point regions
- Build species-specific U2 BP PWM
- Replace default human PWM

**Non-Canonical Boundary Adjustment** (`u12_correction()`, line 2298):
- Scan ±200bp for strong U12-like motifs
- Adjust boundaries if better canonical sites found
- Tag corrected introns with [c:distance]

---

## Comparative Analysis

### Feature Comparison Table

| Feature | U12Scan (Classical) | Sheth et al. | intronIC |
|---------|---------------------|--------------|----------|
| **Scoring Method** | PWM log-odds | PWM raw scores | PWM frequencies → z-scores |
| **Number of Features** | 2 (5'ss, BP) | 2 (5'ss, BP) | **3 (5'ss, BP, 3'ss)** |
| **Classifier** | Statistical outlier (t²>20) | Hard thresholds | **SVM (ML-based)** |
| **Normalization** | Z-scores (manual) | None | **Z-scores (sklearn)** |
| **Threshold Type** | Fixed geometric | Fixed score diffs | **Learned + tunable (default 90%)** |
| **Distance Constraints** | Implicit | **Explicit (8-30 nt)** | Implicit in BP region (-55,-5) |
| **Handles Class Imbalance** | No | No | **Yes (balanced weights)** |
| **Species Adaptation** | No | No | **Yes (recursive training)** |
| **Ensemble Methods** | No | No | **Yes (optional)** |
| **3'ss Scoring** | Limited | No | **Explicit** |
| **Training Data Size** | Small | Small | **Large (100k U2, 500 U12)** |
| **Hyperparameter Tuning** | No | No | **Yes (5-round grid search)** |

### Key Improvements in intronIC

1. **Three-Dimensional Feature Space**
   - Adds explicit 3' splice site scoring
   - More information for classification
   - May capture polypyrimidine tract differences

2. **Machine Learning Framework**
   - Learns decision boundary from data
   - Not dependent on hard-coded thresholds
   - Better handles continuous variation

3. **Automatic Hyperparameter Optimization**
   - Finds optimal C parameter via cross-validation
   - Adapts to specific datasets
   - Reduces manual tuning

4. **Species Adaptability**
   - Recursive training for distant species
   - Dynamic matrix generation
   - Not limited to human-trained models

5. **Class Imbalance Handling**
   - Balanced class weights in SVM
   - Balanced accuracy metric
   - Prevents U2 overwhelming the classifier

6. **Ensemble Voting**
   - Multiple models with different U2 subsamples
   - F1-weighted averaging
   - More robust predictions

### Potential Gaps vs Classical Methods

1. **Distance Constraints Not Explicit**
   - Classical: Sheth et al. enforces 8-30 nt BP-to-3'ss distance
   - intronIC: Uses (-55, -5) search window but doesn't enforce found BP position
   - **Impact:** May miss a valuable filtering criterion
   - **Recommendation:** Add explicit distance check in classification

2. **Hard Threshold Options**
   - Classical: Multiple threshold levels (25, 10 for 5'ss differential)
   - intronIC: Single probability threshold (default 90%)
   - **Impact:** Less interpretability for borderline cases
   - **Recommendation:** Could add score differential reporting

3. **Polypyrimidine Tract (PPT) Analysis**
   - Classical: Explicitly looks for PPT absence in U12
   - intronIC: Implicit via 3'ss scoring
   - **Impact:** May not fully capture PPT differences
   - **Recommendation:** Could add explicit PPT scoring/detection

---

## Recommendations for Development

### High Priority (Immediate Impact)

#### 1. Implement Explicit Distance Constraints

**Current State:**
- BP search region: (-55, -5) from 3'ss
- Optimal BP position found via sliding window
- **No validation** that found BP is 8-30 nt from 3'ss

**Recommendation:**
```python
# In classification stage, add BP distance check
def validate_bp_distance(intron):
    """
    Verify branch point is 8-30 nt upstream of 3' splice site.
    Classical constraint from Sheth et al. method.
    """
    bp_position = intron.bp_stop  # position relative to 3'ss
    distance_from_3ss = abs(bp_position)

    if 8 <= distance_from_3ss <= 30:
        return True
    else:
        return False

# Add to classification criteria:
# High confidence U12: svm_score > 90 AND valid_bp_distance
# Medium confidence U12: svm_score > 90 OR valid_bp_distance
```

**Benefits:**
- Leverages known biological constraint
- May reduce false positives
- Provides additional confidence metric

**Location:** Add to `model_score()` or classification section

#### 2. Add Score Differential Reporting

**Current State:**
- Only outputs final SVM probability
- No U12 vs U2 score comparison

**Recommendation:**
```python
# Add to Intron class or output
def calculate_score_differentials(intron):
    """
    Calculate U12 vs U2 score differentials for each region.
    Useful for interpreting borderline classifications.
    """
    return {
        'five_ss_diff': intron.u12_five_score - intron.u2_five_score,
        'bp_diff': intron.u12_bp_score - intron.u2_bp_score,
        'three_ss_diff': intron.u12_three_score - intron.u2_three_score
    }
```

**Benefits:**
- Better interpretability
- Can compare to classical thresholds (25, 10)
- Useful for debugging borderline cases

**Location:** Add to output generation, `.scores.iic` file

#### 3. Implement Polypyrimidine Tract (PPT) Detection

**Current State:**
- 3'ss region scored but PPT not explicitly identified
- U12 introns lack strong PPT

**Recommendation:**
```python
def detect_polypyrimidine_tract(intron, window=(-40, -10)):
    """
    Calculate pyrimidine (T/C) content in region upstream of 3'ss.
    U12 introns typically have weaker PPT than U2.
    """
    region_seq = intron.get_seq(intron.seq, *window)
    pyrimidine_count = region_seq.count('T') + region_seq.count('C')
    ppt_strength = pyrimidine_count / len(region_seq)

    return ppt_strength

# Add as feature or filter:
# Strong PPT (>0.8) suggests U2
# Weak PPT (<0.6) consistent with U12
```

**Benefits:**
- Explicitly captures key U12 feature
- Could be added as 4th SVM feature
- Useful for filtering

**Location:** Add to feature extraction, optionally to SVM input

### Medium Priority (Enhancement)

#### 4. Multi-Level Confidence Classification

**Current State:**
- Binary classification at single threshold (90%)

**Recommendation:**
```python
# Implement confidence tiers
def classify_with_confidence(intron):
    score = intron.svm_score
    bp_valid = validate_bp_distance(intron)

    if score >= 95 and bp_valid:
        return 'U12_high_confidence'
    elif score >= 90:
        return 'U12_medium_confidence'
    elif score >= 75 and bp_valid:
        return 'U12_borderline'
    elif score < 25:
        return 'U2_high_confidence'
    else:
        return 'U2_medium_confidence'
```

**Benefits:**
- Provides uncertainty quantification
- Users can choose confidence level
- Better for downstream filtering

#### 5. AT-AC Specific Detection

**Current State:**
- Treats all terminal dinucleotides similarly
- AT-AC is strongest U12 signal

**Recommendation:**
```python
# Add AT-AC pre-filter or boosting
def boost_atac_score(intron):
    """
    AT-AC introns are almost exclusively U12.
    Could use as strong prior or separate classification.
    """
    if intron.dnts == 'AT-AC':
        # Option 1: Automatic U12 classification
        # Option 2: Boost SVM score
        intron.svm_score = min(100, intron.svm_score * 1.5)
        intron.dynamic_tag.add('AT-AC_boosted')
```

**Benefits:**
- Leverages strongest U12 signal
- May improve AT-AC recall
- Aligns with biological knowledge

#### 6. Implement U12Scan-Style Geometric Classifier

**Current State:**
- Only SVM classifier

**Recommendation:**
```python
# Add as alternative or ensemble member
def u12scan_classifier(intron, threshold_t2=20):
    """
    Classical U12Scan method: t² = x² + y² > 20
    Requires first quadrant (both positive z-scores).
    """
    x = intron.five_z_score
    y = intron.bp_z_score

    t_squared = x**2 + y**2
    first_quadrant = (x > 0 and y > 0)

    if t_squared > threshold_t2 and first_quadrant:
        return 'U12'
    else:
        return 'U2'

# Could ensemble with SVM:
# U12 if (SVM says U12) OR (U12Scan says U12)
```

**Benefits:**
- Provides comparison to classical method
- May catch different edge cases
- Ensemble could improve recall

### Low Priority (Future Research)

#### 7. Deep Learning Alternative

**Current State:**
- Linear SVM (good but limited)

**Recommendation:**
- Explore CNN or transformer models
- Could learn sequence patterns directly
- May capture long-range dependencies
- Requires larger training set

#### 8. Active Learning Pipeline

**Current State:**
- Fixed training set

**Recommendation:**
- Iteratively add high-confidence predictions to training
- Human-in-the-loop for borderline cases
- Gradually improve species-specific performance

#### 9. Motif Discovery

**Current State:**
- Uses pre-defined PWMs

**Recommendation:**
- Automatically discover motifs from data
- Tools: MEME, DREME
- May find novel U12 subtypes

---

## Field Consensus Sequences Reference

### U2-type (Major Introns)

| Feature | Consensus | Notes |
|---------|-----------|-------|
| **5' SS** | KAG GTRAGT | K=G/T, R=A/G; variable |
| **Branch Point** | CURACU | Degenerate; R=A/G, U=T |
| **PPT** | YYYYYYYYYY | Y=C/T; strong pyrimidine tract |
| **3' SS** | YAG | Preceded by PPT |
| **Terminal Di-nts** | GT-AG | ~99.6% of U2 introns |
| **BP-to-3'SS Distance** | Variable | Often 18-40 nt |

### U12-type (Minor Introns)

| Feature | Consensus | Notes |
|---------|-----------|-------|
| **5' SS** | RTATCCTT | R=A/G; highly conserved |
| **Branch Point** | TTCCTTRACYCY or UUUCCUUAACYCY | Highly conserved; R=A/G, Y=C/T |
| **BPS Adenosine** | Position 9 or 10 | Critical for branch formation |
| **PPT** | Absent or weak | Key distinguishing feature |
| **3' SS** | Goo | G followed by weak sequence; no strong PPT |
| **Terminal Di-nts** | GT-AG (70%), AT-AC (20%), others (10%) | More diverse than U2 |
| **BP-to-3'SS Distance** | **8-30 nt** (strict) | **Most diagnostic feature** |

### Rare U12 Terminal Dinucleotides

Computational methods have confirmed surprising diversity:
- **AT-AG** (relatively common U12 variant)
- **GT-AT** (rare)
- **AT-AT** (very rare)
- **GT-GG** (very rare)
- **AT-AA** (very rare)
- **GT-AA** (very rare)

**Implication:** Cannot rely solely on terminal dinucleotides for classification.

---

## Verification & Filtering Methods

### Classical Approaches

#### EST/cDNA Confirmation
- **Purpose:** Validate predicted introns with expressed sequences
- **Tools:** SSAHA, BLAST
- **Method:** Search for sequences spanning predicted splice junction
- **Benefit:** Confirms introns are actually spliced in vivo

#### RNA-seq Filtering
- **Purpose:** Detect splice junctions from RNA-seq data
- **Tools:** MapSplice, STAR, HISAT2
- **Challenges:**
  - Alignment errors create false junctions
  - Repetitive regions (filter with Repbase)
  - Low complexity sequences (filter with DUST)
  - SNPs/indels create false non-canonical sites

**U2/U12-like Classification from RNA-seq:**
- Apply PWMs to non-canonical junctions detected
- Junction is U2/U12-like if:
  - Score > 70.00, OR
  - Score > 95th percentile of scrambled sequences AND shares splice site with canonical junction

### intronIC's Approach

**Current Filtering:**
- Length filter (default: ≥30 bp)
- Ambiguous sequence filter (non-ACTG)
- Non-canonical dinucleotide handling
- Duplicate detection
- Isoform selection (longest by default)

**Not Implemented:**
- EST/cDNA confirmation
- RNA-seq junction validation
- SNP/indel filtering

**Recommendation:** These are valuable for genome-wide scans of unannotated data but less critical for annotated genomes where introns are already known.

---

## Summary & Strategic Recommendations

### intronIC's Position in the Field

**Strengths:**
1. Most sophisticated ML-based approach
2. Largest training dataset
3. Best handling of class imbalance
4. Species adaptability via recursive training
5. Comprehensive metadata and output

**Opportunities for Enhancement:**
1. Incorporate explicit biological constraints (BP distance)
2. Add classical threshold reporting for interpretability
3. Explicit PPT detection
4. Multi-level confidence classification

### Development Roadmap

**Phase 1: Quick Wins (1-2 weeks)**
- Add BP distance validation (8-30 nt check)
- Report score differentials in `.scores.iic`
- Add confidence tiers to output

**Phase 2: Feature Enhancement (1 month)**
- Implement PPT detection
- Add U12Scan-style geometric classifier as ensemble option
- AT-AC specific handling

**Phase 3: Research & Validation (2-3 months)**
- Benchmark against U12Scan and Sheth methods on published datasets
- Test on diverse species
- Validate distance constraints improve accuracy
- Publish comparison study

**Phase 4: Advanced Methods (future)**
- Deep learning models
- Active learning pipeline
- Automated motif discovery

---

## References & Further Reading

### Key Publications

1. **Burge, Levine, Durbin** - Original U12Scan methodology
2. **Sheth et al.** - Threshold-based classification with distance constraints
3. **Moyer et al. (2020)** - intronIC publication
   - "Comprehensive database and evolutionary dynamics of U12-type introns"
   - *Nucleic Acids Research*, 48(13):7066–7078
   - https://doi.org/10.1093/nar/gkaa464

### Biological Background

- **U12 spliceosome machinery:** Distinct snRNPs (U11, U12, U4atac, U6atac vs U1, U2, U4, U6)
- **Evolutionary conservation:** U12 introns ancient but rare
- **Functional importance:** Often in genes critical for development, neurogenesis

### Computational Resources

- **U12DB:** Database of confirmed U12 introns
- **SpliceRack:** Alternative splicing database
- **ERISdb:** Database of experimentally validated introns
- **Ensembl/GENCODE:** High-quality genome annotations

---

## Appendix: Code Implementation Notes

### Adding BP Distance Constraint

**File:** `intronIC/intronIC.py`
**Location:** After SVM scoring, before final classification

```python
# Around line 5800 in parallel_svm_score or in model_score
def add_bp_distance_feature(intron):
    """
    Calculate branch point distance from 3' splice site.
    Valid range: 8-30 nt (Sheth et al.)
    """
    # intron.bp_stop is position of BP end relative to 3'ss
    # Negative values indicate upstream
    bp_distance = abs(intron.bp_stop) - abs(intron.bp_start)
    intron.bp_distance_from_3ss = bp_distance
    intron.bp_distance_valid = (8 <= bp_distance <= 30)

    return intron
```

### Adding Score Differential Output

**File:** `intronIC/intronIC.py`
**Location:** In output_format() function around line 4200

```python
# Modify .scores.iic output to include differentials
score_diff_header = [
    'five_ss_differential',
    'bp_differential',
    'three_ss_differential'
]

# Calculate during scoring phase:
intron.five_diff = intron.five_raw_score_u12 - intron.five_raw_score_u2
intron.bp_diff = intron.bp_raw_score_u12 - intron.bp_raw_score_u2
intron.three_diff = intron.three_raw_score_u12 - intron.three_raw_score_u2
```

### PPT Strength Calculation

**File:** `intronIC/intronIC.py`
**Location:** In get_sub_seqs() function around line 2650

```python
def calculate_ppt_strength(intron, window=(-40, -10)):
    """
    Calculate polypyrimidine tract strength.
    U12 introns typically have weak PPT (<0.6).
    U2 introns typically have strong PPT (>0.7).
    """
    region = intron.get_seq(intron.seq, *window)
    if len(region) == 0:
        return 0.0

    pyrimidines = region.count('C') + region.count('T')
    ppt_strength = pyrimidines / len(region)

    intron.ppt_strength = ppt_strength
    intron.ppt_weak = (ppt_strength < 0.6)  # Suggests U12

    return intron
```

---

**Document Version:** 1.0
**Maintained by:** intronIC Development Team
**Last Review:** 2025-11-01
