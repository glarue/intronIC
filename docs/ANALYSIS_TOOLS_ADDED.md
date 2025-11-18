# Post-Hoc Analysis Tools - Documentation Complete

**Date:** 2025-11-16
**Purpose:** Document the post-hoc prediction analysis methodology for future use

---

## What Was Created

### 1. Main Analysis Tool (`tools/analyze_predictions.py`)

**Purpose:** Diagnose false positives by manually computing features from z-scores

**Features:**
- ✅ Parses `.score_info.iic` files
- ✅ Manually computes ALL BothEndsStrong features (even those not in model)
- ✅ Identifies problematic patterns (negative scores, one-end-strong, imbalance)
- ✅ Shows both summary tables and detailed analysis
- ✅ Configurable filtering (score range, top N, extremes)

**Usage Examples:**
```bash
# Basic - analyze top 10 predictions
pixi run python tools/analyze_predictions.py species.score_info.iic

# Quick summary only
pixi run python tools/analyze_predictions.py species.score_info.iic --summary

# Top 5 with detailed analysis
pixi run python tools/analyze_predictions.py species.score_info.iic -n 5

# High and low scoring predictions
pixi run python tools/analyze_predictions.py species.score_info.iic --show-extremes

# Filter by score range
pixi run python tools/analyze_predictions.py species.score_info.iic --min-score 5.0
```

**Output - Summary Table:**
```
#     RelScore   SVM%      s5     sBP      s3  min_all Pattern
1        10.00  100.0    0.45    1.37   -1.90    -1.90 NEGATIVE scores: s3
2        10.00  100.0    0.28    1.12    0.41     0.28 ONE-END-STRONG patte
3         7.96   98.0   -0.15    0.30   11.45    -0.15 ONE-END-STRONG patte
```

**Output - Detailed Analysis:**
```
#1: CaeEle-WBGene00021036@W05F2.4a.1-intron_8(35)
  Relative score: 10.00
  SVM probability: 100.0%

  Base z-scores:
    s5:   0.4450
    sBP:  1.3708
    s3:  -1.9010

  Min features (require all/both strong):
    min_all:   -1.9010  ← ALL THREE must be strong
    min_5_bp:   0.4450
    min_5_3:   -1.9010
    min_bp_3:  -1.9010

  Max features (at least one strong - can reward FPs):
    max_all:    1.3708
    max_5_bp:   1.3708
    max_5_3:    0.4450
    max_bp_3:   1.3708

  Imbalance penalties (should be small for real U12s):
    neg_absdiff_5_bp: -0.9258  ← 5'/BP balance
    neg_absdiff_5_3:  -2.3460  ← 5'/3' balance
    neg_absdiff_bp_3: -3.2718  ← BP/3' balance

  ⚠ Pattern warnings:
    • NEGATIVE scores: s3
    • ONE-END-STRONG pattern detected
    • High 5'/3' imbalance (|diff|=2.35)
    • High BP/3' imbalance (|diff|=3.27)
    • min_all is NEGATIVE (should strongly penalize)
```

### 2. Model Inspector (`tools/inspect_model.py`)

**Purpose:** Inspect trained model coefficients and configuration

**Features:**
- ✅ Loads pickled models (handles compression)
- ✅ Shows pipeline structure and transformer config
- ✅ Displays all feature coefficients
- ✅ Shows all models in ensemble

**Usage:**
```bash
# Default location
pixi run python tools/inspect_model.py

# Specific model
pixi run python tools/inspect_model.py path/to/model.pkl
```

**Output:**
```
Transformer: <class 'classification.transformers.BothEndsStrongTransformer'>
  include_max: True
  Features (11): ['s5', 'sBP', 's3', 'min_5_bp', 'min_5_3', 'min_all', ...]

Coefficients (11 features):
  s5                  : +1.249035
  sBP                 : +0.000000  ← ZEROED!
  s3                  : +0.000000  ← ZEROED!
  min_5_bp            : +0.527692
  min_5_3             : +0.171555
  min_all             : +0.000000  ← ZEROED! (Critical feature)
  neg_absdiff_5_bp    : +0.000000  ← ZEROED!
  neg_absdiff_5_3     : +0.000000  ← ZEROED!
  neg_absdiff_bp_3    : +0.000000  ← ZEROED!
  max_5_bp            : +0.000000
  max_5_3             : +0.396311  ← REWARDING ONE-END-STRONG!
```

### 3. Comprehensive Documentation (`docs/POST_HOC_ANALYSIS.md`)

**Contents:**
- ✅ When to use post-hoc analysis (4 scenarios)
- ✅ How to interpret output (summary tables, detailed analysis)
- ✅ Feature computation details (all 14 possible features)
- ✅ Common diagnosis workflows (3 detailed workflows)
- ✅ Case study: L1 penalty investigation (complete example)
- ✅ Advanced: Manual feature contribution calculation
- ✅ Tips and best practices

**Key Sections:**

**When to Use:**
1. Investigating false positives
2. Comparing model versions
3. Understanding training issues
4. Species-specific issues

**Workflows:**
1. Why are there false positives?
2. Did my fix work?
3. Understanding species differences

**Case Study:**
- Complete walkthrough of the L1 penalty issue
- From symptom (FPs) to diagnosis (zeroed features) to fix (L2 + remove redundancy)

### 4. Tools Directory README (`tools/README.md`)

**Contents:**
- ✅ Quick reference for both tools
- ✅ Usage examples with actual commands
- ✅ When to use each tool
- ✅ Red flags to watch for
- ✅ Typical workflow for diagnosing FPs
- ✅ Complete example: L1 penalty investigation
- ✅ Dependencies listed

---

## File Organization

```
intronIC/
├── tools/
│   ├── README.md                    # Quick reference for tools
│   ├── analyze_predictions.py      # Main analysis tool (executable)
│   └── inspect_model.py             # Model coefficient inspector (executable)
├── docs/
│   └── POST_HOC_ANALYSIS.md        # Comprehensive methodology guide
├── FP_ROOT_CAUSE_ANALYSIS.md        # Case study: L1 penalty issue
└── ANALYSIS_TOOLS_ADDED.md          # This file
```

---

## Key Features Documented

### Feature Computation

The tool computes ALL 14 possible BothEndsStrong features:

**Base (3):**
- s5, sBP, s3

**Min features (4):**
- min_5_bp, min_5_3, min_bp_3, min_all

**Max features (4):**
- max_5_bp, max_5_3, max_bp_3, max_all

**Imbalance penalties (3):**
- neg_absdiff_5_bp, neg_absdiff_5_3, neg_absdiff_bp_3

**Note:** Model only uses subset (e.g., 7 features after our fix)

### Pattern Detection

The tool automatically identifies:

1. **Negative scores** - Any base score < 0
2. **One-end-strong** - Exactly one score > 1.0, others weak
3. **High imbalance** - |difference| > 1.0 between any pair
4. **Negative min_all** - Weakest signal is negative (should reject)
5. **All weak** - All scores < 0.5 (weak but balanced)

### Interpretation Guidance

**Red flags (likely FPs):**
- ⚠ min_all < 0 → At least one signal weak/negative
- ⚠ Large imbalance (>2-3) → Inconsistent signals
- ⚠ One-end-strong pattern → One motif matches by chance
- ⚠ Negative base scores → Below background level

**Green flags (likely real U12s):**
- ✓ All z-scores > 1.0 → All signals strong
- ✓ min_all > 1.0 → Consistent strength
- ✓ Small imbalances (<0.5) → Balanced signals
- ✓ No warnings → Looks legitimate

---

## Example Use Cases (Documented)

### 1. Investigating False Positives

```bash
# Check predictions
pixi run python tools/analyze_predictions.py species.score_info.iic --summary

# Count warning types
pixi run python tools/analyze_predictions.py species.score_info.iic --summary | \
    grep -oE "(NEGATIVE|ONE-END|imbalance|min_all)" | sort | uniq -c

# Inspect model
pixi run python tools/inspect_model.py
```

### 2. Comparing Before/After Fix

```bash
# Before
pixi run python tools/analyze_predictions.py old.score_info.iic --summary > before.txt

# After
pixi run python tools/analyze_predictions.py new.score_info.iic --summary > after.txt

# Compare
diff before.txt after.txt
```

### 3. Species Comparison

```bash
for species in human mouse celegans drosophila; do
    echo "=== $species ==="
    pixi run python tools/analyze_predictions.py ${species}.score_info.iic --summary -n 3
done
```

---

## What Makes This Reusable

1. **Command-line tool** - Easy to run on any score_info.iic file
2. **Configurable options** - Adjust to different needs (--summary, -n, --min-score)
3. **Clear output** - Both human-readable tables and parseable formats
4. **Well documented** - Comprehensive guide with examples
5. **Methodology explained** - Not just "what" but "why" and "when"
6. **Case study included** - Real example of using the tools to solve a problem
7. **Executable scripts** - chmod +x, ready to use

---

## Future Use

**Next time you see false positives:**

1. Run: `pixi run python tools/analyze_predictions.py species.score_info.iic --summary`
2. Look for patterns in "Pattern" column
3. Run: `pixi run python tools/inspect_model.py`
4. Check if problematic patterns correspond to missing/zeroed features
5. Refer to: `docs/POST_HOC_ANALYSIS.md` for interpretation
6. Fix model configuration based on findings

**The tools will:**
- ✅ Show you exactly what features are being computed
- ✅ Identify which patterns are problematic
- ✅ Help you diagnose whether issue is in model or data
- ✅ Guide you toward the right fix

---

## Testing

Both tools tested and working:

```bash
$ pixi run python tools/analyze_predictions.py \
    run_tests/caenorhabditis_elegans.pretrained.score_info.iic -n 3 --summary

Found 109830 total introns
After filtering: 8 introns

#     RelScore   SVM%      s5     sBP      s3  min_all Pattern
1        10.00  100.0    0.45    1.37   -1.90    -1.90 NEGATIVE scores: s3
2        10.00  100.0    0.28    1.12    0.41     0.28 ONE-END-STRONG patte
3        10.00  100.0    0.39    0.10    0.74     0.10 OK
```

✓ Works correctly!

---

## Documentation Quality

✅ **Comprehensive:** Covers when, why, how, and what
✅ **Practical:** Real examples and commands
✅ **Educational:** Explains methodology and interpretation
✅ **Actionable:** Clear workflows and decision trees
✅ **Referenced:** Case study showing actual use
✅ **Maintained:** All in version control, easy to update

---

**Created by:** Claude Code
**Date:** 2025-11-16
**Purpose:** Document post-hoc analysis approach for future FP investigations
**Status:** ✅ Complete and tested
