# TODO Items from Systematic Review
**Date:** 2025-11-04
**Source:** Applied PORTING_CHECKLIST.md to all ported modules

## IMMEDIATE (Next Session)

### 1. Add PWM Matrix Name Version Tag Support
**Priority:** LOW  
**Module:** `scoring/pwm.py`  
**Function:** `PWMLoader._parse_matrix_name()`  
**Line:** 319

**What's missing:**
```python
import re

# Add after line 357, before return statement:
# Extract optional version tag (e.g., "v2" from "u12_atac_bp_v2")
matrix_version = re.findall(r'[^A-Za-z]v\.?([^_\s]+)', matrix_name)
if matrix_version:
    name_bits.append(matrix_version[0])
```

**Why needed:**
- Allows users to provide versioned custom matrices
- Original intronIC supports this
- Currently no default matrices use versions, so low impact

**Test to add:**
```python
def test_parse_matrix_name_with_version():
    """Test parsing matrix names with version tags."""
    assert PWMLoader._parse_matrix_name("u12_atac_five_v2") == ('u12', 'atac', 'five', 'v2')
    assert PWMLoader._parse_matrix_name("u2_gtag_bp_v3.1") == ('u2', 'gtag', 'bp', '3.1')
```

**Estimated effort:** 15 minutes

---

## VERIFIED AS COMPLETE ✓

These were double-checked and confirmed as correctly ported:

1. ✅ Branch point window clamping (`_short_bp_adjust`)
2. ✅ Sliding window for exon pairs (`_window`)
3. ✅ Overlap detection
4. ✅ PWM core scoring logic
5. ✅ Pseudocount handling
6. ✅ Ambiguous base handling
7. ✅ Ignore positions support
8. ✅ Branch point search algorithm

---

## NOT APPLICABLE (Future Features)

These nested functions belong to features not yet ported:

1. `_shift_phase` - Part of u12_correction (boundary adjustment)
2. `_correct_exon` - Part of correct_annotation (GFF modification)
3. `_change_coords` - Part of correct_annotation
4. `_iter_bps` - Part of build_u2_bp_matrix (dynamic PWM building)
5. `_get_attrs` - Part of assign_clusters (advanced output)

**Action:** None needed now. Will port when implementing those features.

---

## LESSONS LEARNED

### What Worked
- ✅ Systematic grep for nested functions caught the issue
- ✅ Checking edge case patterns (length < X, if None) useful
- ✅ Cross-referencing line numbers between original and ported code

### What To Improve
- Need better initial dependency mapping before porting
- Should check for ALL nested helpers at start, not just main function
- Integration tests would have caught BP window issue earlier

### Process Update
**Always run before considering module "complete":**
```bash
# 1. Find all nested helpers
grep "^    def " original.py | grep -B5 -A2 "main_function"

# 2. Check for edge cases
grep -n "if.*length\|if.*None\|try:" original.py | grep "feature"

# 3. Verify all are ported
# 4. Add integration test
# 5. Run on real data
```
