# Porting Checklist - Prevent Missing Functionality

When porting a feature from original intronIC, complete ALL steps:

## 1. Identify Primary Function
- [ ] Located main function (e.g., `bp_score()`)
- [ ] Read docstring and understand purpose
- [ ] Note line numbers and file location

## 2. Find All Dependencies
```bash
# Find nested helpers
grep -B5 "def main_function" original.py | grep "def "

# Find callers
grep -n "main_function(" original.py

# Find related conditionals  
grep -n "if.*length\|if.*short" original.py | grep -i "feature_name"
```

- [ ] Listed all nested helper functions
- [ ] Listed all functions that call this function
- [ ] Listed all conditionals that affect this function

## 3. Trace Data Flow
- [ ] Input: Where does data come from?
- [ ] Preprocessing: Any transformations before main function?
- [ ] Main processing: What does function do?
- [ ] Postprocessing: Any transformations after?
- [ ] Output: Where does result go?

## 4. Check for Edge Cases
Search for patterns:
- [ ] `if.*< ` (minimum bounds)
- [ ] `if.*> ` (maximum bounds)
- [ ] `if.*== ` (exact values)
- [ ] `if.*None` (null checks)
- [ ] `try/except` (error handling)
- [ ] Short/long special cases

## 5. Write Tests BEFORE Porting
- [ ] Integration test (end-to-end with real data)
- [ ] Boundary tests (min, max, edge cases)
- [ ] Error tests (invalid inputs)
- [ ] Performance test (if applicable)

## 6. Port Implementation
- [ ] Port main function
- [ ] Port all helpers (including nested)
- [ ] Port all preprocessing logic
- [ ] Port all postprocessing logic
- [ ] Add type hints
- [ ] Add docstrings with "Port from:" references

## 7. Verification
- [ ] All new tests pass
- [ ] Run on test data from original
- [ ] Compare outputs (bit-exact if possible)
- [ ] Check for warnings/errors in logs
- [ ] Profile performance (shouldn't be slower)

## 8. Documentation
- [ ] Update module docstring
- [ ] Add usage examples
- [ ] Document any differences from original
- [ ] Note any TODOs or known limitations

## Example: Branch Point Scoring

### What We Missed:
- ❌ Step 2: Didn't find `_short_bp_adjust()` nested helper
- ❌ Step 3: Didn't trace extraction → scoring data flow  
- ❌ Step 4: Didn't check for length-based edge cases
- ❌ Step 5: No integration test for short introns

### How Checklist Would Prevent:
```bash
# Step 2: Find nested helpers
$ grep -A50 "def assign_seqs" intronIC.py | grep "    def "
    def _short_bp_adjust(intron, bpc, fsl):  # ← Would catch this!

# Step 4: Check edge cases
$ grep -n "if intron.length" intronIC.py | grep bp
2610:    if intron.length < abs(bp_coords[0]) + five_score_length:
         # ← Would catch this!

# Step 5: Integration test
def test_short_intron_bp_scoring():
    """Test 39bp intron scores with clamped window."""
    # Would fail until _short_bp_adjust logic is ported
```

## Lessons Learned

**Nested functions are dangerous:** Easy to miss when scanning top-level functions only.

**Separation of concerns can hide dependencies:** Just because scoring and extraction are separate doesn't mean they don't interact.

**Edge cases reveal gaps:** The first time you run real data, edge cases appear.

**Integration tests > Unit tests for completeness:** Unit tests verify what you implemented; integration tests verify you didn't miss anything.
