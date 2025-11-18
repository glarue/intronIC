# Major Concerns with Existing intronIC Codebase

**Date:** 2025-11-01
**Purpose:** Document critical issues to address during refactoring

---

## Critical Issues (Must Address)

### 1. ⚠️ Coordinate System Complexity - HIGH RISK

**Problem:**
The code juggles multiple coordinate systems with manual conversions:

- **BED input:** 0-based, half-open intervals `[start, stop)`
- **Internal:** 1-based, closed intervals `[start, stop]`
- **Output (BED):** Convert back to 0-based

**Evidence:**
```python
# Line 4187: BED input conversion
start += 1  # BED files are 0-indexed for the start coord

# Line 4277: BED output conversion
attribs[1] -= 1  # revert start to 0-based indexing
```

**Recent Bugs:**
- `99b5aaa BUGFIX: correct off-by-one error in BED output`
- `ca0429b Correctly avoid annot coord adj. for BED input`

**Why This Is Dangerous:**
- Off-by-one errors are subtle and hard to detect
- Can produce introns with wrong boundaries (affects classification)
- Error may not be obvious in output (off by 1 bp looks "close enough")
- Affects every downstream analysis

**Refactoring Priority:** 🔴 **HIGHEST**

**Solution:**
- Create explicit `GenomicCoordinate` class that encodes coordinate system
- Never mix coordinate systems
- All conversions happen at I/O boundaries only
- Internal representation is always consistent (1-based)
- Add validation that start < stop

---

### 2. ⚠️ State Mutation in u12_correction() - HIGH RISK

**Problem:**
The `u12_correction()` function directly modifies intron coordinates in-place:

```python
# Line 2350-2351
intron.start += shift
intron.stop += shift
```

**Known Issue (from TODO comment):**
```python
# Line 2309
#TODO make u12 correction modify a *copy* of the intron object,
# so that the result can be checked for e.g. NC bounds or ambiguous
# sequence. It may be desirable to revert certain cases e.g. intron is
# shifted but the resulting intron is omitted (b/c of ambiguous sequence,
# short, etc) and thus, the corrected boundaries aren't propagated to
# the annotation.iic file)
```

**Why This Is Dangerous:**
- Modifying state makes code hard to reason about
- Can't easily rollback if correction causes other issues
- Strand-dependent logic (line 2348): `if intron.strand == '-': shift *= -1`
- Phase adjustment also happens (line 2346)
- Multiple inter-dependent state changes

**Observed Complexity:**
```python
# Different adjustment for negative strand
if intron.strand == '-':
    shift *= -1
intron.start += shift
intron.stop += shift
```

**Refactoring Priority:** 🔴 **HIGH**

**Solution:**
- Make `u12_correction()` return a new Intron object (immutable approach)
- OR return a validated correction that can be applied/rejected
- Validate corrected intron before accepting changes
- Make coordinate changes atomic (all or nothing)

---

### 3. ⚠️ No Validation of Biological Constraints

**Problem:**
Code doesn't validate biological constraints during processing:

- No check that BP is 8-30 nt from 3'ss (field-proven constraint)
- No validation that intron length is reasonable (>30 bp minimum exists, but no maximum)
- No check that sequences are valid (beyond "contains N")
- No validation of PWM scores (could be NaN, infinite, etc.)

**Current Validation (Line 4174-4182):**
```python
try:
    loc, start, stop, name, _, strand = l.strip().split('\t')[:6]
except:
    continue  # Silently skip - no error message!
try:
    start, stop = sorted(map(int, [start, stop]))
except:
    continue  # Silently skip - no error message!
```

**Why This Is Dangerous:**
- Invalid data propagates through pipeline
- Errors manifest late (hard to debug)
- Silent failures hide problems
- No way to know if input is malformed

**Refactoring Priority:** 🟡 **MEDIUM-HIGH**

**Solution:**
- Fail fast with clear error messages
- Validate at I/O boundaries
- Add biological constraint validation (BP distance, etc.)
- Log skipped entries with reason

---

### 4. ⚠️ Complex State Management (41 Intron Attributes)

**Problem:**
The Intron class has 41 `__slots__` attributes, many mutable, set at different pipeline stages:

```python
__slots__ = [
    '__dict__', 'bp_raw_score', 'bp_region_seq', 'bp_seq',
    'bp_start', 'bp_stop', 'bp_z_score', 'corrected',
    'dnts', 'downstream_exon', 'downstream_flank', 'duplicate',
    # ... 29 more
]
```

**Why This Is Concerning:**
- Hard to track when/where each attribute is set
- Attributes may be `None` at various pipeline stages
- Easy to access attribute before it's set
- No clear state transitions
- Hard to debug: "why is this None?"

**Example Issue (Line 4758):**
```python
#TODO the problem is that for a duplicate short intron, say,
# we need to keep the info for the duplicate map...
# to make sure the intron isn't omitted, but that's hacky.
```

**Refactoring Priority:** 🟡 **MEDIUM**

**Solution:**
- Group related attributes into dataclasses (IntronScores, IntronSequences, etc.)
- Make clear which attributes are required vs optional at each stage
- Use composition instead of flat 41-attribute class
- Consider state machine pattern for pipeline stages

---

### 5. ⚠️ No Automated Test Suite

**Problem:**
- No unit tests
- No integration tests
- No regression tests
- Only example data (Chr19)

**Why This Is Critical for Refactoring:**
- Can't safely refactor without tests
- No way to verify correctness
- Recent bugs suggest things break
- Manual testing is incomplete

**Evidence of Issues:**
Multiple bug fixes in recent commits suggest testing gaps:
- Off-by-one errors
- Coordinate adjustment bugs
- Version string issues (multiple fixes)

**Refactoring Priority:** 🔴 **CRITICAL**

**Solution:**
- Write regression test FIRST: run original on chr19, save results
- Refactored code must produce identical results
- Add unit tests for coordinate handling
- Add unit tests for PWM scoring
- Test boundary conditions

---

## Moderate Issues (Should Address)

### 6. Missing Field-Proven Biological Constraints

**Issue:**
Based on field analysis (FIELD_APPROACHES.md), intronIC is missing validated constraints:

1. **Branch point distance validation** (8-30 nt from 3'ss)
   - Sheth et al. method uses this
   - Strong biological constraint
   - Not validated in intronIC

2. **Explicit polypyrimidine tract detection**
   - U12 introns lack strong PPT
   - Currently implicit via 3'ss scoring
   - Should be explicit feature

3. **Score differentials**
   - Can't compare U12 vs U2 scores directly
   - Sheth et al. uses thresholds (25, 10)
   - Useful for interpretation

**Priority:** 🟡 **MEDIUM** (enhancement, not bug)

**Solution:**
Add during refactoring as new features

---

### 7. Silent Failures and Poor Error Messages

**Issue:**
Many try/except blocks that silently skip:

```python
try:
    # parse line
except:
    continue  # No logging, no error message
```

**Why Problematic:**
- User doesn't know data was skipped
- Can't debug malformed input
- Silent failures hide real problems

**Priority:** 🟡 **MEDIUM**

**Solution:**
- Log skipped entries with reason
- Provide clear error messages
- Fail fast on critical errors

---

### 8. Complex Strand Handling

**Issue:**
Strand handling scattered throughout code with different conventions:

```python
# Line 2348-2349: u12_correction
if intron.strand == '-':
    shift *= -1  # Reverse adjustment

# Line 2413-2416: _correct_exon
if exon.strand == '-':
    shift = shift * -1
    target = next(e for e in ('start', 'stop') if e != target)
```

**Why Concerning:**
- Easy to forget strand check
- Different approaches in different functions
- Hard to verify correctness

**Priority:** 🟢 **LOW-MEDIUM**

**Solution:**
- Centralize strand handling
- Always work with forward strand internally
- Convert only at I/O boundaries

---

## Minor Issues (Nice to Fix)

### 9. Performance: Non-Vectorized Scoring

**Issue:**
PWM scoring is character-by-character loop (line 2114-2140)

**Impact:**
- Acceptable (human genome in 20-35 min)
- Could be 2-10x faster with vectorization

**Priority:** 🟢 **LOW**

---

### 10. Documentation Gaps

**Issue:**
- Many functions lack docstrings
- Complex algorithms not explained
- No type hints

**Priority:** 🟢 **LOW** (but should fix during refactoring)

---

## Refactoring Implications

### What This Means for Our Approach

**1. Test-First Strategy is ESSENTIAL**
   - Must create regression test before refactoring
   - Run original intronIC on chr19
   - Save all outputs as "gold standard"
   - Refactored version must match exactly

**2. Coordinate Handling Requires Extreme Care**
   - This is where bugs have occurred
   - Need explicit coordinate system tracking
   - Validate at every conversion
   - Add extensive tests for boundary conditions

**3. Can't Blindly Copy Logic**
   - Some logic has known issues (TODOs)
   - u12_correction needs redesign
   - Silent failures need to be fixed

**4. Immutability Where Possible**
   - Reduce state mutation
   - Make data flow explicit
   - Easier to test and debug

---

## Recommendations for Refactoring Order

### Phase 0: Foundation (CRITICAL - DO FIRST)
1. ✅ **Create regression test suite**
   - Run original on chr19
   - Save all outputs
   - This is our safety net

2. ✅ **Create coordinate system abstraction**
   - GenomicCoordinate class
   - Explicit conversion functions
   - Validation at creation

### Phase 1: Core (Build on Safe Foundation)
3. Data models with immutable coordinates
4. I/O with explicit coordinate conversion
5. Test against regression suite

### Phase 2: Processing (Careful with State)
6. Extraction pipeline
7. Sequence handling
8. Filtering (careful with duplicates)

### Phase 3: Scoring (Most Stable Part)
9. PWM scoring (mostly straightforward)
10. Feature extraction
11. Normalization

### Phase 4: Classification (Add Enhancements)
12. SVM (well-tested, low risk)
13. **NEW:** Validators (BP distance, PPT)
14. Ensemble

### Phase 5: Advanced (Careful with Corrections)
15. u12_correction (REDESIGN, don't just copy)
16. Recursive training
17. Plots

---

## Risk Assessment

### Highest Risk Areas (Test Extensively)
1. 🔴 Coordinate conversions (BED ↔ internal)
2. 🔴 u12_correction boundary adjustment
3. 🔴 Strand handling in coordinate calculations
4. 🔴 Duplicate detection logic
5. 🟡 Phase adjustment (complex, subtle)

### Medium Risk Areas (Test Thoroughly)
6. 🟡 Non-canonical intron handling
7. 🟡 Isoform selection (longest)
8. 🟡 Overlap detection

### Low Risk Areas (Well-Understood)
9. 🟢 PWM scoring (well-defined math)
10. 🟢 SVM training (sklearn handles it)
11. 🟢 Z-score normalization (standard)
12. 🟢 File I/O (straightforward)

---

## Bottom Line

**The existing code works and produces good scientific results**, but has several areas where bugs can creep in:

1. **Coordinate handling is the #1 risk** - Recent bugs prove this
2. **State mutation makes code hard to reason about**
3. **No tests makes refactoring risky**
4. **Some design decisions have known issues** (see TODOs)

**Our refactoring must:**
- ✅ Create regression tests FIRST
- ✅ Fix coordinate handling with explicit types
- ✅ Reduce state mutation
- ✅ Add validation and fail-fast
- ✅ Not blindly copy problematic patterns

**We should NOT:**
- ❌ Assume all existing logic is correct (check TODOs)
- ❌ Refactor without regression tests
- ❌ Ignore the coordinate system complexity
- ❌ Skip validation "to save time"

---

**Next Step:** Create regression test suite from chr19 before touching any code.
