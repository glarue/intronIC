# intronIC Development TODO

This document tracks planned features, enhancements, and known issues for future development.

---

## Priority Levels

- 🔴 **HIGH** - Critical for next release or blocks other work
- 🟡 **MEDIUM** - Important but not blocking
- 🟢 **LOW** - Nice to have, implement when convenient
- 💭 **RESEARCH** - Requires investigation or proof of concept

---

## Feature Requests

### 🟢 Flexible Feature Selection for Training
**Status:** Deferred (post-v2.0)
**Complexity:** Moderate-High (~19 hours)
**Issue:** Currently requires all three motifs (5'SS, BP, 3'SS). Cannot train with subset.

**Use Cases:**
- Train using only 5'SS scores (1D model)
- Train using 5'SS + BP (2D model, skip 3'SS)
- Train using 5'SS + 3'SS (2D model, skip BP)
- Feature ablation studies to understand motif importance

**Implementation Plan:** See `FEATURE_SELECTION_IMPLEMENTATION_PLAN.md`

**Effort Estimate:**
- Configuration changes: 2 hours
- Scoring/normalization: 4 hours
- Transformer refactor: 3 hours
- Training/prediction updates: 2 hours
- Testing: 6 hours
- **Total:** ~19 hours (2-3 days)

**Files Affected:** 16 files, ~600-800 LOC changes

**Blocking Issues:** None

**Decision Criteria to Proceed:**
- User explicitly requests single-feature training
- Research question requires feature ablation
- Performance issues with BP scoring in certain species
- After v2.0 release and architecture stabilization

**Related Files:**
- Implementation plan: `FEATURE_SELECTION_IMPLEMENTATION_PLAN.md`
- Config: `config/config.yaml`, `src/intronIC/cli/config.py`
- Core logic: `src/intronIC/classification/transformers.py`

---

## Known Issues

### 🔴 [Add any critical bugs here]

---

## Enhancements

### 🟡 [Add medium priority enhancements here]

---

## Documentation

### 🟢 [Add documentation TODOs here]

---

## Performance Optimizations

### 🟢 [Add optimization ideas here]

---

## Testing

### 🟡 [Add testing improvements here]

---

## Technical Debt

### 🟢 [Add refactoring needs here]

---

## Research Ideas

### 💭 [Add research questions here]

---

## Completed Items

Items are moved here after implementation.

---

## Notes

- Keep this file updated as new issues/features are identified
- Link to detailed planning documents for complex features
- Update priority levels as project needs change
- Mark items as completed with date and commit hash
