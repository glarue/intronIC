# intronIC Refactoring: Execution Plan

**Realistic Timeline:** Work in focused chunks, ~2-3 weeks total
**Approach:** Build incrementally, test as we go, maintain working state

---

## Execution Strategy

### Principles
1. **Work in logical chunks** - Complete one module at a time
2. **Keep it working** - Each chunk should be runnable
3. **Test immediately** - Validate each piece before moving on
4. **Iterate fast** - No bureaucracy, just execute

---

## Chunk Breakdown

### Chunk 1: Foundation & Core Models (Start Here)
**Time:** ~2-3 hours of focused work
**Goal:** Basic structure and data models

**Tasks:**
```
1. Create directory structure
2. Build core/models.py (GenomeFeature, Gene, Transcript, Exon)
3. Build core/intron.py (refactored Intron class)
4. Build utils/sequences.py (rev_comp, basic helpers)
5. Build utils/coordinates.py (coordinate handling)
6. Write basic tests for models
```

**Deliverable:** Can create Intron objects, test passes

---

### Chunk 2: I/O Layer
**Time:** ~3-4 hours
**Goal:** Read files, write outputs

**Tasks:**
```
1. Build io/genome.py (GenomeReader with caching)
2. Build io/parsers.py (AnnotationParser, BEDParser, SequenceParser)
3. Build io/writers.py (output file generation)
4. Test with actual chr19 data
```

**Deliverable:** Can read genome + annotation, extract sequences

---

### Chunk 3: Extraction Pipeline
**Time:** ~4-5 hours
**Goal:** Annotation → Introns

**Tasks:**
```
1. Build extraction/annotator.py (hierarchy parsing)
2. Build extraction/intronator.py (exon pairs → introns)
3. Build extraction/sequences.py (sequence extraction)
4. Build extraction/filters.py (duplicates, longest isoform)
5. Integration test: annotation → filtered introns
```

**Deliverable:** Can extract introns from chr19, matches original count

---

### Chunk 4: Scoring System
**Time:** ~4-5 hours
**Goal:** PWM scoring + new features

**Tasks:**
```
1. Build scoring/pwm.py (matrix loading)
2. Build scoring/scorer.py (seq_score, bp_score)
3. Build scoring/features.py (extract 5'ss, BP, 3'ss)
4. Build scoring/normalization.py (z-scores)
5. Build scoring/ppt.py (NEW: PPT detection)
6. Test scoring matches original
```

**Deliverable:** Can score introns, z-scores match original

---

### Chunk 5: Classification & New Features
**Time:** ~4-5 hours
**Goal:** SVM + field improvements

**Tasks:**
```
1. Build classification/svm.py (training, prediction)
2. Build classification/optimization.py (hyperparameter tuning)
3. Build classification/validators.py (NEW: BP distance, confidence)
4. Build classification/ensemble.py (multi-model)
5. Build classification/recursive.py (recursive training)
6. Integration test: full pipeline
```

**Deliverable:** Can classify introns, new validators working

---

### Chunk 6: CLI & Configuration
**Time:** ~2-3 hours
**Goal:** User interface

**Tasks:**
```
1. Build config/defaults.py (configuration dataclasses)
2. Build cli/parser.py (argparse wrapper)
3. Build cli/main.py (entry point)
4. Test command-line interface
```

**Deliverable:** Can run from command line like original

---

### Chunk 7: Analysis & Visualization
**Time:** ~2-3 hours
**Goal:** Plots and statistics

**Tasks:**
```
1. Build analysis/plotting.py (scatter, PR-AUC)
2. Build analysis/statistics.py (summary stats)
3. Test plot generation
```

**Deliverable:** Generates all output files and plots

---

### Chunk 8: Testing & Validation
**Time:** ~3-4 hours
**Goal:** Comprehensive test suite

**Tasks:**
```
1. Write unit tests for each module
2. Write integration test (full pipeline)
3. Write regression test (compare to original)
4. Verify >85% coverage
```

**Deliverable:** Test suite passes, validates against original

---

### Chunk 9: Documentation & Polish
**Time:** ~2-3 hours
**Goal:** Finalize for use

**Tasks:**
```
1. Add docstrings where missing
2. Type check with mypy
3. Format with black
4. Write API documentation
5. Create migration guide
```

**Deliverable:** Production-ready v2.0

---

## Actual Timeline

**If we work through chunks sequentially:**
- 1 chunk per session
- ~3-5 hours per chunk
- 9 chunks total
- **= ~30-40 hours of actual work**
- **= Can complete in 1-2 weeks of focused effort**

**If we parallelize where possible:**
- Some chunks can be done simultaneously
- **Could complete in as little as 1 week**

---

## Let's Start: Chunk 1 Proposal

I can begin **right now** with Chunk 1 if you'd like. Here's what I'll create:

```
intronIC_refactored/
├── __init__.py
├── core/
│   ├── __init__.py
│   ├── models.py          # GenomeFeature, Gene, Transcript, Exon
│   └── intron.py          # Intron with clean structure
└── utils/
    ├── __init__.py
    ├── sequences.py       # rev_comp, sequence helpers
    └── coordinates.py     # Coordinate handling
```

**This takes ~30 minutes to write, then you can review.**

Shall I start with Chunk 1?

---

## Working Agreement

**My role:**
- Write complete, working modules
- Add type hints and docstrings
- Create tests
- Explain design decisions

**Your role:**
- Review chunks as completed
- Test with real data if desired
- Provide feedback/direction
- Approve before moving to next chunk

**Communication:**
- I complete a chunk
- Share code for review
- You validate/test
- Provide feedback
- I iterate if needed
- Move to next chunk

This way we maintain quality while moving quickly.

Ready to start?
