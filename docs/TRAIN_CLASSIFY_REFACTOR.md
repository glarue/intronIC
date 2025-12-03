# Train/Classify Architecture Refactor

**Date:** December 3, 2025
**Status:** Design Phase

## Problem Statement

The current intronIC architecture has two confusing and overlapping training paradigms:

1. **`intronIC train`** - Pure training on reference sequences (no genome/annotation)
2. **`intronIC classify --train`** - Extract from genome/annotation + train + classify in one operation

### Issues with Current Design

1. **Conceptual confusion**: Training on a genome/annotation doesn't make semantic sense. Training should use curated reference sequences, not arbitrary genomic introns.

2. **Lack of separation**: `intronIC classify --train` conflates three distinct operations:
   - Intron extraction from genome/annotation
   - Model training on extracted introns
   - Classification of extracted introns

3. **No intermediate checkpoints**: Users can't inspect extracted introns before training, or train a model and then apply it to different data.

4. **Inconsistent workflow**: The `train` subcommand does one thing (pure training), while `--train` flag does something completely different (extract + train + classify).

## Current Architecture

```
┌─────────────────────────────────────────────────────────────┐
│ intronIC train                                              │
│ ────────────────────────────────────────────────────────── │
│ Input: Reference U12/U2 sequences (.iic files)             │
│ Output: Trained model (.model.pkl)                         │
│                                                             │
│ Flow: Load references → Score → Normalize → Train → Save   │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│ intronIC classify --train                                   │
│ ────────────────────────────────────────────────────────── │
│ Input: Genome + Annotation                                  │
│ Output: Classified introns + Trained model                  │
│                                                             │
│ Flow:                                                       │
│   1. Extract introns from genome/annotation                 │
│   2. Train model on extracted introns (uses as "reference") │
│   3. Classify same extracted introns                        │
│   4. Write outputs                                          │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│ intronIC classify --model MODEL.pkl                         │
│ ────────────────────────────────────────────────────────── │
│ Input: Genome + Annotation + Pretrained Model               │
│ Output: Classified introns                                  │
│                                                             │
│ Flow: Extract → Score → Normalize → Classify → Write       │
└─────────────────────────────────────────────────────────────┘
```

## Proposed Architecture

### Core Principle
**Training and classification should be completely separate operations.** Training happens on curated reference sequences. Classification happens on genomic data using a trained model.

### Proposed Commands

```
┌─────────────────────────────────────────────────────────────┐
│ intronIC train                                              │
│ ────────────────────────────────────────────────────────── │
│ Input: Reference U12/U2 sequences (.iic files)             │
│ Output: Trained model (.model.pkl)                         │
│                                                             │
│ Purpose: Train a model on curated reference sequences       │
│ Usage: intronIC train -n species_name                       │
│        [--reference-u12s U12.iic] [--reference-u2s U2.iic]  │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│ intronIC classify                                           │
│ ────────────────────────────────────────────────────────── │
│ Input: Genome + Annotation + Pretrained Model               │
│ Output: Classified introns                                  │
│                                                             │
│ Purpose: Classify introns using a pretrained model          │
│ Usage: intronIC classify -g genome.fa -a annotation.gff     │
│        -n species --model species.model.pkl                 │
│                                                             │
│ Note: --model is now REQUIRED (no --train flag)            │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│ intronIC extract (NEW)                                      │
│ ────────────────────────────────────────────────────────── │
│ Input: Genome + Annotation                                  │
│ Output: Extracted intron sequences (.introns.iic)           │
│                                                             │
│ Purpose: Extract introns without classification             │
│ Usage: intronIC extract -g genome.fa -a annotation.gff      │
│        -n species                                           │
│                                                             │
│ Equivalent to current: intronIC classify -s (sequences-only)│
└─────────────────────────────────────────────────────────────┘
```

### Migration Path for `--train` Users

For users who currently use `intronIC classify --train`, we provide a clear two-step workflow:

**Old workflow (deprecated):**
```bash
# Extract + train + classify in one step (confusing!)
intronIC classify -g genome.fa -a annotation.gff -n species --train
```

**New workflow (clear separation):**
```bash
# Step 1: Train a model on built-in references
intronIC train -n species

# Step 2: Classify genomic introns with trained model
intronIC classify -g genome.fa -a annotation.gff -n species --model species.model.pkl
```

**Alternative workflow (custom references):**
```bash
# Step 1: Extract introns from genome
intronIC extract -g genome.fa -a annotation.gff -n species

# Step 2: (Manual) Curate introns into U12/U2 reference sets
# User manually creates curated_u12.iic and curated_u2.iic

# Step 3: Train model on curated references
intronIC train -n species --reference-u12s curated_u12.iic --reference-u2s curated_u2.iic

# Step 4: Classify introns with new model
intronIC classify -g genome.fa -a annotation.gff -n species --model species.model.pkl
```

## Implementation Plan

### Phase 1: Add `extract` Subcommand (Low Risk)

This is essentially renaming/exposing the existing `-s/--sequences-only` functionality:

1. Add `extract` subcommand to argument parser
2. Route to extraction-only pipeline (already exists)
3. Update documentation

**Backward compatibility:** Keep `-s` flag working in `classify` mode

### Phase 2: Deprecate `--train` Flag (Breaking Change)

1. Add deprecation warning when `--train` is used
2. Update error message to suggest new workflow:
   ```
   ⚠️  The --train flag is deprecated and will be removed in v3.0.0

   Please use the new two-step workflow:
     1. intronIC train -n species
     2. intronIC classify -g genome.fa -a annotation.gff --model species.model.pkl
   ```
3. Update all documentation and examples

**Backward compatibility:** Keep `--train` working with deprecation warning in v2.x

### Phase 3: Make `--model` Required in `classify` (Breaking Change)

After deprecation period (v3.0.0):

1. Remove `--train` flag entirely
2. Make `--model` required for `classify` command
3. Clean up training logic from `main_classify()`

## Benefits of Refactor

1. **Clear semantics**: Training = reference sequences, Classification = genomic data
2. **Modularity**: Each command does one thing well
3. **Inspectable**: Users can examine extracted introns before training
4. **Flexible workflows**: Train once, classify many times on different data
5. **Better testing**: Each operation can be tested independently
6. **Reduced complexity**: `main_classify()` only handles classification

## Edge Cases and Considerations

### 1. Default Model
**Question:** Should we ship a default model so users don't need to train?

**Proposal:**
- Ship with a default "universal" model trained on diverse species
- Users can optionally train species-specific models for better accuracy
- `--model` remains required, but we provide a built-in model path

### 2. Quick Prototyping
**Question:** What if users want to quickly test on new data without two steps?

**Proposal:** Consider a convenience command:
```bash
intronIC quick -g genome.fa -a annotation.gff -n species
# Internally runs: train (on defaults) → classify
```

This makes the two-step process explicit even in convenience mode.

### 3. Custom Training from Genome
**Question:** What if users truly want to train on genomic introns as "references"?

**Answer:** They can, but it's explicit:
```bash
# Extract
intronIC extract -g genome.fa -a annotation.gff -n species

# Manually label introns as U12/U2 (external tools/manual curation)
# Creates: u12_labeled.iic, u2_labeled.iic

# Train on labeled set
intronIC train -n species --reference-u12s u12_labeled.iic --reference-u2s u2_labeled.iic
```

This makes it clear that training requires curated labels, not raw genomic extractions.

## Timeline

- **v2.1.0 (current)**: Add `extract` subcommand, add deprecation warning for `--train`
- **v2.2.0**: Update all documentation to show new workflows
- **v3.0.0**: Remove `--train` flag, make `--model` required

## Related Files

- [src/intronIC/cli/args.py](../src/intronIC/cli/args.py) - Argument parsing
- [src/intronIC/cli/main.py](../src/intronIC/cli/main.py) - Main entry points
- [src/intronIC/cli/config.py](../src/intronIC/cli/config.py) - Configuration dataclasses

## Questions for Discussion

1. Should we provide a default/universal model, or always require users to train?
2. Should `extract` be a subcommand or keep as `-s` flag in `classify`?
3. What's an appropriate deprecation timeline for `--train`?
4. Should we add a `quick` convenience command for one-step train+classify?
