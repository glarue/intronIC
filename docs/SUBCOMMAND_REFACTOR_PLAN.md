# Subcommand Refactor Plan: Separate Train and Classify

**Date:** 2025-11-15
**Goal:** Split `--train` flag into separate `train` and `classify` subcommands
**Rationale:** Training is slow and should be decoupled from extraction/classification

---

## Current Problem

Currently, to train a model you must provide genome/annotation files even though training only uses reference data:

```bash
# Current (wasteful - need dummy genome/annotation just to train)
intronIC -g dummy.fa -a dummy.gff -n species --train

# After training, classify with pretrained model
intronIC -g real.fa -a real.gff -n species --pretrained_model species.model.pkl
```

**Issues:**
1. Must provide genome/annotation even for pure training
2. `--train` flag is ambiguous (train-then-classify vs. train-only)
3. Can't easily train models without dummy input files
4. Mixes two distinct workflows (training vs. inference)

---

## Proposed Solution

Create two subcommands:

### 1. `intronIC train` - Train model only

```bash
intronIC train -n species_name [training_options]
```

**Input:**
- Reference data (U12/U2 sequences) - built-in or custom
- Training parameters (C, n_models, CV settings, etc.)
- Species name for output

**Output:**
- `{species_name}.model.pkl` - Trained model
- `{species_name}.model.metadata.json` - Training metadata
- `{species_name}.log` - Training log

**No genome/annotation needed!**

### 2. `intronIC classify` (or default) - Extract and classify

```bash
# With pretrained model
intronIC classify -g genome.fa -a annotation.gff -n species --model species.model.pkl

# Or train on-the-fly (current --train behavior)
intronIC classify -g genome.fa -a annotation.gff -n species --train

# Default (no subcommand) = classify
intronIC -g genome.fa -a annotation.gff -n species --model species.model.pkl
```

**Input:**
- Genome/annotation/BED/sequences (at least one)
- Species name
- Model source: `--model path.pkl` OR `--train` (train on-the-fly)

**Output:**
- All current outputs (.bed, .meta, .introns, .score_info, etc.)
- Model file (if `--train` used)

---

## Implementation Plan

### Phase 1: Analysis & Design

**1.1 Analyze current argument structure**
- Review `cli/args.py` (or wherever args are defined)
- Identify training-specific args vs. extraction/classification args
- Document dependencies between args

**1.2 Design subcommand structure**
```
intronIC
├── train
│   ├── Required: -n/--species_name
│   ├── Optional: --reference_u12s, --reference_u2s
│   ├── Training: -C, --n_models, --eval_mode, --n_cv_folds, etc.
│   └── Output: -o/--output_dir
├── classify (or default)
│   ├── Required: -n/--species_name AND (genome+annotation OR bed OR sequences)
│   ├── Input: -g, -a, -b, -q
│   ├── Model: --model path.pkl OR --train
│   ├── Extraction: -f, --min_intron_len, -i, -v, -d, etc.
│   ├── Scoring: -t, --no_nc, --pseudocount, etc.
│   └── Output: -o, --clean_names, etc.
└── Global args: --quiet, --debug, --config, --version, --help
```

**1.3 Plan backward compatibility**
- Support old CLI without subcommands (default to `classify`)
- Support `--train` flag (maps to `--train` in `classify` mode)
- Deprecation warnings for old usage patterns

---

### Phase 2: Argument Parser Refactoring

**2.1 Create subparser structure**

**File:** `cli/args.py` (or new file `cli/subcommands.py`)

```python
def create_argument_parser() -> argparse.ArgumentParser:
    """Create main parser with subcommands."""

    parser = argparse.ArgumentParser(
        prog='intronIC',
        description='Intron classification pipeline'
    )

    # Global arguments (apply to all subcommands)
    parser.add_argument('--version', action='version', version='2.0.0')
    parser.add_argument('--quiet', action='store_true')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--config', type=str, help='Config file')

    # Create subparsers
    subparsers = parser.add_subparsers(
        dest='command',
        help='Command to run'
    )

    # === TRAIN SUBCOMMAND ===
    train_parser = subparsers.add_parser(
        'train',
        help='Train a classifier on reference data'
    )
    add_train_arguments(train_parser)

    # === CLASSIFY SUBCOMMAND ===
    classify_parser = subparsers.add_parser(
        'classify',
        help='Extract and classify introns'
    )
    add_classify_arguments(classify_parser)

    # === BACKWARD COMPATIBILITY ===
    # If no subcommand given, add all classify args to main parser
    # This allows: intronIC -g genome.fa -a annotation.gff -n species
    add_classify_arguments(parser)

    return parser
```

**2.2 Define train-specific arguments**

```python
def add_train_arguments(parser: argparse.ArgumentParser):
    """Add arguments for train subcommand."""

    # Required
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-n', '--species_name',
        required=True,
        help='Species name for output files'
    )

    # Reference data
    reference = parser.add_argument_group('reference data')
    reference.add_argument(
        '--reference_u12s',
        type=str,
        help='Custom U12 reference sequences (default: built-in)'
    )
    reference.add_argument(
        '--reference_u2s',
        type=str,
        help='Custom U2 reference sequences (default: built-in)'
    )

    # Training parameters
    training = parser.add_argument_group('training parameters')
    training.add_argument('-C', type=float, help='Fixed C value (skip optimization)')
    training.add_argument('--n_models', type=int, default=4, help='Number of ensemble models')
    training.add_argument('--eval_mode', choices=['nested_cv', 'split', 'none'], default='nested_cv')
    training.add_argument('--n_cv_folds', type=int, default=5)
    training.add_argument('--n_optimization_rounds', type=int, default=5)
    training.add_argument('--max_iter', type=int, default=10000)
    training.add_argument('--cv_processes', type=int, default=1)
    training.add_argument('--seed', type=int, help='Random seed')

    # Scoring regions (for scoring reference data during training)
    scoring = parser.add_argument_group('scoring parameters')
    scoring.add_argument('--five_score_coords', nargs=2, type=int, metavar=('START', 'END'))
    scoring.add_argument('--bp_region_coords', nargs=2, type=int, metavar=('START', 'END'))
    scoring.add_argument('--three_score_coords', nargs=2, type=int, metavar=('START', 'END'))
    scoring.add_argument('--pseudocount', type=float, default=0.0)

    # Output
    output = parser.add_argument_group('output options')
    output.add_argument('-o', '--output_dir', type=str, default='.')
    output.add_argument('--optimizer-config', type=str, help='Custom optimizer config')
```

**2.3 Define classify-specific arguments**

```python
def add_classify_arguments(parser: argparse.ArgumentParser):
    """Add arguments for classify subcommand (or main parser for backward compat)."""

    # Required (at least one input source)
    required = parser.add_argument_group('required arguments')
    required.add_argument('-n', '--species_name', required=True)

    # Input sources (at least one required)
    input_group = parser.add_argument_group('input sources (at least one required)')
    input_group.add_argument('-g', '--genome', type=str)
    input_group.add_argument('-a', '--annotation', type=str)
    input_group.add_argument('-b', '--bed', type=str)
    input_group.add_argument('-q', '--sequence_file', type=str)

    # Model source (choose one)
    model_group = parser.add_argument_group('model source (choose one)')
    model_mutex = model_group.add_mutually_exclusive_group()
    model_mutex.add_argument('--model', type=str, help='Path to pretrained model')
    model_mutex.add_argument('--train', action='store_true', help='Train model on-the-fly')

    # Extraction parameters
    extraction = parser.add_argument_group('extraction parameters')
    extraction.add_argument('-f', '--feature', choices=['cds', 'exon', 'both'], default='both')
    extraction.add_argument('--min_intron_len', type=int, default=30)
    extraction.add_argument('-i', '--allow_multiple_isoforms', action='store_true')
    extraction.add_argument('-v', '--no_intron_overlap', action='store_true')
    extraction.add_argument('-d', '--include_duplicates', action='store_true')
    extraction.add_argument('--flank_len', type=int, default=400)
    extraction.add_argument('--no_nc_ss_adjustment', action='store_false', dest='u12_boundary_correction')

    # Scoring parameters
    scoring = parser.add_argument_group('scoring parameters')
    scoring.add_argument('-s', '--sequences_only', action='store_true')
    scoring.add_argument('-t', '--threshold', type=int, default=90)
    scoring.add_argument('--no_nc', action='store_true')
    scoring.add_argument('--pseudocount', type=float, default=0.0)
    scoring.add_argument('--no_ignore_nc_dnts', action='store_false', dest='ignore_nc_dnts')
    scoring.add_argument('--five_score_coords', nargs=2, type=int, metavar=('START', 'END'))
    scoring.add_argument('--bp_region_coords', nargs=2, type=int, metavar=('START', 'END'))
    scoring.add_argument('--three_score_coords', nargs=2, type=int, metavar=('START', 'END'))

    # Training (only if --train flag used)
    training = parser.add_argument_group('training parameters (only with --train)')
    training.add_argument('--reference_u12s', type=str)
    training.add_argument('--reference_u2s', type=str)
    training.add_argument('-C', type=float)
    training.add_argument('--n_models', type=int, default=4)
    training.add_argument('--eval_mode', choices=['nested_cv', 'split', 'none'], default='nested_cv')
    training.add_argument('--n_cv_folds', type=int, default=5)
    # ... etc.

    # Performance
    performance = parser.add_argument_group('performance')
    performance.add_argument('-p', '--processes', type=int, default=1)
    performance.add_argument('--cv_processes', type=int, default=1)

    # Output
    output = parser.add_argument_group('output options')
    output.add_argument('-o', '--output_dir', type=str, default='.')
    output.add_argument('--clean_names', action='store_true')
    output.add_argument('--no_abbreviate', action='store_true')
    # ... etc.
```

**2.4 Backward compatibility helper**

```python
def parse_args_with_backward_compat(args=None):
    """Parse args with backward compatibility for old CLI."""

    parser = create_argument_parser()
    parsed = parser.parse_args(args)

    # If no subcommand given, default to 'classify'
    if not hasattr(parsed, 'command') or parsed.command is None:
        parsed.command = 'classify'

        # Check for old --pretrained_model flag (rename to --model)
        if hasattr(parsed, 'pretrained_model') and parsed.pretrained_model:
            parsed.model = parsed.pretrained_model

    # Validate inputs based on command
    if parsed.command == 'classify':
        # Must have at least one input source
        if not any([parsed.genome, parsed.annotation, parsed.bed, parsed.sequence_file]):
            parser.error("classify: at least one input source required (-g, -a, -b, or -q)")

        # Must have model source
        if not parsed.model and not parsed.train:
            parser.error("classify: must specify --model or --train")

    return parsed
```

---

### Phase 3: Main Entry Point Refactoring

**3.1 Create separate main functions**

**File:** `cli/main.py`

```python
def main():
    """Main entry point - route to subcommand."""

    args = parse_args_with_backward_compat()

    if args.command == 'train':
        return main_train(args)
    elif args.command == 'classify':
        return main_classify(args)
    else:
        raise ValueError(f"Unknown command: {args.command}")


def main_train(args):
    """Train a model on reference data only.

    This is a pure training workflow:
    1. Load reference U12/U2 sequences
    2. Extract scoring regions from sequences
    3. Calculate PWM scores
    4. Normalize scores
    5. Optimize hyperparameters
    6. Train ensemble of models
    7. Save model to disk

    No genome/annotation needed!
    """

    # Setup
    config = create_config_from_args(args, mode='train')
    messenger, reporter, logger = setup_logging_and_reporting(config)

    messenger.header("intronIC TRAINING MODE")
    messenger.info(f"Training model: {config.output.species_name}")

    # Step 1: Load reference sequences
    messenger.step(1, "Load Reference Data")
    u12_reference, u2_reference = load_reference_sequences(config, messenger)
    messenger.success(f"Loaded {len(u12_reference)} U12 and {len(u2_reference)} U2 reference introns")

    # Step 2: Score reference sequences
    messenger.step(2, "Score Reference Sequences")
    scored_u12 = score_introns(u12_reference, config, messenger, reporter)
    scored_u2 = score_introns(u2_reference, config, messenger, reporter)
    messenger.success("Reference sequences scored")

    # Step 3: Normalize scores
    messenger.step(3, "Normalize Scores")
    combined = scored_u12 + scored_u2
    normalized, _, _, normalizer = normalize_scores(combined, config, messenger, reporter)
    messenger.success("Scores normalized")

    # Step 4: Train model
    messenger.step(4, "Train Classifier")
    model, metrics = train_classifier(
        normalized, u12_reference, u2_reference, normalizer, config, messenger, reporter
    )
    messenger.success("Model trained")

    # Step 5: Save model
    messenger.step(5, "Save Model")
    model_path = save_model(model, config, metrics, messenger)
    messenger.success(f"Model saved: {model_path}")

    # Done
    messenger.success("Training complete!")
    return 0


def main_classify(args):
    """Extract and classify introns.

    This is the standard workflow:
    1. Load input data (genome/annotation/bed/sequences)
    2. Extract introns
    3. Filter introns
    4. Score introns
    5. Load or train model
    6. Classify introns
    7. Write outputs

    Requires genome/annotation OR bed OR sequences.
    Requires model (--model path.pkl OR --train).
    """

    # Setup
    config = create_config_from_args(args, mode='classify')
    messenger, reporter, logger = setup_logging_and_reporting(config)

    messenger.header("intronIC CLASSIFICATION MODE")

    # ... existing main() logic ...
    # (current pipeline code goes here)

    return 0
```

**3.2 Refactor config creation**

```python
def create_config_from_args(args, mode='classify'):
    """Create config from parsed args.

    Args:
        args: Parsed arguments
        mode: 'train' or 'classify'

    Returns:
        IntronICConfig object
    """

    if mode == 'train':
        # Training-specific config
        return IntronICConfig(
            input=InputConfig(
                mode='train',
                # No genome/annotation needed
            ),
            training=TrainingConfig(
                reference_u12s=args.reference_u12s,
                reference_u2s=args.reference_u2s,
                n_models=args.n_models,
                eval_mode=args.eval_mode,
                # ... etc.
            ),
            output=OutputConfig(
                species_name=args.species_name,
                output_dir=args.output_dir,
                # ...
            ),
            # ... other configs
        )

    elif mode == 'classify':
        # Classification-specific config
        return IntronICConfig(
            input=InputConfig(
                mode=determine_input_mode(args),
                genome=args.genome,
                annotation=args.annotation,
                bed=args.bed,
                sequence_file=args.sequence_file,
            ),
            training=TrainingConfig(
                pretrained_model_path=args.model if hasattr(args, 'model') else None,
                train_on_fly=args.train if hasattr(args, 'train') else False,
                # ... etc.
            ),
            # ... other configs
        )
```

---

### Phase 4: Model Save/Load Functions

**4.1 Standardize model file format**

**File:** `utils/model_io.py` (new file)

```python
from pathlib import Path
import joblib
from typing import Dict, Any
from datetime import datetime

def save_model(
    model,
    config: IntronICConfig,
    metrics: Dict[str, Any],
    messenger
) -> Path:
    """Save trained model to disk.

    Args:
        model: Trained classifier (or ensemble)
        config: Pipeline configuration
        metrics: Training metrics (F1, precision, recall, etc.)
        messenger: Messenger for logging

    Returns:
        Path to saved model file

    Saves:
        - {species_name}.model.pkl - Pickled model
        - {species_name}.model.metadata.json - Training metadata
    """

    output_dir = config.output.output_dir
    base_name = config.output.base_filename

    # Save model
    model_path = output_dir / f"{base_name}.model.pkl"
    joblib.dump(model, model_path)
    messenger.log_only(f"Saved model: {model_path}")

    # Save metadata
    metadata = {
        'species_name': config.output.species_name,
        'trained_date': datetime.now().isoformat(),
        'intronIC_version': '2.0.0',
        'metrics': metrics,
        'training_config': {
            'n_models': config.training.n_models,
            'eval_mode': config.training.eval_mode,
            'n_cv_folds': config.training.n_cv_folds,
            'reference_u12_count': metrics.get('n_u12_reference', 0),
            'reference_u2_count': metrics.get('n_u2_reference', 0),
        },
        'scoring_config': {
            'five_coords': (config.scoring.scoring_regions.five_start,
                           config.scoring.scoring_regions.five_end),
            'bp_coords': (config.scoring.scoring_regions.bp_start,
                         config.scoring.scoring_regions.bp_end),
            'three_coords': (config.scoring.scoring_regions.three_start,
                            config.scoring.scoring_regions.three_end),
            'pseudocount': config.scoring.pseudocount,
        }
    }

    metadata_path = output_dir / f"{base_name}.model.metadata.json"
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    messenger.log_only(f"Saved metadata: {metadata_path}")

    return model_path


def load_model(model_path: Path, messenger) -> Any:
    """Load trained model from disk.

    Args:
        model_path: Path to .model.pkl file
        messenger: Messenger for logging

    Returns:
        Loaded model

    Raises:
        FileNotFoundError: If model file not found
        Exception: If model loading fails
    """

    if not model_path.exists():
        raise FileNotFoundError(f"Model file not found: {model_path}")

    try:
        model = joblib.load(model_path)
        messenger.log_only(f"Loaded model: {model_path}")

        # Load metadata if available
        metadata_path = model_path.with_suffix('.metadata.json')
        if metadata_path.exists():
            with open(metadata_path) as f:
                metadata = json.load(f)
            messenger.log_only(f"Model trained on: {metadata.get('trained_date', 'unknown')}")
            messenger.log_only(f"Training metrics: F1={metadata.get('metrics', {}).get('f1', 'N/A')}")

        return model

    except Exception as e:
        raise Exception(f"Failed to load model from {model_path}: {str(e)}")
```

---

### Phase 5: Testing Strategy

**5.1 Test train subcommand**

```bash
# Basic training
intronIC train -n homo_sapiens

# Custom reference data
intronIC train -n species \
  --reference_u12s custom_u12.iic \
  --reference_u2s custom_u2.iic

# Custom training parameters
intronIC train -n species \
  --n_models 8 \
  --eval_mode split \
  --n_cv_folds 10

# Fixed C value (skip optimization)
intronIC train -n species -C 0.1
```

**Expected outputs:**
- `homo_sapiens.model.pkl`
- `homo_sapiens.model.metadata.json`
- `homo_sapiens.log`

**5.2 Test classify subcommand**

```bash
# With pretrained model
intronIC classify \
  -g genome.fa -a annotation.gff \
  -n species \
  --model species.model.pkl

# Train on-the-fly (old behavior)
intronIC classify \
  -g genome.fa -a annotation.gff \
  -n species \
  --train

# Default (no subcommand) - backward compatible
intronIC \
  -g genome.fa -a annotation.gff \
  -n species \
  --model species.model.pkl
```

**5.3 Test backward compatibility**

```bash
# Old CLI (should still work)
intronIC -g genome.fa -a annotation.gff -n species --train

# Old pretrained model flag (should map to --model)
intronIC -g genome.fa -a annotation.gff -n species --pretrained_model model.pkl
```

---

### Phase 6: Documentation Updates

**6.1 Update README.md**

Add subcommand examples:

```markdown
## Quick Start

### Train a model
```bash
intronIC train -n homo_sapiens
```

### Classify introns with pretrained model
```bash
intronIC classify -g genome.fa.gz -a annotation.gff3.gz -n species --model homo_sapiens.model.pkl
```

### Train and classify in one step
```bash
intronIC classify -g genome.fa.gz -a annotation.gff3.gz -n species --train
```
```

**6.2 Update CLAUDE.md**

Document new CLI structure and subcommands.

**6.3 Create MIGRATION_GUIDE.md**

Help users migrate from old CLI to new subcommand structure.

---

## Implementation Timeline

### Sprint 1: Planning & Design (30-60 min)
- [x] Review current argument structure
- [x] Design subcommand hierarchy
- [x] Plan backward compatibility strategy
- [x] Document implementation plan

### Sprint 2: Argument Parser (1-2 hours)
- [ ] Create subparser structure
- [ ] Define train arguments
- [ ] Define classify arguments
- [ ] Implement backward compatibility
- [ ] Add validation logic
- [ ] Test argument parsing

### Sprint 3: Main Entry Point (1-2 hours)
- [ ] Create `main_train()` function
- [ ] Create `main_classify()` function
- [ ] Refactor current `main()` logic into `main_classify()`
- [ ] Implement subcommand routing
- [ ] Test both workflows

### Sprint 4: Model I/O (30-60 min)
- [ ] Create `save_model()` function
- [ ] Create `load_model()` function
- [ ] Standardize model metadata format
- [ ] Test save/load round-trip

### Sprint 5: Testing (1-2 hours)
- [ ] Test `train` subcommand
- [ ] Test `classify` subcommand with `--model`
- [ ] Test `classify` subcommand with `--train`
- [ ] Test backward compatibility (no subcommand)
- [ ] Test error handling

### Sprint 6: Documentation (30-60 min)
- [ ] Update README.md
- [ ] Update CLAUDE.md
- [ ] Create migration guide
- [ ] Update help text

**Total Estimated Time:** 5-8 hours

---

## Benefits

### For Users

1. **Clearer workflow:** Separate train vs. classify
2. **Faster iteration:** Train once, classify many datasets
3. **No dummy files:** Don't need genome/annotation just to train
4. **Better help:** Subcommand-specific help text
5. **Backward compatible:** Old CLI still works

### For Developers

1. **Separation of concerns:** Training logic separate from extraction
2. **Easier testing:** Can test training without full pipeline
3. **Better organization:** Clear entry points for each workflow
4. **Future extensibility:** Easy to add new subcommands (e.g., `validate`, `compare`)

---

## Future Enhancements

### Additional Subcommands

```bash
# Validate a trained model
intronIC validate --model species.model.pkl --test_data test.iic

# Compare two models
intronIC compare --model1 v1.pkl --model2 v2.pkl --test_data test.iic

# Generate PWM matrices from sequences
intronIC generate-pwms --u12s u12.iic --u2s u2.iic -o matrices.fasta.iic

# Convert between file formats
intronIC convert --input data.iic --output data.bed --format bed
```

---

## Risks & Mitigation

### Risk 1: Breaking backward compatibility
**Mitigation:**
- Support old CLI without subcommands (default to `classify`)
- Add deprecation warnings, don't remove old behavior
- Thorough testing of old CLI patterns

### Risk 2: Increased complexity
**Mitigation:**
- Clear documentation and examples
- Helpful error messages
- Migration guide for users

### Risk 3: Incomplete separation of concerns
**Mitigation:**
- Careful code review
- Ensure `main_train()` doesn't touch extraction code
- Ensure `main_classify()` can work with pretrained models

---

## Success Criteria

- [ ] Can train model without genome/annotation
- [ ] Can classify with pretrained model
- [ ] Can classify with on-the-fly training
- [ ] Old CLI still works (backward compatible)
- [ ] All tests pass
- [ ] Documentation updated
- [ ] Help text clear and accurate

---

## Next Steps

Ready to implement? Let's start with:

1. **Phase 2:** Argument parser refactoring
   - Create subparser structure
   - Define arguments for each subcommand
   - Implement backward compatibility

2. **Phase 3:** Main entry point
   - Create `main_train()` skeleton
   - Create `main_classify()` skeleton
   - Route based on subcommand

3. **Phase 4:** Extract training logic
   - Move training code into `main_train()`
   - Test training without genome/annotation

Let me know when you're ready to proceed!
