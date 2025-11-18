# intronIC Refactoring & Improvement Plan

**Version:** 1.0
**Date:** 2025-11-01
**Target:** Improved clarity, maintainability, and technical approach
**Strategy:** Create refactored version in new subdirectory while preserving original

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [Current State Analysis](#current-state-analysis)
3. [Proposed Architecture](#proposed-architecture)
4. [Module Breakdown](#module-breakdown)
5. [Technical Improvements](#technical-improvements)
6. [Testing Strategy](#testing-strategy)
7. [Implementation Roadmap](#implementation-roadmap)
8. [Migration Strategy](#migration-strategy)
9. [Risk Mitigation](#risk-mitigation)

---

## Executive Summary

### Goals

**Primary Objectives:**
1. **Modularization**: Split 6,093-line monolith into logical, testable modules
2. **Maintainability**: Add type hints, docstrings, and clear interfaces
3. **Technical Enhancement**: Implement field-proven improvements (BP distance, PPT detection)
4. **Testing**: Create comprehensive test suite
5. **Performance**: Optimize bottlenecks through vectorization and caching

**Non-Goals:**
- Complete algorithmic rewrite (preserve proven classification approach)
- Breaking backward compatibility (maintain CLI interface)
- Removing features (additive improvements only)

### Success Metrics

- **Code Quality**: <500 lines per module, 100% type coverage, >80% docstring coverage
- **Performance**: No regression, ideally 10-20% speedup
- **Accuracy**: Match or exceed current classification performance
- **Testing**: >85% code coverage, all critical paths tested
- **Usability**: Identical CLI experience, improved error messages

---

## Current State Analysis

### Pain Points

#### 1. Code Organization
- **Single file**: 6,093 lines in `intronIC.py`
- **~107 functions** in global namespace
- **6 classes** mixed with utility functions
- **Deep nesting**: Some functions 4-5 levels deep
- **No clear separation** of concerns

#### 2. Documentation
- **Inconsistent docstrings**: Many functions lack documentation
- **No type hints**: Makes IDE support and maintenance harder
- **Complex logic**: Undocumented algorithmic decisions

#### 3. Testing
- **No automated tests**: Only example data provided
- **No CI/CD**: Manual testing only
- **No regression tests**: Hard to verify changes don't break functionality

#### 4. Technical Gaps (from Field Analysis)
- **No BP distance validation**: Missing 8-30 nt constraint
- **No explicit PPT detection**: Implicit via 3'ss scoring
- **No score differential reporting**: Hard to interpret borderline cases
- **Single confidence level**: Binary classification only

#### 5. Performance
- **Not fully vectorized**: PWM scoring could be faster
- **No caching**: Repeated calculations (e.g., genome lookups)
- **Sequential bottlenecks**: Some parallelizable sections remain serial

### Strengths to Preserve

✅ **Proven SVM classification** (80/20 split, balanced weights)
✅ **Comprehensive metadata tracking**
✅ **Recursive training** for species adaptation
✅ **Non-canonical boundary adjustment**
✅ **Multi-format input/output**
✅ **Parallel processing** support

---

## Proposed Architecture

### New Directory Structure

```
intronIC/
├── intronIC/                      # Original (preserved)
│   ├── intronIC.py
│   ├── ...
│
├── intronIC_refactored/           # NEW: Refactored version
│   ├── __init__.py
│   ├── __main__.py                # Entry point for `python -m intronIC_refactored`
│   │
│   ├── cli/                       # Command-line interface
│   │   ├── __init__.py
│   │   ├── parser.py              # Argument parsing
│   │   └── main.py                # Main entry point
│   │
│   ├── core/                      # Core data models
│   │   ├── __init__.py
│   │   ├── models.py              # GenomeFeature, Gene, Transcript, Exon, Intron
│   │   ├── intron.py              # Intron class (separate due to size)
│   │   └── validators.py          # Validation functions (omit_check, etc.)
│   │
│   ├── io/                        # Input/Output handling
│   │   ├── __init__.py
│   │   ├── parsers.py             # Annotation, BED, sequence parsers
│   │   ├── genome.py              # Genome file handling
│   │   ├── writers.py             # Output file generation
│   │   └── formats.py             # File format specifications
│   │
│   ├── extraction/                # Intron extraction logic
│   │   ├── __init__.py
│   │   ├── annotator.py           # Annotation hierarchy parsing
│   │   ├── intronator.py          # Intron generation from exons
│   │   ├── sequences.py           # Sequence extraction
│   │   └── filters.py             # Filtering and tagging (duplicates, overlaps)
│   │
│   ├── scoring/                   # PWM and scoring
│   │   ├── __init__.py
│   │   ├── pwm.py                 # PWM loading, matrix operations
│   │   ├── scorer.py              # seq_score, bp_score functions
│   │   ├── features.py            # Feature extraction (5'ss, BP, 3'ss)
│   │   ├── normalization.py       # Z-score normalization
│   │   └── ppt.py                 # NEW: Polypyrimidine tract detection
│   │
│   ├── classification/            # ML classification
│   │   ├── __init__.py
│   │   ├── svm.py                 # SVM training and prediction
│   │   ├── optimization.py        # Hyperparameter optimization
│   │   ├── ensemble.py            # Ensemble voting
│   │   ├── recursive.py           # Recursive training
│   │   └── validators.py          # NEW: BP distance, confidence levels
│   │
│   ├── analysis/                  # Analysis and visualization
│   │   ├── __init__.py
│   │   ├── plotting.py            # Scatter plots, PR-AUC
│   │   └── statistics.py          # Summary statistics
│   │
│   ├── utils/                     # Utility functions
│   │   ├── __init__.py
│   │   ├── sequences.py           # Sequence manipulation (rev_comp, etc.)
│   │   ├── coordinates.py         # Coordinate system handling
│   │   ├── logging_config.py      # Logging configuration
│   │   └── helpers.py             # General utilities
│   │
│   ├── data/                      # Training data (symlink to original)
│   │   └── [symlink to ../intronIC/data/]
│   │
│   └── config/                    # NEW: Configuration management
│       ├── __init__.py
│       ├── defaults.py            # Default parameters
│       └── schemas.py             # Configuration validation
│
├── tests/                         # NEW: Test suite
│   ├── __init__.py
│   ├── conftest.py                # pytest configuration
│   ├── unit/                      # Unit tests
│   │   ├── test_models.py
│   │   ├── test_parsers.py
│   │   ├── test_scorer.py
│   │   ├── test_svm.py
│   │   └── ...
│   ├── integration/               # Integration tests
│   │   ├── test_extraction_pipeline.py
│   │   ├── test_scoring_pipeline.py
│   │   └── test_classification_pipeline.py
│   ├── regression/                # Regression tests
│   │   └── test_chr19_benchmark.py
│   └── fixtures/                  # Test data
│       ├── small_genome.fa
│       ├── small_annotation.gff3
│       └── expected_outputs/
│
├── docs/                          # NEW: Enhanced documentation
│   ├── API.md                     # API reference
│   ├── ALGORITHMS.md              # Algorithm details
│   ├── TUTORIAL.md                # User tutorial
│   └── CONTRIBUTING.md            # Development guide
│
├── benchmarks/                    # NEW: Performance benchmarks
│   ├── benchmark_scoring.py
│   ├── benchmark_extraction.py
│   └── results/
│
├── CLAUDE.md                      # Existing documentation
├── FIELD_APPROACHES.md
├── REFACTORING_PLAN.md            # This document
├── README.md
├── setup.py
├── setup.cfg
├── pyproject.toml                 # NEW: Modern Python packaging
└── .github/                       # NEW: CI/CD
    └── workflows/
        ├── tests.yml
        └── lint.yml
```

### Design Principles

1. **Single Responsibility**: Each module has one clear purpose
2. **Dependency Injection**: Pass dependencies explicitly (no global state)
3. **Interfaces over Implementation**: Use protocols/abstract base classes
4. **Immutability**: Prefer immutable data structures where possible
5. **Fail Fast**: Validate inputs early, provide clear error messages
6. **Type Safety**: Full type hint coverage with mypy validation

---

## Module Breakdown

### 1. Core Models (`core/`)

#### `core/models.py`
```python
"""Base genomic feature classes."""
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional, List

@dataclass(frozen=True)
class GenomicCoordinate:
    """Immutable genomic coordinate."""
    region: str
    start: int
    stop: int
    strand: str

    def __post_init__(self):
        if self.start > self.stop:
            raise ValueError(f"Start {self.start} > stop {self.stop}")

class GenomeFeature(ABC):
    """Base class for genomic features."""

    def __init__(
        self,
        coords: GenomicCoordinate,
        name: str,
        parent: Optional[str] = None
    ):
        self.coords = coords
        self.name = name
        self.parent = parent

    @property
    @abstractmethod
    def length(self) -> int:
        """Feature length."""
        pass

    @abstractmethod
    def get_seq(self, genome: 'GenomeReader', flank: int = 0) -> str:
        """Extract sequence from genome."""
        pass

class Gene(GenomeFeature):
    """Gene feature."""
    def __init__(self, coords, name, transcripts: List['Transcript']):
        super().__init__(coords, name)
        self.transcripts = transcripts

    @property
    def length(self) -> int:
        return self.coords.stop - self.coords.start

# ... similar for Transcript, Exon
```

#### `core/intron.py`
```python
"""Intron class with comprehensive attributes."""
from dataclasses import dataclass, field
from typing import Optional, Set, Tuple

@dataclass
class IntronScores:
    """Scoring data for an intron."""
    five_raw: Optional[float] = None
    five_z: Optional[float] = None
    bp_raw: Optional[float] = None
    bp_z: Optional[float] = None
    three_raw: Optional[float] = None
    three_z: Optional[float] = None
    svm_score: Optional[float] = None
    relative_score: Optional[float] = None

@dataclass
class IntronSequences:
    """Sequence data for an intron."""
    full: str = ""
    five_ss: str = ""
    bp: str = ""
    bp_region: str = ""
    three_ss: str = ""
    upstream_flank: str = ""
    downstream_flank: str = ""

@dataclass
class IntronMetadata:
    """Metadata for an intron."""
    index: Optional[int] = None
    family_size: Optional[int] = None
    phase: Optional[str] = None
    dnts: Optional[str] = None
    is_longest_isoform: bool = False
    is_noncanonical: bool = False
    is_duplicate: bool = False
    omission_reason: Optional[str] = None
    tags: Set[str] = field(default_factory=set)

class Intron(GenomeFeature):
    """Intron with comprehensive scoring and metadata."""

    def __init__(
        self,
        coords: GenomicCoordinate,
        name: str,
        parent: Optional[str] = None,
        grandparent: Optional[str] = None
    ):
        super().__init__(coords, name, parent)
        self.grandparent = grandparent

        # Use composition for related attributes
        self.scores = IntronScores()
        self.sequences = IntronSequences()
        self.metadata = IntronMetadata()

        # Coordinates for scoring regions
        self.five_coords: Optional[Tuple[int, int]] = None
        self.bp_coords: Optional[Tuple[int, int]] = None
        self.three_coords: Optional[Tuple[int, int]] = None

    @classmethod
    def from_exon_pair(
        cls,
        exon1: 'Exon',
        exon2: 'Exon',
        parent: Optional[str] = None
    ) -> 'Intron':
        """Factory method to create intron from two exons."""
        # Implementation...
        pass

    @property
    def length(self) -> int:
        return self.coords.stop - self.coords.start

    @property
    def type_id(self) -> str:
        """Classify as 'u2' or 'u12' based on score."""
        if self.scores.svm_score is None:
            return 'unknown'
        threshold = 90  # TODO: Make configurable
        return 'u12' if self.scores.svm_score > threshold else 'u2'
```

**Benefits:**
- Clear separation of concerns (coordinates, scores, sequences, metadata)
- Type safety with dataclasses
- Immutable coordinates prevent accidental modification
- Factory methods for complex construction

---

### 2. I/O Module (`io/`)

#### `io/parsers.py`
```python
"""Input file parsers."""
from typing import Iterator, Dict
from pathlib import Path
from biogl import GxfParse, fasta_parse

class AnnotationParser:
    """Parse GFF3/GTF annotation files."""

    def __init__(self, filepath: Path, feature_types: Tuple[str, ...] = ('cds', 'exon')):
        self.filepath = filepath
        self.feature_types = feature_types

    def parse(self) -> Iterator[Dict[str, Any]]:
        """Parse annotation into feature dictionaries."""
        with GxfParse(self.filepath) as parser:
            for feature in parser:
                if feature.type.lower() in self.feature_types:
                    yield self._feature_to_dict(feature)

    def _feature_to_dict(self, feature) -> Dict[str, Any]:
        """Convert GxfParse feature to dictionary."""
        # Implementation...
        pass

class BEDParser:
    """Parse BED format intron coordinates."""

    def __init__(self, filepath: Path):
        self.filepath = filepath

    def parse(self) -> Iterator[GenomicCoordinate]:
        """Parse BED file into coordinates."""
        with open(self.filepath) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                # BED is 0-based, half-open; convert to 1-based
                yield GenomicCoordinate(
                    region=fields[0],
                    start=int(fields[1]) + 1,  # Convert to 1-based
                    stop=int(fields[2]),
                    strand=fields[5] if len(fields) > 5 else '+'
                )

class SequenceParser:
    """Parse .iic sequence files."""

    def __init__(self, filepath: Path):
        self.filepath = filepath

    def parse(self) -> Iterator[Tuple[str, str, str, str]]:
        """Parse sequence file."""
        # name, 5'flank, intron, 3'flank
        # Implementation...
        pass
```

#### `io/genome.py`
```python
"""Genome file handling with caching."""
from typing import Dict, Optional
from pathlib import Path
from biogl import fasta_parse
from functools import lru_cache

class GenomeReader:
    """Cached genome sequence reader."""

    def __init__(self, filepath: Path, cache_size: int = 128):
        self.filepath = filepath
        self._sequences: Dict[str, str] = {}
        self._load_genome()

    def _load_genome(self):
        """Load genome into memory."""
        for header, seq in fasta_parse(self.filepath):
            chrom = header.split()[0]  # First word is chromosome
            self._sequences[chrom] = seq

    @lru_cache(maxsize=1024)
    def get_sequence(
        self,
        region: str,
        start: int,
        stop: int,
        strand: str = '+'
    ) -> str:
        """Get sequence with caching."""
        if region not in self._sequences:
            raise ValueError(f"Region {region} not in genome")

        seq = self._sequences[region][start-1:stop]  # Convert to 0-based

        if strand == '-':
            from intronIC_refactored.utils.sequences import reverse_complement
            seq = reverse_complement(seq)

        return seq

    def regions(self) -> List[str]:
        """Get list of chromosome/region names."""
        return list(self._sequences.keys())
```

**Benefits:**
- Clear separation of parsing logic
- Cached genome access reduces I/O
- Explicit coordinate system handling
- Easy to test with mock genomes

---

### 3. Scoring Module (`scoring/`)

#### `scoring/scorer.py`
```python
"""PWM scoring functions."""
from typing import Dict, Tuple, Optional, List
import numpy as np

PWMMatrix = Dict[str, Dict[int, float]]

def seq_score(
    seq: str,
    matrix: PWMMatrix,
    start_index: int = 0,
    ignore_positions: Optional[List[int]] = None,
    pseudocount: float = 0.0001
) -> float:
    """
    Score sequence using PWM.

    Args:
        seq: Sequence to score
        matrix: PWM matrix {base: {position: frequency}}
        start_index: Starting position in matrix
        ignore_positions: Positions to ignore (set to 1.0)
        pseudocount: Value added to avoid zero frequencies

    Returns:
        Product of position frequencies (raw score)

    Example:
        >>> matrix = {'A': {0: 0.8, 1: 0.2}, 'T': {0: 0.1, 1: 0.7}}
        >>> seq_score('AT', matrix)
        0.56
    """
    if ignore_positions is None:
        ignore_positions = []

    score = 1.0
    for i, base in enumerate(seq, start=start_index):
        if i in ignore_positions:
            continue

        try:
            freq = matrix[base][i]
        except KeyError:
            freq = pseudocount  # Handle N or other ambiguous bases

        score *= freq

    return score


def bp_score(
    seq: str,
    matrix: PWMMatrix,
    pseudocount: float = 0.0001
) -> Tuple[float, Tuple[int, int], str]:
    """
    Find best-scoring branch point motif in sequence.

    Uses sliding window to find optimal position.

    Args:
        seq: Branch point region sequence
        matrix: Branch point PWM
        pseudocount: Pseudocount for scoring

    Returns:
        Tuple of (best_score, (start, stop), best_sequence)

    Example:
        >>> seq = 'ACTTTCCTTAACTGAG'
        >>> score, coords, bp_seq = bp_score(seq, u12_bp_matrix)
        >>> bp_seq
        'TTCCTTAAC'
    """
    # Get matrix length (window size)
    matrix_length = len(next(iter(matrix.values())))

    best_score = -np.inf
    best_coords = (0, 0)
    best_seq = ""

    # Sliding window
    for start in range(len(seq) - matrix_length + 1):
        stop = start + matrix_length
        window = seq[start:stop]

        score = seq_score(window, matrix, pseudocount=pseudocount)

        if score > best_score:
            best_score = score
            best_coords = (start, stop)
            best_seq = window

    return best_score, best_coords, best_seq


def vectorized_score_batch(
    sequences: List[str],
    matrix: PWMMatrix,
    pseudocount: float = 0.0001
) -> np.ndarray:
    """
    Score multiple sequences efficiently using vectorization.

    Args:
        sequences: List of sequences (all same length)
        matrix: PWM matrix
        pseudocount: Pseudocount for scoring

    Returns:
        Array of scores

    Note:
        Much faster than looping for large batches.
    """
    # Convert sequences to numpy array of integers
    # A=0, C=1, G=2, T=3
    base_to_int = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    # Vectorized implementation
    # TODO: Implement for performance improvement
    pass
```

#### `scoring/ppt.py` (NEW)
```python
"""Polypyrimidine tract detection."""
from typing import Tuple

def calculate_ppt_strength(
    sequence: str,
    window: Tuple[int, int] = (-40, -10)
) -> float:
    """
    Calculate polypyrimidine tract strength.

    U12 introns typically have weak PPT (<0.6).
    U2 introns typically have strong PPT (>0.7).

    Args:
        sequence: Intron sequence
        window: Region to analyze (relative to 3'ss)

    Returns:
        Fraction of pyrimidines (C/T) in window

    Example:
        >>> seq = 'ACTG' * 20 + 'CCCTTTTCCCT' + 'AG'
        >>> calculate_ppt_strength(seq, window=(-15, -2))
        0.92  # Strong PPT (U2-like)
    """
    start, stop = window

    # Handle relative coordinates (negative from end)
    if start < 0:
        start = len(sequence) + start
    if stop < 0:
        stop = len(sequence) + stop

    region = sequence[start:stop]

    if len(region) == 0:
        return 0.0

    pyrimidines = region.count('C') + region.count('T')
    ppt_strength = pyrimidines / len(region)

    return ppt_strength


def classify_ppt(ppt_strength: float) -> str:
    """
    Classify PPT strength.

    Args:
        ppt_strength: PPT strength (0-1)

    Returns:
        'strong' (>0.7), 'weak' (<0.6), or 'moderate'
    """
    if ppt_strength > 0.7:
        return 'strong'  # U2-like
    elif ppt_strength < 0.6:
        return 'weak'    # U12-like
    else:
        return 'moderate'
```

**Benefits:**
- Clean, testable scoring functions
- Docstrings with examples
- Type hints for clarity
- NEW: PPT detection from field analysis
- Future: Vectorized batch scoring for performance

---

### 4. Classification Module (`classification/`)

#### `classification/validators.py` (NEW)
```python
"""Classification validation functions."""
from typing import Tuple, Optional
from ..core.intron import Intron

def validate_bp_distance(
    intron: Intron,
    valid_range: Tuple[int, int] = (8, 30)
) -> bool:
    """
    Validate branch point distance from 3' splice site.

    This is a key biological constraint from Sheth et al.
    U12 introns have BP 8-30 nt upstream of 3'ss.

    Args:
        intron: Intron object with BP coordinates
        valid_range: Acceptable distance range (nt)

    Returns:
        True if BP distance is valid

    Reference:
        Sheth et al. method - see FIELD_APPROACHES.md
    """
    if intron.bp_coords is None:
        return False

    # BP coords are relative to 3' splice site
    # Both should be negative (upstream)
    bp_start, bp_stop = intron.bp_coords

    # Distance is from end of BP to 3'ss
    distance_from_3ss = abs(bp_stop)

    min_dist, max_dist = valid_range
    return min_dist <= distance_from_3ss <= max_dist


def calculate_confidence_level(
    intron: Intron,
    threshold: float = 90.0,
    bp_distance_valid: Optional[bool] = None
) -> str:
    """
    Calculate multi-level confidence for classification.

    Args:
        intron: Classified intron
        threshold: Base threshold for classification
        bp_distance_valid: Whether BP distance is valid

    Returns:
        Confidence level: 'high', 'medium', 'low', 'borderline'

    Classification tiers:
        - high: score ≥ threshold+5 AND valid BP distance
        - medium: score ≥ threshold
        - borderline: score ≥ threshold-15 AND valid BP distance
        - low: score < threshold-15
    """
    score = intron.scores.svm_score

    if score is None:
        return 'unknown'

    if bp_distance_valid is None:
        bp_distance_valid = validate_bp_distance(intron)

    # High confidence: strong score + biological constraint
    if score >= threshold + 5 and bp_distance_valid:
        return 'high'

    # Medium confidence: exceeds threshold
    elif score >= threshold:
        return 'medium'

    # Borderline: below threshold but has valid BP distance
    elif score >= threshold - 15 and bp_distance_valid:
        return 'borderline'

    # Low confidence
    else:
        return 'low'


def calculate_score_differentials(intron: Intron) -> Dict[str, float]:
    """
    Calculate U12 vs U2 score differentials for each region.

    Useful for interpreting borderline classifications and
    comparing to classical thresholds (Sheth et al.: 25, 10).

    Args:
        intron: Intron with U12 and U2 scores

    Returns:
        Dictionary of differentials for each region

    Example:
        >>> diffs = calculate_score_differentials(intron)
        >>> diffs['five_ss']
        32.5  # Strong preference for U12 5'ss (> 25 threshold)
    """
    # TODO: Requires storing both U12 and U2 raw scores
    # Currently only stores final scores

    return {
        'five_ss_diff': 0.0,  # u12_five - u2_five
        'bp_diff': 0.0,
        'three_ss_diff': 0.0
    }
```

#### `classification/svm.py`
```python
"""SVM classification module."""
from typing import List, Tuple, Optional
import numpy as np
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score, precision_recall_curve, auc
from dataclasses import dataclass

@dataclass
class SVMTrainingResult:
    """Results from SVM training."""
    model: svm.SVC
    f1_score: float
    pr_auc: float
    train_accuracy: float
    test_accuracy: float
    optimal_c: float

class IntronClassifier:
    """SVM-based intron classifier."""

    def __init__(
        self,
        kernel: str = 'linear',
        class_weight: str = 'balanced',
        probability: bool = True,
        random_state: int = 42
    ):
        self.kernel = kernel
        self.class_weight = class_weight
        self.probability = probability
        self.random_state = random_state
        self.model: Optional[svm.SVC] = None
        self.scaler = None

    def train(
        self,
        X_train: np.ndarray,
        y_train: np.ndarray,
        C: float = 1.0,
        test_size: float = 0.2
    ) -> SVMTrainingResult:
        """
        Train SVM classifier.

        Args:
            X_train: Feature matrix (n_samples, n_features)
            y_train: Labels (0=U2, 1=U12)
            C: Regularization parameter
            test_size: Fraction for test set

        Returns:
            Training results with metrics
        """
        # Split train/test
        X_tr, X_te, y_tr, y_te = train_test_split(
            X_train, y_train,
            test_size=test_size,
            stratify=y_train,
            random_state=self.random_state
        )

        # Train model
        self.model = svm.SVC(
            kernel=self.kernel,
            C=C,
            class_weight=self.class_weight,
            probability=self.probability,
            random_state=self.random_state
        )
        self.model.fit(X_tr, y_tr)

        # Evaluate
        y_pred_train = self.model.predict(X_tr)
        y_pred_test = self.model.predict(X_te)

        train_acc = (y_pred_train == y_tr).mean()
        test_acc = (y_pred_test == y_te).mean()
        f1 = f1_score(y_te, y_pred_test)

        # PR-AUC
        y_proba = self.model.predict_proba(X_te)[:, 1]
        precision, recall, _ = precision_recall_curve(y_te, y_proba)
        pr_auc_score = auc(recall, precision)

        return SVMTrainingResult(
            model=self.model,
            f1_score=f1,
            pr_auc=pr_auc_score,
            train_accuracy=train_acc,
            test_accuracy=test_acc,
            optimal_c=C
        )

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        """Predict U12 probabilities."""
        if self.model is None:
            raise ValueError("Model not trained")
        return self.model.predict_proba(X)[:, 1]  # P(U12)

    def predict(self, X: np.ndarray, threshold: float = 0.9) -> np.ndarray:
        """Predict class labels with custom threshold."""
        probas = self.predict_proba(X)
        return (probas > threshold).astype(int)
```

**Benefits:**
- NEW: BP distance validation (field-proven constraint)
- NEW: Multi-level confidence classification
- NEW: Score differentials for interpretability
- Clean SVM interface with dataclass results
- Easy to test and extend

---

### 5. Configuration Module (`config/`) (NEW)

#### `config/defaults.py`
```python
"""Default configuration values."""
from dataclasses import dataclass, field
from typing import Tuple, Optional

@dataclass
class ScoringRegions:
    """Coordinates for scoring regions."""
    five_ss: Tuple[int, int] = (-3, 9)
    branch_point: Tuple[int, int] = (-55, -5)
    three_ss: Tuple[int, int] = (-6, 4)

@dataclass
class ClassificationConfig:
    """Classification parameters."""
    threshold: float = 90.0
    validate_bp_distance: bool = True
    bp_distance_range: Tuple[int, int] = (8, 30)
    detect_ppt: bool = True
    ppt_window: Tuple[int, int] = (-40, -10)

@dataclass
class SVMConfig:
    """SVM training configuration."""
    kernel: str = 'linear'
    class_weight: str = 'balanced'
    n_optimization_rounds: int = 5
    cv_folds: int = 5
    test_size: float = 0.2
    random_state: int = 42

@dataclass
class FilteringConfig:
    """Intron filtering parameters."""
    min_length: int = 30
    allow_noncanonical: bool = True
    longest_isoform_only: bool = True
    allow_overlap: bool = False
    adjust_nc_boundaries: bool = True

@dataclass
class IntronICConfig:
    """Complete intronIC configuration."""
    scoring_regions: ScoringRegions = field(default_factory=ScoringRegions)
    classification: ClassificationConfig = field(default_factory=ClassificationConfig)
    svm: SVMConfig = field(default_factory=SVMConfig)
    filtering: FilteringConfig = field(default_factory=FilteringConfig)

    # I/O
    intron_flank_size: int = 200
    pseudocount: float = 0.0001

    # Performance
    n_processes: int = 1
    cv_processes: Optional[int] = None

    # Output
    generate_plots: bool = True
    output_sequences: bool = True
    include_duplicates: bool = False
    fig_dpi: int = 300
```

**Benefits:**
- Type-safe configuration
- Clear defaults in code
- Easy to serialize to/from YAML/JSON
- Validation via dataclass

---

## Technical Improvements

### 1. Branch Point Distance Validation ⭐⭐⭐

**Priority:** HIGH
**Effort:** Low (1-2 days)
**Impact:** High (biological constraint)

**Implementation:**
- Module: `classification/validators.py`
- Function: `validate_bp_distance()`
- Add to classification pipeline after SVM scoring
- Use as additional filter or confidence boost

**Code:**
```python
# In classification pipeline
for intron in introns:
    # SVM classification
    intron.scores.svm_score = classifier.predict_proba(features)[0] * 100

    # NEW: Validate BP distance
    bp_valid = validate_bp_distance(intron)
    intron.metadata.tags.add('bp_distance_valid' if bp_valid else 'bp_distance_invalid')

    # NEW: Multi-level confidence
    confidence = calculate_confidence_level(intron, threshold=90.0, bp_distance_valid=bp_valid)
    intron.metadata.tags.add(f'confidence_{confidence}')
```

### 2. Polypyrimidine Tract Detection ⭐⭐

**Priority:** HIGH
**Effort:** Low (1-2 days)
**Impact:** Medium (additional feature)

**Implementation:**
- Module: `scoring/ppt.py`
- Function: `calculate_ppt_strength()`
- Add during feature extraction
- Option 1: Use as 4th SVM feature
- Option 2: Use as filter/validator

**Code:**
```python
# In feature extraction
from intronIC_refactored.scoring.ppt import calculate_ppt_strength

for intron in introns:
    ppt_strength = calculate_ppt_strength(intron.sequences.full)
    intron.metadata.ppt_strength = ppt_strength

    # Optional: Add as SVM feature
    features.append([
        intron.scores.five_z,
        intron.scores.bp_z,
        intron.scores.three_z,
        ppt_strength  # NEW: 4th dimension
    ])
```

### 3. Score Differential Reporting ⭐⭐⭐

**Priority:** HIGH
**Effort:** Low (1 day)
**Impact:** Medium (interpretability)

**Implementation:**
- Requires storing both U12 and U2 raw scores (not just final)
- Add to `IntronScores` dataclass
- Report in `.scores.iic` output
- Allows comparison to classical thresholds (25, 10)

**Code:**
```python
@dataclass
class IntronScores:
    # Existing
    five_raw: Optional[float] = None
    # ...

    # NEW: Store U12 and U2 separately
    five_raw_u12: Optional[float] = None
    five_raw_u2: Optional[float] = None
    bp_raw_u12: Optional[float] = None
    bp_raw_u2: Optional[float] = None
    three_raw_u12: Optional[float] = None
    three_raw_u2: Optional[float] = None

    @property
    def five_differential(self) -> float:
        """U12 - U2 score for 5'ss."""
        if self.five_raw_u12 is None or self.five_raw_u2 is None:
            return 0.0
        return self.five_raw_u12 - self.five_raw_u2
```

### 4. Vectorized PWM Scoring ⭐⭐

**Priority:** MEDIUM
**Effort:** Medium (3-5 days)
**Impact:** High (performance)

**Implementation:**
- Module: `scoring/scorer.py`
- Function: `vectorized_score_batch()`
- Use NumPy for batch processing
- Expected speedup: 2-10x for scoring phase

**Pseudocode:**
```python
def vectorized_score_batch(sequences: List[str], matrix: PWMMatrix) -> np.ndarray:
    """Score batch of sequences efficiently."""
    # Convert sequences to integer array
    # Shape: (n_sequences, seq_length)
    seq_array = sequences_to_array(sequences)

    # Convert matrix to array
    # Shape: (4, seq_length) for ACGT
    matrix_array = matrix_to_array(matrix)

    # Fancy indexing to get frequencies
    # Shape: (n_sequences, seq_length)
    frequencies = matrix_array[seq_array, np.arange(seq_array.shape[1])]

    # Product along sequence dimension
    scores = np.prod(frequencies, axis=1)

    return scores
```

### 5. Configuration File Support ⭐

**Priority:** MEDIUM
**Effort:** Medium (2-3 days)
**Impact:** Medium (usability)

**Implementation:**
```python
# config.yaml
classification:
  threshold: 95.0  # Stricter than default
  validate_bp_distance: true

scoring_regions:
  five_ss: [-5, 10]  # Custom coordinates

svm:
  n_optimization_rounds: 7
  cv_folds: 10
```

**Usage:**
```bash
intronIC_refactored --config config.yaml -g genome.fa -a annot.gff3 -n species
```

### 6. AT-AC Specific Handling ⭐

**Priority:** LOW
**Effort:** Low (1 day)
**Impact:** Low (edge case improvement)

**Implementation:**
```python
def classify_atac_intron(intron: Intron) -> str:
    """AT-AC introns are almost exclusively U12."""
    if intron.metadata.dnts == 'AT-AC':
        # Strong prior for U12
        if intron.scores.svm_score > 50:  # Lower threshold
            intron.metadata.tags.add('AT-AC_boosted')
            return 'u12'
    return classify_standard(intron)
```

---

## Testing Strategy

### Test Structure

```
tests/
├── unit/                          # Fast, isolated tests
│   ├── test_models.py             # Data model tests
│   ├── test_parsers.py            # Parser tests
│   ├── test_scorer.py             # PWM scoring tests
│   ├── test_validators.py         # NEW: BP distance, PPT tests
│   └── test_svm.py                # Classifier tests
│
├── integration/                   # Multi-component tests
│   ├── test_extraction_pipeline.py
│   ├── test_scoring_pipeline.py
│   └── test_full_pipeline.py
│
├── regression/                    # Accuracy preservation
│   ├── test_chr19_benchmark.py   # Match original results
│   └── test_known_u12s.py        # Validate known U12 detection
│
└── fixtures/                      # Test data
    ├── small_genome.fa           # 100kb test genome
    ├── small_annotation.gff3     # 50 genes
    ├── known_u12_introns.bed     # Validated U12s
    └── expected_outputs/         # Expected results
```

### Key Test Cases

#### Unit Tests

```python
# tests/unit/test_scorer.py
def test_seq_score_basic():
    """Test basic PWM scoring."""
    matrix = {
        'A': {0: 0.8, 1: 0.1},
        'T': {0: 0.1, 1: 0.8},
        'G': {0: 0.05, 1: 0.05},
        'C': {0: 0.05, 1: 0.05}
    }

    score = seq_score('AT', matrix)
    assert score == pytest.approx(0.64, rel=1e-5)

def test_bp_score_finds_best_motif():
    """Test branch point motif search."""
    # TTCCTTAAC is U12 BP consensus
    seq = 'ACTGACTGTTCCTTAACTGAG'

    score, coords, motif = bp_score(seq, u12_bp_matrix)

    assert 'TTCCTTAAC' in motif
    assert coords[0] >= 8
    assert coords[1] <= 30

# tests/unit/test_validators.py (NEW)
def test_bp_distance_validation():
    """Test branch point distance constraint."""
    intron = create_test_intron()
    intron.bp_coords = (0, -25)  # 25 nt upstream

    assert validate_bp_distance(intron) == True

    intron.bp_coords = (0, -5)  # Too close
    assert validate_bp_distance(intron) == False

def test_ppt_detection():
    """Test polypyrimidine tract detection."""
    # Strong PPT (U2-like)
    u2_seq = 'ACTG' * 10 + 'CCCTTTTCCCTTT' + 'AG'
    assert calculate_ppt_strength(u2_seq) > 0.7

    # Weak PPT (U12-like)
    u12_seq = 'ACTG' * 10 + 'ACGTACGT' + 'AG'
    assert calculate_ppt_strength(u12_seq) < 0.6
```

#### Integration Tests

```python
# tests/integration/test_extraction_pipeline.py
def test_annotation_to_introns():
    """Test full extraction pipeline."""
    genome = GenomeReader('tests/fixtures/small_genome.fa')
    parser = AnnotationParser('tests/fixtures/small_annotation.gff3')

    hierarchy = build_annotation_hierarchy(parser.parse())
    introns = extract_introns(hierarchy, genome)

    assert len(introns) > 0
    assert all(i.length >= 30 for i in introns)
    assert all(i.sequences.full != "" for i in introns)
```

#### Regression Tests

```python
# tests/regression/test_chr19_benchmark.py
def test_chr19_matches_original():
    """Verify refactored version matches original results."""
    original_results = load_original_results('chr19_original.iic')

    # Run refactored version
    config = IntronICConfig()
    results = run_intronic_pipeline(
        genome='test_data/Homo_sapiens.Chr19.fa.gz',
        annotation='test_data/Homo_sapiens.Chr19.gff3.gz',
        config=config
    )

    # Compare
    assert len(results) == len(original_results)

    for orig, new in zip(original_results, results):
        # Allow small numerical differences
        assert orig.coords == new.coords
        assert orig.type_id == new.type_id
        assert abs(orig.svm_score - new.svm_score) < 1.0
```

### Test Coverage Goals

- **Overall:** >85% line coverage
- **Core modules:** >95% coverage
- **Critical paths:** 100% coverage
  - Coordinate handling
  - PWM scoring
  - SVM classification
  - BP distance validation

### Continuous Integration

```yaml
# .github/workflows/tests.yml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: [3.8, 3.9, 3.10, 3.11]

    steps:
    - uses: actions/checkout@v2
    - name: Set up Python ${{ matrix.python-version }}
      uses: actions/setup-python@v2
      with:
        python-version: ${{ matrix.python-version }}

    - name: Install dependencies
      run: |
        pip install -e .[dev]
        pip install pytest pytest-cov mypy black

    - name: Run tests
      run: pytest tests/ --cov=intronIC_refactored --cov-report=xml

    - name: Type check
      run: mypy intronIC_refactored/

    - name: Code style
      run: black --check intronIC_refactored/
```

---

## Implementation Roadmap

### Phase 1: Foundation (2-3 weeks)

**Week 1: Core Infrastructure**
- [ ] Set up directory structure
- [ ] Create `core/models.py` with base classes
- [ ] Create `core/intron.py` with refactored Intron class
- [ ] Create `utils/` module with coordinate handling
- [ ] Set up pytest framework
- [ ] Write unit tests for models

**Week 2: I/O Layer**
- [ ] Create `io/parsers.py` (annotation, BED, sequence)
- [ ] Create `io/genome.py` with caching
- [ ] Create `io/writers.py` for output
- [ ] Write unit tests for parsers
- [ ] Test with chr19 data

**Week 3: Extraction Pipeline**
- [ ] Create `extraction/annotator.py`
- [ ] Create `extraction/intronator.py`
- [ ] Create `extraction/sequences.py`
- [ ] Create `extraction/filters.py`
- [ ] Integration tests for extraction
- [ ] Verify matches original output

### Phase 2: Scoring & Classification (2-3 weeks)

**Week 4: Scoring**
- [ ] Create `scoring/pwm.py` (matrix loading)
- [ ] Create `scoring/scorer.py` (seq_score, bp_score)
- [ ] Create `scoring/features.py` (feature extraction)
- [ ] Create `scoring/normalization.py` (z-scores)
- [ ] **NEW:** Create `scoring/ppt.py` (PPT detection)
- [ ] Unit tests for all scoring functions

**Week 5: Classification**
- [ ] Create `classification/svm.py` (training & prediction)
- [ ] Create `classification/optimization.py` (hyperparameter tuning)
- [ ] **NEW:** Create `classification/validators.py` (BP distance, confidence)
- [ ] Create `classification/ensemble.py` (multi-model voting)
- [ ] Unit tests for SVM

**Week 6: Advanced Features**
- [ ] Create `classification/recursive.py` (recursive training)
- [ ] Create `analysis/plotting.py` (visualizations)
- [ ] Create `analysis/statistics.py` (summary stats)
- [ ] Integration tests for full pipeline
- [ ] Regression tests vs original

### Phase 3: CLI & Configuration (1 week)

**Week 7: User Interface**
- [ ] Create `cli/parser.py` (argument parsing)
- [ ] Create `cli/main.py` (entry point)
- [ ] Create `config/` module (defaults, schemas)
- [ ] Add YAML config file support
- [ ] Update documentation
- [ ] End-to-end testing

### Phase 4: Optimization & Polish (2 weeks)

**Week 8: Performance**
- [ ] Implement vectorized PWM scoring
- [ ] Profile and optimize bottlenecks
- [ ] Add caching where beneficial
- [ ] Parallel processing verification
- [ ] Benchmark vs original

**Week 9: Documentation & Release**
- [ ] Write API documentation
- [ ] Create tutorial
- [ ] Write migration guide
- [ ] Set up CI/CD
- [ ] Create release checklist
- [ ] Alpha release for testing

### Total Timeline: 8-10 weeks

---

## Migration Strategy

### Backward Compatibility

**Preserve CLI Interface:**
```bash
# Original command continues to work
intronIC -g genome.fa -a annot.gff3 -n species -t 90

# New command uses refactored version
intronIC_refactored -g genome.fa -a annot.gff3 -n species -t 90

# Or via Python module
python -m intronIC_refactored -g genome.fa -a annot.gff3 -n species
```

**Output Compatibility:**
- Generate same output file formats
- Validate numerical agreement with original
- Document any intentional changes (e.g., added columns)

### Phased Rollout

**Phase 1: Parallel Development**
- Keep original `intronIC/` untouched
- Develop refactored version in `intronIC_refactored/`
- Both versions coexist

**Phase 2: Alpha Testing**
- Internal testing on diverse datasets
- Compare results to original
- Fix discrepancies
- Performance benchmarking

**Phase 3: Beta Release**
- Public beta release
- Community testing
- Bug fixes
- Documentation improvements

**Phase 4: Stable Release**
- Version 2.0 release
- Original becomes "legacy"
- Support both for transition period
- Eventually deprecate original

### Data Migration

**No data migration needed** - both versions use same:
- Reference sequences (`data/`)
- Test data (`test_data/`)
- Input file formats

Simply symlink data directory:
```bash
cd intronIC_refactored/
ln -s ../intronIC/data/ data
```

---

## Risk Mitigation

### Technical Risks

**Risk 1: Performance Regression**
- **Mitigation:** Benchmark each module during development
- **Fallback:** Profile and optimize before release
- **Acceptance:** No >10% slowdown

**Risk 2: Numerical Differences**
- **Mitigation:** Regression tests comparing to original
- **Fallback:** Investigate and document differences
- **Acceptance:** Differences must be explained and justified

**Risk 3: Bug Introduction**
- **Mitigation:** Comprehensive test suite (>85% coverage)
- **Fallback:** Beta testing period
- **Acceptance:** No critical bugs in core functionality

### Project Risks

**Risk 1: Scope Creep**
- **Mitigation:** Stick to defined roadmap
- **Fallback:** Defer nice-to-have features to v2.1
- **Acceptance:** Complete core features first

**Risk 2: Timeline Overrun**
- **Mitigation:** Regular progress tracking
- **Fallback:** Reduce scope if needed
- **Acceptance:** Minimum viable product: working pipeline with tests

**Risk 3: Breaking Changes**
- **Mitigation:** Maintain backward compatibility
- **Fallback:** Document migration path
- **Acceptance:** Users can use original until ready to migrate

---

## Success Criteria

### Must Have (v2.0)
- ✅ Modular architecture (<500 lines/module)
- ✅ Full type hints
- ✅ >85% test coverage
- ✅ BP distance validation
- ✅ PPT detection
- ✅ Multi-level confidence
- ✅ Matches original accuracy
- ✅ No performance regression

### Should Have (v2.0)
- ✅ Score differential reporting
- ✅ Configuration file support
- ✅ Improved documentation
- ✅ CI/CD pipeline
- ✅ Vectorized scoring (optional)

### Nice to Have (v2.1+)
- ⭐ Deep learning classifier option
- ⭐ Interactive visualization
- ⭐ REST API
- ⭐ Web interface
- ⭐ Database backend

---

## Appendix: Code Comparison

### Before (Original)

```python
# intronIC.py - Line 2114, monolithic file
def seq_score(seq, matrix, start_index=0, ignore=None):
    """
    Score {seq} using values from {matrix}.
    """
    score = None
    for i, e in enumerate(seq, start=start_index):
        if ignore is not None and i in ignore:
            base_freq = 1.0
        else:
            try:
                base_freq = matrix[e][i]
            except KeyError:
                base_freq = 0.0001
        if score is None:
            score = base_freq
        else:
            score *= base_freq
    return score
```

### After (Refactored)

```python
# intronIC_refactored/scoring/scorer.py - focused module
from typing import Dict, Optional, List

PWMMatrix = Dict[str, Dict[int, float]]

def seq_score(
    seq: str,
    matrix: PWMMatrix,
    start_index: int = 0,
    ignore_positions: Optional[List[int]] = None,
    pseudocount: float = 0.0001
) -> float:
    """
    Score sequence using position weight matrix.

    Computes the product of position-specific base frequencies
    for each base in the sequence. Used for evaluating splice
    site strength against U2 or U12 consensus sequences.

    Args:
        seq: Nucleotide sequence to score
        matrix: PWM as {base: {position: frequency}}
        start_index: Starting position in matrix (default: 0)
        ignore_positions: Positions to ignore in scoring
        pseudocount: Value for missing/ambiguous bases

    Returns:
        Raw PWM score (product of frequencies)

    Example:
        >>> matrix = {'A': {0: 0.9}, 'T': {0: 0.1}}
        >>> seq_score('A', matrix)
        0.9

    Note:
        Scores are multiplicative; longer sequences will have
        lower absolute scores. Use z-score normalization for
        comparison across different sequence lengths.
    """
    if ignore_positions is None:
        ignore_positions = []

    score = 1.0

    for i, base in enumerate(seq, start=start_index):
        if i in ignore_positions:
            continue

        try:
            freq = matrix[base][i]
        except KeyError:
            # Handle N or other ambiguous bases
            freq = pseudocount

        score *= freq

    return score
```

**Improvements:**
- Type hints for all parameters and return value
- Comprehensive docstring with examples
- Explicit parameter naming
- Better variable names
- Module organization

---

## Next Steps

1. **Review this plan** with stakeholders
2. **Prioritize features** based on user needs
3. **Set up development environment** (fork, branch strategy)
4. **Begin Phase 1** (Foundation) with core models
5. **Regular progress reviews** (weekly syncs)

---

**Document Version:** 1.0
**Author:** Claude Code Development Team
**Status:** Proposal - Awaiting Approval
**Last Updated:** 2025-11-01
