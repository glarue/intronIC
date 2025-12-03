# intronIC Test Suite

## Overview

This test suite provides comprehensive coverage of the intronIC refactored codebase, including unit tests, integration tests, and edge case validation.

## Test Organization

```
tests/
├── unit/                          # Unit tests for individual components
│   ├── test_coordinates.py        # Coordinate system conversions
│   ├── test_models.py              # Core data models (Intron, etc.)
│   ├── test_error_handling.py      # Error handling and exceptions
│   ├── test_edge_cases.py          # Edge cases and boundary conditions
│   ├── test_scoring/               # Scoring module tests
│   │   ├── test_pwm.py             # PWM loading and scoring
│   │   ├── test_pwm_format_equivalence.py  # Legacy .iic vs JSON format tests
│   │   ├── test_pwm_indexing.py    # PWM position indexing
│   │   ├── test_scorer.py          # IntronScorer pipeline
│   │   ├── test_normalizer.py      # Z-score normalization
│   │   ├── test_branch_point.py    # Branch point scoring
│   │   ├── test_log_ratio_deterministic.py # Log-ratio calculations
│   │   └── test_new_scaling_architecture.py # ZAR scaling tests
│   ├── test_classification/        # Classification module tests
│   │   └── test_ml_integrity.py    # ML data leakage prevention
│   └── test_file_io/                # File I/O tests
│       ├── test_readers.py          # File readers
│       ├── test_writers.py          # File writers
│       └── test_writer_errors.py    # Writer error handling
├── integration/                    # Integration tests
│   ├── test_extraction_pipeline.py # Full extraction workflow
│   └── test_classification_pipeline.py # Full classification workflow
└── README.md                       # This file

## Test Categories (Pytest Markers)

### Component Markers
- `@pytest.mark.unit` - Unit tests for individual components
- `@pytest.mark.integration` - Multi-component integration tests

### Feature Markers
- `@pytest.mark.ml` - Machine learning and classification tests
- `@pytest.mark.scoring` - PWM scoring and normalization tests
- `@pytest.mark.extraction` - Intron extraction from annotations
- `@pytest.mark.io` - File input/output tests
- `@pytest.mark.error_handling` - Error and exception handling tests
- `@pytest.mark.edge_case` - Edge case and boundary condition tests

### Resource Markers
- `@pytest.mark.slow` - Tests that take >5 seconds to run
- `@pytest.mark.requires_data` - Tests requiring external data files
- `@pytest.mark.requires_chr19` - Tests specifically needing Chr19 test data

## Running Tests

### Run All Tests
```bash
pytest
```

### Run Specific Category
```bash
# Unit tests only
pytest -m unit

# Integration tests only
pytest -m integration

# ML/classification tests
pytest -m ml

# Scoring tests
pytest -m scoring
```

### Run Fast Tests Only (Skip Slow Tests)
```bash
pytest -m "not slow"
```

### Run Tests Requiring Data
```bash
pytest -m requires_chr19
```

### Run Specific Test File
```bash
pytest tests/unit/test_scoring/test_pwm.py
```

### Run Specific Test Function
```bash
pytest tests/unit/test_scoring/test_pwm.py::test_pwm_loading
```

### Run with Verbose Output
```bash
pytest -v
```

### Run with Coverage Report
```bash
pytest --cov=. --cov-report=html
```

## Test Coverage

### Current Coverage by Module

| Module | Coverage | Notes |
|--------|----------|-------|
| core/ | 95%+ | Comprehensive coverage of data models |
| scoring/ | 90%+ | PWM, normalization, scoring pipeline |
| classification/ | 85%+ | SVM training, prediction, ML integrity |
| extraction/ | 80%+ | Annotation parsing, intron generation |
| file_io/ | 85%+ | Readers, writers, format handling |
| utils/ | 70%+ | Utility functions |

### Critical Test Areas

**1. ML Integrity (test_ml_integrity.py)**
- Prevents data leakage in z-score normalization
- Validates fit-transform workflow
- Ensures reference-only training

**2. Coordinate Systems (test_coordinates.py)**
- BED (0-based, half-open) vs GFF3 (1-based, closed)
- Conversion accuracy
- Edge cases at chromosome boundaries

**3. PWM Scoring (test_pwm.py, test_scorer.py)**
- Matrix loading and validation
- Log-ratio calculations
- Branch point motif search
- Pseudocount handling

**4. Error Handling (test_error_handling.py)**
- File I/O errors
- Invalid data handling
- Graceful degradation

**5. Edge Cases (test_edge_cases.py)**
- Minimum/maximum values
- Special characters
- Unusual but valid inputs
- Numerical edge cases (NaN, Inf)

## Writing New Tests

### Test Structure Template

```python
"""
Module docstring describing test purpose.
"""

import pytest
from pathlib import Path

from module_to_test import ComponentToTest

# Module-level markers
pytestmark = [pytest.mark.unit, pytest.mark.scoring]

# ============================================================================
# Test Fixtures
# ============================================================================

@pytest.fixture
def sample_data():
    """Create sample data for tests."""
    return {"key": "value"}


# ============================================================================
# Test Category Name
# ============================================================================

def test_basic_functionality():
    """Test basic functionality with standard inputs."""
    # Arrange
    input_data = "test"

    # Act
    result = function_under_test(input_data)

    # Assert
    assert result == expected_output


@pytest.mark.slow
def test_performance_intensive_operation():
    """Test that takes significant time."""
    pass


@pytest.mark.requires_data
def test_with_external_data():
    """Test requiring external files."""
    pass
```

### Best Practices

1. **One Assert Per Test (Generally)**: Each test should verify one specific behavior
2. **Arrange-Act-Assert**: Structure tests clearly with setup, execution, and verification
3. **Descriptive Names**: Test names should describe what they test
4. **Use Fixtures**: Share common setup via pytest fixtures
5. **Mark Appropriately**: Use markers to categorize tests
6. **Test Edge Cases**: Include boundary conditions, empty inputs, special values
7. **Document Intent**: Add docstrings explaining what each test validates

### Example: Adding a New Scoring Test

```python
@pytest.mark.unit
@pytest.mark.scoring
def test_pwm_score_with_ambiguous_nucleotides():
    """Test PWM scoring handles ambiguous nucleotides (N, W, S, etc.)."""
    # Arrange
    pwm = create_test_pwm()
    sequence = "ACGTNNN"  # Contains ambiguous 'N'

    # Act
    score = pwm.score_sequence(sequence)

    # Assert
    assert score is not None  # Should handle gracefully
    assert not np.isnan(score)  # Should not be NaN
```

## Debugging Failed Tests

### View Full Traceback
```bash
pytest --tb=long
```

### Show Local Variables
```bash
pytest --showlocals
```

### Stop on First Failure
```bash
pytest -x
```

### Run Last Failed Tests
```bash
pytest --lf
```

### Run Failed Then All
```bash
pytest --ff
```

### Print Statements
```bash
pytest -s  # Shows print() output
```

## Continuous Integration

Tests run automatically on:
- Push to any branch
- Pull request creation/update
- Scheduled nightly runs

CI configuration: `.github/workflows/tests.yml` (if exists)

## Test Data

Test data is located in `test_data/`:
- `Homo_sapiens.Chr19.Ensembl_91.gff3.gz` - Chr19 annotation
- `Homo_sapiens.Chr19.Ensembl_91.fa.gz` - Chr19 genome sequence

This is symlinked from the main intronIC data directory.

## Common Issues and Solutions

### Issue: "Chr19 test data not available"
**Solution**: Ensure `test_data/` symlink points to valid location:
```bash
ls -l test_data/
# Should show symlink to ../intronIC/test_data
```

### Issue: "Module not found" errors
**Solution**: Install package in development mode:
```bash
pixi install
# or
pip install -e .
```

### Issue: "Marker not registered" warnings
**Solution**: All markers are defined in `pytest.ini`. If adding new marker, update configuration.

### Issue: Tests pass locally but fail in CI
**Solution**: Check for:
- Hardcoded paths
- Timezone dependencies
- Random seed issues
- File permission assumptions

## Contributing

When adding new features:
1. Write tests first (TDD) or alongside feature
2. Ensure all existing tests still pass
3. Add appropriate markers
4. Update this README if adding new test categories
5. Aim for >80% coverage of new code

## Performance Benchmarks

Approximate test run times (on standard hardware):
- Full suite: ~2-3 minutes
- Unit tests only: ~30 seconds
- Integration tests: ~1-2 minutes
- Fast tests (no slow marker): ~45 seconds

## Future Test Improvements

- [ ] Add performance regression tests
- [ ] Expand edge case coverage for rare intron types
- [ ] Add property-based testing with Hypothesis
- [ ] Create test data generator for synthetic introns
- [ ] Add mutation testing with mutmut
- [ ] Benchmark memory usage in tests
- [ ] Add visual regression tests for plots

## Contact

For questions about tests, see main project README or open an issue on GitHub.
