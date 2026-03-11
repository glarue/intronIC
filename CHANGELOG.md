# Changelog

All notable changes to intronIC will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.1.3] - 2026-03-10

### Added
- New `NO_SEQUENCE` omission reason (`[x]` tag) for introns on genome regions missing from the FASTA (e.g. organellar genes referencing contigs not included in nuclear assemblies)

### Fixed
- Pipeline no longer crashes when introns reference contigs absent from the genome FASTA
  - `SequenceExtractor` now yields these introns marked as `omitted_no_sequence` instead of silently dropping them
  - `SequenceWriter` skips introns without sequence data instead of raising `ValueError`
  - Affected species include those with mitochondrial/plastid gene annotations (e.g. *Manihot esculenta*, *Protopterus annectens*, *Spinacia oleracea*)

### Dependencies
- Bump pillow from 12.0.0 to 12.1.1

## [2.1.2] - 2026-02-15

### Fixed
- Log files no longer contain ANSI color codes
- Disabled line wrapping in log file output

## [2.1.0] - 2024-12-15

### Added
- New `intronIC test` subcommand for quick installation validation
  - Runs classification on bundled Chr19 test data
  - Displays intron counts and classification results
  - `--show-only` flag to show test data location without running
  - `--output-dir` flag to specify output location

### Changed
- **Default model switched to isotonic calibration** for improved cross-species performance
  - Better classification accuracy in C. elegans and other non-human species
  - Previous sigmoid-calibrated model preserved as `default_pretrained.model.sigmoid.pkl`
- Removed Python version upper bounds (now compatible with Python 3.14+)
- Removed numpy <2.0 constraint (verified compatibility in commit 0525903)

### Fixed
- `intronIC test` now displays correct intron counts from metrics file
  - Reads from `.metrics.iic.json` instead of looking for non-existent JSON
  - Shows `total_scored` (classified introns) instead of `total_introns_generated`
  - Shows `high_confidence_u12` to match classification results table
  - Added thousand separators for better readability
- Fixed UnboundLocalError in test command by initializing variables before loop
- Completely suppressed sklearn InconsistentVersionWarning
  - All `joblib.load()` calls replaced with `load_model()` wrapper
  - Warning filter now catches all sklearn warnings, not just UserWarning
  - Applies to pretrained model, normalizer, and ensemble loading paths
  - Suppresses version mismatch warnings for LinearSVC across sklearn 1.7.2 → 1.8.0

## [2.0.10] - 2024-12-15

### Changed
- Default model switched to isotonic calibration
- Previous sigmoid model saved as `default_pretrained.model.sigmoid.pkl`

### Fixed
- Sklearn version warnings suppressed in model loading

## [2.0.9] - 2024-12-15

### Fixed
- `intronIC test` command now reads correct metrics from `.metrics.iic.json`
- Displays accurate intron counts with thousand separators

## [2.0.8] - 2024-12-15

### Fixed
- UnboundLocalError in test command variable initialization
- Initial attempt at suppressing sklearn warnings

## [2.0.7] - 2024-12-14

### Added
- `intronIC test` subcommand for installation testing
- Bundled Chr19 test data in PyPI package (~20MB)

### Changed
- Updated README with new test command

## [2.0.6] - 2024-12-08

### Changed
- Remove Python version upper bounds
- Remove numpy <2.0 constraint
- Update pixi.toml to match pyproject.toml constraints

## [2.0.0] - 2024-11-26

### Changed
- Complete rewrite with modular architecture
- src/ layout for better packaging isolation
- Streaming mode for memory-efficient processing
- Domain adaptation with adaptive/frozen normalizer modes
- Improved CLI with classify/train/extract subcommands

### Added
- Per-contig streaming classification
- JSON metrics output (.metrics.iic.json)
- Rich console output with progress tracking
- Comprehensive logging system
- Model metadata tracking

### Fixed
- Memory optimization for large genomes
- Coordinate system handling
- Feature extraction pipeline

## [1.5.1] - Earlier

Previous monolithic version. See git history for details.

---

[2.1.3]: https://github.com/glarue/intronIC/compare/v2.1.2...v2.1.3
[2.1.2]: https://github.com/glarue/intronIC/compare/v2.1.0...v2.1.2
[2.1.0]: https://github.com/glarue/intronIC/compare/v2.0.10...v2.1.0
[2.0.10]: https://github.com/glarue/intronIC/compare/v2.0.9...v2.0.10
[2.0.9]: https://github.com/glarue/intronIC/compare/v2.0.8...v2.0.9
[2.0.8]: https://github.com/glarue/intronIC/compare/v2.0.7...v2.0.8
[2.0.7]: https://github.com/glarue/intronIC/compare/v2.0.6...v2.0.7
[2.0.6]: https://github.com/glarue/intronIC/compare/v2.0.0...v2.0.6
[2.0.0]: https://github.com/glarue/intronIC/compare/v1.5.1...v2.0.0
