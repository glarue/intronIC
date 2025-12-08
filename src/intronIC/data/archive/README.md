# Archived PWM Scoring Matrix Files

This directory contains legacy and supplementary PWM (Position Weight Matrix) scoring files that have been superseded by the unified JSON format (`intronIC_scoring_PWMs.json`).

## Current Production File

The active PWM file used by intronIC is:
- **`../intronIC_scoring_PWMs.json`** - Unified JSON format containing all scoring matrices

## Archived Files

### `scoring_matrices.fasta.iic`
The original legacy PWM file in FASTA-like `.iic` format. Contains:
- U12-type AT-AC five' and three' splice site matrices
- U12-type GT-AG five' and three' splice site matrices  
- U12-type branch point matrices (A9 and A10 variants)
- U2-type GT-AG and GC-AG five' and three' splice site matrices

**Note:** This file does NOT contain U2 branch point matrices. Those were added separately (see below).

### `u2.conserved_empirical_bp_pwm.iic`
U2-type branch point PWM derived from human U2-type introns conserved across mouse, zebrafish, and horseshoe crab. Branchpoint sequences were extracted from supplemental data in:

> Mercer et al. 2015 (https://pubmed.ncbi.nlm.nih.gov/25561518/)

Sequences with an adenosine at position +9 or +10 (of 12) were used to construct the PWM (n=2044 sequences).

This matrix is now included in the production JSON file as `u2-gtag-bp`.

### `scoring_matrices.2+.pub_required.fasta.iic`
Alternative scoring matrices requiring additional publication citations.

### `scoring_matrices.2+_pubs.fasta.iic`
Scoring matrices with expanded publication-derived training data.

### `empirical.u12s.2+.pub_required.matrices`
Empirical U12 matrix data requiring additional publication citations.

## Format Migration

The JSON format (`intronIC_scoring_PWMs.json`) was created to:
1. Consolidate all matrices into a single file
2. Include the U2 branch point matrix from Mercer et al. 2015
3. Provide machine-readable metadata (sample sizes, start indices, bases)
4. Enable easier validation and tooling

## Usage

The archived `.iic` files should not be used directly. They are retained for:
- Historical reference
- Format equivalence testing
- Regeneration of the JSON file if needed
