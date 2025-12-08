"""
Scoring module for intron PWM scoring and normalization.

This module provides the scoring pipeline for intronIC, including:
- PWM (Position Weight Matrix) scoring
- Branch point detection
- Z-score normalization (with ML integrity guarantees)
- Score orchestration
"""

from intronIC.scoring.normalizer import ScoreNormalizer, DatasetType
from intronIC.scoring.pwm import PWM, PWMSet, PWMLoader
from intronIC.scoring.branch_point import BranchPointMatch, BranchPointScorer
from intronIC.scoring.scorer import IntronScorer

__all__ = [
    'ScoreNormalizer',
    'DatasetType',
    'PWM',
    'PWMSet',
    'PWMLoader',
    'BranchPointMatch',
    'BranchPointScorer',
    'IntronScorer',
]
