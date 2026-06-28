"""
Scoring module for intron PWM scoring.

This module provides the scoring pipeline for intronIC, including:
- PWM (Position Weight Matrix) scoring
- Branch point detection
- Score orchestration

(Per-species z-normalization was removed in the supplant 2b refactor; scoring now
operates directly on the background-corrected raw motif log-odds.)
"""

from intronIC.scoring.pwm import PWM, PWMSet, PWMLoader
from intronIC.scoring.branch_point import BranchPointMatch, BranchPointScorer
from intronIC.scoring.scorer import IntronScorer

__all__ = [
    'PWM',
    'PWMSet',
    'PWMLoader',
    'BranchPointMatch',
    'BranchPointScorer',
    'IntronScorer',
]
