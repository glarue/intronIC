"""
Classification module - SVM-based intron classification pipeline.

This module provides components for U2/U12 intron classification using
Support Vector Machines (SVM) with ensemble training.

Components:
    - SVMOptimizer: Hyperparameter optimization via geometric grid search
    - SVMTrainer: Ensemble training with U2 subsampling
    - SVMPredictor: Apply trained models to classify introns
    - IntronClassifier: High-level orchestrator for complete pipeline
"""

from classification.optimizer import (
    SVMOptimizer,
    SVMParameters,
    OptimizationRound,
)
from classification.trainer import (
    SVMTrainer,
    SVMModel,
    SVMEnsemble,
)
from classification.predictor import (
    SVMPredictor,
)
from classification.classifier import (
    IntronClassifier,
    ClassificationResult,
)

__all__ = [
    "SVMOptimizer",
    "SVMParameters",
    "OptimizationRound",
    "SVMTrainer",
    "SVMModel",
    "SVMEnsemble",
    "SVMPredictor",
    "IntronClassifier",
    "ClassificationResult",
]
