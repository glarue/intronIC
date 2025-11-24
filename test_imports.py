#!/usr/bin/env python3
"""Quick test of imports after src layout migration."""

print("Testing imports...")

from intronIC.cli.main import main
print("✓ intronIC.cli.main")

from intronIC.scoring.pwm import PWMLoader
print("✓ intronIC.scoring.pwm")

from intronIC.core.intron import Intron
print("✓ intronIC.core.intron")

from intronIC.extraction.annotator import AnnotationHierarchyBuilder
print("✓ intronIC.extraction.annotator")

from intronIC.file_io.parsers import BioGLAnnotationParser
print("✓ intronIC.file_io.parsers")

from intronIC.classification.classifier import IntronClassifier
print("✓ intronIC.classification.classifier")

from intronIC.visualization.plots import plot_training_results
print("✓ intronIC.visualization.plots")

print("\n✅ All imports successful!")
