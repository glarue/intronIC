#!/usr/bin/env python3
"""
Convert a trained intronIC model's calibration method.

Extracts the fitted LinearSVC from each ensemble model and refits
only the calibration layer (sigmoid or isotonic) on reference data.
Much faster than full retraining since hyperparameters are preserved.

Usage:
    python scripts/convert_calibration.py input.model.pkl sigmoid output.model.pkl
    python scripts/convert_calibration.py input.model.pkl isotonic output.model.pkl
"""

import argparse
import sys
import warnings
from pathlib import Path

import joblib
import numpy as np
from sklearn.calibration import CalibratedClassifierCV

# Suppress FutureWarning about cv='prefit' deprecation
# This warning is about sklearn 1.8+ which will use FrozenEstimator instead
warnings.filterwarnings(
    "ignore",
    message=".*cv='prefit'.*",
    category=FutureWarning,
)


def get_data_dir() -> Path:
    """Get the intronIC data directory."""
    return Path(__file__).parent.parent / "src" / "intronIC" / "data"


def load_and_score_reference_data(normalizer):
    """Load reference data, score it, and normalize using the model's normalizer.

    Args:
        normalizer: The ScoreNormalizer from the model bundle

    Returns:
        Tuple of (X features array, y labels array)
    """
    from intronIC.cli.main import load_reference_sequences
    from intronIC.scoring.pwm import PWMLoader
    from intronIC.scoring.scorer import IntronScorer

    data_dir = get_data_dir()

    # Load reference sequences
    u12_file = data_dir / "u12_reference.introns.iic.gz"
    u2_file = data_dir / "u2_reference.introns.iic.gz"

    if not u12_file.exists() or not u2_file.exists():
        raise FileNotFoundError(
            f"Reference data not found. U12: {u12_file}, U2: {u2_file}"
        )

    u12_introns = load_reference_sequences(u12_file)
    u2_introns = load_reference_sequences(u2_file)
    print(f"Loaded {len(u12_introns)} U12 and {len(u2_introns)} U2 reference sequences")

    # Load PWM matrices
    pwm_file = data_dir / "intronIC_scoring_PWMs.json"
    pwm_sets = PWMLoader.load_from_file(pwm_file, pseudocount=0.0001)

    # Score reference introns
    scorer = IntronScorer(
        pwm_sets=pwm_sets,
        five_coords=(-3, 9),
        bp_coords=(-55, -5),
        three_coords=(-6, 4),
        ignore_nc_dnts=True,
    )

    u12_scored = [scorer.score_intron(i) for i in u12_introns]
    u2_scored = [scorer.score_intron(i) for i in u2_introns]

    # Normalize using the model's normalizer
    u12_normalized = list(normalizer.transform(u12_scored, dataset_type="reference"))
    u2_normalized = list(normalizer.transform(u2_scored, dataset_type="reference"))

    # Extract z-score features
    def extract_features(introns):
        features = []
        for intron in introns:
            if (intron.scores.five_z_score is not None and
                intron.scores.bp_z_score is not None and
                intron.scores.three_z_score is not None):
                features.append([
                    intron.scores.five_z_score,
                    intron.scores.bp_z_score,
                    intron.scores.three_z_score,
                ])
        return np.array(features)

    X_u12 = extract_features(u12_normalized)
    X_u2 = extract_features(u2_normalized)

    print(f"Valid features: {len(X_u12)} U12, {len(X_u2)} U2")

    # Combine and create labels
    X = np.vstack([X_u12, X_u2])
    y = np.array([1] * len(X_u12) + [0] * len(X_u2))

    return X, y


def convert_calibration(input_path: Path, method: str, output_path: Path):
    """Convert model calibration method."""

    if method not in ('sigmoid', 'isotonic'):
        raise ValueError(f"Method must be 'sigmoid' or 'isotonic', got '{method}'")

    print(f"Loading model from {input_path}")
    model_data = joblib.load(input_path)

    ensemble = model_data['ensemble']
    normalizer = model_data['normalizer']
    print(f"Ensemble has {len(ensemble.models)} models")

    # Check current calibration method
    first_model = ensemble.models[0]
    current_method = first_model.parameters.calibration_method if hasattr(first_model, 'parameters') else 'unknown'
    print(f"Current calibration method: {current_method}")

    if current_method == method:
        print(f"Model already uses {method} calibration. Nothing to do.")
        return

    # Load reference data for refitting calibration
    print("Loading and scoring reference data...")
    X, y = load_and_score_reference_data(normalizer)
    print(f"Reference data shape: {X.shape}")

    # Convert each model in ensemble
    print(f"Converting to {method} calibration...")

    # We need to modify the models in the ensemble
    # SVMEnsemble and SVMModel are frozen dataclasses, so we need to recreate them
    from intronIC.classification.trainer import SVMModel, SVMEnsemble
    from intronIC.classification.optimizer import SVMParameters

    new_models = []
    for i, svm_model in enumerate(ensemble.models):
        old_calibrated = svm_model.model

        # Extract the fitted base estimator
        # CalibratedClassifierCV stores calibrated_classifiers_ after fit
        # Each has .estimator which is the fitted base model (a Pipeline)
        if hasattr(old_calibrated, 'calibrated_classifiers_'):
            # Get the base estimator from the first calibrator
            base_estimator = old_calibrated.calibrated_classifiers_[0].estimator
        elif hasattr(old_calibrated, 'estimator'):
            base_estimator = old_calibrated.estimator
        else:
            raise ValueError(f"Cannot extract base estimator from model {i}")

        # Create new calibrated classifier with desired method
        new_calibrated = CalibratedClassifierCV(
            estimator=base_estimator,
            method=method,
            cv='prefit'  # Use already-fitted estimator
        )

        # Fit only the calibration layer
        new_calibrated.fit(X, y)

        # Update parameters metadata
        old_params = svm_model.parameters
        new_params = SVMParameters(
            C=old_params.C,
            calibration_method=method,
            saturate_enabled=old_params.saturate_enabled,
            include_max=old_params.include_max,
            include_pairwise_mins=old_params.include_pairwise_mins,
            penalty=old_params.penalty,
            class_weight_multiplier=old_params.class_weight_multiplier,
            loss=old_params.loss,
            gamma_imbalance=old_params.gamma_imbalance,
            dual=old_params.dual,
            intercept_scaling=old_params.intercept_scaling,
            cv_score=old_params.cv_score,
            round_found=old_params.round_found,
        )

        # Create new SVMModel (frozen dataclass, so we make a new one)
        new_model = SVMModel(
            model=new_calibrated,
            train_size=svm_model.train_size,
            u12_count=svm_model.u12_count,
            u2_count=svm_model.u2_count,
            parameters=new_params,
        )
        new_models.append(new_model)

        print(f"  Converted model {i+1}/{len(ensemble.models)}")

    # Create new ensemble
    new_ensemble = SVMEnsemble(
        models=new_models,
        subsample_ratio=ensemble.subsample_ratio,
    )

    # Update model data
    model_data['ensemble'] = new_ensemble

    # Save converted model
    print(f"Saving converted model to {output_path}")
    joblib.dump(model_data, output_path, compress=3)
    print("Done!")


def main():
    parser = argparse.ArgumentParser(
        description="Convert intronIC model calibration method",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        "input_model",
        type=Path,
        help="Path to input model (.model.pkl)"
    )
    parser.add_argument(
        "method",
        choices=['sigmoid', 'isotonic'],
        help="Calibration method to use"
    )
    parser.add_argument(
        "output_model",
        type=Path,
        help="Path for output model (.model.pkl)"
    )

    args = parser.parse_args()

    if not args.input_model.exists():
        print(f"Error: Input model not found: {args.input_model}", file=sys.stderr)
        sys.exit(1)

    if args.output_model.exists():
        print(f"Warning: Output path exists, will overwrite: {args.output_model}")

    convert_calibration(args.input_model, args.method, args.output_model)


if __name__ == '__main__':
    main()
