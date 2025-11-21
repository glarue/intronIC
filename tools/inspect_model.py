#!/usr/bin/env python3
"""
Inspect a trained model to see its features and coefficients.
"""
import joblib
import sys
from pathlib import Path

def inspect_model(model_path):
    """Load and inspect a pickled model."""
    print(f"Loading model from {model_path}...")

    ensemble = joblib.load(model_path)

    print(f"\nModel type: {type(ensemble)}")
    if isinstance(ensemble, dict):
        print(f"Keys: {list(ensemble.keys())}")
        print(f"Number of items: {len(ensemble)}")
    else:
        print(f"Number of models in ensemble: {len(ensemble)}")

    # Get models from dict or list
    if isinstance(ensemble, dict):
        if 'models' in ensemble:
            models_obj = ensemble['models']
        elif 'ensemble' in ensemble:
            models_obj = ensemble['ensemble']
        else:
            print("Unknown dict structure!")
            return
    else:
        models_obj = ensemble

    print(f"\nEnsemble object type: {type(models_obj)}")

    # If it's an SVMEnsemble object, get the models attribute
    if hasattr(models_obj, 'models'):
        models = models_obj.models
        print(f"Number of models: {len(models)}")
    elif hasattr(models_obj, '__len__'):
        models = models_obj
        print(f"Number of models: {len(models)}")
    else:
        print(f"Object attributes: {dir(models_obj)}")
        return

    # Inspect first model in ensemble
    print("\n" + "="*80)
    print("FIRST MODEL IN ENSEMBLE:")
    print("="*80)

    model = models[0]
    print(f"\nModel structure: {type(model)}")
    print(f"Model attributes: {[a for a in dir(model) if not a.startswith('_')]}")

    # Check if it's an SVMModel object
    if hasattr(model, 'model'):
        print("\nFound model attribute")
        pipeline = model.model
        print(f"Pipeline type: {type(pipeline)}")
        if hasattr(model, 'parameters'):
            try:
                print(f"Parameters: {model.parameters}")
            except AttributeError as e:
                # Handle old pickled models with missing fields
                print(f"Parameters (partial - old model format):")
                params = model.parameters
                for attr in ['C', 'penalty', 'loss', 'calibration_method', 'class_weight_multiplier']:
                    if hasattr(params, attr):
                        print(f"  {attr}: {getattr(params, attr)}")
                # Try new fields
                if hasattr(params, 'gamma_imbalance'):
                    print(f"  gamma_imbalance: {params.gamma_imbalance}")
                print(f"  (Note: Some fields may be missing in old model format)")

        # Check if it's a CalibratedClassifierCV
        if hasattr(pipeline, 'calibrated_classifiers_'):
            print("\nCalibratedClassifierCV found")
            base_estimator = pipeline.calibrated_classifiers_[0].estimator
            print(f"Base estimator: {type(base_estimator)}")

            # Check if base estimator is a Pipeline
            if hasattr(base_estimator, 'named_steps'):
                print("\nBase estimator is a Pipeline")
                print("Pipeline steps:")
                for name, step in base_estimator.named_steps.items():
                    print(f"  {name}: {type(step)}")

                # Get the transformer (check both 'transform' and 'augment' names)
                transformer = None
                if 'transform' in base_estimator.named_steps:
                    transformer = base_estimator.named_steps['transform']
                elif 'augment' in base_estimator.named_steps:
                    transformer = base_estimator.named_steps['augment']

                if transformer:
                    print(f"\nTransformer: {type(transformer)}")
                    # Check for explicit features list (new models)
                    if hasattr(transformer, 'features') and transformer.features is not None:
                        print(f"  features: {transformer.features}")
                    # Check for legacy flags (old models)
                    elif hasattr(transformer, 'include_max'):
                        print(f"  include_max: {transformer.include_max}")
                        if hasattr(transformer, 'include_pairwise_mins'):
                            print(f"  include_pairwise_mins: {transformer.include_pairwise_mins}")

                    if hasattr(transformer, 'get_feature_names_out'):
                        feature_names = transformer.get_feature_names_out()
                        print(f"  Feature names ({len(feature_names)}): {list(feature_names)}")
                    else:
                        feature_names = None
                else:
                    feature_names = None

                # Get the final classifier
                if 'svc' in base_estimator.named_steps:
                    final_svc = base_estimator.named_steps['svc']
                    print(f"\nFinal classifier: {type(final_svc)}")

                    if hasattr(final_svc, 'coef_'):
                        coefs = final_svc.coef_[0]
                        print(f"\n  Coefficients ({len(coefs)} features):")

                        if feature_names is None:
                            feature_names = [f"feature_{i}" for i in range(len(coefs))]

                        for name, coef in zip(feature_names, coefs):
                            sign = '+' if coef >= 0 else ''
                            print(f"    {name:20s}: {sign}{coef:.6f}")

                    if hasattr(final_svc, 'intercept_'):
                        print(f"\n  Intercept: {final_svc.intercept_[0]:.6f}")

    # OLD code for sklearn Pipeline
    elif hasattr(model, 'named_steps'):
        print("\nPipeline steps:")
        for name, step in model.named_steps.items():
            print(f"  {name}: {type(step)}")

        # Get the transformer
        if 'augment' in model.named_steps:
            transformer = model.named_steps['augment']
            print(f"\nTransformer: {type(transformer)}")
            if hasattr(transformer, 'include_max'):
                print(f"  include_max: {transformer.include_max}")
            if hasattr(transformer, 'get_feature_names_out'):
                feature_names = transformer.get_feature_names_out()
                print(f"  Features ({len(feature_names)}): {list(feature_names)}")

        # Get the classifier
        if 'svc' in model.named_steps:
            svc = model.named_steps['svc']
            print(f"\nClassifier: {type(svc)}")

            # Check if it's a CalibratedClassifierCV
            if hasattr(svc, 'calibrated_classifiers_'):
                print("  Calibrated classifier found")
                base_estimator = svc.calibrated_classifiers_[0].estimator
                print(f"  Base estimator: {type(base_estimator)}")

                if hasattr(base_estimator, 'coef_'):
                    coefs = base_estimator.coef_[0]
                    print(f"\n  Coefficients ({len(coefs)} features):")

                    # Get feature names from transformer
                    if 'augment' in model.named_steps:
                        feature_names = model.named_steps['augment'].get_feature_names_out()
                    else:
                        feature_names = [f"feature_{i}" for i in range(len(coefs))]

                    for name, coef in zip(feature_names, coefs):
                        sign = '+' if coef >= 0 else ''
                        print(f"    {name:20s}: {sign}{coef:.6f}")

                if hasattr(base_estimator, 'intercept_'):
                    print(f"\n  Intercept: {base_estimator.intercept_[0]:.6f}")

    # Check metadata
    print("\n" + "="*80)
    print("ENSEMBLE METADATA:")
    print("="*80)

    # Try to get metadata if available
    if hasattr(ensemble, 'metadata'):
        print(f"\nMetadata: {ensemble.metadata}")

    # Check all models in ensemble
    print(f"\nAll {len(models)} models:")
    print("="*80)
    for i, mdl in enumerate(models, 1):
        print(f"\nModel {i}:")
        if hasattr(mdl, 'parameters'):
            try:
                print(f"  Parameters: {mdl.parameters}")
            except AttributeError:
                # Handle old pickled models with missing fields
                params = mdl.parameters
                if hasattr(params, 'C'):
                    print(f"  C: {params.C}")
                if hasattr(params, 'penalty'):
                    print(f"  penalty: {params.penalty}")
                if hasattr(params, 'calibration_method'):
                    print(f"  calibration_method: {params.calibration_method}")

        if hasattr(mdl, 'model'):
            cal_clf = mdl.model
            if hasattr(cal_clf, 'calibrated_classifiers_'):
                base_est = cal_clf.calibrated_classifiers_[0].estimator
                if hasattr(base_est, 'named_steps'):
                    # Check for both 'transform' (new) and 'augment' (old) step names
                    transformer = None
                    if 'transform' in base_est.named_steps:
                        transformer = base_est.named_steps['transform']
                    elif 'augment' in base_est.named_steps:
                        transformer = base_est.named_steps['augment']

                    if transformer and 'svc' in base_est.named_steps:
                        svc = base_est.named_steps['svc']

                        if hasattr(svc, 'coef_'):
                            coefs = svc.coef_[0]
                            if hasattr(transformer, 'get_feature_names_out'):
                                feature_names = transformer.get_feature_names_out()
                            else:
                                feature_names = [f"feature_{j}" for j in range(len(coefs))]

                            print(f"  Coefficients:")
                            for name, coef in zip(feature_names, coefs):
                                if abs(coef) > 0.0001:  # Only non-zero
                                    sign = '+' if coef >= 0 else ''
                                    print(f"    {name:20s}: {sign}{coef:.6f}")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python inspect_model.py <model_path.pkl>")
        print("\nExample:")
        print("  python inspect_model.py homo_sapiens.model.pkl")
        sys.exit(1)

    model_path = Path(sys.argv[1])

    if not model_path.exists():
        print(f"Error: {model_path} not found")
        sys.exit(1)

    inspect_model(model_path)
