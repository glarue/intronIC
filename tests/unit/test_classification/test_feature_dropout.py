"""
Smoke tests for feature-dropout ensemble training.

Validates the feature dropout mechanism on a small synthetic dataset
before integrating into the full training pipeline. Tests:

1. Feature masking produces different predictions than unmasked
2. Models trained with different dropped features disagree on
   feature-dependent introns but agree on robust ones
3. The ensemble mean is still accurate despite dropout
4. Ensemble σ is higher for feature-dependent introns
"""

import numpy as np
import pytest
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.calibration import CalibratedClassifierCV


def _make_synthetic_data(n_u2=500, n_u12=50, n_features=6, seed=42):
    """Create synthetic intron-like data with known properties.

    U12 introns: high on all features (multivariate normal, mean=3)
    U2 introns: centered at 0 (multivariate normal, mean=0)
    FP introns: high on feature 0 only (mean=3), low elsewhere (mean=0)

    Returns:
        X_train, y_train: training data (U12 + U2)
        X_fp: false positive introns (strong on one feature only)
        X_real_u12: additional real U12 introns for testing
    """
    rng = np.random.RandomState(seed)

    # U2 introns: centered at origin
    X_u2 = rng.randn(n_u2, n_features) * 0.5

    # U12 introns: all features elevated
    X_u12 = rng.randn(n_u12, n_features) * 0.5 + 3.0

    X_train = np.vstack([X_u2, X_u12])
    y_train = np.array([0] * n_u2 + [1] * n_u12)

    # FP introns: feature 0 is high (like U12), rest are U2-like
    n_fp = 20
    X_fp = rng.randn(n_fp, n_features) * 0.5
    X_fp[:, 0] = 3.0 + rng.randn(n_fp) * 0.5  # strong feature 0

    # Additional real U12 for testing
    X_real = rng.randn(20, n_features) * 0.5 + 3.0

    return X_train, y_train, X_fp, X_real


def _train_model_with_dropout(X_train, y_train, drop_feature=None, seed=42):
    """Train a single SVM with one feature dropped (masked to median)."""
    X = X_train.copy()
    if drop_feature is not None:
        median_val = np.median(X[:, drop_feature])
        X[:, drop_feature] = median_val

    pipe = Pipeline([
        ('scale', StandardScaler()),
        ('svc', SVC(kernel='rbf', C=10, gamma='scale', probability=True,
                    random_state=seed)),
    ])
    pipe.fit(X, y_train)
    return pipe, drop_feature, np.median(X_train[:, drop_feature]) if drop_feature is not None else None


def _predict_with_mask(model, X, drop_feature, median_val):
    """Predict with the same feature mask used during training."""
    X_masked = X.copy()
    if drop_feature is not None:
        X_masked[:, drop_feature] = median_val
    return model.predict_proba(X_masked)[:, 1]


class TestFeatureDropoutMechanism:
    """Basic mechanism tests."""

    def test_masking_changes_predictions(self):
        """Masking a feature produces different predictions than unmasked."""
        X_train, y_train, X_fp, _ = _make_synthetic_data()

        model_full, _, _ = _train_model_with_dropout(X_train, y_train, drop_feature=None)
        model_drop0, _, med0 = _train_model_with_dropout(X_train, y_train, drop_feature=0)

        probs_full = model_full.predict_proba(X_fp)[:, 1]
        probs_drop = _predict_with_mask(model_drop0, X_fp, 0, med0)

        # FPs have strong feature 0 — dropping it should change their scores
        assert not np.allclose(probs_full, probs_drop, atol=0.1), \
            "Dropping feature 0 should change FP predictions"

    def test_real_u12_robust_to_dropout(self):
        """Real U12 introns (strong on all features) are robust to any single dropout."""
        X_train, y_train, _, X_real = _make_synthetic_data()

        model_full, _, _ = _train_model_with_dropout(X_train, y_train, drop_feature=None)
        probs_full = model_full.predict_proba(X_real)[:, 1]

        for feat in range(6):
            model_drop, _, med = _train_model_with_dropout(X_train, y_train, drop_feature=feat)
            probs_drop = _predict_with_mask(model_drop, X_real, feat, med)

            # Real U12s should still be classified as U12 (>0.5) even with dropout
            assert np.mean(probs_drop > 0.5) >= 0.8, \
                f"Dropping feature {feat} should not destroy real U12 classification"

    def test_fps_sensitive_to_their_strong_feature(self):
        """FPs with one strong feature are sensitive to dropping that feature."""
        X_train, y_train, X_fp, _ = _make_synthetic_data()

        # FPs are strong on feature 0
        model_drop0, _, med0 = _train_model_with_dropout(X_train, y_train, drop_feature=0)
        probs_drop0 = _predict_with_mask(model_drop0, X_fp, 0, med0)

        # Dropping feature 1 (which FPs are NOT strong on) should matter less
        model_drop1, _, med1 = _train_model_with_dropout(X_train, y_train, drop_feature=1)
        probs_drop1 = _predict_with_mask(model_drop1, X_fp, 1, med1)

        # FP scores should drop more when feature 0 is dropped than feature 1
        assert np.median(probs_drop0) < np.median(probs_drop1), \
            "FPs should be more affected by dropping their strong feature"


class TestDropoutEnsemble:
    """Test the full dropout ensemble behavior."""

    def test_ensemble_variance_separates_fps_from_real(self):
        """Dropout ensemble produces higher σ for FPs than real U12s."""
        X_train, y_train, X_fp, X_real = _make_synthetic_data()

        # Train 6 models, each dropping a different feature
        models = []
        for feat in range(6):
            model, drop_feat, med = _train_model_with_dropout(
                X_train, y_train, drop_feature=feat, seed=42 + feat
            )
            models.append((model, drop_feat, med))

        # Get ensemble predictions for FPs and real U12s
        def ensemble_stats(X):
            all_probs = np.zeros((len(X), len(models)))
            for i, (model, drop_feat, med) in enumerate(models):
                all_probs[:, i] = _predict_with_mask(model, X, drop_feat, med)
            means = np.mean(all_probs, axis=1)
            stds = np.std(all_probs, axis=1)
            return means, stds

        fp_means, fp_stds = ensemble_stats(X_fp)
        real_means, real_stds = ensemble_stats(X_real)

        # FPs should have higher ensemble variance than real U12s
        assert np.median(fp_stds) > np.median(real_stds), \
            "FP ensemble σ ({:.3f}) should exceed real U12 σ ({:.3f})".format(
                np.median(fp_stds), np.median(real_stds))

    def test_ensemble_mean_still_accurate(self):
        """Dropout ensemble mean is still accurate for real U12s and U2s."""
        X_train, y_train, _, X_real = _make_synthetic_data()

        models = []
        for feat in range(6):
            model, drop_feat, med = _train_model_with_dropout(
                X_train, y_train, drop_feature=feat, seed=42 + feat
            )
            models.append((model, drop_feat, med))

        # Real U12s should still have high mean probability
        all_probs = np.zeros((len(X_real), len(models)))
        for i, (model, drop_feat, med) in enumerate(models):
            all_probs[:, i] = _predict_with_mask(model, X_real, drop_feat, med)
        means = np.mean(all_probs, axis=1)

        assert np.median(means) > 0.7, \
            "Dropout ensemble should still correctly classify real U12s"

    def test_sigma_penalty_rejects_fps_keeps_real(self):
        """The σ penalty formula rejects FPs while keeping real U12s."""
        X_train, y_train, X_fp, X_real = _make_synthetic_data()

        models = []
        for feat in range(6):
            model, drop_feat, med = _train_model_with_dropout(
                X_train, y_train, drop_feature=feat, seed=42 + feat
            )
            models.append((model, drop_feat, med))

        def get_adjusted_scores(X, sigma_ref=0.05, k=5):
            all_probs = np.zeros((len(X), len(models)))
            for i, (model, drop_feat, med) in enumerate(models):
                all_probs[:, i] = _predict_with_mask(model, X, drop_feat, med)
            means = np.mean(all_probs, axis=1)
            stds = np.std(all_probs, axis=1)
            excess = np.maximum(0, stds - sigma_ref)
            adjusted = means - k * excess
            return adjusted

        fp_adjusted = get_adjusted_scores(X_fp)
        real_adjusted = get_adjusted_scores(X_real)

        # Real U12s should mostly survive (adjusted > 0.5)
        real_survive = np.mean(real_adjusted > 0.5)
        # FPs should mostly be rejected (adjusted < 0.5)
        fp_rejected = np.mean(fp_adjusted < 0.5)

        assert real_survive > 0.7, \
            "Most real U12s should survive σ penalty (got {:.0%})".format(real_survive)
        assert fp_rejected > 0.5, \
            "Most FPs should be rejected by σ penalty (got {:.0%})".format(fp_rejected)
