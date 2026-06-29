"""
Factory for creating SVM estimators with unified interface.

Handles the LinearSVC vs SVC distinction:
- Linear kernel: uses LinearSVC (primal formulation, supports L1/L2, faster)
- RBF/other kernels: uses SVC (dual formulation, supports kernels)

Both are wrapped in the same Pipeline structure, so downstream code
(CalibratedClassifierCV, GridSearchCV) doesn't need to branch.
"""

from sklearn.svm import LinearSVC, SVC
from sklearn.preprocessing import StandardScaler, RobustScaler


#: Accepted ``scaler`` names for the in-pipeline feature scaler (bundle-stamped on SVMParameters).
SCALER_CHOICES = ("standard", "robust", "none")


def make_scaler_step(scaler: str = "standard"):
    """Return the in-pipeline feature scaler estimator for ``scaler``, or ``None`` to drop the scale step.

    The scaler is **global** — fit once on the training corpus *inside* the estimator (so it never touches
    the interpretable ``score_info.iic`` columns), and frozen at inference. It is NEVER refit per-species;
    per-species re-anchoring is the z-inflation removed in supplant 2b (see docs/raw_gated_scoring.md).

    - ``"standard"`` -> :class:`StandardScaler` (mean/std; the historical default).
    - ``"robust"``   -> :class:`RobustScaler` (median/IQR; robust to the fat negative motif tails).
    - ``"none"``     -> ``None`` (caller drops the ``'scale'`` step; raw features feed the kernel directly).
    """
    s = (scaler or "standard").lower()
    if s in ("standard", "standardscaler"):
        return StandardScaler()
    if s in ("robust", "robustscaler"):
        return RobustScaler()
    if s in ("none", "off", "identity", ""):
        return None
    raise ValueError(f"unknown scaler {scaler!r}; expected one of {SCALER_CHOICES}")


def create_svm(
    kernel: str = "rbf",
    C: float = 1.0,
    gamma="scale",
    class_weight=None,
    penalty: str = "l2",
    loss: str = "squared_hinge",
    max_iter: int = 50000,
    tol: float = 1e-4,
    random_state: int = 42,
):
    """Create an SVM estimator based on kernel type.

    For linear kernel: returns LinearSVC (faster, supports L1 penalty).
    For other kernels: returns SVC (supports RBF, poly, etc.).

    Both have decision_function(). Neither has predict_proba() —
    that's handled by CalibratedClassifierCV wrapping the pipeline.

    Args:
        kernel: 'linear', 'rbf', 'poly', etc.
        C: Regularization parameter
        gamma: Kernel coefficient (ignored for linear)
        class_weight: Class weights dict or 'balanced'
        penalty: 'l1' or 'l2' (only used for linear kernel)
        loss: 'hinge' or 'squared_hinge' (only used for linear kernel)
        max_iter: Maximum iterations
        tol: Convergence tolerance
        random_state: Random seed

    Returns:
        Configured LinearSVC or SVC instance
    """
    if kernel == "linear":
        return LinearSVC(
            C=C,
            penalty=penalty,
            loss=loss,
            dual=False,
            class_weight=class_weight,
            max_iter=max_iter,
            tol=tol,
            random_state=random_state,
        )
    else:
        return SVC(
            kernel=kernel,
            C=C,
            gamma=gamma,
            class_weight=class_weight,
            max_iter=max_iter,
            tol=tol,
            random_state=random_state,
        )


def get_grid_params(kernel: str, svc_key: str = "svc__",
                    penalty_options=None, loss_options=None,
                    gamma_search=None):
    """Get kernel-specific parameter grid entries.

    Args:
        kernel: 'linear' or 'rbf'
        svc_key: Parameter prefix ('svc__' or 'estimator__svc__')
        penalty_options: Penalty types for linear kernel (e.g. ['l1', 'l2'])
        loss_options: Loss functions for linear kernel
        gamma_search: Gamma values for RBF kernel

    Returns:
        Dict of parameter grid entries to merge into the main grid
    """
    params = {}

    if kernel == "linear":
        if penalty_options:
            params[f"{svc_key}penalty"] = penalty_options
        if loss_options:
            # Filter to valid losses for dual=False
            valid = [l for l in loss_options if l == "squared_hinge"]
            if valid:
                params[f"{svc_key}loss"] = valid
    else:
        # RBF or other kernel — add gamma to grid
        if gamma_search:
            resolved = []
            for g in gamma_search:
                if isinstance(g, str) and g in ("scale", "auto"):
                    resolved.append(g)
                else:
                    resolved.append(float(g))
            params[f"{svc_key}gamma"] = resolved

    return params
