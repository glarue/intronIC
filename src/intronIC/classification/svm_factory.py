"""
Factory for creating SVM estimators with unified interface.

Handles the LinearSVC vs SVC distinction:
- Linear kernel: uses LinearSVC (primal formulation, supports L1/L2, faster)
- RBF/other kernels: uses SVC (dual formulation, supports kernels)

Both are wrapped in the same Pipeline structure, so downstream code
(CalibratedClassifierCV, GridSearchCV) doesn't need to branch.
"""

from sklearn.svm import LinearSVC, SVC


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
