"""
Conservative margin-space domain adaptation for cross-species classification.

This module implements the expert-recommended approach for cross-species adaptation:
instead of re-normalizing features (which changes the coordinate system), we align
target species margins to the human negative margin distribution.

Key principles:
1. Preserve human coordinate system (keep human scaler fixed)
2. Align in margin space (after feature transformation)
3. Apply guardrails (shrinkage, clamping) to prevent extreme adaptations
4. Provide rich diagnostics for troubleshooting

Algorithm:
    Given target species margins f_target and human U2 margin statistics (μ_h, σ_h):
    1. Compute target statistics: μ_t = median(f_target), σ_t = IQR(f_target)
    2. Raw alignment: a = σ_h / σ_t, c = μ_h - a * μ_t
    3. Apply shrinkage toward identity: λ = k / (k + N)
       a' = (1-λ)*a + λ*1, c' = (1-λ)*c + λ*0
    4. Clamp scale to reasonable range [0.33, 3.0]
    5. Return: f_aligned = a' * f_target + c'

This is safer than feature re-normalization because:
- The bulk of introns in any species are U2 (99.5%+)
- Aligning bulk margins to human U2 margins is statistically sound
- No change to feature coordinate system → SVM decision boundary stays calibrated
- Guardrails prevent pathological adaptations

Reference:
    Expert guidance on robust cross-species domain adaptation (2025)
    "If you add an 'adaptive' mode, make it: margin alignment + (optionally)
    prior-shift, with clear logging and guardrails."
"""

import numpy as np
from typing import Tuple, Optional
from dataclasses import dataclass


@dataclass(frozen=True)
class MarginAlignmentStats:
    """
    Diagnostics and statistics for margin alignment.

    Attributes:
        mu_human: Median of human U2 margins (from training)
        sigma_human: IQR of human U2 margins (from training)
        n_human: Number of human U2 samples used
        mu_target: Median of target species margins
        sigma_target: IQR of target species margins
        n_target: Number of target species samples
        scale: Final scale factor (a) after shrinkage/clamping
        shift: Final shift factor (c) after shrinkage/clamping
        scale_raw: Raw scale before guardrails
        shift_raw: Raw shift before guardrails
        shrinkage_lambda: Shrinkage parameter λ = k/(k+N)
        scale_clamped: Whether scale was clamped to valid range
    """
    # Human statistics
    mu_human: float
    sigma_human: float
    n_human: int

    # Target statistics
    mu_target: float
    sigma_target: float
    n_target: int

    # Final alignment parameters
    scale: float
    shift: float

    # Diagnostic info
    scale_raw: float
    shift_raw: float
    shrinkage_lambda: float
    scale_clamped: bool

    def format_summary(self) -> str:
        """
        Generate human-readable summary for logging.

        Returns:
            Multi-line string with alignment diagnostics
        """
        lines = [
            "Margin Alignment Diagnostics:",
            f"  Human U2:       μ={self.mu_human:>7.3f}, σ={self.sigma_human:>6.3f} (N={self.n_human:,})",
            f"  Target species: μ={self.mu_target:>7.3f}, σ={self.sigma_target:>6.3f} (N={self.n_target:,})",
            f"  Raw alignment:  scale={self.scale_raw:>6.3f}, shift={self.shift_raw:>7.3f}",
        ]

        # Add shrinkage info if applied
        if self.shrinkage_lambda > 0.001:  # Only show if meaningful
            lines.append(f"  Shrinkage:      λ={self.shrinkage_lambda:.3f} (toward identity)")

        # Warn if scale was clamped
        if self.scale_clamped:
            lines.append(f"  ⚠️  Scale clamped to valid range [0.33, 3.0]")

        # Final transformation
        lines.append(f"  Final:          f' = {self.scale:.3f} × f + {self.shift:.3f}")

        return "\n".join(lines)


def compute_margin_alignment(
    target_margins: np.ndarray,
    human_median: float,
    human_iqr: float,
    human_n: int,
    shrinkage_k: float = 5000.0,
    min_scale: float = 0.33,
    max_scale: float = 3.0,
    epsilon: float = 1e-6
) -> Tuple[float, float, MarginAlignmentStats]:
    """
    Compute conservative margin alignment parameters.

    Aligns target species margin distribution to human U2 margin distribution
    using affine transformation with guardrails.

    Args:
        target_margins: SVM decision function values for target species (N,)
        human_median: Median of human U2 training margins
        human_iqr: IQR (Q75 - Q25) of human U2 training margins
        human_n: Number of human U2 samples (for context)
        shrinkage_k: Shrinkage hyperparameter (default: 5000)
                    Controls λ = k/(k+N). Larger k → more shrinkage.
                    At N=5000: λ=0.5 (moderate shrinkage)
                    At N=50000: λ=0.09 (minimal shrinkage)
        min_scale: Minimum allowed scale factor (default: 0.33)
        max_scale: Maximum allowed scale factor (default: 3.0)
        epsilon: Small value to prevent division by zero

    Returns:
        Tuple of (scale, shift, stats)
        - scale (float): Multiply target margins by this
        - shift (float): Add this after scaling
        - stats (MarginAlignmentStats): Full diagnostics

    Algorithm:
        1. Compute target distribution statistics (median, IQR)
        2. Compute raw alignment: a = σ_h / σ_t, c = μ_h - a*μ_t
        3. Apply shrinkage toward identity transformation
        4. Clamp scale to [min_scale, max_scale]
        5. If clamped, recompute shift to maintain median alignment

    Example:
        >>> # Human U2s have median=-2.0, IQR=1.5
        >>> # Target has median=-1.5, IQR=2.0 (wider spread)
        >>> target = np.random.randn(10000) * 2.0 - 1.5
        >>> a, c, stats = compute_margin_alignment(
        ...     target, human_median=-2.0, human_iqr=1.5, human_n=50000
        ... )
        >>> # Result: a ≈ 0.75 (compress wider target distribution)
        >>> #         c ≈ -0.875 (shift median to match human)
        >>> aligned = a * target + c
        >>> # aligned now has similar distribution to human U2s
    """
    # Validate inputs
    if len(target_margins) == 0:
        raise ValueError("target_margins cannot be empty")
    if human_iqr <= 0:
        raise ValueError(f"human_iqr must be positive, got {human_iqr}")

    # Compute target distribution statistics
    # Use same robust estimators as human (median, IQR)
    mu_target = float(np.median(target_margins))
    q25, q75 = np.percentile(target_margins, [25, 75])
    sigma_target = float(q75 - q25)
    n_target = len(target_margins)

    # Compute raw alignment parameters
    # Goal: Transform target to match human distribution
    # a * (mu_target) + c = mu_human  →  c = mu_human - a * mu_target
    # a * sigma_target = sigma_human  →  a = sigma_human / sigma_target
    a_raw = human_iqr / max(sigma_target, epsilon)
    c_raw = human_median - a_raw * mu_target

    # Guardrail 1: Shrinkage toward identity transformation
    # Rationale: Don't fully trust alignment when:
    # - Sample size is small (high variance in estimates)
    # - Distribution shapes differ significantly
    # Shrink toward (a=1, c=0) which is identity (no adaptation)
    lambda_shrink = shrinkage_k / (shrinkage_k + n_target)

    a_shrunk = (1.0 - lambda_shrink) * a_raw + lambda_shrink * 1.0
    c_shrunk = (1.0 - lambda_shrink) * c_raw + lambda_shrink * 0.0

    # Guardrail 2: Clamp scale to reasonable range
    # Rationale: Extreme scales indicate distribution mismatch
    # - scale < 0.33: Target is 3× wider than human (suspicious)
    # - scale > 3.0: Target is 3× narrower than human (suspicious)
    scale_clamped = False
    if a_shrunk < min_scale or a_shrunk > max_scale:
        a_final = float(np.clip(a_shrunk, min_scale, max_scale))
        scale_clamped = True

        # If we clamped scale, recompute shift to maintain median alignment
        # This preserves central tendency even if spread doesn't match
        c_final = human_median - a_final * mu_target
    else:
        a_final = a_shrunk
        c_final = c_shrunk

    # Build diagnostics object
    stats = MarginAlignmentStats(
        mu_human=human_median,
        sigma_human=human_iqr,
        n_human=human_n,
        mu_target=mu_target,
        sigma_target=sigma_target,
        n_target=n_target,
        scale=a_final,
        shift=c_final,
        scale_raw=a_raw,
        shift_raw=c_raw,
        shrinkage_lambda=lambda_shrink,
        scale_clamped=scale_clamped
    )

    return a_final, c_final, stats


def align_margins(
    margins: np.ndarray,
    scale: float,
    shift: float
) -> np.ndarray:
    """
    Apply affine alignment transformation to margins.

    Args:
        margins: Raw SVM decision function values (N,)
        scale: Scale factor (a)
        shift: Shift factor (c)

    Returns:
        Aligned margins: a * margins + c

    Example:
        >>> margins = np.array([-2.5, -1.0, 0.5, 2.0])
        >>> aligned = align_margins(margins, scale=0.8, shift=-0.5)
        >>> # Result: [-2.5, -1.3, -0.1, 1.1]
    """
    return scale * margins + shift


def fold_alignment_into_platt(
    platt_A: float,
    platt_B: float,
    scale: float,
    shift: float
) -> Tuple[float, float]:
    """
    Fold margin alignment into Platt calibration parameters.

    Instead of computing:
        p = sigmoid(A * (a*f + c) + B)

    We can compute equivalent:
        p = sigmoid(A' * f + B')
        where A' = A * a, B' = A * c + B

    This is more efficient (one fewer operation per prediction) and
    numerically equivalent.

    Args:
        platt_A: Original Platt scale parameter
        platt_B: Original Platt offset parameter
        scale: Margin alignment scale (a)
        shift: Margin alignment shift (c)

    Returns:
        Tuple of (A_new, B_new) for aligned predictions

    Derivation:
        p = sigmoid(A * f_aligned + B)
          = sigmoid(A * (a*f + c) + B)
          = sigmoid(A*a*f + A*c + B)
          = sigmoid(A'*f + B')  where A'=A*a, B'=A*c+B

    Example:
        >>> # Original Platt: A=2.0, B=-1.0
        >>> # Alignment: scale=0.8, shift=-0.5
        >>> A_new, B_new = fold_alignment_into_platt(2.0, -1.0, 0.8, -0.5)
        >>> # Result: A_new=1.6, B_new=-2.0
        >>> # Now: sigmoid(1.6*f - 2.0) ≡ sigmoid(2.0*(0.8*f - 0.5) - 1.0)
    """
    A_new = platt_A * scale
    B_new = platt_A * shift + platt_B

    return A_new, B_new
