"""
Calibration quality validation tools for Phase 2A.

This module provides metrics and visualization tools to compare calibration
quality between sigmoid and isotonic calibration methods.

Metrics:
- Log-loss: Measures probability calibration quality
- Expected Calibration Error (ECE): Measures calibration accuracy
- Brier score: Combines calibration + discrimination

Usage:
    from intronIC.classification.calibration_validator import CalibrationValidator

    validator = CalibrationValidator()
    metrics = validator.compare_calibration_methods(
        X_val, y_val, model_sigmoid, model_isotonic
    )
"""

from dataclasses import dataclass
from typing import Tuple, Dict, Optional
import numpy as np
from sklearn.metrics import log_loss, brier_score_loss
from sklearn.calibration import calibration_curve
import matplotlib.pyplot as plt


@dataclass(frozen=True, slots=True)
class CalibrationMetrics:
    """Calibration quality metrics for a single model."""

    log_loss: float  # Lower = better calibrated probabilities
    brier_score: float  # Lower = better (combines calibration + discrimination)
    ece: float  # Expected Calibration Error (lower = better)
    mce: float  # Maximum Calibration Error (lower = better)

    def __str__(self) -> str:
        return (
            f"CalibrationMetrics(\n"
            f"  log_loss={self.log_loss:.6f},\n"
            f"  brier_score={self.brier_score:.6f},\n"
            f"  ECE={self.ece:.6f},\n"
            f"  MCE={self.mce:.6f}\n"
            f")"
        )


class CalibrationValidator:
    """
    Validate and compare calibration quality of different models.

    Provides comprehensive calibration metrics and visualizations for
    comparing sigmoid vs isotonic calibration methods.
    """

    def __init__(self, n_bins: int = 10):
        """
        Initialize validator.

        Args:
            n_bins: Number of bins for ECE calculation and reliability diagrams
        """
        self.n_bins = n_bins

    def compute_ece(
        self,
        y_true: np.ndarray,
        y_prob: np.ndarray,
        n_bins: Optional[int] = None
    ) -> Tuple[float, float]:
        """
        Compute Expected Calibration Error (ECE) and Maximum Calibration Error (MCE).

        ECE measures average calibration error across bins.
        MCE measures worst-case calibration error.

        Args:
            y_true: True binary labels (0/1)
            y_prob: Predicted probabilities for class 1
            n_bins: Number of bins (defaults to self.n_bins)

        Returns:
            (ece, mce): Expected and maximum calibration errors

        Reference:
            Naeini et al. (2015): "Obtaining Well Calibrated Probabilities Using
            Bayesian Binning" (AAAI)
        """
        if n_bins is None:
            n_bins = self.n_bins

        # Get calibration curve
        prob_true, prob_pred = calibration_curve(
            y_true, y_prob, n_bins=n_bins, strategy='uniform'
        )

        # Bin probabilities
        bins = np.linspace(0, 1, n_bins + 1)
        bin_indices = np.digitize(y_prob, bins) - 1
        bin_indices = np.clip(bin_indices, 0, n_bins - 1)

        # Compute calibration error per bin
        bin_errors = []
        bin_weights = []

        for i in range(n_bins):
            mask = (bin_indices == i)
            if mask.sum() > 0:
                # Average predicted probability in this bin
                avg_pred = y_prob[mask].mean()
                # Actual fraction of positives in this bin
                avg_true = y_true[mask].mean()
                # Calibration error for this bin
                error = abs(avg_pred - avg_true)
                weight = mask.sum() / len(y_prob)

                bin_errors.append(error)
                bin_weights.append(weight)

        # ECE: weighted average of bin errors
        if bin_errors:
            ece = np.average(bin_errors, weights=bin_weights)
            mce = max(bin_errors)  # MCE: maximum bin error
        else:
            ece = mce = 0.0

        return float(ece), float(mce)

    def evaluate_model(
        self,
        y_true: np.ndarray,
        y_prob: np.ndarray
    ) -> CalibrationMetrics:
        """
        Compute all calibration metrics for a single model.

        Args:
            y_true: True binary labels (0/1)
            y_prob: Predicted probabilities for class 1

        Returns:
            CalibrationMetrics with all metrics
        """
        # Ensure probabilities are for class 1
        if len(y_prob.shape) > 1 and y_prob.shape[1] == 2:
            y_prob = y_prob[:, 1]  # Extract class 1 probabilities

        # Log-loss (cross-entropy)
        loss = log_loss(y_true, y_prob)

        # Brier score
        brier = brier_score_loss(y_true, y_prob)

        # ECE and MCE
        ece, mce = self.compute_ece(y_true, y_prob)

        return CalibrationMetrics(
            log_loss=loss,
            brier_score=brier,
            ece=ece,
            mce=mce
        )

    def compare_methods(
        self,
        y_true: np.ndarray,
        y_prob_sigmoid: np.ndarray,
        y_prob_isotonic: np.ndarray
    ) -> Dict[str, CalibrationMetrics]:
        """
        Compare calibration quality of sigmoid vs isotonic methods.

        Args:
            y_true: True binary labels
            y_prob_sigmoid: Probabilities from sigmoid calibration
            y_prob_isotonic: Probabilities from isotonic calibration

        Returns:
            Dictionary with metrics for each method
        """
        results = {
            'sigmoid': self.evaluate_model(y_true, y_prob_sigmoid),
            'isotonic': self.evaluate_model(y_true, y_prob_isotonic)
        }

        return results

    def plot_reliability_diagram(
        self,
        y_true: np.ndarray,
        y_prob_sigmoid: np.ndarray,
        y_prob_isotonic: np.ndarray,
        save_path: Optional[str] = None,
        title: str = "Calibration Reliability Diagram"
    ) -> None:
        """
        Plot reliability diagram comparing sigmoid vs isotonic calibration.

        A reliability diagram plots predicted probability vs observed frequency.
        Perfect calibration = diagonal line.

        Args:
            y_true: True binary labels
            y_prob_sigmoid: Probabilities from sigmoid calibration
            y_prob_isotonic: Probabilities from isotonic calibration
            save_path: Optional path to save figure
            title: Plot title
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Sigmoid calibration
        prob_true_sig, prob_pred_sig = calibration_curve(
            y_true, y_prob_sigmoid, n_bins=self.n_bins
        )
        ax1.plot(prob_pred_sig, prob_true_sig, marker='o', linewidth=2, label='Sigmoid')
        ax1.plot([0, 1], [0, 1], 'k--', label='Perfect calibration')
        ax1.set_xlabel('Predicted probability', fontsize=12)
        ax1.set_ylabel('Observed frequency', fontsize=12)
        ax1.set_title('Sigmoid Calibration', fontsize=14)
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Isotonic calibration
        prob_true_iso, prob_pred_iso = calibration_curve(
            y_true, y_prob_isotonic, n_bins=self.n_bins
        )
        ax2.plot(prob_pred_iso, prob_true_iso, marker='o', linewidth=2,
                 label='Isotonic', color='orange')
        ax2.plot([0, 1], [0, 1], 'k--', label='Perfect calibration')
        ax2.set_xlabel('Predicted probability', fontsize=12)
        ax2.set_ylabel('Observed frequency', fontsize=12)
        ax2.set_title('Isotonic Calibration', fontsize=14)
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        fig.suptitle(title, fontsize=16, fontweight='bold')
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Reliability diagram saved to: {save_path}")
        else:
            plt.show()

        plt.close()

    def print_comparison(
        self,
        metrics: Dict[str, CalibrationMetrics],
        verbose: bool = True
    ) -> None:
        """
        Print formatted comparison of calibration metrics.

        Args:
            metrics: Dictionary with 'sigmoid' and 'isotonic' metrics
            verbose: Print detailed breakdown
        """
        print("\n" + "="*80)
        print("CALIBRATION QUALITY COMPARISON")
        print("="*80)

        if verbose:
            print("\nSigmoid Calibration:")
            print(metrics['sigmoid'])
            print("\nIsotonic Calibration:")
            print(metrics['isotonic'])

        print("\n" + "-"*80)
        print("SUMMARY (lower = better for all metrics)")
        print("-"*80)
        print(f"{'Metric':<20} {'Sigmoid':<15} {'Isotonic':<15} {'Winner':<10}")
        print("-"*80)

        # Log-loss
        sig_loss = metrics['sigmoid'].log_loss
        iso_loss = metrics['isotonic'].log_loss
        winner_loss = 'Sigmoid' if sig_loss < iso_loss else 'Isotonic'
        print(f"{'Log-loss':<20} {sig_loss:<15.6f} {iso_loss:<15.6f} {winner_loss:<10}")

        # Brier score
        sig_brier = metrics['sigmoid'].brier_score
        iso_brier = metrics['isotonic'].brier_score
        winner_brier = 'Sigmoid' if sig_brier < iso_brier else 'Isotonic'
        print(f"{'Brier score':<20} {sig_brier:<15.6f} {iso_brier:<15.6f} {winner_brier:<10}")

        # ECE
        sig_ece = metrics['sigmoid'].ece
        iso_ece = metrics['isotonic'].ece
        winner_ece = 'Sigmoid' if sig_ece < iso_ece else 'Isotonic'
        print(f"{'ECE':<20} {sig_ece:<15.6f} {iso_ece:<15.6f} {winner_ece:<10}")

        # MCE
        sig_mce = metrics['sigmoid'].mce
        iso_mce = metrics['isotonic'].mce
        winner_mce = 'Sigmoid' if sig_mce < iso_mce else 'Isotonic'
        print(f"{'MCE':<20} {sig_mce:<15.6f} {iso_mce:<15.6f} {winner_mce:<10}")

        print("-"*80)

        # Determine overall winner (by majority vote)
        winners = [winner_loss, winner_brier, winner_ece, winner_mce]
        sig_wins = winners.count('Sigmoid')
        iso_wins = winners.count('Isotonic')

        if sig_wins > iso_wins:
            overall = 'Sigmoid'
        elif iso_wins > sig_wins:
            overall = 'Isotonic'
        else:
            overall = 'Tie'

        print(f"\n✓ Overall Winner: {overall} ({sig_wins} vs {iso_wins} metrics)")
        print("="*80 + "\n")
