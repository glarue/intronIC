"""
Global progress tracking for training pipeline.

Provides unified step counting across nested CV, optimization, and ensemble training.
"""

from dataclasses import dataclass
from typing import Optional


@dataclass
class ProgressTracker:
    """
    Track global progress across entire training pipeline.

    Attributes:
        total_steps: Total number of steps in pipeline
        current_step: Current step number (0-indexed internally, 1-indexed for display)
        verbose: Whether to print progress updates
    """

    total_steps: int
    current_step: int = 0
    verbose: bool = True

    def increment(self, description: str = "") -> None:
        """
        Increment step counter and optionally print progress.

        Args:
            description: Optional description of current step
        """
        self.current_step += 1
        if self.verbose and description:
            self._print_progress(description)

    def _print_progress(self, description: str) -> None:
        """Print current progress with description."""
        percent = int((self.current_step / self.total_steps) * 100)
        progress_str = f"[{percent}% complete]"
        print(f"{progress_str} {description}", flush=True)

    def format_step(self, local_step: Optional[int] = None, local_total: Optional[int] = None) -> str:
        """
        Format step indicator for progress messages.

        Args:
            local_step: Current step within local phase (e.g., round 3 of optimization)
            local_total: Total steps in local phase (e.g., 5 optimization rounds)

        Returns:
            Formatted string like "Round 3/5 (15% complete)"
        """
        # Calculate percentage based on the NEXT step (current_step + 1)
        # since this is typically called before the step is actually incremented
        percent = int(((self.current_step + 1) / self.total_steps) * 100)
        return f"({percent}% complete)"

    @staticmethod
    def calculate_total_steps(
        eval_mode: str,
        n_cv_folds: int,
        n_optimization_rounds: int,
        n_ensemble_models: int,
        skip_final_optimization: bool = False
    ) -> int:
        """
        Calculate total number of trackable steps in pipeline.

        Args:
            eval_mode: Evaluation mode ('nested_cv', 'split', 'none')
            n_cv_folds: Number of CV folds (for nested_cv)
            n_optimization_rounds: Optimization rounds
            n_ensemble_models: Models in ensemble (for final production model)
            skip_final_optimization: Whether final optimization is skipped (fold-averaged params)

        Returns:
            Total number of steps

        Note:
            During evaluation (nested CV or split), we only train 1 model per fold
            for speed, but the final production model uses the full n_ensemble_models.
        """
        total = 0

        # Phase 1: Evaluation (if enabled)
        # IMPORTANT: Evaluation uses n_ensemble_models=1 for speed (see classifier.py)
        if eval_mode == 'nested_cv':
            # Per fold: n_rounds (optimization) + 1 model (ensemble) + 1 (evaluation)
            steps_per_fold = n_optimization_rounds + 1 + 1
            total += n_cv_folds * steps_per_fold
        elif eval_mode == 'split':
            # Single split: optimization + 1 model training + evaluation
            total += n_optimization_rounds + 1 + 1
        # eval_mode == 'none': no evaluation steps

        # Phase 2: Production model (uses full n_ensemble_models)
        if not skip_final_optimization:
            total += n_optimization_rounds  # Final optimization
        total += n_ensemble_models  # Final ensemble training

        return total
