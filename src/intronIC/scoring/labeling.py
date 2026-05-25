"""Unified U12/U2 labeling on three orthogonal axes.

Replaces the prior inconsistent type_id assignment paths (predictor.py:129
used a 50% raw-classifier boundary; main.py:3197/3301 used the configurable
threshold default 90%). All three are superseded by the labels emitted here.

Three label axes:
    type_id    in {"u12", "u2"}                   binary call on final score
    confidence in {"strong", "borderline"}        within-call gradation
    history    in {"stable", "promoted", "demoted"}  score path through pipeline

Decision logic (v2.7.1 defaults):
    type_id = "u12" if adjusted_score >= 50 else "u2"
    confidence:
        type_id=u12 → "strong" if adjusted_score >= 90 else "borderline"
        type_id=u2  → "strong" if first_pass_svm <  10 else "borderline"
    history:
        u12 + first_pass_svm <  50 → "promoted"   (mode-sep rescued)
        u2  + first_pass_svm >= 50 → "demoted"    (discount suppressed)
        otherwise                  → "stable"

Thresholds are derived empirically from a 73-species IPA-gold panel
(see WtMTA_db/v2.7_runs/diptera_ipa_check/labeling_thresholds.tsv):
    u2_strong (<10) captures 99.92% of true U2s with 0.02% U12 loss.
    Demoted introns have a median gap of ~70 score-units — no minimum
    gap filter needed; all observed demotions are meaningful.

For Tier F species (gate-fail, no second-pass), first_pass_svm is absent
from score_info.iic; callers should pass svm_score as the first_pass_svm
fallback. In that regime mode-sep didn't run, so promoted is impossible
by construction.
"""
from dataclasses import dataclass
from typing import Literal

TypeId = Literal["u12", "u2"]
Confidence = Literal["strong", "borderline"]
History = Literal["stable", "promoted", "demoted"]


# Default thresholds (v2.7.1) — overridable per call if calibration changes.
U12_STRONG_THRESHOLD = 90.0     # adjusted_score for high-confidence U12
U12_CALL_THRESHOLD = 50.0       # adjusted_score for binary U12/U2 boundary
U2_STRONG_THRESHOLD = 10.0      # first_pass_svm below which U2 is "strong"
PROMOTED_DEMOTED_THRESHOLD = 50.0  # first_pass_svm boundary for "did first-pass call it U12?"


@dataclass(frozen=True)
class IntronLabels:
    type_id: TypeId
    confidence: Confidence
    history: History


def assign_labels(
    first_pass_svm: float,
    adjusted_score: float,
    u12_strong_threshold: float = U12_STRONG_THRESHOLD,
    u12_call_threshold: float = U12_CALL_THRESHOLD,
    u2_strong_threshold: float = U2_STRONG_THRESHOLD,
    promoted_demoted_threshold: float = PROMOTED_DEMOTED_THRESHOLD,
) -> IntronLabels:
    """Compute (type_id, confidence, history) for a single intron.

    Args:
        first_pass_svm: first-pass SVM probability (0-100). For Tier F
            species (no second-pass), pass svm_score as fallback — promoted
            will then never fire (correctly: there was no mode-sep to do it).
        adjusted_score: final post-discount score (0-100).
        u12_strong_threshold: adjusted_score floor for "strong U12" (default 90).
        u12_call_threshold: adjusted_score floor for "is U12" (default 50).
        u2_strong_threshold: first_pass_svm ceiling for "strong U2" (default 10).
        promoted_demoted_threshold: first_pass_svm boundary for the history
            axis (default 50, matching the original raw-classifier boundary).

    Returns:
        IntronLabels with all three fields set.
    """
    if adjusted_score >= u12_call_threshold:
        type_id: TypeId = "u12"
        confidence: Confidence = (
            "strong" if adjusted_score >= u12_strong_threshold else "borderline"
        )
        history: History = (
            "promoted" if first_pass_svm < promoted_demoted_threshold else "stable"
        )
    else:
        type_id = "u2"
        confidence = (
            "strong" if first_pass_svm < u2_strong_threshold else "borderline"
        )
        history = (
            "demoted" if first_pass_svm >= promoted_demoted_threshold else "stable"
        )
    return IntronLabels(type_id=type_id, confidence=confidence, history=history)


def count_by_label(intron_labels) -> dict:
    """Aggregate counts across an iterable of IntronLabels.

    Returns:
        Dict with keys: u12_count, u12_strong_count, u12_borderline_count,
            u12_promoted_count, u2_count, u2_strong_count, u2_borderline_count,
            u2_demoted_count. All zero by default; populated as the iterable
            is consumed.
    """
    counts = {
        "u12_count": 0, "u12_strong_count": 0, "u12_borderline_count": 0,
        "u12_promoted_count": 0,
        "u2_count": 0, "u2_strong_count": 0, "u2_borderline_count": 0,
        "u2_demoted_count": 0,
    }
    for lab in intron_labels:
        if lab.type_id == "u12":
            counts["u12_count"] += 1
            counts[f"u12_{lab.confidence}_count"] += 1
            if lab.history == "promoted":
                counts["u12_promoted_count"] += 1
        else:
            counts["u2_count"] += 1
            counts[f"u2_{lab.confidence}_count"] += 1
            if lab.history == "demoted":
                counts["u2_demoted_count"] += 1
    return counts
