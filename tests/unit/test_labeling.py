"""Unit tests for intronIC.scoring.labeling — three-axis label derivation."""
from intronIC.scoring.labeling import assign_labels, count_by_label


def test_strong_u12():
    L = assign_labels(first_pass_svm=99, adjusted_score=98)
    assert L.type_id == "u12"
    assert L.confidence == "strong"
    assert L.history == "stable"


def test_borderline_u12():
    L = assign_labels(first_pass_svm=80, adjusted_score=70)
    assert L.type_id == "u12"
    assert L.confidence == "borderline"
    assert L.history == "stable"


def test_promoted_u12_borderline():
    # mode-sep rescue: first-pass below 50, adjusted between 50 and 90
    L = assign_labels(first_pass_svm=20, adjusted_score=70)
    assert L.type_id == "u12"
    assert L.confidence == "borderline"
    assert L.history == "promoted"


def test_promoted_u12_strong():
    # mode-sep aggressive rescue: first-pass below 50, adjusted reaches strong
    L = assign_labels(first_pass_svm=20, adjusted_score=95)
    assert L.type_id == "u12"
    assert L.confidence == "strong"
    assert L.history == "promoted"


def test_demoted_u2():
    # discount caught overcall: first-pass said U12, adjusted is U2
    L = assign_labels(first_pass_svm=80, adjusted_score=40)
    assert L.type_id == "u2"
    assert L.confidence == "borderline"
    assert L.history == "demoted"


def test_borderline_u2():
    # first-pass moderate, adjusted slightly below threshold
    L = assign_labels(first_pass_svm=30, adjusted_score=30)
    assert L.type_id == "u2"
    assert L.confidence == "borderline"
    assert L.history == "stable"


def test_strong_u2():
    # confidently U2 throughout
    L = assign_labels(first_pass_svm=0, adjusted_score=0)
    assert L.type_id == "u2"
    assert L.confidence == "strong"
    assert L.history == "stable"


def test_u2_boundary_at_strong_threshold():
    # first_pass_svm exactly at u2_strong threshold (10) → borderline
    L = assign_labels(first_pass_svm=10, adjusted_score=10)
    assert L.type_id == "u2"
    assert L.confidence == "borderline"


def test_u12_call_boundary():
    # adjusted_score exactly at 50 → u12 (>= boundary)
    L = assign_labels(first_pass_svm=50, adjusted_score=50)
    assert L.type_id == "u12"
    assert L.confidence == "borderline"
    assert L.history == "stable"


def test_count_by_label_aggregation():
    labels = [
        assign_labels(99, 98),   # u12 strong stable
        assign_labels(20, 70),   # u12 borderline promoted
        assign_labels(80, 40),   # u2 borderline demoted
        assign_labels(0, 0),     # u2 strong stable
    ]
    counts = count_by_label(labels)
    assert counts["u12_count"] == 2
    assert counts["u12_strong_count"] == 1
    assert counts["u12_borderline_count"] == 1
    assert counts["u12_promoted_count"] == 1
    assert counts["u2_count"] == 2
    assert counts["u2_strong_count"] == 1
    assert counts["u2_borderline_count"] == 1
    assert counts["u2_demoted_count"] == 1


def test_tier_f_fallback_with_first_pass_eq_adjusted():
    """For Tier F species, caller passes svm_score as first_pass_svm fallback.
    Promoted is impossible (would require first_pass < adjusted, but they're equal)."""
    # Strong u12 in Tier F (rare)
    L = assign_labels(first_pass_svm=95, adjusted_score=95)
    assert L.type_id == "u12" and L.confidence == "strong" and L.history == "stable"
    # Strong u2 in Tier F
    L = assign_labels(first_pass_svm=5, adjusted_score=5)
    assert L.type_id == "u2" and L.confidence == "strong" and L.history == "stable"
    # Demoted in Tier F (only possible via discount, since mode-sep didn't run)
    L = assign_labels(first_pass_svm=70, adjusted_score=30)
    assert L.type_id == "u2" and L.history == "demoted"
