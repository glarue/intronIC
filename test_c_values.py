#!/usr/bin/env python3
"""
Direct test of C values to determine if optimizer is finding wrong C
or if LinearSVC+calibration fundamentally doesn't work.

Tests:
- C=9.3e-05 (current optimal found by optimizer)
- C=1000 (approximately where original SVC finds optimal)
"""

import numpy as np
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import f1_score, average_precision_score, balanced_accuracy_score
import gzip
import sys

def load_reference_sequences(filepath):
    """Load reference sequences from .iic.gz file."""
    introns = []
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                introns.append(parts[1])  # Intron sequence
    return introns

def extract_splice_sites(seq):
    """Extract 5' splice site, branch point region, and 3' splice site."""
    if len(seq) < 50:
        return None

    # 5' splice site: -3 to +9 relative to donor (12 bp total)
    five_seq = seq[:12] if len(seq) >= 12 else seq

    # 3' splice site: -6 to +4 relative to acceptor (10 bp total)
    three_seq = seq[-10:] if len(seq) >= 10 else seq

    # Branch point: -55 to -5 relative to acceptor (50 bp region)
    bp_start = max(0, len(seq) - 55)
    bp_end = max(0, len(seq) - 5)
    bp_seq = seq[bp_start:bp_end] if bp_end > bp_start else seq

    return five_seq, bp_seq, three_seq

def simple_pwm_score(seq, motif_type='u12'):
    """
    Simple PWM-like scoring based on known U12/U2 motifs.
    This is a placeholder - actual scoring would use the PWM matrices.
    """
    if motif_type == 'u12':
        # U12 5' donor: ATATCCTT
        # U12 3' acceptor: strong AG
        score = 0
        if 'ATAT' in seq[:10]:
            score += 2
        if 'CCTT' in seq[:10]:
            score += 2
        if 'AG' in seq[-10:]:
            score += 1
    else:
        # U2 5' donor: GT
        # U2 3' acceptor: AG
        score = 0
        if seq.startswith('GT'):
            score += 1
        if seq.endswith('AG'):
            score += 1

    return score

def create_features_from_sequences(u12_seqs, u2_seqs):
    """
    Create feature matrix from sequences.
    Using simplified scoring since we don't have PWM matrices loaded.
    """
    features = []
    labels = []

    # Process U12
    for seq in u12_seqs:
        sites = extract_splice_sites(seq)
        if sites:
            five_seq, bp_seq, three_seq = sites
            # Use sequence composition as proxy features
            five_score = five_seq.count('A') + five_seq.count('T')
            bp_score = bp_seq.count('CTAAC') * 10 + bp_seq.count('CTAA') * 5
            three_score = three_seq.count('AG') * 5
            features.append([five_score, bp_score, three_score])
            labels.append(1)  # U12

    # Process U2
    for seq in u2_seqs:
        sites = extract_splice_sites(seq)
        if sites:
            five_seq, bp_seq, three_seq = sites
            five_score = five_seq.count('G') + five_seq.count('T')
            bp_score = bp_seq.count('YNYURAY') * 2  # Simplified
            three_score = three_seq.count('AG') * 3
            features.append([five_score, bp_score, three_score])
            labels.append(0)  # U2

    return np.array(features), np.array(labels)

def test_c_value(C, X_train, X_test, y_train, y_test, X_full, y_full):
    """Test a specific C value and return metrics."""
    print(f"\n{'='*60}")
    print(f"Testing C = {C}")
    print(f"{'='*60}")

    # Create model
    base_svm = LinearSVC(
        C=C,
        class_weight='balanced',
        loss='squared_hinge',
        penalty='l2',
        dual=False,
        max_iter=10000,
        tol=1e-4,
        random_state=42
    )

    # Calibrate
    cal_svm = CalibratedClassifierCV(
        base_svm,
        method='sigmoid',
        cv=5
    )

    # Train
    print("Training model...")
    cal_svm.fit(X_train, y_train)

    # Predict
    print("Making predictions...")
    y_pred = cal_svm.predict(X_test)
    y_proba = cal_svm.predict_proba(X_test)[:, 1]

    # Metrics
    f1 = f1_score(y_test, y_pred)
    pr_auc = average_precision_score(y_test, y_proba)

    print(f"\nTest Set Results:")
    print(f"  F1 Score: {f1:.4f}")
    print(f"  PR-AUC: {pr_auc:.4f}")
    print(f"  U12s predicted: {y_pred.sum()}/{(y_test == 1).sum()} actual")

    # Cross-validation on full dataset
    print("\nCross-validation (5-fold) on full dataset...")
    cv_scores = cross_val_score(
        cal_svm, X_full, y_full,
        cv=5,
        scoring='balanced_accuracy',
        n_jobs=1
    )
    print(f"  Mean CV Balanced Accuracy: {cv_scores.mean():.4f} (+/- {cv_scores.std():.4f})")

    return {
        'C': C,
        'f1': f1,
        'pr_auc': pr_auc,
        'cv_mean': cv_scores.mean(),
        'cv_std': cv_scores.std()
    }

def main():
    print("Loading reference sequences...")

    # Load U12 and U2 references
    u12_path = 'intronIC/data/u12_reference.introns.iic.gz'
    u2_path = 'intronIC/data/u2_reference.introns.iic.gz'

    try:
        u12_seqs = load_reference_sequences(u12_path)
        u2_seqs = load_reference_sequences(u2_path)
    except FileNotFoundError as e:
        print(f"Error: Could not find reference files: {e}")
        print("Please ensure you're running from /home/user/intronIC")
        return

    print(f"Loaded {len(u12_seqs)} U12 and {len(u2_seqs)} U2 sequences")

    # Create features
    print("\nExtracting features...")
    X, y = create_features_from_sequences(u12_seqs[:387], u2_seqs[:20690])

    print(f"Feature matrix shape: {X.shape}")
    print(f"Class distribution: {(y == 0).sum()} U2, {(y == 1).sum()} U12")

    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y,
        test_size=0.2,
        stratify=y,
        random_state=42
    )

    print(f"\nTrain set: {len(y_train)} samples")
    print(f"Test set: {len(y_test)} samples")

    # Test both C values
    results = []

    # C found by current optimizer
    results.append(test_c_value(9.321564e-05, X_train, X_test, y_train, y_test, X, y))

    # C approximately where original finds optimal
    results.append(test_c_value(1000, X_train, X_test, y_train, y_test, X, y))

    # Also test a few other values for context
    for C in [0.001, 0.01, 0.1, 1.0, 10, 100]:
        results.append(test_c_value(C, X_train, X_test, y_train, y_test, X, y))

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"{'C Value':<15} {'F1':<8} {'PR-AUC':<8} {'CV Bal. Acc.':<15}")
    print("-"*60)
    for r in results:
        print(f"{r['C']:<15.2e} {r['f1']:<8.4f} {r['pr_auc']:<8.4f} {r['cv_mean']:<15.4f}")

    # Find best
    best_by_cv = max(results, key=lambda x: x['cv_mean'])
    print(f"\nBest C by CV score: {best_by_cv['C']:.2e} (CV={best_by_cv['cv_mean']:.4f})")

    best_by_f1 = max(results, key=lambda x: x['f1'])
    print(f"Best C by F1 score: {best_by_f1['C']:.2e} (F1={best_by_f1['f1']:.4f})")

if __name__ == '__main__':
    main()
