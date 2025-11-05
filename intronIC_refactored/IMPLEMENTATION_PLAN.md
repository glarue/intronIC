# Scoring & Classification Implementation Plan

**Based on:** SCORING_ANALYSIS.md
**Goal:** Implement M1.4 (Scoring) and M1.5 (Classification) while fixing identified issues
**Approach:** Build it right from the start, don't replicate bad patterns

---

## Table of Contents

1. [Overall Strategy](#overall-strategy)
2. [Phase 1: M1.4 Scoring System](#phase-1-m14-scoring-system)
3. [Phase 2: M1.5 Classification System](#phase-2-m15-classification-system)
4. [Phase 3: Integration & Validation](#phase-3-integration--validation)

---

## Overall Strategy

### Key Principles

1. **Build clean from the start** - Don't replicate issues, even temporarily
2. **Test-first for ML integrity** - Write tests for data leakage prevention BEFORE implementation
3. **Port the good, fix the bad** - Keep excellent ensemble approach, fix normalization issues
4. **Document the rationale** - Explain WHY we do things (geometric mean, ensemble, etc.)
5. **Incremental validation** - Compare to original at each step

### What to Port vs. What to Fix

**Port Directly (These are excellent):**
- ✅ PWM scoring logic (`seq_score`, `bp_score`)
- ✅ Ensemble training with U2 subsampling
- ✅ Stratified train/test split
- ✅ Hyperparameter optimization with geometric refinement
- ✅ Rank-1 parameter averaging
- ✅ F1-weighted ensemble predictions
- ✅ Evaluation metrics (F1, PR-AUC)

**Fix During Port:**
- ⚠️ Remove post-classification re-normalization
- ⚠️ Simplify recursive scoring (or skip for now)
- ⚠️ Use StandardScaler consistently
- ⚠️ Add data leakage prevention tests

**Skip for MVP (M1.4 + M1.5):**
- ❌ Recursive scoring implementation (save for M1.6 - see decision below)
- ❌ Clustering features (not used in main pipeline)
- ❌ Performance optimizations (do later)

**Recursive Scoring Decision: Option 2 ✅**

**Decision:** Design API now, implement in M1.6 (future phase)

**Rationale:**
- User hasn't used it much in practice (human PWMs work well)
- Current implementation has 3 issues that need careful fixes
- Adds significant complexity to initial validation
- Easy to add later with proper API hooks
- Can prioritize based on actual need after core pipeline validated

**What We'll Do Now:**
- ✅ Design swappable PWM interface
- ✅ Add `PWMBuilder` stub for building matrices from sequences
- ✅ Add `RecursiveScorer` stub with clean two-pass interface
- ✅ Document TODOs for M1.6

**What We'll Do Later (M1.6):**
- Build PWMs from high-confidence predictions
- Implement clean two-pass classification (fixes Issues #2, #3)
- Validate on non-model organisms

---

## Phase 1: M1.4 Scoring System

**Goal:** PWM scoring + normalization (WITHOUT the data leakage issue)

### Step 1: Set Up Module Structure

**Commands:**
```bash
cd intronIC_refactored
mkdir -p scoring
touch scoring/__init__.py
touch scoring/pwm.py
touch scoring/branch_point.py
touch scoring/normalizer.py
touch scoring/scorer.py

mkdir -p tests/unit/scoring
touch tests/unit/scoring/__init__.py
touch tests/unit/scoring/test_pwm.py
touch tests/unit/scoring/test_normalizer.py
```

**Time:** 5 minutes

---

### Step 2: Write ML Integrity Tests FIRST

**File:** `tests/unit/scoring/test_normalizer.py`

**Why First?** These tests ensure we DON'T replicate Issue #1.

**Critical Tests:**
```python
def test_normalizer_only_sees_training_data():
    """Normalizer should ONLY be fit on reference/training sequences."""
    # This test will fail if we accidentally fit on experimental data
    pass

def test_zscore_consistency_through_pipeline():
    """Z-scores should not change after classification."""
    # This test will fail if we do post-classification re-normalization
    pass

def test_normalizer_rejects_experimental_data_in_fit():
    """Normalizer.fit() should raise error if given non-reference data."""
    # Make it impossible to accidentally fit on wrong data
    pass
```

**Reference:** SCORING_ANALYSIS.md, "Test Strategy Recommendations"

**Time:** 30 minutes

**Action:** Let me help you write these tests first, then we implement to pass them.

---

### Step 3: Implement PWM Scoring

**File:** `scoring/pwm.py`

**Reference Original:** `intronIC.py:2114-2200`

**Classes:**
```python
@dataclass(frozen=True, slots=True)
class PWM:
    """Position Weight Matrix for sequence scoring."""
    name: str
    matrix: np.ndarray  # Shape: (4, length) for ACGT
    length: int
    pseudocount: float = 0.0001

    def score_sequence(
        self,
        seq: str,
        ignore_positions: Optional[Set[int]] = None
    ) -> float:
        """
        Score sequence with PWM (log-odds ratio).

        Port from: intronIC.py:2114-2142 (seq_score)
        """
        pass

@dataclass(frozen=True, slots=True)
class PWMSet:
    """U2 and U12 PWMs for a splice site type."""
    u2_canonical: PWM
    u2_noncanonical: Optional[PWM]
    u12_canonical: PWM
    u12_noncanonical: Optional[PWM]

class PWMLoader:
    """Load PWMs from scoring_matrices.fasta.iic format."""

    @staticmethod
    def load_from_file(filepath: Path) -> Dict[str, PWMSet]:
        """
        Load all PWMs from file.

        Returns: {'five': PWMSet, 'bp': PWMSet, 'three': PWMSet}
        """
        pass
```

**Port Strategy:**
- Copy `seq_score` logic line-by-line
- Verify log-odds calculation is correct
- Add type hints
- Write unit tests

**Time:** 3-4 hours

---

### Step 4: Implement Branch Point Detection

**File:** `scoring/branch_point.py`

**Reference Original:** `intronIC.py:2143-2200`

**Classes:**
```python
@dataclass(frozen=True, slots=True)
class BranchPointMatch:
    """Result from branch point search."""
    sequence: str
    score: float
    position: int  # Relative to intron 3' end
    coordinate: GenomicCoordinate  # Absolute position

class BranchPointScorer:
    """Find optimal branch point location."""

    def __init__(self, u12_pwm: PWM, u2_pwm: PWM):
        self.u12_pwm = u12_pwm
        self.u2_pwm = u2_pwm

    def find_best_match(
        self,
        intron: Intron,
        search_window: Tuple[int, int] = (-55, -5)
    ) -> BranchPointMatch:
        """
        Sliding window search for best branch point.

        Port from: intronIC.py:2143-2200 (bp_score)
        """
        pass
```

**Port Strategy:**
- Copy sliding window logic
- Return structured result (not just score)
- Track position and coordinates
- Write unit tests

**Time:** 2-3 hours

---

### Step 5: Implement Normalizer (CORRECTLY)

**File:** `scoring/normalizer.py`

**Reference Original:** `intronIC.py:3725-3732` (GOOD part only, NOT 5247-5251)

**Design to Prevent Issue #1:**

```python
from typing import Literal

DatasetType = Literal["reference", "experimental"]

class ScoreNormalizer:
    """
    Normalize PWM scores to z-scores.

    CRITICAL: This class enforces ML best practices by tracking
    which data is used for fitting vs. transforming.
    """

    def __init__(self):
        self._scaler: Optional[StandardScaler] = None
        self._fitted_on: Optional[DatasetType] = None
        self._is_fitted: bool = False

    def fit(
        self,
        introns: Iterable[Intron],
        dataset_type: DatasetType = "reference"
    ) -> "ScoreNormalizer":
        """
        Fit scaler on intron scores.

        Args:
            introns: Should be REFERENCE introns only for proper ML practice
            dataset_type: Must be "reference" to prevent accidental misuse

        Raises:
            ValueError: If dataset_type is "experimental" (prevents Issue #1)

        Port from: intronIC.py:3727-3731
        DO NOT PORT: intronIC.py:5247-5251 (this is the bad re-normalization)
        """
        if dataset_type == "experimental":
            raise ValueError(
                "Cannot fit normalizer on experimental data! "
                "This would cause data leakage. Use dataset_type='reference'."
            )

        # Extract raw scores into matrix
        score_matrix = self._extract_scores(introns)

        # Fit StandardScaler (NOT RobustScaler - Issue #3)
        self._scaler = StandardScaler().fit(score_matrix)
        self._fitted_on = dataset_type
        self._is_fitted = True

        return self

    def transform(
        self,
        introns: Iterable[Intron],
        dataset_type: DatasetType
    ) -> Iterator[Intron]:
        """
        Transform intron scores to z-scores.

        Can be called on either reference or experimental data,
        but scaler must have been fit on reference data first.

        Port from: intronIC.py:3655-3665 (scale_scores)
        """
        if not self._is_fitted:
            raise RuntimeError("Must call fit() before transform()")

        score_matrix = self._extract_scores(introns)
        z_score_matrix = self._scaler.transform(score_matrix)

        # Yield introns with updated z-scores
        for intron, z_scores in zip(introns, z_score_matrix):
            yield self._update_intron_scores(intron, z_scores)

    def fit_transform(
        self,
        introns: Iterable[Intron],
        dataset_type: DatasetType = "reference"
    ) -> Iterator[Intron]:
        """Convenience method: fit and transform in one step."""
        self.fit(introns, dataset_type)
        return self.transform(introns, dataset_type)
```

**Key Design Decisions:**

1. **Explicit `dataset_type` parameter** - Makes it impossible to accidentally fit on experimental data
2. **Raises error if fitting on experimental** - Prevents Issue #1 at runtime
3. **Always uses StandardScaler** - Fixes Issue #3 (consistency)
4. **Tracks what it was fitted on** - Helps with debugging and validation
5. **Immutable introns** - Returns new introns with updated scores (functional style)

**Tests to Write:**
```python
def test_cannot_fit_on_experimental_data():
    """Should raise ValueError if trying to fit on experimental data."""
    normalizer = ScoreNormalizer()
    with pytest.raises(ValueError, match="data leakage"):
        normalizer.fit(experimental_introns, dataset_type="experimental")

def test_can_transform_experimental_after_fitting_reference():
    """Can transform experimental data after fitting on reference."""
    normalizer = ScoreNormalizer()
    normalizer.fit(reference_introns, dataset_type="reference")
    transformed = list(normalizer.transform(
        experimental_introns,
        dataset_type="experimental"
    ))
    assert len(transformed) > 0
```

**Time:** 2-3 hours

---

### Step 6: Implement Scoring Pipeline Orchestrator

**File:** `scoring/scorer.py`

**Reference Original:** `intronIC.py:3115-3400` (get_raw_scores)

**Class:**
```python
class IntronScorer:
    """
    Main scoring pipeline: PWM scoring + normalization.

    This class orchestrates the full scoring process while
    maintaining proper train/test separation.
    """

    def __init__(
        self,
        pwm_sets: Dict[str, PWMSet],
        five_coords: Tuple[int, int] = (-3, 9),
        bp_coords: Tuple[int, int] = (-55, -5),
        three_coords: Tuple[int, int] = (-6, 4),
        pseudocount: float = 0.0001
    ):
        self.pwm_sets = pwm_sets
        self.five_coords = five_coords
        self.bp_coords = bp_coords
        self.three_coords = three_coords
        self.pseudocount = pseudocount

        self.bp_scorer = BranchPointScorer(
            pwm_sets['bp'].u12_canonical,
            pwm_sets['bp'].u2_canonical
        )

    def score_introns(
        self,
        introns: Iterable[Intron],
        ignore_nc_dnts: bool = True
    ) -> Iterator[Intron]:
        """
        Score introns with PWMs (raw scores only, no normalization).

        Port from: intronIC.py:3115-3400

        Returns: Introns with .scores.five_raw_score, .bp_raw_score,
                 .three_raw_score populated
        """
        for intron in introns:
            # Score 5' splice site
            five_score = self._score_five_site(intron, ignore_nc_dnts)

            # Find and score branch point
            bp_match = self.bp_scorer.find_best_match(intron, self.bp_coords)

            # Score 3' splice site
            three_score = self._score_three_site(intron, ignore_nc_dnts)

            # Update intron scores
            updated_scores = replace(
                intron.scores,
                five_raw_score=five_score,
                bp_raw_score=bp_match.score,
                three_raw_score=three_score,
                bp_sequence_u12=bp_match.sequence,
                bp_position=bp_match.position
            )

            yield replace(intron, scores=updated_scores)
```

**Time:** 3-4 hours

---

### M1.4 Validation

**Integration Test:**
```python
def test_scoring_pipeline_chr19():
    """Score chr19 introns and verify distributions."""
    # Load chr19 introns (from M1.3)
    # Load PWMs
    # Create scorer
    # Score all introns
    # Assert: raw scores are reasonable
    # Assert: distributions match expected ranges
```

**Gold Standard Comparison:**
```python
def test_raw_scores_match_original():
    """Raw scores should match original intronIC."""
    # Load same introns scored by original intronIC
    # Score with new implementation
    # Assert: raw scores are nearly identical (within floating point tolerance)
```

**Time:** 2 hours

**Total M1.4 Time:** ~15-18 hours

---

## Phase 2: M1.5 Classification System

**Goal:** SVM training + prediction (with all the GOOD patterns)

### Step 7: Implement Hyperparameter Optimizer

**File:** `classification/optimizer.py`

**Reference Original:** `intronIC.py:5431-5528`

**Port Strategy:** This is EXCELLENT code - port almost exactly

**Classes:**
```python
@dataclass(frozen=True, slots=True)
class SVMParameters:
    """SVM hyperparameters."""
    C: float
    kernel: str = 'linear'
    class_weight: str = 'balanced'
    probability: bool = True
    cache_size: int = 1000

class SVMOptimizer:
    """
    Optimize SVM C parameter via geometric grid search.

    This is a PORT of the excellent optimization logic from
    intronIC.py:5431-5528. The geometric refinement strategy
    is sophisticated and well-designed.
    """

    def optimize(
        self,
        u2_introns: Sequence[Intron],
        u12_introns: Sequence[Intron],
        n_rounds: int = 5,
        cv_folds: int = 5,
        scorer: str = 'balanced_accuracy',
        cv_jobs: int = 1,
        seed: int = 42
    ) -> SVMParameters:
        """
        Multi-round geometric grid search for C parameter.

        Process (from original):
        1. Round 1: Coarse grid from 10^-6 to 10^6
        2. Rounds 2-5: Refine around best value
        3. Use geometric mean of rank-1 parameters

        Port from: intronIC.py:5431-5528
        Keep: Geometric refinement, rank-1 averaging
        """
        # Initial coarse grid
        param_grid = self._create_initial_grid()

        for round_num in range(1, n_rounds + 1):
            logger.info(f"Optimization round {round_num}/{n_rounds}")

            # Grid search with cross-validation
            best_params = self._grid_search_round(
                u2_introns, u12_introns, param_grid,
                cv_folds, scorer, cv_jobs, seed + round_num
            )

            # Get ALL rank-1 parameters (not just first)
            rank1_params = self._extract_rank1_params(best_params)

            # Geometric mean (conservative for log-scale parameters)
            avg_C = gmean(rank1_params)

            # Refine grid around best value
            param_grid = self._refine_grid(param_grid, avg_C)

        return SVMParameters(C=avg_C)

    def _extract_rank1_params(self, cv_results: dict) -> List[float]:
        """
        Extract ALL rank-1 C values, not just first.

        Port from: intronIC.py:5290-5322 (rank1_param_avg, rank_ones)

        Why: sklearn's best_params_ returns FIRST rank-1 value.
        For C parameter, we want GEOMETRIC MEAN of all tied values
        to be conservative (higher C = more conservative).
        """
        ranks = cv_results['rank_test_score']
        params = cv_results['params']
        rank1_C_values = [
            p['C'] for p, rank in zip(params, ranks) if rank == 1
        ]
        return rank1_C_values
```

**Documentation Note:** Add detailed docstring explaining WHY geometric mean of rank-1 params

**Time:** 3-4 hours

---

### Step 8: Implement SVM Trainer

**File:** `classification/trainer.py`

**Reference Original:** `intronIC.py:5345-5428`

**Port Strategy:** The ensemble approach is excellent - port carefully

**Classes:**
```python
@dataclass(frozen=True, slots=True)
class SVMModel:
    """Trained SVM model with performance metrics."""
    classifier: SVC  # The actual sklearn model
    parameters: SVMParameters
    f1_score: float
    pr_auc: float
    precision_recall_curve: Tuple[np.ndarray, np.ndarray, np.ndarray]
    n_training_samples: int
    feature_names: Tuple[str, ...] = ('five_z_score', 'bp_z_score', 'three_z_score')

class SVMTrainer:
    """
    Train SVM ensemble with U2 subsampling.

    Port from: intronIC.py:5345-5428 (train_svm)

    Why ensemble? U2 introns outnumber U12 by ~200:1.
    Subsampling U2 to match U12 count creates balanced training,
    and training multiple models with different subsamples creates
    a robust ensemble that reduces variance.
    """

    def train_ensemble(
        self,
        u2_introns: Sequence[Intron],
        u12_introns: Sequence[Intron],
        parameters: SVMParameters,
        n_models: int = 10,
        train_fraction: float = 0.80,
        seed: int = 42
    ) -> List[SVMModel]:
        """
        Train ensemble of SVM models with U2 subsampling.

        Process (from original):
        1. For each model iteration:
           a. Randomly sample U2 introns (size = len(U12))
           b. Combine with all U12 introns
           c. Stratified 80/20 train/test split
           d. Train SVM
           e. Evaluate on test set (F1, PR-AUC)
        2. Return all trained models

        Port from: intronIC.py:5356-5428
        """
        models = []
        np.random.seed(seed)

        # Create U2 subsamples
        subset_size = len(u12_introns)
        u2_indices = self._create_subsamples(
            len(u2_introns), subset_size, n_models
        )

        for iteration, u2_idx in enumerate(u2_indices, start=1):
            logger.info(f"Training model {iteration}/{n_models}")

            # Subsample U2
            u2_subset = [u2_introns[i] for i in u2_idx]

            # Train single model
            model = self._train_single_model(
                u2_subset, u12_introns, parameters,
                train_fraction, seed
            )

            models.append(model)

        return models

    def _train_single_model(
        self,
        u2_introns: Sequence[Intron],
        u12_introns: Sequence[Intron],
        parameters: SVMParameters,
        train_fraction: float,
        seed: int
    ) -> SVMModel:
        """
        Train single SVM model.

        Port from: intronIC.py:5395-5423
        """
        # Extract feature matrix
        X, y = self._create_feature_matrix(u2_introns, u12_introns)

        # Stratified split (CRITICAL for imbalanced data)
        X_train, X_test, y_train, y_test = train_test_split(
            X, y,
            train_size=train_fraction,
            stratify=y,  # Equal class proportions in train/test
            random_state=seed
        )

        # Train
        clf = SVC(**asdict(parameters))
        clf.fit(X_train, y_train)

        # Evaluate
        y_pred = clf.predict(X_test)
        f1 = f1_score(y_test, y_pred)

        y_prob = clf.predict_proba(X_test)[:, 1]
        precision, recall, thresholds = precision_recall_curve(y_test, y_prob)
        pr_auc = auc(recall, precision)

        return SVMModel(
            classifier=clf,
            parameters=parameters,
            f1_score=f1,
            pr_auc=pr_auc,
            precision_recall_curve=(precision, recall, thresholds),
            n_training_samples=len(X_train)
        )
```

**Time:** 4-5 hours

---

### Step 9: Implement Ensemble Predictor

**File:** `classification/predictor.py`

**Reference Original:** `intronIC.py:5651-5687, 5802-5853`

**Classes:**
```python
class SVMPredictor:
    """
    Apply ensemble of SVM models to classify introns.

    Port from: intronIC.py:5651-5853
    """

    def __init__(self, models: Sequence[SVMModel], use_f1_weights: bool = True):
        self.models = models
        self.use_f1_weights = use_f1_weights

        if use_f1_weights:
            self.weights = np.array([m.f1_score for m in models])
        else:
            self.weights = np.ones(len(models))

    def classify_introns(
        self,
        introns: Iterable[Intron],
        threshold: float = 90.0
    ) -> Iterator[Intron]:
        """
        Classify introns using ensemble predictions.

        Process (from original):
        1. Extract z-score features
        2. Apply each model to get U12 probabilities
        3. Weight by F1 scores
        4. Average probabilities
        5. Compare to threshold
        6. Assign type_id

        Port from: intronIC.py:5651-5687 (average_svm_score_info)
                   intronIC.py:5802-5853 (assign_svm_score, parallel_svm_score)
        """
        for intron in introns:
            # Extract features
            features = self._extract_features(intron)

            # Get predictions from all models
            probabilities = []
            distances = []
            labels = []

            for model in self.models:
                prob = model.classifier.predict_proba([features])[0, 1]
                dist = model.classifier.decision_function([features])[0]
                label = model.classifier.predict([features])[0]

                probabilities.append(prob)
                distances.append(dist)
                labels.append(label)

            # Weighted average
            avg_prob = self._weighted_average(probabilities, self.weights)
            avg_dist = self._weighted_average(distances, self.weights)

            # Convert to percentage
            svm_score = avg_prob * 100

            # Relative to threshold
            relative_score = svm_score - threshold

            # Classify
            type_id: IntronType = 'u12' if svm_score >= threshold else 'u2'

            # Update intron
            updated_scores = replace(
                intron.scores,
                svm_score=svm_score,
                relative_score=relative_score,
                decision_distance=avg_dist
            )

            updated_metadata = replace(
                intron.metadata,
                type_id=type_id
            )

            yield replace(
                intron,
                scores=updated_scores,
                metadata=updated_metadata
            )
```

**Time:** 3-4 hours

---

### Step 10: High-Level Pipeline

**File:** `classification/classifier.py`

**Classes:**
```python
class IntronClassifier:
    """
    Complete classification pipeline.

    Combines scoring, normalization, training, and prediction.
    """

    def __init__(
        self,
        scorer: IntronScorer,
        normalizer: ScoreNormalizer,
        trainer: SVMTrainer,
        optimizer: Optional[SVMOptimizer] = None
    ):
        self.scorer = scorer
        self.normalizer = normalizer
        self.trainer = trainer
        self.optimizer = optimizer

    def train_from_references(
        self,
        u2_ref_file: Path,
        u12_ref_file: Path,
        optimize_hyperparameters: bool = True,
        n_models: int = 10
    ) -> SVMPredictor:
        """
        Train classifier from reference sequence files.

        Process:
        1. Load reference sequences
        2. Score with PWMs
        3. Fit normalizer on references ONLY (prevents Issue #1)
        4. Normalize z-scores
        5. Optimize hyperparameters (if requested)
        6. Train ensemble
        7. Return predictor
        """
        # Load
        u2_introns = self._load_references(u2_ref_file, type_id='u2')
        u12_introns = self._load_references(u12_ref_file, type_id='u12')

        # Score
        u2_scored = list(self.scorer.score_introns(u2_introns))
        u12_scored = list(self.scorer.score_introns(u12_introns))

        # Normalize (fit on references ONLY)
        all_refs = u2_scored + u12_scored
        self.normalizer.fit(all_refs, dataset_type="reference")

        u2_normalized = list(self.normalizer.transform(
            u2_scored, dataset_type="reference"
        ))
        u12_normalized = list(self.normalizer.transform(
            u12_scored, dataset_type="reference"
        ))

        # Optimize
        if optimize_hyperparameters and self.optimizer:
            params = self.optimizer.optimize(u2_normalized, u12_normalized)
        else:
            params = SVMParameters(C=1.0)  # Default

        # Train
        models = self.trainer.train_ensemble(
            u2_normalized, u12_normalized, params, n_models
        )

        return SVMPredictor(models)

    def classify(
        self,
        introns: Iterable[Intron],
        predictor: SVMPredictor,
        threshold: float = 90.0
    ) -> Iterator[Intron]:
        """
        Full classification pipeline for experimental introns.

        Process:
        1. Score with PWMs (if not already scored)
        2. Normalize z-scores (using fitted normalizer)
        3. Classify with predictor
        4. Return classified introns

        IMPORTANT: Z-scores are NOT re-normalized after classification.
        This prevents Issue #1 (data leakage).
        """
        # Score (if needed)
        scored = self.scorer.score_introns(introns)

        # Normalize (transform only, using reference-fitted scaler)
        normalized = self.normalizer.transform(
            scored,
            dataset_type="experimental"
        )

        # Classify
        classified = predictor.classify_introns(normalized, threshold)

        # NO re-normalization here! (Issue #1 fix)
        yield from classified
```

**Time:** 2-3 hours

---

### M1.5 Validation

**Integration Tests:**
```python
def test_full_classification_pipeline_chr19():
    """Train and classify chr19 introns."""
    # Train on references
    # Classify chr19
    # Assert: known U12s are identified (>95% recall)

def test_no_data_leakage_in_normalization():
    """Experimental data should never be used for fitting normalizer."""
    # This is caught by design (ValueError in normalizer.fit)
    # Test that pipeline respects this

def test_zscore_consistency():
    """Z-scores should not change after classification."""
    # Score and normalize introns
    # Save z-scores
    # Classify
    # Assert: z-scores unchanged
```

**Gold Standard:**
```python
def test_classification_matches_original():
    """Classification should match original intronIC (after fixes)."""
    # Note: Will differ slightly due to Issue #1 fix
    # But classifications should be nearly identical
```

**Time:** 3 hours

**Total M1.5 Time:** ~18-22 hours

---

## Phase 3: Integration & Validation

### Step 11: End-to-End Pipeline Test

```python
def test_complete_pipeline():
    """
    Full pipeline from annotation to classification.

    Steps:
    1. Load chr19 annotation + genome (M1.3)
    2. Extract introns (M1.3)
    3. Filter introns (M1.3)
    4. Load PWMs (M1.4)
    5. Score introns (M1.4)
    6. Train classifier on references (M1.5)
    7. Classify chr19 introns (M1.5)
    8. Write outputs (M1.2)
    9. Validate outputs match expected format
    """
```

**Time:** 2 hours

---

### Step 12: Comparison with Original

**Gold Standard Validation:**
```bash
# Run original intronIC on chr19
cd /home/glarue/code/intronIC
python -m intronIC.intronIC \
    -g test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
    -a test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz \
    -n homo_sapiens_original \
    --no_recursive

# Run refactored intronIC on chr19
cd /home/glarue/code/intronIC/intronIC_refactored
python -m main \
    -g test_data/Homo_sapiens.Chr19.Ensembl_91.fa.gz \
    -a test_data/Homo_sapiens.Chr19.Ensembl_91.gff3.gz \
    -n homo_sapiens_refactored

# Compare outputs
python compare_outputs.py \
    homo_sapiens_original.meta.iic \
    homo_sapiens_refactored.meta.iic
```

**Expected Differences:**
- Z-scores will differ (we fixed Issue #1 - no post-classification re-normalization)
- Classifications should be nearly identical (same models, same logic)
- SVM scores should be identical (same ensemble averaging)

**Time:** 2-3 hours

---

## Implementation Timeline

### Week 1: M1.4 Scoring
- **Day 1:** Setup, ML integrity tests, PWM scoring (8 hours)
- **Day 2:** Branch point, normalizer (8 hours)
- **Day 3:** Scorer orchestrator, validation (8 hours)

### Week 2: M1.5 Classification
- **Day 4:** Optimizer, trainer (8 hours)
- **Day 5:** Predictor, high-level pipeline (8 hours)
- **Day 6:** Integration tests, validation (8 hours)

### Week 3: Polish
- **Day 7:** Gold standard comparison, documentation (8 hours)

**Total:** ~40-50 hours

---

## Key Success Criteria

### M1.4 Complete When:
- [ ] All ML integrity tests pass
- [ ] Raw PWM scores match original (±0.001)
- [ ] Normalizer CANNOT be fit on experimental data (enforced by API)
- [ ] Z-scores have μ≈0, σ≈1 for reference data
- [ ] 268 + new tests passing

### M1.5 Complete When:
- [ ] Hyperparameter optimization converges
- [ ] Ensemble training succeeds (F1 >0.90 per model)
- [ ] Classifications match original (>98% agreement)
- [ ] Known U12s identified (>95% recall on chr19)
- [ ] Z-scores unchanged after classification (Issue #1 fixed)
- [ ] All tests passing

### Integration Complete When:
- [ ] Full pipeline runs on chr19 (<5 minutes)
- [ ] All output formats match expected
- [ ] Gold standard validation passes
- [ ] Documentation complete

---

## Commands Reference

```bash
# Setup
cd intronIC_refactored
mkdir -p scoring classification tests/unit/scoring tests/unit/classification

# Run tests during development
PYTHONPATH=. pixi run pytest tests/unit/scoring/ -v
PYTHONPATH=. pixi run pytest tests/unit/classification/ -v
PYTHONPATH=. pixi run pytest tests/integration/ -v

# Run specific test
PYTHONPATH=. pixi run pytest tests/unit/scoring/test_normalizer.py::test_cannot_fit_on_experimental_data -v

# Check coverage
PYTHONPATH=. pixi run pytest --cov=scoring --cov=classification --cov-report=html

# Gold standard comparison
PYTHONPATH=. pixi run python scripts/compare_to_original.py
```

---

## Notes

**What Makes This Different from Original:**
1. ✅ No post-classification re-normalization (Issue #1 fixed)
2. ✅ Normalizer API prevents data leakage by design
3. ✅ Consistent scaler choice (StandardScaler)
4. ✅ Modular, testable architecture
5. ✅ Full type hints
6. ✅ Comprehensive tests
7. ✅ Documentation of ML rationale

**What Stays the Same:**
1. ✅ PWM scoring logic (excellent)
2. ✅ Ensemble approach (excellent)
3. ✅ Hyperparameter optimization (excellent)
4. ✅ Geometric mean for rank-1 params (excellent)
5. ✅ F1-weighted predictions (excellent)
6. ✅ Evaluation metrics (excellent)

**What We Skip (For Now):**
1. ❌ Recursive scoring (has issues, complexity not worth it for MVP)
2. ❌ Clustering (not used in main pipeline)
3. ❌ Parallel processing (add after correctness validated)

---

## Next Session Checklist

Ready to start when you have:
- [ ] Read SCORING_ANALYSIS.md
- [ ] Read this implementation plan
- [ ] Understand the Issues #1, #2, #3
- [ ] Understand why the normalizer API is designed to prevent Issue #1
- [ ] Ready to write ML integrity tests first
- [ ] Coffee ☕

**Let's build it right!** 🚀
