# CLI Improvements Summary (2025-01-17)

This document summarizes the CLI usability improvements made during this session.

## Issues Fixed

### 1. Train Subcommand Not Working
**Issue:** Running `intronIC train -n species_name` incorrectly tried to load a pretrained model instead of training a new one.

**Error:**
```
Error: Default pretrained model not found at data/default_pretrained.model.pkl.
Use --train to train a new model instead.
```

**Root Cause:** `cli/config.py:231` only checked for `args.train` flag (classify subcommand), but didn't account for `args.command == 'train'` (train subcommand).

**Fix:** Updated `cli/config.py:231-242` to check both:
```python
is_training = (
    getattr(args, 'command', None) == 'train' or
    getattr(args, 'train', False)
)
```

**Files Modified:**
- `cli/config.py:231-242`

---

### 2. Raw Rich Markup in Log Files
**Issue:** Log files contained literal Rich markup strings like `[bold blue]` instead of proper ANSI color codes.

**Root Cause:** `RichHandler` in `cli/main.py:365` had `markup=False`, causing Rich markup to be written as literal strings instead of being converted to ANSI codes.

**Fix:** Changed `markup=False` to `markup=True` in `cli/main.py:365`:
```python
file_handler = RichHandler(
    console=log_console,
    show_time=True,
    show_level=True,
    show_path=False,
    markup=True,  # Enable Rich markup interpretation for proper ANSI formatting
    rich_tracebacks=True,
    tracebacks_show_locals=config.output.debug,
    level=logging.DEBUG
)
```

**Files Modified:**
- `cli/main.py:365`

**Verification:** ANSI codes now appear correctly in log files (e.g., `[1;34m` for bold blue instead of literal `[bold blue]`).

---

### 3. Confusing Optimization Round Logging
**Issue:** Optimization round logging was confusing with redundant messages:
```
ROUND 2 SUMMARY
...
  Range too narrow (1.05×), expanding to 2× around best_C
Optimization round 3/5...

ROUND 3/5 - Grid Search
```

**Root Cause:**
1. Line 349 in `optimizer.py` printed "Optimization round X/Y..." before `_grid_search_round()` also printed its own header
2. Range expansion messages didn't clearly indicate they were preparing for the next round

**Fix:**
1. Removed redundant line 349 message
2. Added "[Next round preparation]" prefix to range expansion message (line 773):
```python
if self.verbose:
    print(f"\n[Next round preparation] Range too narrow ({current_ratio:.2f}×), expanding to {min_ratio:.0f}× around best_C", flush=True)
```

**Files Modified:**
- `classification/optimizer.py:349` (removed)
- `classification/optimizer.py:773` (added prefix)

**Expected Result:**
```
ROUND 2 SUMMARY
...

[Next round preparation] Range too narrow (1.05×), expanding to 2× around best_C

ROUND 3/5 - Grid Search
```

---

### 4. Progress Bars Not Updating During Training
**Issue:** Progress bars showed 0% and didn't update during grid search operations, appearing to "jump" from one step to the next.

**Root Cause:**
1. `leave=False` parameter caused progress bars to be erased when complete, making it look like they were jumping
2. Default `mininterval` (1.0s) meant infrequent refreshes
3. No fixed width caused inconsistent display

**Fix:** Updated `classification/optimizer.py:635-642` with better progress bar settings:
```python
with tqdm_joblib(tqdm(
    total=total_tasks,
    desc=desc,
    unit="fit",
    leave=True,  # Keep bar visible after completion
    mininterval=0.1,  # Refresh every 0.1s for smoother updates
    ncols=100  # Fixed width for consistent display
)):
    grid_search.fit(X, y)
```

**Files Modified:**
- `classification/optimizer.py:635-642`

**Expected Result:**
- Progress bars persist after completion (stack vertically)
- More frequent updates during long operations (every 0.1s)
- Consistent bar width across rounds

**Note:** With high parallelism (`-p 10`), joblib batches work such that progress updates occur in chunks rather than smoothly. This is expected behavior with parallel processing.

**Additional improvements for high parallelism (2025-01-17):**
1. Added `file=sys.stdout` to force tqdm output to stdout (prevent logging system capture)
2. Added `miniters=1` for more granular update thresholds
3. Reduced `mininterval` from 0.1s to 0.05s for faster refresh
4. Added background monitoring thread that prints time-based updates every 30s for long tasks (>50 fits)
5. Added informational note in output when using high parallelism explaining expected behavior

With these changes, users get feedback even if the tqdm bar doesn't update frequently:
```
Note: With high parallelism (n_jobs=10), progress bar may update
      in large increments as worker batches complete.
[Still running... 30s elapsed]
[Still running... 60s elapsed]
```

---

## Summary of Changes

**Files Modified:**
1. `cli/config.py` - Fixed train/classify mode detection
2. `cli/main.py` - Enabled Rich markup in log files
3. `classification/optimizer.py` - Improved logging clarity and progress bar visibility

**Testing:**
- Train subcommand now works correctly
- Log files show proper ANSI formatting
- Optimization round logging is clearer
- Progress bars persist and update more frequently

**Impact:**
These changes significantly improve the user experience during training by:
- Making the CLI more intuitive (train subcommand works as expected)
- Improving log file readability (proper formatting)
- Clarifying training progress (better round organization)
- Providing better feedback during long-running operations (visible progress bars)
