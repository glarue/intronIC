"""
Debug script to check if BP and 3'SS raw scores are being calculated.
"""
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

from intronIC.cli.main import load_reference_sequences, score_introns
from config.config import IntronICConfig
from intronIC.cli.progress import IntronICProgressReporter
from intronIC.cli.messenger import UnifiedMessenger
import logging

# Setup minimal config
config = IntronICConfig()
config.output.quiet = False
config.performance.processes = 1

# Setup logging
logger = logging.getLogger('intronIC')
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter('%(message)s'))
logger.addHandler(handler)

# Setup messenger
reporter = IntronICProgressReporter(quiet=False)
log_console = reporter.console  # Dummy
messenger = UnifiedMessenger(
    console=reporter.console,
    log_console=log_console,
    logger=logger,
    quiet=False
)

# Load reference sequences (from src/intronIC/data/)
data_dir = Path(__file__).parent.parent / "src" / "intronIC" / "data"
u12_file = data_dir / "u12_reference.introns.iic.gz"
u2_file = data_dir / "u2_reference.introns.iic.gz"

print("Loading reference sequences...")
u12_ref = load_reference_sequences(u12_file, messenger=messenger)
u2_ref = load_reference_sequences(u2_file, messenger=messenger)

print(f"Loaded {len(u12_ref)} U12 and {len(u2_ref)} U2 references")

# Score a sample
print("\nScoring first 10 U12 introns...")
sample = u12_ref[:10]
scored = score_introns(sample, config, messenger, reporter)

# Check scores
print("\nChecking raw scores:")
print(f"{'Intron':30s} {'5_raw':>10s} {'BP_raw':>10s} {'3_raw':>10s}")
print("-" * 65)

for intron in scored[:10]:
    if intron.scores:
        print(f"{intron.intron_id[:30]:30s} "
              f"{intron.scores.five_raw_score:>10.4f} "
              f"{intron.scores.bp_raw_score:>10.4f} "
              f"{intron.scores.three_raw_score:>10.4f}")
    else:
        print(f"{intron.intron_id[:30]:30s} NO SCORES")

# Check if any are zero or None
zeros = sum(1 for i in scored if i.scores and (
    i.scores.bp_raw_score == 0 or i.scores.three_raw_score == 0
))
nones = sum(1 for i in scored if i.scores and (
    i.scores.bp_raw_score is None or i.scores.three_raw_score is None
))

print(f"\nStatistics:")
print(f"  Introns with BP or 3' score = 0: {zeros}/{len(scored)}")
print(f"  Introns with BP or 3' score = None: {nones}/{len(scored)}")

# Check variance
import numpy as np
bp_scores = [i.scores.bp_raw_score for i in scored if i.scores and i.scores.bp_raw_score is not None]
three_scores = [i.scores.three_raw_score for i in scored if i.scores and i.scores.three_raw_score is not None]

if bp_scores:
    print(f"\nBP raw scores: mean={np.mean(bp_scores):.4f}, std={np.std(bp_scores):.4f}, min={np.min(bp_scores):.4f}, max={np.max(bp_scores):.4f}")
if three_scores:
    print(f"3' raw scores: mean={np.mean(three_scores):.4f}, std={np.std(three_scores):.4f}, min={np.min(three_scores):.4f}, max={np.max(three_scores):.4f}")
