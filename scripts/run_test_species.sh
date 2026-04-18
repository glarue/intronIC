#!/usr/bin/env bash
#
# Run intronIC classify on standard test species with a given model/config.
#
# Usage:
#   bash scripts/run_test_species.sh <model.pkl> <config.yaml> <output_dir> [processes]
#
# Example:
#   bash scripts/run_test_species.sh \
#     models/linear_3d_rebuilt_pwms/homo_sapiens.model.pkl \
#     config/config_linear_3d.yaml \
#     models/test_runs_linear_3d_rebuilt \
#     6
#
# After completion, summarize with:
#   python scripts/summarize_run.py <output_dir>

set -euo pipefail

MODEL="${1:?Usage: $0 <model.pkl> <config.yaml> <output_dir> [processes]}"
CONFIG="${2:?Usage: $0 <model.pkl> <config.yaml> <output_dir> [processes]}"
OUTDIR="${3:?Usage: $0 <model.pkl> <config.yaml> <output_dir> [processes]}"
PROCS="${4:-6}"

INTRONIC="${INTRONIC:-/mnt/data/u12/intronIC_v2_env/bin/intronIC}"

# ── Species data paths ──
# Format: name|genome|annotation
SPECIES=(
"homo_sapiens|/mnt/data/u12/refseq/vertebrates/homo_sapiens-GRCh38.p13/homo_sapiens.GCF_000001405.39_GRCh38.p13_genomic.fna.gz|/mnt/data/u12/refseq/vertebrates/homo_sapiens-GRCh38.p13/homo_sapiens.GCF_000001405.39_GRCh38.p13_genomic.gff.gz"
"drosophila_melanogaster|/mnt/data/u12/refseq/invertebrates/drosophila_melanogaster-Release_6_plus_ISO1_MT/drosophila_melanogaster.GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz|/mnt/data/u12/refseq/invertebrates/drosophila_melanogaster-Release_6_plus_ISO1_MT/drosophila_melanogaster.GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff.gz"
"arabidopsis_thaliana|/mnt/data/u12/refseq/plants/arabidopsis_thaliana-TAIR10.1/arabidopsis_thaliana.GCF_000001735.4_TAIR10.1_genomic.fna.gz|/mnt/data/u12/refseq/plants/arabidopsis_thaliana-TAIR10.1/arabidopsis_thaliana.GCF_000001735.4_TAIR10.1_genomic.gff.gz"
"caenorhabditis_elegans|/mnt/data/u12/refseq/invertebrates/caenorhabditis_elegans-WBcel235/caenorhabditis_elegans.GCF_000002985.6_WBcel235_genomic.fna.gz|/mnt/data/u12/refseq/invertebrates/caenorhabditis_elegans-WBcel235/caenorhabditis_elegans.GCF_000002985.6_WBcel235_genomic.gff.gz"
"ascaris_suum|/mnt/data/u12/new_species/ascaris_suum__pig_roundworm__ag01_-ASM18702v3/Ascaris_suum.ASM18702v3.dna.toplevel.fa.gz|/mnt/data/u12/new_species/ascaris_suum__pig_roundworm__ag01_-ASM18702v3/Ascaris_suum.ASM18702v3.62.gtf.gz"
"symbiodinium_microadriaticum|/mnt/data/u12/genbank/protozoa/symbiodinium_microadriaticum-ASM193914v1/symbiodinium_microadriaticum.GCA_001939145.1_ASM193914v1_genomic.fna.gz|/mnt/data/u12/genbank/protozoa/symbiodinium_microadriaticum-ASM193914v1/symbiodinium_microadriaticum.GCA_001939145.1_ASM193914v1_genomic.gff.gz"
)

# ── Validate inputs ──
if [ ! -f "$MODEL" ]; then
    echo "Error: model not found: $MODEL" >&2
    exit 1
fi
if [ ! -f "$CONFIG" ]; then
    echo "Error: config not found: $CONFIG" >&2
    exit 1
fi

mkdir -p "$OUTDIR"

echo "Model:  $MODEL"
echo "Config: $CONFIG"
echo "Output: $OUTDIR"
echo "Procs:  $PROCS"
echo ""

# ── Run all species in parallel ──
PIDS=()
for entry in "${SPECIES[@]}"; do
    IFS='|' read -r name genome annot <<< "$entry"

    if [ ! -f "$genome" ]; then
        echo "Warning: genome not found for $name, skipping" >&2
        continue
    fi
    if [ ! -f "$annot" ]; then
        echo "Warning: annotation not found for $name, skipping" >&2
        continue
    fi

    echo "Starting: $name"
    "$INTRONIC" --config "$CONFIG" classify \
        -g "$genome" -a "$annot" -n "$name" \
        --model "$MODEL" -o "$OUTDIR" -p "$PROCS" \
        > /dev/null 2>&1 &
    PIDS+=($!)
done

echo ""
echo "Waiting for ${#PIDS[@]} species to finish..."

# Wait for all and report
FAILED=0
for pid in "${PIDS[@]}"; do
    if ! wait "$pid"; then
        FAILED=$((FAILED + 1))
    fi
done

echo ""
if [ "$FAILED" -gt 0 ]; then
    echo "Done ($FAILED failures). Check .iic.log files for details."
else
    echo "Done (all succeeded)."
fi

echo ""
echo "Summary:"
python3 "$(dirname "$0")/summarize_run.py" "$OUTDIR"
