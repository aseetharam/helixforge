#!/bin/bash
set -euo pipefail

# HelixForge chunk processing wrapper
# Usage: bash chunk_wrapper.sh <chunk_id>
#
# This script processes a single chunk from a chunk plan.
# It reads chunk details from the CHUNK_PLAN environment variable
# or from a default chunks.json file.
#
# Environment variables:
#   CHUNK_PLAN - Path to chunk plan JSON (default: chunks.json)
#   OUTPUT_DIR - Output directory (optional)
#
# The script exports these variables for use in commands:
#   $CHUNK_ID  - Chunk identifier
#   $SEQID     - Scaffold/chromosome name
#   $START     - Start coordinate (0-based)
#   $END       - End coordinate (0-based, exclusive)
#   $SIZE      - Chunk size in bases

# =============================================================================
# Configuration
# =============================================================================

CHUNK_ID="${1:?Chunk ID required}"
CHUNK_PLAN="${CHUNK_PLAN:-chunks.json}"
OUTPUT_DIR="${OUTPUT_DIR:-outputs}"

# =============================================================================
# Parse Chunk Information
# =============================================================================

if [ ! -f "$CHUNK_PLAN" ]; then
    echo "ERROR: Chunk plan not found: $CHUNK_PLAN" >&2
    exit 1
fi

# Parse chunk info from JSON using Python
read SEQID START END < <(python3 -c "
import json
import sys

with open('${CHUNK_PLAN}') as f:
    plan = json.load(f)

for chunk in plan['chunks']:
    if chunk['chunk_id'] == '${CHUNK_ID}':
        print(chunk['seqid'], chunk['start'], chunk['end'])
        sys.exit(0)

print('ERROR: Chunk ${CHUNK_ID} not found in plan', file=sys.stderr)
sys.exit(1)
")

SIZE=$((END - START))

# Export for use in subprocesses
export CHUNK_ID SEQID START END SIZE CHUNK_PLAN OUTPUT_DIR

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "HelixForge Chunk Processing"
echo "Chunk ID: ${CHUNK_ID}"
echo "Region: ${SEQID}:${START}-${END}"
echo "Size: ${SIZE} bp"
echo "Output: ${OUTPUT_DIR}"
echo "Start time: $(date)"
echo "=========================================="

# =============================================================================
# Setup Commands (customize for your environment)
# =============================================================================

# Uncomment and modify as needed for your HPC environment:

# module load python/3.10
# module load conda
# source $(conda info --base)/etc/profile.d/conda.sh
# conda activate helixforge

# =============================================================================
# Main Commands (customize for your workflow)
# =============================================================================

# Example: Confidence scoring
# helixforge confidence \
#     --helixer-h5 predictions.h5 \
#     --helixer-gff predictions.gff3 \
#     --genome genome.fa \
#     --chunk-id "${CHUNK_ID}" \
#     --region "${SEQID}:${START}-${END}" \
#     --output "${OUTPUT_DIR}/${CHUNK_ID}_confidence.tsv"

# Example: Splice refinement
# helixforge splice \
#     --helixer-gff predictions.gff3 \
#     --genome genome.fa \
#     --rnaseq-bam rnaseq.bam \
#     --region "${SEQID}:${START}-${END}" \
#     --output-gff "${OUTPUT_DIR}/${CHUNK_ID}.gff3" \
#     --report "${OUTPUT_DIR}/${CHUNK_ID}_splice.tsv"

# Placeholder - replace with your commands
echo "Processing region ${SEQID}:${START}-${END}..."
echo "Add your helixforge commands here"

# =============================================================================
# Completion
# =============================================================================

# Mark chunk as complete
touch "${OUTPUT_DIR}/${CHUNK_ID}.done"

echo "=========================================="
echo "Chunk ${CHUNK_ID} complete"
echo "End time: $(date)"
echo "=========================================="
