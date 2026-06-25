#!/bin/bash
# == Information ==============================================================
#  Create bigWigs files from strand BAMa


# == Configuration ============================================================
set -euo pipefail

SAMPLE=$1
THREADS=$2
BASE="/home/josvarcas/Archivos/master/c0_TFMEXP/TFM_DRIPc-seq"
CLEANDIR="${BASE}/data/bam_cleaned"
STRANDDIR="${BASE}/data/bam_stranded"
OUTDIR="${BASE}/data/bigwig"

mkdir -p ${OUTDIR}

FWD_BAM="${STRANDDIR}/${SAMPLE}.fwd.bam"
REV_BAM="${STRANDDIR}/${SAMPLE}.rev.bam"
FULL_BAM="${CLEANDIR}/${SAMPLE}.clean.bam"


# == bigWigs ==================================================================
echo "=== bigWig: ${SAMPLE} ==="
echo "  Start: $(date)"

# Calculate scale factor for each strand
TOTAL=$(samtools idxstats ${FULL_BAM} \
  | grep -v "^MT" | awk '{sum += $3} END {print sum}')
SCALE=$(LC_ALL=C awk "BEGIN {printf \"%.10f\", 1000000 / ${TOTAL}}")

echo "  Total reads (w/o MT):   ${TOTAL}"
echo "  ScaleFactor:            ${SCALE}"

# Create bigWig for forward strand
bamCoverage \
    -b ${FWD_BAM} \
    -o "${OUTDIR}/${SAMPLE}.fwd.bw" \
    --scaleFactor ${SCALE} \
    --normalizeUsing None \
    --binSize 10 \
    --extendReads \
    --numberOfProcessors ${THREADS}

# Create bigWig for reverse strand
bamCoverage \
    -b ${REV_BAM} \
    -o "${OUTDIR}/${SAMPLE}.rev.bw" \
    --scaleFactor ${SCALE} \
    --normalizeUsing None \
    --binSize 10 \
    --extendReads \
    --numberOfProcessors ${THREADS}

echo "  End: $(date)"
echo "✓ ${SAMPLE} → ${SAMPLE}.fwd.bw + ${SAMPLE}.rev.bw"
echo ""

