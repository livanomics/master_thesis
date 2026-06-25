#!/bin/bash
# == Information ==============================================================
#  Flag the duplicate reads in the .bam files


# == Configuration ============================================================
set -euo pipefail

SAMPLE=$1
THREADS=$2
BASE="/home/josvarcas/Archivos/master/c0_TFMEXP/TFM_DRIPc-seq"

cd "${BASE}/data/bam"


# == Flagging =================================================================
echo "=== samtools markdup: ${SAMPLE} ==="
echo "  Start: $(date)"

# Sort .bam file by name
samtools sort -@ ${THREADS} -n -o ${SAMPLE}.nsort.bam ${SAMPLE}.bam

# Invoque fixmate
samtools fixmate -@ ${THREADS} -m ${SAMPLE}.nsort.bam ${SAMPLE}.fixmate.bam

# Sort fixmate by coordinates
samtools sort -@ ${THREADS} -o ${SAMPLE}.csort.bam ${SAMPLE}.fixmate.bam

# Mark duplicates in .bam
samtools markdup -@ ${THREADS} -s \
    ${SAMPLE}.csort.bam \
    ${SAMPLE}.markdup.bam \
    2>> "${BASE}/bash/logs/bam_markdup_stats.log"

# Index .bam that has duplicates flagged
samtools index ${SAMPLE}.markdup.bam

# Delete intermediary files
rm ${SAMPLE}.nsort.bam ${SAMPLE}.fixmate.bam ${SAMPLE}.csort.bam

echo "  End: $(date)"
echo "✓ ${SAMPLE} markdup finished"


# == Quick stats ==============================================================
# Count the reads
TOTAL=$(samtools view -c ${SAMPLE}.markdup.bam)

# Count reads marked as duplicates
DUP=$(samtools view -c -f 1024 ${SAMPLE}.markdup.bam)

# Calculate percentage of duplicates
PCT=$(awk "BEGIN {printf \"%.2f\", ${DUP}*100/${TOTAL}}")

# Show stats
echo "  Total reads: ${TOTAL}"
echo "  Duplicates:    ${DUP} (${PCT}%)"
echo ""

