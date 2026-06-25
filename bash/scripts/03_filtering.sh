#!/bin/bash
# == Information ==============================================================
#  Filter .bam files to have cleaned ones for other analysis


# == Configuration ============================================================
set -euo pipefail

SAMPLE=$1
THREADS=$2
BASE="/home/josvarcas/Archivos/master/c0_TFMEXP/TFM_DRIPc-seq"
BLACKLIST="${BASE}/data/genome/hg38-blacklist.v2.ensembl.bed"

cd "${BASE}/data/bam"
mkdir -p "${BASE}/data/bam_cleaned"


# == Filtering ================================================================
echo "=== Filtering: ${SAMPLE} ==="
echo "  Start: $(date)"

# Filter flagged, MAPQ and names of chromosomes
samtools view -@ ${THREADS} -b \
    -f 2 -F 1804 -q 20 \
    "${SAMPLE}.markdup.bam" \
    1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT \
    -o "${SAMPLE}.filt.bam"

# Filter blacklist
bedtools intersect -v \
    -abam "${SAMPLE}.filt.bam" \
    -b ${BLACKLIST} \
    > "${BASE}/data/bam_cleaned/${SAMPLE}.clean.bam"

# Index .bam that is cleaned
samtools index "${BASE}/data/bam_cleaned/${SAMPLE}.clean.bam"

# Delete intermediary files
rm "${SAMPLE}.filt.bam"

echo "  End: $(date)"
echo "✓ ${SAMPLE} filtering finished"


# == Quick stats ==============================================================
# Count reads before filtering
BEFORE=$(samtools view -c ${SAMPLE}.markdup.bam)

# Count reads after filtering
AFTER=$(samtools view -c ${BASE}/bam_cleaned/${SAMPLE}.clean.bam)

# Calculate percentage of saved reads
PCT=$(awk "BEGIN {printf \"%.2f\", ${AFTER}*100/${BEFORE}}")

# Counts reads in MT chromosome
CHRM=$(samtools idxstats ${BASE}/bam_cleaned/${SAMPLE}.clean.bam | grep -P "^MT\t" | cut -f3 || echo "0")

# Show stats
echo "${SAMPLE}: ${BEFORE} → ${AFTER} (${PCT}%)  | reads MT: ${CHRM}"
echo ""

