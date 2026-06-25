#!/bin/bash
# == Information ==============================================================
#  Map trimmed data to the reference genome to get .bam files


# == Configuration ============================================================
set -euo pipefail

SAMPLE=$1
THREADS=$2
BASE="/home/josvarcas/Archivos/master/c0_TFMEXP/TFM_DRIPc-seq"

cd "${BASE}/data/trimmed"
mkdir -p "${BASE}/data/bam"


# == Mapping ==================================================================
echo "=== bowtie2 mapping: ${SAMPLE} ==="
echo "  Start: $(date)"

# Map trimmed data to the genome outputting a .sam file
bowtie2 -x "${BASE}/data/genome/index_bowtie2/index" \
  -1 "${SAMPLE}_R1.trim.fastq.gz" \
  -2 "${SAMPLE}_R2.trim.fastq.gz" \
  --no-unal \
  -p ${THREADS} \
  -S "${BASE}/data/bam/${SAMPLE}.sam" \
  2>> "${BASE}/bash/logs/bowtie2.log"

# Sort and convert .sam into .bam file
samtools sort -@ ${THREADS} -O bam \
  -o "${BASE}/data/bam/${SAMPLE}.bam" \
  "${BASE}/data/bam/${SAMPLE}.sam"

# Delete .sam files to free up space
rm "${BASE}/data/bam/${SAMPLE}.sam"

echo "  End: $(date)"
echo "✓ ${SAMPLE} map finished"

