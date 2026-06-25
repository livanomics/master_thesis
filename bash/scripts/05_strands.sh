#!/bin/bash
# == Information ==============================================================
#  Separate BAMs into strand specific BAMs


# == Configuration ============================================================
set -euo pipefail

SAMPLE=$1
THREADS=$2
BASE="/home/josvarcas/Archivos/master/c0_TFMEXP/TFM_DRIPc-seq"
BAMDIR="${BASE}/data/bam_cleaned"
OUTDIR="${BASE}/data/bam_stranded"

mkdir -p ${OUTDIR}

IN="${BAMDIR}/${SAMPLE}.clean.bam"


# == Separation ===============================================================
echo "=== Separating strands: ${SAMPLE} ==="
echo "  Start: $(date)"

# Create forward strand from R1-reverse and R2-forward
samtools view -@ ${THREADS} -b -f 80 ${IN} \
  -o "${OUTDIR}/${SAMPLE}.fwd1.bam"

samtools view -@ ${THREADS} -b -f 128 -F 16 ${IN} \
  -o "${OUTDIR}/${SAMPLE}.fwd2.bam"

# Merge, re-sort and index forward strand
samtools merge -@ ${THREADS} -f -u - \
    "${OUTDIR}/${SAMPLE}.fwd1.bam" \
    "${OUTDIR}/${SAMPLE}.fwd2.bam" \
    | samtools sort -@ ${THREADS} -o "${OUTDIR}/${SAMPLE}.fwd.bam" -

samtools index "${OUTDIR}/${SAMPLE}.fwd.bam"

# Create reverse strand from R1-forward and R2-reverse
samtools view -@ ${THREADS} -b -f 64 -F 16 ${IN} \
  -o "${OUTDIR}/${SAMPLE}.rev1.bam"

samtools view -@ ${THREADS} -b -f 144 ${IN} \
  -o "${OUTDIR}/${SAMPLE}.rev2.bam"

# Merge, re-sort and index reverse strand
samtools merge -@ ${THREADS} -f -u - \
    "${OUTDIR}/${SAMPLE}.rev1.bam" \
    "${OUTDIR}/${SAMPLE}.rev2.bam" \
    | samtools sort -@ ${THREADS} -o "${OUTDIR}/${SAMPLE}.rev.bam" -

samtools index "${OUTDIR}/${SAMPLE}.rev.bam"

# Delete intermediary files
rm "${OUTDIR}/${SAMPLE}.fwd1.bam" "${OUTDIR}/${SAMPLE}.fwd2.bam"
rm "${OUTDIR}/${SAMPLE}.rev1.bam" "${OUTDIR}/${SAMPLE}.rev2.bam"

echo "  End: $(date)"
echo "✓ ${SAMPLE} separation finished"


# == Verification =============================================================
# Count reads in cleaned BAM
TOTAL=$(samtools view -c ${IN})

# Count reads in forward BAM
FWD=$(samtools view -c ${OUTDIR}/${SAMPLE}.fwd.bam)

# Count reads in reverse BAM
REV=$(samtools view -c ${OUTDIR}/${SAMPLE}.rev.bam)

# Sum reads from strands
SUM=$((FWD + REV))

# Show stats
echo ""
echo "  Total: ${TOTAL} | Fwd: ${FWD} | Rev: ${REV} | Suma: ${SUM}"
echo ""

