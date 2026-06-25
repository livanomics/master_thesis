#!/bin/bash
# == Information ==============================================================
#  Call peaks using macs3


# == Configuration ============================================================
set -euo pipefail

SAMPLE=$1
STRAND=$2
BASE="/home/josvarcas/Archivos/master/c0_TFMEXP/TFM_DRIPc-seq"
STRANDDIR="${BASE}/data/bam_stranded"
OUTDIR="${BASE}/data/peaks"

mkdir -p ${OUTDIR}

TREAT="${STRANDDIR}/${SAMPLE}.${STRAND}.bam"
CONTROL="${STRANDDIR}/RNH_pooled.${STRAND}.bam"
NAME="${SAMPLE}_${STRAND}"


# == Peak calling =============================================================
echo "=== MACS3 (broad + control RNH): ${NAME} ==="
echo "  Start: $(date)"

# Call peaks
macs3 callpeak \
    -t ${TREAT} \
    -c ${CONTROL} \
    -f BAMPE \
    -g hs \
    -n ${NAME} \
    --outdir ${OUTDIR} \
    --keep-dup all \
    --broad \
    -q 0.05 \
    --broad-cutoff 0.1 \
    2>> "${BASE}/bash/logs/macs3.log"

echo "  End: $(date)"
echo "✓ ${NAME} peak calling finished"


# == Quick stats ==============================================================
# Count number of peaks
NPEAKS=$(wc -l < ${OUTDIR}/${NAME}_peaks.broadPeak)
echo "✓ ${NAME}: ${NPEAKS} peaks minus RNAse H"
echo ""

