#!/bin/bash
# == Information ==============================================================
#  Merge RNAse H conditions into a single RNAse H input file by strand


# == Configuration ============================================================
set -euo pipefail

SAMPLE=$1
THREADS=$2
BASE="/home/josvarcas/Archivos/master/c0_TFMEXP/TFM_DRIPc-seq"
STRANDDIR="${BASE}/data/bam_stranded"


# == Merger ===================================================================
echo "=== RNAse H merger ==="
echo "  Start: $(date)"

if [[ ! -f "${STRANDDIR}/RNH_pooled.fwd.bam" ]]; then
    echo "=== Creating RNH pooled (fwd) ==="
    samtools merge -@ ${THREADS} -f \
        "${STRANDDIR}/RNH_pooled.fwd.bam" \
        "${STRANDDIR}/siC_RNH.fwd.bam" \
        "${STRANDDIR}/siUAP56_RNH.fwd.bam"

    samtools index "${STRANDDIR}/RNH_pooled.fwd.bam"

    echo "=== Creating RNH pooled (rev) ==="
    samtools merge -@ ${THREADS} -f \
        "${STRANDDIR}/RNH_pooled.rev.bam" \
        "${STRANDDIR}/siC_RNH.rev.bam" \
        "${STRANDDIR}/siUAP56_RNH.rev.bam"

    samtools index "${STRANDDIR}/RNH_pooled.rev.bam"

    echo "✓ RNH pooled created"
    echo "  fwd: $(samtools view -c ${STRANDDIR}/RNH_pooled.fwd.bam) reads"
    echo "  rev: $(samtools view -c ${STRANDDIR}/RNH_pooled.rev.bam) reads"
fi

echo "   End: $(date)"
echo "✓ RNAse H merger finished"

