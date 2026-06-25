#!/bin/bash
# == Information ==============================================================
#  Perform QC on the cleaned .bam files


# == Configuration ============================================================
set -euo pipefail

SAMPLE=$1
THREADS=$2
BASE="/home/josvarcas/Archivos/master/c0_TFMEXP/TFM_DRIPc-seq"
BAMDIR="${BASE}/data/bam_cleaned"
QCDIR="${BASE}/results/qc/bam_cleaned_qc"

mkdir -p ${QCDIR}


# == QC analysis ==============================================================
# Creates paths to .bam and labels in the same order
BAMFILES=""
LABELS=""
for S in "${SAMPLE[@]}"; do
    BAMFILES="${BAMFILES} ${BAMDIR}/${S}.clean.bam"
    LABELS="${LABELS} ${S}"
done

echo "=== Cleaned BAMs QC ==="
echo "  Start: $(date)"

# Show stats per sample
echo ""
echo "=== [1/4] flagstat ==="

for S in "${SAMPLE[@]}"; do
    samtools flagstat "${BAMDIR}/${S}.clean.bam" > "${QCDIR}/${S}.flagstat.txt"
    echo "  ✓ ${S}"
done

# Create multiBamSummary with 10kb bins
echo ""
echo "=== [2/4] multiBamSummary ==="

multiBamSummary bins \
    --bamfiles ${BAMFILES} \
    --labels ${LABELS} \
    --binSize 10000 \
    --numberOfProcessors ${THREADS} \
    --outRawCounts "${QCDIR}/counts.tab" \
    -o "${QCDIR}/multibam.npz"

echo "  ✓ Matrix created: ${QCDIR}/multibam.npz"

# Create a plotCorrelation heatmap with Spearman coefficient
echo ""
echo "=== [3/4] plotCorrelation ==="

plotCorrelation \
    -in "${QCDIR}/multibam.npz" \
    --corMethod spearman \
    --whatToPlot heatmap \
    --skipZeros \
    --plotNumbers \
    --colorMap RdYlBu_r \
    --plotTitle "Spearman correlation - DRIPc-seq (10kb bins)" \
    -o "${QCDIR}/correlation_heatmap.png" \
    --outFileCorMatrix "${QCDIR}/correlation_matrix.tab"

echo "  ✓ Heatmap: ${QCDIR}/correlation_heatmap.png"

# Create a PCA of the data
echo ""
echo "=== [4/4] plotPCA ==="

plotPCA \
    -in "${QCDIR}/multibam.npz" \
    --labels ${LABELS} \
    --plotTitle "PCA - DRIPc-seq (10kb bins)" \
    --plotFileFormat png \
    -o "${QCDIR}/pca.png" \
    --outFileNameData "${QCDIR}/pca_data.tab"

echo "  ✓ PCA: ${QCDIR}/pca.png"

echo ""
echo "  End: $(date)"
echo "✓ Cleaned BAM QC finished"

