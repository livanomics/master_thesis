# == 01. Information ==========================================================
#' Physical mapping of the R-loops: fraction of affected genes per chromosome,
#' correlation with chromosomal GC content, long-gene feature distribution and
#' the intronic profile of siUAP56_UP

# Load source files
source("src/config.R")
source("src/differential_annotation_functions.R")
source("src/mapping_functions.R")

# Define paths for the script
DIFDIR <- file.path(results, "differential")
OUTDIR <- file.path(results, "physical_mapping")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)


# == 02. Libraries ============================================================
library(ChIPseeker)
library(GenomicFeatures)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

genome <- BSgenome.Hsapiens.UCSC.hg38
canon  <- paste0("chr", c(1:22, "X", "Y"))


# == 03. Gene Lists per Condition + Direction =================================
gene_lists <- readRDS(file.path(DIFDIR, "annotation_analysis.rds"))

# Combine sense + antisense symbols, then translate to ENTREZID (for the TxDb)
lists_sym <- list(
  siUAP56_UP   = genes_of(gene_lists, "siUAP56_UP"),
  siUAP56_DOWN = genes_of(gene_lists, "siUAP56_DOWN"),
  siBRG1_UP    = genes_of(gene_lists, "siBRG1_UP"),
  siBRG1_DOWN  = genes_of(gene_lists, "siBRG1_DOWN")
)
lists_entrez <- lapply(lists_sym, sym_to_entrez)


# == 04. Fraction of Genes per Chromosome =====================================
chrom_data <- build_chrom_tab(lists_entrez, canon)

# Split 'list' into condition and direction (used by the figures)
chrom_data$condition <- sub("_(UP|DOWN)$", "", chrom_data$list)
chrom_data$direction <- factor(sub("^.*_", "", chrom_data$list),
                               levels = c("UP", "DOWN"))
saveRDS(chrom_data, file.path(OUTDIR, "chrom_data.rds"))


# == 05. Correlation with Chromosomal GC Content ==============================
prop_chr  <- chrom_gc_content(genome, canon)
corr_data <- merge(chrom_data, prop_chr, by = "chrom")
saveRDS(corr_data, file.path(OUTDIR, "corr_data.rds"))

# Spearman correlation GC vs fraction, per list
correlations <- data.frame()
for (l in unique(corr_data$list)) {
  sub <- corr_data[corr_data$list == l, ]
  ct  <- cor.test(sub$GC, sub$fraction, method = "spearman")
  correlations <- rbind(correlations, data.frame(
    list = l, rho = round(ct$estimate, 3), p_value = signif(ct$p.value, 3)))
}
# FDR correction across the 4 correlations
correlations$p_adj       <- signif(p.adjust(correlations$p_value,
                                            method = "BH"), 3)
correlations$significant <- correlations$p_adj < 0.05

write.csv(correlations, file.path(OUTDIR, "GC_correlations.csv"),
          row.names = FALSE)


# == 06. Long-gene Feature Distribution =======================================
orient_anno <- readRDS(file.path(DIFDIR, "orientation_annotation.rds"))

# Long-gene threshold: 75th percentile of all gene lengths
long_threshold <- quantile(width(genes(txdb)), 0.75)

# Re-annotate only the regions inside long genes
anno_long <- list(
  siC        = long_genes_anno(orient_anno[["siC"]],        long_threshold),
  siUAP56_UP = long_genes_anno(orient_anno[["siUAP56_UP"]], long_threshold)
)
anno_long <- anno_long[!sapply(anno_long, is.null)]

saveRDS(anno_long, file.path(OUTDIR, "anno_long.rds"))


# == 07. siUAP56_UP: Intronic vs Non-intronic Gene Size =======================
df_int <- orient_anno[["siUAP56_UP"]]$df
df_int$is_intron <- grepl("Intron", df_int$annotation)
df_int$gene_kb   <- (df_int$geneEnd - df_int$geneStart) / 1000
df_int <- df_int[!is.na(df_int$gene_kb), ]

# Wilcoxon: gene size intronic vs non-intronic
w <- wilcox.test(gene_kb ~ is_intron, data = df_int)
test_res <- data.frame(
  median_intron    = round(median(df_int$gene_kb[df_int$is_intron]), 1),
  median_no_intron = round(median(df_int$gene_kb[!df_int$is_intron]), 1),
  n_intron         = sum(df_int$is_intron),
  n_no_intron      = sum(!df_int$is_intron),
  p_value          = signif(w$p.value, 3))

write.csv(test_res, file.path(OUTDIR, "siUAP56_UP_intron_size_test.csv"),
          row.names = FALSE)

# Save the table needed for the plot
saveRDS(df_int[, c("is_intron", "gene_kb")],
        file.path(OUTDIR, "siUAP56_UP_intron.rds"))
