# == 01. Information ==========================================================
#' Perform differential analysis of the peaks in each condition and annotate
#' the coordinates to genes and the orientation (sense/antisense) of the
#' R-loops within them

# Load source files
source("src/config.R")
source("src/differential_annotation_functions.R")

# Define paths for the script
BAMDIR  <- file.path(data, "bam_stranded")
BEDDIR  <- file.path(data, "bed")
PEAKDIR <- file.path(data, "peaks")
OUTDIR  <- file.path(results, "differential")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)


# == 02. Libraries ============================================================
library(DiffBind)
library(ChIPseeker)


# == 03. Differential Analysis ================================================
# Run differential analysis of peaks in different samples and strands
db_uap56_fwd <- run_diffbind(c("siUAP56_1", "siUAP56_2"), "siUAP56", "fwd")
db_uap56_rev <- run_diffbind(c("siUAP56_1", "siUAP56_2"), "siUAP56", "rev")
db_brg1_fwd  <- run_diffbind(c("siBRG1_1", "siBRG1_2"), "siBRG1", "fwd")
db_brg1_rev  <- run_diffbind(c("siBRG1_1", "siBRG1_2"), "siBRG1", "rev")


# == 04. Gene Annotation ======================================================
# Define lists and data.frame that store the data
gene_lists    <- list()
orient_anno   <- list()
sense_summary <- data.frame()

# Define lists of data to process
to_process <- list(
  list(cond = "siC",     source = "bed",  direction = NULL),
  list(cond = "siUAP56", source = "diff", direction = "UP"),
  list(cond = "siUAP56", source = "diff", direction = "DOWN"),
  list(cond = "siBRG1",  source = "diff", direction = "UP"),
  list(cond = "siBRG1",  source = "diff", direction = "DOWN")
)

# Loop that runs the combination of data and annotation
for (item in to_process) {
  res <- annotate_combined(item$cond, item$source, item$direction)
  
  # Create label for each run
  label <- if (is.null(item$direction)) item$cond
  else paste0(item$cond, "_", item$direction)
  
  orient_anno[[label]] <- res
  
  # Gene lists by orientation
  for (ori in c("sense", "antisense")) {
    genes <- genes_by_orientation(res$df, ori)
    gene_lists[[paste0(label, "_", ori)]] <- genes
    writeLines(genes, file.path(OUTDIR,
                                paste0("genes_", label, "_", ori, ".txt")))
  }
  
  # Save the full annotation table with orientation column
  write.csv(res$df, file.path(OUTDIR, paste0("anno_", label, ".csv")),
            row.names = FALSE)
  
  # Create a summary of each iteration
  sense_summary <- rbind(sense_summary,
                         data.frame(list = label, pct_sense = res$pct_sense))
  cat(sprintf("  [ %s ] %.1f%% sense  (%d sense / %d antisense genes)\n",
              label, res$pct_sense,
              length(gene_lists[[paste0(label, "_sense")]]),
              length(gene_lists[[paste0(label, "_antisense")]])))
}


# == 05. Save Results =========================================================
saveRDS(gene_lists,  file.path(OUTDIR, "annotation_analysis.rds"))
saveRDS(orient_anno, file.path(OUTDIR, "orientation_annotation.rds"))
write.csv(sense_summary, file.path(OUTDIR, "sense_summary.csv"),
          row.names = FALSE)
