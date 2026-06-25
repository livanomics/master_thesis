# == 01. Information ==========================================================
#' Perform feature analysis (GC content, gene size, exon number, biotype) of
#' the annotated gene lists and the normality / Wilcoxon tests

# Load source files
source("src/config.R")
source("src/features_functions.R")

# Define paths for the script
DIFDIR  <- file.path(results, "differential")
OUTDIR  <- file.path(results, "features")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)


# == 02. Libraries ============================================================
library(biomaRt)
library(GenomicFeatures)


# == 03. Master Feature Table =================================================
# Load differential analysis and get the full list of symbols without NA
gene_lists <- readRDS(file.path(DIFDIR, "annotation_analysis.rds"))
all_genes  <- unique(unlist(gene_lists))
all_genes  <- all_genes[!is.na(all_genes)]

# Translate SYMBOL to ENSEMBL and ENTREZID and filter by ENSEMBL
map <- AnnotationDbi::select(org.Hs.eg.db,
                             keys    = all_genes,
                             keytype = "SYMBOL",
                             columns = c("ENSEMBL", "ENTREZID", "SYMBOL"))
map <- map[!is.na(map$ENSEMBL), ]
ensembl_ids <- unique(map$ENSEMBL)

# Connect to the Ensembl server and download attributes for each gene
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
bm <- getBM(attributes = c("ensembl_gene_id",
                           "hgnc_symbol",
                           "percentage_gene_gc_content",
                           "start_position",
                           "end_position",
                           "gene_biotype",
                           "chromosome_name"),
            filters = "ensembl_gene_id",
            values  = ensembl_ids,
            mart    = mart)

# Calculate gene length
bm$gene_length <- bm$end_position - bm$start_position + 1

# Number of exons per gene (overlapping exons reduced first)
exons_by_gene <- exonsBy(txdb, by = "gene")
exons_n <- sapply(exons_by_gene, function(g) length(reduce(g)))
exon_df <- data.frame(ENTREZID = names(exons_n),
                      exons_n  = as.integer(exons_n),
                      stringsAsFactors = FALSE)

# Combine all the information into a single table
master <- merge(map, bm,
                by.x = "ENSEMBL", by.y = "ensembl_gene_id", all.x = TRUE)
master <- merge(master, exon_df, by = "ENTREZID", all.x = TRUE)


# == 04. Per-List Feature Tables and Summary ==================================
features_summary <- data.frame()
for (l in names(gene_lists)) {
  # Filter the master table by the symbols of the list, remove duplicates
  sub <- master[master$SYMBOL %in% gene_lists[[l]], ]
  sub <- sub[!duplicated(sub$ENSEMBL), ]
  
  # Write per-list results
  write.csv(sub, file.path(OUTDIR, paste0("features_", l, ".csv")),
            row.names = FALSE)
  
  # Add a row to the summary table
  features_summary <- rbind(features_summary, data.frame(
    list           = l,
    gene_n         = nrow(sub),
    mean_GC        = round(mean(sub$percentage_gene_gc_content,
                                na.rm = TRUE), 2),
    mean_size_kb   = round(mean(sub$gene_length, na.rm = TRUE) / 1000, 1),
    median_size_kb = round(median(sub$gene_length, na.rm = TRUE) / 1000, 1),
    mean_exons     = round(mean(sub$exons_n, na.rm = TRUE), 1),
    pct_protcoding = round(mean(sub$gene_biotype == "protein_coding",
                                na.rm = TRUE) * 100, 1)))
}

write.csv(features_summary, file.path(OUTDIR, "features_summary.csv"),
          row.names = FALSE)


# == 05. Complete Long Table (sample / direction / orientation) ===============
# Parse the list name: last part = orientation, first = sample,
# middle = direction (NA for siC)
complete <- data.frame()
for (l in names(gene_lists)) {
  df <- read.csv(file.path(OUTDIR, paste0("features_", l, ".csv")),
                 stringsAsFactors = FALSE)
  if (nrow(df) == 0) next
  
  parts <- strsplit(l, "_")[[1]]
  df$sample      <- parts[1]
  df$direction   <- if (length(parts) == 3) parts[2] else NA
  df$orientation <- parts[length(parts)]
  
  complete <- rbind(complete,
                    df[, c("sample", "direction", "orientation",
                           "percentage_gene_gc_content",
                           "gene_length", "exons_n")])
}

# Derived columns: gene length in kb and group label (siC for control)
complete$length_kb <- complete$gene_length / 1000
complete$group <- ifelse(complete$sample == "siC", "siC",
                         paste(complete$sample, complete$direction, sep = "_"))
complete$group <- factor(complete$group, levels = group_levels)

# Save the complete table for the figures script
saveRDS(complete, file.path(OUTDIR, "features_complete.rds"))


# == 06. Statistical Analysis (per Orientation) ===============================
# Variables evaluated in every test
vars <- c(percentage_gene_gc_content = "GC content (%)",
          length_kb                  = "Gene size (kb)",
          exons_n                    = "Number of exons")

# Split data by orientation
sense_data     <- complete[complete$orientation == "sense", ]
antisense_data <- complete[complete$orientation == "antisense", ]

# Normality test (Shapiro-Wilk)
write.csv(run_normality(sense_data, vars),
          file.path(OUTDIR, "shapiro_sense.csv"), row.names = FALSE)
write.csv(run_normality(antisense_data, vars),
          file.path(OUTDIR, "shapiro_antisense.csv"), row.names = FALSE)

# Wilcoxon test (UP vs DOWN per condition, and each list vs siC)
write.csv(run_wilcoxon(sense_data, vars),
          file.path(OUTDIR, "wilcoxon_sense.csv"), row.names = FALSE)
write.csv(run_wilcoxon(antisense_data, vars),
          file.path(OUTDIR, "wilcoxon_antisense.csv"), row.names = FALSE)
