# == 01. Information ==========================================================
#' Functions for the differential analysis of peaks (DiffBind) and the gene
#' annotation of the resulting regions (ChIPseeker), including the sense /
#' antisense assignment and helpers to extract gene lists


# == 02. Helpers ==============================================================
# Translate chromosome names from Ensembl to UCSC notation
to_ucsc <- function(seqnames) {
  seqnames <- ifelse(grepl("^chr", seqnames), seqnames,
                     paste0("chr", seqnames))
  gsub("^chrMT$", "chrM", seqnames)
}


# == 03. Differential Analysis ================================================
# Build the DiffBind sample sheet for one condition and strand
build_sheet <- function(treat_samples, treat_cond, strand) {
  # Basic structure of the experiment (controls + treatment)
  ctrl_samples <- c("siC_2", "siC_3")
  all_samples  <- c(ctrl_samples, treat_samples)
  all_conds    <- c(rep("siC", length(ctrl_samples)),
                    rep(treat_cond, length(treat_samples)))
  
  data.frame(
    SampleID   = paste0(all_samples, "_", strand),
    Tissue     = "K562",
    Factor     = "Rloop",
    Condition  = all_conds,
    Replicate  = c(seq_along(ctrl_samples), seq_along(treat_samples)),
    bamReads   = file.path(BAMDIR, paste0(all_samples, ".", strand, ".bam")),
    ControlID  = paste0("RNH_", strand),
    bamControl = file.path(BAMDIR, paste0("RNH_pooled.", strand, ".bam")),
    Peaks      = file.path(PEAKDIR, paste0(all_samples, "_",
                                           strand, "_peaks.broadPeak")),
    PeakCaller = "bed",
    stringsAsFactors = FALSE
  )
}

# Run the full DiffBind workflow and save the results (CSV + RDS)
run_diffbind <- function(treat_samples, treat_cond, strand) {
  sheet <- build_sheet(treat_samples, treat_cond, strand)
  
  # Load samples and count reads over the consensus regions
  dba_obj <- dba(sampleSheet = sheet)
  dba_obj <- dba.count(dba_obj, summits = FALSE)
  
  # Normalize samples (RLE method from DESeq2)
  dba_obj <- dba.normalize(dba_obj)
  
  # Define the statistical contrast and run the DESeq2 test
  dba_obj <- dba.contrast(dba_obj, categories = DBA_CONDITION, minMembers = 2)
  dba_obj <- dba.analyze(dba_obj, method = DBA_DESEQ2)
  
  # Extract and save the differential regions (CSV + RDS)
  db_results <- dba.report(dba_obj)
  cat("\nDifferential regions (FDR<0.05):", length(db_results), "\n")
  
  write.csv(as.data.frame(db_results),
            file.path(OUTDIR, paste0(treat_cond, "_vs_siC_", strand, ".csv")),
            row.names = FALSE)
  saveRDS(dba_obj,
          file.path(OUTDIR, paste0(treat_cond, "_vs_siC_", strand, ".rds")))
  
  dba_obj
}


# == 04. Region Loading =======================================================
# Read a BED file and build a GRanges object (UCSC notation)
load_bed <- function(bedfile) {
  df <- read.table(bedfile, sep = "\t", stringsAsFactors = FALSE)[, 1:3]
  colnames(df) <- c("seqnames", "start", "end")
  df$seqnames  <- to_ucsc(df$seqnames)
  
  GRanges(seqnames = df$seqnames,
          ranges   = IRanges(df$start, df$end))
}

# Read DiffBind results, filter by FDR / direction and build a GRanges object
load_diff <- function(csvfile, direction, fdr = 0.05) {
  df <- read.csv(csvfile, stringsAsFactors = FALSE)
  df <- df[df$FDR < fdr, ]
  
  if (direction == "UP")   df <- df[df$Fold > 0, ]
  if (direction == "DOWN") df <- df[df$Fold < 0, ]
  if (nrow(df) == 0) return(GRanges())
  
  df$seqnames <- to_ucsc(df$seqnames)
  
  GRanges(seqnames = df$seqnames,
          ranges   = IRanges(start = df$start, end = df$end),
          Fold     = df$Fold,
          FDR      = df$FDR)
}


# == 05. Gene Annotation ======================================================
# Combine fwd/rev regions, annotate to genes and assign sense / antisense
annotate_combined <- function(cond, source = "diff", direction = NULL) {
  # Load both strands according to the source
  if (source == "bed") {
    fwd <- load_bed(file.path(BEDDIR, paste0(cond, "_fwd.bed")))
    rev <- load_bed(file.path(BEDDIR, paste0(cond, "_rev.bed")))
  } else {
    fwd <- load_diff(file.path(OUTDIR, paste0(cond, "_vs_siC_fwd.csv")),
                     direction)
    rev <- load_diff(file.path(OUTDIR, paste0(cond, "_vs_siC_rev.csv")),
                     direction)
  }
  
  # Assign genomic strand before combining
  if (length(fwd) > 0) strand(fwd) <- "+"
  if (length(rev) > 0) strand(rev) <- "-"
  
  # Unify seqlevels so c() does not warn about mismatched seqinfo
  common_levels  <- union(seqlevels(fwd), seqlevels(rev))
  seqlevels(fwd) <- common_levels
  seqlevels(rev) <- common_levels
  
  gr <- c(fwd, rev)
  if (length(gr) == 0) return(NULL)
  
  # Annotate regions to genes
  anno <- annotatePeak(gr, TxDb = txdb, tssRegion = c(-500, 500),
                       annoDb = "org.Hs.eg.db", verbose = FALSE)
  df <- as.data.frame(anno)
  
  # Strand of the host gene
  gene_str <- ifelse(df$geneStrand == 1, "+",
                     ifelse(df$geneStrand == 2, "-", NA))
  
  # Sense / antisense relative to the host gene
  df$orientation <- ifelse(is.na(gene_str), NA,
                           ifelse(df$strand == gene_str, "sense", "antisense"))
  
  # Return the table, the csAnno object and the % of sense regions
  list(df = df, anno = anno,
       pct_sense = round(100 * mean(df$orientation == "sense",
                                    na.rm = TRUE), 1))
}


# == 06. Gene List Helpers ====================================================
# Unique gene symbols of a given orientation from an annotation data.frame
genes_by_orientation <- function(df, orientation) {
  unique(df$SYMBOL[df$orientation == orientation & !is.na(df$SYMBOL)])
}

# Combine sense + antisense symbols of a label
genes_of <- function(gene_lists, prefix) {
  idx <- grepl(paste0("^", prefix, "_(sense|antisense)$"), names(gene_lists))
  unique(unlist(gene_lists[idx]))
}
