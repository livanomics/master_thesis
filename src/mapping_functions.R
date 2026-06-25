# == 01. Information ==========================================================
#' Functions for the physical mapping analysis: SYMBOL/ENTREZ conversion,
#' per-chromosome gene fractions, chromosomal GC content
#' and long-gene annotation


# == 02. ID Conversion ========================================================
# SYMBOL -> ENTREZID (to map gene lists onto the TxDb)
sym_to_entrez <- function(syms) {
  m <- AnnotationDbi::select(org.Hs.eg.db, keys = syms,
                             keytype = "SYMBOL", columns = "ENTREZID")
  unique(m$ENTREZID[!is.na(m$ENTREZID)])
}


# == 03. Per-chromosome Fractions =============================================
# For every list: affected genes per chromosome/total genes of that chromosome
build_chrom_tab <- function(lists_entrez, canon) {
  # Chromosome of every gene in the TxDb (canonical chromosomes only)
  g <- genes(txdb)
  g <- g[as.character(seqnames(g)) %in% canon]
  gene_chrom <- data.frame(gene  = names(g),
                           chrom = as.character(seqnames(g)),
                           stringsAsFactors = FALSE)

  # Denominator: total genes per chromosome
  total <- table(factor(gene_chrom$chrom, levels = canon))

  # Numerator: affected genes per chromosome, for each list
  out <- data.frame()
  for (l in names(lists_entrez)) {
    aff   <- gene_chrom$chrom[gene_chrom$gene %in% lists_entrez[[l]]]
    n_aff <- table(factor(aff, levels = canon))
    out <- rbind(out, data.frame(
      list       = l,
      chrom      = canon,
      n_total    = as.integer(total),
      n_affected = as.integer(n_aff),
      fraction   = as.numeric(n_aff) / as.integer(total),
      stringsAsFactors = FALSE))
  }
  out$chrom <- factor(out$chrom, levels = canon)
  out
}

# GC content (%) of each canonical chromosome
chrom_gc_content <- function(genome, canon) {
  prop <- data.frame(chrom = canon, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(prop))) {
    af <- alphabetFrequency(getSeq(genome, prop$chrom[i]), baseOnly = TRUE)
    prop$GC[i] <- 100 * sum(af[c("G", "C")]) / sum(af[c("A", "C", "G", "T")])
  }
  prop
}


# == 04. Long-gene Annotation =================================================
# Re-annotate only the regions whose host gene is "long" (>= threshold)
long_genes_anno <- function(res, threshold) {
  df <- res$df
  df$gene_len <- df$geneEnd - df$geneStart
  gr <- res$anno@anno[which(df$gene_len >= threshold)]
  if (length(gr) == 0) return(NULL)
  annotatePeak(gr, TxDb = txdb, tssRegion = c(-500, 500), verbose = FALSE)
}
