# == 01. Information ==========================================================
#' Functions for creating plots and panels


# == 02. Features Plots =======================================================
# Pairs compared on the plots: each list vs siC, and UP vs DOWN per condition
make_comparisons <- function(data) {
  groups <- as.character(unique(data$group))
  comps  <- list()
  
  # Each list vs siC
  if ("siC" %in% groups) {
    diff_groups <- setdiff(groups, "siC")
    comps <- c(comps, lapply(diff_groups, function(g) c("siC", g)))
  }
  # UP vs DOWN within each condition
  for (smp in c("siUAP56", "siBRG1")) {
    a <- paste0(smp, "_UP"); b <- paste0(smp, "_DOWN")
    if (a %in% groups && b %in% groups)
      comps <- c(comps, list(c(a, b)))
  }
  comps
}

# Violin + boxplot of one variable across groups, with significance brackets
plot_feature <- function(data, var, ylab, title, comparisons, logy = FALSE) {
  data$group <- factor(data$group, levels = group_levels)
  p <- ggplot(data, aes(x = group, y = .data[[var]], fill = group)) +
    geom_violin(alpha = 0.30, color = NA, scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.18, outlier.size = 0.3, alpha = 0.95,
                 color = "grey20") +
    scale_fill_manual(values = group_palette) +
    labs(title = title, x = NULL, y = ylab) +
    theme_tfm(base_size = 13) +
    theme(legend.position = "none")
  if (logy) p <- p + scale_y_log10()
  p + stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                         label = "p.signif", size = 3.5)
}

# Three feature panels (GC, size, exons) for sense and antisense
make_compositional <- function(sense, antisense) {
  comps_sense <- make_comparisons(sense)
  comps_antis <- make_comparisons(antisense)
  A <- plot_feature(sense, "percentage_gene_gc_content", "GC content (%)",
                    "GC content", comps_sense)
  B <- plot_feature(sense, "length_kb", "Gene size (kb, log10)",
                    "Gene size", comps_sense, logy = TRUE)
  C <- plot_feature(sense, "exons_n", "Number of exons (log10)",
                    "Number of exons", comps_sense, logy = TRUE)
  D <- plot_feature(antisense, "percentage_gene_gc_content", "GC content (%)",
                    "GC content", comps_antis)
  E <- plot_feature(antisense, "length_kb", "Gene size (kb, log10)",
                    "Gene size", comps_antis, logy = TRUE)
  G <- plot_feature(antisense, "exons_n", "Number of exons (log10)",
                    "Number of exons", comps_antis, logy = TRUE)
  (A | B | C) / (D | E | G) + plot_annotation(tag_levels = "A",
                                              tag_prefix = "(",
                                              tag_suffix = ")")
}


# == 03. Metagene Plots =======================================================
# Density metagene for one orientation + direction, from the position table
metagene_panel <- function(metagene_tab, orientation, direction,
                           flank_rel = 0.3, min_n = 20) {
  df_plot <- metagene_tab[metagene_tab$orientation == orientation &
                            metagene_tab$direction == direction, ]
  if (nrow(df_plot) == 0) return(NULL)
  
  # Drop curves with too few points (unreliable density estimate)
  keep    <- names(which(table(df_plot$curve) >= min_n))
  df_plot <- df_plot[df_plot$curve %in% keep, ]
  if (nrow(df_plot) == 0) return(NULL)
  
  # Curve colors from the shared palette
  cols <- c(
    setNames(group_palette[[paste0("siUAP56_", direction)]],
             paste0("siUAP56 ", direction)),
    setNames(group_palette[[paste0("siBRG1_",  direction)]],
             paste0("siBRG1 ",  direction)),
    setNames(group_palette[["siC"]], "siC"))
  
  ggplot(df_plot, aes(x = pos, color = curve)) +
    geom_density(linewidth = 1, bounds = c(-flank_rel, 1 + flank_rel)) +
    geom_vline(xintercept = c(0, 1), linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    scale_color_manual(values = cols) +
    scale_x_continuous(
      breaks = c(-flank_rel / 2, 0, 0.5, 1, 1 + flank_rel / 2),
      labels = c("-500", "TSS", "body", "TTS", "+500"),
      limits = c(-flank_rel, 1 + flank_rel), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 3)) +
    labs(title = paste0("R-loop position — ", orientation, " ", direction),
         x = NULL, y = "Density", color = NULL) +
    theme_tfm()
}


# == 04. Chromosome Location Plots ============================================
# Chromosome fraction for one condition (UP vs DOWN)
plot_chrom <- function(chrom_data, cond) {
  sub <- chrom_data[chrom_data$condition == cond, ]
  ggplot(sub, aes(x = chrom, y = fraction, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = direction_palette) +
    labs(title = paste0("Fraction of genes with R-loops per chromosome (", cond, ")"),
         subtitle = "affected genes / total genes of the chromosome",
         x = NULL, y = "Fraction of chromosome genes", fill = NULL) +
    theme_tfm() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Scatter of chromosomal GC vs fraction of affected genes (faceted by list)
plot_gc_scatter <- function(corr_data) {
  corr_data$lab <- ifelse(corr_data$chrom %in% paste0("chr", c(19, 16, 17, 22)),
                          as.character(corr_data$chrom), "")
  ggplot(corr_data, aes(x = GC, y = fraction)) +
    geom_smooth(method = "lm", se = TRUE, color = "grey50",
                linetype = "dashed", linewidth = 0.5) +
    geom_point(aes(color = list), size = 2.5, alpha = 0.8) +
    ggrepel::geom_text_repel(aes(label = lab), size = 3, max.overlaps = 20) +
    facet_wrap(~list, scales = "free_y", ncol = 2) +
    labs(title = "Chromosomal GC vs fraction of genes with R-loops",
         x = "Chromosomal GC content (%)", y = "Fraction of genes affected") +
    theme_tfm() +
    theme(legend.position = "none")
}


# == 05. Long-genes Control Plot ==============================================
# Gene size by intronic vs non-intronic R-loop localization
plot_intron_size <- function(df, title = NULL) {
  ggplot(df, aes(x = is_intron, y = gene_kb, fill = is_intron)) +
    geom_violin(alpha = 0.3, color = NA, scale = "width") +
    geom_boxplot(width = 0.2, outlier.size = 0.3) +
    scale_y_log10() +
    scale_x_discrete(labels = c("FALSE" = "Not intron", "TRUE" = "Intron")) +
    scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#cb181d")) +
    labs(title = title, x = NULL, y = "Gene size (kb, log10)") +
    theme_tfm() +
    theme(legend.position = "none") +
    stat_compare_means(method = "wilcox.test", label = "p.format")
}
