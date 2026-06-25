# == 01. Information ==========================================================
#' Build every figure and composed panel of the TFM from the data saved by the
#' analysis scripts. This is the only script that produces plots

# Load source files
source("src/config.R")
source("src/differential_annotation_functions.R")
source("src/figures_functions.R")

# Define paths for the script
DIFDIR <- file.path(results, "differential")
FEATDIR <- file.path(results, "features")
MAPDIR  <- file.path(results, "physical_mapping")
METADIR <- file.path(results, "metagene")


# == 02. Libraries ============================================================
library(ChIPseeker)
library(ggpubr)
library(patchwork)
library(ggVennDiagram)


# == 03. Load Data ============================================================
gene_lists   <- readRDS(file.path(DIFDIR,  "annotation_analysis.rds"))
orient_anno  <- readRDS(file.path(DIFDIR,  "orientation_annotation.rds"))
complete     <- readRDS(file.path(FEATDIR, "features_complete.rds"))
chrom_data   <- readRDS(file.path(MAPDIR,  "chrom_data.rds"))
corr_data    <- readRDS(file.path(MAPDIR,  "corr_data.rds"))
anno_long    <- readRDS(file.path(MAPDIR,  "anno_long.rds"))
intron_df    <- readRDS(file.path(MAPDIR,  "siUAP56_UP_intron.rds"))
metagene_tab <- readRDS(file.path(METADIR, "metagene_positions.rds"))

# Split feature table by orientation
sense_data     <- complete[complete$orientation == "sense", ]
antisense_data <- complete[complete$orientation == "antisense", ]


# == 04. Venn Diagrams ========================================================
# Combine sense + antisense per condition/direction
venn_UP <- list(siC        = genes_of(gene_lists, "siC"),
                siUAP56_UP = genes_of(gene_lists, "siUAP56_UP"),
                siBRG1_UP  = genes_of(gene_lists, "siBRG1_UP"))

p_UP <- ggVennDiagram(venn_UP, label_alpha = 0, label = "count",
                      color = "black", lwd = 0.8, lty = 1) +
  scale_fill_gradient(low = "#FFFFFF", high = "#FF8C00") +
  theme(legend.position = "right")

ggsave(file.path(images, "figure1.jpg"),
       width = 7, height = 6, dpi = 300)


# == 05. Features Panel =======================================================
features <- make_compositional(sense_data, antisense_data)

ggsave(file.path(images, "figure2.jpg"),
       features, width = 14, height = 14, dpi = 300)


# == 06. Genomic Annotation and Metagene ======================================
# Feature distribution (plotAnnoBar) in the fixed group order
bar_labels <- c("siC", "siUAP56 UP", "siUAP56 DOWN", "siBRG1 UP", "siBRG1 DOWN")
anno_bars  <- lapply(orient_anno[group_levels], function(x) x$anno)
names(anno_bars) <- bar_labels

A <- plotAnnoBar(anno_bars) +
  xlab(NULL) + ylab("Percentage (%)") +
  theme_tfm() +
  theme(legend.position = "right", legend.title = element_blank())

# Metagenes plots
B <- metagene_panel(metagene_tab, "sense", "UP")
C <- metagene_panel(metagene_tab, "sense", "DOWN")
D <- metagene_panel(metagene_tab, "antisense", "UP")
E <- metagene_panel(metagene_tab, "antisense", "DOWN")

anno_meta_sense <- A / (B | C) / (D | E) + plot_annotation(tag_levels = "A",
                                                           tag_prefix = "(",
                                                           tag_suffix = ")")
ggsave(file.path(images, "figure3.jpg"),
       anno_meta_sense, width = 12, height = 15, dpi = 300)


# == 07. Chromosomal Distribution and GC Correlation ==========================
# Chromosome fraction, double panel (siUAP56 + siBRG1), UP/DOWN colored
A <- ggplot(chrom_data, aes(x = chrom, y = fraction, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~condition, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = direction_palette) +
  labs(x = NULL, y = "Fraction of chromosome genes", fill = NULL) +
  theme_tfm() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# GC vs fraction scatter
B <- plot_gc_scatter(corr_data)

chrom_distr_gc <- A / B + plot_annotation(tag_levels = "A",
                                          tag_prefix = "(",
                                          tag_suffix = ")")

ggsave(file.path(images, "figure4.jpg"),
       chrom_distr_gc, width = 10, height = 14, dpi = 300)


# == 08. UAP56 Long-gene Specificity Control ==================================
# Feature distribution restricted to long genes (siC vs siUAP56_UP)
names(anno_long)[1] <- "siC_75p"

A <- plotAnnoBar(anno_long) +
  xlab(NULL) +
  theme_tfm() +
  theme(legend.position = "right")

# Gene size by intronic localization in siUAP56_UP
B <- plot_intron_size(intron_df)

intron_control <- A / B + plot_annotation(tag_levels = "A",
                                          tag_prefix = "(",
                                          tag_suffix = ")")
ggsave(file.path(images, "figure5.jpg"), intron_control,
       width = 8, height = 10, dpi = 300)
