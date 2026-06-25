# == 01. Information ==========================================================
#' Main configuration for the entire program: paths, genome annotation and the
#' shared palette/theme/order used across every figure

author  <- "José Livan Vargas Castro"
version <- "1.0.0"


# == 02. Global Libraries =====================================================
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)


# == 03. Paths ================================================================
root    <- "/home/josvarcas/Archivos/master/c0_TFMEXP/TFM_DRIPc-seq"
data    <- file.path(root, "data")
results <- file.path(root, "results")
images  <- file.path(results, "images")


# == 04. Genome Annotation ====================================================
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


# == 05. Shared Palette and Order =============================================
# Color per group
group_palette <- c(
  "siC"          = "#6a51a3",
  "siUAP56_UP"   = "#cb181d",
  "siUAP56_DOWN" = "#fb6a4a",
  "siBRG1_UP"    = "#238b45",
  "siBRG1_DOWN"  = "#74c476"
)

# Color per direction
direction_palette <- c(UP = "#cb181d", DOWN = "#2171b5")

# Fixed group order for the x axis
group_levels <- c("siC", "siUAP56_UP", "siUAP56_DOWN",
                  "siBRG1_UP", "siBRG1_DOWN")


# == 06. Shared Theme =========================================================
# Clean, paper-style theme used by all figures
theme_tfm <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      legend.position  = "bottom",
      plot.title       = element_text(face = "bold", size = base_size + 1),
      axis.text.x      = element_text(angle = 30, hjust = 1),
      panel.grid.minor = element_blank()
    )
}
