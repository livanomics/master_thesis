# == 01. Information ==========================================================
#' Build the metagene position table: relative position of every R-loop within
#' its host gene, for each orientation x direction combination

# Load source files
source("src/config.R")
source("src/metagene_functions.R")

# Define paths for the script
DIFDIR <- file.path(results, "differential")
OUTDIR <- file.path(results, "metagene")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)


# == 02. Position Table =======================================================
orient_anno <- readRDS(file.path(DIFDIR, "orientation_annotation.rds"))

# Long table with all curves (siUAP56, siBRG1, siC) for every orientation x
# direction. FLANK = 500 bp mapped to a 0.3 fraction of the axis outside genes
metagene_tab <- build_metagene_table(orient_anno, flank_bp = 500,
                                     flank_rel = 0.3)

saveRDS(metagene_tab, file.path(OUTDIR, "metagene_positions.rds"))
