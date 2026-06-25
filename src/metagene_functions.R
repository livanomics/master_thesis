# == 01. Information ==========================================================
#' Functions for the metagene analysis: relative position of each region within
#' its host gene, a long position table for all combinations, and the density
#' (metagene) plot per orientation / direction


# == 02. Position Profile =====================================================
# Relative position of each region within its gene (0 = TSS, 1 = TTS)
position_from_df <- function(df, flank_bp = 0, flank_rel = 0.3) {
  df <- df[!is.na(df$geneStart) & !is.na(df$geneEnd), ]
  if (nrow(df) == 0) return(numeric(0))

  gs <- as.character(df$geneStrand)
  gs[gs %in% c("1", "+")] <- "+"; gs[gs %in% c("2", "-")] <- "-"

  center <- (df$start + df$end) / 2
  glen   <- df$geneEnd - df$geneStart
  pos    <- numeric(nrow(df)); keep <- rep(TRUE, nrow(df))

  for (i in seq_len(nrow(df))) {
    plus   <- gs[i] == "+"
    inside <- center[i] >= df$geneStart[i] & center[i] <= df$geneEnd[i]
    if (inside) {
      pos[i] <- if (plus) (center[i] - df$geneStart[i]) / glen[i]
                else      (df$geneEnd[i] - center[i]) / glen[i]
    } else if (flank_bp > 0) {
      if (plus) {
        if (center[i] < df$geneStart[i]) {
          d <- df$geneStart[i] - center[i]
          if (d > flank_bp) keep[i] <- FALSE else pos[i] <- -d / flank_bp * flank_rel
        } else {
          d <- center[i] - df$geneEnd[i]
          if (d > flank_bp) keep[i] <- FALSE else pos[i] <- 1 + d / flank_bp * flank_rel
        }
      } else {
        if (center[i] > df$geneEnd[i]) {
          d <- center[i] - df$geneEnd[i]
          if (d > flank_bp) keep[i] <- FALSE else pos[i] <- -d / flank_bp * flank_rel
        } else {
          d <- df$geneStart[i] - center[i]
          if (d > flank_bp) keep[i] <- FALSE else pos[i] <- 1 + d / flank_bp * flank_rel
        }
      }
    } else keep[i] <- FALSE
  }
  pos[keep]
}

# Positions of one list filtered by orientation, tagged with a curve name
positions_of <- function(orient_anno, label, orientation, curve_name,
                         flank_bp = 500, flank_rel = 0.3) {
  if (is.null(orient_anno[[label]])) return(NULL)
  df  <- orient_anno[[label]]$df
  sub <- df[!is.na(df$orientation) & df$orientation == orientation, ]
  if (nrow(sub) == 0) return(NULL)
  p <- position_from_df(sub, flank_bp = flank_bp, flank_rel = flank_rel)
  if (length(p) == 0) return(NULL)
  data.frame(pos = p, curve = curve_name)
}


# == 03. Position Table =======================================================
# Long table of relative positions for every orientation x direction, with the
# three curves (siUAP56, siBRG1, siC)
build_metagene_table <- function(orient_anno, flank_bp = 500, flank_rel = 0.3) {
  out <- data.frame()
  for (ori in c("sense", "antisense")) {
    for (dir in c("UP", "DOWN")) {
      uap <- paste0("siUAP56_", dir)
      brg <- paste0("siBRG1_",  dir)
      part <- rbind(
        positions_of(orient_anno, uap,   ori, paste0("siUAP56 ", dir),
                     flank_bp, flank_rel),
        positions_of(orient_anno, brg,   ori, paste0("siBRG1 ",  dir),
                     flank_bp, flank_rel),
        positions_of(orient_anno, "siC", ori, "siC", flank_bp, flank_rel))
      if (!is.null(part) && nrow(part) > 0) {
        part$orientation <- ori
        part$direction   <- dir
        out <- rbind(out, part)
      }
    }
  }
  out
}
