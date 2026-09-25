# ═══════════════════════════════════════════════════════════════════════════
# One set of colours for every figure, so figures can be read side by side.
# ═══════════════════════════════════════════════════════════════════════════
suppressMessages({library(ComplexHeatmap); library(circlize); library(grid)})
ht_opt$message <- FALSE

endotype_colours <- function(levels) {
  setNames(c("#D55E7A", "#4E9A06", "#3465A4", "#C4A000", "#75507B", "#CE5C00")[seq_along(levels)], levels)
}
SITE_COL  <- c(A = "#1b9e77", B = "#d95f02", C = "#7570b3")
SPLIT_COL <- c(discovery = "grey30", validation = "grey70", healthy = "white")
STAGE_COL <- setNames(hcl.colors(5, "YlOrRd", rev = TRUE), c("none", "mild", "moderate", "high", "very high"))
NEPH_COL  <- c(NoLN = "grey85", Mesan = "#c6dbef", Membr = "#6baed6",
               Prolif = "#08519c", "Proli+Membr" = "#54278f")
MDG_COL   <- setNames(hcl.colors(9, "Set 2"), as.character(0:8))
IFN_COL   <- colorRamp2(c(-1, 0, 2, 4), c("#2166AC", "white", "#F4A582", "#B2182B"))
VISIT_COL <- colorRamp2(c(1, 22), c("white", "#3f007d"))
Z_COL     <- colorRamp2(c(-2.5, 0, 2.5), c("#2166AC", "white", "#B2182B"))
