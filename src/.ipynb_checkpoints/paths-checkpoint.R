# Paths.R
#   Sets the locations for all the files:
#
#   data/<GSE>/                one folder per GEO study
#   data/ncbi/                 NCBI gene_info 
#   data/run_artifacts/<GSE>/  the output of all noteboos, organized by GEO study, reenerable
#   figures/<GSE>/             figures are both inline and saved, regenerable
#   genes/                     curated gene lists
#
# `..` refers to the parent to a current directory
# `.`  refers to the current directory
#  
# ipynb/ or from the repository root.

.find_root <- function(start = getwd()) {
  d <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(d, "endotypes-transcriptomics.yml"))) return(d)
    up <- dirname(d)
    if (identical(up, d)) stop("endotypes-transcriptomics.yml not found above ", start, call. = FALSE)
    d <- up
  }
}

ROOT <- .find_root()
DATA <- file.path(ROOT, "data")
NCBI <- file.path(DATA, "ncbi")

# raw("GSE65391", "file.gz")  ->  data/GSE65391/file.gz
raw <- function(gse, file = "") {
  d <- file.path(DATA, gse); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  if (nzchar(file)) file.path(d, file) else d
}

# art("GSE65391", "file.rds")  ->  data/run_artifacts/GSE65391/file.rds
art <- function(gse, file = "") {
  d <- file.path(DATA, "run_artifacts", gse); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  if (nzchar(file)) file.path(d, file) else d
}

# fig("GSE65391", "07_heatmap_hclust.png")  ->  figures/GSE65391/07_heatmap_hclust.png
fig <- function(gse, file) {
  d <- file.path(ROOT, "figures", gse); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  file.path(d, file)
}

# genes/<name>: one gene symbol per line; lines starting with # are notes
read_gene_set <- function(name) {
  x <- trimws(readLines(file.path(ROOT, "genes", name), warn = FALSE))
  unique(x[nzchar(x) & !startsWith(x, "#")])
}

SEED <- 20260928L
