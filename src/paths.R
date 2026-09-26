# Where things live. Sourced by every notebook. Locations only; no analysis.
#
#   data/<GSE>/                one folder per GEO study (downloads)
#   data/ncbi/                 NCBI gene_info
#   data/run_artifacts/<GSE>/  everything a notebook writes; regenerable
#   genes/                     curated gene lists
#
# The root is found by walking up to the environment file, so this works from
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

# genes/<name>: one gene symbol per line; lines starting with # are notes
read_gene_set <- function(name) {
  x <- trimws(readLines(file.path(ROOT, "genes", name), warn = FALSE))
  unique(x[nzchar(x) & !startsWith(x, "#")])
}

SEED <- 20260928L
