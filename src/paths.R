# ═══════════════════════════════════════════════════════════════════════════
# WHERE THINGS LIVE -- one definition, sourced by every notebook and script.
#
#   cohorts/            COMMITTED. The frozen site assignment and the frozen
#                       discovery/validation split. Redrawing either moves every
#                       downstream number, so a run reads them back, never redraws.
#   genes/              COMMITTED. Curated gene sets the analysis is interpreted
#                       against, not derived from (interferon scores).
#   figures/            COMMITTED. The rendered heatmaps, kept for sharing.
#   data/<accession>/   GITIGNORED. One folder per GEO (Gene Expression Omnibus)
#                       series: GSE65391 (RNA_array), GSE232381 (bulk_RNA_seq),
#                       GSE135779 (scRNA_seq). data/ncbi/ holds NCBI gene_info.
#   data/run_artifacts/ GITIGNORED. Everything a run can regenerate: .rds and
#                       derived .csv tables. Safe to delete.
#   src/                R that is not a notebook.
#   ipynb/              notebooks, and nothing else.
#   reference/          the earlier Nextflow scaffold. Logic source, not run.
#
# The root is found by walking up for the environment file, so this works
# from ipynb/, from src/, or from the repository root:
#
#   source("../src/paths.R")     # from a notebook
#   source("src/paths.R")        # from the root
# ═══════════════════════════════════════════════════════════════════════════

.find_root <- function(start = getwd()) {
  d <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(d, "endotypes-transcriptomics.yml"))) return(d)
    up <- dirname(d)
    if (identical(up, d))
      stop("endotypes-transcriptomics.yml not found above '", start,
           "' -- are you inside the repository?", call. = FALSE)
    d <- up
  }
}

ROOT      <- .find_root()
DATA      <- file.path(ROOT, "data")
NCBI      <- file.path(DATA, "ncbi")
COHORTS   <- file.path(ROOT, "cohorts")
GENES     <- file.path(ROOT, "genes")
FIGURES   <- file.path(ROOT, "figures")
ARTIFACTS <- file.path(ROOT, "data", "run_artifacts")

for (d in c(NCBI, COHORTS, FIGURES, ARTIFACTS))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# raw("GSE65391", "GSE65391_series_matrix.txt.gz") -- the study's own folder,
# created on first use.
raw <- function(accession, fmt = "", ...) {
  d <- file.path(DATA, accession)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  if (nzchar(fmt)) file.path(d, sprintf(fmt, ...)) else d
}
coh <- function(fmt, ...) file.path(COHORTS,   sprintf(fmt, ...))
gen <- function(fmt, ...) file.path(GENES,     sprintf(fmt, ...))
fig <- function(fmt, ...) file.path(FIGURES,   sprintf(fmt, ...))
art <- function(fmt, ...) file.path(ARTIFACTS, sprintf(fmt, ...))

# A gene set file is one HGNC symbol per line. Lines starting with # are notes.
read_gene_set <- function(name) {
  x <- trimws(readLines(gen(name), warn = FALSE))
  x <- x[nzchar(x) & !startsWith(x, "#")]
  unique(x)
}

# Every notebook seeds from here, so one number controls every random draw.
SEED <- 20260928L

# The three simulated sites of GSE65391.
SITES <- c("A", "B", "C")

# Where one site's data lives after step 04. A notebook touching site data
# reads one of these at a time.
site_file <- function(step, site) art("step%s_site_%s.rds", step, site)

# Download once; a second run touches no network.
fetch <- function(url, dest) {
  if (!file.exists(dest)) {
    options(timeout = max(3600, getOption("timeout")))
    download.file(url, dest, mode = "wb", quiet = TRUE)
  }
  stopifnot(file.exists(dest))
  invisible(dest)
}
