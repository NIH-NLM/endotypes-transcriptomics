# ═══════════════════════════════════════════════════════════════════════════
# Platforms: each study is normalised with its own platform's method, and
# studies meet only in shared units.
#
#   RNA_array     GSE65391   log2 intensities as deposited -> limma
#   bulk_RNA_seq  GSE232381  raw counts -> edgeR TMM -> voom log2 CPM
#   scRNA_seq     GSE135779  UMI counts summed per donor (pseudobulk) -> as bulk
#
# Values from different platforms are never put in one matrix. What crosses
# between studies is a gene list, gene-wise standardised units, centroids and
# scale-free scores.
#
# Genes are keyed by NCBI Entrez GeneID, the one identifier all three studies
# can reach: the array annotation carries it, the NCBI RNA-seq counts are keyed
# by it, and NCBI gene_info maps Ensembl IDs (single-cell) to it. Each gene is
# labelled with its current NCBI symbol.
# ═══════════════════════════════════════════════════════════════════════════

GENE_INFO_URL <- "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"

read_gene_info <- function() {
  f <- fetch(GENE_INFO_URL, file.path(NCBI, "Homo_sapiens.gene_info.gz"))
  gi <- utils::read.delim(f, quote = "", comment.char = "", check.names = FALSE,
                          colClasses = "character")
  names(gi)[1] <- "tax_id"
  ens <- regmatches(gi$dbXrefs, regexpr("Ensembl:ENSG[0-9]+", gi$dbXrefs))
  has <- grepl("Ensembl:ENSG[0-9]+", gi$dbXrefs)
  gi$ensembl <- NA_character_
  gi$ensembl[has] <- sub("Ensembl:", "", ens)
  gi[, c("GeneID", "Symbol", "ensembl", "type_of_gene")]
}

# Collapse a matrix with duplicated row keys to one row per key, keeping the
# row with the highest mean (the brightest probe, or the most-counted gene).
collapse_max_mean <- function(X, key) {
  ok  <- !is.na(key) & nzchar(key)
  X   <- X[ok, , drop = FALSE]; key <- key[ok]
  pick <- tapply(seq_len(nrow(X)), key, function(i) i[which.max(rowMeans(X[i, , drop = FALSE]))])
  out <- X[unlist(pick), , drop = FALSE]
  rownames(out) <- names(pick)
  out
}

# Counts -> log2 counts per million, with trimmed-mean-of-M-values (TMM) library
# normalisation and voom's precision weights discarded (we keep the values for
# scoring and clustering; the weights matter for differential expression).
counts_to_logcpm <- function(counts, group = NULL) {
  suppressMessages({library(edgeR); library(limma)})
  d <- DGEList(counts = counts, group = group)
  keep <- filterByExpr(d, group = group)
  d <- calcNormFactors(d[keep, , keep.lib.sizes = FALSE], method = "TMM")
  design <- if (is.null(group)) matrix(1, ncol(d), 1) else model.matrix(~ group)
  list(E = voom(d, design)$E, genes_kept = sum(keep), genes_in = nrow(counts))
}

# Sum a 10x-style sparse counts matrix (genes x cells) into one column.
pseudobulk_one <- function(mtx_file) {
  suppressMessages(library(Matrix))
  m <- readMM(gzfile(mtx_file))
  Matrix::rowSums(m)
}

# Gene-wise standardisation against a reference set of samples of the SAME
# study. This is the shared unit: "reference SDs above the reference mean".
standardise_to_reference <- function(E, ref) {
  mu <- rowMeans(E[, ref, drop = FALSE])
  sd <- apply(E[, ref, drop = FALSE], 1, sd)
  sd[sd == 0] <- NA
  list(Z = (E - mu) / sd, mean = mu, sd = sd)
}
