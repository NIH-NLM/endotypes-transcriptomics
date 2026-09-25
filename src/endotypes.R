# ═══════════════════════════════════════════════════════════════════════════
# Endotype discovery helpers, ported from reference/packages/endotypes/core.py.
#
# Every matrix here is GENES x SAMPLES, the Bioconductor convention, unless a
# function says otherwise.
#
# Rules that protect against circular inference:
#   * patients, not samples, are split into discovery and validation, so the
#     repeated visits of one patient never sit on both sides;
#   * gene selection, centring and clustering use discovery patients only;
#   * every other sample is labelled by nearest centroid, never by the
#     clustering itself.
# ═══════════════════════════════════════════════════════════════════════════

# ── disease stage ─────────────────────────────────────────────────────────
# SLEDAI (Systemic Lupus Erythematosus Disease Activity Index) category per
# visit. The cut points are the usual ones: 0, 1-5, 6-10, 11-19, 20 and above.
SLEDAI_LEVELS <- c("none", "mild", "moderate", "high", "very high")

sledai_category <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  factor(cut(x, c(-Inf, 0, 5, 10, 19, Inf), labels = SLEDAI_LEVELS),
         levels = SLEDAI_LEVELS, ordered = TRUE)
}

# ── interferon score ──────────────────────────────────────────────────────
# Each gene is z-scored against the healthy controls' mean and standard
# deviation, and the score is the median z across the genes. Without at least
# three controls, genes are z-scored across all samples instead and the basis
# column says so.
ifn_score <- function(E, genes, reference) {
  present <- intersect(genes, rownames(E))
  if (!length(present)) stop("none of the interferon genes are in the matrix", call. = FALSE)
  sub <- E[present, , drop = FALSE]
  ref <- if (sum(reference) >= 3) sub[, reference, drop = FALSE] else sub
  mu  <- rowMeans(ref)
  sd  <- apply(ref, 1, sd); sd[sd == 0] <- NA
  z   <- (sub - mu) / sd
  data.frame(sample       = colnames(E),
             ifn_score    = apply(z, 2, median, na.rm = TRUE),
             n_genes_used = colSums(!is.na(z)),
             basis        = if (sum(reference) >= 3) "healthy controls" else "all samples",
             row.names    = colnames(E))
}

# ── consensus clustering (Monti et al. 2003) ──────────────────────────────
# x is SAMPLES x FEATURES here, because k-means clusters rows.
# Each resample draws a fraction of the samples without replacement and runs
# k-means. The consensus value for a pair is the fraction of the resamples that
# drew both in which they landed in the same cluster.
consensus_matrix <- function(x, k, n_resamples = 100, frac = 0.8) {
  n <- nrow(x)
  together <- sampled <- matrix(0, n, n, dimnames = list(rownames(x), rownames(x)))
  m <- max(k + 1, round(frac * n))
  for (r in seq_len(n_resamples)) {
    idx <- sample.int(n, m)
    lab <- kmeans(x[idx, , drop = FALSE], centers = k, nstart = 5, iter.max = 50)$cluster
    together[idx, idx] <- together[idx, idx] + outer(lab, lab, "==")
    sampled[idx, idx]  <- sampled[idx, idx] + 1
  }
  C <- ifelse(sampled > 0, together / sampled, 0)
  diag(C) <- 1
  C
}

# PAC, the proportion of ambiguous clustering (Senbabaoglu et al. 2014): the
# share of pairs whose consensus value lies strictly between 0.1 and 0.9.
# A clean partition puts every pair near 0 or near 1. Lower is better.
pac <- function(C, lo = 0.1, hi = 0.9) {
  v <- C[upper.tri(C)]
  mean(v > lo & v < hi)
}

# Lowest PAC wins. When k values are within `tie` of the minimum they are
# equally stable, and the larger k is taken because it resolves more structure
# at no cost in stability.
choose_k <- function(pac_table, tie = 0.005) {
  tied <- pac_table$k[pac_table$pac <= min(pac_table$pac) + tie]
  max(tied)
}

# Final labels: average-linkage hierarchical clustering on 1 - consensus, cut
# at k, renamed E1..Ek by descending size so the labels are stable and readable.
labels_from_consensus <- function(C, k) {
  lab <- cutree(hclust(as.dist(1 - C), method = "average"), k = k)
  by_size <- names(sort(table(lab), decreasing = TRUE))
  setNames(factor(paste0("E", match(lab, by_size)), levels = paste0("E", seq_len(k))),
           rownames(C))
}

# ── nearest centroid ──────────────────────────────────────────────────────
# z is SAMPLES x GENES, centroids is ENDOTYPES x GENES, both on the same
# centring. Assignment is by Pearson correlation. The margin is the best
# correlation minus the second best: a small margin means the sample sits
# between two endotypes.
nearest_centroid <- function(z, centroids) {
  r    <- cor(t(z), t(centroids))
  ord  <- t(apply(r, 1, order, decreasing = TRUE))
  best <- ord[, 1]
  second <- if (ncol(r) > 1) r[cbind(seq_len(nrow(r)), ord[, 2])] else 0
  data.frame(endotype = factor(rownames(centroids)[best], levels = rownames(centroids)),
             margin   = r[cbind(seq_len(nrow(r)), best)] - second,
             row.names = rownames(z))
}

# ── pseudo-sites ──────────────────────────────────────────────────────────
# Adds a known site effect of the kind ComBat models: for every gene and site,
# an additive shift and a multiplicative change of spread around the site's own
# mean. Because we planted it, we can measure how much a correction removes.
add_site_effects <- function(E, site, shift_sd = 0.4, scale_sdlog = 0.25) {
  sites <- sort(unique(site))
  shift <- matrix(rnorm(nrow(E) * length(sites), 0, shift_sd), nrow(E),
                  dimnames = list(rownames(E), sites))
  scale <- matrix(rlnorm(nrow(E) * length(sites), 0, scale_sdlog), nrow(E),
                  dimnames = list(rownames(E), sites))
  out <- E
  for (s in sites) {
    j  <- site == s
    mu <- rowMeans(E[, j, drop = FALSE])
    out[, j] <- mu + (E[, j, drop = FALSE] - mu) * scale[, s] + shift[, s]
  }
  list(E = out, shift = shift, scale = scale)
}

# Share of each gene's variance explained by a grouping factor, averaged over
# genes. One number for "how much of this matrix is site?".
variance_explained <- function(E, group) {
  group <- factor(group)
  tot <- rowSums((E - rowMeans(E))^2)
  gm  <- sapply(levels(group), function(g) rowMeans(E[, group == g, drop = FALSE]))
  n   <- as.numeric(table(group))
  between <- rowSums(sweep((gm - rowMeans(E))^2, 2, n, "*"))
  mean(between / tot, na.rm = TRUE)
}

# ── federated location and scale correction ───────────────────────────────
# Each site computes three numbers per gene -- sample count, mean and variance --
# and sends only those. The coordinator pools them into one mean and one
# within-site variance per gene and sends those back. Each site then rescales
# its own data locally. No sample-level value leaves a site.
site_summary <- function(E_site)
  list(n = ncol(E_site), mean = rowMeans(E_site), var = apply(E_site, 1, var))

pool_summaries <- function(summaries) {
  n  <- sapply(summaries, `[[`, "n")
  mu <- sapply(summaries, `[[`, "mean")
  v  <- sapply(summaries, `[[`, "var")
  list(mean = drop(mu %*% n) / sum(n),
       var  = drop(v %*% (n - 1)) / (sum(n) - length(n)))
}

apply_pooled <- function(E_site, own, pooled)
  (E_site - own$mean) / sqrt(own$var) * sqrt(pooled$var) + pooled$mean

# ── agreement between two labellings ──────────────────────────────────────
# Adjusted Rand index (Hubert and Arabie 1985): 1 for identical partitions
# (whatever the label names), about 0 for agreement no better than chance.
ari <- function(a, b) {
  tab <- table(a, b)
  comb2 <- function(x) x * (x - 1) / 2
  s_ij <- sum(comb2(tab)); s_a <- sum(comb2(rowSums(tab))); s_b <- sum(comb2(colSums(tab)))
  expected <- s_a * s_b / comb2(sum(tab))
  (s_ij - expected) / ((s_a + s_b) / 2 - expected)
}

# ── applying a model at a site that did not build it ──────────────────────
# z is SAMPLES x GENES in the receiving site's own standardised units; the
# centroids are ENDOTYPES x GENES in the sender's. Only shared genes are used,
# and assignment is refused when fewer than `min_coverage` of the centroid
# genes are present (rule from reference/packages/endotypes/.../core.py).
assign_external <- function(z, centroids, min_coverage = 0.7) {
  shared <- intersect(colnames(centroids), colnames(z))
  coverage <- length(shared) / ncol(centroids)
  if (coverage < min_coverage)
    stop(sprintf("only %.0f%% of centroid genes present; refusing to assign", 100 * coverage),
         call. = FALSE)
  out <- nearest_centroid(z[, shared, drop = FALSE], centroids[, shared, drop = FALSE])
  attr(out, "coverage") <- coverage
  out
}
