# ═══════════════════════════════════════════════════════════════════════════
# Federation: what a site computes locally, and what a coordinator does with
# the little that is sent.
#
# The rule: no sample-level value leaves a site. A site may send
#   * per-gene counts, sums and sums of squares (one number per gene);
#   * per-cluster counts and per-cluster sums of a gene vector;
#   * model parameters (centroids, loadings, correction coefficients);
#   * histograms and contingency counts.
# A gene x gene matrix (a Gram matrix) is never sent: with more genes than
# patients it can be inverted back towards the patients' data.
#
# Matrices are GENES x SAMPLES, as in src/endotypes.R, unless a function says
# otherwise. Every function whose name starts with `site_` runs at one site and
# sees only that site's data. Every function whose name starts with `pool_` or
# `coordinator_` sees only what the sites sent.
#
# Every message a site sends is recorded by `send()`, so step 17 can audit
# exactly what left each site, and in what shape.
# ═══════════════════════════════════════════════════════════════════════════

# ── the message log ───────────────────────────────────────────────────────
# send() returns its payload unchanged and appends a line to the log file for
# this step: which step, which site, what, and the dimensions of every piece.
MESSAGE_LOG <- NULL          # set by start_log() at the top of a notebook

start_log <- function(step) {
  MESSAGE_LOG <<- art("messages_%s.csv", step)
  if (file.exists(MESSAGE_LOG)) file.remove(MESSAGE_LOG)
  invisible(MESSAGE_LOG)
}

send <- function(payload, site, what, n_patients_at_site) {
  if (!is.null(MESSAGE_LOG)) {
    pieces <- if (is.list(payload)) payload else list(value = payload)
    rows <- do.call(rbind, lapply(names(pieces), function(nm) {
      d <- dim(pieces[[nm]]); if (is.null(d)) d <- length(pieces[[nm]])
      data.frame(site = site, what = what, piece = nm,
                 dims = paste(d, collapse = " x "),
                 numbers = prod(d),
                 # a piece "has a patient dimension" if one of its dimensions
                 # equals the number of patients at the site; step 17 checks it
                 patient_dimension = any(d == n_patients_at_site) && n_patients_at_site > 1)
    }))
    write.table(rows, MESSAGE_LOG, sep = ",", row.names = FALSE,
                col.names = !file.exists(MESSAGE_LOG), append = file.exists(MESSAGE_LOG))
  }
  payload
}

# ── sufficient statistics: n, sum, sum of squares per gene ────────────────
site_suff <- function(X) list(n = ncol(X), sum = rowSums(X), sumsq = rowSums(X^2))

pool_suff <- function(msgs) {
  n   <- sum(vapply(msgs, `[[`, 0, "n"))
  s   <- Reduce(`+`, lapply(msgs, `[[`, "sum"))
  ss  <- Reduce(`+`, lapply(msgs, `[[`, "sumsq"))
  mean <- s / n
  list(n = n, mean = mean,
       var = (ss - n * mean^2) / (n - 1))      # total variance, the same as var() on all samples
}

# Pooled WITHIN-site variance: what ComBat calls var.pooled. Each site's sum of
# squared deviations from its own mean, added up, divided by the total count.
pool_within_var <- function(msgs) {
  n  <- sum(vapply(msgs, `[[`, 0, "n"))
  within <- Reduce(`+`, lapply(msgs, function(m) m$sumsq - m$sum^2 / m$n))
  within / n
}

# ── federated ComBat (Johnson et al. 2007), exact ─────────────────────────
# ComBat has two kinds of quantity:
#   * global: the grand mean and the pooled within-batch variance of each gene.
#     These need every batch, but only through n, sum and sum of squares.
#   * local: each batch's shift (gamma) and scale (delta), and the empirical
#     Bayes priors, which ComBat estimates across genes WITHIN one batch.
# So ComBat federates exactly: one round of summaries, then everything else
# at the site. The code below follows sva::ComBat(mod = NULL, par.prior = TRUE)
# line for line, so the gate in step 06 can demand equality.

combat_aprior <- function(d) { m <- mean(d); s2 <- var(d); (2 * s2 + m^2) / s2 }
combat_bprior <- function(d) { m <- mean(d); s2 <- var(d); (m * s2 + m^3) / s2 }

site_combat_adjust <- function(X, grand_mean, var_pooled, conv = 1e-4) {
  s   <- (X - grand_mean) / sqrt(var_pooled)             # standardise with the global values
  g_hat <- rowMeans(s)                                   # this site's shift, per gene
  d_hat <- apply(s, 1, var)                              # this site's scale, per gene
  g_bar <- mean(g_hat); t2 <- var(g_hat)                 # priors across genes, local
  a <- combat_aprior(d_hat); b <- combat_bprior(d_hat)
  n <- ncol(s)
  g_old <- g_hat; d_old <- d_hat; change <- 1
  while (change > conv) {                                # sva:::it.sol
    g_new  <- (t2 * n * g_hat + d_old * g_bar) / (t2 * n + d_old)
    sum2   <- rowSums((s - g_new)^2)
    d_new  <- (0.5 * sum2 + b) / (n / 2 + a - 1)
    change <- max(abs(g_new - g_old) / g_old, abs(d_new - d_old) / d_old)
    g_old <- g_new; d_old <- d_new
  }
  adj <- (s - g_new) / sqrt(d_new) * sqrt(var_pooled) + grand_mean
  list(X = adj, gamma_star = g_new, delta_star = d_new)
}

# ── federated principal components ────────────────────────────────────────
# Subspace iteration. The coordinator holds a genes x q basis V. Each site
# centres its own patients on the pooled mean and returns Xc' (Xc V), a
# genes x q matrix: never a genes x genes matrix, never a patient row. The sum
# over sites is the pooled X'X V. Orthonormalise, repeat. q is larger than the
# number of components kept so the leading ones converge quickly.
site_xtxv <- function(X, mean, V) { Xc <- t(X - mean); crossprod(Xc, Xc %*% V) }

federated_pca <- function(sites, mean, n_total, k = 10, oversample = 10,
                          tol = 1e-12, max_iter = 2000, seed = 1L, step = NULL) {
  set.seed(seed)
  q <- k + oversample
  V <- qr.Q(qr(matrix(rnorm(length(mean) * q), length(mean), q)))
  for (it in seq_len(max_iter)) {
    msgs <- lapply(names(sites), function(s)
      send(site_xtxv(sites[[s]], mean, V), s, "X'XV", ncol(sites[[s]])))
    W    <- Reduce(`+`, msgs)
    Vnew <- qr.Q(qr(W))
    # the change in the leading k-dimensional subspace
    delta <- 1 - min(svd(crossprod(Vnew[, 1:k], V[, 1:k]))$d)
    V <- Vnew
    if (delta < tol) break
  }
  # Rayleigh-Ritz: rotate the basis to the eigenvectors of V'X'XV
  msgs <- lapply(names(sites), function(s)
    send(site_xtxv(sites[[s]], mean, V), s, "X'XV", ncol(sites[[s]])))
  H  <- crossprod(V, Reduce(`+`, msgs))
  e  <- eigen((H + t(H)) / 2, symmetric = TRUE)
  rot <- (V %*% e$vectors)[, 1:k]
  rot <- sweep(rot, 2, sign(rot[which.max(abs(rot[, 1])), ] + 1e-300), "*")
  list(rotation = rot, sdev = sqrt(e$values[1:k] / (n_total - 1)), iterations = it)
}

# ── federated k-means (Lloyd) ─────────────────────────────────────────────
# The coordinator sends k centroids. Each site assigns its own points to the
# nearest centroid and returns, per cluster, a count and the sum of the points,
# and its within-cluster sum of squares. The coordinator averages. Repeat until
# no site changes an assignment.
#
# `sites` is a list of SAMPLES x FEATURES matrices (rows are points).
site_lloyd_step <- function(P, C) {
  d   <- as.matrix(dist(rbind(C, P)))[-(1:nrow(C)), 1:nrow(C), drop = FALSE]
  lab <- max.col(-d, ties.method = "first")
  k   <- nrow(C)
  list(count = tabulate(lab, k),
       sum   = t(sapply(seq_len(k), function(j) colSums(P[lab == j, , drop = FALSE]))),
       wss   = sum(d[cbind(seq_along(lab), lab)]^2),
       labels = lab)                                    # stays at the site
}

federated_kmeans <- function(sites, init, max_iter = 100, log = FALSE) {
  C <- init; prev <- NULL
  for (it in seq_len(max_iter)) {
    local <- lapply(sites, site_lloyd_step, C = C)
    if (log) for (s in names(local))
      send(local[[s]][c("count", "sum", "wss")], s, "k-means cluster sums", nrow(sites[[s]]))
    cnt <- Reduce(`+`, lapply(local, `[[`, "count"))
    sm  <- Reduce(`+`, lapply(local, `[[`, "sum"))
    keep <- cnt > 0
    C[keep, ] <- sm[keep, , drop = FALSE] / cnt[keep]
    labs <- lapply(local, `[[`, "labels")
    if (!is.null(prev) && identical(labs, prev)) break
    prev <- labs
  }
  final <- lapply(sites, site_lloyd_step, C = C)
  list(centers = C, labels = lapply(final, `[[`, "labels"),
       wss = sum(vapply(final, `[[`, 0, "wss")),
       size = Reduce(`+`, lapply(final, `[[`, "count")), iterations = it)
}

# Starting centroids without seeing any patient: random points drawn from the
# pooled distribution of each feature (pooled mean and SD, from summaries).
random_init <- function(k, mean, sd) t(replicate(k, rnorm(length(mean), mean, sd)))

# Several random starts; keep the run with the smallest total within-cluster
# sum of squares (each site reports its own as one number).
federated_kmeans_best <- function(sites, k, mean, sd, nstart = 5, log = FALSE) {
  best <- NULL
  for (r in seq_len(nstart)) {
    fit <- federated_kmeans(sites, random_init(k, mean, sd), log = log)
    if (any(fit$size == 0)) next
    if (is.null(best) || fit$wss < best$wss) best <- fit
  }
  best
}

# ── federated consensus clustering ────────────────────────────────────────
# Monti consensus, with each resample drawn WITHIN each site (80% of that
# site's points). Each site keeps co-clustering counts for its own pairs of
# patients only; pairs across sites are never observed. For PAC a site sends
# only a histogram of its consensus values.
BINS <- seq(0, 1, by = 0.01)

federated_consensus <- function(sites, k, mean, sd, n_resamples = 100, frac = 0.8, nstart = 5) {
  together <- lapply(sites, function(P) matrix(0, nrow(P), nrow(P), dimnames = list(rownames(P), rownames(P))))
  sampled  <- together
  for (r in seq_len(n_resamples)) {
    idx <- lapply(sites, function(P) sort(sample.int(nrow(P), max(2, round(frac * nrow(P))))))
    sub <- Map(function(P, i) P[i, , drop = FALSE], sites, idx)
    fit <- federated_kmeans_best(sub, k, mean, sd, nstart = nstart)
    for (s in names(sites)) {
      i <- idx[[s]]; lab <- fit$labels[[s]]
      together[[s]][i, i] <- together[[s]][i, i] + outer(lab, lab, "==")
      sampled[[s]][i, i]  <- sampled[[s]][i, i] + 1
    }
  }
  cons <- Map(function(t, s) { C <- ifelse(s > 0, t / s, 0); diag(C) <- 1; C }, together, sampled)
  hist <- lapply(names(cons), function(s) {
    v <- cons[[s]][upper.tri(cons[[s]])]
    send(list(counts = as.numeric(table(cut(v, BINS, include.lowest = TRUE)))),
         s, sprintf("consensus histogram k=%d", k), nrow(sites[[s]]))$counts
  })
  names(hist) <- names(cons)
  list(consensus = cons, histogram = hist)
}

# PAC from pooled histograms: the share of pairs with consensus strictly
# between 0.1 and 0.9. Bins are 0.01 wide, and the bin edges sit exactly on 0.1
# and 0.9, so this equals the pair-level PAC up to values on the edges.
pac_from_hist <- function(h) {
  mids <- (head(BINS, -1) + tail(BINS, -1)) / 2
  sum(h[mids > 0.1 & mids < 0.9]) / sum(h)
}

# ── a site's picture of its own samples, without the samples ──────────────
# A smoothed two-dimensional density on a fixed grid (a Gaussian kernel, as in
# MASS::kde2d). A site sends the grid; the coordinator draws contours. No
# sample's coordinates leave the site.
density_grid <- function(x, y, lims, n = 60, h = c(diff(lims[1:2]), diff(lims[3:4])) / 12) {
  gx <- seq(lims[1], lims[2], length.out = n); gy <- seq(lims[3], lims[4], length.out = n)
  ax <- outer(gx, x, function(g, v) dnorm((g - v) / h[1]))
  ay <- outer(gy, y, function(g, v) dnorm((g - v) / h[2]))
  list(x = gx, y = gy, z = tcrossprod(ax, ay) / (length(x) * h[1] * h[2]))
}
