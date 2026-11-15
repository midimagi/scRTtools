# IZ heterogeneity phases 4-7: clustering, nearby co-usage, non-negative
# tensor factorisation and epigenomic integration. Migrated from
# jupyter/iz_heterogeneity_phase4_7.R; dplyr steps are rewritten in base R and
# set.seed() calls replaced by withr::with_seed() (same random streams).

#' Winsorised, median-imputed and scaled metric matrix
#' @keywords internal
#' @noRd
.metric_matrix <- function(x, metrics, winsor = c(0.01, 0.99)) {
  m <- as.matrix(x[, metrics, drop = FALSE])
  storage.mode(m) <- "double"
  for (j in seq_len(ncol(m))) {
    q <- stats::quantile(m[, j], winsor, na.rm = TRUE, names = FALSE)
    m[, j] <- pmin(pmax(m[, j], q[1]), q[2])
    m[!is.finite(m[, j]), j] <- stats::median(m[, j], na.rm = TRUE)
  }
  scale(m)
}

#' Index of the nearest centre for each row
#' @keywords internal
#' @noRd
.nearest_center <- function(x, centers) {
  d <- vapply(seq_len(nrow(centers)), function(k)
    rowSums((x - matrix(centers[k, ], nrow(x), ncol(x), byrow = TRUE))^2),
    numeric(nrow(x)))
  max.col(-d, ties.method = "first")
}

#' Cluster initiation zones by heterogeneity metrics
#'
#' Metrics are winsorised (1%-99%), median-imputed and scaled, then clustered
#' with k-means for every `k` in `k_range`; `k` defaults to the value with the
#' highest mean silhouette width (computed on at most `silhouette_n` IZs).
#' Ward hierarchical clustering cut at `k` and a PCA are returned as well.
#'
#' @param iz_summary Per-IZ table with `IZ_id` and the metric columns (e.g.
#'   [assemble_iz_summary()]).
#' @param metrics Metric column names.
#' @param k_range Candidate numbers of clusters.
#' @param k `NULL` or the number of clusters to use (must be in `k_range`).
#' @param nstart `nstart` of `stats::kmeans()`.
#' @param silhouette_n Maximum number of IZs used for silhouettes.
#' @param seed Seed for the subsampling and k-means.
#' @return A list with `assignments`, `profiles` (metric means per cluster),
#'   `diagnostics`, `z`, `centers`, `hclust`, `pca`, `k` and `metrics`.
#' @export
cluster_iz_heterogeneity <- function(iz_summary, metrics, k_range = 2:8,
                                     k = NULL, nstart = 50L,
                                     silhouette_n = 2000L, seed = 1L) {
  rlang::check_installed("cluster", reason = "to compute silhouette widths.")
  z <- .metric_matrix(iz_summary, metrics)
  rng <- withr::with_seed(seed, {
    take <- if (nrow(z) > silhouette_n) sort(sample(seq_len(nrow(z)), silhouette_n))
            else seq_len(nrow(z))
    list(take = take,
         fits = lapply(k_range, function(kk) stats::kmeans(z, kk, nstart = nstart)))
  })
  take <- rng$take
  fits <- rng$fits
  ds <- stats::dist(z[take, , drop = FALSE])
  diag <- data.frame(
    k = k_range,
    total_withinss = vapply(fits, function(f) f$tot.withinss, numeric(1)),
    silhouette = vapply(fits, function(f) {
      lab <- .nearest_center(z[take, , drop = FALSE], f$centers)
      mean(cluster::silhouette(lab, ds)[, "sil_width"])
    }, numeric(1))
  )
  if (is.null(k)) k <- diag$k[which.max(diag$silhouette)]
  fit <- fits[[match(k, k_range)]]
  hc <- stats::hclust(stats::dist(z), method = "ward.D2")
  pc <- stats::prcomp(z, center = FALSE, scale. = FALSE)
  assignment <- data.frame(
    IZ_id = iz_summary$IZ_id,
    cluster_kmeans = paste0("cluster", fit$cluster),
    cluster_hierarchical = paste0("cluster", stats::cutree(hc, k)),
    PC1 = pc$x[, 1], PC2 = pc$x[, 2]
  )
  # mean of each metric per k-means cluster (clusters in C-locale order)
  cl <- assignment$cluster_kmeans[match(iz_summary$IZ_id, assignment$IZ_id)]
  lev <- sort(unique(cl), method = "radix")
  profile <- data.frame(cluster_kmeans = lev, stringsAsFactors = FALSE)
  for (m in metrics) {
    profile[[m]] <- vapply(lev, function(g) mean(iz_summary[[m]][cl == g], na.rm = TRUE),
                           numeric(1), USE.NAMES = FALSE)
  }
  profile$n_IZ <- vapply(lev, function(g) sum(cl == g), integer(1), USE.NAMES = FALSE)
  list(assignments = assignment, profiles = profile, diagnostics = diag,
       z = z, centers = fit$centers, hclust = hc, pca = pc, k = k,
       metrics = metrics)
}

#' IZ pairs on the same chromosome within max_distance
#' @keywords internal
#' @noRd
.nearby_pairs <- function(iz_metadata, max_distance = 1e6) {
  ans <- vector("list", length(unique(iz_metadata$chr)))
  chrs <- unique(iz_metadata$chr)
  for (cc in seq_along(chrs)) {
    d <- iz_metadata[iz_metadata$chr == chrs[cc], , drop = FALSE]
    d <- d[order(d$center), , drop = FALSE]
    out <- vector("list", nrow(d))
    for (i in seq_len(max(0L, nrow(d) - 1L))) {
      jj <- which(d$center > d$center[i] &
                    d$center - d$center[i] <= max_distance)
      if (length(jj)) out[[i]] <- data.frame(
        chr = chrs[cc], IZ1 = d$IZ_id[i], IZ2 = d$IZ_id[jj],
        distance = d$center[jj] - d$center[i]
      )
    }
    ans[[cc]] <- .bind_rows_df(out)
  }
  .bind_rows_df(ans)
}

#' Co-usage of nearby initiation zones
#'
#' For every pair of IZs on the same chromosome with centres at most
#' `max_distance` apart, correlate their usage across cells after removing the
#' replication-progress trend ([correct_signal_matrix()]); pairs with fewer
#' than `min_complete` complete cells or constant residuals get `NA`.
#' P-values (t-test of the correlation) are BH-adjusted; significant pairs are
#' `"co_usage"` (positive) or `"competition"` (negative), others `"neutral"`.
#'
#' @param usage Cells x IZs usage matrix with IZ ids as column names.
#' @param progress Replication progress of the cells.
#' @param iz_metadata Table with `IZ_id`, `chr` and `center`.
#' @param max_distance Maximum centre distance (bp).
#' @param correction,n_progress_bins,spline_df Trend-correction settings.
#' @param min_complete Minimum number of complete cells.
#' @param fdr_threshold FDR threshold.
#' @param BPPARAM A `BiocParallel::BiocParallelParam` used over chromosomes.
#' @return A `data.frame` of pairs: `chr`, `IZ1`, `IZ2`, `distance`,
#'   `association`, `p_value`, `n_cells`, `FDR`, `association_type` and
#'   `distance_bin`.
#' @export
compute_nearby_iz_associations <- function(usage, progress, iz_metadata,
                                           max_distance = 1e6,
                                           correction = c("bin", "spline"),
                                           n_progress_bins = 10L,
                                           spline_df = 4L,
                                           min_complete = 30L,
                                           fdr_threshold = 0.05,
                                           BPPARAM = BiocParallel::SerialParam()) {
  correction <- match.arg(correction)
  pairs <- .nearby_pairs(iz_metadata, max_distance)
  split_pairs <- split(pairs, pairs$chr)
  one_chr <- function(p) {
    ids <- sort(unique(c(p$IZ1, p$IZ2)))
    j <- match(as.character(ids), colnames(usage))
    x <- as.matrix(usage[, j, drop = FALSE])
    r <- correct_signal_matrix(x, progress, correction,
                               n_progress_bins, spline_df)$residual
    colnames(r) <- as.character(ids)
    stat <- lapply(seq_len(nrow(p)), function(i) {
      a <- r[, as.character(p$IZ1[i])]
      b <- r[, as.character(p$IZ2[i])]
      ok <- is.finite(a) & is.finite(b)
      n <- sum(ok)
      rho <- if (n >= min_complete && stats::sd(a[ok]) > 0 &&
                   stats::sd(b[ok]) > 0) stats::cor(a[ok], b[ok]) else NA_real_
      tval <- if (is.finite(rho) && abs(rho) < 1)
        rho * sqrt((n - 2) / (1 - rho^2)) else NA_real_
      data.frame(association = rho, p_value = if (is.finite(tval))
        2 * stats::pt(-abs(tval), df = n - 2) else NA_real_, n_cells = n)
    })
    cbind(p, .bind_rows_df(stat))
  }
  edges <- .bind_rows_df(BiocParallel::bplapply(split_pairs, one_chr,
                                                   BPPARAM = BPPARAM))
  edges$FDR <- stats::p.adjust(edges$p_value, method = "BH")
  edges$association_type <- ifelse(
    is.finite(edges$FDR) & edges$FDR <= fdr_threshold,
    ifelse(edges$association > 0, "co_usage", "competition"), "neutral"
  )
  edges$distance_bin <- cut(
    edges$distance, breaks = c(0, 1e5, 2.5e5, 5e5, 1e6, Inf),
    include.lowest = TRUE, right = TRUE,
    labels = c("<100kb", "100-250kb", "250-500kb", "500kb-1Mb", ">1Mb")
  )
  edges
}

#' Column-wise Khatri-Rao product
#' @keywords internal
#' @noRd
.khatri_rao <- function(x, y) {
  do.call(cbind, lapply(seq_len(ncol(x)), function(r)
    kronecker(x[, r], y[, r])))
}

#' Relative reconstruction error of a CP model
#' @keywords internal
#' @noRd
.cp_error <- function(y, a, b, cc) {
  y1 <- matrix(y, nrow = dim(y)[1])
  inner <- sum(a * (y1 %*% .khatri_rao(cc, b)))
  model_ss <- sum(crossprod(a) * crossprod(b) * crossprod(cc))
  sse <- max(0, sum(y^2) - 2 * inner + model_ss)
  sqrt(sse / max(sum(y^2), .Machine$double.eps))
}

#' Non-negative CP decomposition of a three-way tensor
#'
#' Multiplicative-update non-negative CANDECOMP/PARAFAC with random uniform
#' initialisation. The relative error is checked every `check_every`
#' iterations; the fit stops when it changes by less than `tol`.
#'
#' @param y Non-negative array with three dimensions and no `NA`.
#' @param rank Number of components.
#' @param seed Seed of the initialisation.
#' @param max_iter Maximum number of iterations.
#' @param tol Convergence tolerance.
#' @param check_every Error-check interval.
#' @param epsilon Lower bound of denominators.
#' @return A list with unit-norm factor matrices `A`, `B`, `C`, `lambda`,
#'   `error` (relative), `explained`, `iterations`, `converged`, `history` and
#'   `rank`.
#' @export
fit_nonnegative_cp <- function(y, rank, seed = 1L, max_iter = 50L,
                               tol = 1e-4, check_every = 5L,
                               epsilon = 1e-10) {
  if (any(y < 0, na.rm = TRUE) || anyNA(y)) stop("y must be finite and non-negative.")
  dd <- dim(y)
  init <- withr::with_seed(seed, list(
    a = matrix(stats::runif(dd[1] * rank) + 0.1, dd[1], rank),
    b = matrix(stats::runif(dd[2] * rank) + 0.1, dd[2], rank),
    cc = matrix(stats::runif(dd[3] * rank) + 0.1, dd[3], rank)
  ))
  a <- init$a
  b <- init$b
  cc <- init$cc
  y1 <- matrix(y, nrow = dd[1])
  y2 <- matrix(aperm(y, c(2, 1, 3)), nrow = dd[2])
  y3 <- matrix(aperm(y, c(3, 1, 2)), nrow = dd[3])
  history <- numeric()
  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    a <- a * (y1 %*% .khatri_rao(cc, b)) /
      pmax(a %*% (crossprod(b) * crossprod(cc)), epsilon)
    b <- b * (y2 %*% .khatri_rao(cc, a)) /
      pmax(b %*% (crossprod(a) * crossprod(cc)), epsilon)
    cc <- cc * (y3 %*% .khatri_rao(b, a)) /
      pmax(cc %*% (crossprod(a) * crossprod(b)), epsilon)
    sb <- pmax(sqrt(colSums(b^2)), epsilon)
    b <- sweep(b, 2, sb, "/"); a <- sweep(a, 2, sb, "*")
    sc <- pmax(sqrt(colSums(cc^2)), epsilon)
    cc <- sweep(cc, 2, sc, "/"); a <- sweep(a, 2, sc, "*")
    if (iter %% check_every == 0L || iter == max_iter) {
      err <- .cp_error(y, a, b, cc)
      history <- c(history, err)
      if (length(history) > 1L &&
          abs(history[length(history) - 1L] - err) < tol * max(1, history[length(history) - 1L])) {
        converged <- TRUE
        break
      }
    }
  }
  sa <- pmax(sqrt(colSums(a^2)), epsilon)
  sb <- pmax(sqrt(colSums(b^2)), epsilon)
  sc <- pmax(sqrt(colSums(cc^2)), epsilon)
  lambda <- sa * sb * sc
  a <- sweep(a, 2, sa, "/")
  b <- sweep(b, 2, sb, "/")
  cc <- sweep(cc, 2, sc, "/")
  list(A = a, B = b, C = cc, lambda = lambda,
       error = utils::tail(history, 1), explained = 1 - utils::tail(history, 1)^2,
       iterations = iter, converged = converged, history = history, rank = rank)
}

#' Mean best-match similarity of the B factors across runs
#' @keywords internal
#' @noRd
.factor_stability <- function(fits) {
  if (length(fits) < 2L) return(NA_real_)
  cmb <- utils::combn(seq_along(fits), 2)
  mean(apply(cmb, 2, function(ii) {
    s <- abs(crossprod(fits[[ii[1]]]$B, fits[[ii[2]]]$B))
    mean(c(apply(s, 1, max), apply(s, 2, max)))
  }))
}

#' Choose the rank of a non-negative CP decomposition
#'
#' Fit [fit_nonnegative_cp()] `n_init` times for every rank (seed
#' `seed + 1000 * rank + run`). Among ranks whose factor stability reaches
#' `stability_threshold` (all ranks if none does), the smallest rank whose
#' mean error is within 10% of the error range from the minimum is chosen.
#'
#' @param y Non-negative three-way array.
#' @param ranks Candidate ranks.
#' @param n_init Number of initialisations per rank.
#' @param seed Base seed.
#' @param max_iter,tol Passed to [fit_nonnegative_cp()].
#' @param stability_threshold Minimum stability.
#' @param BPPARAM A `BiocParallel::BiocParallelParam` used over fits.
#' @return A list with `summary`, `chosen_rank`, `best_fits` (per rank),
#'   `all_fits` and `jobs`.
#' @export
select_cp_rank <- function(y, ranks = 2:6, n_init = 2L, seed = 1L,
                           max_iter = 30L, tol = 1e-4,
                           stability_threshold = 0.80,
                           BPPARAM = BiocParallel::SerialParam()) {
  jobs <- expand.grid(rank = ranks, run = seq_len(n_init))
  fits <- BiocParallel::bplapply(seq_len(nrow(jobs)), function(i)
    fit_nonnegative_cp(y, jobs$rank[i], seed + 1000L * jobs$rank[i] + jobs$run[i],
                       max_iter, tol), BPPARAM = BPPARAM)
  by_rank <- split(seq_along(fits), jobs$rank)
  summary <- .bind_rows_df(lapply(names(by_rank), function(rr) {
    ff <- fits[by_rank[[rr]]]
    data.frame(rank = as.integer(rr),
               reconstruction_error = mean(vapply(ff, `[[`, numeric(1), "error")),
               best_error = min(vapply(ff, `[[`, numeric(1), "error")),
               explained_variance = max(vapply(ff, `[[`, numeric(1), "explained")),
               stability = .factor_stability(ff))
  }))
  stable <- summary$rank[is.finite(summary$stability) &
                           summary$stability >= stability_threshold]
  pool <- if (length(stable)) summary[summary$rank %in% stable, ] else summary
  cutoff <- min(pool$reconstruction_error) +
    0.10 * diff(range(pool$reconstruction_error))
  chosen <- min(pool$rank[pool$reconstruction_error <= cutoff])
  best <- lapply(by_rank, function(ii) fits[[ii[which.min(vapply(fits[ii], `[[`,
                                                                  numeric(1), "error"))]]])
  list(summary = summary, chosen_rank = chosen, best_fits = best,
       all_fits = fits, jobs = jobs)
}

#' Split a signed residual tensor into non-negative parts
#'
#' @param residual Array, cells x IZs x positions (e.g.
#'   `correct_replication_tensor()$residual`).
#' @return An array with the third dimension doubled: positive parts followed
#'   by negative parts.
#' @export
split_signed_residual <- function(residual) {
  out <- array(c(pmax(residual, 0), pmax(-residual, 0)),
               dim = c(dim(residual)[1:2], 2L * dim(residual)[3]))
  dn <- dimnames(residual)
  dimnames(out) <- list(dn[[1]], dn[[2]],
                        c(paste0("positive_", dn[[3]]),
                          paste0("negative_", dn[[3]])))
  out
}

#' Sample cells evenly across replication progress
#'
#' @param cell_metadata Table with `cell_id` and `replication_progress`.
#' @param n_cells Number of cells to keep (all cells if there are fewer).
#' @param n_bins Number of progress bins sampled from.
#' @param seed Seed.
#' @return Character vector of cell ids.
#' @export
sample_cells_by_progress <- function(cell_metadata, n_cells = 250L,
                                     n_bins = 10L, seed = 1L) {
  if (nrow(cell_metadata) <= n_cells) return(cell_metadata$cell_id)
  pb <- make_progress_bins(cell_metadata$replication_progress, n_bins)
  per <- ceiling(n_cells / n_bins)
  ids <- withr::with_seed(seed, {
    ids <- unlist(lapply(levels(pb), function(b) {
      x <- which(pb == b)
      sample(x, min(per, length(x)))
    }))
    if (length(ids) > n_cells) ids <- sort(sample(ids, n_cells))
    ids
  })
  cell_metadata$cell_id[ids]
}

#' Select informative initiation zones for the tensor
#'
#' Score IZs by their mean absolute scaled metric plus relative
#' `log1p(n_active_cells)`, keep those with at least 10 active cells, and take
#' the top IZs per k-means cluster, then the top `n_iz` overall.
#'
#' @param iz_summary Per-IZ table with the metrics and `n_active_cells`.
#' @param assignments `cluster_iz_heterogeneity()$assignments`.
#' @param metrics Metric column names.
#' @param n_iz Number of IZs.
#' @return Integer vector of `IZ_id`.
#' @export
select_iz_for_tensor <- function(iz_summary, assignments, metrics,
                                 n_iz = 80L) {
  d <- .left_join_by_id(as.data.frame(iz_summary),
                        assignments[, c("IZ_id", "cluster_kmeans")], "IZ_id")
  d$tensor_score <- rowMeans(abs(.metric_matrix(d, metrics)), na.rm = TRUE) +
    log1p(d$n_active_cells) / max(log1p(d$n_active_cells), na.rm = TRUE)
  per <- ceiling(n_iz / length(unique(d$cluster_kmeans)))
  d <- d[!is.na(d$n_active_cells) & d$n_active_cells >= 10 & is.finite(d$tensor_score), ,
         drop = FALSE]
  # top `per` IZs per cluster (clusters in C-locale order, NA last), then the
  # top n_iz overall; ties keep the earlier row
  groups <- sort(unique(d$cluster_kmeans), method = "radix", na.last = TRUE)
  rows <- unlist(lapply(groups, function(g) {
    ii <- if (is.na(g)) which(is.na(d$cluster_kmeans)) else which(d$cluster_kmeans %in% g)
    utils::head(ii[order(-d$tensor_score[ii])], per)
  }))
  rows <- utils::head(rows[order(-d$tensor_score[rows])], n_iz)
  ids <- d$IZ_id[rows]
  as.integer(ids)
}

#' Tables of CP factors
#'
#' @param fit Output of [fit_nonnegative_cp()].
#' @param cell_ids,iz_ids,relative_bp Labels of the three tensor modes.
#' @param cell_metadata Table with `cell_id` and `replication_progress`.
#' @param iz_summary Per-IZ table with `IZ_id`.
#' @param signed_track Logical; the tensor came from
#'   [split_signed_residual()] (track factors = positive minus negative part).
#' @return A list with `cell`, `iz`, `track`, `progress` (Spearman correlation
#'   of each cell factor with progress) and `top_iz` (top 5% IZs per factor).
#' @export
cp_factor_tables <- function(fit, cell_ids, iz_ids, relative_bp,
                             cell_metadata, iz_summary,
                             signed_track = TRUE) {
  factors <- paste0("factor", seq_len(fit$rank))
  a <- sweep(fit$A, 2, fit$lambda, "*")
  cell <- data.frame(cell_id = cell_ids, a, check.names = FALSE)
  names(cell)[-(1)] <- factors
  cell <- .left_join_by_id(cell, as.data.frame(cell_metadata), "cell_id")
  iz <- data.frame(IZ_id = iz_ids, fit$B, check.names = FALSE)
  names(iz)[-(1)] <- factors
  iz <- .left_join_by_id(iz, as.data.frame(iz_summary), "IZ_id")
  if (signed_track) {
    k <- length(relative_bp)
    tr <- fit$C[seq_len(k), , drop = FALSE] -
      fit$C[k + seq_len(k), , drop = FALSE]
  } else tr <- fit$C
  track <- data.frame(relative_bp = relative_bp, tr, check.names = FALSE)
  names(track)[-(1)] <- factors
  progress <- data.frame(
    factor = factors,
    progress_rho = vapply(factors, function(f)
      stats::cor(cell[[f]], cell$replication_progress,
                 method = "spearman", use = "complete.obs"), numeric(1))
  )
  top_n <- max(1L, ceiling(0.05 * nrow(iz)))
  top_iz <- .bind_rows_df(lapply(factors, function(f) {
    top <- iz[utils::head(order(-iz[[f]]), top_n), , drop = FALSE]
    top$factor <- f
    top$loading <- top[[f]]
    top
  }))
  list(cell = cell, iz = iz, track = track, progress = progress,
       top_iz = top_iz)
}

#' Mean bigWig signal around IZ centres
#'
#' @param iz_metadata Table with `IZ_id`, `chr` and `center`.
#' @param features Named character vector of bigWig files.
#' @param window_bp Window width centred on the IZ centre.
#' @return A `data.frame` with `IZ_id` and one column per feature.
#' @export
annotate_iz_bigwig <- function(iz_metadata, features, window_bp = 100000L) {
  rlang::check_installed("rtracklayer", reason = "to read bigWig files.")
  stopifnot(all(file.exists(features)), !is.null(names(features)))
  gr <- GenomicRanges::GRanges(
    iz_metadata$chr,
    IRanges::IRanges(start = pmax(1, round(iz_metadata$center)), width = 1)
  )
  gr <- GenomicRanges::trim(GenomicRanges::resize(gr, window_bp, fix = "center"))
  out <- data.frame(IZ_id = iz_metadata$IZ_id)
  for (nm in names(features)) {
    message("summarize ", nm)
    z <- rtracklayer::summary(rtracklayer::BigWigFile(features[[nm]]), gr,
                              size = 1L, type = "mean", as = "matrix")
    out[[nm]] <- as.numeric(z[, 1])
  }
  out
}

#' Correlate heterogeneity metrics with epigenomic features
#'
#' @param d Per-IZ table with metric and feature columns.
#' @param metrics,features Column names.
#' @param methods Correlation methods.
#' @return A `data.frame` with one row per metric, feature and method
#'   (`correlation`, `p_value`, `n`, BH `FDR`; `NA` with fewer than 10 values).
#' @export
correlate_epigenomic_features <- function(d, metrics, features,
                                          methods = c("spearman", "pearson")) {
  ans <- expand.grid(metric = metrics, feature = features, method = methods,
                     stringsAsFactors = FALSE)
  stat <- lapply(seq_len(nrow(ans)), function(i) {
    x <- d[[ans$metric[i]]]; y <- d[[ans$feature[i]]]
    ok <- is.finite(x) & is.finite(y)
    z <- if (sum(ok) >= 10L) suppressWarnings(stats::cor.test(
      x[ok], y[ok], method = ans$method[i], exact = FALSE)) else NULL
    data.frame(correlation = if (is.null(z)) NA_real_ else unname(z$estimate),
               p_value = if (is.null(z)) NA_real_ else z$p.value,
               n = sum(ok))
  })
  ans <- cbind(ans, .bind_rows_df(stat))
  ans$FDR <- stats::p.adjust(ans$p_value, "BH")
  ans
}

#' sign(x) * log1p(|x|)
#' @keywords internal
#' @noRd
.signed_log1p <- function(x) sign(x) * log1p(abs(x))

#' Linear models of heterogeneity metrics on epigenomic features
#'
#' For each metric, fit `scale(metric) ~ features + covariates` with features
#' transformed by signed `log1p` and scaled, and `width` (log1p), `S50` and
#' `mean_usage` scaled.
#'
#' @param d Per-IZ table with the metrics, features, `width`, `S50`,
#'   `mean_usage` and the covariates.
#' @param metrics,features Column names.
#' @param covariates Covariate column names.
#' @return A `data.frame` of feature coefficients per metric with BH `FDR`.
#' @export
fit_epigenomic_models <- function(d, metrics, features,
                                  covariates = c("width", "S50", "mean_usage", "chr")) {
  z <- d
  for (f in features) z[[f]] <- as.numeric(scale(.signed_log1p(z[[f]])))
  z$width <- as.numeric(scale(log1p(z$width)))
  z$S50 <- as.numeric(scale(z$S50))
  z$mean_usage <- as.numeric(scale(z$mean_usage))
  rhs <- paste(c(features, covariates), collapse = " + ")
  ans <- lapply(metrics, function(metric) {
    fit <- stats::lm(stats::as.formula(paste0("scale(", metric, ") ~ ", rhs)),
                     data = z, na.action = stats::na.omit)
    co <- summary(fit)$coefficients
    keep <- intersect(features, rownames(co))
    data.frame(metric = metric, feature = keep, beta = co[keep, 1],
               std_error = co[keep, 2], p_value = co[keep, 4],
               model_r2 = summary(fit)$adj.r.squared,
               n = stats::nobs(fit))
  })
  ans <- .bind_rows_df(ans)
  ans$FDR <- stats::p.adjust(ans$p_value, "BH")
  ans
}

#' Logistic classifier of high versus low heterogeneity
#'
#' IZs in the bottom / top quartile of the composite heterogeneity (mean of
#' the scaled metrics) are labelled 0 / 1 and modelled by logistic regression
#' on the transformed features plus `width`, `S50`, `mean_usage` and `chr`.
#'
#' @inheritParams fit_epigenomic_models
#' @return A list with `model`, `coefficients` (with odds ratios, training
#'   AUC and BH `FDR`) and `labeled_data`.
#' @export
fit_epigenomic_classifier <- function(d, metrics, features) {
  z <- d
  m <- .metric_matrix(z, metrics)
  z$heterogeneity_composite <- rowMeans(m)
  q <- stats::quantile(z$heterogeneity_composite, c(0.25, 0.75), na.rm = TRUE)
  z$class <- ifelse(z$heterogeneity_composite <= q[1], 0L,
                    ifelse(z$heterogeneity_composite >= q[2], 1L, NA_integer_))
  for (f in features) z[[f]] <- as.numeric(scale(.signed_log1p(z[[f]])))
  z$width <- as.numeric(scale(log1p(z$width)))
  z$S50 <- as.numeric(scale(z$S50))
  z$mean_usage <- as.numeric(scale(z$mean_usage))
  fit <- stats::glm(stats::as.formula(paste(
    "class ~", paste(c(features, "width", "S50", "mean_usage", "chr"), collapse = " + "))),
    data = z, family = stats::binomial(), na.action = stats::na.omit)
  co <- summary(fit)$coefficients
  keep <- intersect(features, rownames(co))
  pred <- stats::fitted(fit)
  yy <- stats::model.response(stats::model.frame(fit))
  auc <- if (length(unique(yy)) == 2L) {
    r <- rank(pred)
    n1 <- sum(yy == 1); n0 <- sum(yy == 0)
    (sum(r[yy == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  } else NA_real_
  table <- data.frame(feature = keep, beta = co[keep, 1],
                      odds_ratio = exp(co[keep, 1]),
                      std_error = co[keep, 2], p_value = co[keep, 4],
                      n = stats::nobs(fit), training_auc = auc)
  table$FDR <- stats::p.adjust(table$p_value, "BH")
  list(model = fit, coefficients = table, labeled_data = z)
}
