# Single-cell initiation-zone heterogeneity (phases 1-3).
# Migrated from jupyter/iz_heterogeneity_phase1_3.R together with
# jupyter/iz_heterogeneity_state_usage.R, which is sourced after it and
# replaces compute_usage_heterogeneity() by an NA-aware version and wraps
# fit_timing_models() (both merged here). The dplyr steps are rewritten in
# base R. tensor_heatmap_data() (plot preparation) is not included.

#' Stop if columns are missing
#' @keywords internal
#' @noRd
.need_cols <- function(x, cols, label) {
  miss <- setdiff(cols, names(x))
  if (length(miss)) stop(label, " is missing: ", paste(miss, collapse = ", "))
}

#' Replication-progress bins
#'
#' Clamp `progress` to `[0, 1]` and cut it into `n_bins` equal bins labelled
#' `P01`, `P02`, ...
#'
#' @param progress Numeric replication progress (fraction).
#' @param n_bins Number of bins.
#' @return A factor.
#' @export
make_progress_bins <- function(progress, n_bins = 10L) {
  p <- pmin(pmax(as.numeric(progress), 0), 1)
  cut(p, breaks = seq(0, 1, length.out = n_bins + 1L),
      include.lowest = TRUE, labels = sprintf("P%02d", seq_len(n_bins)))
}

#' Prepare track, cell and IZ tables for heterogeneity analysis
#'
#' Collect the IZ-assigned tracks (one row per cell and IZ; duplicated
#' cell-IZ assignments are merged: outer start/end, width-weighted centre,
#' maximum flag), the IZ metadata and the cell replication progress, and
#' derive per-track features relative to the IZ centre (`q50`): relative
#' initiation position, left/right/total track length and fork asymmetry.
#'
#' @param tracks Named `GRangesList` of tracks (every cell with at least one
#'   track) with columns `cell_id`, `rep_percent`, `iz_id` and `flag` (e.g.
#'   `define_initiation_zones()$tracks`).
#' @param gr_izones `GRanges` of IZs with `iz_id`, `q5`-`q95`, `n_tracks` and
#'   optionally the firing range (e.g. [calculate_firing_range()]).
#' @param cell_metadata `NULL` or a `data.frame` with `cell_id` and
#'   `progress_col`; if it has `NA`s, the track `rep_percent` is used.
#' @param progress_col Progress column of `cell_metadata`.
#'
#' @return A list with `track_intervals`, `cell_metadata` (with
#'   `replication_progress`), `iz_metadata` and `diagnostics`.
#' @export
prepare_iz_data <- function(tracks, gr_izones, cell_metadata = NULL,
                            progress_col = "rep_percent") {
  stopifnot(inherits(tracks, "GRangesList"), inherits(gr_izones, "GRanges"))
  gr <- BiocGenerics::unlist(tracks, use.names = FALSE)
  mc <- S4Vectors::mcols(gr)
  .need_cols(mc, c("cell_id", "rep_percent", "iz_id", "flag"), "track metadata")

  keep <- !is.na(mc$iz_id)
  g <- gr[keep]
  gm <- S4Vectors::mcols(g)
  raw <- data.frame(
    cell_id = as.character(gm$cell_id), IZ_id = as.integer(gm$iz_id),
    chr = as.character(GenomeInfoDb::seqnames(g)),
    start = BiocGenerics::start(g), end = BiocGenerics::end(g),
    segment_width = BiocGenerics::width(g), flag = as.integer(gm$flag),
    replication_progress = as.numeric(gm$rep_percent)
  )
  raw$segment_center <- (raw$start + raw$end) / 2

  # Merge the few duplicated cell-IZ assignments; the weighted center is the
  # center of mass of their replicated intervals.
  track_intervals <- .merge_cell_iz_tracks(raw)

  im <- S4Vectors::mcols(gr_izones)
  .need_cols(im, c("iz_id", "q5", "q25", "q50", "q75", "q95", "n_tracks"),
             "IZ metadata")
  iz_metadata <- data.frame(
    IZ_id = as.integer(im$iz_id),
    chr = as.character(GenomeInfoDb::seqnames(gr_izones)),
    start = BiocGenerics::start(gr_izones), end = BiocGenerics::end(gr_izones),
    center = as.numeric(im$q50), width = BiocGenerics::width(gr_izones),
    q5 = as.numeric(im$q5), q25 = as.numeric(im$q25), q50 = as.numeric(im$q50),
    q75 = as.numeric(im$q75), q95 = as.numeric(im$q95),
    n_tracks = as.integer(im$n_tracks),
    earliest_rep_percent = as.numeric(im$earliest_rep_percent),
    latest_unrep_percent = as.numeric(im$latest_unrep_percent),
    corroborated = as.logical(im$corroborated)
  )

  lens <- lengths(tracks)
  if (any(lens == 0L)) stop("Every cell must contain at least one track.")
  first <- cumsum(c(1L, utils::head(lens, -1L)))
  track_cells <- data.frame(
    cell_id = names(tracks),
    replication_progress = as.numeric(mc$rep_percent[first])
  )
  if (is.null(cell_metadata)) {
    cell_metadata <- track_cells
  } else {
    .need_cols(cell_metadata, c("cell_id", progress_col), "cell metadata")
    cell_metadata <- as.data.frame(cell_metadata)
    cell_metadata <- cell_metadata[match(track_cells$cell_id, cell_metadata$cell_id), , drop = FALSE]
    if (anyNA(cell_metadata$cell_id)) stop("Some track cells are absent from cell metadata.")
    cell_metadata$replication_progress <- as.numeric(cell_metadata[[progress_col]])
    if (anyNA(cell_metadata$replication_progress)) {
      cell_metadata$replication_progress <- track_cells$replication_progress
    }
  }
  if (anyNA(cell_metadata$replication_progress)) stop("replication_progress contains NA.")

  z <- match(track_intervals$IZ_id, iz_metadata$IZ_id)
  track_intervals$IZ_center <- iz_metadata$center[z]
  track_intervals$IZ_width <- iz_metadata$width[z]
  track_intervals$relative_init_position <-
    (track_intervals$track_center_bp - track_intervals$IZ_center) /
    pmax(track_intervals$IZ_width / 2, 1)
  track_intervals$left_track_length <-
    pmax(track_intervals$IZ_center - track_intervals$start, 0)
  track_intervals$right_track_length <-
    pmax(track_intervals$end - track_intervals$IZ_center, 0)
  track_intervals$total_track_length <-
    track_intervals$left_track_length + track_intervals$right_track_length
  den <- track_intervals$total_track_length
  track_intervals$fork_asymmetry <- ifelse(
    den > 0,
    (track_intervals$right_track_length - track_intervals$left_track_length) / den,
    NA_real_
  )

  list(
    track_intervals = track_intervals,
    cell_metadata = cell_metadata,
    iz_metadata = iz_metadata,
    diagnostics = data.frame(
      n_cells = nrow(cell_metadata), n_IZ = nrow(iz_metadata),
      n_assigned_tracks = sum(keep), n_active_cell_IZ = nrow(track_intervals),
      duplicated_cell_IZ_merged = sum(track_intervals$n_segments - 1L),
      bin_size_bp = 20000L
    )
  )
}

#' Merge duplicated cell-IZ track assignments (base-R version of the
#' original dplyr group_by / summarise)
#' @keywords internal
#' @noRd
.merge_cell_iz_tracks <- function(raw) {
  keys <- c("cell_id", "IZ_id", "chr", "replication_progress")
  ord <- do.call(order, c(unname(as.list(raw[keys])), list(method = "radix")))
  raw <- raw[ord, , drop = FALSE]
  grp <- do.call(paste, c(unname(as.list(raw[keys])), sep = "\r"))
  first <- !duplicated(grp)
  idx <- split(seq_len(nrow(raw)), factor(grp, levels = unique(grp)))
  out <- raw[first, keys, drop = FALSE]
  out$start <- vapply(idx, function(i) min(raw$start[i]), raw$start[1], USE.NAMES = FALSE)
  out$end <- vapply(idx, function(i) max(raw$end[i]), raw$end[1], USE.NAMES = FALSE)
  out$track_center_bp <- vapply(idx, function(i) {
    stats::weighted.mean(raw$segment_center[i], raw$segment_width[i])
  }, numeric(1), USE.NAMES = FALSE)
  out$n_segments <- lengths(idx, use.names = FALSE)
  out$flag <- vapply(idx, function(i) max(raw$flag[i], na.rm = TRUE), raw$flag[1],
                     USE.NAMES = FALSE)
  rownames(out) <- NULL
  out
}

#' Cell x IZ usage matrix from assigned tracks
#'
#' @param track_intervals `prepare_iz_data()$track_intervals`.
#' @param cell_ids,iz_ids Row and column axes.
#' @param method `"assigned_track"` (1 where a cell has a track assigned to
#'   the IZ) or `"existing"` (`existing_usage != 0`).
#' @param existing_usage Matrix used with `method = "existing"`.
#' @return A sparse matrix, cells x IZs.
#' @export
call_iz_usage <- function(track_intervals = NULL, cell_ids = NULL, iz_ids = NULL,
                          method = c("assigned_track", "existing"),
                          existing_usage = NULL) {
  method <- match.arg(method)
  if (method == "existing") {
    if (is.null(existing_usage)) stop("Supply existing_usage for method='existing'.")
    return(Matrix::Matrix(existing_usage != 0, sparse = TRUE))
  }
  .need_cols(track_intervals, c("cell_id", "IZ_id"), "track_intervals")
  i <- match(track_intervals$cell_id, cell_ids)
  j <- match(track_intervals$IZ_id, iz_ids)
  if (anyNA(i) || anyNA(j)) stop("Track cell/IZ IDs do not match the requested axes.")
  Matrix::sparseMatrix(
    i = i, j = j, x = 1, dims = c(length(cell_ids), length(iz_ids)),
    dimnames = list(cell_id = cell_ids, IZ_id = as.character(iz_ids))
  )
}

#' Cell x IZ usage matrix from replication states
#'
#' Read the state of the genome bin containing each IZ centre: `1` if
#' replicated, `0` if unreplicated, `NA` otherwise.
#'
#' @param genomic_states Matrix of states, bins x cells, with row names
#'   `"chr-start-end"` of `bin_size` bins (1-based starts).
#' @param cell_ids Cells to use (columns of `genomic_states`).
#' @param iz_metadata `prepare_iz_data()$iz_metadata` (`IZ_id`, `chr`,
#'   `center`).
#' @param replicated_state,unreplicated_state State codes.
#' @param bin_size Bin size (bp).
#' @return An integer matrix, cells x IZs.
#' @export
call_iz_usage_from_states <- function(genomic_states, cell_ids, iz_metadata,
                                      replicated_state = 4L, unreplicated_state = 2L,
                                      bin_size = 20000L) {
  if (is.null(rownames(genomic_states)) || is.null(colnames(genomic_states)))
    stop("genomic_states needs genomic-bin row names and cell-ID column names.")
  cidx <- match(cell_ids, colnames(genomic_states))
  if (anyNA(cidx)) stop("Some requested cells are absent from genomic_states.")
  bin_start <- floor((iz_metadata$center - 1) / bin_size) * bin_size + 1
  keys <- sprintf("%s-%.0f-%.0f", iz_metadata$chr, bin_start,
                  bin_start + bin_size - 1)
  ridx <- match(keys, rownames(genomic_states))
  if (anyNA(ridx)) stop("Some IZ q50 positions are absent from genomic_states.")

  state <- t(genomic_states[ridx, cidx, drop = FALSE])
  usage <- matrix(NA_integer_, nrow(state), ncol(state),
                  dimnames = list(cell_id = cell_ids,
                                  IZ_id = as.character(iz_metadata$IZ_id)))
  usage[state == replicated_state] <- 1L
  usage[state == unreplicated_state] <- 0L
  usage
}

#' Dense cell x IZ x position replication tensor
#'
#' `1` where a cell's merged track covers position `center + rel_bp` of the
#' IZ (`rel_bp` from `-window_bp` to `window_bp` by `bin_size`), else `0`.
#'
#' @param track_intervals `prepare_iz_data()$track_intervals`.
#' @param cell_ids Cells (first dimension).
#' @param iz_metadata `prepare_iz_data()$iz_metadata`.
#' @param iz_ids IZs (second dimension).
#' @param window_bp,bin_size Window half-width and step (bp).
#' @param max_elements Maximum number of tensor elements.
#' @return A numeric array, cells x IZs x positions.
#' @export
build_replication_tensor <- function(track_intervals, cell_ids, iz_metadata, iz_ids,
                                     window_bp = 500000L, bin_size = 20000L,
                                     max_elements = 5e7) {
  iz_ids <- as.integer(iz_ids)
  iz <- iz_metadata[match(iz_ids, iz_metadata$IZ_id), , drop = FALSE]
  if (anyNA(iz$IZ_id)) stop("Unknown IZ_id requested.")
  rel_bp <- seq(-window_bp, window_bp, by = bin_size)
  n_elem <- length(cell_ids) * length(iz_ids) * length(rel_bp)
  if (n_elem > max_elements) {
    stop("Dense tensor would contain ", format(n_elem, big.mark = ","),
         " elements. Select fewer IZs; track_intervals is the full sparse/RLE tensor.")
  }
  x <- array(
    0, dim = c(length(cell_ids), length(iz_ids), length(rel_bp)),
    dimnames = list(cell_id = cell_ids, IZ_id = as.character(iz_ids),
                    relative_bp = as.character(rel_bp))
  )
  d <- track_intervals[track_intervals$IZ_id %in% iz_ids &
                         track_intervals$cell_id %in% cell_ids, , drop = FALSE]
  ci <- match(d$cell_id, cell_ids)
  ji <- match(d$IZ_id, iz_ids)
  centers <- iz$center[ji]
  for (r in seq_len(nrow(d))) {
    k <- which(rel_bp >= d$start[r] - centers[r] & rel_bp <= d$end[r] - centers[r])
    if (length(k)) x[ci[r], ji[r], k] <- 1
  }
  x
}

#' Remove the replication-progress trend from a signal matrix
#'
#' Fit the expected value of each column as a function of replication
#' progress, either as the mean within progress bins (`"bin"`) or with a
#' natural spline of `spline_df` degrees of freedom (`"spline"`, falling back
#' to `"bin"` with fewer than 4 cells or 3 distinct progress values), and
#' return the residuals.
#'
#' @param x Numeric matrix, cells x features.
#' @param progress Replication progress of the cells.
#' @param method `"bin"` or `"spline"`.
#' @param n_progress_bins Number of progress bins.
#' @param spline_df Spline degrees of freedom.
#' @return A list with `residual`, `fitted` and `method`.
#' @export
correct_signal_matrix <- function(x, progress, method = c("bin", "spline"),
                                  n_progress_bins = 10L, spline_df = 4L) {
  method <- match.arg(method)
  x <- as.matrix(x)
  progress <- as.numeric(progress)
  fitted <- matrix(NA_real_, nrow(x), ncol(x), dimnames = dimnames(x))
  okp <- is.finite(progress)

  if (method == "bin") {
    pb <- make_progress_bins(progress, n_progress_bins)
    for (b in levels(pb)) {
      ii <- which(pb == b)
      if (!length(ii)) next
      mu <- colMeans(x[ii, , drop = FALSE], na.rm = TRUE)
      mu[is.nan(mu)] <- NA_real_
      fitted[ii, ] <- matrix(mu, length(ii), ncol(x), byrow = TRUE)
    }
  } else {
    n_unique <- length(unique(progress[okp]))
    if (sum(okp) < 4L || n_unique < 3L) {
      return(correct_signal_matrix(x, progress, "bin", n_progress_bins, spline_df))
    }
    df <- min(as.integer(spline_df), n_unique - 1L, sum(okp) - 2L)
    design <- if (df >= 2L) {
      stats::model.matrix(~ splines::ns(progress[okp], df = df))
    } else cbind(1, progress[okp])
    y <- x[okp, , drop = FALSE]
    if (!anyNA(y)) {
      fitted[okp, ] <- design %*% qr.coef(qr(design), y)
    } else {
      for (j in seq_len(ncol(y))) {
        ok <- is.finite(y[, j])
        if (sum(ok) >= ncol(design))
          fitted[which(okp)[ok], j] <- stats::lm.fit(design[ok, , drop = FALSE], y[ok, j])$fitted.values
      }
    }
  }
  list(residual = x - fitted, fitted = fitted, method = method)
}

#' Remove the replication-progress trend from a replication tensor
#'
#' Apply [correct_signal_matrix()] to the tensor flattened to cells x
#' (IZ, position).
#'
#' @inheritParams correct_signal_matrix
#' @param x Array, cells x IZs x positions (e.g.
#'   [build_replication_tensor()]).
#' @param keep_fitted Logical; also return the fitted tensor.
#' @return A list with `residual`, `fitted` (or `NULL`) and `method`.
#' @export
correct_replication_tensor <- function(x, progress, method = c("bin", "spline"),
                                       n_progress_bins = 10L, spline_df = 4L,
                                       keep_fitted = FALSE) {
  method <- match.arg(method)
  d <- dim(x)
  if (length(d) != 3L || length(progress) != d[1]) stop("x must be cell x IZ x bin.")
  ans <- correct_signal_matrix(matrix(x, nrow = d[1]), progress, method,
                               n_progress_bins, spline_df)
  residual <- array(ans$residual, d, dimnames = dimnames(x))
  list(residual = residual,
       fitted = if (keep_fitted) array(ans$fitted, d, dimnames = dimnames(x)) else NULL,
       method = method)
}

#' Usage heterogeneity of initiation zones
#'
#' Per progress bin, the usage probability `p` of each IZ is the mean usage of
#' the cells in the bin (ignoring `NA`). The heterogeneity is the
#' cell-count-weighted mean of `4 p (1 - p)` over the bins with a finite `p`
#' (weights renormalised over those bins); without `NA` all bins are used.
#'
#' @param usage Cells x IZs usage matrix (0/1/`NA`; base or sparse), e.g.
#'   [call_iz_usage()] or [call_iz_usage_from_states()].
#' @param progress Replication progress of the cells.
#' @param n_progress_bins Number of progress bins.
#' @return A list with `summary` (`IZ_id`, `mean_usage`, `usage_by_progress`,
#'   `usage_heterogeneity`, `n_usage_observed`), `by_progress`,
#'   `probability` (bins x IZs) and `progress_bin`.
#' @export
compute_usage_heterogeneity <- function(usage, progress, n_progress_bins = 10L) {
  pb <- make_progress_bins(progress, n_progress_bins)
  lev <- levels(pb)
  n_by <- tabulate(pb, nbins = length(lev))
  p <- vapply(lev, function(b) {
    ii <- which(pb == b)
    if (length(ii)) .col_means(usage[ii, , drop = FALSE], na.rm = TRUE)
    else rep(NA_real_, ncol(usage))
  }, numeric(ncol(usage)))
  p <- t(p)
  p[is.nan(p)] <- NA_real_
  rownames(p) <- lev
  colnames(p) <- colnames(usage)
  w <- n_by / sum(n_by)
  h <- vapply(seq_len(ncol(p)), function(j) {
    ok <- is.finite(p[, j])
    if (any(ok)) stats::weighted.mean(4 * p[ok, j] * (1 - p[ok, j]), w[ok])
    else NA_real_
  }, numeric(1))
  mean_usage <- .col_means(usage, na.rm = TRUE)
  mean_usage[is.nan(mean_usage)] <- NA_real_
  p_text <- vapply(seq_len(ncol(p)), function(j)
    paste(sprintf("%s=%.4f", lev, p[, j]), collapse = ";"), character(1))
  summary <- data.frame(
    IZ_id = as.integer(colnames(usage)), mean_usage = mean_usage,
    usage_by_progress = p_text, usage_heterogeneity = h,
    n_usage_observed = .col_sums(is.finite(usage))
  )
  mids <- (seq_len(length(lev)) - 0.5) / length(lev)
  by_progress <- data.frame(
    IZ_id = rep(as.integer(colnames(usage)), each = length(lev)),
    progress_bin = rep(lev, ncol(usage)), progress_mid = rep(mids, ncol(usage)),
    n_cells = rep(n_by, ncol(usage)), usage_probability = as.vector(p)
  )
  list(summary = summary, by_progress = by_progress, probability = p,
       progress_bin = pb)
}

#' Column means / sums for base or Matrix objects
#' @keywords internal
#' @noRd
.col_means <- function(x, na.rm = FALSE) {
  if (methods::is(x, "Matrix")) Matrix::colMeans(x, na.rm = na.rm) else colMeans(x, na.rm = na.rm)
}

#' @keywords internal
#' @noRd
.col_sums <- function(x, na.rm = FALSE) {
  if (methods::is(x, "Matrix")) Matrix::colSums(x, na.rm = na.rm) else colSums(x, na.rm = na.rm)
}

#' Row-bind data frames (dplyr::bind_rows for identical columns)
#' @keywords internal
#' @noRd
.bind_rows_df <- function(x) {
  x <- Filter(Negate(is.null), x)
  if (!length(x)) return(data.frame())
  out <- do.call(rbind, unname(x))
  rownames(out) <- NULL
  out
}

#' Logistic timing model of each initiation zone
#'
#' For each IZ, fit `usage ~ progress` by logistic regression
#' (`stats::glm.fit`) on cells with finite values and report the intercept,
#' `S50 = -intercept / slope` (progress at 50% usage), slope,
#' `timing_heterogeneity = 1 / (1 + |slope|)` and the deviance-explained
#' quality. Flags: `always_0` / `always_1`, `no_observations` (no finite
#' usage), `fit_failed`, `non_converged`, `flat_slope`, `decreasing`,
#' `S50_out_of_range`, `separation` (|coefficient| > 30) and `non_monotonic`
#' (a drop larger than `nonmonotonic_tolerance` in `usage_probability`).
#'
#' @param usage Cells x IZs usage matrix.
#' @param progress Replication progress of the cells.
#' @param usage_probability `NULL` or bins x IZs probabilities
#'   (`compute_usage_heterogeneity()$probability`).
#' @param nonmonotonic_tolerance Tolerated drop between consecutive bins.
#' @param BPPARAM A `BiocParallel::BiocParallelParam` used over IZ chunks.
#' @return A `data.frame`, one row per IZ.
#' @export
fit_timing_models <- function(usage, progress, usage_probability = NULL,
                              nonmonotonic_tolerance = 0.10,
                              BPPARAM = BiocParallel::SerialParam()) {
  ans <- .fit_timing_models_impl(usage, progress, usage_probability,
                                 nonmonotonic_tolerance, BPPARAM)
  n_observed <- .col_sums(is.finite(usage))
  ans$timing_flag[n_observed[match(ans$IZ_id, as.integer(colnames(usage)))] == 0L] <-
    "no_observations"
  ans
}

#' Timing models without the no_observations flag (phase1_3 version)
#' @keywords internal
#' @noRd
.fit_timing_models_impl <- function(usage, progress, usage_probability = NULL,
                              nonmonotonic_tolerance = 0.10,
                              BPPARAM = BiocParallel::SerialParam()) {
  progress <- as.numeric(progress)
  one <- function(j) {
    y <- as.numeric(usage[, j])
    ok <- is.finite(progress) & is.finite(y)
    y <- y[ok]; s <- progress[ok]
    flag <- character()
    if (!length(y) || all(y == 0) || all(y == 1)) {
      flag <- if (!length(y) || all(y == 0)) "always_0" else "always_1"
      return(data.frame(IZ_id = as.integer(colnames(usage)[j]), intercept = NA_real_,
                        S50 = NA_real_, slope = NA_real_, timing_heterogeneity = NA_real_,
                        model_quality = NA_real_, n_cells = length(y), timing_flag = flag))
    }
    fit <- try(suppressWarnings(stats::glm.fit(cbind(1, s), y,
                                               family = stats::binomial())), silent = TRUE)
    if (inherits(fit, "try-error")) {
      return(data.frame(IZ_id = as.integer(colnames(usage)[j]), intercept = NA_real_,
                        S50 = NA_real_, slope = NA_real_, timing_heterogeneity = NA_real_,
                        model_quality = NA_real_, n_cells = length(y), timing_flag = "fit_failed"))
    }
    b <- fit$coefficients
    slope <- unname(b[2]); s50 <- if (is.finite(slope) && abs(slope) > 1e-8) -b[1] / slope else NA_real_
    if (!isTRUE(fit$converged)) flag <- c(flag, "non_converged")
    if (!is.finite(slope) || abs(slope) <= 1e-8) flag <- c(flag, "flat_slope")
    if (is.finite(slope) && slope < 0) flag <- c(flag, "decreasing")
    if (is.finite(s50) && (s50 < 0 || s50 > 1)) flag <- c(flag, "S50_out_of_range")
    if (any(abs(b) > 30, na.rm = TRUE)) flag <- c(flag, "separation")
    if (!is.null(usage_probability)) {
      dp <- diff(usage_probability[, j])
      if (any(dp < -nonmonotonic_tolerance, na.rm = TRUE)) flag <- c(flag, "non_monotonic")
    }
    quality <- if (is.finite(fit$null.deviance) && fit$null.deviance > 0)
      1 - fit$deviance / fit$null.deviance else NA_real_
    data.frame(
      IZ_id = as.integer(colnames(usage)[j]), intercept = unname(b[1]), S50 = s50,
      slope = slope,
      timing_heterogeneity = if (is.finite(slope)) 1 / (1 + abs(slope)) else NA_real_,
      model_quality = quality, n_cells = length(y),
      timing_flag = if (length(flag)) paste(unique(flag), collapse = ";") else "ok"
    )
  }
  n_worker <- max(1L, BiocParallel::bpworkers(BPPARAM))
  chunks <- split(seq_len(ncol(usage)),
                  cut(seq_len(ncol(usage)), min(ncol(usage), 4L * n_worker), labels = FALSE))
  ans <- BiocParallel::bplapply(chunks, function(ii) do.call(rbind, lapply(ii, one)),
                                BPPARAM = BPPARAM)
  .bind_rows_df(ans)
}

#' Progress-corrected feature residuals of one IZ
#' @keywords internal
#' @noRd
.correct_feature <- function(y, progress, method, n_progress_bins, spline_df) {
  ok <- is.finite(y) & is.finite(progress)
  out <- rep(NA_real_, length(y))
  if (sum(ok) < 4L) {
    out[ok] <- y[ok] - mean(y[ok])
  } else {
    out[ok] <- correct_signal_matrix(matrix(y[ok], ncol = 1), progress[ok], method,
                                     n_progress_bins, spline_df)$residual[, 1]
  }
  out
}

#' Initiation-position and fork-asymmetry heterogeneity per IZ
#'
#' Per IZ, remove the replication-progress trend from the relative initiation
#' position and the fork asymmetry of its tracks ([correct_signal_matrix()];
#' plain centring with fewer than 4 tracks) and summarise means, medians and
#' the spread (SD, MAD, IQR, variance) of the residuals.
#'
#' @param track_intervals `prepare_iz_data()$track_intervals`.
#' @param method,n_progress_bins,spline_df Correction settings.
#' @param BPPARAM A `BiocParallel::BiocParallelParam` used over IZs.
#' @return A list with `summary` (one row per IZ), `track_intervals` (with
#'   `position_residual` and `fork_asymmetry_residual`) and
#'   `correction_method`.
#' @export
summarize_track_features <- function(track_intervals, method = c("spline", "bin"),
                                     n_progress_bins = 10L, spline_df = 4L,
                                     BPPARAM = BiocParallel::SerialParam()) {
  method <- match.arg(method)
  idx <- split(seq_len(nrow(track_intervals)), track_intervals$IZ_id)
  one <- function(ii) {
    d <- track_intervals[ii, , drop = FALSE]
    rp <- .correct_feature(d$relative_init_position, d$replication_progress,
                           method, n_progress_bins, spline_df)
    ra <- .correct_feature(d$fork_asymmetry, d$replication_progress,
                           method, n_progress_bins, spline_df)
    sm <- data.frame(
      IZ_id = d$IZ_id[1], mean_init_position = mean(d$relative_init_position, na.rm = TRUE),
      median_init_position = stats::median(d$relative_init_position, na.rm = TRUE),
      position_SD = stats::sd(rp, na.rm = TRUE), position_MAD = stats::mad(rp, na.rm = TRUE),
      position_IQR = stats::IQR(rp, na.rm = TRUE), position_variance = stats::var(rp, na.rm = TRUE),
      position_heterogeneity = stats::sd(rp, na.rm = TRUE),
      mean_left_track_length = mean(d$left_track_length, na.rm = TRUE),
      mean_right_track_length = mean(d$right_track_length, na.rm = TRUE),
      mean_total_track_length = mean(d$total_track_length, na.rm = TRUE),
      mean_fork_asymmetry = mean(d$fork_asymmetry, na.rm = TRUE),
      median_fork_asymmetry = stats::median(d$fork_asymmetry, na.rm = TRUE),
      fork_SD = stats::sd(ra, na.rm = TRUE), fork_MAD = stats::mad(ra, na.rm = TRUE),
      fork_heterogeneity = stats::sd(ra, na.rm = TRUE),
      n_active_cells = length(unique(d$cell_id))
    )
    list(index = ii, position_residual = rp, fork_residual = ra, summary = sm)
  }
  ans <- BiocParallel::bplapply(idx, one, BPPARAM = BPPARAM)
  track_intervals$position_residual <- NA_real_
  track_intervals$fork_asymmetry_residual <- NA_real_
  for (z in ans) {
    track_intervals$position_residual[z$index] <- z$position_residual
    track_intervals$fork_asymmetry_residual[z$index] <- z$fork_residual
  }
  list(summary = .bind_rows_df(lapply(ans, `[[`, "summary")),
       track_intervals = track_intervals, correction_method = method)
}

#' Track-shape heterogeneity per IZ
#'
#' For IZs with at least `min_active_cells` tracks, the binary coverage
#' profile of each track around the IZ centre is compared with its
#' progress-expected profile ([correct_signal_matrix()]): mean correlation
#' (`mean_track_similarity`; `1 - RMSE` where undefined), residual variance
#' and the residual error after a rank-`pca_rank` SVD.
#'
#' @param track_intervals,iz_metadata From [prepare_iz_data()].
#' @param method,n_progress_bins,spline_df Correction settings.
#' @param window_bp,bin_size Profile window half-width and step (bp).
#' @param min_active_cells Minimum number of tracks.
#' @param pca_rank,compute_pca SVD settings.
#' @param BPPARAM A `BiocParallel::BiocParallelParam` used over IZs.
#' @return A `data.frame`, one row per IZ.
#' @export
summarize_track_shape <- function(track_intervals, iz_metadata,
                                  method = c("spline", "bin"),
                                  n_progress_bins = 10L, spline_df = 4L,
                                  window_bp = 500000L, bin_size = 20000L,
                                  min_active_cells = 5L, pca_rank = 2L,
                                  compute_pca = TRUE,
                                  BPPARAM = BiocParallel::SerialParam()) {
  method <- match.arg(method)
  rel_bp <- seq(-window_bp, window_bp, by = bin_size)
  idx <- split(seq_len(nrow(track_intervals)), track_intervals$IZ_id)
  centers <- stats::setNames(iz_metadata$center, iz_metadata$IZ_id)
  one <- function(ii) {
    d <- track_intervals[ii, , drop = FALSE]
    id <- d$IZ_id[1]
    if (nrow(d) < min_active_cells) {
      return(data.frame(IZ_id = id, mean_track_similarity = NA_real_,
                        track_variance = NA_real_, track_pca_error = NA_real_,
                        track_shape_heterogeneity = NA_real_, shape_flag = "too_few_active"))
    }
    pos <- as.numeric(centers[as.character(id)]) + rel_bp
    observed <- (outer(d$start, pos, `<=`) & outer(d$end, pos, `>=`)) * 1
    corrected <- correct_signal_matrix(observed, d$replication_progress, method,
                                       n_progress_bins, spline_df)
    expected <- corrected$fitted
    a <- observed - rowMeans(observed)
    b <- expected - rowMeans(expected)
    den <- sqrt(rowSums(a * a) * rowSums(b * b))
    similarity <- rowSums(a * b) / den
    rmse <- sqrt(rowMeans(corrected$residual^2, na.rm = TRUE))
    similarity[!is.finite(similarity)] <- 1 - pmin(1, rmse[!is.finite(similarity)])
    track_var <- mean(corrected$residual^2, na.rm = TRUE)
    pca_error <- NA_real_
    if (isTRUE(compute_pca)) {
      r <- scale(corrected$residual, center = TRUE, scale = FALSE)
      ss <- sum(r^2, na.rm = TRUE)
      if (is.finite(ss) && ss == 0) pca_error <- 0
      else if (is.finite(ss)) {
        sv <- try(svd(r, nu = 0L, nv = 0L)$d, silent = TRUE)
        if (!inherits(sv, "try-error")) {
          k <- min(as.integer(pca_rank), length(sv))
          pca_error <- max(0, (ss - sum(utils::head(sv^2, k))) / length(r))
        }
      }
    }
    data.frame(
      IZ_id = id, mean_track_similarity = mean(similarity, na.rm = TRUE),
      track_variance = track_var, track_pca_error = pca_error,
      track_shape_heterogeneity = mean(1 - similarity, na.rm = TRUE), shape_flag = "ok"
    )
  }
  ans <- BiocParallel::bplapply(idx, one, BPPARAM = BPPARAM)
  .bind_rows_df(ans)
}

#' Join IZ metadata and heterogeneity summaries
#'
#' Left-join the summaries to `iz_metadata` by `IZ_id` (clashing column names
#' get `.x` / `.y` suffixes, as `dplyr::left_join()`) and order by `IZ_id`.
#'
#' @param iz_metadata `prepare_iz_data()$iz_metadata`.
#' @param usage_summary,timing_summary,feature_summary,shape_summary
#'   Per-IZ tables with an `IZ_id` column.
#' @return A `data.frame`.
#' @export
assemble_iz_summary <- function(iz_metadata, usage_summary, timing_summary,
                                feature_summary, shape_summary) {
  out <- as.data.frame(iz_metadata)
  for (y in list(usage_summary, timing_summary, feature_summary, shape_summary)) {
    out <- .left_join_by_id(out, as.data.frame(y), "IZ_id")
  }
  out <- out[order(out$IZ_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Left join on a unique key column (dplyr naming of clashing columns)
#' @keywords internal
#' @noRd
.left_join_by_id <- function(x, y, by) {
  if (anyDuplicated(y[[by]])) stop("Duplicated '", by, "' values in a summary table.")
  j <- match(x[[by]], y[[by]])
  add <- setdiff(names(y), by)
  clash <- intersect(add, names(x))
  if (length(clash)) {
    names(x)[match(clash, names(x))] <- paste0(clash, ".x")
  }
  for (nm in add) {
    new_nm <- if (nm %in% clash) paste0(nm, ".y") else nm
    x[[new_nm]] <- y[[nm]][j]
  }
  x
}
