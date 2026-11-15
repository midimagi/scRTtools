# Cell QC: MAPD and G1-cell identification.
# Migrated from rcode/processing/calculate_mapd.R and find_g1_cells.R
# (notebooks 1.2 create_Seurat_object_*).

#' Calculate per-cell coverage and MAPD from a tile count matrix
#'
#' Compute per-cell coverage, MAPD (median absolute pairwise difference) and
#' scaled MAPD from a tile-by-cell count matrix.
#'
#' Steps:
#' \enumerate{
#'   \item A per-tile GC correction value is computed from the pseudo-bulk
#'     (row sums): tiles are grouped by `round(gc_frequency, 1)`, and the value
#'     is `median(all reads/mappability) / median(group reads/mappability)`.
#'     Only mappable tiles with finite GC, positive mappability and finite
#'     correction are used further.
#'   \item For each cell, corrected reads are `reads / mappability * gc_corr`;
#'     `coverage = 1e6 * sum(corrected) / effective genome size` (sum of the
#'     widths of the used tiles).
#'   \item Corrected reads are summed in groups of `aggr_bins` consecutive used
#'     tiles per chromosome; MAPD is the median absolute deviation (around the
#'     median, unscaled) of the differences between adjacent groups on the same
#'     chromosome, and `scaled_mapd = round(mapd / sqrt(coverage), 3)`.
#' }
#'
#' @param tile_matrix A tile-by-cell count matrix (base or \pkg{Matrix}
#'   sparse matrix) with row names equal to `names(genome_granges)` and column
#'   names equal to `rownames(meta_data)`, in the same order.
#' @param genome_granges A named `GRanges` with one range per row of
#'   `tile_matrix` and metadata columns `gc_frequency`, `mappability` and
#'   `is_mappable` (see [compute_bin_mappability()]).
#' @param meta_data A `data.frame` of cell metadata, row names = cell IDs.
#' @param aggr_bins Positive integer; number of consecutive tiles aggregated
#'   before computing MAPD. We use `50`.
#' @param verbose Logical; print progress messages. Default `FALSE`.
#' @param BPPARAM A `BiocParallel::BiocParallelParam`; cells are processed in
#'   chunks in parallel. Results do not depend on the number of workers.
#'
#' @return `meta_data` with added or updated columns `coverage`, `mapd`,
#'   `scaled_mapd` (numeric) and `mapd_error` (character, `NA` when the cell
#'   succeeded).
#'
#' @examples
#' set.seed(1)
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(seq(1, by = 1000, length.out = 200),
#'                                                       width = 1000))
#' names(gr) <- paste0("bin", seq_along(gr))
#' gr$gc_frequency <- runif(200, 0.35, 0.5)
#' gr$mappability <- 1
#' gr$is_mappable <- TRUE
#' m <- matrix(rpois(200 * 5, 5), 200, 5, dimnames = list(names(gr), paste0("c", 1:5)))
#' m <- Matrix::Matrix(m, sparse = TRUE)
#' md <- data.frame(row.names = colnames(m), cell_id = colnames(m))
#' calculate_mapd(m, gr, md, aggr_bins = 10, verbose = FALSE)
#'
#' @export
calculate_mapd <- function(tile_matrix,
                           genome_granges,
                           meta_data,
                           aggr_bins = 50,
                           verbose = FALSE,
                           BPPARAM = BiocParallel::SerialParam()) {
  # -----------------------------
  # 1. Parameter verification
  # -----------------------------
  if (!is.matrix(tile_matrix) && !inherits(tile_matrix, "Matrix")) {
    stop("'tile_matrix' must be a base matrix or Matrix sparse matrix.")
  }

  .validate_bins(genome_granges, "genome_granges")

  if (nrow(tile_matrix) != length(genome_granges)) {
    stop("Number of rows in 'tile_matrix' must equal length of 'genome_granges'.")
  }

  if (!is.data.frame(meta_data)) {
    stop("'meta_data' must be a data.frame.")
  }

  cell_ids <- colnames(tile_matrix)
  if (is.null(cell_ids) || length(cell_ids) == 0L) {
    stop("'tile_matrix' must contain column names, i.e. cell IDs.")
  }

  if (is.null(rownames(meta_data))) {
    stop("'meta_data' must have row names matching cell IDs.")
  }

  required_cols <- c("gc_frequency", "mappability", "is_mappable")
  range_mcols <- names(S4Vectors::mcols(genome_granges))
  missing_cols <- setdiff(required_cols, range_mcols)
  if (length(missing_cols) > 0L) {
    stop(
      "Missing required columns in mcols(genome_granges): ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # Check tile_matrix row names and genome_granges names are identical
  tile_bin_ids <- rownames(tile_matrix)
  grange_bin_ids <- names(genome_granges)

  if (is.null(tile_bin_ids) || length(tile_bin_ids) == 0L) {
    stop("'tile_matrix' must contain row names")
  }

  if (is.null(grange_bin_ids) || length(grange_bin_ids) == 0L) {
    stop("'genome_granges' must contain range names, ie names(genome_granges)")
  }

  if (length(tile_bin_ids) != length(grange_bin_ids)) {
    stop("Length mismatch between tile_matrix row names and names(genome_granges).")
  }

  if (!identical(as.character(tile_bin_ids), grange_bin_ids)) {
    stop(
      "'tile_matrix' row names and names(genome_granges) are not identical ",
      "or not in the same order."
    )
  }

  # Check tile_matrix column names and meta_data row names are identical.
  # (The original read tile_matrix@Dimnames[[2]], which fails for base
  # matrices; colnames() returns the same value for Matrix objects.)
  tile_cell_ids <- colnames(tile_matrix)
  meta_cell_ids <- rownames(meta_data)

  if (length(tile_cell_ids) != length(meta_cell_ids)) {
    stop("Length mismatch between tile_matrix column names and rownames(meta_data).")
  }

  if (!identical(as.character(tile_cell_ids), as.character(meta_cell_ids))) {
    stop(
      "'tile_matrix' column names and rownames(meta_data) are not identical ",
      "or not in the same order."
    )
  }

  if (!is.numeric(aggr_bins) || length(aggr_bins) != 1L ||
      is.na(aggr_bins) || aggr_bins < 1) {
    stop("'aggr_bins' must be a single positive integer.")
  }
  aggr_bins <- as.integer(aggr_bins)

  if (!inherits(BPPARAM, "BiocParallelParam")) {
    stop("'BPPARAM' must be a BiocParallelParam object.")
  }

  # -----------------------------
  # 2. Compute GC correction values
  # -----------------------------
  .vmsg(verbose, "\n1. Computing GC correction values...")

  gc_frequency <- as.numeric(S4Vectors::mcols(genome_granges)$gc_frequency)
  mappability <- as.numeric(S4Vectors::mcols(genome_granges)$mappability)
  is_mappable <- as.logical(S4Vectors::mcols(genome_granges)$is_mappable)

  bin_sum <- as.numeric(Matrix::rowSums(tile_matrix))

  gc_corr <- rep(NA_real_, length(genome_granges))

  gc_idx <- which(
    is_mappable &
      is.finite(gc_frequency) &
      is.finite(mappability) &
      (mappability > 0) &
      is.finite(bin_sum)
  )

  if (length(gc_idx) > 0L) {
    gc_round <- round(gc_frequency[gc_idx], 1)
    reads_mappability <- bin_sum[gc_idx] / mappability[gc_idx]
    global_median <- stats::median(reads_mappability, na.rm = TRUE)

    if (is.finite(global_median)) {
      gc_groups <- split(seq_along(gc_idx), gc_round)
      for (idx_local in gc_groups) {
        mgc <- stats::median(reads_mappability[idx_local], na.rm = TRUE)
        corr_val <- if (is.na(mgc) || mgc == 0) {
          NA_real_
        } else {
          global_median / mgc
        }
        gc_corr[gc_idx[idx_local]] <- corr_val
      }
    }
  }

  # -----------------------------
  # 3. Prepare static layout for MAPD
  # -----------------------------
  .vmsg(verbose, "\n2. Preparing data for MAPD calculation...")

  valid_bins <- which(
    is_mappable &
      is.finite(mappability) &
      (mappability > 0) &
      is.finite(gc_corr)
  )

  if (length(valid_bins) == 0L) {
    warning("No valid bins available for MAPD after filtering.", call. = FALSE)

    out_cols <- c("coverage", "mapd", "scaled_mapd", "mapd_error")
    for (nm in out_cols) {
      if (!nm %in% colnames(meta_data)) {
        meta_data[[nm]] <- if (nm == "mapd_error") NA_character_ else NA_real_
      }
    }

    # NOTE (kept from original): indexing an unnamed column by row names
    # extends the vector instead of filling it, so this line errors.
    # Reported as a suspected bug.
    meta_data$mapd_error[rownames(meta_data)] <- "No valid bins available for MAPD."
    return(meta_data)
  }

  seqlevels_all <- GenomeInfoDb::seqlevels(genome_granges)
  bin_chr <- as.character(GenomicRanges::seqnames(genome_granges)[valid_bins])
  bin_start <- GenomicRanges::start(genome_granges)[valid_bins]
  bin_end <- GenomicRanges::end(genome_granges)[valid_bins]
  bin_width <- IRanges::width(genome_granges)[valid_bins]

  chr_factor <- factor(bin_chr, levels = seqlevels_all)
  ord <- order(chr_factor, bin_start, bin_end, na.last = NA)

  valid_bins <- valid_bins[ord]
  bin_chr <- bin_chr[ord]
  bin_width <- bin_width[ord]
  map_valid <- mappability[valid_bins]
  gc_corr_valid <- gc_corr[valid_bins]

  effective_genome_size <- sum(bin_width, na.rm = TRUE)
  if (!is.finite(effective_genome_size) || effective_genome_size <= 0) {
    stop("Computed effective genome size is not valid.")
  }

  chr_rle <- rle(bin_chr)
  n_aggr <- sum((chr_rle$lengths + aggr_bins - 1L) %/% aggr_bins)

  agg_group <- integer(length(valid_bins))
  agg_chr <- character(n_aggr)

  pos_bin <- 1L
  pos_group <- 1L
  group_offset <- 0L

  for (k in seq_along(chr_rle$lengths)) {
    n_chr_bins <- chr_rle$lengths[k]
    chr_name <- chr_rle$values[k]
    n_chr_groups <- (n_chr_bins + aggr_bins - 1L) %/% aggr_bins

    idx_bin <- pos_bin:(pos_bin + n_chr_bins - 1L)
    local_group <- (seq_len(n_chr_bins) - 1L) %/% aggr_bins + 1L
    agg_group[idx_bin] <- group_offset + local_group

    idx_group <- pos_group:(pos_group + n_chr_groups - 1L)
    agg_chr[idx_group] <- chr_name

    pos_bin <- pos_bin + n_chr_bins
    pos_group <- pos_group + n_chr_groups
    group_offset <- group_offset + n_chr_groups
  }

  same_chr_adj <- if (length(agg_chr) >= 2L) {
    agg_chr[-1L] == agg_chr[-length(agg_chr)]
  } else {
    logical(0)
  }

  # -----------------------------
  # 4. Chunked parallelism with BiocParallel
  # -----------------------------
  .vmsg(verbose, "\n3. Calculating MAPD with BiocParallel...")

  n_cells <- length(cell_ids)
  n_workers <- BiocParallel::bpworkers(BPPARAM)

  .vmsg(verbose, "   Workers: ", n_workers)
  .vmsg(verbose, "   Cells: ", n_cells)

  chunk_size <- max(1L, ceiling(n_cells / max(1L, n_workers * 2L)))
  chunk_id <- ceiling(seq_len(n_cells) / chunk_size)
  cell_chunks <- split(seq_len(n_cells), chunk_id)

  .vmsg(verbose, "   Chunk size: ", chunk_size)
  .vmsg(verbose, "   Number of chunks: ", length(cell_chunks))

  .mapd_worker <- function(idx_chunk) {
    chunk_cell_ids <- cell_ids[idx_chunk]
    chunk_counts <- tile_matrix[valid_bins, idx_chunk, drop = FALSE]

    n_chunk <- ncol(chunk_counts)
    res <- data.frame(
      cell_id = chunk_cell_ids,
      coverage = rep(NA_real_, n_chunk),
      mapd = rep(NA_real_, n_chunk),
      scaled_mapd = rep(NA_real_, n_chunk),
      mapd_error = rep(NA_character_, n_chunk),
      stringsAsFactors = FALSE
    )

    for (j in seq_len(n_chunk)) {
      tryCatch({
        reads <- as.numeric(chunk_counts[, j, drop = TRUE])

        reads_mpgc <- (reads / map_valid) * gc_corr_valid
        if (!any(is.finite(reads_mpgc))) {
          stop("No finite corrected read values for this cell.")
        }

        coverage <- 1e6 * sum(reads_mpgc, na.rm = TRUE) / effective_genome_size

        aggr_sum <- rowsum(
          x = reads_mpgc,
          group = agg_group,
          reorder = FALSE
        )
        aggr_sum <- as.numeric(aggr_sum[, 1])

        if (length(aggr_sum) < 2L || !any(same_chr_adj)) {
          mapd <- NA_real_
          scaled_mapd <- NA_real_
        } else {
          diff_reads <- aggr_sum[-1L] - aggr_sum[-length(aggr_sum)]
          diff_reads <- diff_reads[same_chr_adj]
          diff_reads <- diff_reads[is.finite(diff_reads)]

          if (length(diff_reads) == 0L) {
            mapd <- NA_real_
            scaled_mapd <- NA_real_
          } else {
            med_diff <- stats::median(diff_reads, na.rm = TRUE)
            mapd <- stats::median(abs(diff_reads - med_diff), na.rm = TRUE)

            scaled_mapd <- if (is.na(mapd) || is.na(coverage) || coverage <= 0) {
              NA_real_
            } else {
              round(mapd / sqrt(coverage), 3)
            }
          }
        }

        res$coverage[j] <- coverage
        res$mapd[j] <- mapd
        res$scaled_mapd[j] <- scaled_mapd
      }, error = function(e) {
        # NOTE (kept from original): this assigns to a local copy of 'res'
        # inside the handler, so the error message is never recorded and
        # the failure warning below never fires. Reported as a suspected bug.
        res$mapd_error[j] <- conditionMessage(e)
      })
    }

    res
  }

  mapd_list <- BiocParallel::bplapply(
    X = cell_chunks,
    FUN = .mapd_worker,
    BPPARAM = BPPARAM
  )

  mapd_results <- do.call(rbind, mapd_list)
  rownames(mapd_results) <- mapd_results$cell_id

  # -----------------------------
  # 5. Write results back to metadata
  # -----------------------------
  .vmsg(verbose, "\n4. Adding MAPD results to metadata...")

  common_cells <- intersect(rownames(meta_data), rownames(mapd_results))

  if (!"coverage" %in% colnames(meta_data)) meta_data$coverage <- NA_real_
  if (!"mapd" %in% colnames(meta_data)) meta_data$mapd <- NA_real_
  if (!"scaled_mapd" %in% colnames(meta_data)) meta_data$scaled_mapd <- NA_real_
  if (!"mapd_error" %in% colnames(meta_data)) meta_data$mapd_error <- NA_character_

  meta_data[common_cells, "coverage"] <- mapd_results[common_cells, "coverage"]
  meta_data[common_cells, "mapd"] <- mapd_results[common_cells, "mapd"]
  meta_data[common_cells, "scaled_mapd"] <- mapd_results[common_cells, "scaled_mapd"]
  meta_data[common_cells, "mapd_error"] <- mapd_results[common_cells, "mapd_error"]

  # -----------------------------
  # 6. Summary
  # -----------------------------
  n_failed <- sum(!is.na(mapd_results$mapd_error) & nzchar(mapd_results$mapd_error))

  if (n_failed > 0L) {
    warning(
      n_failed, " cell(s) failed during MAPD calculation. ",
      "See 'meta_data$mapd_error' for details.",
      call. = FALSE
    )
  }

  meta_data
}

#' Identify G1-phase cells from coverage and scaled MAPD
#'
#' Classify cells as G1 using the relationship between coverage and scaled
#' MAPD. A two-component linear model is fitted with a hard-assignment EM
#' (internal `.find_g1_coefficients_em()`); the component with the lower predicted
#' MAPD at the maximum coverage is taken as G1, and cells whose residual from
#' that line is `< 0.05` are classified as G1. With fewer than 4 cells a
#' single linear regression is used instead (with a warning).
#'
#' @param meta_data A `data.frame` with columns `cell_id`, `coverage`,
#'   `scaled_mapd` and `sample` (e.g. output of [calculate_mapd()]); row names
#'   must equal `cell_id`.
#' @param by_sample Logical; fit the model separately for each `sample`.
#'   We use `FALSE`.
#' @param seed Optional integer. The EM breaks likelihood ties and re-seeds
#'   degenerate assignments at random; if `seed` is given, this runs inside
#'   `withr::with_seed(seed, ...)`. With `NULL` (default) the current RNG state
#'   is used.
#' @param verbose Logical; print EM iteration messages. Default `FALSE`.
#'
#' @return `meta_data` **restricted to rows with non-`NA` `coverage` and
#'   `scaled_mapd`**, with an added/updated logical column `is_G1`.
#'
#' @examples
#' set.seed(1)
#' cov <- runif(200, 0.5, 3)
#' g1 <- rep(c(TRUE, FALSE), c(120, 80))
#' smapd <- ifelse(g1, 0.2 - 0.02 * cov, 0.5 - 0.02 * cov) + rnorm(200, 0, 0.01)
#' md <- data.frame(cell_id = paste0("c", 1:200), coverage = cov,
#'                  scaled_mapd = smapd, sample = "s1")
#' rownames(md) <- md$cell_id
#' res <- find_g1_cells(md, seed = 1, verbose = FALSE)
#' table(res$is_G1, g1)
#'
#' @export
find_g1_cells <- function(meta_data,
                          by_sample = FALSE,
                          seed = NULL,
                          verbose = FALSE) {
  if (is.null(seed)) {
    .find_g1_cells_impl(meta_data, by_sample, verbose)
  } else {
    withr::with_seed(seed, .find_g1_cells_impl(meta_data, by_sample, verbose))
  }
}

.find_g1_cells_impl <- function(meta_data, by_sample, verbose) {
  if (!is.data.frame(meta_data)) {
    stop("'meta_data' must be a data.frame.")
  }
  required_cols <- c("cell_id", "coverage", "scaled_mapd", "sample")
  if (!all(required_cols %in% names(meta_data))) {
    missing_cols <- setdiff(required_cols, names(meta_data))
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  # Drop cells with NA coverage / scaled MAPD
  valid <- !is.na(meta_data$coverage) & !is.na(meta_data$scaled_mapd)
  meta_data <- meta_data[valid, ]

  if (by_sample) {
    result_list <- list()
    samples <- unique(meta_data$sample)

    for (samp in samples) {
      sample_df <- meta_data[meta_data$sample == samp, ]

      if (nrow(sample_df) >= 4) {
        g1_coefficients <- .find_g1_coefficients_em(
          coverage = sample_df$coverage,
          scaled_mapd = sample_df$scaled_mapd,
          verbose = verbose
        )
      } else {
        warning(paste("Sample", samp, "has too few cells (<4) for EM. Using simple linear regression."))
        fit <- stats::lm(scaled_mapd ~ coverage, data = sample_df)
        g1_coefficients <- stats::coef(fit)
      }

      predicted_mapd <- g1_coefficients[1] + g1_coefficients[2] * sample_df$coverage
      residual_mapd <- sample_df$scaled_mapd - predicted_mapd
      sample_df$is_G1 <- residual_mapd < 0.05

      result_list[[samp]] <- sample_df
    }

    result_df <- do.call(rbind, result_list)
    rownames(result_df) <- NULL
  } else {
    if (nrow(meta_data) >= 4) {
      g1_coefficients <- .find_g1_coefficients_em(
        coverage = meta_data$coverage,
        scaled_mapd = meta_data$scaled_mapd,
        verbose = verbose
      )
    } else {
      warning("Too few cells (<4) for EM. Using simple linear regression.")
      fit <- stats::lm(scaled_mapd ~ coverage, data = meta_data)
      g1_coefficients <- stats::coef(fit)
    }

    predicted_mapd <- g1_coefficients[1] + g1_coefficients[2] * meta_data$coverage
    residual_mapd <- meta_data$scaled_mapd - predicted_mapd
    meta_data$is_G1 <- residual_mapd < 0.05

    result_df <- meta_data
  }

  # Write is_G1 back by cell_id (requires rownames(meta_data) == cell_id)
  result_df <- result_df[, c("cell_id", "is_G1")]
  rownames(result_df) <- result_df$cell_id
  common_cells <- intersect(rownames(meta_data), rownames(result_df))
  meta_data[common_cells, "is_G1"] <- result_df[common_cells, "is_G1"]

  meta_data
}

#' Estimate G1 linear-model coefficients with a two-component EM
#'
#' Fit two lines of scaled MAPD versus coverage with a hard-assignment EM:
#' cells start in group 2 if their scaled MAPD is at or above the 10th
#' percentile (group 1 otherwise); each iteration fits a line per group,
#' computes normal likelihoods of all residuals (SD from the group's own
#' residuals), and reassigns each cell to the most likely group (ties broken
#' at random). If a group empties, assignments are re-drawn at random. The
#' loop stops when assignments no longer change or after 100 iterations. The
#' group with the lower predicted MAPD at the maximum coverage is G1.
#'
#' @param coverage,scaled_mapd Numeric vectors of equal length.
#' @param verbose Logical; print the iteration counter. Default `FALSE`.
#'
#' @return Named numeric vector `c(intercept, slope)` of the G1 line (unnamed
#'   `lm` coefficients `(Intercept)`, `coverage` for the < 4 cells fallback).
#'
#' @keywords internal
#' @noRd
.find_g1_coefficients_em <- function(coverage, scaled_mapd, verbose = FALSE) {
  num_cells <- length(coverage)
  if (length(scaled_mapd) != num_cells) {
    stop("coverage and scaled_mapd must have the same length")
  }

  valid <- !is.na(coverage) & !is.na(scaled_mapd)
  coverage <- coverage[valid]
  scaled_mapd <- scaled_mapd[valid]
  num_cells <- length(coverage)

  if (num_cells < 4) {
    warning("Too few cells for two-component model fitting. Using simple linear regression.")
    fit <- stats::lm(scaled_mapd ~ coverage)
    return(stats::coef(fit))
  }

  # 1. Initial grouping by the 10th percentile of scaled MAPD
  group <- rep(1, num_cells)
  threshold <- stats::quantile(scaled_mapd, probs = 0.10, na.rm = TRUE)
  is_high_mapd <- scaled_mapd >= threshold
  group[is_high_mapd] <- 2

  previous_assignments <- rep(1, num_cells)

  coefficients <- matrix(0, nrow = 2, ncol = 2)
  residuals <- matrix(0, nrow = 2, ncol = num_cells)
  likelihoods <- matrix(0, nrow = 2, ncol = num_cells)

  # 2. EM iterations
  iteration <- 0
  max_iterations <- 100

  while (any(previous_assignments != group) && iteration < max_iterations) {
    iteration <- iteration + 1

    .vmsg(verbose, "current iteration: ", iteration)

    for (group_i in 1:2) {
      group_idx <- which(group == group_i)

      if (length(group_idx) >= 2) {
        fit <- stats::lm(scaled_mapd[group_idx] ~ coverage[group_idx])
        coefs <- stats::coef(fit)
        coefficients[group_i, ] <- coefs

        predicted <- coefficients[group_i, 1] + coefficients[group_i, 2] * coverage
        residuals[group_i, ] <- scaled_mapd - predicted

        resid_sd <- stats::sd(residuals[group_i, group_idx], na.rm = TRUE)
        if (is.na(resid_sd) || resid_sd == 0) {
          resid_sd <- 1e-6
        }
        likelihoods[group_i, ] <- stats::dnorm(residuals[group_i, ], mean = 0, sd = resid_sd)
      } else {
        coefficients[group_i, ] <- c(0, 0)
        residuals[group_i, ] <- 0
        likelihoods[group_i, ] <- 1e-6
      }
    }

    previous_assignments <- group

    for (barcode in 1:num_cells) {
      max_likelihood <- max(likelihoods[, barcode], na.rm = TRUE)
      if (!is.na(max_likelihood)) {
        max_indices <- which(likelihoods[, barcode] == max_likelihood)
        group[barcode] <- if (length(max_indices) > 1) {
          sample(max_indices, 1)
        } else {
          max_indices
        }
      }
    }

    if (length(unique(group)) < 2) {
      group <- sample(1:2, num_cells, replace = TRUE, prob = c(0.5, 0.5))
    }
  }

  # 3. Prediction of each line at the maximum coverage
  max_coverage <- max(coverage, na.rm = TRUE)
  predictions <- numeric(2)

  for (c in 1:2) {
    predictions[c] <- coefficients[c, 1] + coefficients[c, 2] * max_coverage
  }

  # 4. The lower line is G1
  g1_group <- which.min(predictions)
  g1_coefficients <- coefficients[g1_group, ]

  if (any(is.na(g1_coefficients))) {
    warning("G1 population detection failed. Using linear regression on all cells.")
    fit <- stats::lm(scaled_mapd ~ coverage)
    g1_coefficients <- stats::coef(fit)
  }

  names(g1_coefficients) <- c("intercept", "slope")

  g1_coefficients
}
