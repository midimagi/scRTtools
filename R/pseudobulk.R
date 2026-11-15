# Pseudobulk replication-timing signal.
# calculate_pseudobulk_rt() is migrated from
# rcode/processing/calculate_pseudobulk_rt.R (notebooks 2.x / 4.x); the S/G1
# ratio functions at the end of the file from calculate_pseudobulk_rt_g1_ratio.R.

#' Calculate a pseudobulk replication-timing signal
#'
#' Aggregate QC-passed S-phase cells (`pass_qc & !is_G1`) and G1 cells
#' (`is_G1`), correct the S-phase counts for GC content and mappability using
#' the G1 pseudobulk, normalise by the genome-wide mean, z-score, set bins
#' with `|z| > extreme_cutoff` to `NA`, smooth with `stats::smooth.spline()`
#' within mappable contigs, then z-score and filter the smoothed signal again.
#'
#' GC/mappability correction (on mappable bins): with
#' `g1m = g1_counts / mappability`, `mGC` = median of `g1m` within bins of the
#' same `round(gc_frequency, 2)` and `m` = median of all `g1m`, the corrected
#' S signal is `(m / mGC) * s_counts / mappability` (`NA` where `mGC == 0`).
#' Contigs are mappable bins merged with `GenomicRanges::reduce(min.gapwidth =
#' min.gapwidth)`; contigs with fewer than `min_bins_per_contig` bins or fewer
#' than `min_non_na_points` distinct positions are not smoothed.
#'
#' @param genobins_matrix `dgCMatrix` of raw counts, genome bins x cells, with
#'   row names equal to `names(genome_granges)` and column names equal to
#'   `rownames(meta_data)`.
#' @param genome_granges Named `GRanges` with metadata columns `is_mappable`,
#'   `gc_frequency` and `mappability` (e.g. from [compute_bin_mappability()]).
#' @param meta_data A `data.frame` with logical columns `is_G1` and `pass_qc`.
#' @param extreme_cutoff Absolute z-score above which bins are set to `NA`.
#' @param min.gapwidth Passed to `GenomicRanges::reduce()` to build contigs.
#' @param min_bins_per_contig Minimum number of bins of a smoothed contig.
#' @param min_non_na_points Minimum number of non-`NA` points for the spline.
#' @param smooth_spar `NULL` (automatic) or `spar` of
#'   `stats::smooth.spline()`.
#'
#' @return A sorted `GRanges` of the bins with a non-`NA` smoothed signal,
#'   with the single metadata column `score`.
#'
#' @export
calculate_pseudobulk_rt <- function(
  genobins_matrix, 
  genome_granges, 
  meta_data,
  extreme_cutoff = 5.0,
  min.gapwidth = 40001L,
  min_bins_per_contig = 5,
  min_non_na_points = 2,
  smooth_spar = NULL
) {
    
  stopifnot(
    inherits(genobins_matrix, "dgCMatrix"),
    inherits(genome_granges, "GRanges"),
    all(rownames(genobins_matrix) == names(genome_granges)),
    all(colnames(genobins_matrix) == rownames(meta_data))
  )

  if (!is.null(smooth_spar)) {
    if (
      !is.numeric(smooth_spar) ||
        length(smooth_spar) != 1 ||
        !is.finite(smooth_spar)
    ) {
      stop("smooth_spar must be NULL or a finite numeric scalar.")
    }
  }

  req_meta <- c("is_G1", "pass_qc")
  req_rng <- c("is_mappable", "gc_frequency", "mappability")

  if (length(m <- setdiff(req_meta, colnames(meta_data)))) {
    stop("meta_data is missing columns: ", paste(m, collapse = ", "))
  }

  if (length(m <- setdiff(req_rng, colnames(S4Vectors::mcols(genome_granges))))) {
    stop("mcols(genome_granges) is missing columns: ", paste(m, collapse = ", "))
  }

  zscore <- function(x, label = "signal") {
    ok <- !is.na(x)

    if (sum(ok) < 2) {
      warning(label, ": fewer than two non-NA values; returning all NA.")
      return(rep(NA_real_, length(x)))
    }

    mu <- mean(x[ok])
    s <- stats::sd(x[ok])

    if (!is.finite(s) || s == 0) {
      warning(label, ": standard deviation is zero or non-finite; returning all NA.")
      return(rep(NA_real_, length(x)))
    }

    out <- rep(NA_real_, length(x))
    out[ok] <- (x[ok] - mu) / s
    out
  }

  remove_extreme <- function(x) {
    i <- !is.na(x) & abs(x) > extreme_cutoff
    x[i] <- NA_real_
    x
  }

  # ------------------------------------------------------------
  # Build mappable contigs from fixed 20 kb mappable bins
  # ------------------------------------------------------------

  is_map <- as.logical(S4Vectors::mcols(genome_granges)$is_mappable)
  is_map[is.na(is_map)] <- FALSE

  mappable_gr <- genome_granges[is_map]
  mappable_gr <- GenomicRanges::reduce(mappable_gr, min.gapwidth = min.gapwidth)
  mappable_gr <- BiocGenerics::sort(mappable_gr)

  if (!length(mappable_gr)) {
    stop("No mappable regions found from genome_granges$is_mappable.")
  }

  mappable_gr$contig_id <- seq_along(mappable_gr)

  hit <- GenomicRanges::findOverlaps(
    genome_granges,
    mappable_gr,
    type = "within",
    ignore.strand = TRUE
  )

  contig_id <- rep(NA_integer_, length(genome_granges))
  contig_id[S4Vectors::queryHits(hit)] <- mappable_gr$contig_id[S4Vectors::subjectHits(hit)]

  midpoint <- floor((BiocGenerics::start(genome_granges) + BiocGenerics::end(genome_granges)) / 2)
  is_in_mappable_contig <- !is.na(contig_id)

  # ------------------------------------------------------------
  # Select S-phase cells
  # ------------------------------------------------------------

  s_cells <- which(meta_data$pass_qc & !meta_data$is_G1)
  if (!length(s_cells)) {
    stop("No S-phase cells found: meta_data$pass_qc == TRUE & meta_data$is_G1 == FALSE.")
  }
  s_counts <- as.numeric(Matrix::rowSums(genobins_matrix[, s_cells, drop = FALSE]))

  g1_cells <- which(meta_data$is_G1)
  g1_counts <- as.numeric(Matrix::rowSums(genobins_matrix[, g1_cells, drop = FALSE]))

  # GC / mappability correction from the G1 pseudobulk. Base-R version of the
  # original dplyr chain (left_join by bin_id, filter(is_mappable),
  # group_by(round(gc_frequency, 2)), median ratios); same arithmetic.
  gc_cov <- .pseudobulk_gc_correct(
    bin_id = rownames(genobins_matrix),
    g1_counts = g1_counts,
    s_counts = s_counts,
    genome_granges = genome_granges
  )

  s_counts <- rep(NA_real_, nrow(genobins_matrix))
  s_counts[gc_cov$bin_idx] <- gc_cov$norm_s

  s_counts[!is_in_mappable_contig] <- NA_real_

  # ------------------------------------------------------------
  # Normalize by genome-wide mean coverage
  # ------------------------------------------------------------

  mean_coverage <- mean(s_counts, na.rm = TRUE)

  if (!is.finite(mean_coverage) || mean_coverage == 0) {
    stop("Mean coverage is NA, non-finite, or zero. Cannot normalize.")
  }

  normalized_counts <- s_counts / mean_coverage

  # ------------------------------------------------------------
  # Z-score and remove extreme bins
  # ------------------------------------------------------------
  unsmoothed <- remove_extreme(zscore(normalized_counts, "normalized_counts"))

  # ------------------------------------------------------------
  # Smooth within chromosome-specific mappable contigs
  # ------------------------------------------------------------

  group_id <- paste0(
    as.character(GenomicRanges::seqnames(genome_granges)),
    "-",
    contig_id
  )

  group_id[is.na(contig_id)] <- NA_character_

  smoothed <- rep(NA_real_, length(unsmoothed))

  groups <- base::split(which(!is.na(group_id)), group_id[!is.na(group_id)])

  for (idx in groups) {
    if (length(idx) < min_bins_per_contig) next

    ok <- !is.na(midpoint[idx]) & !is.na(unsmoothed[idx])

    if (sum(ok) < min_non_na_points) next

    df <- data.frame(
      x = midpoint[idx][ok],
      y = unsmoothed[idx][ok]
    )

    df <- stats::aggregate(y ~ x, df, mean)
    df <- df[order(df$x), , drop = FALSE]

    if (nrow(df) < min_non_na_points) next

    fit <- tryCatch(
      {
        if (is.null(smooth_spar)) {
          stats::smooth.spline(df$x, df$y)
        } else {
          stats::smooth.spline(df$x, df$y, spar = smooth_spar)
        }
      },
      error = function(e) NULL
    )

    if (is.null(fit)) next

    pred_idx <- idx[!is.na(midpoint[idx])]
    smoothed[pred_idx] <- stats::predict(fit, x = midpoint[pred_idx])$y
  }

  # ------------------------------------------------------------
  # Re-standardize smoothed signal and remove extreme bins
  # ------------------------------------------------------------

  smoothed_z <- remove_extreme(zscore(smoothed, "smoothed"))

  # ------------------------------------------------------------
  # return GRanges
  # ------------------------------------------------------------

  bw_gr <- genome_granges
  S4Vectors::mcols(bw_gr) <- S4Vectors::DataFrame(score = smoothed_z)
  bw_gr <- bw_gr[!is.na(bw_gr$score)]
  bw_gr <- BiocGenerics::sort(bw_gr)

  return(bw_gr)
}

#' GC / mappability-corrected S counts on mappable bins
#' @keywords internal
#' @noRd
.pseudobulk_gc_correct <- function(bin_id, g1_counts, s_counts, genome_granges) {
  mc <- .mcols_df(S4Vectors::mcols(genome_granges))
  # left_join(df_counts, mcols, by = "bin_id") with unique bin ids
  j <- match(bin_id, names(genome_granges))
  keep <- which(mc$is_mappable[j] %in% TRUE)        # filter(is_mappable)

  gc_frequency <- round(mc$gc_frequency[j][keep], 2)
  mappability <- mc$mappability[j][keep]
  g1_mappability <- g1_counts[keep] / mappability
  # group_by(gc_frequency) keeps NA as its own group
  mGC <- stats::ave(g1_mappability, factor(gc_frequency, exclude = NULL),
                    FUN = stats::median)
  m <- stats::median(g1_mappability)
  gc_cov <- ifelse(mGC == 0, NA, m / mGC)

  data.frame(
    bin_id = bin_id[keep],
    bin_idx = keep,
    norm_s = gc_cov * s_counts[keep] / mappability,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# S/G1 ratio pseudobulk (per-cell CP10K). Migrated from
# rcode/processing/calculate_pseudobulk_rt_g1_ratio.R and the glue function
# export_pseudobulk_g1_ratio() of retreat_figure6_helpers.R.
# ---------------------------------------------------------------------------

#' Aggregate per-cell CP10K counts of S and G1 cells over several samples
#'
#' Reads one raw genome-bin count matrix per sample, normalises every selected
#' cell to counts per 10,000 (the denominator is the cell's total over all raw
#' bins, before any mappability filter) and sums the normalised counts per bin
#' over the QC-passed S cells (`is_G1` is `FALSE` and `pass_qc` is `TRUE`) and
#' over the G1 cells (`is_G1` is `TRUE`). The matrices are read one at a time.
#'
#' Cell names of a matrix that do not start with `"<sample>-"` get that
#' prefix before they are matched to `rownames(meta_data)`; every cell of
#' every matrix must be present in `meta_data`. Each sample must contribute at
#' least one S cell and one G1 cell.
#'
#' @param genobins_files Named character vector of RDS files, one sparse
#'   bins x cells count matrix per sample; the names are the sample ids.
#' @param meta_data A `data.frame` with cell ids as row names and logical
#'   columns `is_G1` and `pass_qc`.
#' @param bin_ids Character vector of bin ids. Each matrix must have exactly
#'   these row names (rows are reordered if needed).
#' @param cache_file `NULL` or the path of an RDS cache. The cache is reused
#'   only when the files (path, size, modification time), the selected cells
#'   and `bin_ids` are unchanged; otherwise it is overwritten.
#' @param verbose Print one message per sample. Default `FALSE`.
#'
#' @return A `data.frame` with columns `bin_id`, `s_counts` and `g1_counts`
#'   (summed CP10K, not raw reads), and attributes `normalization`
#'   (`"cp10k_per_cell_v1"`) and `selected_cells` (cells used per sample).
#'
#' @seealso [calculate_pseudobulk_rt_g1_ratio()],
#'   [export_pseudobulk_g1_ratio()]
#' @export
aggregate_pseudobulk_s_g1_counts <- function(
    genobins_files, meta_data, bin_ids, cache_file = NULL, verbose = FALSE
) {
    stopifnot(
        length(genobins_files) > 0L,
        !is.null(names(genobins_files)),
        !anyDuplicated(names(genobins_files)),
        all(file.exists(genobins_files)),
        !anyDuplicated(rownames(meta_data)),
        !anyDuplicated(bin_ids),
        all(c("is_G1", "pass_qc") %in% colnames(meta_data))
    )
    s_cells <- rownames(meta_data)[
        meta_data$is_G1 %in% FALSE & meta_data$pass_qc %in% TRUE
    ]
    g1_cells <- rownames(meta_data)[meta_data$is_G1 %in% TRUE]
    signature <- list(
        normalization = "cp10k_per_cell_v1",
        files = normalizePath(genobins_files),
        file_info = file.info(genobins_files)[, c("size", "mtime")],
        s_cells = s_cells, g1_cells = g1_cells, bin_ids = bin_ids
    )
    if (!is.null(cache_file) && file.exists(cache_file)) {
        cached <- readRDS(cache_file)
        if (identical(cached$signature, signature)) return(cached$counts)
    }

    counts <- data.frame(
        bin_id = bin_ids, s_counts = numeric(length(bin_ids)),
        g1_counts = numeric(length(bin_ids)), stringsAsFactors = FALSE
    )
    selected <- vector("list", length(genobins_files))
    for (i in seq_along(genobins_files)) {
        sample_id <- names(genobins_files)[i]
        .vmsg(verbose, "Aggregating per-cell S/G1 CP10K: ", sample_id)
        mat <- readRDS(genobins_files[[i]])
        stopifnot(inherits(mat, "sparseMatrix"))
        if (!identical(rownames(mat), bin_ids)) {
            stopifnot(setequal(rownames(mat), bin_ids))
            mat <- mat[bin_ids, , drop = FALSE]
        }
        ids <- colnames(mat)
        if (is.null(ids)) stop("Raw matrix has no cell IDs: ", sample_id)
        ids <- ifelse(startsWith(ids, paste0(sample_id, "-")),
                      ids, paste0(sample_id, "-", ids))
        stopifnot(!anyDuplicated(ids), all(ids %in% rownames(meta_data)))
        s_idx <- which(ids %in% s_cells)
        g1_idx <- which(ids %in% g1_cells)
        if (!length(s_idx) || !length(g1_idx)) {
            stop("Both QC-passed S cells and G1 cells are required: ", sample_id)
        }
        library_size <- Matrix::colSums(mat)
        selected_idx <- c(s_idx, g1_idx)
        invalid <- selected_idx[
            !is.finite(library_size[selected_idx]) | library_size[selected_idx] <= 0
        ]
        if (length(invalid)) {
            stop("Selected cells have zero or non-finite total reads: ",
                 paste(utils::head(ids[invalid], 10L), collapse = ", "))
        }
        # rowSums(mat %*% Diagonal(x = 1e4 / colSums(mat))), evaluated as
        # a sparse matrix-vector product without materializing normalized cells.
        counts$s_counts <- counts$s_counts + as.numeric(
            mat[, s_idx, drop = FALSE] %*% (1e4 / library_size[s_idx])
        )
        counts$g1_counts <- counts$g1_counts + as.numeric(
            mat[, g1_idx, drop = FALSE] %*% (1e4 / library_size[g1_idx])
        )
        selected[[i]] <- data.frame(
            sample = sample_id, n_s = length(s_idx), n_g1 = length(g1_idx)
        )
        rm(mat)
        invisible(gc())
    }
    attr(counts, "normalization") <- "cp10k_per_cell_v1"
    attr(counts, "selected_cells") <- do.call(rbind, selected)
    if (!is.null(cache_file)) {
        saveRDS(list(signature = signature, counts = counts), cache_file)
    }
    counts
}

#' Pseudobulk replication timing as a z-scored S/G1 ratio
#'
#' On the bins that are mappable and have a finite, non-negative S signal and
#' a finite, positive G1 signal, computes
#' `(s / sum(s)) / (g1 / sum(g1))` and z-scores it. Unlike
#' [calculate_pseudobulk_rt()] there is no GC / mappability model, no removal
#' of extreme bins and no smoothing: the G1 pseudobulk is the only baseline.
#'
#' @param counts A `data.frame` with columns `bin_id`, `s_counts` and
#'   `g1_counts`, as returned by [aggregate_pseudobulk_s_g1_counts()].
#'   `bin_id` must be identical to `names(genome_granges)`.
#' @param genome_granges Named `GRanges` with the logical metadata column
#'   `is_mappable`.
#'
#' @return A sorted `GRanges` of the bins used, with metadata columns `score`
#'   (z-score of the ratio), `s_g1_ratio`, `s_counts` and `g1_counts`.
#'
#' @export
calculate_pseudobulk_rt_g1_ratio <- function(counts, genome_granges) {
    stopifnot(
        inherits(genome_granges, "GRanges"),
        identical(as.character(counts$bin_id), names(genome_granges)),
        all(c("s_counts", "g1_counts") %in% colnames(counts)),
        "is_mappable" %in% colnames(S4Vectors::mcols(genome_granges))
    )
    valid <- genome_granges$is_mappable %in% TRUE &
        is.finite(counts$s_counts) & counts$s_counts >= 0 &
        is.finite(counts$g1_counts) & counts$g1_counts > 0
    if (sum(valid) < 2L || sum(counts$s_counts[valid]) <= 0) {
        stop("Insufficient mappable bins with S reads and a positive G1 baseline.")
    }
    ratio <- (counts$s_counts[valid] / sum(counts$s_counts[valid])) /
        (counts$g1_counts[valid] / sum(counts$g1_counts[valid]))
    if (!is.finite(stats::sd(ratio)) || stats::sd(ratio) == 0) {
        stop("S/G1 ratio has no finite variation.")
    }
    gr <- genome_granges[valid]
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
        score = as.numeric(scale(ratio)), s_g1_ratio = ratio,
        s_counts = counts$s_counts[valid], g1_counts = counts$g1_counts[valid]
    )
    BiocGenerics::sort(gr)
}

#' S/G1 ratio pseudobulk for a group of samples, optionally written to bigWig
#'
#' Runs [aggregate_pseudobulk_s_g1_counts()] and
#' [calculate_pseudobulk_rt_g1_ratio()] for one group of samples (for example
#' the replicates of one condition). The result is returned, and a bigWig is
#' written only when `file` is given.
#'
#' @inheritParams aggregate_pseudobulk_s_g1_counts
#' @param genome_granges Named `GRanges` with the metadata column
#'   `is_mappable`, or the path of an RDS file holding one. Its names are the
#'   bin ids. Writing a bigWig needs chromosome lengths in its `seqinfo`.
#' @param file `NULL` (default, nothing is written) or the path of the bigWig
#'   that receives the `score` column (needs the rtracklayer package).
#' @param ... Passed to [aggregate_pseudobulk_s_g1_counts()], e.g.
#'   `cache_file`.
#'
#' @return The `GRanges` of [calculate_pseudobulk_rt_g1_ratio()], with the
#'   per-sample cell numbers in `S4Vectors::metadata(x)$selected_cells`;
#'   returned invisibly when `file` is given.
#'
#' @export
export_pseudobulk_g1_ratio <- function(
    genobins_files, meta_data, genome_granges, file = NULL, ...
) {
    if (is.character(genome_granges)) {
        stopifnot(length(genome_granges) == 1L, file.exists(genome_granges))
        genome_granges <- readRDS(genome_granges)
    }
    stopifnot(
        inherits(genome_granges, "GRanges"),
        !is.null(names(genome_granges))
    )
    if (!is.null(file)) {
        stopifnot(is.character(file), length(file) == 1L)
        rlang::check_installed("rtracklayer", reason = "to write bigWig files.")
    }
    counts <- aggregate_pseudobulk_s_g1_counts(
        genobins_files, meta_data, names(genome_granges), ...
    )
    gr <- calculate_pseudobulk_rt_g1_ratio(counts, genome_granges)
    S4Vectors::metadata(gr)$selected_cells <- attr(counts, "selected_cells")
    if (is.null(file)) return(gr)
    rtracklayer::export(gr[, "score"], file, format = "bigWig")
    invisible(gr)
}
