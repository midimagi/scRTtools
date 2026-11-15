# Peak calling on binned signals and peak annotation with feature tracks.
# find_rt_peaks() is migrated from rcode/processing/find_rt_peaks.R,
# annotate_peak_scores() from rcode/analysis/annotate_peak_scores.R.

#' Find local maxima in a binned genomic signal
#'
#' Within each chromosome, bins are split into blocks at gaps between
#' consecutive starts larger than `max_gap` (default: 1.5 times the median bin
#' step). Runs of equal scores (within `plateau_tol`) form plateaus. A
#' plateau of at most `max_plateau_bins` bins is a peak if it has at least
#' `max(1, ceiling(flank_bins / 2))` non-`NA` bins among the `flank_bins` bins
#' on each side, is higher than all of them by more than `plateau_tol`, is at
#' least `min_score` (unless `NULL`) and exceeds the median of the flanking
#' bins by at least `min_delta`. The middle bin of the plateau is reported.
#'
#' @param granges `GRanges` of bins with a numeric score column (e.g. from
#'   [calculate_pseudobulk_rt()]).
#' @param score_col Name of the score metadata column.
#' @param flank_bins Positive integer; flank size on each side.
#' @param min_score `NULL` or the minimum peak score. We use `0.2` on
#'   RT scores and `10` on track-centre counts.
#' @param min_delta Minimum difference between the peak and the flank median.
#' @param max_plateau_bins Maximum plateau size.
#' @param max_gap `NULL` or the maximum start-to-start distance within a
#'   block.
#' @param plateau_tol Tolerance for equal scores.
#'
#' @return A sorted `GRanges` of peak bins with additional metadata columns
#'   `RT_value` (the score), `local_median`, `local_delta` and
#'   `plateau_size`.
#'
#' @export
find_rt_peaks <- function(
  granges,
  score_col = "score",
  flank_bins = 5,
  min_score = 0.2,
  min_delta = 0.02,
  max_plateau_bins = 5,
  max_gap = NULL,
  plateau_tol = 1e-8
) {
  stopifnot(inherits(granges, "GRanges"))
  if (!score_col %in% colnames(S4Vectors::mcols(granges))) {
    stop("score_col not found in mcols(granges): ", score_col)
  }

  if (!is.null(min_score) && (!is.numeric(min_score) || length(min_score) != 1)) {
    stop("min_score must be NULL or a numeric scalar.")
  }

  if (!is.numeric(flank_bins) || length(flank_bins) != 1 || flank_bins < 1) {
    stop("flank_bins must be a positive integer scalar.")
  }

  flank_bins <- as.integer(flank_bins)
  min_flank_bins_each_side <- max(1L, ceiling(flank_bins / 2))

  granges <- BiocGenerics::sort(granges)
  split_gr <- S4Vectors::split(granges, GenomicRanges::seqnames(granges))

  peak_list <- lapply(split_gr, function(chr_gr) {
    chr_gr <- BiocGenerics::sort(chr_gr)
    y_all <- S4Vectors::mcols(chr_gr)[[score_col]]
    n_all <- length(chr_gr)

    if (n_all == 0) return(chr_gr[FALSE])

    starts <- BiocGenerics::start(chr_gr)

    if (is.null(max_gap)) {
      bin_step <- stats::median(diff(starts), na.rm = TRUE)
      if (!is.finite(bin_step) || bin_step <= 0) {
        bin_step <- stats::median(BiocGenerics::width(chr_gr), na.rm = TRUE)
      }
      max_gap_use <- bin_step * 1.5
    } else {
      max_gap_use <- max_gap
    }

    if (n_all == 1) {
      block_id <- 1L
    } else {
      gap_break <- c(TRUE, diff(starts) > max_gap_use)
      block_id <- cumsum(gap_break)
    }

    keep_all <- rep(FALSE, n_all)
    local_median_all <- rep(NA_real_, n_all)
    local_delta_all <- rep(NA_real_, n_all)
    plateau_size_all <- rep(NA_integer_, n_all)

    block_indices <- split(seq_len(n_all), block_id)

    for (block_idx in block_indices) {
      y <- y_all[block_idx]
      n <- length(y)

      if (n < 2 * min_flank_bins_each_side + 1) next

      i <- 1L

      while (i <= n) {
        if (is.na(y[i])) {
          i <- i + 1L
          next
        }

        # ------------------------------------------------------
        # Identify plateau around current bin
        # ------------------------------------------------------

        left_plateau <- i
        while (
          left_plateau > 1L &&
            !is.na(y[left_plateau - 1L]) &&
            abs(y[left_plateau - 1L] - y[i]) <= plateau_tol
        ) {
          left_plateau <- left_plateau - 1L
        }

        right_plateau <- i
        while (
          right_plateau < n &&
            !is.na(y[right_plateau + 1L]) &&
            abs(y[right_plateau + 1L] - y[i]) <= plateau_tol
        ) {
          right_plateau <- right_plateau + 1L
        }

        plateau_idx <- left_plateau:right_plateau
        plateau_size <- length(plateau_idx)
        plateau_value <- y[i]

        # Move i to the next position after this plateau
        i <- right_plateau + 1L

        if (plateau_size > max_plateau_bins) next

        # ------------------------------------------------------
        # Require enough bins on both sides within same continuous block
        # ------------------------------------------------------

        left_start <- max(1L, left_plateau - flank_bins)
        left_end <- left_plateau - 1L

        right_start <- right_plateau + 1L
        right_end <- min(n, right_plateau + flank_bins)

        if (left_end < left_start || right_end < right_start) next

        left_idx <- left_start:left_end
        right_idx <- right_start:right_end

        left_vals <- y[left_idx]
        right_vals <- y[right_idx]

        n_left_valid <- sum(!is.na(left_vals))
        n_right_valid <- sum(!is.na(right_vals))

        if (
          n_left_valid < min_flank_bins_each_side ||
            n_right_valid < min_flank_bins_each_side
        ) {
          next
        }

        bg_vals <- c(left_vals, right_vals)
        bg_vals <- bg_vals[!is.na(bg_vals)]

        if (!length(bg_vals)) next

        # ------------------------------------------------------
        # Plateau local maximum:
        # plateau value must be strictly higher than surrounding bins
        # ------------------------------------------------------

        is_local_max <- plateau_value > max(bg_vals, na.rm = TRUE) + plateau_tol

        pass_score <- TRUE
        if (!is.null(min_score)) {
          pass_score <- plateau_value >= min_score
        }

        local_median <- stats::median(bg_vals, na.rm = TRUE)
        local_delta <- plateau_value - local_median
        pass_delta <- is.finite(local_delta) && local_delta >= min_delta

        if (is_local_max && pass_score && pass_delta) {
          selected <- plateau_idx[ceiling(length(plateau_idx) / 2)]
          global_selected <- block_idx[selected]

          keep_all[global_selected] <- TRUE
          local_median_all[global_selected] <- local_median
          local_delta_all[global_selected] <- local_delta
          plateau_size_all[global_selected] <- plateau_size
        }
      }
    }

    chr_gr$RT_value <- y_all
    chr_gr$local_median <- local_median_all
    chr_gr$local_delta <- local_delta_all
    chr_gr$plateau_size <- plateau_size_all

    chr_gr[keep_all]
  })

  peak_list <- peak_list[lengths(peak_list) > 0]

  if (!length(peak_list)) {
    out <- granges[FALSE]
    out$RT_value <- numeric(0)
    out$local_median <- numeric(0)
    out$local_delta <- numeric(0)
    out$plateau_size <- integer(0)
    return(out)
  }

  peak_gr <- BiocGenerics::unlist(GenomicRanges::GRangesList(peak_list), use.names = FALSE)
  peak_gr <- BiocGenerics::sort(peak_gr)

  peak_gr
}

#' Annotate regions with scores from feature tracks
#'
#' For each feature file (imported with `rtracklayer::import()`, must have a
#' `score` column), aggregate the scores of the feature ranges overlapping
#' each region: `"weighted_mean"` (overlap-width weighted mean), `"sum"`
#' (sum of score x overlap width) or `"max"`. Regions without overlap get
#' `NA`.
#'
#' @param regions A `GRanges`, or the path of a BED/TSV file with columns
#'   `chr`, `start` (0-based) and `end` (read with `data.table::fread()`).
#' @param features Named character vector (or list) of feature file paths;
#'   names become column names (file basenames are used for unnamed entries,
#'   duplicated names are made unique).
#' @param window `NULL` or a width (bp) to resize each region to, around its
#'   centre (then trimmed).
#' @param method Aggregation method.
#'
#' @return `as.data.frame()` of the regions with one column per feature.
#'
#' @export
annotate_peak_scores <- function(
  regions = NULL,
  features,
  window = NULL,
  method = c("weighted_mean", "sum", "max")
) {

  rlang::check_installed("rtracklayer", reason = "to import feature files.")

  ## -------------------------------
  ## 1. argument check
  ## -------------------------------
  if (is.list(features)) {
    features <- unlist(features, use.names = TRUE)
  }

  if (!is.character(features) || length(features) == 0) {
    stop("features must be a non-empty named character vector of feature files")
  }

  if (!all(file.exists(features))) {
    missing_files <- features[!file.exists(features)]
    stop(
      "Some feature files do not exist: ",
      paste(missing_files, collapse = ", ")
    )
  }

  feature_names <- names(features)

  if (is.null(feature_names)) {
    feature_names <- rep("", length(features))
  }

  empty_names <- is.na(feature_names) | feature_names == ""

  if (any(empty_names)) {
    feature_names[empty_names] <- tools::file_path_sans_ext(
      basename(features[empty_names])
    )
    warning(
      "Some features are unnamed; use file basename as feature column name: ",
      paste(feature_names[empty_names], collapse = ", ")
    )
  }

  method <- match.arg(method)

  ## -------------------------------
  ## 2. read regions -> GRanges
  ## -------------------------------
  if (inherits(regions, "GRanges")) {

    gr <- regions

  } else if (is.character(regions) && length(regions) == 1) {

    if (!file.exists(regions)) {
      stop("regions file not found: ", regions)
    }

    df <- data.table::fread(regions)

    required <- c("chr", "start", "end")
    if (!all(required %in% colnames(df))) {
      stop("BED/TSV file must contain columns: chr, start, end")
    }

    gr <- GenomicRanges::makeGRangesFromDataFrame(
      df,
      seqnames.field = "chr",
      start.field = "start",
      end.field = "end",
      keep.extra.columns = TRUE,
      starts.in.df.are.0based = TRUE
    )

  } else {
    stop("`regions` must be either a BED/TSV file path or a GRanges object")
  }

  ## -------------------------------
  ## 3. resize to window optional
  ## -------------------------------
  if (!is.null(window)) {
    if (
      !is.numeric(window) ||
        length(window) != 1 ||
        !is.finite(window) ||
        window <= 0
    ) {
      stop("window must be a positive number")
    }

    window_int <- as.integer(round(window))

    if (!isTRUE(all.equal(window, window_int))) {
      message("coerce window to integer: ", window_int, " bp")
    } else {
      message("resize window to ", window_int, " bp")
    }

    gr <- GenomicRanges::resize(gr, width = window_int, fix = "center")
    gr <- GenomicRanges::trim(gr)
  }

  ## -------------------------------
  ## 4. annotate each feature
  ## -------------------------------
  for (i in seq_along(features)) {

    feature <- features[[i]]
    raw_name <- feature_names[[i]]

    message("annotate ", raw_name, ": ", feature)

    gr_feat <- rtracklayer::import(feature)

    if (!"score" %in% colnames(S4Vectors::mcols(gr_feat))) {
      stop("feature file must contain 'score' metadata: ", feature)
    }

    existing_cols <- colnames(S4Vectors::mcols(gr))
    if (is.null(existing_cols)) {
      existing_cols <- character(0)
    }

    feat_name <- make.unique(
      c(existing_cols, raw_name),
      sep = "_"
    )[length(existing_cols) + 1]

    if (!identical(raw_name, feat_name)) {
      message(
        "rename duplicated feature column '",
        raw_name,
        "' to '",
        feat_name,
        "'"
      )
    }

    S4Vectors::mcols(gr)[[feat_name]] <- rep(NA_real_, length(gr))

    hits <- GenomicRanges::findOverlaps(
      gr,
      gr_feat,
      ignore.strand = TRUE
    )

    if (length(hits) == 0) {
      warning("No overlaps found for feature: ", feat_name)
      next
    }

    ov <- GenomicRanges::pintersect(
      gr[S4Vectors::queryHits(hits)],
      gr_feat[S4Vectors::subjectHits(hits)]
    )

    w <- IRanges::width(ov)
    s <- gr_feat$score[S4Vectors::subjectHits(hits)]

    dt <- data.table::data.table(
      idx = S4Vectors::queryHits(hits),
      score = s,
      weight = w
    )

    if (method == "weighted_mean") {

      dt_sum <- dt[
        ,
        list(
          value = {
            total_weight <- sum(weight, na.rm = TRUE)
            if (total_weight > 0) {
              sum(score * weight, na.rm = TRUE) / total_weight
            } else {
              NA_real_
            }
          }
        ),
        by = idx
      ]

    } else if (method == "sum") {

      dt_sum <- dt[
        ,
        list(
          value = sum(score * weight, na.rm = TRUE)
        ),
        by = idx
      ]

    } else if (method == "max") {

      dt_sum <- dt[
        ,
        list(
          value = max(score, na.rm = TRUE)
        ),
        by = idx
      ]

    }

    S4Vectors::mcols(gr)[[feat_name]][dt_sum$idx] <- dt_sum$value
  }

  ## -------------------------------
  ## 5. return as data.frame
  ## -------------------------------
  .mcols_df(gr)
}
