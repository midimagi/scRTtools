# Read counting, window filtering, GC correction, smoothing and export.
# Migrated from rcode/processing/calculate_rt.R (identical to
# scRT_learning rcode/tigeR/calculate_rt.R).

#' Count reads per genomic window
#'
#' Count reads from per-chromosome read-position files and add raw read counts
#' and raw copy-number estimates to a `GRangesList`.
#'
#' Reads are assigned to windows with `findInterval()` on the cut points
#' `c(start(windows), end(last window))`, so windows are assumed to be sorted,
#' contiguous and non-overlapping within each chromosome. Reads outside the
#' covered span are dropped. `raw_cn = 2 * raw_counts / mean(raw_counts)`,
#' where the mean is taken over all windows of all chromosomes; values above
#' `raw_cn_th` or equal to 0 are set to `NA`.
#'
#' @param wins_grl A named `GRangesList`; each element holds the windows of one
#'   chromosome, named like `"chr1"`.
#' @param data_dir Directory containing one whitespace-separated read-position
#'   file per chromosome, named `readsChr<N>.txt` where `<N>` is the element
#'   name without the `"chr"` prefix.
#' @param bad_coords Optional named list (names = chromosome names, as in
#'   `wins_grl`) of read positions to drop before counting, or `NULL`.
#' @param raw_cn_th Positive numeric threshold for `raw_cn`.
#' @param verbose Logical; print one progress message per chromosome. Default `FALSE`.
#'
#' @return A list with
#'   \describe{
#'     \item{`wins_count`}{`wins_grl` with metadata columns `raw_counts`
#'       (integer) and `raw_cn` (numeric).}
#'     \item{`reads_list`}{named list of numeric read positions per chromosome
#'       after removing `bad_coords`.}
#'   }
#'
#' @examples
#' d <- tempfile(); dir.create(d)
#' writeLines(as.character(c(5, 15, 15, 25, 35)), file.path(d, "readsChr1.txt"))
#' wins <- GenomicRanges::GRangesList(
#'   chr1 = GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 11, 21, 31), width = 10))
#' )
#' res <- count_reads(wins, d, verbose = FALSE)
#' res$wins_count$chr1
#'
#' @export
count_reads <- function(wins_grl,
                        data_dir,
                        bad_coords = NULL,
                        raw_cn_th = 6,
                        verbose = FALSE) {
  ## ------------------------------------------------------------
  ## Parameter checks
  ## ------------------------------------------------------------
  if (missing(wins_grl) || is.null(wins_grl)) {
    stop("'wins_grl' must be supplied and cannot be NULL.", call. = FALSE)
  }
  if (!methods::is(wins_grl, "GRangesList")) {
    stop("'wins_grl' must be a GenomicRanges::GRangesList object.", call. = FALSE)
  }
  if (length(wins_grl) == 0L) {
    stop("'wins_grl' must contain at least one chromosome/window set.", call. = FALSE)
  }
  if (missing(data_dir) || is.null(data_dir)) {
    stop("'data_dir' must be supplied and cannot be NULL.", call. = FALSE)
  }
  if (!is.character(data_dir) || length(data_dir) != 1L || is.na(data_dir)) {
    stop("'data_dir' must be a single non-NA character string.", call. = FALSE)
  }
  if (!dir.exists(data_dir)) {
    stop("'data_dir' does not exist: ", data_dir, call. = FALSE)
  }
  if (!is.numeric(raw_cn_th) || length(raw_cn_th) != 1L || is.na(raw_cn_th)) {
    stop("'raw_cn_th' must be a single non-NA numeric value.", call. = FALSE)
  }
  if (!is.finite(raw_cn_th) || raw_cn_th <= 0) {
    stop("'raw_cn_th' must be a finite positive numeric value.", call. = FALSE)
  }
  if (!is.null(bad_coords)) {
    if (!is.list(bad_coords)) {
      stop("'bad_coords' must be a list.", call. = FALSE)
    }
    if (length(bad_coords) != length(wins_grl)) {
      stop("'bad_coords' must have the same length as 'wins_grl'. ")
    }
  }

  # NOTE (kept from original): this pre-check builds file names from the
  # element *index* (readsChr1..N), while the loop below uses the element
  # *name* (e.g. readsChrX). They disagree for non-numeric or non-consecutive
  # chromosomes. Reported as a suspected bug.
  read_files <- file.path(data_dir, paste0("readsChr", seq_along(wins_grl), ".txt"))
  missing_files <- read_files[!file.exists(read_files)]

  if (length(missing_files) > 0L) {
    stop(
      "The following read-position files are missing:\n",
      paste(missing_files, collapse = "\n"),
      call. = FALSE
    )
  }

  n_chr <- length(wins_grl)
  chrs <- names(wins_grl)
  reads_list <- stats::setNames(vector("list", n_chr), chrs)

  ## Count reads per window
  for (chr in chrs) {
    .vmsg(verbose, "processing ", chr)

    chr_num <- gsub("chr", "", chr)
    filename <- file.path(data_dir, paste0("readsChr", chr_num, ".txt"))
    if (!file.exists(filename)) {
      stop("File not found: ", filename)
    }

    reads_chr <- scan(filename, what = numeric(), quiet = TRUE)
    reads_chr <- reads_chr[!(reads_chr %in% bad_coords[[chr]])]

    reads_list[[chr]] <- reads_chr

    gr_chr <- wins_grl[[chr]]
    if (length(gr_chr) == 0L) {
      S4Vectors::mcols(gr_chr)$raw_counts <- integer(0)
      wins_grl[[chr]] <- gr_chr
      next
    }

    wins_chr_cut <- c(
      BiocGenerics::start(gr_chr),
      BiocGenerics::end(gr_chr)[length(gr_chr)]
    )

    n_bins <- length(wins_chr_cut) - 1L

    bin_per_read <- findInterval(
      reads_chr,
      wins_chr_cut,
      rightmost.closed = TRUE,
      all.inside = FALSE
    )

    valid <- bin_per_read >= 1L & bin_per_read <= n_bins
    bin_valid <- bin_per_read[valid]
    raw_counts <- tabulate(bin_valid, nbins = n_bins)
    S4Vectors::mcols(gr_chr)$raw_counts <- raw_counts
    wins_grl[[chr]] <- gr_chr
  }

  ## Mean coverage based on raw_counts
  mean_coverage <- mean(
    unlist(
      lapply(chrs, function(chr) {
        S4Vectors::mcols(wins_grl[[chr]])$raw_counts
      }),
      use.names = FALSE
    ),
    na.rm = TRUE
  )

  ## Calculate raw copy number
  for (chr in chrs) {
    gr_chr <- wins_grl[[chr]]
    if (length(gr_chr) == 0L) {
      S4Vectors::mcols(gr_chr)$raw_cn <- numeric(0)
      wins_grl[[chr]] <- gr_chr
      next
    }
    raw_counts <- S4Vectors::mcols(gr_chr)$raw_counts
    raw_cn <- raw_counts / mean_coverage * 2
    raw_cn[raw_cn > raw_cn_th | raw_cn == 0] <- NA_real_
    S4Vectors::mcols(gr_chr)$raw_cn <- raw_cn
    wins_grl[[chr]] <- gr_chr
  }

  list(wins_count = wins_grl,
       reads_list = reads_list)
}

#' Segment and filter copy-number windows
#'
#' For each chromosome, split the windows at genome gaps, segment the values of
#' `mcols_id` within each gap-free fragment with `changepoint::cpt.mean()`, and
#' replace each window value by its segment median (`seg_score`). Windows are
#' then removed if any of the following holds:
#' \itemize{
#'   \item the window is more than twice as long as the shortest window of the
#'     chromosome (`end - start` used as length);
#'   \item `|seg_score - genome median| > std_th * genome SD` (median and SD of
#'     the raw `mcols_id` values over all chromosomes);
#'   \item `|seg_score - chromosome median of retained seg_score| > std_th * genome SD`;
#'   \item `seg_score == 0` or `seg_score` is `NA`.
#' }
#' Only windows with non-`NA` values `<= 20` are segmented, and only fragments
#' with more than 20 such windows; all other windows get `NA` and are removed.
#' If `cpt.mean()` fails on a fragment, all its windows get the fragment median.
#'
#' @param wins_grl A `GRangesList` (one element per chromosome) containing the
#'   numeric metadata column `mcols_id`.
#' @param gaps A `GRanges` of genome gap regions.
#' @param mcols_id Character scalar. Metadata column to segment, e.g.
#'   `"raw_cn"`.
#' @param std_th Positive numeric scalar. SD multiplier for filtering. We
#'   use `0.5` (sometimes `0.35`).
#' @param cpt_arg List of arguments passed to `changepoint::cpt.mean()`.
#'
#' @return A `GRangesList` with the same names as `wins_grl` holding only the
#'   retained windows; metadata columns are inherited, plus a numeric column
#'   named `paste0(mcols_id, "_seg_score")`.
#'
#' @examples
#' set.seed(1)
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(seq(1, by = 1000, length.out = 60),
#'                                                       width = 1000))
#' gr$raw_cn <- c(rnorm(30, 2, 0.1), rnorm(30, 4, 0.1))
#' gaps <- GenomicRanges::GRanges("chr1", IRanges::IRanges(10001, 12000))
#' res <- filter_windows(GenomicRanges::GRangesList(chr1 = gr), gaps, "raw_cn",
#'                       std_th = 2)
#' length(res$chr1)
#'
#' @export
filter_windows <- function(wins_grl,
                           gaps,
                           mcols_id,
                           std_th = 2.0,
                           cpt_arg = list(method = "PELT",
                                          penalty = "MBIC",
                                          test.stat = "Normal",
                                          minseglen = 5,
                                          class = TRUE)) {
  ## ------------------------------------------------------------
  ## Parameter checks
  ## ------------------------------------------------------------
  if (missing(wins_grl)) {
    stop("'wins_grl' must be a GenomicRanges::GRangesList object.", call. = FALSE)
  }
  .validate_grl(wins_grl, "wins_grl")
  if (missing(gaps) || !methods::is(gaps, "GRanges")) {
    stop("'gaps' must be a GenomicRanges::GRanges object.", call. = FALSE)
  }
  if (missing(mcols_id) ||
      !is.character(mcols_id) ||
      length(mcols_id) != 1L ||
      is.na(mcols_id) ||
      !nzchar(mcols_id)) {
    stop("'mcols_id' must be a single non-empty character string.", call. = FALSE)
  }
  .validate_grl(wins_grl, "wins_grl", mcol = mcols_id)
  if (!is.numeric(std_th) ||
      length(std_th) != 1L ||
      is.na(std_th) ||
      !is.finite(std_th) ||
      std_th <= 0) {
    stop("'std_th' must be a single finite positive numeric value.", call. = FALSE)
  }
  if (!is.list(cpt_arg)) {
    stop("'cpt_arg' must be a list.", call. = FALSE)
  }

  mcols_id_all <- unlist(lapply(seq_along(wins_grl), function(i) {
    S4Vectors::mcols(wins_grl[[i]])[[mcols_id]]
  }), use.names = FALSE)

  depth_mean <- stats::median(mcols_id_all, na.rm = TRUE)
  depth_std <- stats::sd(mcols_id_all, na.rm = TRUE)
  std_th <- as.numeric(std_th)

  .seg_fun <- function(data_in) {
    n <- length(data_in)
    if (n == 0L) return(numeric(0))
    if (all(is.na(data_in))) return(rep(NA_real_, n))

    data_out <- rep(NA_real_, n)
    keep <- which(!is.na(data_in))

    if (length(keep) <= 1L) {
      data_out[keep] <- data_in[keep]
      return(data_out)
    }

    y <- data_in[keep]
    arg <- c(list(data = y), cpt_arg)

    cpt_res <- tryCatch(
      do.call(changepoint::cpt.mean, arg),
      error = function(e) NULL
    )

    if (is.null(cpt_res)) {
      data_out[keep] <- stats::median(y, na.rm = TRUE)
      return(data_out)
    }

    cps <- changepoint::cpts(cpt_res)
    bounds <- c(0L, cps, length(y))
    seg_median <- rep(NA_real_, length(y))

    for (i in seq_len(length(bounds) - 1L)) {
      idx <- seq.int(bounds[i] + 1L, bounds[i + 1L])
      seg_median[idx] <- stats::median(y[idx], na.rm = TRUE)
    }

    data_out[keep] <- seg_median
    data_out
  }

  mcols_id_seg_score <- paste0(mcols_id, "_seg_score")
  wins_filtered <- GenomicRanges::GRangesList(
    lapply(seq_along(wins_grl), function(Chr) {
      gr <- wins_grl[[Chr]]
      n <- length(gr)

      if (n == 0L) {
        S4Vectors::mcols(gr)[[mcols_id_seg_score]] <- numeric(0)
        return(gr)
      }

      chr_name <- as.character(GenomeInfoDb::seqnames(gr)[1])
      mcols_id_values <- S4Vectors::mcols(gr)[[mcols_id]]
      seg_score <- rep(NA_real_, n)

      window_length <- GenomicRanges::end(gr) - GenomicRanges::start(gr)
      min_window_length <- min(window_length, na.rm = TRUE)
      filter_long_window <- which(window_length > min_window_length * 2)

      gap_chr <- gaps[as.character(GenomeInfoDb::seqnames(gaps)) == chr_name]
      gap_chr <- BiocGenerics::sort(gap_chr)

      chr_start <- GenomicRanges::start(gr)[1L] - 1L
      chr_end <- GenomicRanges::end(gr)[n]

      if (length(gap_chr) > 0L) {
        G <- cbind(
          c(chr_start, GenomicRanges::end(gap_chr)),
          c(GenomicRanges::start(gap_chr), chr_end)
        )
        G <- G[G[, 1L] < G[, 2L], , drop = FALSE]
      } else {
        G <- matrix(c(chr_start, chr_end), ncol = 2L, byrow = TRUE)
      }

      for (frag in seq_len(nrow(G))) {
        # NOTE (kept from original): hard-coded value cap (<= 20) and minimum
        # fragment size (> 20 windows).
        used_wins <- which(
          GenomicRanges::start(gr) > G[frag, 1L] &
            GenomicRanges::end(gr) <= G[frag, 2L] &
            !is.na(mcols_id_values) &
            mcols_id_values <= 20
        )

        if (length(used_wins) > 20L) {
          seg_score[used_wins] <- .seg_fun(mcols_id_values[used_wins])
        }
      }

      filter_2 <- which(abs(seg_score - depth_mean) > std_th * depth_std)
      remain <- setdiff(seq_len(n), filter_2)
      chr_median <- stats::median(seg_score[remain], na.rm = TRUE)

      filter_3 <- which(abs(seg_score - chr_median) > std_th * depth_std)
      filter_4 <- which(seg_score == 0)
      filter_5 <- which(is.na(seg_score))

      filtered_all <- sort(unique(c(
        filter_long_window,
        filter_2,
        filter_3,
        filter_4,
        filter_5
      )))

      keep_wins <- setdiff(seq_len(n), filtered_all)
      S4Vectors::mcols(gr)[[mcols_id_seg_score]] <- seg_score

      gr[keep_wins]
    })
  )

  names(wins_filtered) <- names(wins_grl)
  wins_filtered
}

#' GC-correct read-count windows for one sample
#'
#' Estimate per-GC-bin bias factors from the reads that fall inside the
#' retained windows `wins_filtered`, then compute expected counts and the
#' observed/expected ratio for every retained window.
#'
#' Reads are assigned greedily to the 401 GC bins (GC 0.000-0.400 step 0.001
#' bins as produced by the tigeR genome processing), visiting the GC bins with
#' the most retained positions first. The bias factor of a GC bin is
#' `(reads fraction) / (genome fraction)`. The expected count of a window is
#' `sum_{gc in gc_bins_use} bases_in_gc_bin * factor * median(raw_counts) / win_size / 2`,
#' where `win_size` is the median of `row_sum` in `wins_gc`.
#'
#' @param wins_filtered Named `GRangesList` of read-count windows with numeric
#'   metadata column `raw_counts` (e.g. output of [filter_windows()]).
#'   **Ranges must be named**: windows are matched to `wins_gc` by
#'   `names()`.
#' @param wins_gc Named `GRangesList` (same chromosome names) with named ranges
#'   and metadata columns `GC_bin_000` ... `GC_bin_400` (bases per GC bin) and
#'   `row_sum`.
#' @param reads_list Named list of numeric read-position vectors, one per
#'   chromosome (e.g. `count_reads()$reads_list`).
#' @param coords_list Named character vector of RDS file paths, one per
#'   chromosome; each RDS holds a list of 401 numeric vectors of genome
#'   positions per GC bin.
#' @param gc_bins_use Integer indices (1-401) of GC bins used for correction.
#' @param verbose Logical; print progress messages. Default `FALSE`.
#'
#' @return A list with
#'   \describe{
#'     \item{`wins_ratio`}{`wins_filtered` with added numeric metadata columns
#'       `expected_counts` and `oe_ratio`. `oe_ratio` is `NA` where
#'       `row_sum != win_size` or the expected count is 0.}
#'     \item{`gc_factors_per_bin`}{named numeric vector (length 401) of GC
#'       bias factors, `NA` for GC bins without genome coverage.}
#'   }
#'
#' @export
calculate_OE_ratio <- function(wins_filtered,
                               wins_gc,
                               reads_list,
                               coords_list,
                               gc_bins_use = 81:321,
                               verbose = FALSE) {
  ## ------------------------------------------------------------
  ## Parameter checks
  ## ------------------------------------------------------------
  if (!methods::is(wins_filtered, "GRangesList")) {
    stop("'wins_filtered' must be a GRangesList.", call. = FALSE)
  }
  if (!methods::is(wins_gc, "GRangesList")) {
    stop("'wins_gc' must be a GRangesList.", call. = FALSE)
  }
  if (!is.list(reads_list) || is.null(names(reads_list))) {
    stop("'reads_list' must be a named list.", call. = FALSE)
  }
  if (!is.character(coords_list) || is.null(names(coords_list))) {
    stop("'coords_list' must be a named vector of RDS file paths.", call. = FALSE)
  }

  chrs <- names(wins_filtered)
  n_chr <- length(chrs)

  if (length(wins_filtered) == 0L ||
      is.null(chrs) || anyNA(chrs) || any(chrs == "")) {
    stop("'wins_filtered' must be a non-empty named GRangesList.", call. = FALSE)
  }

  if (!all(chrs %in% names(wins_gc))) {
    stop("'wins_gc' must contain all chromosomes in names(wins_filtered).", call. = FALSE)
  }
  if (!all(chrs %in% names(reads_list))) {
    stop("'reads_list' must contain all chromosomes in names(wins_filtered).", call. = FALSE)
  }
  if (!all(chrs %in% names(coords_list))) {
    stop("'coords_list' must contain all chromosomes in names(wins_filtered).", call. = FALSE)
  }

  if (!all(vapply(reads_list[chrs], is.numeric, logical(1)))) {
    stop("'reads_list' elements must be numeric read-position vectors.", call. = FALSE)
  }

  if (!all(vapply(coords_list[chrs], function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && file.exists(x)
  }, logical(1)))) {
    stop("'coords_list' elements must be existing RDS file paths.", call. = FALSE)
  }

  if (!is.numeric(gc_bins_use) ||
      anyNA(gc_bins_use) ||
      any(gc_bins_use != as.integer(gc_bins_use)) ||
      any(gc_bins_use < 1L | gc_bins_use > 401L)) {
    stop("'gc_bins_use' must contain integer indices between 1 and 401.", call. = FALSE)
  }

  gc_cols <- sprintf("GC_bin_%03d", 0:400)

  has_raw_counts <- vapply(chrs, function(chr) {
    gr <- wins_filtered[[chr]]
    "raw_counts" %in% colnames(S4Vectors::mcols(gr)) &&
      is.numeric(S4Vectors::mcols(gr)$raw_counts)
  }, logical(1))

  if (!all(has_raw_counts)) {
    stop("'wins_filtered' must contain numeric metadata column 'raw_counts' for all chromosomes.", call. = FALSE)
  }

  has_gc_cols <- vapply(chrs, function(chr) {
    all(c(gc_cols, "row_sum") %in% colnames(S4Vectors::mcols(wins_gc[[chr]])))
  }, logical(1))

  if (!all(has_gc_cols)) {
    stop("'wins_gc' must contain metadata columns GC_bin_000 to GC_bin_400 and 'row_sum'.", call. = FALSE)
  }

  row_sum_all <- unlist(lapply(chrs, function(chr) {
    S4Vectors::mcols(wins_gc[[chr]])$row_sum
  }), use.names = FALSE)
  win_size <- stats::median(row_sum_all, na.rm = TRUE)

  reads_filtered <- reads_list
  total_bps_per_gc_bin <- numeric(401)
  n_reads_gc_bins_by_chr <- matrix(0, nrow = 401, ncol = n_chr,
                                   dimnames = list(sprintf("GC_bin_%03d", 0:400), chrs))

  ## 1. Estimate GC factors from retained ranges
  for (i in seq_along(chrs)) {
    chr <- chrs[i]
    .vmsg(verbose, "finding reads in GC bins, ", chr)

    reads_chr <- reads_filtered[[chr]]
    filt_gr <- wins_filtered[[chr]]

    if (length(filt_gr) > 0L && length(reads_chr) > 0L) {
      in_filter <- IRanges::overlapsAny(
        IRanges::IRanges(reads_chr, reads_chr),
        IRanges::ranges(filt_gr)
      )
      reads_chr <- reads_chr[in_filter]
    } else {
      reads_chr <- numeric(0)
    }

    reads_filtered[[chr]] <- reads_chr

    coords <- readRDS(coords_list[[chr]])
    if (length(coords) != 401L) {
      stop("`coords_list[[", chr, "]]` should contain 401 GC-bin vectors.")
    }

    if (length(filt_gr) > 0L) {
      coords_filtered <- lapply(coords, function(x) {
        if (length(x) == 0L) return(x)
        keep <- IRanges::overlapsAny(
          IRanges::IRanges(x, x),
          IRanges::ranges(filt_gr)
        )
        x[keep]
      })
    } else {
      coords_filtered <- lapply(coords, function(x) x[0])
    }

    bin_len <- vapply(coords_filtered, length, integer(1))
    bin_order <- order(bin_len, decreasing = TRUE)

    for (j in bin_order) {
      if (length(reads_chr) == 0L || length(coords_filtered[[j]]) == 0L) next
      in_bin <- reads_chr %in% coords_filtered[[j]]
      n_reads_gc_bins_by_chr[j, chr] <- sum(in_bin)
      reads_chr <- reads_chr[!in_bin]
    }

    total_bps_per_gc_bin <- total_bps_per_gc_bin + bin_len
  }

  n_reads_per_gc_bin <- rowSums(n_reads_gc_bins_by_chr, na.rm = TRUE)
  genome_fraction_per_gc_bin <- total_bps_per_gc_bin / sum(total_bps_per_gc_bin, na.rm = TRUE)

  gc_factors_per_bin <- rep(NA_real_, 401)
  names(gc_factors_per_bin) <- sprintf("GC_bin_%03d", 0:400)

  ok <- genome_fraction_per_gc_bin > 0 &
    !is.na(genome_fraction_per_gc_bin) &
    !is.na(n_reads_per_gc_bin)

  total_reads_ok <- sum(n_reads_per_gc_bin[ok], na.rm = TRUE)
  total_fraction_ok <- sum(genome_fraction_per_gc_bin[ok], na.rm = TRUE)

  if (total_reads_ok > 0 && total_fraction_ok > 0) {
    reads_fractions_ok <- n_reads_per_gc_bin[ok] / total_reads_ok
    genome_fractions_ok <- genome_fraction_per_gc_bin[ok] / total_fraction_ok
    gc_factors_per_bin[ok] <- reads_fractions_ok / genome_fractions_ok
  }

  ## 2. Library coverage correction factor
  wins_raw_count <- unlist(lapply(chrs, function(chr) {
    S4Vectors::mcols(wins_filtered[[chr]])$raw_counts
  }), use.names = FALSE)

  expected_reads_per_base <- stats::median(wins_raw_count, na.rm = TRUE) / win_size / 2

  ## 3. Expected reads and normalized counts
  wins_ratio <- wins_filtered
  for (i in seq_along(chrs)) {
    chr <- chrs[i]
    .vmsg(verbose, "calculate observed/expected, ", chr)

    gr <- wins_ratio[[chr]]
    wins_raw_counts <- S4Vectors::mcols(gr)$raw_counts

    gr_gc <- wins_gc[[chr]]
    keep <- names(gr_gc) %in% names(gr)
    gr_gc <- gr_gc[keep]

    wins_by_gc_bins <- as.matrix(S4Vectors::mcols(gr_gc)[, sprintf("GC_bin_%03d", 0:400)])
    wins_expected_nominal <- matrix(NA_real_, nrow = nrow(wins_by_gc_bins), ncol = 401)
    wins_expected_nominal[, gc_bins_use] <- sweep(
      wins_by_gc_bins[, gc_bins_use, drop = FALSE],
      2,
      gc_factors_per_bin[gc_bins_use],
      `*`
    )

    wins_expected_bases <- rowSums(wins_expected_nominal[, gc_bins_use, drop = FALSE], na.rm = TRUE)
    wins_expected_counts <- wins_expected_bases * expected_reads_per_base
    wins_oe_ratio <- wins_raw_counts / wins_expected_counts

    in_bad_unexpected <- S4Vectors::mcols(gr_gc)$row_sum != win_size
    in_zero_expected <- wins_expected_counts == 0
    wins_oe_ratio[in_bad_unexpected | in_zero_expected] <- NA_real_

    S4Vectors::mcols(gr)[["expected_counts"]] <- wins_expected_counts
    S4Vectors::mcols(gr)[["oe_ratio"]] <- wins_oe_ratio
    wins_ratio[[chr]] <- gr
  }

  list(wins_ratio = wins_ratio,
       gc_factors_per_bin = gc_factors_per_bin)
}

#' Smooth window values and compute z-scores
#'
#' Smooth one metadata column of a named `GRangesList` without crossing genome
#' gaps or large data gaps. Within each gap-free fragment holding at least
#' `min_wins` non-`NA` windows, values are first averaged with a centred
#' rolling mean of width 10 (`zoo::rollapply(..., partial = TRUE)`), then a
#' smoothing spline is fitted to the rolling means against window midpoints.
#' Genome-wide z-scores of both the rolling mean and the spline fit are added.
#'
#' Only genome gaps of width `>= gap_th` are used as breaks; a data gap is
#' declared between consecutive non-`NA` windows whose midpoints are
#' `>= gap_th` apart.
#'
#' @param wins_grl Named `GRangesList` of genomic windows (names must equal the
#'   seqnames, e.g. `"chr1"`).
#' @param mcols_id Character scalar. Numeric metadata column to smooth, e.g.
#'   `"raw_cn"`.
#' @param gaps `GRanges` of genome gap regions.
#' @param min_wins Minimum number of non-`NA` windows required to smooth a
#'   fragment (integer `>= 2`).
#' @param gap_th Minimum gap size (bp) used to split smoothing fragments.
#' @param smooth_arg List of extra arguments for `stats::smooth.spline()`,
#'   e.g. `list(spar = 0.6)`, `list(cv = TRUE)` or `list(df = 10)`. We
#'   use `list(spar = 0.1)`.
#'
#' @return `wins_grl` with four added numeric metadata columns:
#'   `{mcols_id}_ave` (rolling mean), `{mcols_id}_smoothed` (spline fit),
#'   `{mcols_id}_ave_z` and `{mcols_id}_smoothed_z` (genome-wide z-scores;
#'   all `NA` if the SD is `NA` or 0). Windows outside smoothed fragments are
#'   `NA`.
#'
#' @examples
#' set.seed(1)
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(seq(1, by = 1000, length.out = 50),
#'                                                       width = 1000))
#' gr$raw_cn <- rnorm(50, 2, 0.2)
#' gaps <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 1))
#' res <- smooth_windows(GenomicRanges::GRangesList(chr1 = gr), "raw_cn", gaps)
#' head(res$chr1)
#'
#' @export
smooth_windows <- function(wins_grl,
                           mcols_id,
                           gaps,
                           min_wins = 20L,
                           gap_th = 100000L,
                           smooth_arg = list(spar = 0.2)) {
  ## ------------------------------------------------------------
  ## Parameter checks
  ## ------------------------------------------------------------
  if (!methods::is(wins_grl, "GRangesList")) {
    stop("'wins_grl' must be a GRangesList.", call. = FALSE)
  }
  if (!methods::is(gaps, "GRanges")) {
    stop("'gaps' must be a GRanges object.", call. = FALSE)
  }
  if (!is.character(mcols_id) ||
      length(mcols_id) != 1L ||
      is.na(mcols_id) ||
      !nzchar(mcols_id)) {
    stop("'mcols_id' must be a non-empty character scalar.", call. = FALSE)
  }

  .validate_grl(wins_grl, "wins_grl", named = TRUE)
  chrs <- names(wins_grl)

  if (!is.numeric(min_wins) ||
      length(min_wins) != 1L ||
      is.na(min_wins) ||
      !is.finite(min_wins) ||
      min_wins < 2L) {
    stop("'min_wins' must be a numeric scalar >= 2.", call. = FALSE)
  }
  min_wins <- as.integer(min_wins)

  if (!is.numeric(gap_th) ||
      length(gap_th) != 1L ||
      is.na(gap_th) ||
      !is.finite(gap_th) ||
      gap_th <= 0) {
    stop("'gap_th' must be a positive numeric scalar.", call. = FALSE)
  }

  if (!is.list(smooth_arg)) {
    stop("'smooth_arg' must be a list of arguments passed to stats::smooth.spline().", call. = FALSE)
  }

  .validate_grl(wins_grl, "wins_grl", mcol = mcols_id, numeric_mcol = TRUE)

  ## ------------------------------------------------------------
  ## Run
  ## ------------------------------------------------------------
  mcols_id_ave <- paste0(mcols_id, "_ave")
  mcols_id_ave_z <- paste0(mcols_id, "_ave_z")
  mcols_id_smoothed <- paste0(mcols_id, "_smoothed")
  mcols_id_smoothed_z <- paste0(mcols_id, "_smoothed_z")

  wins_out <- wins_grl

  ## Keep only genome gaps >= gap_th, using the same threshold as data gaps.
  gaps <- gaps[GenomicRanges::width(gaps) >= gap_th]

  for (chr in chrs) {
    gr <- wins_out[[chr]]

    y <- S4Vectors::mcols(gr)[[mcols_id]]
    y_ave <- rep(NA_real_, length(y))
    y_smoothed <- rep(NA_real_, length(y))

    if (length(gr) == 0L) {
      S4Vectors::mcols(gr)[[mcols_id_ave]] <- y_ave
      S4Vectors::mcols(gr)[[mcols_id_smoothed]] <- y_smoothed
      wins_out[[chr]] <- gr
      next
    }

    win_start <- GenomicRanges::start(gr)
    win_end <- GenomicRanges::end(gr)
    win_mid <- floor((win_start + win_end) / 2)

    ## 1. Large data gaps between non-NA windows
    idx_not_nan <- which(!is.na(y))

    data_gaps_chr <- GenomicRanges::GRanges()

    if (length(idx_not_nan) >= 2L) {
      data_gap_idx <- which(diff(win_mid[idx_not_nan]) >= gap_th)

      if (length(data_gap_idx) > 0L) {
        data_gap_start <- win_end[idx_not_nan[data_gap_idx]] + 1L
        data_gap_end <- win_start[idx_not_nan[data_gap_idx + 1L]] - 1L

        keep_data_gap <- data_gap_start <= data_gap_end

        if (any(keep_data_gap)) {
          data_gaps_chr <- GenomicRanges::GRanges(
            seqnames = chr,
            ranges = IRanges::IRanges(
              start = data_gap_start[keep_data_gap],
              end = data_gap_end[keep_data_gap]))
        }
      }
    }

    ## 2. Genome gaps for this chromosome
    genome_gaps_chr <- gaps[as.character(GenomicRanges::seqnames(gaps)) == chr]
    break_gr <- c(genome_gaps_chr, data_gaps_chr)

    chr_start <- min(win_start, na.rm = TRUE)
    chr_end <- max(win_end, na.rm = TRUE)

    if (length(break_gr) > 0L) {
      break_ir <- IRanges::reduce(GenomicRanges::ranges(break_gr))
      frag_ir <- IRanges::setdiff(
        IRanges::IRanges(start = chr_start, end = chr_end),
        break_ir
      )
    } else {
      frag_ir <- IRanges::IRanges(start = chr_start, end = chr_end)
    }

    ## 3. Smooth each fragment separately
    for (frag in seq_along(frag_ir)) {
      frag_start <- IRanges::start(frag_ir)[frag]
      frag_end <- IRanges::end(frag_ir)[frag]

      # NOTE (kept from original): 'win_start > frag_start' excludes the
      # first window of a chromosome / fragment whose start equals frag_start.
      wins_in_frag <- which(
        win_start > frag_start &
          win_end <= frag_end &
          !is.na(y)
      )

      if (length(wins_in_frag) >= min_wins) {
        x_in <- win_mid[wins_in_frag]
        y_in <- y[wins_in_frag]

        y_in_ave <- zoo::rollapply(data = y_in, width = 10, FUN = mean,
                                   fill = NA_real_, align = "center", partial = TRUE)
        y_ave[wins_in_frag] <- y_in_ave

        fit <- tryCatch(
          do.call(stats::smooth.spline, c(list(x = x_in, y = y_in_ave), smooth_arg)),
          error = function(e) NULL
        )

        if (!is.null(fit)) {
          y_smoothed[wins_in_frag] <- stats::predict(fit, x = x_in)$y
        }
      }
    }
    S4Vectors::mcols(gr)[[mcols_id_ave]] <- y_ave
    S4Vectors::mcols(gr)[[mcols_id_smoothed]] <- y_smoothed
    wins_out[[chr]] <- gr
  }

  ## 4. Z-score normalization
  ave_all <- unlist(lapply(chrs, function(chr) {
    S4Vectors::mcols(wins_out[[chr]])[[mcols_id_ave]]
  }), use.names = FALSE)
  ave_mean <- mean(ave_all, na.rm = TRUE)
  ave_sd <- stats::sd(ave_all, na.rm = TRUE)

  smoothed_all <- unlist(lapply(chrs, function(chr) {
    S4Vectors::mcols(wins_out[[chr]])[[mcols_id_smoothed]]
  }), use.names = FALSE)
  smoothed_mean <- mean(smoothed_all, na.rm = TRUE)
  smoothed_sd <- stats::sd(smoothed_all, na.rm = TRUE)

  for (chr in chrs) {
    gr <- wins_out[[chr]]

    y_ave <- S4Vectors::mcols(gr)[[mcols_id_ave]]
    if (is.na(ave_sd) || ave_sd == 0) {
      ave_z <- rep(NA_real_, length(y_ave))
    } else {
      ave_z <- (y_ave - ave_mean) / ave_sd
    }

    y_smoothed <- S4Vectors::mcols(gr)[[mcols_id_smoothed]]
    if (is.na(smoothed_sd) || smoothed_sd == 0) {
      smoothed_z <- rep(NA_real_, length(y_smoothed))
    } else {
      smoothed_z <- (y_smoothed - smoothed_mean) / smoothed_sd
    }

    S4Vectors::mcols(gr)[[mcols_id_ave_z]] <- ave_z
    S4Vectors::mcols(gr)[[mcols_id_smoothed_z]] <- smoothed_z
    wins_out[[chr]] <- gr
  }

  wins_out
}

#' Export a GRangesList metadata column to BigWig
#'
#' Unlist a `GRangesList`, use one metadata column as the BigWig score, drop
#' `NA` values, attach the supplied `Seqinfo`, sort, and write a BigWig file
#' with \pkg{rtracklayer}.
#'
#' @param grl A `GRangesList`.
#' @param mcol_name Character scalar. Metadata column to export (coerced with
#'   `as.numeric()`).
#' @param file_path Character scalar. Output BigWig path; its directory must
#'   exist.
#' @param seqinf A `Seqinfo` with lengths for every chromosome present in
#'   `grl`.
#'
#' @return Invisibly, `file_path`.
#'
#' @examples
#' if (requireNamespace("rtracklayer", quietly = TRUE)) {
#'   si <- GenomeInfoDb::Seqinfo("chr1", 10000L)
#'   gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 1001), width = 1000),
#'                                val = c(1.5, NA))
#'   f <- tempfile(fileext = ".bw")
#'   export_mcol(GenomicRanges::GRangesList(chr1 = gr), "val", f, si)
#' }
#'
#' @export
export_mcol <- function(grl,
                        mcol_name,
                        file_path,
                        seqinf) {
  rlang::check_installed("rtracklayer", reason = "to write BigWig files.")

  if (!methods::is(grl, "GRangesList")) {
    stop("'grl' must be a GRangesList.", call. = FALSE)
  }
  if (length(grl) == 0L) {
    stop("'grl' must not be empty.", call. = FALSE)
  }

  if (!is.character(mcol_name) ||
      length(mcol_name) != 1L ||
      is.na(mcol_name) ||
      !nzchar(mcol_name)) {
    stop("'mcol_name' must be a non-empty character scalar.", call. = FALSE)
  }

  if (!is.character(file_path) ||
      length(file_path) != 1L ||
      is.na(file_path) ||
      !nzchar(file_path)) {
    stop("'file_path' must be a non-empty character scalar.", call. = FALSE)
  }

  out_dir <- dirname(file_path)
  if (!dir.exists(out_dir)) {
    stop("Output directory does not exist: ", out_dir, call. = FALSE)
  }

  if (!methods::is(seqinf, "Seqinfo")) {
    stop("'seqinf' must be a Seqinfo object.", call. = FALSE)
  }

  gr <- BiocGenerics::unlist(grl, use.names = FALSE)

  if (length(gr) == 0L) {
    stop("'grl' contains no ranges.", call. = FALSE)
  }

  if (!mcol_name %in% colnames(S4Vectors::mcols(gr))) {
    stop("Column '", mcol_name, "' not found in mcols(grl).", call. = FALSE)
  }

  score <- suppressWarnings(as.numeric(S4Vectors::mcols(gr)[[mcol_name]]))
  keep <- !is.na(score)
  gr <- gr[keep]
  score <- score[keep]
  if (length(gr) == 0L) {
    stop("No non-NA numeric values found in column '", mcol_name, "'.", call. = FALSE)
  }

  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(score = score)
  used_chrs <- as.character(unique(GenomeInfoDb::seqnames(gr)))
  missing_chrs <- setdiff(used_chrs, GenomeInfoDb::seqlevels(seqinf))
  if (length(missing_chrs) > 0L) {
    stop(
      "Chromosomes in 'gr' but not in 'seqinf': ",
      paste(missing_chrs, collapse = ", "),
      call. = FALSE
    )
  }

  seq_len_used <- GenomeInfoDb::seqlengths(seqinf)[used_chrs]
  if (anyNA(seq_len_used)) {
    stop(
      "Missing seqlengths for: ",
      paste(used_chrs[is.na(seq_len_used)], collapse = ", "),
      call. = FALSE
    )
  }

  GenomeInfoDb::seqlevels(gr, pruning.mode = "coarse") <- used_chrs
  GenomeInfoDb::seqinfo(gr) <- seqinf[used_chrs]
  gr <- GenomeInfoDb::sortSeqlevels(gr)
  gr <- BiocGenerics::sort(gr)

  if (anyNA(GenomeInfoDb::seqlengths(gr))) {
    stop("Some seqlengths in final 'gr' are NA.", call. = FALSE)
  }
  rtracklayer::export(gr, file_path, format = "BigWig")
  invisible(file_path)
}
