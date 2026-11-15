# G1-derived variable-width windows and multi-sample bin synchronisation.
# define_G1_windows() is migrated from rcode/processing/define_G1_windows.R;
# synchronize_bins() is extracted from notebook 1.3.synchronize_gr_filtered_*.

#' Define G1-derived windows from a fragment file
#'
#' Build variable-width genomic windows that each hold `target_reads` G1
#' fragments, and count fragments of all cells in those windows.
#'
#' Steps:
#' \enumerate{
#'   \item `genome_granges` is merged with `GenomicRanges::reduce(min.gapwidth)`
#'     into contigs; contigs shorter than `target_window_length` and contigs on
#'     chromosomes absent from the tabix index are dropped.
#'   \item G1 cells (`meta_data[[g1_col]]` is `TRUE`, finite coverage) are
#'     ranked by decreasing coverage.
#'   \item On a pilot chromosome (largest total contig width), a binary search
#'     chooses the number `n` of top-ranked G1 cells whose windows have a
#'     median length closest to `target_window_length`.
#'   \item For every chromosome, fragment midpoints
#'     (`((start0 + 1) + end) %/% 2`) of the `n` selected G1 cells are sorted
#'     per contig and cut into consecutive blocks of `target_reads`; each
#'     window runs from the first fragment of its block to just before the
#'     first fragment of the next block (the last block ends at its last
#'     fragment; the incomplete remainder is dropped). If a chromosome has
#'     `>= 100` windows, windows outside the 1st-99th length percentiles are
#'     removed.
#'   \item Fragments of all cells in `rownames(meta_data)` are counted per
#'     window (by midpoint).
#' }
#'
#' @param genome_granges A `GRanges` of retained genome bins (already filtered
#'   for gaps / CNV / mappability, e.g. the output of [filter_windows()],
#'   unlisted, or of [synchronize_bins()]).
#' @param frag_path Path to a bgzipped, tabix-indexed 10x fragment file
#'   (`<frag_path>.tbi` must exist). Columns: chr, start (0-based), end,
#'   barcode, ...
#' @param meta_data A `data.frame` whose unique row names are cell barcodes.
#' @param g1_col Column of `meta_data` flagging G1 cells (coerced to logical).
#' @param cov_col Numeric column of `meta_data` with per-cell coverage.
#' @param target_reads Positive integer; G1 fragments per window. We use
#'   `200`.
#' @param target_window_length Positive integer; target median window length
#'   (bp) and minimum contig width, e.g. `20000`.
#' @param min.gapwidth Positive integer passed to `GenomicRanges::reduce()`.
#'   Default `100000L`.
#' @param verbose Logical; report progress. Default `FALSE`.
#'
#' @return A list with
#'   \describe{
#'     \item{`g1bins_counts`}{`dgCMatrix` windows x cells (column names =
#'       `rownames(meta_data)`, row names = window IDs `"<chr>-<start>-<end>"`).}
#'     \item{`g1bins`}{named `GRanges` of windows with metadata columns
#'       `midpoint`, `g1_read_count`, `is_masked_window` (G1 count differs
#'       from `target_reads` by more than 1) and `contig_id`. Its
#'       `S4Vectors::metadata()` holds `target_reads`, `target_window_length`,
#'       `pilot_chr`, `selected_g1_n`, `selected_g1_barcodes`,
#'       `selected_cov_cutoff`, `pilot_median_length`, `pilot_n_windows`.}
#'   }
#'
#' @export
define_G1_windows <- function(genome_granges,
                              frag_path,
                              meta_data,
                              g1_col = "is_G1",
                              cov_col = "coverage",
                              target_reads = 200,
                              target_window_length = 20000,
                              min.gapwidth = 100000L,
                              verbose = FALSE) {
  vmsg <- function(...) .vmsg(verbose, ...)

  # =========================================================================
  # 1-2. Input validation
  # =========================================================================
  vmsg("[2/10] Validating inputs...")

  check_pos_int <- function(x, name) {
    if (!is.numeric(x) ||
        length(x) != 1L ||
        is.na(x) ||
        !is.finite(x) ||
        x <= 0 ||
        x != as.integer(x)) {
      stop("'", name, "' must be one positive integer.")
    }
    as.integer(x)
  }

  .validate_bins(genome_granges, "genome_granges")
  if (!is.data.frame(meta_data)) {
    stop("'meta_data' must be a data.frame.")
  }
  if (!g1_col %in% colnames(meta_data)) {
    stop("Missing column '", g1_col, "' in 'meta_data'.")
  }
  if (!cov_col %in% colnames(meta_data)) {
    stop("Missing column '", cov_col, "' in 'meta_data'.")
  }
  if (!is.numeric(meta_data[[cov_col]])) {
    stop("'", cov_col, "' must be numeric.")
  }
  if (!is.character(frag_path) || length(frag_path) != 1L || is.na(frag_path)) {
    stop("'frag_path' must be a single non-missing character string.")
  }
  if (!file.exists(frag_path)) {
    stop("Fragment file does not exist: ", frag_path)
  }
  if (!file.exists(paste0(frag_path, ".tbi"))) {
    stop("Tabix index not found: ", frag_path, ".tbi")
  }

  target_reads <- check_pos_int(target_reads, "target_reads")
  target_window_length <- check_pos_int(target_window_length, "target_window_length")
  min.gapwidth <- check_pos_int(min.gapwidth, "min.gapwidth")

  cell_ids <- rownames(meta_data)
  if (is.null(cell_ids) || anyNA(cell_ids) || anyDuplicated(cell_ids)) {
    stop("Cell names in rownames(meta_data) must be non-missing and unique.")
  }

  # =========================================================================
  # 4. Prepare G1 cells ranked by coverage
  # =========================================================================
  vmsg("[4/10] Selecting and ranking G1 cells...")
  g1_flag <- as.logical(meta_data[[g1_col]])
  g1_idx <- !is.na(g1_flag) & g1_flag
  if (!any(g1_idx)) {
    stop("No G1 cells found.")
  }

  g1_cov <- meta_data[g1_idx, cov_col, drop = TRUE]
  g1_ids <- rownames(meta_data)[g1_idx]

  keep_cov <- !is.na(g1_cov) & is.finite(g1_cov)
  g1_cov <- as.numeric(g1_cov[keep_cov])
  g1_ids <- g1_ids[keep_cov]

  if (length(g1_ids) == 0L) {
    stop("No G1 cells with finite coverage.")
  }

  ord_cov <- order(-g1_cov, g1_ids)
  g1_ranked <- g1_ids[ord_cov]
  g1_cov_ranked <- g1_cov[ord_cov]

  # =========================================================================
  # 6. Mappable regions
  # =========================================================================
  vmsg("[6/10] Preparing mappable regions...")

  mappable_gr <- GenomicRanges::reduce(genome_granges, min.gapwidth = min.gapwidth)
  mappable_gr <- mappable_gr[IRanges::width(mappable_gr) >= target_window_length]

  if (length(mappable_gr) == 0L) {
    stop("No mappable regions remain.")
  }

  tabix_seqnames <- Rsamtools::seqnamesTabix(frag_path)
  keep_seq <- as.character(GenomicRanges::seqnames(mappable_gr)) %in% tabix_seqnames
  mappable_gr <- mappable_gr[keep_seq]

  if (length(mappable_gr) == 0L) {
    stop(
      "No overlapping chromosome names between genome_granges and fragment file. ",
      "genome_granges examples: ",
      paste(utils::head(unique(as.character(GenomicRanges::seqnames(genome_granges)))),
            collapse = ", "),
      "; fragment file examples: ",
      paste(utils::head(tabix_seqnames), collapse = ", ")
    )
  }

  map_dt <- data.table::data.table(
    chr = as.character(GenomicRanges::seqnames(mappable_gr)),
    chr_norm = .g1w_strip_chr(GenomicRanges::seqnames(mappable_gr)),
    start = BiocGenerics::start(mappable_gr),
    end = BiocGenerics::end(mappable_gr)
  )

  chr_order <- unique(map_dt$chr_norm)
  ord <- order(match(map_dt$chr_norm, chr_order), map_dt$start, map_dt$end)

  map_dt <- map_dt[ord]
  mappable_gr <- mappable_gr[ord]

  # =========================================================================
  # 7. Pilot chromosome selection
  # =========================================================================
  vmsg("[7/10] Selecting G1 cells using a pilot chromosome...")

  chr_width_dt <- map_dt[, list(total_width = sum(end - start + 1L)), by = chr_norm]
  data.table::setorder(chr_width_dt, -total_width, chr_norm)

  pilot_chr <- chr_width_dt$chr_norm[1]
  vmsg("  - Pilot chromosome: ", pilot_chr)

  pilot_fetched <- .g1w_fetch_chr_valid(
    chr_i = pilot_chr,
    map_dt = map_dt,
    mappable_gr = mappable_gr,
    frag_path = frag_path,
    cell_ids = cell_ids,
    verbose = verbose
  )

  if (is.null(pilot_fetched) || nrow(pilot_fetched$chr_valid) == 0L) {
    stop("Pilot chromosome '", pilot_chr, "' has no valid fragments.")
  }

  eval_cache <- new.env(parent = emptyenv())

  eval_n_cells <- function(n_use) {
    key <- as.character(n_use)

    if (exists(key, envir = eval_cache, inherits = FALSE)) {
      return(get(key, envir = eval_cache, inherits = FALSE))
    }

    g1_use <- g1_ranked[seq_len(n_use)]

    windows_dt <- .g1w_build_windows_dt(
      chr_valid = pilot_fetched$chr_valid,
      chr_contigs = pilot_fetched$chr_contigs,
      g1_barcodes = g1_use,
      target_reads = target_reads
    )

    if (is.null(windows_dt) || nrow(windows_dt) == 0L) {
      stat <- list(
        n_cells = n_use,
        median_length = NA_real_,
        n_windows = 0L
      )
    } else {
      lens <- windows_dt$end - windows_dt$start + 1L
      stat <- list(
        n_cells = n_use,
        median_length = stats::median(lens),
        n_windows = nrow(windows_dt)
      )
    }

    assign(key, stat, envir = eval_cache)

    if (isTRUE(verbose)) {
      med_str <- if (is.finite(stat$median_length)) {
        format(round(stat$median_length, 2), trim = TRUE)
      } else {
        "NA"
      }
      message(
        "    [pilot eval] n_cells=", stat$n_cells,
        ", median_length=", med_str,
        ", n_windows=", stat$n_windows
      )
    }

    stat
  }

  n_g1 <- length(g1_ranked)

  if (n_g1 == 1L) {
    best_stat <- eval_n_cells(1L)
  } else {
    st_lo <- eval_n_cells(1L)
    st_hi <- eval_n_cells(n_g1)

    if (!is.finite(st_hi$median_length)) {
      stop("Not enough G1 fragments even when using all ranked G1 cells.")
    }

    best_stat <- st_hi
    best_diff <- abs(st_hi$median_length - target_window_length)

    if (is.finite(st_lo$median_length) &&
        abs(st_lo$median_length - target_window_length) < best_diff) {
      best_stat <- st_lo
      best_diff <- abs(st_lo$median_length - target_window_length)
    }

    low <- 1L
    high <- n_g1
    max_iter <- ceiling(log2(n_g1)) + 6L
    iter <- 0L

    while (low <= high && iter < max_iter) {
      iter <- iter + 1L
      mid <- low + (high - low) %/% 2L

      st_mid <- eval_n_cells(mid)

      if (is.finite(st_mid$median_length)) {
        d_mid <- abs(st_mid$median_length - target_window_length)
        if (d_mid < best_diff) {
          best_stat <- st_mid
          best_diff <- d_mid
        }
      }

      if (!is.finite(st_mid$median_length)) {
        low <- mid + 1L
      } else if (st_mid$median_length > target_window_length) {
        low <- mid + 1L
      } else if (st_mid$median_length < target_window_length) {
        high <- mid - 1L
      } else {
        best_stat <- st_mid
        break
      }
    }
  }

  selected_n <- best_stat$n_cells
  selected_g1_barcodes <- g1_ranked[seq_len(selected_n)]
  selected_cov_cutoff <- g1_cov_ranked[selected_n]
  barcode_index <- stats::setNames(seq_along(cell_ids), cell_ids)
  seqinfo_ref <- GenomeInfoDb::seqinfo(mappable_gr)

  vmsg("[8/10] Releasing temporary pilot objects...")
  rm(pilot_fetched, eval_cache)

  # =========================================================================
  # 8. Final per-chromosome processing
  # =========================================================================
  vmsg("[9/10] Running genome-wide window construction and counting...")

  process_one_chr <- function(chr_i) {
    fetched <- .g1w_fetch_chr_valid(
      chr_i = chr_i,
      map_dt = map_dt,
      mappable_gr = mappable_gr,
      frag_path = frag_path,
      cell_ids = cell_ids,
      verbose = verbose
    )

    if (is.null(fetched)) return(NULL)

    windows_chr <- .g1w_build_windows_dt(
      chr_valid = fetched$chr_valid,
      chr_contigs = fetched$chr_contigs,
      g1_barcodes = selected_g1_barcodes,
      target_reads = target_reads
    )

    if (is.null(windows_chr) || nrow(windows_chr) == 0L) {
      return(NULL)
    }

    n_win_chr <- nrow(windows_chr)

    idx <- findInterval(fetched$chr_valid$pos, windows_chr$start)
    valid <- idx > 0L & idx <= n_win_chr
    keep <- rep(FALSE, length(idx))
    keep[valid] <- fetched$chr_valid$pos[valid] <= windows_chr$end[idx[valid]]

    if (any(keep)) {
      local_i <- idx[keep]
      local_j <- unname(barcode_index[fetched$chr_valid$barcode[keep]])

      chr_counts <- Matrix::sparseMatrix(
        i = local_i,
        j = local_j,
        x = 1L,
        dims = c(n_win_chr, length(cell_ids))
      )
    } else {
      chr_counts <- Matrix::sparseMatrix(
        i = integer(0),
        j = integer(0),
        dims = c(n_win_chr, length(cell_ids))
      )
    }

    chr_gr <- GenomicRanges::GRanges(
      seqnames = windows_chr$seqnames,
      ranges = IRanges::IRanges(
        start = windows_chr$start,
        end = windows_chr$end
      ),
      seqinfo = seqinfo_ref
    )

    S4Vectors::mcols(chr_gr) <- S4Vectors::DataFrame(
      midpoint = windows_chr$midpoint,
      g1_read_count = windows_chr$g1_read_count,
      is_masked_window = windows_chr$is_masked_window,
      contig_id = windows_chr$contig_id
    )

    list(gr = chr_gr, counts = chr_counts)
  }

  res_list <- lapply(chr_order, process_one_chr)
  res_list <- res_list[!vapply(res_list, is.null, logical(1))]

  if (length(res_list) == 0L) {
    vmsg("  - No windows were produced genome-wide; creating empty result.")
    g1bins <- GenomicRanges::GRanges()
    g1bins_counts <- Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      dims = c(0L, length(cell_ids))
    )
  } else {
    gr_list <- lapply(res_list, `[[`, "gr")
    count_list <- lapply(res_list, `[[`, "counts")

    g1bins <- do.call(c, gr_list)
    g1bins_counts <- do.call(rbind, count_list)

    rm(gr_list, count_list, res_list)
  }

  # =========================================================================
  # 9. Finalize and return
  # =========================================================================
  vmsg("[10/10] Finalizing output...")

  if (length(g1bins) > 0L) {
    window_ids <- paste0(
      as.character(GenomicRanges::seqnames(g1bins)), "-",
      BiocGenerics::start(g1bins), "-",
      BiocGenerics::end(g1bins)
    )
    names(g1bins) <- window_ids
    rownames(g1bins_counts) <- window_ids
  } else {
    rownames(g1bins_counts) <- character(0)
  }

  colnames(g1bins_counts) <- cell_ids

  md <- S4Vectors::metadata(g1bins)
  md$target_reads <- target_reads
  md$target_window_length <- target_window_length
  md$pilot_chr <- pilot_chr
  md$selected_g1_n <- selected_n
  md$selected_g1_barcodes <- selected_g1_barcodes
  md$selected_cov_cutoff <- selected_cov_cutoff
  md$pilot_median_length <- best_stat$median_length
  md$pilot_n_windows <- best_stat$n_windows
  S4Vectors::metadata(g1bins) <- md

  list(
    g1bins_counts = g1bins_counts,
    g1bins = g1bins
  )
}

# ---------------------------------------------------------------------------
# Internal helpers of define_G1_windows() (were closures in the original).
# ---------------------------------------------------------------------------

.g1w_strip_chr <- function(x) {
  sub("^chr", "", as.character(x), ignore.case = TRUE)
}

.g1w_parse_tabix_lines <- function(lines_list) {
  n_per_contig <- lengths(lines_list)
  keep_contig <- which(n_per_contig > 0L)
  if (length(keep_contig) == 0L) return(NULL)
  out <- vector("list", length(keep_contig))
  kk <- 0L

  for (cid in keep_contig) {
    lines <- lines_list[[cid]]
    keep <- !is.na(lines) & nzchar(lines) & substr(lines, 1L, 1L) != "#"
    lines <- lines[keep]
    if (!length(lines)) next
    parts <- data.table::tstrsplit(lines, "\t", fixed = TRUE, keep = 1:4)
    dt <- data.table::data.table(
      chr = parts[[1]],
      start0 = suppressWarnings(as.integer(parts[[2]])),
      end = suppressWarnings(as.integer(parts[[3]])),
      barcode = parts[[4]],
      contig_id = cid)
    dt <- dt[!is.na(start0) & !is.na(end) & !is.na(barcode)]
    if (nrow(dt) == 0L) next
    dt[, pos := as.integer(((start0 + 1L) + end) %/% 2L)]
    kk <- kk + 1L
    out[[kk]] <- dt[, list(pos, barcode, contig_id)]
    rm(lines, parts, dt)
  }
  if (kk == 0L) return(NULL)
  data.table::rbindlist(out[seq_len(kk)], use.names = TRUE)
}

.g1w_build_windows_dt <- function(chr_valid, chr_contigs, g1_barcodes, target_reads) {
  if (is.null(chr_valid) || nrow(chr_valid) == 0L) return(NULL)
  g1_chr <- chr_valid[barcode %in% g1_barcodes]
  if (nrow(g1_chr) == 0L) return(NULL)
  data.table::setorder(g1_chr, contig_id, pos)
  out_list <- vector("list", length(chr_contigs$contig_id))
  kk <- 0L
  for (cid in chr_contigs$contig_id) {
    g1_pos <- sort(g1_chr[contig_id == cid, pos])
    n_complete <- length(g1_pos) %/% target_reads
    if (n_complete < 1L) next
    start_idx <- seq.int(1L, by = target_reads, length.out = n_complete)
    end_idx <- start_idx + target_reads - 1L
    win_start <- g1_pos[start_idx]
    if (n_complete == 1L) {
      win_end <- g1_pos[end_idx]
    } else {
      win_end <- c(
        g1_pos[start_idx[-1L]] - 1L,
        g1_pos[end_idx[n_complete]]
      )
    }
    # NOTE (kept from original): with duplicated midpoints at a block
    # boundary, pmax() can make a window end at the next window's start,
    # producing overlapping windows. Reported as a suspected edge case.
    win_end <- pmax(win_end, win_start)
    actual_g1_reads <-
      findInterval(win_end, g1_pos, rightmost.closed = TRUE) -
      findInterval(win_start - 1L, g1_pos, rightmost.closed = TRUE)

    is_masked <- actual_g1_reads > (target_reads + 1L) |
      actual_g1_reads < (target_reads - 1L)

    kk <- kk + 1L
    out_list[[kk]] <- data.table::data.table(
      seqnames = chr_contigs[contig_id == cid, chr][1],
      start = win_start,
      end = win_end,
      midpoint = floor((win_start + win_end) / 2),
      g1_read_count = actual_g1_reads,
      is_masked_window = is_masked,
      contig_id = cid
    )
  }

  if (kk == 0L) return(NULL)
  windows_chr <- data.table::rbindlist(out_list[seq_len(kk)])
  data.table::setorder(windows_chr, start, end)
  window_length <- windows_chr$end - windows_chr$start + 1L
  if (length(window_length) >= 100L) {
    q <- stats::quantile(
      window_length,
      probs = c(0.01, 0.99),
      na.rm = TRUE,
      names = FALSE
    )
    keep_len <- window_length >= q[1] & window_length <= q[2]
    windows_chr <- windows_chr[keep_len]
  }
  if (nrow(windows_chr) == 0L) return(NULL)
  windows_chr
}

.g1w_fetch_chr_valid <- function(chr_i, map_dt, mappable_gr, frag_path, cell_ids,
                                 verbose = FALSE) {
  idx_chr <- which(map_dt$chr_norm == chr_i)
  if (length(idx_chr) == 0L) return(NULL)

  chr_contigs_gr <- mappable_gr[idx_chr]
  chr_contigs <- data.table::copy(map_dt[idx_chr])
  chr_contigs[, contig_id := seq_len(.N)]

  .vmsg(verbose, "    [fetch fragments] scanTabix start for chr", chr_i)
  lines_list <- Rsamtools::scanTabix(file = frag_path, param = chr_contigs_gr)
  .vmsg(verbose, "    [fetch fragments] scanTabix done: ", sum(lengths(lines_list)), " lines")

  .vmsg(verbose, "    [fetch fragments] parse start")
  chr_valid <- .g1w_parse_tabix_lines(lines_list)
  .vmsg(verbose, "    [fetch fragments] parse done: ",
        if (is.null(chr_valid)) 0L else nrow(chr_valid), " rows")

  if (is.null(chr_valid) || nrow(chr_valid) == 0L) return(NULL)

  chr_valid <- chr_valid[barcode %in% cell_ids]
  if (nrow(chr_valid) == 0L) return(NULL)

  data.table::setorder(chr_valid, contig_id, pos)

  list(chr_valid = chr_valid, chr_contigs = chr_contigs)
}

#' Synchronise retained bins across samples
#'
#' Apply a shared bin filter to several samples: for each sample the *removed* bins are the
#' bins of `gr_cn_list[[s]]` with no exactly matching bin (same seqnames,
#' start, end, strand) in `gr_filtered_list[[s]]`; the removed bins of the
#' samples in `union_samples` are pooled (unique coordinates), and every
#' sample's retained bins that exactly match a pooled removed bin are dropped.
#' A bin is therefore kept only if it survived filtering in all
#' `union_samples` (among the bins present in each sample's unfiltered set).
#'
#' @param gr_cn_list Named list (or `GRangesList`) of the unfiltered per-sample
#'   bins.
#' @param gr_filtered_list Named list (or `GRangesList`) of the per-sample
#'   retained bins (`gr_filtered.rds`, e.g. unlisted [filter_windows()]
#'   output), with the same names as `gr_cn_list`.
#' @param union_samples Names of the samples whose removed bins are pooled.
#'   Default: all samples.
#'
#' @return A list with
#'   \describe{
#'     \item{`filtered`}{named list of `GRanges`: `gr_filtered_list` with the
#'       pooled removed bins dropped; metadata columns are kept.}
#'     \item{`removed`}{named list of `GRanges`: per-sample removed bins.}
#'     \item{`removed_union`}{`GRanges` of pooled removed bins (no metadata
#'       columns).}
#'   }
#'
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(seq(1, by = 10, length.out = 6),
#'                                                       width = 10))
#' res <- synchronize_bins(
#'   list(a = gr, b = gr),
#'   list(a = gr[-2], b = gr[-5])
#' )
#' res$filtered$a   # bins 2 and 5 removed in both samples
#'
#' @export
synchronize_bins <- function(gr_cn_list, gr_filtered_list,
                             union_samples = names(gr_cn_list)) {
  if (methods::is(gr_cn_list, "GRangesList")) gr_cn_list <- as.list(gr_cn_list)
  if (methods::is(gr_filtered_list, "GRangesList")) {
    gr_filtered_list <- as.list(gr_filtered_list)
  }
  if (!is.list(gr_cn_list) || !is.list(gr_filtered_list)) {
    stop("'gr_cn_list' and 'gr_filtered_list' must be lists of GRanges.")
  }
  samples <- names(gr_cn_list)
  if (is.null(samples) || anyNA(samples) || any(samples == "") ||
      anyDuplicated(samples)) {
    stop("'gr_cn_list' must have unique, non-empty names.")
  }
  if (!setequal(samples, names(gr_filtered_list))) {
    stop("'gr_cn_list' and 'gr_filtered_list' must have the same names.")
  }
  for (s in samples) {
    .validate_bins(gr_cn_list[[s]], paste0("gr_cn_list[['", s, "']]"),
                   allow_empty = TRUE)
    .validate_bins(gr_filtered_list[[s]], paste0("gr_filtered_list[['", s, "']]"),
                   allow_empty = TRUE)
  }
  if (!is.character(union_samples) || length(union_samples) == 0L ||
      !all(union_samples %in% samples)) {
    stop("'union_samples' must be a non-empty subset of names(gr_cn_list).")
  }

  removed <- lapply(stats::setNames(samples, samples), function(s) {
    gr_cn <- gr_cn_list[[s]]
    hits <- GenomicRanges::findOverlaps(gr_cn, gr_filtered_list[[s]], type = "equal")
    kept_idx <- S4Vectors::queryHits(hits)
    # NOTE (kept from notebook): gr_cn[-integer(0)] is empty, so a sample with
    # no retained bins contributes no removed bins. Reported as an edge case.
    gr_cn[-kept_idx]
  })

  pooled <- do.call(c, unname(lapply(removed[union_samples], function(gr) {
    S4Vectors::mcols(gr) <- NULL
    gr
  })))
  removed_union <- unique(pooled)

  filtered <- lapply(stats::setNames(samples, samples), function(s) {
    gr <- gr_filtered_list[[s]]
    hits <- GenomicRanges::findOverlaps(gr, removed_union, type = "equal")
    remove_indices <- S4Vectors::queryHits(hits)
    gr[!(seq_along(gr) %in% remove_indices)]
  })

  list(filtered = filtered, removed = removed, removed_union = removed_union)
}
