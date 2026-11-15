# Replication-state QC and percent replication.
# compute_rep_percent() is migrated from rcode/processing/compute_rep_percent.R
# (notebooks 1.4 / 1.5).

#' Compute percent replication and QC cells from bin-level replication states
#'
#' Run two QC passes on a bins x cells replication-state matrix (states `2` =
#' unreplicated, `4` = replicated, `NA` = unknown), compute the per-cell
#' replicated fraction and return the filtered matrix with updated metadata.
#'
#' \describe{
#'   \item{First pass (cell level)}{`rep_percent` = fraction of non-`NA` bins
#'     in state 4; `sigma_chr` = population SD of the per-chromosome mean
#'     states. A cell fails if `rep_percent < min_rep_percent` or is `NA`,
#'     `sigma_chr > max_sigma`, or at least one chromosome is entirely 2 while
#'     another is entirely 4 (`has_cnv`). Passing cells are ordered by
#'     increasing `rep_percent`.}
#'   \item{Second pass (per chromosome)}{Pairwise Pearson correlations of the
#'     cells' states on the chromosome are computed (values equal to 1 are set
#'     to `NA`); cells with a negative mean correlation are masked. For each
#'     cell with both states present on the chromosome, among its
#'     `neighbor_window` neighbours on each side in `rep_percent` order, the
#'     cell is masked if no neighbour correlates better with any column than
#'     it does, or if the fraction of MAD outliers (|r - median| > 3 MAD, IQR
#'     rule if MAD = 0) in its row exceeds `outlier_threshold`. Masked cells get
#'     `NA` on that chromosome; a cell masked on at least
#'     `n_chromosomes - 1` chromosomes fails the second pass.}
#' }
#'
#' @param genobins_rep_states Matrix-like bins x cells replication states;
#'   row names = `names(genome_granges)`, column names = `rownames(meta_data)`
#'   (same order).
#' @param genome_granges Named `GRanges` with one range per row.
#' @param meta_data A `data.frame` whose unique row names are cell IDs.
#' @param min_rep_percent Minimum replicated fraction (first pass). We
#'   use `0.01`.
#' @param max_sigma Maximum `sigma_chr` (first pass). We use `0.3`.
#' @param neighbor_window Neighbourhood half-width for the second pass.
#'   We use `5`.
#' @param outlier_threshold Maximum fraction of local correlation outliers
#'   (second pass). We use `0.2`.
#' @param BPPARAM A `BiocParallel::BiocParallelParam`; chromosomes are
#'   processed in parallel in the second pass.
#' @param verbose Logical; print progress messages. Default `FALSE`.
#'
#' @return A list with
#'   \describe{
#'     \item{`genobins_rep_states_filtered`}{numeric matrix, all bins x cells
#'       passing both QC passes (columns ordered by increasing `rep_percent`),
#'       with second-pass masked chromosome blocks set to `NA`.}
#'     \item{`meta_data`}{input `meta_data` with added columns `rep_percent`,
#'       `sigma_chr` (numeric, all cells) and `pass_qc` (logical).}
#'   }
#'
#' @examples
#' set.seed(1)
#' gr <- GenomicRanges::GRanges(rep(c("chr1", "chr2", "chr3"), each = 100),
#'                              IRanges::IRanges(rep(seq(1, by = 1000, length.out = 100), 3),
#'                                               width = 1000))
#' names(gr) <- paste0("bin", seq_along(gr))
#' rp <- seq(0.1, 0.9, length.out = 20)
#' st <- sapply(rp, function(p) ifelse(rep(seq(0, 1, length.out = 100), 3) < p, 4, 2))
#' dimnames(st) <- list(names(gr), paste0("c", 1:20))
#' md <- data.frame(row.names = colnames(st), cell_id = colnames(st))
#' res <- compute_rep_percent(st, gr, md, min_rep_percent = 0.01, max_sigma = 0.3,
#'                            verbose = FALSE)
#' head(res$meta_data)
#'
#' @export
compute_rep_percent <- function(genobins_rep_states,
                                genome_granges,
                                meta_data,
                                min_rep_percent = 0.001,
                                max_sigma = 0.4,
                                neighbor_window = 5,
                                outlier_threshold = 0.2,
                                BPPARAM = BiocParallel::SerialParam(),
                                verbose = FALSE) {
  vmsg <- function(...) .vmsg(verbose, ...)

  # =========================================================================
  # 2. Input validation
  # =========================================================================
  vmsg("[2/7] Validating inputs...")

  .validate_bins(genome_granges, "genome_granges")
  genome_granges <- GenomeInfoDb::keepSeqlevels(
    genome_granges,
    value = GenomeInfoDb::seqlevelsInUse(genome_granges),
    pruning.mode = "coarse"
  )
  if (!is.data.frame(meta_data)) {
    stop("'meta_data' must be a data.frame.")
  }

  meta_cells <- rownames(meta_data)
  if (is.null(meta_cells) || anyNA(meta_cells) || anyDuplicated(meta_cells)) {
    stop("rownames(meta_data) must be non-missing and unique.")
  }

  if (!is.numeric(min_rep_percent) || length(min_rep_percent) != 1L ||
      is.na(min_rep_percent) || min_rep_percent < 0 || min_rep_percent > 1) {
    stop("'min_rep_percent' must be a single numeric value between 0 and 1.")
  }
  if (!is.numeric(max_sigma) || length(max_sigma) != 1L ||
      is.na(max_sigma) || max_sigma < 0) {
    stop("'max_sigma' must be a single non-negative numeric value.")
  }
  if (!is.numeric(neighbor_window) || length(neighbor_window) != 1L ||
      is.na(neighbor_window) || neighbor_window < 0) {
    stop("'neighbor_window' must be a single non-negative numeric value.")
  }
  neighbor_window <- as.integer(neighbor_window)
  if (!is.numeric(outlier_threshold) || length(outlier_threshold) != 1L ||
      is.na(outlier_threshold) || outlier_threshold < 0 || outlier_threshold > 1) {
    stop("'outlier_threshold' must be a single numeric value between 0 and 1.")
  }
  if (!inherits(BPPARAM, "BiocParallelParam")) {
    stop("'BPPARAM' must be a BiocParallelParam object.")
  }

  # =========================================================================
  # 4. Prepare replication-state matrix
  # =========================================================================
  vmsg("[4/7] Preparing replication-state matrix...")
  state <- .crp_as_state_matrix(genobins_rep_states, genome_granges)
  barcode_df <- .crp_make_barcodes(state, meta_bc = meta_cells)

  chr <- as.character(GenomicRanges::seqnames(genome_granges))
  chr <- factor(chr, levels = GenomeInfoDb::seqlevels(genome_granges))
  idx <- split(seq_len(length(genome_granges)), chr)
  state_chr <- lapply(idx, function(i) state[i, , drop = FALSE])

  vmsg("  - ", nrow(state), " bins x ", ncol(state), " cells across ",
       length(idx), " chromosomes.")

  # =========================================================================
  # 5. First QC
  # =========================================================================
  vmsg("[5/7] Running first cell-level QC...")
  rep_percent <- .crp_calc_percent(state_chr)

  mu_chr <- do.call(rbind, lapply(state_chr, .crp_col_means_na))
  sigma_chr <- apply(mu_chr, 2, .crp_sd_pop_na)

  has_cnv <- apply(mu_chr == 2, 2, any, na.rm = TRUE) &
    apply(mu_chr == 4, 2, any, na.rm = TRUE)

  bad <- rep_percent < min_rep_percent | is.na(rep_percent) |
    sigma_chr > max_sigma | has_cnv
  pass_1st_qc <- !bad
  names(pass_1st_qc) <- names(rep_percent)

  keep <- which(pass_1st_qc)
  keep <- keep[order(rep_percent[keep])]

  barcode_df <- barcode_df[keep, , drop = FALSE]
  state_chr <- lapply(state_chr, function(m) m[, keep, drop = FALSE])

  vmsg("  - Kept ", length(keep), "/", length(rep_percent), " cells after QC.")

  # =========================================================================
  # 6. Second-pass QC
  # =========================================================================
  vmsg("[6/7] Running second chromosome-wise QC...")
  nc <- length(keep)
  chr_names <- names(state_chr)

  if (nc >= 3L && length(chr_names) > 0L) {
    .qc_one_chr <- function(chr) {
      matrix_chr <- state_chr[[chr]]
      r <- suppressWarnings(stats::cor(matrix_chr, use = "pairwise.complete.obs"))
      if (is.null(dim(r))) {
        r <- matrix(r, nrow = nc, ncol = nc)
      }
      # NOTE (kept from original): this removes the diagonal but also any
      # off-diagonal correlation of exactly 1.
      r[r == 1] <- NA_real_

      mean_r <- colMeans(r, na.rm = TRUE)
      mean_r[is.nan(mean_r)] <- NA_real_
      neg <- !is.na(mean_r) & mean_r < 0

      if (any(neg)) {
        r[neg, ] <- NA_real_
        r[, neg] <- NA_real_
      }

      bad_local <- rep(FALSE, nc)
      for (i in seq_len(nc)) {
        if (sum(matrix_chr[, i] == 4, na.rm = TRUE) == 0 ||
            sum(matrix_chr[, i] == 2, na.rm = TRUE) == 0) {
          next
        }

        nb <- seq.int(max(1L, i - neighbor_window), min(nc, i + neighbor_window))
        r1 <- r[nb, nb, drop = FALSE]
        pos <- which(nb == i)

        comp <- sweep(r1, 2, r1[pos, ], `>`)
        comp[is.na(comp)] <- FALSE

        out <- matrix(
          .crp_is_outlier_mad(as.vector(r1)),
          nrow = nrow(r1),
          ncol = ncol(r1)
        )

        denom_local <- nrow(r1) - 1L
        frac_outlier <- if (denom_local <= 0L) {
          0
        } else {
          sum(out[pos, ], na.rm = TRUE) / denom_local
        }

        bad_local[i] <- (sum(comp) == 0) || (frac_outlier > outlier_threshold)
      }

      keep_vec <- !(neg | bad_local)
      if (any(!keep_vec)) {
        matrix_chr[, !keep_vec] <- NA_real_
      }

      list(chr = chr, keep = keep_vec, state = matrix_chr)
    }

    qc_res <- BiocParallel::bplapply(
      X = chr_names,
      FUN = .qc_one_chr,
      BPPARAM = BPPARAM
    )
    names(qc_res) <- chr_names

    keep_chr <- do.call(rbind, lapply(qc_res, function(x) x$keep))
    rownames(keep_chr) <- chr_names
    colnames(keep_chr) <- barcode_df$cell_id

    state_chr <- stats::setNames(lapply(qc_res, function(x) x$state), chr_names)

    bad2 <- if (length(chr_names) > 1L) {
      colSums(!keep_chr) >= (length(chr_names) - 1L)
    } else {
      rep(FALSE, nc)
    }

    keep2 <- which(!bad2)
    pass_2nd_qc <- stats::setNames(rep(FALSE, nc), barcode_df$cell_id)
    if (length(keep2) > 0L) {
      pass_2nd_qc[barcode_df$cell_id[keep2]] <- TRUE
    }

    barcode_df <- barcode_df[keep2, , drop = FALSE]
    state_chr <- lapply(state_chr, function(m) m[, keep2, drop = FALSE])

    vmsg("  - Kept ", length(keep2), "/", nc, " cells after second QC.")
  } else {
    pass_2nd_qc <- stats::setNames(rep(TRUE, nc), barcode_df$cell_id)
    vmsg("  - Skipped second QC (fewer than 3 cells retained).")
  }

  # =========================================================================
  # 7. Finalize output
  # =========================================================================
  vmsg("[7/7] Computing final filtered matrix and metadata...")

  row_ids <- names(genome_granges)

  clean_states <- matrix(
    NA_real_,
    nrow = length(genome_granges),
    ncol = nrow(barcode_df),
    dimnames = list(row_ids, barcode_df$cell_id)
  )

  for (chr in names(idx)) {
    clean_states[idx[[chr]], ] <- state_chr[[chr]]
  }

  meta_data$rep_percent <- NA_real_
  common_cells <- intersect(rownames(meta_data), names(rep_percent))
  meta_data[common_cells, "rep_percent"] <- rep_percent[common_cells]

  meta_data$sigma_chr <- NA_real_
  common_cells <- intersect(rownames(meta_data), names(sigma_chr))
  meta_data[common_cells, "sigma_chr"] <- sigma_chr[common_cells]

  meta_data$pass_qc <- FALSE
  common_cells <- intersect(rownames(meta_data), names(pass_2nd_qc))
  meta_data[common_cells, "pass_qc"] <- pass_2nd_qc[common_cells]

  list(
    genobins_rep_states_filtered = clean_states,
    meta_data = meta_data
  )
}

# ---------------------------------------------------------------------------
# Internal helpers of compute_rep_percent() (were closures in the original).
# ---------------------------------------------------------------------------

.crp_as_state_matrix <- function(x, win) {
  n_win <- length(win)
  if (is.null(x)) {
    stop("'genobins_rep_states' is missing.")
  }
  if (inherits(x, "DataFrame") || is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x)) {
    stop("'genobins_rep_states' must be matrix-like.")
  }
  if (nrow(x) != n_win) {
    stop("nrow(genobins_rep_states) must equal length(genome_granges).")
  }
  if (is.null(rownames(x)) || is.null(names(win))) {
    stop("rownames(genobins_rep_states) or rownames(genome_granges) is NULL.")
  }
  if (!all(rownames(x) == names(win))) {
    stop("rownames(genobins_rep_states) must be identical to rownames(genome_granges).")
  }
  storage.mode(x) <- "numeric"
  x
}

.crp_col_means_na <- function(mat) {
  if (ncol(mat) == 0L) return(numeric(0))
  out <- colMeans(mat, na.rm = TRUE)
  out[colSums(!is.na(mat)) == 0L] <- NA_real_
  out
}

.crp_sd_pop_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_real_)
  mu <- mean(x)
  sqrt(mean((x - mu)^2))
}

.crp_is_outlier_mad <- function(x) {
  out <- rep(FALSE, length(x))
  ok <- is.finite(x)
  if (sum(ok) < 3L) return(out)

  med <- stats::median(x[ok])
  mad0 <- stats::mad(x[ok], center = med, constant = 1)

  if (is.na(mad0) || mad0 == 0) {
    q <- stats::quantile(x[ok], c(0.25, 0.75), na.rm = TRUE, names = FALSE)
    iqr <- q[2] - q[1]
    if (is.na(iqr) || iqr == 0) return(out)
    out[ok] <- x[ok] < (q[1] - 1.5 * iqr) | x[ok] > (q[2] + 1.5 * iqr)
    return(out)
  }

  out[ok] <- abs(x[ok] - med) > 3 * mad0
  out
}

.crp_make_barcodes <- function(state, meta_bc = NULL) {
  nc <- ncol(state)
  state_bc <- colnames(state)

  if (is.null(state_bc) || anyNA(state_bc) || anyDuplicated(state_bc)) {
    stop("colnames(genobins_rep_states) are missing/invalid ")
  }

  if (!identical(state_bc, meta_bc)) {
    stop("colnames(genobins_rep_states) must be identical to rownames(meta_data).")
  }

  data.frame(
    cell_index = seq_len(nc),
    cell_id = state_bc,
    stringsAsFactors = FALSE
  )
}

.crp_calc_percent <- function(lst) {
  if (length(lst) == 0L) return(numeric(0))
  if (ncol(lst[[1]]) == 0L) return(stats::setNames(numeric(0), character(0)))

  denom <- Reduce(`+`, lapply(lst, function(m) colSums(!is.na(m))))
  numer <- Reduce(`+`, lapply(lst, function(m) colSums(m == 4, na.rm = TRUE)))

  out <- numer / denom
  out[denom == 0L] <- NA_real_
  names(out) <- colnames(lst[[1]])
  out
}

# assign_replication_states() is migrated from
# rcode/processing/assign_replication_states.R (notebooks 1.4 / DS 1.4).

#' Assign replication states to G1 and genome windows
#'
#' Run [infer_single_cell_states()] on every cell of a G1-window count matrix
#' (in chunks of `chunk_size` cells, parallelised with `BPPARAM`), using the
#' chromosome starts of `g1_granges` as HMM breakpoints, then project the
#' states onto genome windows.
#'
#' States of masked G1 windows (`mcols(g1_granges)$is_masked_window`, `NA`
#' counts as masked) are set to `NA`. Each genome window takes the state of the
#' unmasked G1 window on the same chromosome with the largest overlap;
#' unmappable genome windows (`mcols(genome_granges)$is_mappable` `FALSE` or
#' `NA`) are set to `NA`, and runs of unmappable windows of at most
#' `gap_max_bp` bp that are not at a chromosome end are filled when both
#' flanking windows have the same state.
#'
#' @param g1bins_matrix A `dgCMatrix` of counts, G1 windows x cells, with row
#'   names identical to `names(g1_granges)` and non-`NULL` column names (e.g.
#'   `define_G1_windows()$g1bins_counts`).
#' @param g1_granges `GRanges` of the G1 windows, grouped by chromosome, with
#'   a logical `is_masked_window` metadata column.
#' @param genome_granges Named `GRanges` of genome windows, grouped by
#'   chromosome, with a logical `is_mappable` metadata column (e.g. from
#'   [compute_bin_mappability()]).
#' @param gap_max_bp Maximum total width (bp) of an unmappable gap to fill.
#' @param BPPARAM A `BiocParallel::BiocParallelParam` used over cell chunks.
#' @param chunk_size Number of cells per chunk.
#' @param infer_args Named list of extra arguments for
#'   [infer_single_cell_states()]. A list given here replaces the default as
#'   a whole; tune the values with [evaluate_state_parameters()].
#' @param verbose Logical; print progress messages. Default `FALSE`.
#'
#' @return A list with
#' \describe{
#'   \item{`genobins_rep_states`}{Integer matrix, genome windows x cells,
#'     values `2`, `4` or `NA`; input for [compute_rep_percent()].}
#'   \item{`g1bins_rep_states`}{Integer matrix, G1 windows x cells.}
#' }
#'
#' @export
assign_replication_states <- function(
    g1bins_matrix,  
    g1_granges,
    genome_granges,
    gap_max_bp = 50000L,
    BPPARAM = BiocParallel::SerialParam(),
    chunk_size = 100L,
    infer_args = list(hmm_aggr_width = 15L, init_ratio_lower = 1.2,
                      final_ratio_lower = 1.4),
    verbose = FALSE
) {
  # ------------------------------------------------------------------------
  # 1. Input validation
  # ------------------------------------------------------------------------
  .vmsg(verbose, "[1/6] Validating input objects...")
  
  raw_counts <- g1bins_matrix
  if (is.null(raw_counts) || !inherits(raw_counts, "dgCMatrix")) {
    stop("g1bins_matrix must be a dgCMatrix.")
  }
  if (nrow(raw_counts) != length(g1_granges)) {
    stop("nrow(g1bins_matrix) must equal length(g1_granges).")
  }
  if (is.null(colnames(raw_counts))) {
    stop("colnames(g1bins_matrix) must not be NULL.")
  }
  if (is.null(rownames(raw_counts))) {
    stop("rownames(g1bins_matrix) must not be NULL.")
  }

  if (!inherits(g1_granges, "GRanges")) {
    stop("g1_ranges must be a GRanges object.")
  }
  if (!"is_masked_window" %in% colnames(S4Vectors::mcols(g1_granges))) {
    stop("mcols(g1_granges)$is_masked_window must exist.")
  }
  is_masked_window <- as.logical(S4Vectors::mcols(g1_granges)$is_masked_window)
  is_masked_window[is.na(is_masked_window)] <- TRUE
  
  g1_granges_names <- names(g1_granges)
  if (is.null(g1_granges_names)) {
    g1_granges_names <- paste0(
      as.character(GenomeInfoDb::seqnames(g1_granges)),
      "-",
      BiocGenerics::start(g1_granges),
      "-",
      BiocGenerics::end(g1_granges)
    )
  }
  if (!identical(as.character(rownames(raw_counts)), as.character(g1_granges_names))) {
    stop("rownames(g1bins_matrix) must be identical to names(g1_granges).")
  }

  if (!inherits(genome_granges, "GRanges")) {
    stop("genome_granges must be a GRanges object.")
  }
  if (!"is_mappable" %in% colnames(S4Vectors::mcols(genome_granges))) {
    stop("mcols(genome_granges)$is_mappable must exist.")
  }
  genome_granges_names <- names(genome_granges)
  # NOTE: the original tests g1_granges_names here (not genome_granges_names);
  # genome_granges_names is not used afterwards. Kept as is.
  if (is.null(g1_granges_names)) {
    genome_granges_names <- paste0(
      as.character(GenomeInfoDb::seqnames(genome_granges)),
      "-",
      BiocGenerics::start(genome_granges),
      "-",
      BiocGenerics::end(genome_granges)
    )
  }
    
  is_mappable <- as.logical(S4Vectors::mcols(genome_granges)$is_mappable)
  is_mappable[is.na(is_mappable)] <- FALSE

  if (!inherits(BPPARAM, "BiocParallelParam")) {
    stop("'BPPARAM' must be a BiocParallelParam object.")
  }
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L ||
      is.na(chunk_size) || chunk_size < 1L) {
    stop("'chunk_size' must be a positive integer.")
  }
  chunk_size <- as.integer(chunk_size)
  
  if (!is.numeric(gap_max_bp) || length(gap_max_bp) != 1L ||
    is.na(gap_max_bp) || gap_max_bp <= 0L) {
    stop("'gap_max_bp' must be a positive integer.")
  }
  gap_max_bp <- as.integer(gap_max_bp)

  # ------------------------------------------------------------------------
  # 2. Helper functions (now internal: .ars_chr_blocks() etc.)
  # ------------------------------------------------------------------------
  .vmsg(verbose, "[2/6] Building helper functions...")

  # ------------------------------------------------------------------------
  # 3. Prepare chromosome blocks and cell chunks
  # ------------------------------------------------------------------------
  .vmsg(verbose, "[3/6] Preparing chromosome blocks and cell chunks...")

  .ars_check_chr_contiguous(g1_granges, "g1_granges")
  .ars_check_chr_contiguous(genome_granges, "genome_granges")

  g1_map <- .ars_chr_blocks(g1_granges)
  genome_map <- .ars_chr_blocks(genome_granges)

  chrom_boundaries <- g1_map$starts
  num_cells <- ncol(raw_counts)
  cell_ids <- seq_len(num_cells)

  chunks <- split(cell_ids, ceiling(seq_along(cell_ids) / chunk_size))

  # ------------------------------------------------------------------------
  # 4. Parallel inference
  # ------------------------------------------------------------------------
  .vmsg(verbose, "[4/6] Running per-cell replication-state inference...")

  chunk_res <- BiocParallel::bplapply(
    X = chunks,
    FUN = .ars_process_chunk,
    raw_counts = raw_counts,
    g1_granges = g1_granges,
    chrom_boundaries = chrom_boundaries,
    infer_args = infer_args,
    BPPARAM = BPPARAM
  )

  # ------------------------------------------------------------------------
  # 5. Flatten results
  # ------------------------------------------------------------------------
  .vmsg(verbose, "[5/6] Flattening chunk results...")

  flat_filtered_counts <- matrix(
    NA_real_,
    nrow = length(g1_granges),
    ncol = num_cells,
    dimnames = list(names(g1_granges), colnames(raw_counts))
  )

  flat_replication_inference <- matrix(
    NA_integer_,
    nrow = length(g1_granges),
    ncol = num_cells,
    dimnames = list(names(g1_granges), colnames(raw_counts))
  )

  for (res in chunk_res) {
    if (!identical(rownames(res$filtered_chunk), rownames(res$state_chunk))) {
      stop("Row names of filtered_chunk and state_chunk are not identical.")
    }
  
    row_idx <- match(rownames(res$filtered_chunk), rownames(flat_filtered_counts))
    keep_row <- !is.na(row_idx)
    if (!any(keep_row)) next
  
    flat_filtered_counts[row_idx[keep_row], res$cell_idx] <- 
      res$filtered_chunk[keep_row, ,drop = FALSE]
  
    flat_replication_inference[row_idx[keep_row], res$cell_idx] <- 
      res$state_chunk[keep_row, , drop = FALSE]
  }

  flat_replication_inference[is_masked_window, ] <- NA_integer_

  # ------------------------------------------------------------------------
  # 6. Project G1 states to genome windows and fill small unmappable gaps
  # ------------------------------------------------------------------------
  .vmsg(verbose, "[6/6] Projecting states to genome windows...")

  genome_replication_state <- matrix(
    NA_integer_,
    nrow = length(genome_granges),
    ncol = num_cells,
    dimnames = list(names(genome_granges), colnames(raw_counts))
  )

  common_chr <- intersect(g1_map$chr, genome_map$chr)

  for (chr in common_chr) {
    g1_idx <- g1_map$blocks[[chr]]
    genome_idx <- genome_map$blocks[[chr]]

    keep_idx <- g1_idx[!is_masked_window[g1_idx]]
    if (length(keep_idx) == 0L || length(genome_idx) == 0L) {
      next
    }

    hits <- GenomicRanges::findOverlaps(
      query = genome_granges[genome_idx],
      subject = g1_granges[keep_idx],
      ignore.strand = TRUE
    )

    if (length(hits) == 0L) {
      next
    }

    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)

    ov_width <- BiocGenerics::width(
      GenomicRanges::pintersect(
        genome_granges[genome_idx[qh]],
        g1_granges[keep_idx[sh]]
      )
    )

    ord <- order(qh, -ov_width)
    qh <- qh[ord]
    sh <- sh[ord]

    best <- !duplicated(qh)
    qh <- qh[best]
    sh <- sh[best]

    genome_replication_state[genome_idx[qh], ] <-
      flat_replication_inference[keep_idx[sh], , drop = FALSE]

    genome_replication_state[genome_idx[!is_mappable[genome_idx]], ] <-
      NA_integer_

    tmp <- genome_replication_state[genome_idx, , drop = FALSE]
    tmp[!(tmp %in% c(2L, 3L, 4L))] <- NA_integer_
    genome_replication_state[genome_idx, ] <- tmp
  }

  for (chr in common_chr) {
    genome_idx <- genome_map$blocks[[chr]]
    missing_local <- which(!is_mappable[genome_idx])

    if (length(missing_local) == 0L) {
      next
    }

    run_start <- c(1L, which(diff(missing_local) > 1L) + 1L)
    run_end <- c(run_start[-1L] - 1L, length(missing_local))

    for (k in seq_along(run_start)) {
      local_start <- missing_local[run_start[k]]
      local_end <- missing_local[run_end[k]]

      gap_bp <- sum(
        BiocGenerics::width(
          genome_granges[genome_idx[local_start:local_end]]
        )
      )

      if (gap_bp > gap_max_bp) {
        next
      }
      if (local_start == 1L || local_end == length(genome_idx)) {
        next
      }

      fill_range <- seq.int(local_start - 1L, local_end + 1L)
      values <- genome_replication_state[genome_idx[fill_range], , drop = FALSE]

      to_fill <-
        !is.na(values[1L, ]) &
        !is.na(values[nrow(values), ]) &
        values[1L, ] == values[nrow(values), ]

      if (nrow(values) > 2L && any(to_fill)) {
        values[2:(nrow(values) - 1L), to_fill] <-
          matrix(
            values[1L, to_fill],
            nrow = nrow(values) - 2L,
            ncol = sum(to_fill),
            byrow = TRUE
          )
      }

      genome_replication_state[genome_idx[fill_range], ] <- values
    }
  }

  list(
   genobins_rep_states = genome_replication_state,
    g1bins_rep_states = flat_replication_inference
  )
}

#' Run-length chromosome blocks of a GRanges
#' @keywords internal
#' @noRd
.ars_chr_blocks <- function(gr) {
  chr <- as.character(GenomeInfoDb::seqnames(gr))
  rr <- rle(chr)
  starts <- cumsum(c(1L, utils::head(rr$lengths, -1L)))
  ends <- cumsum(rr$lengths)
  blocks <- Map(seq.int, starts, ends)
  names(blocks) <- rr$values

  list(
    chr = rr$values,
    starts = starts,
    ends = ends,
    blocks = blocks
  )
}

#' Stop unless each chromosome forms one contiguous block
#' @keywords internal
#' @noRd
.ars_check_chr_contiguous <- function(gr, object_name) {
  chr <- as.character(GenomeInfoDb::seqnames(gr))
  rr <- rle(chr)

  if (anyDuplicated(rr$values)) {
    dup_chr <- unique(rr$values[duplicated(rr$values)])
    stop(
      "'", object_name, "' must be sorted/grouped by chromosome. ",
      "Non-contiguous chromosome blocks found: ",
      paste(dup_chr, collapse = ", ")
    )
  }
}

#' Run infer_single_cell_states() over a chunk of cells
#' @keywords internal
#' @noRd
.ars_process_chunk <- function(idx, raw_counts, g1_granges, chrom_boundaries, infer_args) {
  filtered_chunk <- matrix(
    NA_real_,
    nrow = length(g1_granges),
    ncol = length(idx),
    dimnames = list(names(g1_granges), colnames(raw_counts)[idx])
  )

  state_chunk <- matrix(
    NA_integer_,
    nrow = length(g1_granges),
    ncol = length(idx),
    dimnames = list(names(g1_granges), colnames(raw_counts)[idx])
  )

  filtered_list <- vector("list", length(idx))
  state_list <- vector("list", length(idx))
  
  names(filtered_list) <- colnames(raw_counts)[idx]
  names(state_list) <- colnames(raw_counts)[idx]

  for (j in seq_along(idx)) {
    barcode <- idx[j]

    res <- do.call(infer_single_cell_states,
           c(list(counts = raw_counts[, barcode],
                  breakpoints = chrom_boundaries),
             infer_args))

    filtered_list[[j]] <- res$filtered_counts
    states <- res$states
    valid <- !is.na(states) & states %in% c(2, 3, 4)
    states <- states[valid]
    state_list[[j]] <- states
  }

  # -------------------------------
  # Batch fill filtered_chunk
  # -------------------------------
  filtered_i <- unlist(
    lapply(filtered_list, function(x) {
      match(names(x), rownames(filtered_chunk))
    }),
    use.names = FALSE
  )
  
  filtered_j <- rep.int(
    seq_along(filtered_list),
    lengths(filtered_list)
  )
  
  filtered_x <- unlist(filtered_list, use.names = FALSE)
  keep <- !is.na(filtered_i)    
  filtered_chunk[cbind(filtered_i[keep], filtered_j[keep])] <- filtered_x[keep]

  # -------------------------------
  # Batch fill state_chunk
  # -------------------------------
  state_i <- unlist(
    lapply(state_list, function(x) {
      match(names(x), rownames(state_chunk))
    }),
    use.names = FALSE
  )
  
  state_j <- rep.int(
    seq_along(state_list),
    lengths(state_list)
  )
  
  state_x <- unlist(state_list, use.names = FALSE)
  keep <- !is.na(state_i) & !is.na(state_x)
  state_chunk[cbind(state_i[keep], state_j[keep])] <- as.integer(state_x[keep])
  
  list(
    cell_idx = idx,
    filtered_chunk = filtered_chunk,
    state_chunk = state_chunk
  )
}
