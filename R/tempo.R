# Migrated from rcode/analysis/tempo.R (rep_tempo) and the inline
# subset_tempo helper of notebook "4.8 tempo SSlogis GAM analysis".

#' Regional replication percentage along genome-wide replication progress
#'
#' For every target region, averages the replication state of the overlapping
#' bins within each cell and then averages cells within groups ordered by
#' genome-wide replication progress. Bins contribute equally whenever they
#' overlap a target (strand is ignored). Missing bins are excluded within each
#' cell; cells without an observed bin are excluded from their group's mean.
#'
#' @param rep_matrix Numeric matrix or `Matrix`, bins by cells, holding
#'   `unrep_state`, `rep_state` and `NA` only. Row and column names must be
#'   unique bin and cell IDs.
#' @param granges `GRanges` of target regions. IDs are taken from the
#'   `range_id` column, `names(granges)`, the `iz_id` column, or the
#'   coordinates, in that order.
#' @param ncells Optional positive integer: consecutive cells per
#'   non-overlapping group (cells must then be sorted by progress). Mutually
#'   exclusive with `nbin_rep`. If neither is given each cell is its own group.
#' @param bin_granges Optional `GRanges` with a `bin_id` column or names
#'   matching the matrix row names. Otherwise the row names are parsed as
#'   `chr-start-end` or `chr:start-end`.
#' @param cell_rep_percent Optional numeric genome-wide progress in 0-100,
#'   named by cell ID or in matrix column order. Defaults to 100 times the
#'   column means of the recoded matrix. Do not combine with `meta_data`.
#' @param meta_data Optional data.frame whose row names are cell IDs; its
#'   `var_use` column (finite, in 0-1) defines genome-wide progress.
#' @param var_use Column of `meta_data` holding progress.
#' @param nbin_rep Optional positive integer: number of equal-width progress
#'   intervals spanning 0-1. Cells in an interval are merged with equal
#'   weight; empty intervals stay as `NA` columns. Requires `meta_data`.
#'   Genome-wide analyses typically use `nbin_rep = 100`.
#' @param unrep_state,rep_state Codes of the unreplicated and replicated
#'   states in `rep_matrix`. The defaults (2 and 4) are the codes written by
#'   [assign_replication_states()]. Use 0 and 1 for a matrix that is already binary.
#'
#' @return A data.frame (regions by progress groups) of percentages in 0-100.
#'   Attributes: `rep_bins`, `range_info`, `cell_info`, `n_valid_cells` and
#'   `units`. Save it with [saveRDS()] to keep the attributes, and subset it
#'   with [subset_rep_tempo()].
#' @seealso [fit_rep_tempo()]
#' @export
rep_tempo <- function(rep_matrix, granges, ncells = NULL,
                      bin_granges = NULL, cell_rep_percent = NULL,
                      meta_data = NULL, var_use = "rep_percent",
                      nbin_rep = NULL, unrep_state = 2, rep_state = 4) {
    for (pkg in c("GenomicRanges", "IRanges", "S4Vectors")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop("Package '", pkg, "' is required.")
        }
    }
    if (!(is.matrix(rep_matrix) || inherits(rep_matrix, "Matrix")) ||
        any(dim(rep_matrix) == 0L)) {
        stop("rep_matrix must be a non-empty matrix or Matrix.")
    }
    check_ids <- function(ids, n, label) {
        if (is.null(ids) || length(ids) != n || anyNA(ids) ||
            any(!nzchar(ids)) || anyDuplicated(ids)) {
            stop(label, " must be non-missing, non-empty and unique.")
        }
    }
    bin_ids <- rownames(rep_matrix)
    cell_ids <- colnames(rep_matrix)
    check_ids(bin_ids, nrow(rep_matrix), "Matrix row names")
    check_ids(cell_ids, ncol(rep_matrix), "Matrix column names")
    if (!inherits(granges, "GRanges")) stop("granges must be a GRanges.")
    if (!is.null(ncells) && !is.null(nbin_rep)) {
        stop("ncells and nbin_rep are mutually exclusive; supply only one.")
    }
    if (is.null(ncells) && is.null(nbin_rep)) ncells <- 1L
    if (!is.null(ncells) &&
        (!is.numeric(ncells) || length(ncells) != 1L || !is.finite(ncells) ||
         ncells < 1 || ncells != floor(ncells))) {
        stop("ncells must be NULL or a positive integer.")
    }
    if (!is.null(nbin_rep) &&
        (!is.numeric(nbin_rep) || length(nbin_rep) != 1L ||
         !is.finite(nbin_rep) || nbin_rep < 1 ||
         nbin_rep != floor(nbin_rep))) {
        stop("nbin_rep must be NULL or a positive integer.")
    }
    if (!is.character(var_use) || length(var_use) != 1L ||
        is.na(var_use) || !nzchar(var_use)) {
        stop("var_use must be one non-empty column name.")
    }
    if (!is.null(meta_data) && !is.data.frame(meta_data)) {
        stop("meta_data must be a data.frame.")
    }
    if (!is.null(meta_data) && !is.null(cell_rep_percent)) {
        stop("meta_data and cell_rep_percent are alternative progress inputs; supply only one.")
    }
    if (!is.null(nbin_rep) && is.null(meta_data)) {
        stop("meta_data is required when nbin_rep is supplied.")
    }
    metadata_progress <- NULL
    if (!is.null(meta_data)) {
        meta_ids <- rownames(meta_data)
        check_ids(meta_ids, nrow(meta_data), "meta_data row names")
        if (!var_use %in% colnames(meta_data)) {
            stop("var_use column '", var_use, "' is absent from meta_data.")
        }
        meta_idx <- match(cell_ids, meta_ids)
        if (anyNA(meta_idx)) {
            stop("meta_data is missing matrix cell IDs.")
        }
        metadata_progress <- meta_data[[var_use]][meta_idx]
        if (!is.numeric(metadata_progress) ||
            any(!is.finite(metadata_progress)) ||
            any(metadata_progress < 0 | metadata_progress > 1)) {
            stop("meta_data[[var_use]] must provide finite values in 0--1 for every matrix cell.")
        }
        metadata_progress <- stats::setNames(100 * metadata_progress, cell_ids)
    }
    if (is.null(bin_granges)) {
        coords <- utils::strcapture(
            "^(.*)[:-]([0-9]+)-([0-9]+)$", bin_ids,
            data.frame(chr = character(), start = integer(), end = integer())
        )
        if (anyNA(coords) || any(!nzchar(coords$chr)) ||
            any(coords$start < 1L) || any(coords$end < coords$start)) {
            stop("Cannot parse bin IDs as chr-start-end or chr:start-end; ",
                 "supply bin_granges with matching IDs.")
        }
        bin_granges <- GenomicRanges::GRanges(
            coords$chr, IRanges::IRanges(coords$start, coords$end)
        )
    } else {
        if (!inherits(bin_granges, "GRanges")) {
            stop("bin_granges must be a GRanges.")
        }
        ids <- if ("bin_id" %in% names(S4Vectors::mcols(bin_granges))) {
            as.character(S4Vectors::mcols(bin_granges)$bin_id)
        } else names(bin_granges)
        check_ids(ids, length(bin_granges), "bin_granges IDs")
        idx <- match(bin_ids, ids)
        if (anyNA(idx)) stop("Some matrix bin IDs are absent from bin_granges.")
        bin_granges <- bin_granges[idx]
    }
    if (any(GenomicRanges::width(bin_granges) < 1L) ||
        any(GenomicRanges::width(granges) < 1L)) {
        stop("Bins and target regions must have positive width.")
    }
    meta <- S4Vectors::mcols(granges)
    range_ids <- if ("range_id" %in% names(meta)) {
        as.character(meta$range_id)
    } else if (!is.null(names(granges))) {
        names(granges)
    } else if ("iz_id" %in% names(meta)) {
        as.character(meta$iz_id)
    } else {
        make.unique(paste0(GenomicRanges::seqnames(granges), ":",
                          GenomicRanges::start(granges), "-",
                          GenomicRanges::end(granges)))
    }
    check_ids(range_ids, length(granges), "Target range IDs")

    # Validate all values in bounded chunks rather than copying a genome-wide
    # dense matrix (the project's merged matrix contains >450 million entries).
    
    # Recode the state matrix to 0 (unreplicated) / 1 (replicated) / NA.
    .rt_state <- function(x, label) {
        if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
            stop(label, " must be a single finite number.")
        }
    }
    .rt_state(unrep_state, "unrep_state")
    .rt_state(rep_state, "rep_state")
    if (unrep_state == rep_state) stop("unrep_state and rep_state must differ.")
    if (unrep_state != 0 || rep_state != 1) {
        is_rep <- rep_matrix == rep_state
        rep_matrix[rep_matrix == unrep_state] <- 0
        rep_matrix[is_rep] <- 1
        rm(is_rep)
    }
    
    
    computed_progress <- numeric(ncol(rep_matrix))
    for (first in seq.int(1L, ncol(rep_matrix), by = 64L)) {
        cols <- seq.int(first, min(first + 63L, ncol(rep_matrix)))
        block <- as.matrix(rep_matrix[, cols, drop = FALSE])
        if (!(is.numeric(block) || is.logical(block)) ||
            any(!is.na(block) & block != 0 & block != 1)) {
            stop("rep_matrix must contain only 0, 1 and NA; convert 2/4 states first.")
        }
        computed_progress[cols] <- 100 * colMeans(block, na.rm = TRUE)
    }
    if (!is.null(metadata_progress)) {
        cell_rep_percent <- metadata_progress
    } else if (is.null(cell_rep_percent)) {
        cell_rep_percent <- computed_progress
    } else if (!is.null(names(cell_rep_percent))) {
        check_ids(names(cell_rep_percent), length(cell_rep_percent),
                  "cell_rep_percent names")
        if (anyNA(match(cell_ids, names(cell_rep_percent)))) {
            stop("cell_rep_percent is missing matrix cell IDs.")
        }
        cell_rep_percent <- cell_rep_percent[cell_ids]
    }
    if (!is.numeric(cell_rep_percent) || length(cell_rep_percent) != ncol(rep_matrix) ||
        any(!is.finite(cell_rep_percent)) ||
        any(cell_rep_percent < 0 | cell_rep_percent > 100)) {
        stop("cell_rep_percent must provide finite progress in 0--100 for every cell.")
    }
    if (!is.null(nbin_rep)) {
        n_groups <- as.integer(nbin_rep)
        breaks <- seq(0, 100, length.out = n_groups + 1L)
        groups <- as.integer(cut(
            unname(cell_rep_percent), breaks = breaks, labels = FALSE,
            include.lowest = TRUE, right = TRUE
        ))
        group_indices <- lapply(seq_len(n_groups), function(i) which(groups == i))
        group_sizes <- lengths(group_indices)
        midpoints <- (breaks[-1L] + breaks[-length(breaks)]) / 2
        group_progress <- vapply(seq_len(n_groups), function(i) {
            if (group_sizes[i]) mean(cell_rep_percent[group_indices[[i]]]) else midpoints[i]
        }, numeric(1))
        rep_bins <- data.frame(
            rep_bin = sprintf("rep_bin_%04d", seq_len(n_groups)),
            rep_percent = group_progress,
            rep_min = breaks[-length(breaks)],
            rep_max = breaks[-1L],
            n_cells = group_sizes, row.names = NULL
        )
    } else {
        if (is.unsorted(unname(cell_rep_percent))) {
            stop("Matrix cells must be sorted by increasing progress when ncells is used.")
        }
        groups <- as.integer((seq_len(ncol(rep_matrix)) - 1) %/% ncells + 1)
        group_indices <- split(seq_along(groups), groups)
        n_groups <- length(group_indices)
        rep_bins <- data.frame(
            rep_bin = sprintf("rep_bin_%04d", seq_len(n_groups)),
            rep_percent = vapply(group_indices, function(i) mean(cell_rep_percent[i]), numeric(1)),
            rep_min = vapply(group_indices, function(i) min(cell_rep_percent[i]), numeric(1)),
            rep_max = vapply(group_indices, function(i) max(cell_rep_percent[i]), numeric(1)),
            n_cells = lengths(group_indices), row.names = NULL
        )
    }
    dim_names <- list(range_ids, rep_bins$rep_bin)
    tempo <- matrix(NA_real_, length(granges), n_groups, dimnames = dim_names)
    n_valid <- matrix(0L, length(granges), n_groups, dimnames = dim_names)
    hits <- GenomicRanges::findOverlaps(granges, bin_granges, ignore.strand = TRUE)
    overlap_bins <- split(S4Vectors::subjectHits(hits),
                          factor(S4Vectors::queryHits(hits),
                                 levels = seq_along(granges)))
    for (i in seq_along(granges)) {
        idx <- overlap_bins[[i]]
        if (!length(idx)) next
        per_cell <- 100 * colMeans(as.matrix(rep_matrix[idx, , drop = FALSE]),
                                  na.rm = TRUE)
        valid <- is.finite(per_cell)
        per_cell[!valid] <- 0
        observed_counts <- rowsum(as.integer(valid), groups, reorder = FALSE)
        observed_sums <- rowsum(per_cell, groups, reorder = FALSE)
        observed_groups <- as.integer(rownames(observed_counts))
        counts <- integer(n_groups)
        sums <- numeric(n_groups)
        counts[observed_groups] <- as.integer(observed_counts[, 1L])
        sums[observed_groups] <- as.numeric(observed_sums[, 1L])
        means <- sums / counts
        means[counts == 0L] <- NA_real_
        tempo[i, ] <- means
        n_valid[i, ] <- counts
    }
    out <- as.data.frame(tempo, check.names = FALSE)
    attr(out, "rep_bins") <- rep_bins
    attr(out, "range_info") <- data.frame(
        range_id = range_ids, chr = as.character(GenomicRanges::seqnames(granges)),
        start = GenomicRanges::start(granges), end = GenomicRanges::end(granges),
        n_overlap_bins = lengths(overlap_bins), row.names = NULL
    )
    attr(out, "cell_info") <- data.frame(
        cell_id = cell_ids, rep_percent = unname(cell_rep_percent),
        rep_bin = rep_bins$rep_bin[groups], row.names = NULL
    )
    attr(out, "n_valid_cells") <- n_valid
    attr(out, "units") <- "percent (0--100)"
    out
}

#' Subset a `rep_tempo` table and keep its attributes
#'
#' Row-subsetting a data.frame drops the custom attributes that
#' [fit_rep_tempo()] needs; this helper carries them along.
#'
#' @param x Output of [rep_tempo()].
#' @param idx Row index (integer, logical or region IDs).
#' @return The subset of `x` with `rep_bins`, `n_valid_cells` and
#'   `range_info` subset to match.
#' @export
subset_rep_tempo <- function(x, idx) {
    if (is.character(idx)) {
        idx <- match(idx, rownames(x))
        if (anyNA(idx)) stop("Some region IDs in 'idx' are absent from 'x'.")
    }
    z <- x[idx, , drop = FALSE]
    attr(z, "rep_bins") <- attr(x, "rep_bins")
    if (!is.null(attr(x, "n_valid_cells"))) {
        attr(z, "n_valid_cells") <- attr(x, "n_valid_cells")[idx, , drop = FALSE]
    }
    if (!is.null(attr(x, "range_info"))) {
        attr(z, "range_info") <- attr(x, "range_info")[idx, , drop = FALSE]
    }
    z
}
