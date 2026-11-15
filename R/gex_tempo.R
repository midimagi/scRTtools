# Gene expression along replication progress
# (from rcode/analysis/gex_tempo_with_cp10k.R).

#' @keywords internal
#' @noRd
.gex_tempo_check_ids <- function(ids, n, label) {
    if (is.null(ids) || length(ids) != n || anyNA(ids) ||
        any(!nzchar(ids)) || anyDuplicated(ids)) {
        stop(label, " must be non-missing, non-empty and unique.",
             call. = FALSE)
    }
    invisible(TRUE)
}

#' @keywords internal
#' @noRd
.gex_tempo_progress <- function(meta_data, var_use, cell_ids) {
    if (!is.data.frame(meta_data)) {
        stop("meta_data must be a data.frame.", call. = FALSE)
    }
    if (!is.character(var_use) || length(var_use) != 1L ||
        is.na(var_use) || !nzchar(var_use)) {
        stop("var_use must be one non-empty column name.", call. = FALSE)
    }
    .gex_tempo_check_ids(
        rownames(meta_data), nrow(meta_data), "meta_data row names"
    )
    if (!var_use %in% colnames(meta_data)) {
        stop("var_use column '", var_use, "' is absent from meta_data.",
             call. = FALSE)
    }
    matched <- match(cell_ids, rownames(meta_data))
    if (anyNA(matched)) {
        stop("meta_data is missing matrix cell IDs.", call. = FALSE)
    }
    progress <- meta_data[[var_use]][matched]
    if (!is.numeric(progress)) {
        stop("meta_data[[var_use]] must be numeric.", call. = FALSE)
    }
    finite <- is.finite(progress)
    if (any(!finite & !is.na(progress)) ||
        any(progress[finite] < 0 | progress[finite] > 1)) {
        stop("meta_data[[var_use]] must contain values in 0--1 or NA.",
             call. = FALSE)
    }
    if (!any(finite)) {
        stop("No cells have a finite replication-progress value.",
             call. = FALSE)
    }
    unname(progress)
}

#' Summarise gene-expression dynamics along replication progress
#'
#' Cells with finite progress are first ordered by increasing replication
#' progress. They are then grouped either into consecutive groups of `ncells`
#' cells or into `nbin_rep` equal-width progress intervals. Within each
#' non-empty group, method = "sc2p" applies [SC2P_gene_unit()] with warnings
#' suppressed and returns Phase I ratios and Phase II means. Alternatively,
#' method = "cp10k" normalises the entire matrix once by each cell's total
#' counts across all input genes and multiplies by 10,000. It then takes the
#' arithmetic mean across cells in each replication bin, including zero gene
#' expression. No log transformation or phase selection is applied. Supply
#' the full gene count matrix before selecting genes for downstream analysis
#' so library sizes reflect all genes. Sparse input remains sparse during
#' normalisation.
#'
#' Replication progress is read only from `meta_data[[var_use]]`, matched by
#' metadata row names to matrix column names. Values must be in 0--1, as in
#' [rep_tempo()]. Cells with NA/NaN progress are excluded; missing metadata cell IDs,
#' infinite/out-of-range progress, or no usable progress cause an error.
#' Matrix attributes and column ranks are never used as progress.
#'
#' @param gex_matrix Numeric matrix, data.frame, or Matrix-like object with
#'   genes (`gene_id`) in rows and cells (`cell_id`) in columns. Values must be
#'   finite, non-negative integer counts. Row and column names must be unique.
#' @param meta_data Required data.frame with unique cell IDs as row names,
#'   containing every matrix cell and the replication-progress column.
#' @param var_use Column of meta_data containing 0--1 replication progress.
#'   Defaults to "rep_percent".
#' @param ncells Optional positive integer giving the number of consecutive
#'   progress-sorted cells per group. The final smaller group is retained.
#' @param nbin_rep Optional positive integer giving the number of equal-width
#'   replication-progress intervals. Empty intervals are retained as `NA`
#'   columns. `ncells` and `nbin_rep` are mutually exclusive. If no grouping
#'   parameter is supplied, each cell forms its own group.
#'
#' @param yy_thres Count threshold passed to [SC2P_gene_unit()]; ignored for cp10k.
#' @param method Either "sc2p" (default) or "cp10k".
#'
#' @return For sc2p, a named list with `df_phase_I` and `df_phase_II`. Both are data
#'   frames with genes in rows and replication bins in columns. The list has
#'   `rep_bins`, `cell_info`, `progress_source`, and `grouping` attributes.
#'   In non-empty bins, genes with all-zero counts have Phase I ratio 1 and
#'   Phase II mean NA. Empty bins retain NA for both summaries.
#'   For cp10k, a list containing `df_cp10k_mean`, a gene-by-bin data.frame,
#'   with the same grouping and cell metadata plus a `method` attribute.
#'   Zero-library cells are excluded from CP10K means; bins with no positive
#'   libraries retain NA. Otherwise all-zero genes have mean 0. The cp10k
#'   rep_bins attribute also includes `n_cells_used`, the mean denominator.
#'   All input counts are validated for cp10k, including cells with NA progress.
#'   SC2P is not run when method = "cp10k".
#'
#' @examples
#' set.seed(1)
#' counts <- matrix(rpois(30 * 40, 2), nrow = 30,
#'                  dimnames = list(paste0("g", 1:30), paste0("c", 1:40)))
#' meta <- data.frame(rep_percent = runif(40), row.names = colnames(counts))
#' res <- gex_tempo(counts, meta, ncells = 10, method = "cp10k")
#' dim(res$df_cp10k_mean)
#' attr(res, "rep_bins")
#' @export
gex_tempo <- function(gex_matrix, meta_data, var_use = "rep_percent",
                      ncells = NULL, nbin_rep = NULL, yy_thres = 15,
                      method = c("sc2p", "cp10k")) {
    method <- match.arg(method)
    matrix_like <- is.matrix(gex_matrix) || is.data.frame(gex_matrix) ||
        inherits(gex_matrix, "Matrix")
    if (!matrix_like || length(dim(gex_matrix)) != 2L ||
        any(dim(gex_matrix) == 0L)) {
        stop(
            "gex_matrix must be a non-empty numeric matrix, data.frame, ",
            "or Matrix-like object.",
            call. = FALSE
        )
    }
    check_positive_integer <- function(value, label) {
        if (!is.numeric(value) || length(value) != 1L ||
            !is.finite(value) || value < 1 || value != floor(value)) {
            stop(label, " must be NULL or a positive integer.",
                 call. = FALSE)
        }
        as.integer(value)
    }
    if (!is.null(nbin_rep)) {
        nbin_rep <- check_positive_integer(nbin_rep, "nbin_rep")
    }
    if (!is.null(ncells)) {
        ncells <- check_positive_integer(ncells, "ncells")
    }
    if (!is.null(ncells) && !is.null(nbin_rep)) {
        stop("ncells and nbin_rep are mutually exclusive; supply only one.",
             call. = FALSE)
    }
    if (is.null(ncells) && is.null(nbin_rep)) ncells <- 1L

    gene_ids <- rownames(gex_matrix)
    cell_ids <- colnames(gex_matrix)
    .gex_tempo_check_ids(gene_ids, nrow(gex_matrix),
                         "gex_matrix row names")
    .gex_tempo_check_ids(cell_ids, ncol(gex_matrix),
                         "gex_matrix column names")
    if (is.data.frame(gex_matrix) &&
        !all(vapply(gex_matrix, is.numeric, logical(1)))) {
        stop("All gex_matrix columns must be numeric counts.",
             call. = FALSE)
    }
    if (!is.data.frame(gex_matrix) &&
        !(is.numeric(gex_matrix) || inherits(gex_matrix, "dMatrix") ||
          inherits(gex_matrix, "iMatrix"))) {
        stop("All gex_matrix values must be numeric counts.",
             call. = FALSE)
    }

    if (method == "cp10k") {
        # Validate and normalise once, before replication-progress grouping.
        counts <- if (is.data.frame(gex_matrix)) as.matrix(gex_matrix) else gex_matrix
        is_matrix_object <- inherits(counts, "Matrix")
        values <- if (is_matrix_object && "x" %in% methods::slotNames(counts)) {
            methods::slot(counts, "x")
        } else {
            counts
        }
        if (anyNA(values) || any(!is.finite(values)) ||
            any(values < 0) || any(values != floor(values))) {
            stop(
                "gex_matrix must contain finite, non-negative integer counts ",
                "without missing values.", call. = FALSE
            )
        }
        library_sizes <- if (is_matrix_object) Matrix::colSums(counts) else colSums(counts)
        if (any(!is.finite(library_sizes))) {
            stop("Cell library sizes must be finite.", call. = FALSE)
        }
        nonzero_library <- library_sizes > 0
        # Keep zero-library columns at zero for storage, but exclude them
        # from bin means below because their normalised expression is undefined.
        divisors <- library_sizes
        divisors[!nonzero_library] <- 1
        if (is_matrix_object) {
            cp10k <- (counts %*% Matrix::Diagonal(x = 1 / divisors)) * 1e4
            dimnames(cp10k) <- dimnames(counts)
        } else {
            cp10k <- sweep(counts, 2L, divisors, "/") * 1e4
        }
    }

    progress <- .gex_tempo_progress(meta_data, var_use, cell_ids)
    valid_progress <- which(is.finite(progress))
    ordered_cells <- valid_progress[
        order(progress[valid_progress], valid_progress, method = "radix")
    ]
    groups <- rep(NA_integer_, length(progress))

    if (!is.null(nbin_rep)) {
        grouping_mode <- "nbin_rep"
        n_groups <- nbin_rep
        breaks <- seq(0, 1, length.out = n_groups + 1L)
        ordered_groups <- as.integer(cut(
            progress[ordered_cells], breaks = breaks, labels = FALSE,
            include.lowest = TRUE, right = TRUE
        ))
        groups[ordered_cells] <- ordered_groups
        group_indices <- lapply(
            seq_len(n_groups),
            function(i) ordered_cells[ordered_groups == i]
        )
        group_sizes <- lengths(group_indices)
        midpoints <- (breaks[-1L] + breaks[-length(breaks)]) / 2
        group_progress <- vapply(seq_len(n_groups), function(i) {
            if (group_sizes[i]) {
                mean(progress[group_indices[[i]]])
            } else {
                midpoints[i]
            }
        }, numeric(1))
        rep_min <- breaks[-length(breaks)]
        rep_max <- breaks[-1L]
    } else {
        grouping_mode <- "ncells"
        ordered_groups <- as.integer(
            (seq_along(ordered_cells) - 1L) %/% ncells + 1L
        )
        n_groups <- max(ordered_groups)
        groups[ordered_cells] <- ordered_groups
        group_indices <- lapply(
            seq_len(n_groups),
            function(i) ordered_cells[ordered_groups == i]
        )
        group_sizes <- lengths(group_indices)
        group_progress <- vapply(
            group_indices,
            function(i) mean(progress[i]),
            numeric(1)
        )
        rep_min <- vapply(
            group_indices,
            function(i) min(progress[i]),
            numeric(1)
        )
        rep_max <- vapply(
            group_indices,
            function(i) max(progress[i]),
            numeric(1)
        )
    }

    width <- max(4L, nchar(as.character(n_groups)))
    bin_ids <- sprintf(paste0("rep_bin_%0", width, "d"),
                       seq_len(n_groups))
    rep_bins <- data.frame(
        rep_bin = bin_ids,
        rep_percent = 100 * group_progress,
        rep_min = 100 * rep_min,
        rep_max = 100 * rep_max,
        n_cells = group_sizes,
        row.names = NULL,
        check.names = FALSE
    )

    dim_names <- list(gene_ids, bin_ids)
    if (method == "sc2p") {
        phase_I <- matrix(
            NA_real_, nrow = length(gene_ids), ncol = n_groups,
            dimnames = dim_names
        )
        phase_II <- phase_I
        sc2p <- SC2P_gene_unit
    } else {
        cp10k_mean <- matrix(
            NA_real_, nrow = length(gene_ids), ncol = n_groups,
            dimnames = dim_names
        )
        rep_bins$n_cells_used <- integer(n_groups)
    }

    for (i in seq_len(n_groups)) {
        cell_idx <- group_indices[[i]]
        if (!length(cell_idx)) next

        if (method == "cp10k") {
            cell_idx <- cell_idx[nonzero_library[cell_idx]]
            rep_bins$n_cells_used[i] <- length(cell_idx)
            if (length(cell_idx)) {
                bin_cp10k <- cp10k[, cell_idx, drop = FALSE]
                cp10k_mean[, i] <- if (inherits(bin_cp10k, "Matrix")) {
                    Matrix::rowMeans(bin_cp10k)
                } else {
                    rowMeans(bin_cp10k)
                }
            }
            next
        }

        counts <- gex_matrix[, cell_idx, drop = FALSE]
        if (inherits(counts, "Matrix")) counts <- as.matrix(counts)
        if (is.data.frame(counts)) {
            numeric_columns <- vapply(counts, is.numeric, logical(1))
            if (!all(numeric_columns)) {
                stop("All gex_matrix columns must be numeric counts.",
                     call. = FALSE)
            }
        }
        count_values <- as.matrix(counts)
        if (anyNA(count_values) ||
            any(!is.finite(count_values)) ||
            any(count_values < 0) ||
            any(count_values != floor(count_values))) {
            stop(
                "gex_matrix must contain finite, non-negative integer counts ",
                "without missing values.",
                call. = FALSE
            )
        }

        # A gene with zero counts in every cell of a non-empty bin is
        # entirely Phase I. Its Phase II mean remains undefined (NA).
        # Set this before skipping bins with only zero-library cells.
        phase_I[rowSums(count_values) == 0, i] <- 1

        # SC2P_gene_unit() cannot normalise cells whose library size is zero.
        nonzero_cells <- colSums(count_values) > 0
        if (!any(nonzero_cells)) next
        if (!all(nonzero_cells)) {
            counts <- counts[, nonzero_cells, drop = FALSE]
        }

        result <- tryCatch(
            suppressWarnings(sc2p(counts, yy_thres)),
            error = function(e) {
                stop(
                    "SC2P_gene_unit() failed for ", bin_ids[i], ": ",
                    conditionMessage(e),
                    call. = FALSE
                )
            }
        )
        required <- c("gene", "phase_I_ratio", "phase_II_mean")
        if (!is.data.frame(result) || !all(required %in% names(result))) {
            stop(
                "SC2P_gene_unit() must return gene, phase_I_ratio, and ",
                "phase_II_mean columns.",
                call. = FALSE
            )
        }
        result_genes <- as.character(result$gene)
        if (anyNA(result_genes) || any(!nzchar(result_genes)) ||
            anyDuplicated(result_genes)) {
            stop("SC2P_gene_unit() returned invalid gene IDs.",
                 call. = FALSE)
        }
        gene_idx <- match(result_genes, gene_ids)
        if (anyNA(gene_idx)) {
            stop(
                "SC2P_gene_unit() returned genes absent from gex_matrix.",
                call. = FALSE
            )
        }
        phase_I[gene_idx, i] <- result$phase_I_ratio
        phase_II[gene_idx, i] <- result$phase_II_mean
    }

    rep_order <- rep(NA_integer_, length(cell_ids))
    rep_order[ordered_cells] <- seq_along(ordered_cells)
    cell_info <- data.frame(
        cell_id = cell_ids,
        rep_percent = 100 * progress,
        rep_bin = ifelse(is.na(groups), NA_character_, bin_ids[groups]),
        rep_order = rep_order,
        row.names = NULL,
        check.names = FALSE
    )
    if (method == "sc2p") {
        output <- list(
            df_phase_I = as.data.frame(phase_I, check.names = FALSE),
            df_phase_II = as.data.frame(phase_II, check.names = FALSE)
        )
    } else {
        output <- list(
            df_cp10k_mean = as.data.frame(cp10k_mean, check.names = FALSE)
        )
        attr(output, "method") <- method
    }
    attr(output, "rep_bins") <- rep_bins
    attr(output, "cell_info") <- cell_info
    attr(output, "progress_source") <- "meta_data"
    attr(output, "var_use") <- var_use
    attr(output, "grouping") <- list(
        mode = grouping_mode,
        ncells = if (grouping_mode == "ncells") ncells else NULL,
        nbin_rep = if (grouping_mode == "nbin_rep") nbin_rep else NULL
    )
    output
}
