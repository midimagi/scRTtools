# S-phase progression marker genes
# (from scRT_learning rcode/processing/identify_S_phase_progression_markers.R).
# The original `config` list is expanded into explicit arguments, the dplyr
# verbs are written in base R and results are returned instead of written.

# -----------------------------------------------------------------------------
# Data loading and sparse normalisation
# -----------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.add_sample_prefix <- function(cell_ids, sample_name) {
    cell_ids <- as.character(cell_ids)
    sample_prefix <- paste0(sample_name, "-")
    ifelse(
        startsWith(cell_ids, sample_prefix),
        cell_ids,
        paste0(sample_prefix, cell_ids)
    )
}

#' @keywords internal
#' @noRd
.read_sample_rna_counts <- function(sample_name, data_dir, file_name, assay) {
    input_file <- file.path(data_dir, sample_name, file_name)
    if (!file.exists(input_file)) {
        stop("Cannot find RNA object for ", sample_name, ": ", input_file)
    }

    wnn <- readRDS(input_file)
    if (!assay %in% SeuratObject::Assays(wnn)) {
        stop("RNA assay named '", assay, "' is absent from: ", input_file)
    }

    # Seurat v5 API; the fallback keeps compatibility with older objects.
    counts <- tryCatch(
        SeuratObject::GetAssayData(wnn, assay = assay, layer = "counts"),
        error = function(e) {
            SeuratObject::GetAssayData(wnn, assay = assay, slot = "counts")
        }
    )

    if (is.null(rownames(counts)) || is.null(colnames(counts))) {
        stop("RNA count matrix lacks gene or cell names for sample: ", sample_name)
    }

    colnames(counts) <- .add_sample_prefix(colnames(counts), sample_name)
    if (anyDuplicated(colnames(counts))) {
        stop("Duplicated RNA cell IDs after prefixing sample: ", sample_name)
    }

    counts
}

#' Load and merge RNA count matrices of several samples
#'
#' Reads one Seurat object per sample from
#' `file.path(data_dir, sample, file_name)`, extracts the raw counts of
#' `assay`, prefixes cell IDs with `"<sample>-"` (unless already prefixed),
#' keeps the genes common to all samples and column-binds the matrices.
#'
#' @param sample_names Character vector of sample names (sub-directories of
#'   `data_dir`).
#' @param data_dir Directory holding one sub-directory per sample.
#' @param file_name Name of the Seurat `.rds` file inside each sample
#'   directory.
#' @param assay Name of the RNA assay.
#'
#' @return A gene-by-cell count matrix (sparse, as stored in the objects).
#' @export
load_and_merge_rna_counts <- function(sample_names, data_dir,
                                      file_name = "wnn_rna_filtered.rds",
                                      assay = "rna") {
    rlang::check_installed("SeuratObject", reason = "to read Seurat objects.")
    counts_by_sample <- stats::setNames(
        lapply(sample_names, .read_sample_rna_counts, data_dir = data_dir,
               file_name = file_name, assay = assay),
        sample_names
    )

    common_genes <- Reduce(intersect, lapply(counts_by_sample, rownames))
    if (length(common_genes) == 0L) {
        stop("No common RNA genes were found across selected samples.")
    }

    counts_by_sample <- lapply(
        counts_by_sample,
        function(x) x[common_genes, , drop = FALSE]
    )
    merged <- do.call(cbind, counts_by_sample)

    if (anyDuplicated(colnames(merged))) {
        duplicated_ids <- unique(colnames(merged)[duplicated(colnames(merged))])
        stop(
            "Duplicated RNA cell IDs after merging (showing up to 20): ",
            paste(utils::head(duplicated_ids, 20), collapse = ", ")
        )
    }

    merged
}

#' Match cell metadata to RNA cells and flag eligible S-phase cells
#'
#' Aligns the merged cell metadata with the RNA cell IDs (prefixing metadata
#' row names with `"<sample>-"` when that matches more cells), optionally
#' derives a `state` column from sample names, and adds the flags
#' `rna_pass_qc`, `dna_pass_qc` and `s_phase_pass`. A cell passes when it
#' passed DNA QC, is not G1, has a finite `rep_percent` within
#' `[min_rep_percent, max_rep_percent]` and a non-missing group.
#'
#' @param meta_data A data.frame of cell metadata with cell IDs as row names
#'   and columns `sample`, `is_G1`, `rep_percent`, `pass_qc`; or the path of
#'   an `.rds` file holding one.
#' @param rna_cell_ids Cell IDs of the RNA count matrix.
#' @param group_by_col Metadata column defining the analysis groups.
#' @param min_rep_percent,max_rep_percent Range of replication progress kept.
#' @param state_patterns Named character vector of regular expressions used
#'   only when `group_by_col` is absent: a `state` column is created whose
#'   value is the name of the first pattern matching `sample` (`NA` when none
#'   matches). `NULL` disables this.
#'
#' @return `meta_data` restricted to (and ordered like) the cells shared with
#'   `rna_cell_ids`, with the three flag columns added.
#' @export
prepare_metadata <- function(meta_data, rna_cell_ids, group_by_col = "state",
                             min_rep_percent = 0.01, max_rep_percent = 0.99,
                             state_patterns = c(ama = "-ama-", ctr = "mESC-rep")) {
    if (is.character(meta_data) && length(meta_data) == 1L) {
        if (!file.exists(meta_data)) {
            stop("Cannot find merged metadata: ", meta_data)
        }
        meta_data <- readRDS(meta_data)
    }
    meta_data <- as.data.frame(meta_data, stringsAsFactors = FALSE)
    required_columns <- c("sample", "is_G1", "rep_percent", "pass_qc")
    missing_columns <- setdiff(required_columns, names(meta_data))
    if (length(missing_columns) > 0L) {
        stop(
            "Merged metadata lacks required columns: ",
            paste(missing_columns, collapse = ", ")
        )
    }

    current_overlap <- length(intersect(rownames(meta_data), rna_cell_ids))
    prefixed_ids <- .add_sample_prefix(rownames(meta_data), meta_data$sample)
    prefixed_overlap <- length(intersect(prefixed_ids, rna_cell_ids))

    if (prefixed_overlap > current_overlap) {
        rownames(meta_data) <- prefixed_ids
        meta_data$cell_id <- prefixed_ids
    }
    if (anyDuplicated(rownames(meta_data))) {
        stop("Merged metadata contains duplicated cell IDs.")
    }

    if (!group_by_col %in% names(meta_data) && length(state_patterns)) {
        if (is.null(names(state_patterns)) || any(!nzchar(names(state_patterns)))) {
            stop("state_patterns must be a named character vector.")
        }
        state <- rep(NA_character_, nrow(meta_data))
        for (i in seq_along(state_patterns)) {
            hit <- is.na(state) & grepl(state_patterns[[i]], meta_data$sample)
            state[hit] <- names(state_patterns)[i]
        }
        meta_data$state <- state
    }
    if (!group_by_col %in% names(meta_data)) {
        stop("Metadata lacks grouping column: ", group_by_col)
    }

    common_cells <- intersect(rna_cell_ids, rownames(meta_data))
    if (length(common_cells) == 0L) {
        stop("No matching cells were found between RNA counts and metadata.")
    }

    meta_data <- meta_data[common_cells, , drop = FALSE]
    meta_data$rna_pass_qc <- TRUE
    meta_data$dna_pass_qc <- !is.na(meta_data$pass_qc) & meta_data$pass_qc
    meta_data$s_phase_pass <- (
        meta_data$rna_pass_qc &
        meta_data$dna_pass_qc &
        !is.na(meta_data$is_G1) &
        !meta_data$is_G1 &
        is.finite(meta_data$rep_percent) &
        meta_data$rep_percent >= min_rep_percent &
        meta_data$rep_percent <= max_rep_percent &
        !is.na(meta_data[[group_by_col]])
    )

    meta_data
}

#' Library-size normalise a sparse count matrix
#'
#' @param counts Sparse gene-by-cell count matrix (`dgCMatrix`).
#' @param scale_factor Counts are scaled to this total per cell.
#'
#' @return A list with `normalized` (counts per `scale_factor`),
#'   `log_normalized` (`log1p` of it) and `library_size`. Cells with a zero or
#'   non-finite library size are an error.
#'
#' @examples
#' m <- Matrix::rsparsematrix(20, 10, 0.4, rand.x = function(n) rpois(n, 3) + 1)
#' dimnames(m) <- list(paste0("g", 1:20), paste0("c", 1:10))
#' str(normalize_sparse_counts(m), max.level = 1)
#' @export
normalize_sparse_counts <- function(counts, scale_factor = 10000) {
    library_size <- Matrix::colSums(counts)
    if (any(!is.finite(library_size) | library_size <= 0)) {
        stop("RNA count matrix contains cells with zero or invalid library size.")
    }

    normalized <- counts %*% Matrix::Diagonal(x = scale_factor / library_size)
    dimnames(normalized) <- dimnames(counts)
    log_normalized <- normalized
    log_normalized@x <- log1p(log_normalized@x)

    list(
        normalized = normalized,
        log_normalized = log_normalized,
        library_size = library_size
    )
}

# -----------------------------------------------------------------------------
# GAM fitting and model summaries
# -----------------------------------------------------------------------------

#' Fit the progression GAM and its null model for one gene
#'
#' Fits `expression ~ s(rep_percent, k)` (optionally `+ sample`) with
#' `mgcv::gam` (Gaussian, GCV) and the matching null model without the
#' smooth. `k = max(4, min(gam_basis_k, n_unique_rep - 1))`.
#'
#' @param y Numeric expression vector, one value per row of `model_data`.
#' @param model_data data.frame with `rep_percent` and a factor `sample`.
#' @param gam_basis_k Upper bound of the smooth basis dimension.
#' @param use_sample_covariate Add `sample` as a covariate when it has more
#'   than one level.
#'
#' @return A list with `gam`, `null`, `model_data` (with `expression` added)
#'   and `sample_is_used`. Errors when there are fewer than five unique
#'   `rep_percent` values or the expression is constant.
#' @export
build_gene_models <- function(y, model_data, gam_basis_k = 10L,
                              use_sample_covariate = FALSE) {
    rlang::check_installed("mgcv", reason = "to fit progression GAMs.")
    # mgcv::gam() resolves s() from the formula environment; mgcv is only
    # suggested, so bind it here instead of attaching mgcv.
    s <- mgcv::s

    model_data$expression <- as.numeric(y)
    n_unique_rep <- length(unique(model_data$rep_percent))
    if (n_unique_rep < 5L) {
        stop("Fewer than five unique replication percentages.")
    }

    expression_variance <- stats::var(model_data$expression)
    if (!is.finite(expression_variance) || expression_variance <= 0) {
        stop("Expression is constant or invalid.")
    }

    k_use <- max(4L, min(as.integer(gam_basis_k), n_unique_rep - 1L))
    sample_is_used <- (
        isTRUE(use_sample_covariate) &&
        "sample" %in% names(model_data) &&
        nlevels(droplevels(model_data$sample)) > 1L
    )

    if (sample_is_used) {
        model_data$sample <- droplevels(model_data$sample)
        gam_fit <- mgcv::gam(
            expression ~ sample + s(rep_percent, k = k_use),
            data = model_data,
            family = stats::gaussian(link = "identity"),
            method = "GCV.Cp"
        )
        null_fit <- stats::lm(expression ~ sample, data = model_data)
    } else {
        gam_fit <- mgcv::gam(
            expression ~ s(rep_percent, k = k_use),
            data = model_data,
            family = stats::gaussian(link = "identity"),
            method = "GCV.Cp"
        )
        null_fit <- stats::lm(expression ~ 1, data = model_data)
    }

    list(
        gam = gam_fit,
        null = null_fit,
        model_data = model_data,
        sample_is_used = sample_is_used
    )
}

#' Summarise the models of one gene
#'
#' @param model_list Output of [build_gene_models()].
#' @return A one-row data.frame: smooth-term `p_value`, `edf`, `f_statistic`,
#'   `adjusted_r_squared`, `deviance_explained`, `aic_gam`, `aic_null` and
#'   `sample_covariate_used`.
#' @export
summarize_gene_models <- function(model_list) {
    gam_summary <- summary(model_list$gam)
    smooth_table <- gam_summary$s.table
    if (is.null(smooth_table) || nrow(smooth_table) != 1L) {
        stop("Expected exactly one smooth term.")
    }

    data.frame(
        p_value = as.numeric(smooth_table[1, ncol(smooth_table)]),
        edf = as.numeric(smooth_table[1, "edf"]),
        f_statistic = if ("F" %in% colnames(smooth_table)) {
            as.numeric(smooth_table[1, "F"])
        } else {
            NA_real_
        },
        adjusted_r_squared = unname(gam_summary$r.sq),
        deviance_explained = unname(gam_summary$dev.expl),
        aic_gam = stats::AIC(model_list$gam),
        aic_null = stats::AIC(model_list$null),
        sample_covariate_used = model_list$sample_is_used,
        stringsAsFactors = FALSE
    )
}

#' Predict the fitted expression curve of one gene
#'
#' Fits the GAM of [build_gene_models()] and predicts on `rep_grid`. With a
#' sample covariate the prediction is averaged over samples.
#'
#' @inheritParams build_gene_models
#' @param rep_grid Replication-progress values to predict at.
#' @return Numeric vector, one value per element of `rep_grid`.
#' @export
predict_marginal_curve <- function(y, model_data, rep_grid, gam_basis_k = 10L,
                                   use_sample_covariate = FALSE) {
    model_list <- build_gene_models(
        y = y,
        model_data = model_data,
        gam_basis_k = gam_basis_k,
        use_sample_covariate = use_sample_covariate
    )

    if (model_list$sample_is_used) {
        sample_levels <- levels(model_list$model_data$sample)
        new_data <- expand.grid(
            rep_percent = rep_grid,
            sample = sample_levels,
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
        )
        new_data$sample <- factor(new_data$sample, levels = sample_levels)
        predicted <- stats::predict(model_list$gam, newdata = new_data)
        predicted <- matrix(
            predicted,
            nrow = length(rep_grid),
            ncol = length(sample_levels)
        )
        rowMeans(predicted)
    } else {
        as.numeric(stats::predict(
            model_list$gam,
            newdata = data.frame(rep_percent = rep_grid)
        ))
    }
}

#' Row-wise z-scores, dropping constant rows
#'
#' @param x Numeric matrix.
#' @return A list with `matrix` (z-scored rows with a finite, positive
#'   standard deviation) and `valid` (logical, which rows of `x` were kept).
#' @examples
#' scale_rows(rbind(a = c(1, 2, 3), b = c(5, 5, 5)))
#' @export
scale_rows <- function(x) {
    centered <- sweep(x, 1, rowMeans(x), FUN = "-")
    row_sd <- sqrt(rowSums(centered^2) / pmax(ncol(x) - 1L, 1L))
    valid <- is.finite(row_sd) & row_sd > 0
    list(
        matrix = sweep(
            centered[valid, , drop = FALSE],
            1,
            row_sd[valid],
            FUN = "/"
        ),
        valid = valid
    )
}

# -----------------------------------------------------------------------------
# Group analysis and marker classification
# -----------------------------------------------------------------------------

#' Collect the cells and matrices of one analysis group
#'
#' @param group_name Value of `group_by_col` to analyse.
#' @param meta_data Output of [prepare_metadata()].
#' @param log_expr,normalized_expr Gene-by-cell matrices from
#'   [normalize_sparse_counts()].
#' @param group_by_col Metadata column defining the analysis groups.
#' @param min_cells_per_group Minimum number of eligible cells; fewer is an
#'   error.
#'
#' @return A list with `group_meta` (sorted by `rep_percent`), `model_data`
#'   (`rep_percent`, `sample`), `log_expr`, `normalized_expr` and `cells`.
#' @export
prepare_group_inputs <- function(group_name, meta_data, log_expr, normalized_expr,
                                 group_by_col = "state",
                                 min_cells_per_group = 50L) {
    group_col <- group_by_col
    if (!group_col %in% names(meta_data)) {
        stop("Metadata lacks grouping column: ", group_col)
    }

    group_cells <- rownames(meta_data)[
        meta_data$s_phase_pass &
        as.character(meta_data[[group_col]]) == group_name
    ]
    if (length(group_cells) < min_cells_per_group) {
        stop(
            "Group '", group_name, "' has only ", length(group_cells),
            " eligible S-phase cells; minimum is ",
            min_cells_per_group, "."
        )
    }
    if (!all(group_cells %in% colnames(log_expr)) ||
        !all(group_cells %in% colnames(normalized_expr))) {
        stop("Some group cells are absent from the expression matrices.")
    }

    group_meta <- meta_data[group_cells, , drop = FALSE]
    group_meta$sample <- droplevels(factor(group_meta$sample))
    group_meta <- group_meta[order(group_meta$rep_percent), , drop = FALSE]
    group_cells <- rownames(group_meta)

    model_data <- data.frame(
        rep_percent = group_meta$rep_percent,
        sample = group_meta$sample,
        row.names = group_cells
    )
    group_log_expr <- log_expr[, group_cells, drop = FALSE]
    group_normalized_expr <- normalized_expr[, group_cells, drop = FALSE]

    if (!identical(rownames(group_log_expr), rownames(group_normalized_expr))) {
        stop("Normalized and log-normalized matrices have different genes.")
    }

    list(
        group_meta = group_meta,
        model_data = model_data,
        log_expr = group_log_expr,
        normalized_expr = group_normalized_expr,
        cells = group_cells
    )
}

#' Fit progression GAMs for every gene
#'
#' Applies [build_gene_models()] and [summarize_gene_models()] to each row of
#' `group_log_expr`. A gene whose fit fails gets `NA` statistics and the
#' error message in `fit_status`.
#'
#' @param group_log_expr,group_normalized_expr Gene-by-cell matrices of one
#'   group, columns aligned with the rows of `model_data`.
#' @param model_data data.frame from [prepare_group_inputs()].
#' @param min_normalized_expression A gene counts as expressed in a cell when
#'   its normalised expression is at least this value.
#' @inheritParams build_gene_models
#' @param BPPARAM A `BiocParallel` back-end; serial by default.
#' @param verbose Print a progress message. Default `FALSE`.
#'
#' @return A data.frame with one row per gene: the columns of
#'   [summarize_gene_models()], `fit_status`, `gene` and `expressed_fraction`.
#' @export
fit_gene_gams <- function(group_log_expr, group_normalized_expr, model_data,
                          min_normalized_expression = 0.1, gam_basis_k = 10L,
                          use_sample_covariate = FALSE,
                          BPPARAM = BiocParallel::SerialParam(),
                          verbose = FALSE) {
    if (ncol(group_log_expr) != nrow(model_data) ||
        !identical(colnames(group_log_expr), rownames(model_data))) {
        stop("Expression columns and model_data rows are not aligned.")
    }

    expressed_fraction <- Matrix::rowMeans(
        group_normalized_expr >= min_normalized_expression
    )

    fit_one_index <- function(gene_index) {
        tryCatch(
            {
                model_list <- suppressWarnings(build_gene_models(
                    y = group_log_expr[gene_index, ],
                    model_data = model_data,
                    gam_basis_k = gam_basis_k,
                    use_sample_covariate = use_sample_covariate
                ))
                result <- summarize_gene_models(model_list)
                result$fit_status <- "ok"
                result
            },
            error = function(e) {
                data.frame(
                    p_value = NA_real_,
                    edf = NA_real_,
                    f_statistic = NA_real_,
                    adjusted_r_squared = NA_real_,
                    deviance_explained = NA_real_,
                    aic_gam = NA_real_,
                    aic_null = NA_real_,
                    sample_covariate_used = NA,
                    fit_status = conditionMessage(e),
                    stringsAsFactors = FALSE
                )
            }
        )
    }

    gene_indices <- seq_len(nrow(group_log_expr))
    .vmsg(
        verbose,
        "Fitting ", length(gene_indices), " genes in ",
        ncol(group_log_expr), " cells..."
    )

    fit_list <- BiocParallel::bplapply(
        X = gene_indices,
        FUN = fit_one_index,
        BPPARAM = BPPARAM
    )

    fit_table <- .bind_rows_df(fit_list)
    fit_table$gene <- rownames(group_log_expr)
    fit_table$expressed_fraction <- as.numeric(
        expressed_fraction[fit_table$gene]
    )
    fit_table
}

#' Classify genes as S-phase progression markers
#'
#' Adds BH-adjusted `fdr` (non-finite p-values count as 1), the AIC
#' difference `delta_aic_null_minus_gam` and the flags `aic_pass`,
#' `expression_pass` (`expressed_fraction >= min_cell_fraction`), `dynamic`
#' (`fdr < fdr_threshold` and `aic_pass`) and `marker` (`dynamic` and
#' `expression_pass`).
#'
#' @param fit_table Output of [fit_gene_gams()].
#' @param group_name,n_cells Group label and cell number, stored in the
#'   `analysis_group` and `n_cells` columns.
#' @param min_cell_fraction Minimum fraction of cells expressing the gene.
#' @param fdr_threshold FDR cut-off for the smooth term.
#' @param aic_rule `"gam_lower"`: the GAM must have the lower AIC;
#'   `"paper_literal_gam_greater"`: the opposite.
#'
#' @return `fit_table` with the added columns, `gene`, `analysis_group`,
#'   `n_cells` first, sorted by `fdr` then decreasing AIC difference.
#' @export
classify_progression_markers <- function(fit_table, group_name, n_cells,
                                         min_cell_fraction = 0.2,
                                         fdr_threshold = 0.01,
                                         aic_rule = c("gam_lower",
                                                      "paper_literal_gam_greater")) {
    aic_rule <- match.arg(aic_rule)
    p_for_adjustment <- fit_table$p_value
    p_for_adjustment[!is.finite(p_for_adjustment)] <- 1
    fit_table$fdr <- stats::p.adjust(p_for_adjustment, method = "BH")
    fit_table$delta_aic_null_minus_gam <- fit_table$aic_null - fit_table$aic_gam

    if (aic_rule == "gam_lower") {
        fit_table$aic_pass <- fit_table$delta_aic_null_minus_gam > 0
    } else {
        fit_table$aic_pass <- fit_table$delta_aic_null_minus_gam < 0
    }

    fit_table$expression_pass <- (
        fit_table$expressed_fraction >= min_cell_fraction
    )
    fit_table$dynamic <- (
        fit_table$fdr < fdr_threshold & fit_table$aic_pass
    )
    fit_table$marker <- fit_table$dynamic & fit_table$expression_pass
    fit_table$analysis_group <- group_name
    fit_table$n_cells <- n_cells
    fit_table$aic_rule <- aic_rule

    first <- c("gene", "analysis_group", "n_cells")
    fit_table <- fit_table[, c(first, setdiff(names(fit_table), first)), drop = FALSE]
    fit_table <- fit_table[
        order(fit_table$fdr, -fit_table$delta_aic_null_minus_gam), , drop = FALSE
    ]
    rownames(fit_table) <- NULL
    fit_table
}

# -----------------------------------------------------------------------------
# Curve calculation and complete per-group analysis
# -----------------------------------------------------------------------------

#' Fitted expression curves of selected genes and their clustering
#'
#' Predicts each gene's GAM curve on an even grid over the observed
#' replication progress, z-scores the curves, clusters them (Ward's method on
#' `1 - Pearson correlation`) and summarises each curve.
#'
#' @param group_log_expr Gene-by-cell log-normalised matrix of one group.
#' @param model_data data.frame from [prepare_group_inputs()].
#' @param genes Genes to compute curves for; unknown genes are ignored.
#' @inheritParams build_gene_models
#' @param n_curve_points Number of grid points.
#' @param n_clusters Number of curve clusters (capped at the number of
#'   curves).
#'
#' @return A list with `curves` and `curve_z` (gene-by-grid matrices),
#'   `rep_grid`, `clusters`, `hclust` (`NULL` for a single curve),
#'   `gene_order` and `curve_summary` (`gene`, `cluster`, `peak_rep_percent`,
#'   `late_minus_early`). Genes whose curve is non-finite or constant are
#'   dropped. With no usable gene name the elements are `NULL` / empty.
#' @export
compute_gene_curves <- function(group_log_expr, model_data, genes,
                                gam_basis_k = 10L, use_sample_covariate = FALSE,
                                n_curve_points = 101L, n_clusters = 3L) {
    genes <- intersect(as.character(genes), rownames(group_log_expr))
    if (length(genes) == 0L) {
        return(list(
            curves = NULL,
            curve_z = NULL,
            rep_grid = NULL,
            clusters = NULL,
            hclust = NULL,
            gene_order = character(0),
            curve_summary = data.frame()
        ))
    }

    rep_grid <- seq(
        min(model_data$rep_percent),
        max(model_data$rep_percent),
        length.out = n_curve_points
    )
    curve_list <- lapply(genes, function(gene) {
        tryCatch(
            suppressWarnings(predict_marginal_curve(
                y = group_log_expr[gene, ],
                model_data = model_data,
                rep_grid = rep_grid,
                gam_basis_k = gam_basis_k,
                use_sample_covariate = use_sample_covariate
            )),
            error = function(e) rep(NA_real_, length(rep_grid))
        )
    })

    curve_matrix <- do.call(rbind, curve_list)
    rownames(curve_matrix) <- genes
    colnames(curve_matrix) <- sprintf("rep_%0.4f", rep_grid)
    valid_curve <- apply(curve_matrix, 1, function(x) all(is.finite(x)))
    curve_matrix <- curve_matrix[valid_curve, , drop = FALSE]

    if (nrow(curve_matrix) == 0L) {
        stop("All requested GAM curves were invalid.")
    }

    scaled_curve <- scale_rows(curve_matrix)
    curve_matrix <- curve_matrix[scaled_curve$valid, , drop = FALSE]
    curve_z <- scaled_curve$matrix
    if (nrow(curve_z) == 0L) {
        stop("All requested GAM curves were constant.")
    }

    if (nrow(curve_z) == 1L) {
        hc <- NULL
        clusters <- stats::setNames(1L, rownames(curve_z))
        gene_order <- rownames(curve_z)
    } else {
        gene_correlation <- stats::cor(t(curve_z), method = "pearson")
        # Do not use pmin/pmax here: with a scalar first argument they can
        # drop matrix dimensions and make as.dist() receive a vector.
        gene_correlation[gene_correlation > 1] <- 1
        gene_correlation[gene_correlation < -1] <- -1
        correlation_distance <- stats::as.dist(1 - gene_correlation)
        hc <- stats::hclust(correlation_distance, method = "ward.D2")
        clusters <- stats::cutree(
            hc,
            k = min(as.integer(n_clusters), nrow(curve_z))
        )
        gene_order <- rownames(curve_z)[hc$order]
    }

    peak_index <- max.col(curve_matrix, ties.method = "first")
    curve_summary <- data.frame(
        gene = rownames(curve_matrix),
        cluster = as.integer(clusters[rownames(curve_matrix)]),
        peak_rep_percent = rep_grid[peak_index],
        late_minus_early = curve_matrix[, ncol(curve_matrix)] - curve_matrix[, 1],
        stringsAsFactors = FALSE
    )

    list(
        curves = curve_matrix,
        curve_z = curve_z,
        rep_grid = rep_grid,
        clusters = clusters,
        hclust = hc,
        gene_order = gene_order,
        curve_summary = curve_summary
    )
}

#' Identify S-phase progression markers in one group of cells
#'
#' Runs [prepare_group_inputs()], [fit_gene_gams()],
#' [classify_progression_markers()] and [compute_gene_curves()] (on the
#' marker genes) for one analysis group.
#'
#' @inheritParams prepare_group_inputs
#' @inheritParams fit_gene_gams
#' @inheritParams classify_progression_markers
#' @inheritParams compute_gene_curves
#'
#' @return A list with `analysis_group`, `gam_results` (all genes),
#'   `marker_table` (markers with a valid curve, with `cluster`,
#'   `peak_rep_percent`, `late_minus_early`, sorted by them and `fdr`),
#'   `group_meta`, `cells`, and the `curves`, `curve_z`, `rep_grid`,
#'   `clusters`, `hclust`, `gene_order` elements of [compute_gene_curves()].
#' @export
analyze_one_group <- function(group_name, meta_data, log_expr, normalized_expr,
                              group_by_col = "state", min_cells_per_group = 50L,
                              min_normalized_expression = 0.1,
                              min_cell_fraction = 0.2, fdr_threshold = 0.01,
                              aic_rule = c("gam_lower", "paper_literal_gam_greater"),
                              gam_basis_k = 10L, use_sample_covariate = FALSE,
                              n_curve_points = 101L, n_clusters = 3L,
                              BPPARAM = BiocParallel::SerialParam(),
                              verbose = FALSE) {
    aic_rule <- match.arg(aic_rule)
    .vmsg(verbose, "\n=== Compute group: ", group_name, " ===")
    inputs <- prepare_group_inputs(
        group_name = group_name,
        meta_data = meta_data,
        log_expr = log_expr,
        normalized_expr = normalized_expr,
        group_by_col = group_by_col,
        min_cells_per_group = min_cells_per_group
    )

    fit_table <- fit_gene_gams(
        group_log_expr = inputs$log_expr,
        group_normalized_expr = inputs$normalized_expr,
        model_data = inputs$model_data,
        min_normalized_expression = min_normalized_expression,
        gam_basis_k = gam_basis_k,
        use_sample_covariate = use_sample_covariate,
        BPPARAM = BPPARAM,
        verbose = verbose
    )
    gam_results <- classify_progression_markers(
        fit_table = fit_table,
        group_name = group_name,
        n_cells = length(inputs$cells),
        min_cell_fraction = min_cell_fraction,
        fdr_threshold = fdr_threshold,
        aic_rule = aic_rule
    )

    marker_genes <- gam_results$gene[gam_results$marker]
    .vmsg(
        verbose,
        "Dynamic genes before expression filter: ", sum(gam_results$dynamic),
        "; final markers: ", length(marker_genes)
    )
    curve_result <- compute_gene_curves(
        group_log_expr = inputs$log_expr,
        model_data = inputs$model_data,
        genes = marker_genes,
        gam_basis_k = gam_basis_k,
        use_sample_covariate = use_sample_covariate,
        n_curve_points = n_curve_points,
        n_clusters = n_clusters
    )

    if (is.null(curve_result$curves)) {
        marker_table <- gam_results[gam_results$marker, , drop = FALSE]
        marker_table$cluster <- integer(0)
        marker_table$peak_rep_percent <- numeric(0)
        marker_table$late_minus_early <- numeric(0)
    } else {
        marker_table <- gam_results[
            gam_results$gene %in% rownames(curve_result$curves), , drop = FALSE
        ]
        marker_table <- .left_join_by_id(
            marker_table, curve_result$curve_summary, by = "gene"
        )
        marker_table <- marker_table[
            order(marker_table$cluster, marker_table$peak_rep_percent,
                  marker_table$fdr), , drop = FALSE
        ]
        rownames(marker_table) <- NULL
    }

    c(
        list(
            analysis_group = group_name,
            gam_results = gam_results,
            marker_table = marker_table,
            group_meta = inputs$group_meta,
            cells = inputs$cells
        ),
        curve_result[c(
            "curves", "curve_z", "rep_grid", "clusters",
            "hclust", "gene_order"
        )]
    )
}
