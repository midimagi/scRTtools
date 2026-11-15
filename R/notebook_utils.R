# Small utilities shared by the analysis notebooks
# (from scRT_src/scripts/jupyter/pkg_notebooks/notebook_helpers.R, where they
# were collected from inline definitions in notebooks 1.5, 3.1, 4.4 and 6.4).
# The code is unchanged apart from namespaced base/stats calls.

#' NA-tolerant summaries for heatmap tables
#'
#' Helpers used to build the iZone velocity and gene tempo heatmaps. All
#' three drop non-finite values before computing.
#'
#' - `cor_safe()` returns the Pearson correlation of the finite pairs of `x`
#'   and `y`, or `NA_real_` when fewer than 2 pairs are left.
#' - `zscore_safe()` centres and scales the finite values of `x`. Non-finite
#'   entries stay `NA`. The whole result is `NA` when fewer than 2 values are
#'   finite or their standard deviation is 0 or not finite.
#' - `peak_bin()` returns the `position` at which `value` is largest, using only
#'   pairs where both are finite (first one on ties), or `NA_real_` when there
#'   are none.
#'
#' @param x,y Numeric vectors of the same length.
#' @param position Numeric vector of positions (e.g. replication progress bins).
#' @param value Numeric vector of values at `position`.
#' @return `cor_safe()` and `peak_bin()`: a single number. `zscore_safe()`: a
#'   numeric vector as long as `x`.
#' @examples
#' cor_safe(c(1, 2, NA, 4), c(2, 4, 6, Inf))
#' zscore_safe(c(1, 2, 3, NA))
#' peak_bin(c(0, 25, 50, 75), c(0.1, 0.8, NaN, 0.3))
#' @name safe_summaries
NULL

#' @rdname safe_summaries
#' @export
cor_safe <- function(x, y) {
    keep <- is.finite(x) & is.finite(y)
    if (sum(keep) < 2L) return(NA_real_)
    stats::cor(x[keep], y[keep])
}

#' @rdname safe_summaries
#' @export
zscore_safe <- function(x) {
    keep <- is.finite(x)
    result <- rep(NA_real_, length(x))
    if (sum(keep) < 2L) return(result)

    spread <- stats::sd(x[keep])
    if (!is.finite(spread) || spread == 0) return(result)

    result[keep] <- (x[keep] - mean(x[keep])) / spread
    result
}

#' @rdname safe_summaries
#' @export
peak_bin <- function(position, value) {
    keep <- is.finite(position) & is.finite(value)
    if (!any(keep)) return(NA_real_)
    position[keep][which.max(value[keep])]
}

#' Per-sample cell IDs and input files
#'
#' Helpers for merging several samples.
#'
#' `add_sample_prefix()` prepends `"<sample_name>-"` to cell IDs that do not
#' already start with it, so merged samples get unique cell IDs.
#'
#' `find_sample_file()` searches `data_dir` recursively for a file called
#' `filename` whose immediate parent directory is exactly `sample_name`. Files in
#' other directories (e.g. a `merge/` output folder) and partial name matches
#' are ignored. It is an error if there is no match or more than one.
#'
#' @param cell_ids Character vector of cell IDs (coerced with `as.character()`).
#'   `NULL` is an error.
#' @param sample_name Sample name (a single string).
#' @param data_dir Directory to search.
#' @param filename File name to look for (exact match, case-sensitive).
#' @return `add_sample_prefix()`: a character vector as long as `cell_ids`.
#'   `find_sample_file()`: the normalised absolute path of the single match.
#' @examples
#' add_sample_prefix(c("AAAC-1", "s1-AAAG-1"), "s1")
#'
#' d <- tempfile()
#' dir.create(file.path(d, "s1"), recursive = TRUE)
#' dir.create(file.path(d, "merge"))
#' invisible(file.create(file.path(d, "s1", "meta_data.rds"),
#'                       file.path(d, "merge", "meta_data.rds")))
#' find_sample_file(d, "s1", "meta_data.rds")
#' @name sample_files
NULL

#' @rdname sample_files
#' @export
add_sample_prefix <- function(cell_ids, sample_name) {
    if (is.null(cell_ids)) {
        stop("Cell IDs are NULL for sample: ", sample_name)
    }

    cell_ids <- as.character(cell_ids)
    sample_prefix <- paste0(sample_name, "-")

    ifelse(
        startsWith(cell_ids, sample_prefix),
        cell_ids,
        paste0(sample_prefix, cell_ids)
    )
}

#' @rdname sample_files
#' @export
find_sample_file <- function(data_dir, sample_name, filename) {
    candidates <- list.files(
        data_dir,
        full.names = TRUE,
        recursive = TRUE,
        ignore.case = FALSE
    )

    candidates <- candidates[
        basename(candidates) == filename &
        basename(dirname(candidates)) == sample_name
    ]

    if (length(candidates) == 0L) {
        stop(
            "Cannot find ", filename,
            " for sample '", sample_name,
            "' under: ", data_dir
        )
    }

    if (length(candidates) > 1L) {
        stop(
            "Found multiple ", filename,
            " files for sample '", sample_name, "':\n",
            paste(candidates, collapse = "\n")
        )
    }

    normalizePath(candidates, mustWork = TRUE)
}
