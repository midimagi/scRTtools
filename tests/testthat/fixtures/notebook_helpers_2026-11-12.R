# Verbatim copy of the 5 functions as they were in
# scRT_src/scripts/jupyter/pkg_notebooks/notebook_helpers.R on 2026-11-12,
# before they moved into scRTtools (2026-11-13). Used by test-notebook-utils.R.

cor_safe <- function(x, y) {
    keep <- is.finite(x) & is.finite(y)
    if (sum(keep) < 2L) return(NA_real_)
    cor(x[keep], y[keep])
}

zscore_safe <- function(x) {
    keep <- is.finite(x)
    result <- rep(NA_real_, length(x))
    if (sum(keep) < 2L) return(result)

    spread <- sd(x[keep])
    if (!is.finite(spread) || spread == 0) return(result)

    result[keep] <- (x[keep] - mean(x[keep])) / spread
    result
}

peak_bin <- function(position, value) {
    keep <- is.finite(position) & is.finite(value)
    if (!any(keep)) return(NA_real_)
    position[keep][which.max(value[keep])]
}

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

find_sample_file <- function(data_dir, sample_name, filename) {
    # Search recursively, but only accept a file whose immediate parent
    # directory exactly equals sample_name. This excludes merge/ outputs and
    # prevents accidental partial-name matching.
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
