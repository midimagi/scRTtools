# Internal input validators. They only check structure and stop() with an
# informative message; they never modify their input.

#' Validate a GRanges of genome bins
#' @keywords internal
#' @noRd
.validate_bins <- function(x, arg = "genome_bins", allow_empty = FALSE) {
  if (!methods::is(x, "GRanges")) {
    stop("'", arg, "' must be a GRanges object.", call. = FALSE)
  }
  if (!allow_empty && length(x) == 0L) {
    stop("'", arg, "' must be non-empty.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate a per-chromosome GRangesList of windows
#'
#' @param x Object to check.
#' @param arg Argument name used in error messages.
#' @param named Require non-empty, non-NA element names.
#' @param mcol Optional metadata column that must exist in every element.
#' @param numeric_mcol If `TRUE`, `mcol` must also be numeric.
#' @keywords internal
#' @noRd
.validate_grl <- function(x, arg = "wins_grl", named = FALSE,
                          mcol = NULL, numeric_mcol = FALSE) {
  if (!methods::is(x, "GRangesList")) {
    stop("'", arg, "' must be a GenomicRanges::GRangesList object.",
         call. = FALSE)
  }
  if (length(x) == 0L) {
    stop("'", arg, "' must not be empty.", call. = FALSE)
  }
  if (named) {
    nm <- names(x)
    if (is.null(nm) || anyNA(nm) || any(nm == "")) {
      stop("'", arg, "' must be a named GRangesList.", call. = FALSE)
    }
  }
  if (!is.null(mcol)) {
    ok <- vapply(seq_along(x), function(i) {
      mc <- S4Vectors::mcols(x[[i]])
      mcol %in% colnames(mc) && (!numeric_mcol || is.numeric(mc[[mcol]]))
    }, logical(1))
    if (!all(ok)) {
      stop("'", mcol, "' must exist", if (numeric_mcol) " and be numeric",
           " in all elements of '", arg, "'. Missing in element(s): ",
           paste(which(!ok), collapse = ", "), call. = FALSE)
    }
  }
  invisible(TRUE)
}
