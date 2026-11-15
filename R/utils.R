# Internal helpers shared across modules.

#' Normalise chromosome names to "chr"-prefixed style
#'
#' @param x Character (or coercible) vector of chromosome names.
#' @return Character vector with a single leading `"chr"`.
#' @keywords internal
#' @noRd
.normalize_chr <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  paste0("chr", x)
}

#' Message only when verbose
#' @keywords internal
#' @noRd
.vmsg <- function(verbose, ...) {
  if (isTRUE(verbose)) message(...)
  invisible(NULL)
}

#' Parse a length given as a number or a string with a unit ("50kb", "1Mb")
#' @param x Numeric or character scalar.
#' @return Numeric value in bp.
#' @keywords internal
#' @noRd
.parse_coordinate <- function(x) {
  if (is.numeric(x)) {
    if (length(x) != 1L) {
      stop("Coordinate value must be a scalar.")
    }
    return(as.numeric(x))
  }

  if (!is.character(x) || length(x) != 1L) {
    stop("Coordinate value must be numeric or character, e.g. 50000, '50kb', '1Mb'.")
  }

  x <- gsub("\\s+", "", x)
  value <- as.numeric(gsub("[^0-9.]", "", x))
  unit <- tolower(gsub("[0-9.]", "", x))

  if (is.na(value)) {
    stop("Cannot parse coordinate value: ", x)
  }

  if (unit == "") {
    multiplier <- 1
  } else {
    multiplier <- switch(
      unit,
      bp = 1,
      kb = 1e3,
      mb = 1e6,
      stop("Unsupported length unit: ", unit)
    )
  }

  value * multiplier
}

# as.data.frame() for S4Vectors / GenomicRanges objects (DataFrame, GRanges).
# Under `R CMD check --as-cran`, R >= 4.3 reports the per-column conversion
# done inside S4Vectors as a "direct call of as.data.frame.<type>()"; only
# that deprecation warning is muffled, the result is unchanged.
.mcols_df <- function(x, ...) {
  withCallingHandlers(
    as.data.frame(x, ...),
    deprecatedWarning = function(w) invokeRestart("muffleWarning")
  )
}
