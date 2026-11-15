#' scRTtools: Single-Cell Replication Timing Analysis Tools
#'
#' Analysis functions for single-cell replication timing (no plotting code).
#'
#' @keywords internal
#' @importClassesFrom Matrix dgCMatrix
#' @useDynLib scRTtools, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"

# data.table NSE is used inside package code without importing data.table
# into the namespace, so declare awareness explicitly.
.datatable.aware <- TRUE

utils::globalVariables(c(
  # data.table columns used in define_G1_windows() helpers
  "barcode", "chr", "chr_norm", "contig_id", "end", "pos", "start", "start0",
  "total_width",
  # data.table columns used in annotate_peak_scores()
  "idx", "score", "value", "weight",
  ".N", ":="
))
