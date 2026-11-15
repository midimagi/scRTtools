#' Small example dataset
#'
#' Loads the example shipped with the package: chromosomes 1 to 19 of 40 cells
#' (12 G1 and 28 S-phase) of one mouse ESC 10x multiome sample (`scRT1`,
#' mm10, 20 kb bins). It is meant for examples, the vignette and tests, not
#' for biological conclusions. The script that built it is
#' `data-raw/make_example_data.R` in the package sources.
#'
#' @return A list with
#' \describe{
#'   \item{`bins`}{`GRanges` of the 123,148 20 kb bins of chr1-chr19 with
#'     `gc_frequency`, `mappability` and `is_mappable`.}
#'   \item{`bin_counts`}{Integer matrix of raw read counts, bins by cells.}
#'   \item{`g1_windows`}{`GRanges` of the 101,297 variable-width G1 windows
#'     of chr1-chr19 (about 200 G1 reads each) as produced by [define_G1_windows()].}
#'   \item{`window_counts`}{Integer matrix of raw read counts, G1 windows by
#'     cells.}
#'   \item{`rep_states`}{Integer matrix of replication states from the
#'     original analysis (2 = unreplicated, 4 = replicated, `NA` = no call),
#'     bins by the 28 S-phase cells.}
#'   \item{`rna_counts`}{Sparse matrix of RNA counts of the 24,860 chr1-chr19
#'     genes, genes by cells. Use `meta$nCount_rna` as the library size:
#'     these columns hold only the chr1-chr19 genes (about 87% of the
#'     counts).}
#'   \item{`genes`}{data.frame with `gene`, `chr`, `start`, `end`, `strand`.}
#'   \item{`meta`}{data.frame of cell metadata (row names = cell IDs):
#'     `is_G1`, `rep_percent` (genome-wide fraction replicated), `coverage`,
#'     `mapd`, RNA depth and QC flags.}
#'   \item{`info`}{Sample, chromosomes, genome, bin width and creation date.}
#' }
#'
#' @examples
#' ex <- scrt_example_data()
#' dim(ex$bin_counts)
#' table(ex$meta$is_G1)
#' @export
scrt_example_data <- function() {
  path <- system.file("extdata", "scRT1_example.rds", package = "scRTtools")
  if (!nzchar(path)) stop("The example data file is not installed.")
  readRDS(path)
}
