#' Example datasets
#'
#' Loads one of the two example datasets of the package.
#'
#' `which = "small"` (default): chromosomes 1 to 19 of 40 cells
#' (12 G1 and 28 S-phase) of one mouse ESC 10x multiome sample (`scRT1`,
#' mm10, 20 kb bins). It is meant for examples, the vignette and tests, not
#' for biological conclusions. The script that built it is
#' `data-raw/make_example_data.R` in the package sources.
#'
#' `which = "vignette"`: the same sample and bins for 500 cells (100 called
#' G1 and 400 called non-G1 by the original analysis, picked evenly along
#' their coverage and not filtered on QC, all with RNA), used by the
#' `figures` vignette. It has the same elements except `rep_states`, plus
#' `bulk_rt`: a `GRanges` of the 20 kb bins whose `score` is a bulk Repli-seq
#' replication timing (Loess smoothed; higher = earlier). Built by
#' `data-raw/make_vignette_data.R`.
#'
#' The vignette data (49 MB) are not installed with the package. The first
#' call downloads them from the GitHub release of the package (option
#' `scRTtools.vignette_data_url`, see below), checks their MD5 sum and keeps
#' them in `cache_dir`; later calls read the cached copy.
#'
#' @param which `"small"` or `"vignette"`.
#' @param cache_dir Directory holding the downloaded vignette data; default
#'   `tools::R_user_dir("scRTtools", "cache")`.
#' @param timeout Download timeout in seconds (the `timeout` option is raised
#'   to at least this value during the download).
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
scrt_example_data <- function(which = c("small", "vignette"),
                              cache_dir = tools::R_user_dir("scRTtools", "cache"),
                              timeout = 600) {
  which <- match.arg(which)
  # readRDS() does not load the namespaces of the S4 classes it restores;
  # without this, e.g. names() on the GRanges falls back to base and gives NULL
  loadNamespace("GenomicRanges")
  if (which == "small") {
    path <- system.file("extdata", "scRT1_example.rds", package = "scRTtools")
    if (!nzchar(path)) stop("The example data file scRT1_example.rds is not installed.")
    return(readRDS(path))
  }
  readRDS(.vignette_data_path(cache_dir, timeout))
}

# Release asset holding the vignette data; update all three together when the
# data are rebuilt (data-raw/make_vignette_data.R) and uploaded again.
.vignette_data <- list(
  file = "scRT1_vignette.rds",
  url = "https://github.com/midimagi/scRTtools/releases/download/data-v1/scRT1_vignette.rds",
  md5 = "90fb57adccee3844bcb6880d5b3d3474"
)

.vignette_data_path <- function(cache_dir, timeout) {
  path <- file.path(cache_dir, .vignette_data$file)
  if (file.exists(path) && unname(tools::md5sum(path)) == .vignette_data$md5) {
    return(path)
  }
  url <- getOption("scRTtools.vignette_data_url", .vignette_data$url)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = cache_dir, fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  message("Downloading the vignette data (49 MB) from ", url)
  old <- options(timeout = max(timeout, getOption("timeout")))
  on.exit(options(old), add = TRUE)
  status <- tryCatch(
    utils::download.file(url, tmp, mode = "wb", quiet = TRUE),
    error = function(e) conditionMessage(e)
  )
  if (!identical(status, 0L)) {
    stop("Could not download the vignette data from ", url, ": ",
         if (is.character(status)) status else paste("status", status), call. = FALSE)
  }
  if (unname(tools::md5sum(tmp)) != .vignette_data$md5) {
    stop("The downloaded vignette data have an unexpected MD5 sum; ",
         "the file at ", url, " may have changed.", call. = FALSE)
  }
  if (!file.rename(tmp, path)) stop("Could not write ", path, call. = FALSE)
  path
}
