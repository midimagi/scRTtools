# Coverage matching and fragment downsampling between conditions.
# optimize_shrink_factor() is extracted from the notebook
# jupyter_downsample_frag_quantile/1.2.1 extract quantile cells (optimize_shrink);
# downsample_fragments() wraps inst/scripts/downsample_fragments_fast.sh.

#' Find the coverage shrink factor that matches two conditions
#'
#' Search for the factor `s` in `interval` that minimises the two-sample
#' Kolmogorov-Smirnov statistic between the coverage of the `target` cells and
#' `s` times the coverage of the `reference` cells, using
#' `stats::optimize()`. The returned `shrink_factor` is
#' `round(optimum, digits) + offset` and is meant to be used as the keep
#' probability in [downsample_fragments()] for the reference samples.
#'
#' @param meta_data A `data.frame` of cell metadata.
#' @param group_col Column holding the condition labels.
#' @param target Label of the condition whose coverage is kept (e.g.
#'   `"ama"`).
#' @param reference Label of the condition whose coverage is scaled (e.g.
#'   `"ctr"`).
#' @param coverage_col Numeric coverage column (e.g. from [calculate_mapd()]).
#' @param interval Search interval for the factor, passed to
#'   `stats::optimize()`.
#' @param digits Rounding digits applied to the optimum.
#' @param offset Number added after rounding; the default `0` returns the
#'   rounded optimum.
#' @param tol Tolerance passed to `stats::optimize()` (its default).
#'
#' @return A list with `optimum` (unrounded minimiser), `objective` (KS
#'   statistic at the optimum) and `shrink_factor`
#'   (`round(optimum, digits) + offset`).
#'
#' @examples
#' set.seed(1)
#' md <- data.frame(state = rep(c("ama", "ctr"), each = 300),
#'                  coverage = c(rlnorm(300, log(80), 0.3), rlnorm(300, log(240), 0.3)))
#' optimize_shrink_factor(md)
#'
#' @export
optimize_shrink_factor <- function(meta_data,
                                   group_col = "state",
                                   target = "ama",
                                   reference = "ctr",
                                   coverage_col = "coverage",
                                   interval = c(0.1, 1),
                                   digits = 2L,
                                   offset = 0,
                                   tol = .Machine$double.eps^0.25) {
  if (!is.data.frame(meta_data)) {
    stop("'meta_data' must be a data.frame.")
  }
  for (col in c(group_col, coverage_col)) {
    if (!col %in% names(meta_data)) {
      stop("Missing column '", col, "' in 'meta_data'.")
    }
  }
  if (!is.numeric(meta_data[[coverage_col]])) {
    stop("'", coverage_col, "' must be numeric.")
  }
  grp <- as.character(meta_data[[group_col]])
  if (!any(grp %in% target, na.rm = TRUE) || !any(grp %in% reference, na.rm = TRUE)) {
    stop("Both 'target' and 'reference' must occur in meta_data$", group_col, ".")
  }
  if (!is.numeric(interval) || length(interval) != 2L || anyNA(interval) ||
      interval[1] >= interval[2]) {
    stop("'interval' must be an increasing numeric vector of length 2.")
  }
  if (!is.numeric(offset) || length(offset) != 1L || is.na(offset)) {
    stop("'offset' must be a single number.")
  }

  cov_target <- meta_data[[coverage_col]][!is.na(grp) & grp == target]
  cov_ref <- meta_data[[coverage_col]][!is.na(grp) & grp == reference]

  objective <- function(shrink_factor) {
    stats::ks.test(cov_target, cov_ref * shrink_factor)$statistic
  }

  res <- stats::optimize(f = objective, interval = interval, maximum = FALSE,
                         tol = tol)

  list(
    optimum = res$minimum,
    objective = unname(res$objective),
    shrink_factor = round(res$minimum, digits) + offset
  )
}

#' Randomly downsample a 10x fragment file
#'
#' Stream a coordinate-sorted, bgzipped fragment file once and keep each
#' fragment record independently with probability `keep_prob` (SplitMix64 RNG
#' seeded with `seed`; header lines starting with `#` are always kept), then
#' write BGZF output and build a tabix BED index. This is an R wrapper around
#' the packaged shell script
#' `system.file("scripts", "downsample_fragments_fast.sh", package = "scRTtools")`,
#' which compiles a small C sampler with `cc` at run time.
#'
#' Besides `output`, the script writes `<output>.tbi` and
#' `<output>.downsample_info.txt` (parameters and kept/total record counts).
#' Results are reproducible for a given `seed` and input.
#'
#' @param input Path to the input `.tsv.gz` fragment file.
#' @param output Path to the output `.tsv.gz` file (must differ from `input`;
#'   its directory must exist and be writable).
#' @param keep_prob Keep probability in `(0, 1]`, e.g.
#'   `optimize_shrink_factor()$shrink_factor`.
#' @param seed Non-negative integer seed of the sampler.
#' @param threads Positive integer; bgzip threads.
#' @param compression_level Integer 0-9; bgzip compression level.
#' @param overwrite Logical; replace an existing output and index.
#' @param bgzip,tabix,cc Paths to the `bgzip`, `tabix` and C compiler
#'   executables. Default: found on `PATH`.
#' @param verbose Logical; show the script's progress output. Default `FALSE`.
#'
#' @return Invisibly, `output`.
#'
#' @export
downsample_fragments <- function(input,
                                 output,
                                 keep_prob,
                                 seed = 20260821,
                                 threads = 16L,
                                 compression_level = 1L,
                                 overwrite = FALSE,
                                 bgzip = Sys.which("bgzip"),
                                 tabix = Sys.which("tabix"),
                                 cc = Sys.which("cc"),
                                 verbose = FALSE) {
  is_scalar_chr <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  is_whole <- function(x, lo) {
    is.numeric(x) && length(x) == 1L && !is.na(x) && x >= lo && x == round(x)
  }

  if (!is_scalar_chr(input) || !file.exists(input)) {
    stop("'input' must be an existing file.")
  }
  if (!is_scalar_chr(output)) {
    stop("'output' must be a single file path.")
  }
  if (normalizePath(input, mustWork = FALSE) ==
      normalizePath(output, mustWork = FALSE)) {
    stop("'input' and 'output' must differ.")
  }
  if (!dir.exists(dirname(output))) {
    stop("Output directory does not exist: ", dirname(output))
  }
  if (!is.numeric(keep_prob) || length(keep_prob) != 1L || is.na(keep_prob) ||
      keep_prob <= 0 || keep_prob > 1) {
    stop("'keep_prob' must be a single number in (0, 1].")
  }
  if (!is_whole(seed, 0)) stop("'seed' must be a non-negative integer.")
  if (!is_whole(threads, 1)) stop("'threads' must be a positive integer.")
  if (!is_whole(compression_level, 0) || compression_level > 9) {
    stop("'compression_level' must be an integer between 0 and 9.")
  }
  tools <- c(bgzip = unname(bgzip), tabix = unname(tabix), cc = unname(cc))
  for (nm in names(tools)) {
    if (!is_scalar_chr(tools[[nm]]) || !file.exists(tools[[nm]])) {
      stop("'", nm, "' executable not found; install it or pass its path via '",
           nm, "'.")
    }
  }

  script <- system.file("scripts", "downsample_fragments_fast.sh",
                        package = "scRTtools", mustWork = TRUE)

  args <- c(
    shQuote(script),
    "-p", format(keep_prob, digits = 15, scientific = FALSE),
    "-s", format(seed, scientific = FALSE),
    "-t", format(threads, scientific = FALSE),
    "-l", format(compression_level, scientific = FALSE),
    if (isTRUE(overwrite)) "-f",
    shQuote(input), shQuote(output)
  )
  env <- c(
    paste0("BGZIP_BIN=", shQuote(tools[["bgzip"]])),
    paste0("TABIX_BIN=", shQuote(tools[["tabix"]])),
    paste0("CC_BIN=", shQuote(tools[["cc"]]))
  )

  out_stream <- if (isTRUE(verbose)) "" else FALSE
  status <- system2("bash", args = args, env = env,
                    stdout = out_stream, stderr = out_stream)
  if (!identical(as.integer(status), 0L)) {
    stop("downsample_fragments_fast.sh failed with exit status ", status,
         if (!isTRUE(verbose)) " (rerun with verbose = TRUE for details)", ".")
  }
  invisible(output)
}
