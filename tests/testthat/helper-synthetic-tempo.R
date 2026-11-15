# Synthetic inputs and original-code loaders for the tempo tests.

source_orig_tempo <- function(file) {
  path <- file.path(orig_rcode_dir(), "analysis", file)
  testthat::skip_if_not(file.exists(path), paste("original", file, "not available"))
  env <- new.env(parent = globalenv())
  # tempo.R holds one non-UTF-8 comment line, so read it as latin1.
  con <- file(path, encoding = "latin1")
  on.exit(close(con))
  eval(parse(con, keep.source = FALSE), envir = env)
  env
}

synth_tempo <- function(n_bins = 300L, n_cells = 80L, seed = 11L, na_frac = 0.03) {
  withr::with_seed(seed, {
    start <- (seq_len(n_bins) - 1L) * 20000L + 1L
    bin_ids <- paste0("chr1-", start, "-", start + 19999L)
    timing <- 0.5 + 0.35 * sin(seq(0, 4 * pi, length.out = n_bins))
    progress <- sort(stats::runif(n_cells, 0.02, 0.98))
    prob <- stats::plogis(outer(timing, progress, function(t, p) (p - t) / 0.08))
    m <- matrix(ifelse(stats::runif(n_bins * n_cells) < prob, 4, 2), n_bins, n_cells,
                dimnames = list(bin_ids, sprintf("cell%03d", seq_len(n_cells))))
    m[sample(length(m), round(na_frac * length(m)))] <- NA
    meta <- data.frame(rep_percent = colMeans(m == 4, na.rm = TRUE),
                       row.names = colnames(m))
    ord <- order(meta$rep_percent)
    m <- m[, ord]
    meta <- meta[ord, , drop = FALSE]
    r_start <- seq(1L, by = 500000L, length.out = 12L)
    regions <- GenomicRanges::GRanges(
      "chr1", IRanges::IRanges(r_start, width = 100000L))
    names(regions) <- sprintf("region%02d", seq_along(regions))
    list(matrix = m, meta = meta, regions = regions)
  })
}
