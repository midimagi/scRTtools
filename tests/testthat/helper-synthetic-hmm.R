# Synthetic G1-window counts for the HMM / replication-state tests.
synth_g1 <- function(n_chr = 3L, n_win = 600L, n_cells = 30L, depth = 12, seed = 11L) {
  withr::with_seed(seed, {
    chrs <- paste0("chr", seq_len(n_chr))
    gr <- GenomicRanges::GRanges(
      rep(chrs, each = n_win),
      IRanges::IRanges(rep(seq(1L, by = 50000L, length.out = n_win), n_chr), width = 50000L)
    )
    names(gr) <- paste(GenomicRanges::seqnames(gr), BiocGenerics::start(gr),
                       BiocGenerics::end(gr), sep = "-")
    # replication-timing domains of 40 windows
    rt <- unlist(lapply(chrs, function(i) {
      rep(stats::runif(ceiling(n_win / 40)), each = 40)[seq_len(n_win)]
    }))
    p <- seq(0.1, 0.9, length.out = n_cells)
    m <- vapply(seq_len(n_cells), function(j) {
      stats::rpois(length(rt), depth * ifelse(rt < p[j], 2, 1))
    }, numeric(length(gr)))
    m[, 1] <- stats::rpois(nrow(m), depth)                                  # G1-like
    m[1:200, 2] <- stats::rpois(200, depth * 3)                             # CNV gain
    m[, 3] <- stats::rpois(nrow(m), depth * (1 + 3 * (seq_len(nrow(m)) %% 300 < 150)))
  })
  colnames(m) <- sprintf("cell%02d", seq_len(n_cells))
  rownames(m) <- names(gr)
  gr$is_masked_window <- FALSE
  gr$is_masked_window[c(5, 700, 1500)] <- TRUE
  list(gr = gr, m = Matrix::Matrix(m, sparse = TRUE),
       breakpoints = c(1, cumsum(rep(n_win, n_chr - 1L)) + 1))
}

# Genome windows (25 kb) covering the G1 windows, with unmappable runs.
synth_genome_bins <- function(g1_gr) {
  chrs <- unique(as.character(GenomicRanges::seqnames(g1_gr)))
  gr <- unlist(GenomicRanges::GRangesList(lapply(chrs, function(chr) {
    e <- max(BiocGenerics::end(g1_gr[GenomicRanges::seqnames(g1_gr) == chr]))
    GenomicRanges::GRanges(chr, IRanges::IRanges(seq(1L, e, by = 25000L), width = 25000L))
  })), use.names = FALSE)
  names(gr) <- paste(GenomicRanges::seqnames(gr), BiocGenerics::start(gr),
                     BiocGenerics::end(gr), sep = "-")
  gr$is_mappable <- TRUE
  gr$is_mappable[c(1:2, 101, 300:301, 500:510, 1300, 2399:2400)] <- FALSE
  gr$is_mappable[777] <- NA
  gr
}

# The original infer_single_cell_states.R (+ assign_replication_states.R),
# with its Rcpp::sourceCpp() calls (absolute paths that are not readable here)
# replaced by the package's compiled kernels, which are byte-identical copies of
# rcode/src/*.cpp.
source_orig_hmm <- function(with_assign = FALSE) {
  dir <- orig_processing_dir()
  path <- file.path(dir, "infer_single_cell_states.R")
  testthat::skip_if_not(file.exists(path), "original infer_single_cell_states.R not available")
  env <- new.env()
  for (f in c("single_cell_baum_welch_cpp", "compute_baum_welch_expectations",
              "viterbi_decode_cpp", "findchangepts_cpp_max2")) {
    assign(f, utils::getFromNamespace(f, "scRTtools"), envir = env)
  }
  code <- readLines(path)
  code <- code[!grepl("^Rcpp::sourceCpp", code)]
  eval(parse(text = code), envir = env)
  if (with_assign) {
    p2 <- file.path(dir, "assign_replication_states.R")
    testthat::skip_if_not(file.exists(p2), "original assign_replication_states.R not available")
    sys.source(p2, envir = env)
  }
  env
}
