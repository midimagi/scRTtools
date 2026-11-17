# Small synthetic inputs for the bins / read-count tests.

synth_windows <- function(n_win = 300L, width = 1000L, chrs = c("chr1", "chr2")) {
  grl <- lapply(chrs, function(chr) {
    gr <- GenomicRanges::GRanges(
      chr, IRanges::IRanges(start = seq(1L, by = width, length.out = n_win),
                            width = width)
    )
    names(gr) <- paste(chr, BiocGenerics::start(gr), BiocGenerics::end(gr), sep = "-")
    gr
  })
  names(grl) <- chrs
  GenomicRanges::GRangesList(grl)
}

synth_gaps <- function() {
  GenomicRanges::GRanges(
    c("chr1", "chr2"),
    IRanges::IRanges(c(100001L, 150001L), width = c(20000L, 120000L))
  )
}

# Read positions with a copy-number step on chr1 and a few outlier windows.
synth_reads <- function(grl, seed = 1L) {
  withr::with_seed(seed, {
    lapply(names(grl), function(chr) {
      gr <- grl[[chr]]
      n <- length(gr)
      lambda <- rep(20, n)
      if (chr == "chr1") lambda[seq(n / 2, n)] <- 35
      lambda[sample.int(n, 5)] <- 200
      cnt <- stats::rpois(n, lambda)
      unlist(mapply(function(s, e, k) sort(sample(s:e, k, replace = TRUE)),
                    BiocGenerics::start(gr), BiocGenerics::end(gr), cnt,
                    SIMPLIFY = FALSE), use.names = FALSE)
    })
  }) |> stats::setNames(names(grl))
}

write_reads_dir <- function(reads) {
  d <- tempfile("reads_")
  dir.create(d)
  for (chr in names(reads)) {
    writeLines(as.character(reads[[chr]]),
               file.path(d, paste0("readsChr", sub("chr", "", chr), ".txt")))
  }
  d
}

# Per-base GC-bin labels -> wins_gc GRangesList + coords (in memory and as RDS paths).
synth_gc <- function(grl, seed = 2L) {
  gc_cols <- sprintf("GC_bin_%03d", 0:400)
  coords_list <- list()
  coords_files <- character()
  wins_gc <- list()
  withr::with_seed(seed, {
    for (chr in names(grl)) {
      gr <- grl[[chr]]
      pos <- seq.int(min(BiocGenerics::start(gr)), max(BiocGenerics::end(gr)))
      lab <- sample(150:250, length(pos), replace = TRUE)
      coords <- lapply(1:401, function(j) as.numeric(pos[lab == j]))
      f <- tempfile(paste0("coords_", chr, "_"), fileext = ".rds")
      saveRDS(coords, f)
      coords_files[chr] <- f
      coords_list[[chr]] <- coords
      win_idx <- findInterval(pos, BiocGenerics::start(gr))
      tab <- matrix(0, length(gr), 401, dimnames = list(NULL, gc_cols))
      tt <- table(factor(win_idx, levels = seq_along(gr)),
                  factor(lab, levels = 1:401))
      tab[] <- as.numeric(tt)
      g <- gr
      S4Vectors::mcols(g) <- S4Vectors::DataFrame(tab, check.names = FALSE)
      S4Vectors::mcols(g)$row_sum <- rowSums(tab)
      wins_gc[[chr]] <- g
    }
  })
  list(wins_gc = GenomicRanges::GRangesList(wins_gc), coords_list = coords_list,
       coords_files = coords_files)
}

# Path to the original calculate_rt.R / make_fixed_bins.R for regression tests.
orig_processing_dir <- function() {
  local_dir("SCRTTOOLS_ORIG_PROCESSING")
}

source_orig <- function(file) {
  path <- file.path(orig_processing_dir(), file)
  testthat::skip_if_not(file.exists(path), paste("original", file, "not available"))
  # The notebooks attach GenomicRanges before source(); the original scripts
  # rely on that for unqualified S4 generics.
  env <- new.env()
  env$seqnames <- GenomeInfoDb::seqnames
  env$start <- BiocGenerics::start
  env$end <- BiocGenerics::end
  env$sort <- BiocGenerics::sort
  env$unlist <- BiocGenerics::unlist
  env$keepSeqlevels <- GenomeInfoDb::keepSeqlevels
  env$seqlevelsInUse <- GenomeInfoDb::seqlevelsInUse
  env$seqlevels <- GenomeInfoDb::seqlevels
  sys.source(path, envir = env)
  env
}

# Generic loader for an original script given its path relative to rcode/,
# with the unqualified generics (and optionally dplyr verbs) it relies on.
orig_rcode_dir <- function() dirname(orig_processing_dir())

source_orig_file <- function(rel_path, dplyr = FALSE) {
  path <- file.path(orig_rcode_dir(), rel_path)
  testthat::skip_if_not(file.exists(path), paste("original", rel_path, "not available"))
  parent <- if (dplyr) {
    testthat::skip_if_not_installed("dplyr")
    asNamespace("dplyr")
  } else {
    globalenv()
  }
  env <- new.env(parent = parent)
  env$mcols <- S4Vectors::mcols
  env$sort <- BiocGenerics::sort
  env$findOverlaps <- GenomicRanges::findOverlaps
  env$queryHits <- S4Vectors::queryHits
  env$subjectHits <- S4Vectors::subjectHits
  sys.source(path, envir = env)
  env
}

# The original scripts (and notebook code copied into tests) convert
# DataFrame / GRanges with as.data.frame(), which `R CMD check --as-cran`
# reports as deprecated direct calls of as.data.frame.<type>() (a false
# positive from S4Vectors under R 4.3). Silence only that warning class.
quiet_deprecated <- function(expr) suppressWarnings(expr, classes = "deprecatedWarning")
