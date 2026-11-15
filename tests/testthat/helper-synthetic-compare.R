# Loader for retreat_figure6_helpers.R and synthetic two-state tables.

source_orig_fig6 <- function() {
  for (pkg in c("dplyr", "tidyr")) testthat::skip_if_not_installed(pkg)
  path <- file.path(
    local_dir("SCRTTOOLS_NOTEBOOKS_DS"),
    "retreat_figure6_helpers.R"
  )
  testthat::skip_if_not(file.exists(path), "retreat_figure6_helpers.R not available")
  env <- new.env(parent = asNamespace("dplyr"))
  # The notebooks attach GenomicRanges before source().
  env$makeGRangesFromDataFrame <- GenomicRanges::makeGRangesFromDataFrame
  env$GRanges <- GenomicRanges::GRanges
  env$IRanges <- IRanges::IRanges
  env$tileGenome <- GenomicRanges::tileGenome
  env$findOverlaps <- GenomicRanges::findOverlaps
  env$queryHits <- S4Vectors::queryHits
  env$subjectHits <- S4Vectors::subjectHits
  env$seqnames <- GenomeInfoDb::seqnames
  env$seqinfo <- GenomeInfoDb::seqinfo
  env[["%in%"]] <- S4Vectors::`%in%`
  env[["seqlevels<-"]] <- GenomeInfoDb::`seqlevels<-`
  sys.source(path, envir = env)
  env
}

# Plain data.frame without tibble class, row names or na.action.
plain_df <- function(x) {
  x <- as.data.frame(x)
  attr(x, "na.action") <- NULL
  rownames(x) <- NULL
  x
}

synth_state_long <- function(n = 60L, seed = 41L) {
  withr::with_seed(seed, {
    id <- sprintf("chr%d:%d-%d", rep(1:3, length.out = n),
                  seq_len(n) * 100000L + 1L, seq_len(n) * 100000L + 50000L)
    df <- data.frame(
      iz_id = rep(id, 2),
      state = factor(rep(c("ctr", "ama"), each = n), levels = c("ctr", "ama")),
      usage = c(stats::rexp(n, 0.2), stats::rexp(n, 0.2)),
      width = stats::runif(2 * n, 1e4, 1e5),
      stringsAsFactors = FALSE
    )
    df$usage[c(3, 10, n + 10, n + 20)] <- NA
    df$width[c(5, n + 7)] <- NA
    # some regions exist in one state only; shuffle the rows
    df <- df[-c(1, 2, n + 30), ]
    df[sample(nrow(df)), ]
  })
}
