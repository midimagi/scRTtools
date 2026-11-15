track_inputs <- function() {
  x <- synth_states()
  x$gr$is_mappable <- TRUE
  x$gr$is_mappable[c(30:31, 150, 200:215)] <- FALSE
  x$gr$is_mappable[100] <- NA
  x$st[50:60, ] <- NA                                   # global gap
  rp <- colSums(x$st == 4, na.rm = TRUE) / colSums(!is.na(x$st))
  x$rp <- rev(rp)                                       # names, other order
  x
}

test_that("define_replication_tracks extracts state-4 segments", {
  x <- track_inputs()
  tr <- define_replication_tracks(x$st, x$gr, x$rp, verbose = FALSE)
  expect_s4_class(tr, "GRangesList")
  expect_identical(names(tr), colnames(x$st))
  expect_length(tr[["cell01"]], 0L)                     # never replicated
  g <- tr[["cell20"]]
  expect_true(all(g$cn == 4))
  expect_true(all(g$cell_id == "cell20"))
  expect_true(all(g$rep_percent == x$rp[["cell20"]]))
  expect_true(all(g$na_fraction <= 0.5))
  expect_true(all(g$flag %in% 0:1))
  # tracks never cross chromosomes and never contain only NA bins
  expect_true(all(BiocGenerics::width(g) == g$num_bins * 20000L))
  expect_identical(
    define_replication_tracks(x$st, x$gr, x$rp, verbose = FALSE,
                              BPPARAM = BiocParallel::MulticoreParam(2)),
    tr
  )
  expect_error(define_replication_tracks(x$st, x$gr, x$rp[-1], verbose = FALSE))
  bad <- x$gr; S4Vectors::mcols(bad)$is_mappable <- NULL
  expect_error(define_replication_tracks(x$st, bad, x$rp, verbose = FALSE), "is_mappable")
})

test_that("track lengths and gaps", {
  x <- track_inputs()
  tr <- define_replication_tracks(x$st, x$gr, x$rp, verbose = FALSE)
  len <- calc_rep_track_lengths(tr)
  expect_identical(nrow(len), sum(lengths(tr)))
  expect_identical(names(len)[1:4], c("chr", "start", "end", "width_bp"))
  expect_true(all(c("cell_id", "rep_percent", "flag") %in% names(len)))
  gap <- calc_rep_track_gaps(tr)
  expect_named(gap, c("cell_id", "rep_percent", "chr", "gap_bp"))
  expect_true(all(gap$gap_bp >= 0))
  n_exp <- sum(vapply(tr, function(g) {
    if (length(g) < 2) return(0L)
    s <- as.character(GenomicRanges::seqnames(BiocGenerics::sort(g)))
    sum(s[-1] == s[-length(s)])
  }, integer(1)))
  expect_identical(nrow(gap), n_exp)
  empty <- GenomicRanges::GRangesList(a = GenomicRanges::GRanges())
  expect_identical(nrow(calc_rep_track_lengths(empty)), 0L)
  expect_identical(nrow(calc_rep_track_gaps(empty)), 0L)
  expect_error(calc_rep_track_gaps(list()), "GRangesList")
})

test_that("tracks, lengths and gaps match the original scripts", {
  skip_on_cran()
  o1 <- source_orig_file("processing/define_replication_tracks.R")
  o2 <- source_orig_file("analysis/calc_rep_track_lengths.R")
  o3 <- source_orig_file("analysis/calc_rep_track_gaps.R")
  x <- track_inputs()
  for (a in list(list(4, 5L, 0.5), list(2, 0L, 0.2), list(4, 2L, 1))) {
    new <- define_replication_tracks(x$st, x$gr, x$rp, a[[1]], a[[2]], a[[3]], verbose = FALSE)
    old <- o1$define_replication_tracks(x$st, x$gr, x$rp, a[[1]], a[[2]], a[[3]],
                                        bpparam = BiocParallel::SerialParam(), verbose = FALSE)
    expect_identical(new, old)
    expect_identical(calc_rep_track_lengths(new),
                     quiet_deprecated(o2$calc_rep_track_lengths(old)))
    expect_identical(calc_rep_track_gaps(new), o3$calc_rep_track_gaps(old))
  }
})
