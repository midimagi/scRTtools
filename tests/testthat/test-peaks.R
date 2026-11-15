peak_signal <- function() {
  gr <- GenomicRanges::GRanges(
    rep(c("chr2", "chr1"), each = 200),
    IRanges::IRanges(rep(seq(1L, by = 10000L, length.out = 200), 2), width = 10000L)
  )
  y <- c(sin(seq(0, 6 * pi, length.out = 200)), cos(seq(0, 5 * pi, length.out = 200)))
  y[c(30:31)] <- max(y[1:60]) + 0.1                       # 2-bin plateau
  y[c(90, 250:252)] <- NA
  gr$score <- y
  gr[-c(150:160)]                                         # a gap on chr2
}

test_that("find_rt_peaks finds local maxima", {
  gr <- peak_signal()
  pk <- find_rt_peaks(gr, flank_bins = 5, min_score = 0.2, min_delta = 0.02)
  expect_s4_class(pk, "GRanges")
  expect_true(all(c("RT_value", "local_median", "local_delta", "plateau_size") %in%
                    names(S4Vectors::mcols(pk))))
  expect_false(is.unsorted(pk))
  expect_true(all(pk$RT_value >= 0.2))
  expect_true(all(pk$local_delta >= 0.02))
  expect_true(any(pk$plateau_size == 2L))
  none <- find_rt_peaks(gr, min_score = 10)
  expect_length(none, 0L)
  expect_true("plateau_size" %in% names(S4Vectors::mcols(none)))
  expect_error(find_rt_peaks(gr, score_col = "x"), "score_col")
  expect_error(find_rt_peaks(gr, flank_bins = 0), "flank_bins")
})

test_that("find_rt_peaks matches the original script", {
  skip_on_cran()
  orig <- source_orig_file("processing/find_rt_peaks.R")
  gr <- peak_signal()
  x <- synth_pseudobulk()
  rt <- suppressWarnings(calculate_pseudobulk_rt(x$m, x$gr, x$md))
  for (a in list(list(gr), list(gr, flank_bins = 3, min_score = NULL, max_plateau_bins = 1),
                 list(gr, max_gap = 1e6, min_delta = 0.3),
                 list(rt, flank_bins = 5, min_score = 0.2, min_delta = 0.02, max_plateau_bins = 3))) {
    expect_identical(do.call(find_rt_peaks, a), do.call(orig$find_rt_peaks, a))
  }
})

feature_files <- function() {
  d <- tempfile("feat_"); dir.create(d)
  f1 <- file.path(d, "f1.bed")
  writeLines(c("chr1\t0\t1000\ta\t10", "chr1\t500\t1500\tb\t20",
               "chr1\t5000\t6000\tc\t30", "chr2\t0\t100\td\t40"), f1)
  f2 <- file.path(d, "f2.bed")
  writeLines("chr1\t900\t1100\tx\t2", f2)
  reg <- file.path(d, "regions.tsv")
  writeLines(c("chr\tstart\tend\tid", "chr1\t0\t1000\tr1", "chr1\t2000\t3000\tr2",
               "chr1\t4999\t5100\tr3"), reg)
  list(f1 = f1, f2 = f2, reg = reg)
}

test_that("annotate_peak_scores aggregates feature scores", {
  skip_if_not_installed("rtracklayer")
  f <- feature_files()
  res <- suppressMessages(suppressWarnings(
    annotate_peak_scores(f$reg, c(A = f$f1, f$f2))
  ))
  expect_true(all(c("id", "A", "f2") %in% names(res)))
  # r1 = 1..1000 overlaps a (1000 bp, 10) and b (500 bp, 20)
  expect_equal(res$A, c((10 * 1000 + 20 * 500) / 1500, NA, 30))
  expect_equal(res$f2, c(2, NA, NA))
  s <- suppressMessages(annotate_peak_scores(f$reg, c(A = f$f1), method = "sum"))
  expect_equal(s$A[1], 10 * 1000 + 20 * 500)
  m <- suppressMessages(annotate_peak_scores(f$reg, c(A = f$f1), method = "max"))
  expect_equal(m$A[1], 20)
  expect_warning(
    w <- suppressMessages(annotate_peak_scores(f$reg, list(A = f$f1, A = f$f2), window = 100)),
    "No overlaps found for feature: A_1"
  )
  expect_true(all(c("A", "A_1") %in% names(w)))
  expect_true(all(w$width == 100))
  expect_error(annotate_peak_scores(f$reg, c(A = "nope.bed")), "do not exist")
  expect_error(annotate_peak_scores(1, c(A = f$f1)), "regions")
})

test_that("annotate_peak_scores matches the original script", {
  skip_on_cran()
  skip_if_not_installed("rtracklayer")
  orig <- source_orig_file("analysis/annotate_peak_scores.R")
  f <- feature_files()
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 800, 5001), width = 300))
  for (a in list(list(f$reg, c(A = f$f1, f$f2)),
                 list(gr, c(A = f$f1, B = f$f2), method = "sum"),
                 list(gr, list(A = f$f1, A = f$f2), window = 1000, method = "max"))) {
    expect_identical(suppressMessages(suppressWarnings(do.call(annotate_peak_scores, a))),
                     suppressMessages(suppressWarnings(do.call(orig$annotate_peak_scores, a))))
  }
})
