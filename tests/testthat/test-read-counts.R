test_that("count_reads counts reads per window and derives raw_cn", {
  d <- tempfile(); dir.create(d)
  writeLines(as.character(c(5, 15, 15, 25, 35, 40, 99)),
             file.path(d, "readsChr1.txt"))
  wins <- GenomicRanges::GRangesList(
    chr1 = GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 11, 21, 31), width = 10))
  )
  res <- count_reads(wins, d, verbose = FALSE)
  mc <- S4Vectors::mcols(res$wins_count$chr1)
  expect_equal(mc$raw_counts, c(1L, 2L, 1L, 2L))  # 99 lies outside
  expect_equal(mc$raw_cn, c(1, 2, 1, 2) / 1.5 * 2)
  expect_equal(res$reads_list$chr1, c(5, 15, 15, 25, 35, 40, 99))

  res2 <- count_reads(wins, d, bad_coords = list(chr1 = 15), raw_cn_th = 2,
                      verbose = FALSE)
  mc2 <- S4Vectors::mcols(res2$wins_count$chr1)
  expect_equal(mc2$raw_counts, c(1L, 0L, 1L, 2L))
  expect_true(is.na(mc2$raw_cn[2]))           # zero -> NA
  expect_true(is.na(mc2$raw_cn[4]))           # 2/1*2 = 4 > 2 -> NA
  expect_error(count_reads(wins, file.path(d, "nope")), "does not exist")
  expect_message(count_reads(wins, d, verbose = TRUE), "processing chr1")
})

test_that("filter_windows keeps a clean step profile and drops outliers", {
  grl <- synth_windows()
  reads <- synth_reads(grl)
  cnt <- count_reads(grl, write_reads_dir(reads), verbose = FALSE)$wins_count
  out <- filter_windows(cnt, synth_gaps(), "raw_cn", std_th = 2)
  expect_s4_class(out, "GRangesList")
  expect_named(out, names(grl))
  expect_true("raw_cn_seg_score" %in% names(S4Vectors::mcols(out$chr1)))
  expect_true(all(!is.na(S4Vectors::mcols(out$chr1)$raw_cn_seg_score)))
  # no retained window overlaps a gap
  expect_false(any(IRanges::overlapsAny(BiocGenerics::unlist(out), synth_gaps())))
  expect_lt(length(BiocGenerics::unlist(out)), length(BiocGenerics::unlist(cnt)))
  expect_error(filter_windows(cnt, synth_gaps(), "nope"), "nope")
})

test_that("calculate_OE_ratio returns expected counts and GC factors", {
  grl <- synth_windows()
  reads <- synth_reads(grl)
  cnt <- count_reads(grl, write_reads_dir(reads), verbose = FALSE)$wins_count
  filt <- filter_windows(cnt, synth_gaps(), "raw_cn", std_th = 2)
  gc <- synth_gc(grl)
  res <- calculate_OE_ratio(filt, gc$wins_gc, reads, gc$coords_list,
                            gc_bins_use = 81:321, verbose = FALSE)
  expect_named(res, c("wins_ratio", "gc_factors_per_bin"))
  expect_length(res$gc_factors_per_bin, 401L)
  expect_true(all(is.na(res$gc_factors_per_bin[-(150:250)])))
  mc <- S4Vectors::mcols(res$wins_ratio$chr1)
  expect_true(all(c("expected_counts", "oe_ratio") %in% names(mc)))
  expect_equal(length(mc$oe_ratio), length(filt$chr1))
  expect_gt(stats::median(mc$oe_ratio, na.rm = TRUE), 0)
  # a loader function gives the same result as the in-memory list
  res_fun <- calculate_OE_ratio(filt, gc$wins_gc, reads,
                                function(chr) readRDS(gc$coords_files[[chr]]),
                                gc_bins_use = 81:321)
  expect_identical(res_fun, res)
  expect_error(calculate_OE_ratio(filt, gc$wins_gc, reads, gc$coords_files),
               "named list or a function")
})

test_that("smooth_windows adds rolling-mean, spline and z-score columns", {
  grl <- synth_windows()
  reads <- synth_reads(grl)
  cnt <- count_reads(grl, write_reads_dir(reads), verbose = FALSE)$wins_count
  out <- smooth_windows(cnt, "raw_cn", synth_gaps(), smooth_arg = list(spar = 0.1))
  cols <- paste0("raw_cn", c("_ave", "_smoothed", "_ave_z", "_smoothed_z"))
  expect_true(all(cols %in% names(S4Vectors::mcols(out$chr1))))
  z <- unlist(lapply(out, function(g) S4Vectors::mcols(g)$raw_cn_smoothed_z))
  expect_equal(mean(z, na.rm = TRUE), 0, tolerance = 1e-12)
  expect_equal(stats::sd(z, na.rm = TRUE), 1, tolerance = 1e-12)
  # windows inside the 120 kb chr2 gap are not smoothed
  in_gap <- IRanges::overlapsAny(out$chr2, synth_gaps())
  expect_true(all(is.na(S4Vectors::mcols(out$chr2)$raw_cn_smoothed[in_gap])))
})

test_that("export_mcol writes a BigWig without NA scores", {
  skip_if_not_installed("rtracklayer")
  si <- GenomeInfoDb::Seqinfo("chr1", 10000L)
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1001, 1), width = 1000),
                               val = c(2.5, 1.5))
  gr <- c(gr, GenomicRanges::GRanges("chr1", IRanges::IRanges(3001, width = 1000),
                                     val = NA_real_))
  f <- tempfile(fileext = ".bw")
  expect_identical(export_mcol(GenomicRanges::GRangesList(chr1 = gr), "val", f, si), f)
  bw <- rtracklayer::import(f)
  expect_equal(bw$score, c(1.5, 2.5))
  expect_error(export_mcol(GenomicRanges::GRangesList(chr1 = gr), "val",
                           file.path(tempdir(), "no", "x.bw"), si), "directory")
})

test_that("read-count functions match the original calculate_rt.R", {
  skip_on_cran()
  orig <- source_orig("calculate_rt.R")
  grl <- synth_windows()
  reads <- synth_reads(grl)
  d <- write_reads_dir(reads)
  bad <- list(chr1 = reads$chr1[1:20], chr2 = reads$chr2[1:20])

  new_cnt <- count_reads(grl, d, bad, verbose = FALSE)
  old_cnt <- suppressMessages(orig$count_reads(grl, d, bad))
  expect_identical(new_cnt, old_cnt)
  cnt <- new_cnt$wins_count

  for (std_th in c(0.35, 0.5, 2)) {
    expect_identical(filter_windows(cnt, synth_gaps(), "raw_cn", std_th = std_th),
                     orig$filter_windows(cnt, synth_gaps(), "raw_cn", std_th = std_th))
  }
  filt <- filter_windows(cnt, synth_gaps(), "raw_cn", std_th = 2)

  gc <- synth_gc(grl)
  expect_identical(
    calculate_OE_ratio(filt, gc$wins_gc, new_cnt$reads_list, gc$coords_list,
                       verbose = FALSE),
    suppressMessages(orig$calculate_OE_ratio(filt, gc$wins_gc, new_cnt$reads_list,
                                             gc$coords_files))
  )

  for (sa in list(list(spar = 0.1), list(spar = 0.2), list(df = 10))) {
    expect_identical(
      smooth_windows(cnt, "raw_cn", synth_gaps(), smooth_arg = sa),
      orig$smooth_windows(cnt, "raw_cn", synth_gaps(), smooth_arg = sa)
    )
  }

  skip_if_not_installed("rtracklayer")
  si <- GenomeInfoDb::Seqinfo(c("chr1", "chr2"), c(400000L, 400000L))
  f_new <- tempfile(fileext = ".bw")
  f_old <- tempfile(fileext = ".bw")
  export_mcol(cnt, "raw_cn", f_new, si)
  orig$export_mcol(cnt, "raw_cn", f_old, si)
  expect_identical(rtracklayer::import(f_new), rtracklayer::import(f_old))
})
