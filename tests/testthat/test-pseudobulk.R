test_that("calculate_pseudobulk_rt returns a smoothed z-scored signal", {
  x <- synth_pseudobulk()
  rt <- suppressWarnings(calculate_pseudobulk_rt(x$m, x$gr, x$md))
  expect_s4_class(rt, "GRanges")
  expect_identical(names(S4Vectors::mcols(rt)), "score")
  expect_false(anyNA(rt$score))
  expect_true(all(abs(rt$score) <= 5))
  # unmappable gaps shorter than min.gapwidth lie inside a contig and get an
  # interpolated score (original behaviour); long gaps do not
  um <- names(x$gr)[!x$gr$is_mappable]
  expect_true(any(um[um %in% names(x$gr)[c(120:125, 300)]] %in% names(rt)))
  expect_false(any(names(x$gr)[600:640] %in% names(rt)))
  expect_false(is.unsorted(rt))
  rt2 <- calculate_pseudobulk_rt(x$m, x$gr, x$md, smooth_spar = 0.5)
  expect_false(isTRUE(all.equal(rt$score, rt2$score)))
  md <- x$md; md$is_G1 <- NULL
  expect_error(calculate_pseudobulk_rt(x$m, x$gr, md), "is_G1")
  g <- x$gr; S4Vectors::mcols(g)$mappability <- NULL
  expect_error(calculate_pseudobulk_rt(x$m, g, x$md), "mappability")
  expect_error(calculate_pseudobulk_rt(x$m, x$gr, x$md, smooth_spar = NA), "smooth_spar")
})

test_that("calculate_pseudobulk_rt matches the original script", {
  skip_on_cran()
  orig <- source_orig_file("processing/calculate_pseudobulk_rt.R", dplyr = TRUE)
  x <- synth_pseudobulk()
  for (a in list(list(), list(extreme_cutoff = 2, min.gapwidth = 1L, smooth_spar = 0.3),
                 list(min_bins_per_contig = 50, min_non_na_points = 10))) {
    new <- suppressWarnings(do.call(calculate_pseudobulk_rt, c(list(x$m, x$gr, x$md), a)))
    old <- suppressWarnings(do.call(orig$calculate_pseudobulk_rt, c(list(x$m, x$gr, x$md), a)))
    expect_identical(new, old)
  }
})

test_that("aggregate_pseudobulk_s_g1_counts sums per-cell CP10K", {
  f <- synth_pseudobulk_files()
  expect_message(
    cnt <- aggregate_pseudobulk_s_g1_counts(f$files, f$md, names(f$gr), verbose = TRUE),
    "sA"
  )
  expect_silent(
    cnt2 <- aggregate_pseudobulk_s_g1_counts(f$files, f$md, names(f$gr), verbose = FALSE)
  )
  expect_identical(cnt, cnt2)
  expect_identical(cnt$bin_id, names(f$gr))
  cp <- as.matrix(f$m) %*% diag(1e4 / Matrix::colSums(f$m))
  md <- f$md[seq_len(ncol(f$m)), ]
  expect_equal(cnt$s_counts, unname(rowSums(cp[, !md$is_G1 & md$pass_qc])))
  expect_equal(cnt$g1_counts, unname(rowSums(cp[, md$is_G1])))
  sel <- attr(cnt, "selected_cells")
  expect_identical(sel$sample, c("sA", "sB"))
  expect_identical(sum(sel$n_s), sum(!md$is_G1 & md$pass_qc))
  expect_identical(sum(sel$n_g1), sum(md$is_G1))
  expect_identical(attr(cnt, "normalization"), "cp10k_per_cell_v1")

  # cache: reused when the signature matches, rebuilt when it does not
  cache <- tempfile(fileext = ".rds")
  c1 <- aggregate_pseudobulk_s_g1_counts(f$files, f$md, names(f$gr), cache, verbose = FALSE)
  expect_true(file.exists(cache))
  expect_silent(c2 <- aggregate_pseudobulk_s_g1_counts(f$files, f$md, names(f$gr), cache))
  expect_identical(c1, c2)
  md2 <- f$md
  md2$pass_qc[which(md2$is_G1 %in% FALSE & md2$pass_qc)[1:2]] <- FALSE
  expect_message(
    c3 <- aggregate_pseudobulk_s_g1_counts(f$files, md2, names(f$gr), cache, verbose = TRUE),
    "Aggregating"
  )
  expect_false(identical(c1$s_counts, c3$s_counts))
  expect_identical(c1$g1_counts, c3$g1_counts)

  # input checks
  expect_error(aggregate_pseudobulk_s_g1_counts(unname(f$files), f$md, names(f$gr)))
  expect_error(aggregate_pseudobulk_s_g1_counts(f$files, f$md[-1, ], names(f$gr), verbose = FALSE))
  expect_error(aggregate_pseudobulk_s_g1_counts(f$files, f$md, names(f$gr)[-1], verbose = FALSE))
  md3 <- f$md
  md3$is_G1[startsWith(rownames(md3), "sB-")] <- FALSE
  expect_error(
    aggregate_pseudobulk_s_g1_counts(f$files, md3, names(f$gr), verbose = FALSE),
    "Both QC-passed S cells and G1 cells"
  )
})

test_that("calculate_pseudobulk_rt_g1_ratio z-scores the S/G1 ratio", {
  f <- synth_pseudobulk_files()
  cnt <- aggregate_pseudobulk_s_g1_counts(f$files, f$md, names(f$gr), verbose = FALSE)
  cnt$g1_counts[10] <- 0
  cnt$s_counts[11] <- NA
  rt <- calculate_pseudobulk_rt_g1_ratio(cnt, f$gr)
  expect_s4_class(rt, "GRanges")
  expect_identical(names(S4Vectors::mcols(rt)),
                   c("score", "s_g1_ratio", "s_counts", "g1_counts"))
  keep <- setdiff(which(f$gr$is_mappable), 10:11)
  expect_identical(names(rt), names(f$gr)[keep])
  expect_equal(mean(rt$score), 0)
  expect_equal(stats::sd(rt$score), 1)
  r <- (cnt$s_counts[keep] / sum(cnt$s_counts[keep])) /
    (cnt$g1_counts[keep] / sum(cnt$g1_counts[keep]))
  expect_equal(rt$s_g1_ratio, r)
  expect_error(calculate_pseudobulk_rt_g1_ratio(cnt[-1, ], f$gr))
  g <- f$gr
  g$is_mappable <- FALSE
  expect_error(calculate_pseudobulk_rt_g1_ratio(cnt, g), "Insufficient")
  cnt$s_counts <- cnt$g1_counts
  expect_error(calculate_pseudobulk_rt_g1_ratio(cnt, f$gr), "no finite variation")
})

test_that("S/G1 ratio functions match the original script", {
  skip_on_cran()
  orig <- source_orig_file("processing/calculate_pseudobulk_rt_g1_ratio.R")
  f <- synth_pseudobulk_files()
  cache_new <- tempfile(fileext = ".rds")
  cache_old <- tempfile(fileext = ".rds")
  new <- aggregate_pseudobulk_s_g1_counts(f$files, f$md, names(f$gr), cache_new, verbose = FALSE)
  old <- suppressMessages(
    orig$aggregate_pseudobulk_s_g1_counts(f$files, f$md, names(f$gr), cache_old)
  )
  expect_identical(new, old)
  expect_identical(readRDS(cache_new), readRDS(cache_old))
  # each implementation accepts the other's cache
  expect_identical(
    aggregate_pseudobulk_s_g1_counts(f$files, f$md, names(f$gr), cache_old), old
  )
  expect_identical(
    calculate_pseudobulk_rt_g1_ratio(new, f$gr),
    orig$calculate_pseudobulk_rt_g1_ratio(old, f$gr)
  )
  new$g1_counts[c(5, 50)] <- 0
  new$s_counts[70] <- NA
  expect_identical(
    calculate_pseudobulk_rt_g1_ratio(new, f$gr),
    orig$calculate_pseudobulk_rt_g1_ratio(new, f$gr)
  )
})

test_that("export_pseudobulk_g1_ratio returns the ratio and writes a bigWig", {
  f <- synth_pseudobulk_files()
  cnt <- aggregate_pseudobulk_s_g1_counts(f$files, f$md, names(f$gr), verbose = FALSE)
  ref <- calculate_pseudobulk_rt_g1_ratio(cnt, f$gr)
  rt <- export_pseudobulk_g1_ratio(f$files, f$md, f$gr, verbose = FALSE)
  expect_identical(S4Vectors::metadata(rt)$selected_cells, attr(cnt, "selected_cells"))
  S4Vectors::metadata(rt) <- list()
  expect_identical(rt, ref)
  # genome bins given as an RDS path, arguments passed on through ...
  gr_file <- tempfile(fileext = ".rds")
  saveRDS(f$gr, gr_file)
  cache <- tempfile(fileext = ".rds")
  rt2 <- export_pseudobulk_g1_ratio(f$files, f$md, gr_file, cache_file = cache, verbose = FALSE)
  expect_true(file.exists(cache))
  expect_identical(rt2$score, ref$score)

  skip_if_not_installed("rtracklayer")
  # same glue as the original export_pseudobulk_g1_ratio(): the bigWig holds
  # the score column of calculate_pseudobulk_rt_g1_ratio()
  g <- f$gr
  GenomeInfoDb::seqlengths(g) <- c(chr1 = 8e6, chr2 = 8e6)
  bw <- tempfile(fileext = ".bigwig")
  expect_invisible(export_pseudobulk_g1_ratio(f$files, f$md, g, file = bw, verbose = FALSE))
  back <- rtracklayer::import(bw)
  expect_identical(length(back), length(ref))
  expect_identical(BiocGenerics::start(back), BiocGenerics::start(ref))
  expect_equal(back$score, ref$score, tolerance = 1e-6)
  expect_error(export_pseudobulk_g1_ratio(f$files, f$md, f$gr, file = bw, verbose = FALSE))
})
