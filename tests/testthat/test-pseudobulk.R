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
    cnt <- aggregate_pseudobulk_s_g1_counts(f$mats, f$md, names(f$gr), verbose = TRUE),
    "sA"
  )
  expect_silent(
    cnt2 <- aggregate_pseudobulk_s_g1_counts(f$mats, f$md, names(f$gr), verbose = FALSE)
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

  # a single matrix is matched as is; an unnamed list gets no prefix
  m1 <- f$mats$sB
  md1 <- f$md[colnames(m1), ]
  c_one <- aggregate_pseudobulk_s_g1_counts(m1, md1, names(f$gr))
  c_ref <- aggregate_pseudobulk_s_g1_counts(list(sB = m1), f$md, names(f$gr))
  expect_identical(c_one$s_counts, c_ref$s_counts)
  expect_identical(attr(c_one, "selected_cells")$sample, "1")
  expect_error(aggregate_pseudobulk_s_g1_counts(unname(f$mats), f$md, names(f$gr)))

  # input checks
  expect_error(aggregate_pseudobulk_s_g1_counts(f$files, f$md, names(f$gr)))
  expect_error(aggregate_pseudobulk_s_g1_counts(f$mats, f$md[-1, ], names(f$gr), verbose = FALSE))
  expect_error(aggregate_pseudobulk_s_g1_counts(f$mats, f$md, names(f$gr)[-1], verbose = FALSE))
  md3 <- f$md
  md3$is_G1[startsWith(rownames(md3), "sB-")] <- FALSE
  expect_error(
    aggregate_pseudobulk_s_g1_counts(f$mats, md3, names(f$gr), verbose = FALSE),
    "Both QC-passed S cells and G1 cells"
  )
})

test_that("calculate_pseudobulk_rt_g1_ratio z-scores the S/G1 ratio", {
  f <- synth_pseudobulk_files()
  cnt <- aggregate_pseudobulk_s_g1_counts(f$mats, f$md, names(f$gr), verbose = FALSE)
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
  new <- aggregate_pseudobulk_s_g1_counts(f$mats, f$md, names(f$gr), verbose = FALSE)
  old <- suppressMessages(
    orig$aggregate_pseudobulk_s_g1_counts(f$files, f$md, names(f$gr), tempfile(fileext = ".rds"))
  )
  expect_identical(new, old)
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
  cnt <- aggregate_pseudobulk_s_g1_counts(f$mats, f$md, names(f$gr), verbose = FALSE)
  ref <- calculate_pseudobulk_rt_g1_ratio(cnt, f$gr)
  rt <- export_pseudobulk_g1_ratio(f$mats, f$md, f$gr, verbose = FALSE)
  expect_identical(S4Vectors::metadata(rt)$selected_cells, attr(cnt, "selected_cells"))
  S4Vectors::metadata(rt) <- list()
  expect_identical(rt, ref)
  # genome bins must be an object, not an RDS path
  gr_file <- tempfile(fileext = ".rds")
  saveRDS(f$gr, gr_file)
  expect_error(export_pseudobulk_g1_ratio(f$mats, f$md, gr_file, verbose = FALSE))
  rt2 <- export_pseudobulk_g1_ratio(f$mats, f$md, f$gr, verbose = FALSE)
  expect_identical(rt2$score, ref$score)

  skip_if_not_installed("rtracklayer")
  # same glue as the original export_pseudobulk_g1_ratio(): the bigWig holds
  # the score column of calculate_pseudobulk_rt_g1_ratio()
  g <- f$gr
  GenomeInfoDb::seqlengths(g) <- c(chr1 = 8e6, chr2 = 8e6)
  bw <- tempfile(fileext = ".bigwig")
  expect_invisible(export_pseudobulk_g1_ratio(f$mats, f$md, g, file = bw, verbose = FALSE))
  back <- rtracklayer::import(bw)
  expect_identical(length(back), length(ref))
  expect_identical(BiocGenerics::start(back), BiocGenerics::start(ref))
  expect_equal(back$score, ref$score, tolerance = 1e-6)
  expect_error(export_pseudobulk_g1_ratio(f$mats, f$md, f$gr, file = bw, verbose = FALSE))
})
