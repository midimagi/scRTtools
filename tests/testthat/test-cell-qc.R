test_that("calculate_mapd adds coverage and MAPD columns", {
  x <- synth_tiles()
  out <- calculate_mapd(x$m, x$gr, x$md, aggr_bins = 10, verbose = FALSE)
  expect_true(all(c("coverage", "mapd", "scaled_mapd", "mapd_error") %in% names(out)))
  expect_equal(rownames(out), rownames(x$md))
  expect_true(all(out$coverage > 0))
  expect_equal(out$scaled_mapd, round(out$mapd / sqrt(out$coverage), 3))
  # noisy cells (even index) have a higher scaled MAPD
  noisy <- seq_len(nrow(out)) %% 2 == 0
  expect_gt(median(out$scaled_mapd[noisy]), median(out$scaled_mapd[!noisy]))
  # base matrix input gives the same result
  expect_identical(
    calculate_mapd(as.matrix(x$m), x$gr, x$md, aggr_bins = 10, verbose = FALSE),
    out
  )
  # number of workers does not change the result
  expect_identical(
    calculate_mapd(x$m, x$gr, x$md, aggr_bins = 10, verbose = FALSE,
                   BPPARAM = BiocParallel::MulticoreParam(2)),
    out
  )
  bad <- x$m
  rownames(bad) <- rev(rownames(bad))
  expect_error(calculate_mapd(bad, x$gr, x$md, verbose = FALSE), "row names")
})

test_that("find_g1_cells recovers the low-MAPD line", {
  x <- synth_mapd()
  res <- find_g1_cells(x$md, seed = 1, verbose = FALSE)
  expect_type(res$is_G1, "logical")
  expect_gt(mean(res$is_G1 == x$truth), 0.95)
  expect_identical(find_g1_cells(x$md, seed = 1, verbose = FALSE), res)
  res_s <- find_g1_cells(x$md, by_sample = TRUE, seed = 1, verbose = FALSE)
  expect_gt(mean(res_s$is_G1 == x$truth), 0.95)

  md_na <- x$md
  md_na$coverage[1] <- NA
  expect_equal(nrow(find_g1_cells(md_na, seed = 1, verbose = FALSE)), nrow(x$md) - 1L)
  expect_warning(find_g1_cells(x$md[1:3, ], verbose = FALSE), "too few|Too few")
  expect_error(find_g1_cells(x$md[, -1]), "cell_id")
  expect_message(find_g1_cells(x$md, seed = 1, verbose = TRUE), "current iteration")
})

test_that("cell QC functions match the original scripts", {
  skip_on_cran()
  orig_mapd <- source_orig("calculate_mapd.R")
  x <- synth_tiles()
  for (ab in c(10, 50)) {
    expect_identical(
      calculate_mapd(x$m, x$gr, x$md, aggr_bins = ab, verbose = FALSE),
      suppressMessages(orig_mapd$calculate_mapd(x$m, x$gr, x$md, aggr_bins = ab,
                                                bpparam = BiocParallel::SerialParam()))
    )
  }

  orig_g1 <- source_orig("find_g1_cells.R")
  y <- synth_mapd()
  for (bs in c(FALSE, TRUE)) {
    set.seed(11)
    old <- suppressMessages(orig_g1$find_g1_cells(y$md, by_sample = bs))
    expect_identical(find_g1_cells(y$md, by_sample = bs, seed = 11, verbose = FALSE),
                     old)
  }
})
