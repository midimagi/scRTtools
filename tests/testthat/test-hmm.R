test_that("infer_single_cell_states returns states for S-phase cells", {
  x <- synth_g1()
  r <- infer_single_cell_states(x$m[, 10], x$breakpoints)
  expect_named(r, c("filtered_counts", "states"))
  expect_length(r$states, nrow(x$m))
  expect_identical(names(r$states), rownames(x$m))
  expect_true(all(r$states %in% c(2, 4)))
  # CNV gain on the first 200 windows of cell 2 is masked
  r2 <- infer_single_cell_states(x$m[, 2], x$breakpoints)
  expect_gt(mean(is.na(r2$filtered_counts[1:200])), 0.8)
  expect_true(all(is.na(r2$filtered_counts[201:600]) == FALSE))
  # failed initialisation -> all states 2
  r3 <- infer_single_cell_states(x$m[, 10], x$breakpoints,
                                 init_ratio_lower = 10, init_ratio_upper = 20)
  expect_true(all(r3$states == 2))
  # missing counts stay NaN
  v <- as.numeric(x$m[, 12]); v[c(3, 50:60)] <- NA
  r4 <- infer_single_cell_states(v, x$breakpoints)
  expect_true(all(is.na(r4$filtered_counts[c(3, 50:60)])))
  expect_error(infer_single_cell_states(matrix(1:4), 1), "numeric vector")
  expect_error(infer_single_cell_states(1:10, 1, init_ratio_upper = 1), "init_ratio_upper")
})

test_that("internal helpers behave like their MATLAB counterparts", {
  # ties go right; ends filled only with "nearest" (max()/min() of empty
  # sets warn, as in the original)
  fm <- function(...) suppressWarnings(scRTtools:::.fillmissing_nearest(...))
  expect_identical(fm(c(NA, 1, NA, NA, NA, 5, NA), "nearest"), c(1, 1, 1, 5, 5, 5, 5))
  expect_identical(fm(c(NA, 1, NA)), c(NA, 1, NA))
  expect_equal(drop(scRTtools:::.aggregate_counts(c(1, 2, NA, 4, 5), 2)), c(3, 4, 5))
  cp <- scRTtools:::.findchangepts_2cpts(c(rep(0, 20), rep(5, 20), rep(0, 20)))
  expect_identical(cp, c(20L, 40L))
  expect_identical(scRTtools:::.calculate_read_count_ratio(c(1, 1, 2, 2), c(2, 2, 4, 4)), 2)
})

test_that("infer_single_cell_states matches the original script", {
  skip_on_cran()
  orig <- source_orig_hmm()
  x <- synth_g1()
  arg_sets <- list(
    list(),
    list(hmm_aggr_width = 15L, init_ratio_lower = 1.2, final_ratio_lower = 1.4),
    list(cnv_fc_cutoff = 1.3, autocorr_threshold = 0.5, min_switch_prob = 0, maxiter = 5L)
  )
  for (args in arg_sets) {
    for (j in seq_len(ncol(x$m))) {
      a <- c(list(counts = x$m[, j], breakpoints = x$breakpoints), args)
      expect_identical(do.call(infer_single_cell_states, a),
                       do.call(orig$infer_single_cell_states, a))
    }
  }
  v <- as.numeric(x$m[, 12]); v[c(3, 50:60)] <- NA
  expect_identical(infer_single_cell_states(v, x$breakpoints),
                   orig$infer_single_cell_states(v, x$breakpoints))
})
