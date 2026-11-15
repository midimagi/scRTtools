test_that("compute_rep_percent computes rep_percent and filters cells", {
  x <- synth_states()
  res <- compute_rep_percent(x$st, x$gr, x$md, min_rep_percent = 0.01,
                             max_sigma = 0.3, verbose = FALSE)
  md <- res$meta_data
  expect_true(all(c("rep_percent", "sigma_chr", "pass_qc") %in% names(md)))
  expect_equal(md$rep_percent, colSums(x$st == 4, na.rm = TRUE) / colSums(!is.na(x$st)),
               ignore_attr = TRUE)
  expect_false(md["cell01", "pass_qc"])  # rep_percent 0
  expect_false(md["cell10", "pass_qc"])  # whole-chromosome 2 vs 4
  m <- res$genobins_rep_states_filtered
  expect_equal(colnames(m), rownames(md)[md$pass_qc][order(md$rep_percent[md$pass_qc])])
  expect_equal(rownames(m), names(x$gr))
  expect_gt(sum(md$pass_qc), 25)
  # parallel chromosomes give the same result
  expect_identical(
    compute_rep_percent(x$st, x$gr, x$md, min_rep_percent = 0.01, max_sigma = 0.3,
                        verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(2)),
    res
  )
  bad <- x$st
  colnames(bad) <- rev(colnames(bad))
  expect_error(compute_rep_percent(bad, x$gr, x$md, verbose = FALSE), "identical")
})

test_that("compute_rep_percent matches the original script", {
  skip_on_cran()
  orig <- source_orig("compute_rep_percent.R")
  x <- synth_states()
  for (args in list(list(0.01, 0.3, 5, 0.2), list(0.001, 0.4, 3, 0.1))) {
    new <- compute_rep_percent(x$st, x$gr, x$md, args[[1]], args[[2]], args[[3]],
                               args[[4]], verbose = FALSE)
    old <- orig$compute_rep_percent(x$st, x$gr, x$md, args[[1]], args[[2]], args[[3]],
                                    args[[4]], bpparam = BiocParallel::SerialParam(),
                                    verbose = FALSE)
    expect_identical(new, old)
  }
})

test_that("assign_replication_states projects states to genome windows", {
  x <- synth_g1()
  gb <- synth_genome_bins(x$gr)
  ia <- list(hmm_aggr_width = 15L, init_ratio_lower = 1.2, final_ratio_lower = 1.4)
  res <- assign_replication_states(x$m, x$gr, gb, chunk_size = 7L, infer_args = ia,
                                   verbose = FALSE)
  expect_named(res, c("genobins_rep_states", "g1bins_rep_states"))
  g1 <- res$g1bins_rep_states
  expect_identical(dimnames(g1), dimnames(as.matrix(x$m)))
  expect_type(g1, "integer")
  expect_true(all(is.na(g1[x$gr$is_masked_window, ])))
  expect_true(all(g1 %in% c(2L, 4L, NA)))
  st <- infer_single_cell_states(x$m[, 10], x$breakpoints, hmm_aggr_width = 15L,
                                 init_ratio_lower = 1.2, final_ratio_lower = 1.4)$states
  st[x$gr$is_masked_window] <- NA
  expect_equal(g1[, 10], as.integer(st), ignore_attr = TRUE)
  gs <- res$genobins_rep_states
  expect_identical(dim(gs), c(length(gb), ncol(x$m)))
  # windows at chromosome ends and long gaps stay NA; short inner gaps
  # (is_mappable FALSE or NA) are filled when both flanks agree
  expect_true(all(is.na(gs[c(1:2, 500:510, 2399:2400), ])))
  for (i in c(101, 777)) {
    ok <- !is.na(gs[i - 1, ]) & gs[i - 1, ] == gs[i + 1, ]
    expect_true(any(ok))
    expect_identical(gs[i, ok], gs[i - 1, ok])
    expect_true(all(is.na(gs[i, !ok])))
  }
  # chunking and parallelism do not change the result
  expect_identical(
    assign_replication_states(x$m, x$gr, gb, chunk_size = 100L, infer_args = ia,
                              BPPARAM = BiocParallel::MulticoreParam(2), verbose = FALSE),
    res
  )
  expect_error(assign_replication_states(as.matrix(x$m), x$gr, gb, verbose = FALSE),
               "dgCMatrix")
  bad <- x$gr; S4Vectors::mcols(bad)$is_masked_window <- NULL
  expect_error(assign_replication_states(x$m, bad, gb, verbose = FALSE), "is_masked_window")
  expect_error(assign_replication_states(x$m, x$gr, gb[c(1300:1400, 1:1299, 1401:3600)], verbose = FALSE),
               "grouped by chromosome")
})

test_that("assign_replication_states matches the original script", {
  skip_on_cran()
  orig <- source_orig_hmm(with_assign = TRUE)
  x <- synth_g1()
  gb <- synth_genome_bins(x$gr)
  for (ia in list(list(hmm_aggr_width = 15L, init_ratio_lower = 1.2, final_ratio_lower = 1.4),
                  list())) {
    new <- assign_replication_states(x$m, x$gr, gb, gap_max_bp = 60000L, chunk_size = 9L,
                                     infer_args = ia, verbose = FALSE)
    old <- orig$assign_replication_states(x$m, x$gr, gb, gap_max_bp = 60000L,
                                          bpparam = BiocParallel::SerialParam(),
                                          chunk_size = 9L, infer_args = ia, verbose = FALSE)
    expect_identical(new, old)
  }
})
