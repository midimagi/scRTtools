test_that("gex_tempo matches the original for both methods and groupings", {
  orig <- source_orig_gex()
  d <- synth_gex()
  meta <- d$meta
  meta$rep_percent[c(4, 9)] <- NA
  sparse <- methods::as(d$counts, "CsparseMatrix")
  args <- list(list(ncells = 25L), list(nbin_rep = 6L), list(nbin_rep = 40L), list())
  for (a in args) {
    for (method in c("sc2p", "cp10k")) {
      # one cell per group makes SC2P very slow and adds nothing
      if (method == "sc2p" && (!length(a) || identical(a$nbin_rep, 40L))) next
      for (x in list(d$counts, sparse)) {
        call_args <- c(list(x, meta, method = method), a)
        expect_identical(do.call(gex_tempo, call_args),
                         do.call(orig$gex_tempo, call_args))
      }
    }
  }
  expect_identical(
    gex_tempo(as.data.frame(d$counts), meta, ncells = 30, yy_thres = 5),
    orig$gex_tempo(as.data.frame(d$counts), meta, ncells = 30, yy_thres = 5)
  )
})

test_that("gex_tempo handles zero-library cells like the original", {
  orig <- source_orig_gex()
  d <- synth_gex()
  counts <- d$counts
  counts[, c(2, 50, 51)] <- 0L
  for (method in c("sc2p", "cp10k")) {
    expect_identical(gex_tempo(counts, d$meta, ncells = 20, method = method),
                     orig$gex_tempo(counts, d$meta, ncells = 20, method = method))
  }
})

test_that("gex_tempo output structure and input validation", {
  d <- synth_gex()
  res <- gex_tempo(d$counts, d$meta, nbin_rep = 5, method = "cp10k")
  expect_named(res, "df_cp10k_mean")
  expect_identical(dim(res$df_cp10k_mean), c(nrow(d$counts), 5L))
  bins <- attr(res, "rep_bins")
  expect_identical(sum(bins$n_cells), ncol(d$counts))
  expect_identical(attr(res, "grouping")$mode, "nbin_rep")
  # cp10k of a single group is the mean of per-cell CP10K
  one <- gex_tempo(d$counts, d$meta, nbin_rep = 1, method = "cp10k")
  manual <- rowMeans(sweep(d$counts, 2, colSums(d$counts), "/") * 1e4)
  expect_equal(one$df_cp10k_mean[[1]], unname(manual))

  res2 <- gex_tempo(d$counts, d$meta, ncells = 40)
  expect_named(res2, c("df_phase_I", "df_phase_II"))
  expect_identical(ncol(res2$df_phase_I), 3L)

  expect_error(gex_tempo(d$counts, d$meta, ncells = 5, nbin_rep = 5), "mutually exclusive")
  expect_error(gex_tempo(d$counts, d$meta[-1, ], ncells = 5), "missing matrix cell IDs")
  expect_error(gex_tempo(d$counts, d$meta, var_use = "nope"), "absent")
  expect_error(gex_tempo(d$counts + 0.5, d$meta, method = "cp10k"), "integer counts")
  bad <- d$meta
  bad$rep_percent <- bad$rep_percent * 100
  expect_error(gex_tempo(d$counts, bad), "0--1")
})
