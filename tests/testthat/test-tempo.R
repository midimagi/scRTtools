test_that("rep_tempo matches the original for nbin_rep and ncells grouping", {
  orig <- source_orig_tempo("tempo.R")
  d <- synth_tempo()
  expect_identical(
    rep_tempo(d$matrix, d$regions, meta_data = d$meta, nbin_rep = 20),
    orig$rep_tempo(d$matrix, d$regions, meta_data = d$meta, nbin_rep = 20)
  )
  expect_identical(
    rep_tempo(d$matrix, d$regions, ncells = 5L),
    orig$rep_tempo(d$matrix, d$regions, ncells = 5L)
  )
  expect_identical(
    rep_tempo(d$matrix, d$regions, meta_data = d$meta, ncells = 8L),
    orig$rep_tempo(d$matrix, d$regions, meta_data = d$meta, ncells = 8L)
  )
})

test_that("rep_tempo state codes are configurable", {
  d <- synth_tempo()
  ref <- rep_tempo(d$matrix, d$regions, meta_data = d$meta, nbin_rep = 20)
  binary <- (d$matrix == 4) + 0
  expect_identical(
    rep_tempo(binary, d$regions, meta_data = d$meta, nbin_rep = 20,
              unrep_state = 0, rep_state = 1), ref)
  swapped <- ifelse(d$matrix == 4, 0, 1)
  dimnames(swapped) <- dimnames(d$matrix)
  expect_identical(
    rep_tempo(swapped, d$regions, meta_data = d$meta, nbin_rep = 20,
              unrep_state = 1, rep_state = 0), ref)
  expect_error(rep_tempo(d$matrix, d$regions, unrep_state = 2, rep_state = 2),
               "must differ")
  expect_error(rep_tempo(d$matrix, d$regions, meta_data = d$meta, nbin_rep = 20,
                         unrep_state = 0, rep_state = 1))
  expect_error(rep_tempo(d$matrix, d$regions, ncells = 2, nbin_rep = 5,
                         meta_data = d$meta), "mutually exclusive")
})

test_that("subset_rep_tempo keeps and subsets the attributes", {
  d <- synth_tempo()
  x <- rep_tempo(d$matrix, d$regions, meta_data = d$meta, nbin_rep = 20)
  z <- subset_rep_tempo(x, c(2L, 5L, 9L))
  expect_identical(rownames(z), rownames(x)[c(2, 5, 9)])
  expect_identical(attr(z, "rep_bins"), attr(x, "rep_bins"))
  expect_identical(attr(z, "n_valid_cells"),
                   attr(x, "n_valid_cells")[c(2, 5, 9), , drop = FALSE])
  expect_identical(attr(z, "range_info")$range_id, rownames(z))
  expect_identical(subset_rep_tempo(x, rownames(x)[c(2, 5, 9)]), z)
  expect_error(subset_rep_tempo(x, "nope"), "absent")
  expect_identical(
    fit_rep_tempo(z, verbose = FALSE)$t50,
    fit_rep_tempo(x, verbose = FALSE)$t50[c(2, 5, 9)])
})

test_that("fit_rep_tempo matches the original for every model type", {
  skip_if_not_installed("mgcv")
  orig <- source_orig_tempo("tempo_models1.R")
  d <- synth_tempo()
  x <- rep_tempo(d$matrix, d$regions, meta_data = d$meta, nbin_rep = 20)
  for (model in c("sslogis", "gam_low", "gam_fast")) {
    expect_identical(
      fit_rep_tempo(x, model_type = model, verbose = FALSE),
      orig$fit_rep_tempo(x, model_type = model, n_boot = 0L, verbose = FALSE),
      info = model
    )
  }
  new <- fit_rep_tempo(x, verbose = FALSE)
  expect_true(all(c("t10", "t50", "t90", "window_10_90", "max_slope",
                    "auc_normalized") %in% names(new)))
  expect_gt(sum(new$fit_status == "ok"), 0)
  expect_identical(
    fit_rep_tempo(x, model_type = "gam_fast", gam_fixed_df = TRUE, gam_k = 6L,
                  gam_method = "REML", anchors = FALSE, curve_ids = rownames(x)[1:2],
                  verbose = FALSE),
    orig$fit_rep_tempo(x, model_type = "gam_fast", gam_fixed_df = TRUE, gam_k = 6L,
                       gam_method = "REML", anchors = FALSE, n_boot = 0L,
                       curve_ids = rownames(x)[1:2], verbose = FALSE)
  )
  plain <- as.matrix(x)
  attributes(plain) <- list(dim = dim(plain), dimnames = dimnames(plain))
  rp <- attr(x, "rep_bins")$rep_percent
  expect_identical(
    fit_rep_tempo(plain, rep_percent = rp, anchors = FALSE, verbose = FALSE),
    orig$fit_rep_tempo(plain, rep_percent = rp, anchors = FALSE, n_boot = 0L,
                       verbose = FALSE))
})

test_that("fit_rep_tempo bootstrap matches the original and restores the RNG", {
  skip_if_not_installed("mgcv")
  orig <- source_orig_tempo("tempo_models1.R")
  d <- synth_tempo()
  x <- subset_rep_tempo(
    rep_tempo(d$matrix, d$regions, meta_data = d$meta, nbin_rep = 20), 1:4)
  set.seed(99)
  before <- .Random.seed
  for (model in c("sslogis", "gam_fast")) {
    expect_identical(
      fit_rep_tempo(x, model_type = model, n_boot = 25L, seed = 7L, verbose = FALSE),
      orig$fit_rep_tempo(x, model_type = model, n_boot = 25L, seed = 7L,
                         verbose = FALSE),
      info = model
    )
  }
  expect_identical(.Random.seed, before)
})

test_that("fit_rep_tempo gives the same result on a parallel backend", {
  skip_on_os("windows")
  d <- synth_tempo()
  x <- subset_rep_tempo(
    rep_tempo(d$matrix, d$regions, meta_data = d$meta, nbin_rep = 20), 1:4)
  serial <- fit_rep_tempo(x, n_boot = 25L, seed = 7L, verbose = FALSE)
  par <- fit_rep_tempo(x, n_boot = 25L, seed = 7L, verbose = FALSE,
                       BPPARAM = BiocParallel::MulticoreParam(workers = 2L))
  expect_identical(attr(par, "execution")$workers, 2L)
  expect_identical(attr(serial, "execution")$workers, 1L)
  attr(par, "execution") <- attr(serial, "execution") <- NULL
  expect_equal(par, serial, tolerance = 1e-12)
  expect_error(fit_rep_tempo(x, BPPARAM = "multicore"), "BiocParallelParam")
})
