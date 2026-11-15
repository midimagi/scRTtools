downstream_inputs <- function() {
  f <- iz_analysis_fixture()
  d <- prepare_iz_data(f$diz$tracks[lengths(f$diz$tracks) > 0], f$fr)
  cells <- d$cell_metadata$cell_id
  prog <- d$cell_metadata$replication_progress
  us <- call_iz_usage_from_states(f$x$st, cells, d$iz_metadata)
  uh <- compute_usage_heterogeneity(us, prog)
  tm <- fit_timing_models(us, prog, uh$probability)
  tf <- summarize_track_features(d$track_intervals)$summary
  sh <- summarize_track_shape(d$track_intervals, d$iz_metadata, window_bp = 200000L)
  s <- assemble_iz_summary(d$iz_metadata, uh$summary, tm, tf, sh)
  withr::with_seed(3, {
    s$SNS <- stats::rexp(nrow(s)) * (1 + s$usage_heterogeneity)
    s$ATAC <- stats::rnorm(nrow(s))
  })
  list(f = f, d = d, us = us, prog = prog, cells = cells, s = s,
       metrics = c("usage_heterogeneity", "position_heterogeneity", "fork_heterogeneity",
                   "track_shape_heterogeneity"))
}

source_orig_downstream <- function() {
  env <- source_orig_heterogeneity()
  path <- file.path(notebook_dir(), "iz_heterogeneity_phase4_7.R")
  testthat::skip_if_not(file.exists(path), "iz_heterogeneity_phase4_7.R not available")
  sys.source(path, envir = env, toplevel.env = env)
  env
}

test_that("IZ heterogeneity downstream functions run and match the original", {
  skip_on_cran()
  skip_if_not_installed("cluster")
  orig <- source_orig_downstream()
  x <- downstream_inputs()
  df <- function(v) { v <- as.data.frame(v); rownames(v) <- NULL; v }

  # clustering
  cl <- cluster_iz_heterogeneity(x$s, x$metrics, k_range = 2:4, nstart = 5L, seed = 7L)
  old <- orig$cluster_iz_heterogeneity(x$s, x$metrics, k_range = 2:4, nstart = 5L, seed = 7L)
  expect_true(cl$k %in% 2:4)
  expect_identical(cl$assignments, old$assignments)
  expect_identical(cl$diagnostics, old$diagnostics)
  expect_identical(cl$profiles, df(old$profiles))
  expect_identical(cl$centers, old$centers)
  expect_identical(cluster_iz_heterogeneity(x$s, x$metrics, 2:4, k = 3, nstart = 5L)$k, 3)

  # nearby associations
  for (m in c("bin", "spline")) {
    na <- compute_nearby_iz_associations(x$us, x$prog, x$d$iz_metadata, max_distance = 3e6,
                                         correction = m, min_complete = 20L)
    expect_identical(na, df(orig$compute_nearby_iz_associations(
      x$us, x$prog, x$d$iz_metadata, max_distance = 3e6, correction = m, min_complete = 20L)))
  }
  expect_true(all(na$distance > 0 & na$distance <= 3e6))
  expect_true(all(na$association_type %in% c("co_usage", "competition", "neutral")))

  # cell / IZ selection and tensor factorisation
  sc <- sample_cells_by_progress(x$d$cell_metadata, n_cells = 40L, n_bins = 5L, seed = 2L)
  expect_identical(sc, orig$sample_cells_by_progress(x$d$cell_metadata, 40L, 5L, 2L))
  expect_lte(length(sc), 40L)
  expect_identical(sample_cells_by_progress(x$d$cell_metadata, 1e4), x$cells)
  iz_sel <- select_iz_for_tensor(x$s, cl$assignments, x$metrics, n_iz = 8L)
  expect_identical(iz_sel, orig$select_iz_for_tensor(x$s, old$assignments, x$metrics, n_iz = 8L))
  expect_lte(length(iz_sel), 8L)
  tens <- build_replication_tensor(x$d$track_intervals, sc, x$d$iz_metadata, iz_sel,
                                   window_bp = 200000L)
  res <- correct_replication_tensor(tens, x$prog[match(sc, x$cells)], "bin")$residual
  y <- split_signed_residual(res)
  expect_identical(y, orig$split_signed_residual(res))
  expect_true(all(y >= 0))
  expect_equal(y[, , 1:21] - y[, , 22:42], res, ignore_attr = TRUE)
  fit <- fit_nonnegative_cp(y, rank = 2L, seed = 4L, max_iter = 20L)
  expect_identical(fit, orig$fit_nonnegative_cp(y, rank = 2L, seed = 4L, max_iter = 20L))
  expect_equal(colSums(fit$B^2), c(1, 1))
  expect_true(fit$explained >= 0 && fit$explained <= 1)
  expect_error(fit_nonnegative_cp(res, 2L), "non-negative")
  rk <- select_cp_rank(y, ranks = 2:3, n_init = 2L, seed = 5L, max_iter = 10L)
  old_rk <- orig$select_cp_rank(y, ranks = 2:3, n_init = 2L, seed = 5L, max_iter = 10L)
  expect_identical(rk$summary, df(old_rk$summary))
  expect_identical(rk$chosen_rank, old_rk$chosen_rank)
  expect_identical(rk$best_fits, old_rk$best_fits)
  rel <- as.numeric(dimnames(tens)[[3]])
  ft <- cp_factor_tables(fit, sc, iz_sel, rel, x$d$cell_metadata, x$s)
  old_ft <- orig$cp_factor_tables(fit, sc, iz_sel, rel, x$d$cell_metadata, x$s)
  for (nm in names(old_ft)) expect_identical(df(ft[[nm]]), df(old_ft[[nm]]), info = nm)
  expect_identical(nrow(ft$track), 21L)

  # epigenomic integration
  ce <- correlate_epigenomic_features(x$s, x$metrics, c("SNS", "ATAC"))
  expect_identical(ce, df(orig$correlate_epigenomic_features(x$s, x$metrics, c("SNS", "ATAC"))))
  expect_identical(nrow(ce), length(x$metrics) * 2L * 2L)
  em <- fit_epigenomic_models(x$s, x$metrics[1:2], c("SNS", "ATAC"))
  expect_identical(em, df(orig$fit_epigenomic_models(x$s, x$metrics[1:2], c("SNS", "ATAC"))))
  ec <- suppressWarnings(fit_epigenomic_classifier(x$s, x$metrics, c("SNS", "ATAC")))
  old_ec <- suppressWarnings(orig$fit_epigenomic_classifier(x$s, x$metrics, c("SNS", "ATAC")))
  expect_identical(ec$coefficients, old_ec$coefficients)
  expect_identical(ec$labeled_data, old_ec$labeled_data)
})

test_that("annotate_iz_bigwig summarises bigWig signal", {
  skip_on_cran()
  skip_if_not_installed("rtracklayer")
  bw <- tempfile(fileext = ".bw")
  g <- GenomicRanges::GRanges("chr1", IRanges::IRanges(seq(1L, 9001L, by = 1000L), width = 1000L),
                              score = as.numeric(1:10))
  GenomeInfoDb::seqlengths(g) <- c(chr1 = 10000L)
  ok <- tryCatch({ rtracklayer::export(g, bw, format = "bigWig"); TRUE },
                 error = function(e) FALSE)
  skip_if_not(ok, "cannot write bigWig files here")
  izm <- data.frame(IZ_id = 1:2, chr = "chr1", center = c(1500, 8000))
  res <- suppressMessages(annotate_iz_bigwig(izm, c(sig = bw), window_bp = 1000L))
  expect_named(res, c("IZ_id", "sig"))
  expect_equal(res$sig[1], 2, tolerance = 0.01)
  orig <- source_orig_downstream()
  expect_identical(res, suppressMessages(orig$annotate_iz_bigwig(izm, c(sig = bw), 1000L)))
  expect_error(annotate_iz_bigwig(izm, c(sig = "nope.bw")))
})
