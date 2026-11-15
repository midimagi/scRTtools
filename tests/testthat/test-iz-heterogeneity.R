het_inputs <- function() {
  f <- iz_analysis_fixture()
  d <- prepare_iz_data(f$diz$tracks[lengths(f$diz$tracks) > 0], f$fr)
  c(f, list(d = d))
}

test_that("IZ heterogeneity pipeline runs on synthetic data", {
  h <- het_inputs()
  d <- h$d
  expect_named(d, c("track_intervals", "cell_metadata", "iz_metadata", "diagnostics"))
  ti <- d$track_intervals
  expect_false(anyDuplicated(paste(ti$cell_id, ti$IZ_id)) > 0)
  expect_true(all(abs(ti$fork_asymmetry) <= 1, na.rm = TRUE))
  cells <- d$cell_metadata$cell_id
  izs <- d$iz_metadata$IZ_id
  u <- call_iz_usage(ti, cells, izs)
  expect_identical(dim(u), c(length(cells), length(izs)))
  us <- call_iz_usage_from_states(h$x$st, cells, d$iz_metadata)
  expect_true(all(us %in% c(0L, 1L, NA)))
  prog <- d$cell_metadata$replication_progress
  uh <- compute_usage_heterogeneity(u, prog)
  expect_true(all(uh$summary$usage_heterogeneity >= 0 & uh$summary$usage_heterogeneity <= 1))
  uh_s <- compute_usage_heterogeneity(us, prog)
  expect_true("n_usage_observed" %in% names(uh_s$summary))
  tm <- fit_timing_models(us, prog, uh_s$probability)
  expect_identical(tm$IZ_id, izs)
  expect_true(any(tm$timing_flag == "ok" | grepl("separation", tm$timing_flag)))
  us_na <- us; us_na[, 1] <- NA
  expect_identical(fit_timing_models(us_na, prog)$timing_flag[1], "no_observations")
  tf <- summarize_track_features(ti)
  expect_identical(nrow(tf$summary), length(unique(ti$IZ_id)))
  sh <- summarize_track_shape(ti, d$iz_metadata, window_bp = 200000L)
  expect_true(all(sh$shape_flag %in% c("ok", "too_few_active")))
  tens <- build_replication_tensor(ti, cells, d$iz_metadata, izs[1:3], window_bp = 100000L)
  expect_identical(dim(tens), c(length(cells), 3L, 11L))
  ct <- correct_replication_tensor(tens, prog, keep_fitted = TRUE)
  expect_equal(ct$residual + ct$fitted, tens, ignore_attr = TRUE)
  expect_error(build_replication_tensor(ti, cells, d$iz_metadata, izs, max_elements = 10), "Dense")
  all_sum <- assemble_iz_summary(d$iz_metadata, uh$summary, tm, tf$summary, sh)
  expect_identical(all_sum$IZ_id, sort(izs))
  expect_true(all(c("usage_heterogeneity", "S50", "position_SD", "track_shape_heterogeneity") %in%
                    names(all_sum)))
})

test_that("IZ heterogeneity functions match the original scripts", {
  skip_on_cran()
  orig <- source_orig_heterogeneity()
  h <- het_inputs()
  tr <- h$diz$tracks[lengths(h$diz$tracks) > 0]
  md <- h$x$md
  md$rep_percent[2] <- NA
  for (cm in list(NULL, md)) {
    expect_identical(prepare_iz_data(tr, h$fr, cm), orig$prepare_iz_data(tr, h$fr, cm))
  }
  d <- h$d
  ti <- d$track_intervals
  cells <- d$cell_metadata$cell_id
  izs <- d$iz_metadata$IZ_id
  prog <- d$cell_metadata$replication_progress
  expect_identical(call_iz_usage(ti, cells, izs), orig$call_iz_usage(ti, cells, izs))
  us <- call_iz_usage_from_states(h$x$st, cells, d$iz_metadata)
  expect_identical(us, orig$call_iz_usage_from_states(h$x$st, cells, d$iz_metadata))
  expect_identical(make_progress_bins(prog, 7), orig$make_progress_bins(prog, 7))
  uh <- compute_usage_heterogeneity(us, prog)
  expect_identical(uh, orig$compute_usage_heterogeneity(us, prog))
  # track-based (sparse, no NA) usage: equal to the phase1_3 version too
  u <- call_iz_usage(ti, cells, izs)
  ph13 <- new.env(parent = asNamespace("dplyr"))
  sys.source(file.path(notebook_dir(), "iz_heterogeneity_phase1_3.R"), envir = ph13,
             toplevel.env = ph13)
  a <- compute_usage_heterogeneity(u, prog)
  b <- ph13$compute_usage_heterogeneity(u, prog)
  expect_equal(a$summary[names(b$summary)], b$summary, tolerance = 1e-12)
  expect_identical(a$probability, b$probability)
  expect_identical(fit_timing_models(us, prog, uh$probability),
                   orig$fit_timing_models(us, prog, uh$probability))
  for (m in c("spline", "bin")) {
    expect_identical(summarize_track_features(ti, m), orig$summarize_track_features(ti, m))
    expect_identical(summarize_track_shape(ti, d$iz_metadata, m, window_bp = 200000L),
                     orig$summarize_track_shape(ti, d$iz_metadata, m, window_bp = 200000L))
  }
  tens <- build_replication_tensor(ti, cells, d$iz_metadata, izs[1:4], window_bp = 100000L)
  expect_identical(tens, orig$build_replication_tensor(ti, cells, d$iz_metadata, izs[1:4],
                                                       window_bp = 100000L))
  for (m in c("bin", "spline")) {
    expect_identical(correct_replication_tensor(tens, prog, m, keep_fitted = TRUE),
                     orig$correct_replication_tensor(tens, prog, m, keep_fitted = TRUE))
  }
  tm <- fit_timing_models(us, prog, uh$probability)
  tf <- summarize_track_features(ti)$summary
  sh <- summarize_track_shape(ti, d$iz_metadata, window_bp = 200000L)
  expect_identical(assemble_iz_summary(d$iz_metadata, uh$summary, tm, tf, sh),
                   orig$assemble_iz_summary(d$iz_metadata, uh$summary, tm, tf, sh))
})
