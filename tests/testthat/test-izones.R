test_that(".parse_coordinate parses lengths", {
  expect_identical(scRTtools:::.parse_coordinate("500Kb"), 5e5)
  expect_identical(scRTtools:::.parse_coordinate(" 1 Mb"), 1e6)
  expect_identical(scRTtools:::.parse_coordinate("0"), 0)
  expect_identical(scRTtools:::.parse_coordinate(20000L), 20000)
  expect_error(scRTtools:::.parse_coordinate("3gb"), "Unsupported")
})

test_that("define_initiation_zones clusters tracks around origins", {
  f <- izone_fixture()
  iz <- f$diz$izones
  expect_s4_class(iz, "GRanges")
  expect_gt(length(iz), 10L)
  expect_identical(iz$iz_id, seq_along(iz))
  expect_false(is.unsorted(iz))
  expect_true(all(iz$n_tracks >= 5L))
  expect_named(S4Vectors::mcols(iz), c("iz_id", "q5", "q25", "q50", "q75", "q95", "n_tracks"))
  tr <- BiocGenerics::unlist(f$diz$tracks, use.names = FALSE)
  expect_identical(lengths(f$diz$tracks), lengths(f$x$tracks))
  expect_identical(as.vector(table(factor(tr$iz_id, levels = iz$iz_id))), iz$n_tracks)
  # IZ centres lie close to the simulated origins (also unseeded ones)
  org <- f$x$gr[f$x$origins$bin]
  d <- GenomicRanges::distanceToNearest(GenomicRanges::resize(iz, 1L, fix = "center"), org)
  expect_lt(stats::median(S4Vectors::mcols(d)$distance), 60000)
  expect_identical(
    define_initiation_zones(f$x$tracks, f$x$peaks, min_usage = 5L, verbose = FALSE,
                            BPPARAM = BiocParallel::MulticoreParam(2)),
    f$diz
  )
  empty <- define_initiation_zones(f$x$tracks, f$x$peaks, min_usage = 1e6, verbose = FALSE)
  expect_length(empty$izones, 0L)
  expect_true(all(is.na(BiocGenerics::unlist(empty$tracks)$iz_id)))
  bad <- f$x$tracks; names(bad) <- NULL
  expect_error(define_initiation_zones(bad, f$x$peaks, verbose = FALSE), "names")
})

test_that("firing range, firing order and predicted order", {
  f <- izone_fixture()
  fr <- calculate_firing_range(f$iz, f$x$st, f$diz$tracks, verbose = FALSE)
  expect_true(all(c("earliest_rep_percent", "latest_unrep_percent", "corroborated") %in%
                    names(S4Vectors::mcols(fr))))
  ok <- !is.na(fr$earliest_rep_percent) & !is.na(fr$latest_unrep_percent)
  expect_gt(sum(ok), 5L)
  expect_true(all(fr$latest_unrep_percent[ok] >= fr$earliest_rep_percent[ok]))
  fo <- calculate_firing_order(f$iz, f$x$st, f$x$md, rep_range = c(0.2, 0.8))
  expect_identical(dim(fo), c(length(f$iz), 4L))
  expect_equal(unname(rowSums(fo)), rep(100, nrow(fo)))
  expect_identical(rownames(fo)[1], {
    i <- order(f$iz$rt, decreasing = TRUE)[1]
    paste0(GenomicRanges::seqnames(f$iz)[i], ":", BiocGenerics::start(f$iz)[i], "-",
           BiocGenerics::end(f$iz)[i])
  })
  pf <- predict_firing_order(f$iz, f$diz$tracks)
  expect_identical(dimnames(pf), dimnames(calculate_firing_order(f$iz, f$x$st, f$x$md)))
  nb <- 5
  pf2 <- predict_firing_order(f$iz, f$diz$tracks, rep_range = c(0, 0.5),
                              track_filter = num_bins <= nb)
  expect_identical(dim(pf2), dim(pf))
  expect_error(calculate_firing_order(f$iz[, -8], f$x$st, f$x$md), "rt")
  expect_error(predict_firing_order(f$iz, f$diz$tracks, rep_range = c(1, 0)), "smaller")
})

test_that("fit_all_izone_velocities fits GAM velocity curves", {
  skip_if_not_installed("mgcv")
  f <- izone_fixture()
  tr <- BiocGenerics::unlist(f$diz$tracks, use.names = FALSE)
  df <- data.frame(iz_id = tr$iz_id, rep_percent = tr$rep_percent,
                   mid_track_len = BiocGenerics::width(tr))
  df <- df[!is.na(df$iz_id), ]
  res <- fit_all_izone_velocities(df, c(f$iz$iz_id, 999L), min_points_per_izone = 8L)
  expect_named(res, c("summary", "curves", "config"))
  s <- res$summary
  expect_identical(s$iz_id, c(f$iz$iz_id, 999L))
  expect_identical(s$fit_status[s$iz_id == 999L], "insufficient_points")
  expect_true(any(s$fit_status == "ok"))
  expect_true(any(s$velocity_direction == "increasing"))
  expect_identical(nrow(res$curves), 101L * sum(s$fit_status == "ok"))
  expect_identical(
    fit_all_izone_velocities(df, c(f$iz$iz_id, 999L), min_points_per_izone = 8L,
                             BPPARAM = BiocParallel::MulticoreParam(2)),
    res
  )
})

test_that("IZ functions match the original scripts", {
  skip_on_cran()
  skip_if_not_installed("mgcv")
  orig <- source_orig_izones()
  x <- synth_izone_data()
  for (a in list(list(min_usage = 5L), list(min_usage = 3L, short_track_th = "200Kb",
                                            max_track_length = 600000, min_overlap = "20kb"))) {
    new <- do.call(define_initiation_zones, c(list(x$tracks, x$peaks, verbose = FALSE), a))
    old <- do.call(orig$define_initiation_zones,
                   c(list(x$tracks, x$peaks, verbose = FALSE,
                          bpparam = BiocParallel::SerialParam()), a))
    expect_identical(new, old)
  }
  f <- izone_fixture(x)
  for (a in list(list(), list(min_rep_tracks = 6))) {
    expect_identical(
      do.call(calculate_firing_range, c(list(f$iz, x$st, f$diz$tracks, verbose = FALSE), a)),
      do.call(orig$calculate_firing_range, c(list(f$iz, x$st, f$diz$tracks, verbose = FALSE), a))
    )
  }
  for (rr in list(c(0, 1), c(0.25, 0.75))) {
    expect_identical(calculate_firing_order(f$iz, x$st, x$md, rr),
                     orig$calculate_firing_order(f$iz, x$st, x$md, rr))
  }
  expect_identical(predict_firing_order(f$iz, f$diz$tracks),
                   quiet_deprecated(orig$predict_firing_order(f$iz, f$diz$tracks)))
  expect_identical(predict_firing_order(f$iz, f$diz$tracks, c(0.2, 0.9), num_bins <= 6),
                   quiet_deprecated(
                     orig$predict_firing_order(f$iz, f$diz$tracks, c(0.2, 0.9), num_bins <= 6)))
  tr <- BiocGenerics::unlist(f$diz$tracks, use.names = FALSE)
  df <- data.frame(iz_id = tr$iz_id, rep_percent = tr$rep_percent,
                   mid_track_len = BiocGenerics::width(tr))
  df <- df[!is.na(df$iz_id), ]
  cfg <- list(min_points_per_izone = 8L, gam_basis_k = 10L, n_curve_points = 101L,
              fdr_threshold = 0.05, workers = 1L)
  old <- orig$fit_all_izone_velocities(df, f$iz$iz_id, cfg)
  new <- fit_all_izone_velocities(df, f$iz$iz_id, min_points_per_izone = 8L)
  expect_identical(new$summary, old$summary)
  expect_identical(new$curves, old$curves)
})
