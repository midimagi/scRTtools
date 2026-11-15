# Package functions against the stored results of the original scripts on the
# real-data example (chr1-chr19 of 40 cells): cell QC through initiation zones,
# replication tempo, RNA, group comparison and HC3 statistics. The reference
# file is built by
# data-raw/make_regression_reference.R and lives outside the package.

test_that("the pipeline reproduces the original scripts on real data", {
  skip_on_cran()
  skip_if_not_installed("mgcv")
  ref_file <- regression_ref_file()
  skip_if_not(file.exists(ref_file), "regression reference results not available")
  ref <- readRDS(ref_file)
  example_md5 <- unname(tools::md5sum(
    system.file("extdata", "scRT1_example.rds", package = "scRTtools")))
  skip_if_not(identical(attr(ref, "example_md5"), example_md5),
              "reference results were built from another example dataset")

  new <- regression_pipeline(regression_impl_package())
  expect_identical(names(new), names(ref))
  for (step in setdiff(names(new), "gene_bins")) {
    expect_equal(new[[!!step]], ref[[!!step]], tolerance = 1e-8)
  }

  # make_bin: the 6.6 original writes round bin ends in scientific notation
  # (e.g. chr19-3980001-4e+06); everywhere else the two agree, and with scipen
  # raised the original gives the package's IDs (checked in test-hc3-stats.R)
  sci <- grepl("e+", ref$gene_bins, fixed = TRUE)
  expect_identical(new$gene_bins[!sci], ref$gene_bins[!sci])
  expect_gt(sum(sci), 0)
  expect_true(all(grepl("0000$", new$gene_bins[sci])))
  expect_true(all(new$gene_bins %in% names(scrt_example_data()$bins)))
  expect_false(all(ref$gene_bins %in% names(scrt_example_data()$bins)))

  # the data are large enough for the comparison to mean something
  expect_gt(sum(new$g1$is_G1), 0)
  expect_gt(mean(new$states$genobins_rep_states == 4L, na.rm = TRUE), 0.3)
  expect_gt(sum(lengths(new$tracks)), 300)
  expect_gt(length(new$rt_peaks), 50)
  expect_gt(length(new$izones$izones), 5)
  expect_gt(nrow(new$track_gaps), 100)
  expect_identical(nrow(new$rep_tempo), length(new$izones$izones))
  expect_gt(sum(new$tempo_fit$sslogis$fit_status == "ok"), 0)
  expect_gt(sum(is.finite(new$sc2p$phase_II_mean)), 20)
  expect_gt(nrow(new$gex_tempo$cp10k$df_cp10k_mean), 50)
  expect_gt(sum(new$markers$result$gam_results$fit_status == "ok"), 100)
  expect_gt(nrow(new$markers$result$marker_table), 2)
  expect_false(is.null(new$markers$result$curves))

  gc <- new$group_compare
  expect_gt(nrow(gc$wide$drop), 50)
  expect_gt(sum(gc$diff_groups$iz_group != "no_sig"), 0)
  expect_gt(nrow(gc$tempo_feature), 5)
  expect_identical(length(gc$tempo_granges), nrow(gc$tempo_feature))
  expect_setequal(stats::na.omit(gc$segments$rt_class), c("early", "late"))
  expect_false(is.null(gc$intersection))
  expect_identical(length(gc$length_labels), nrow(gc$velocity))

  hc3 <- new$hc3
  expect_gt(hc3$fits$half1$n_genes, 500)
  expect_identical(sum(hc3$fit_table$status == "ok"), 8L)
  expect_identical(sum(hc3$fit_table$status == "rank_or_n"), 2L)
  expect_identical(hc3$combined_keep$status,
                   c("fit_failed", "ok", "ok", "ok", "ok"))
  expect_identical(nrow(hc3$combined_drop), 4L)
  expect_identical(sort(hc3$permutation), seq_along(hc3$permutation))
  expect_true(any(hc3$permutation != seq_along(hc3$permutation)))
})

test_that("recomputed states agree with the states of the original analysis", {
  # The shipped rep_states came from the genome-wide run (all cells, all
  # chromosomes), so only close agreement is expected here.
  ex <- scrt_example_data()
  s_cells <- colnames(ex$rep_states)
  res <- assign_replication_states(
    Matrix::Matrix(ex$window_counts[, s_cells], sparse = TRUE), ex$g1_windows, ex$bins,
    infer_args = regression_infer_args(), verbose = FALSE)
  agree <- mean(res$genobins_rep_states == ex$rep_states, na.rm = TRUE)
  expect_gt(agree, 0.85)
})
