marker_config <- function(...) {
  utils::modifyList(list(
    group_by_col = "state", use_sample_covariate = FALSE,
    min_cells_per_group = 50L, min_normalized_expression = 0.1,
    min_cell_fraction = 0.20, fdr_threshold = 0.01, aic_rule = "gam_lower",
    gam_basis_k = 10L, n_curve_points = 101L, n_clusters = 3L, workers = 1L
  ), list(...))
}

test_that("prepare_metadata and normalize_sparse_counts match the original", {
  orig <- source_orig_markers()
  d <- synth_markers()
  meta_file <- withr::local_tempfile(fileext = ".rds")
  saveRDS(d$meta, meta_file)
  old <- orig$prepare_metadata(meta_file, colnames(d$counts), "state", 0.01, 0.99)
  expect_identical(prepare_metadata(meta_file, colnames(d$counts)), old)
  expect_identical(prepare_metadata(d$meta, colnames(d$counts)), old)
  expect_setequal(unique(old$state), c("ctr", "ama"))
  expect_true(any(old$s_phase_pass) && !all(old$s_phase_pass))

  # an existing grouping column is used as is
  meta2 <- d$meta
  meta2$grp <- rep(c("a", "b", NA), length.out = nrow(meta2))
  saveRDS(meta2, meta_file)
  expect_identical(
    prepare_metadata(meta2, colnames(d$counts), "grp", 0.1, 0.9),
    orig$prepare_metadata(meta_file, colnames(d$counts), "grp", 0.1, 0.9)
  )
  expect_error(prepare_metadata(d$meta, colnames(d$counts), state_patterns = NULL),
               "grouping column")
  expect_error(prepare_metadata(d$meta[, -1], colnames(d$counts)), "required columns")
  expect_error(prepare_metadata(d$meta, c("x", "y")), "No matching cells")

  expect_identical(normalize_sparse_counts(d$counts),
                   orig$normalize_sparse_counts(d$counts))
  expect_identical(normalize_sparse_counts(d$counts, 1e6),
                   orig$normalize_sparse_counts(d$counts, 1e6))
  zero <- d$counts
  zero[, 1] <- 0
  expect_error(normalize_sparse_counts(zero), "zero or invalid library size")
})

test_that("per-gene model helpers match the original", {
  orig <- source_orig_markers()
  d <- synth_markers()
  meta <- prepare_metadata(d$meta, colnames(d$counts))
  norm <- normalize_sparse_counts(d$counts)
  inputs <- prepare_group_inputs("ctr", meta, norm$log_normalized, norm$normalized)
  expect_identical(
    inputs,
    orig$prepare_group_inputs("ctr", meta, norm$log_normalized, norm$normalized,
                              marker_config())
  )
  expect_error(
    prepare_group_inputs("ctr", meta, norm$log_normalized, norm$normalized,
                         min_cells_per_group = 1000L), "minimum is")

  grid <- seq(0.05, 0.95, length.out = 25)
  for (use_sample in c(FALSE, TRUE)) {
    for (g in c(1L, 2L, 3L)) {
      y <- inputs$log_expr[g, ]
      new <- build_gene_models(y, inputs$model_data, 10L, use_sample)
      old <- orig$build_gene_models(y, inputs$model_data, 10L, use_sample)
      expect_identical(new$sample_is_used, use_sample)
      expect_identical(stats::coef(new$gam), stats::coef(old$gam))
      expect_identical(summarize_gene_models(new), orig$summarize_gene_models(old))
      expect_identical(
        predict_marginal_curve(y, inputs$model_data, grid, 10L, use_sample),
        orig$predict_marginal_curve(y, inputs$model_data, grid, 10L, use_sample)
      )
    }
  }
  expect_error(build_gene_models(rep(1, nrow(inputs$model_data)), inputs$model_data),
               "constant")

  m <- rbind(a = c(1, 2, 3, 7), b = rep(5, 4), c = c(4, 1, 0, 2))
  expect_identical(scale_rows(m), orig$scale_rows(m))
  expect_identical(scale_rows(m)$valid, c(a = TRUE, b = FALSE, c = TRUE))
})

test_that("analyze_one_group and its steps match the original", {
  orig <- source_orig_markers()
  d <- synth_markers()
  meta <- prepare_metadata(d$meta, colnames(d$counts))
  norm <- normalize_sparse_counts(d$counts)

  for (use_sample in c(FALSE, TRUE)) {
    cfg <- marker_config(use_sample_covariate = use_sample, n_curve_points = 41L)
    for (grp in c("ctr", "ama")) {
      old <- suppressMessages(orig$analyze_one_group(
        grp, meta, norm$log_normalized, norm$normalized, cfg))
      new <- analyze_one_group(
        grp, meta, norm$log_normalized, norm$normalized,
        use_sample_covariate = use_sample, n_curve_points = 41L, verbose = FALSE)
      expect_identical(new, old)
      expect_gt(nrow(new$marker_table), 3)
      expect_true("Expression is constant or invalid." %in% new$gam_results$fit_status ||
                    all(new$gam_results$fit_status == "ok"))
    }
  }

  # the other AIC rule and a stricter expression filter: no curve is left
  cfg <- marker_config(aic_rule = "paper_literal_gam_greater", fdr_threshold = 1e-300)
  old <- suppressMessages(orig$analyze_one_group(
    "ctr", meta, norm$log_normalized, norm$normalized, cfg))
  new <- analyze_one_group("ctr", meta, norm$log_normalized, norm$normalized,
                           aic_rule = "paper_literal_gam_greater",
                           fdr_threshold = 1e-300, verbose = FALSE)
  expect_identical(new, old)
  expect_null(new$curves)
  expect_identical(nrow(new$marker_table), 0L)

  # parallel back-end gives the same table
  inputs <- prepare_group_inputs("ctr", meta, norm$log_normalized, norm$normalized)
  serial <- fit_gene_gams(inputs$log_expr, inputs$normalized_expr, inputs$model_data,
                          verbose = FALSE)
  expect_identical(
    fit_gene_gams(inputs$log_expr, inputs$normalized_expr, inputs$model_data,
                  BPPARAM = BiocParallel::MulticoreParam(2), verbose = FALSE),
    serial)
  expect_message(
    fit_gene_gams(inputs$log_expr[1:2, ], inputs$normalized_expr[1:2, ],
                  inputs$model_data, verbose = TRUE), "Fitting 2 genes")
  expect_error(classify_progression_markers(serial, "ctr", 10L, aic_rule = "x"))

  # a single curve cannot be clustered
  one <- compute_gene_curves(inputs$log_expr, inputs$model_data, "gene01")
  expect_null(one$hclust)
  expect_identical(
    one, orig$compute_gene_curves(inputs$log_expr, inputs$model_data, "gene01",
                                  marker_config()))
})

# Seurat objects for load_and_merge_rna_counts: one sample keeps the sample
# prefix on its cell names, the others have bare barcodes; one lacks a gene.
write_marker_seurat_objects <- function(d, dir) {
  for (s in d$samples) {
    cells <- startsWith(colnames(d$counts), paste0(s, "-"))
    m <- d$counts[, cells]
    if (s != d$samples[1]) colnames(m) <- sub(paste0(s, "-"), "", colnames(m), fixed = TRUE)
    if (s == d$samples[2]) m <- m[-1, ]
    obj <- suppressWarnings(SeuratObject::CreateSeuratObject(m, assay = "rna"))
    dir.create(file.path(dir, s))
    saveRDS(obj, file.path(dir, s, "wnn_rna_filtered.rds"))
  }
  invisible(dir)
}

test_that("load_and_merge_rna_counts merges per-sample Seurat objects", {
  skip_if_not_installed("SeuratObject")
  d <- synth_markers()
  dir <- write_marker_seurat_objects(d, withr::local_tempdir())
  new <- load_and_merge_rna_counts(d$samples, dir)
  expect_identical(colnames(new), colnames(d$counts))
  expect_identical(nrow(new), nrow(d$counts) - 1L)
  expect_error(load_and_merge_rna_counts("nope", dir), "Cannot find RNA object")
  expect_error(load_and_merge_rna_counts(d$samples[1], dir, assay = "atac"), "absent")
})

test_that("load_and_merge_rna_counts matches the original", {
  skip_if_not_installed("SeuratObject")
  # the original calls Seurat::Assays() / Seurat::GetAssayData(); under
  # `R CMD check --as-cran` Seurat cannot load (undeclared recommended
  # packages such as KernSmooth are hidden), so this comparison skips there.
  skip_if_not_installed("Seurat")
  orig <- source_orig_markers()
  d <- synth_markers()
  dir <- write_marker_seurat_objects(d, withr::local_tempdir())
  expect_identical(load_and_merge_rna_counts(d$samples, dir),
                   orig$load_and_merge_rna_counts(d$samples, dir))
})
