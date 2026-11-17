test_that("the example dataset is consistent and usable", {
  ex <- scrt_example_data()
  expect_named(ex, c("bins", "bin_counts", "g1_windows", "window_counts", "rep_states",
                     "rna_counts", "genes", "meta", "info"))
  cells <- rownames(ex$meta)
  s_cells <- cells[!ex$meta$is_G1]
  expect_identical(length(cells), 40L)
  expect_identical(sum(ex$meta$is_G1), 12L)
  expect_identical(dimnames(ex$bin_counts), list(names(ex$bins), cells))
  expect_identical(dimnames(ex$window_counts), list(names(ex$g1_windows), cells))
  expect_identical(dimnames(ex$rep_states), list(names(ex$bins), s_cells))
  expect_identical(dimnames(ex$rna_counts), list(ex$genes$gene, cells))
  expect_true(all(ex$rep_states %in% c(2L, 4L, NA)))
  chrs <- paste0("chr", 1:19)
  expect_identical(ex$info$chromosomes, chrs)
  expect_identical(length(ex$bins), 123148L)
  expect_setequal(as.character(GenomicRanges::seqnames(ex$bins)), chrs)
  expect_setequal(as.character(GenomicRanges::seqnames(ex$g1_windows)), chrs)
  expect_true(all(ex$genes$chr %in% chrs))
  expect_true(all(BiocGenerics::width(ex$bins) <= 20000L))
  expect_false(is.unsorted(ex$meta$rep_percent[!ex$meta$is_G1]))

  # the replicated fraction of the shipped states follows rep_percent
  rep_fraction <- colMeans(ex$rep_states == 4L, na.rm = TRUE)
  expect_false(anyNA(rep_fraction))
  expect_gt(stats::cor(rep_fraction, ex$meta[s_cells, "rep_percent"]), 0.99)

  # package functions run on it
  res <- gex_tempo(ex$rna_counts[, s_cells], ex$meta, nbin_rep = 4, method = "cp10k")
  expect_identical(dim(res$df_cp10k_mean), c(nrow(ex$genes), 4L))
  expect_identical(sum(attr(res, "rep_bins")$n_cells), length(s_cells))
  mapd <- calculate_mapd(ex$bin_counts, ex$bins, ex$meta["cell_id"], verbose = FALSE)
  expect_true(all(is.finite(mapd$mapd)))
  expect_true(all(is.na(mapd$mapd_error)))
})

test_that("the vignette data are downloaded, checked and cached", {
  cache <- withr::local_tempdir()
  # a file with the wrong content is rejected and not cached
  bad <- withr::local_tempfile(fileext = ".rds")
  saveRDS(1:3, bad)
  withr::local_options(scRTtools.vignette_data_url = paste0("file://", bad))
  expect_error(scrt_example_data("vignette", cache_dir = cache), "MD5")
  expect_length(list.files(cache), 0L)
  expect_error(
    withr::with_options(list(scRTtools.vignette_data_url = "file:///no/such/file.rds"),
                        suppressWarnings(scrt_example_data("vignette", cache_dir = cache))),
    "Could not download"
  )

  # the real file, when a local copy is available (SCRTTOOLS_VIGNETTE_DATA)
  local_copy <- Sys.getenv("SCRTTOOLS_VIGNETTE_DATA")
  skip_if_not(nzchar(local_copy) && file.exists(local_copy), "no local copy of the vignette data")
  withr::local_options(scRTtools.vignette_data_url = paste0("file://", normalizePath(local_copy)))
  ex <- suppressMessages(scrt_example_data("vignette", cache_dir = cache))
  expect_named(ex, c("bins", "bin_counts", "g1_windows", "window_counts", "rna_counts",
                     "genes", "bulk_rt", "meta", "info"))
  expect_identical(dim(ex$bin_counts), c(123148L, 500L))
  expect_identical(list.files(cache), "scRT1_vignette.rds")
  # the second call reads the cache, so a broken URL does not matter
  withr::local_options(scRTtools.vignette_data_url = "file:///no/such/file.rds")
  expect_identical(dim(scrt_example_data("vignette", cache_dir = cache)$meta), c(500L, 13L))
})
