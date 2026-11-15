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
