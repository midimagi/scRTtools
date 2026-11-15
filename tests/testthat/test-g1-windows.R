test_that("define_G1_windows builds windows of ~target_reads G1 fragments", {
  x <- synth_fragments()
  res <- define_G1_windows(x$bins, x$frag_path, x$md, target_reads = 20,
                           target_window_length = 5000, min.gapwidth = 1001L,
                           verbose = FALSE)
  g <- res$g1bins
  m <- res$g1bins_counts
  expect_s4_class(g, "GRanges")
  expect_s4_class(m, "dgCMatrix")
  expect_equal(dim(m), c(length(g), nrow(x$md)))
  expect_equal(rownames(m), names(g))
  expect_equal(colnames(m), rownames(x$md))
  expect_true(all(c("midpoint", "g1_read_count", "is_masked_window", "contig_id") %in%
                    names(S4Vectors::mcols(g))))
  md <- S4Vectors::metadata(g)
  expect_true(md$selected_g1_n >= 1 && md$selected_g1_n <= sum(x$md$is_G1))
  expect_equal(md$selected_g1_barcodes,
               rownames(x$md)[order(-x$md$coverage)][seq_len(md$selected_g1_n)])
  # no window overlaps the removed chr1 gap
  expect_false(any(IRanges::overlapsAny(
    g, GenomicRanges::GRanges("chr1", IRanges::IRanges(100001, 150000)))))
  # G1 counts of the selected cells in unmasked windows are target_reads (+-1)
  sel <- md$selected_g1_barcodes
  g1_counts <- Matrix::rowSums(m[, sel, drop = FALSE])
  unmasked <- !g$is_masked_window
  expect_true(all(abs(g1_counts[unmasked] - 20) <= 1))
  expect_error(define_G1_windows(x$bins, "nope.tsv.gz", x$md), "does not exist")
})

test_that("define_G1_windows matches the original script", {
  skip_on_cran()
  orig <- source_orig("define_G1_windows.R")
  x <- synth_fragments()
  for (tr in c(20, 50)) {
    new <- define_G1_windows(x$bins, x$frag_path, x$md, target_reads = tr,
                             target_window_length = 5000, min.gapwidth = 1001L,
                             verbose = FALSE)
    old <- suppressMessages(orig$define_G1_windows(
      x$bins, x$frag_path, x$md, target_reads = tr,
      target_window_length = 5000, min.gapwidth = 1001L, verbose = FALSE))
    expect_identical(new, old)
  }
})

# Notebook 1.3 logic verbatim (2-sample version), used as the reference.
notebook_sync <- function(gr_cn_list, gr_filtered_list) {
  gr_removed_list <- list()
  for (s in names(gr_cn_list)) {
    gr_cn <- gr_cn_list[[s]]
    hits <- GenomicRanges::findOverlaps(gr_cn, gr_filtered_list[[s]], type = "equal")
    bad_indices <- S4Vectors::queryHits(hits)
    gr_removed_list[[s]] <- gr_cn[-bad_indices]
  }
  gr_temp <- c(gr_removed_list[[1]], gr_removed_list[[2]])
  S4Vectors::mcols(gr_temp)$raw_cn <- NULL
  gr_removed_union <- unique(gr_temp)
  out <- list()
  for (s in names(gr_cn_list)) {
    gr_cn <- gr_filtered_list[[s]]
    hits <- GenomicRanges::findOverlaps(gr_cn, gr_removed_union, type = "equal")
    remove_indices <- S4Vectors::queryHits(hits)
    out[[s]] <- gr_cn[!(seq_along(gr_cn) %in% remove_indices)]
  }
  list(filtered = out, removed = gr_removed_list, removed_union = gr_removed_union)
}

test_that("synchronize_bins keeps only bins retained in all samples", {
  grl <- synth_windows()
  reads <- synth_reads(grl)
  reads2 <- synth_reads(grl, seed = 9L)
  cn1 <- BiocGenerics::unlist(count_reads(grl, write_reads_dir(reads), verbose = FALSE)$wins_count)
  cn2 <- BiocGenerics::unlist(count_reads(grl, write_reads_dir(reads2), verbose = FALSE)$wins_count)
  f1 <- BiocGenerics::unlist(filter_windows(S4Vectors::split(cn1, GenomicRanges::seqnames(cn1)),
                                            synth_gaps(), "raw_cn", std_th = 2))
  f2 <- BiocGenerics::unlist(filter_windows(S4Vectors::split(cn2, GenomicRanges::seqnames(cn2)),
                                            synth_gaps(), "raw_cn", std_th = 2))
  cn <- list(s1 = cn1, s2 = cn2)
  fl <- list(s1 = f1, s2 = f2)
  res <- synchronize_bins(cn, fl)

  key <- function(gr) paste(GenomicRanges::seqnames(gr), BiocGenerics::start(gr),
                            BiocGenerics::end(gr))
  shared <- intersect(key(f1), key(f2))
  expect_setequal(key(res$filtered$s1), shared)
  expect_setequal(key(res$filtered$s2), shared)
  expect_true("raw_cn" %in% names(S4Vectors::mcols(res$filtered$s1)))
  expect_equal(ncol(S4Vectors::mcols(res$removed_union)), 0L)

  ref <- notebook_sync(cn, fl)
  expect_identical(res$filtered, ref$filtered)
  expect_identical(res$removed, ref$removed)
  expect_identical(GenomicRanges::granges(res$removed_union), GenomicRanges::granges(ref$removed_union))

  # pooling only s1 leaves s1's own filter applied to both
  res1 <- synchronize_bins(cn, fl, union_samples = "s1")
  expect_setequal(key(res1$filtered$s1), key(f1))
  expect_error(synchronize_bins(cn, list(s1 = f1)), "same names")
  expect_error(synchronize_bins(cn, fl, union_samples = "s3"), "union_samples")
})
