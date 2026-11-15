test_that("make_fixed_bins tiles chromosomes and normalises names", {
  si <- GenomeInfoDb::Seqinfo(c("1", "2", "X"), c(105000L, 60000L, 30000L))
  gr <- make_fixed_bins(si, 20000L, use_chr = c("chr1", 2))
  expect_s4_class(gr, "GRanges")
  expect_equal(as.character(GenomeInfoDb::seqnames(gr)),
               c(rep("chr1", 6), rep("chr2", 3)))
  expect_equal(BiocGenerics::end(gr)[6], 105000L)
  expect_equal(BiocGenerics::width(gr)[6], 5000L)
  expect_true(all(BiocGenerics::width(gr)[-c(6, 9)] == 20000L))
  expect_equal(GenomeInfoDb::seqlevels(gr), c("chr1", "chr2", "chrX"))
  expect_error(make_fixed_bins(si, 0), "bin_width")
  expect_error(make_fixed_bins(si, 100, use_chr = "chr5"), "chr5")
  expect_error(make_fixed_bins("x", 100), "Seqinfo")
})

make_map_inputs <- function() {
  fa <- tempfile(fileext = ".fa")
  seqs <- Biostrings::DNAStringSet(c(
    "1" = paste0(strrep("GC", 250), strrep("AT", 250))
  ))
  Biostrings::writeXStringSet(seqs, fa)
  si <- GenomeInfoDb::Seqinfo("chr1", 1000L)
  bins <- make_fixed_bins(si, 250L, use_chr = 1)
  map <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(c(1, 201, 601), c(200, 600, 1000)),
    score = c(1, 0.5, 0.9)
  )
  GenomeInfoDb::seqinfo(map) <- si
  list(fa = fa, bins = bins, map = map)
}

test_that("compute_bin_mappability computes GC and mappable fractions", {
  skip_if_not_installed("Biostrings")
  x <- make_map_inputs()
  out <- compute_bin_mappability(x$bins, x$map, th = 0.8, window_prop = 0.7,
                                 ref_path = x$fa)
  mc <- S4Vectors::mcols(out)
  expect_named(mc, c("bin_id", "gc_frequency", "mappability", "is_mappable"))
  expect_equal(mc$gc_frequency, c(1, 1, 0, 0))
  # bin1 1-250: 200 bp mappable; bin2 251-500: 0; bin3 501-750: 150; bin4: 250
  expect_equal(mc$mappability, c(200, 0, 150, 250) / 250)
  expect_equal(mc$is_mappable, c(TRUE, FALSE, FALSE, TRUE))
  expect_equal(mc$bin_id[1], "chr1-1-250")
})

test_that("compute_bin_mappability falls back when nothing passes th", {
  skip_if_not_installed("Biostrings")
  x <- make_map_inputs()
  out2 <- compute_bin_mappability(x$bins, x$map[2], th = 0.8, ref_path = x$fa)
  expect_true(all(!S4Vectors::mcols(out2)$is_mappable))
  expect_true("frac_mappable" %in% names(S4Vectors::mcols(out2)))
  expect_error(compute_bin_mappability(x$bins, x$map), "ref_path")
})

test_that("bins functions match the original make_fixed_bins.R", {
  skip_on_cran()
  skip_if_not_installed("Biostrings")
  orig <- source_orig("make_fixed_bins.R")
  si <- GenomeInfoDb::Seqinfo(c("chr1", "chr2"), c(105000L, 60000L))
  expect_identical(make_fixed_bins(si, 20000L, use_chr = 1:2),
                   orig$make_fixed_bins(si, 20000L, use_chr = 1:2))
  x <- make_map_inputs()
  for (th in c(0.8, 0.5, 1)) {
    expect_identical(
      compute_bin_mappability(x$bins, x$map, th = th, ref_path = x$fa),
      orig$compute_bin_mappability(x$bins, x$map, th = th, ref_path = x$fa,
                                   bpparam = BiocParallel::SerialParam())
    )
  }
})
