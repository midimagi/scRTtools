# Small synthetic inputs for the cell QC / G1 window tests.

# Tile matrix + annotated bins for calculate_mapd().
synth_tiles <- function(n_bins = 400L, n_cells = 30L, seed = 3L) {
  withr::with_seed(seed, {
    chr <- rep(c("chr1", "chr2"), each = n_bins / 2)
    start <- c(seq(1L, by = 1000L, length.out = n_bins / 2),
               seq(1L, by = 1000L, length.out = n_bins / 2))
    gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(start, width = 1000L))
    names(gr) <- paste(chr, start, start + 999L, sep = "-")
    gr$gc_frequency <- stats::runif(n_bins, 0.3, 0.6)
    gr$mappability <- stats::runif(n_bins, 0.6, 1)
    gr$is_mappable <- gr$mappability > 0.7
    depth <- stats::runif(n_cells, 2, 10)
    # half the cells carry a noisy (S-phase like) profile
    noise <- rep(c(0, 1), length.out = n_cells)
    m <- vapply(seq_len(n_cells), function(j) {
      prof <- 1 + noise[j] * rep(c(0, 1), each = 25, length.out = n_bins)
      stats::rpois(n_bins, depth[j] * prof * gr$mappability)
    }, numeric(n_bins))
    dimnames(m) <- list(names(gr), sprintf("cell%02d", seq_len(n_cells)))
  })
  md <- data.frame(cell_id = colnames(m), sample = rep(c("s1", "s2"), length.out = n_cells))
  rownames(md) <- md$cell_id
  list(gr = gr, m = Matrix::Matrix(m, sparse = TRUE), md = md)
}

# coverage / scaled MAPD for find_g1_cells(): two lines.
synth_mapd <- function(n = 200L, seed = 4L) {
  withr::with_seed(seed, {
    cov <- stats::runif(n, 0.5, 3)
    g1 <- rep(c(TRUE, FALSE), c(round(n * 0.6), n - round(n * 0.6)))
    smapd <- ifelse(g1, 0.2 - 0.02 * cov, 0.5 - 0.02 * cov) + stats::rnorm(n, 0, 0.01)
  })
  md <- data.frame(cell_id = sprintf("c%03d", seq_len(n)), coverage = cov,
                   scaled_mapd = smapd, sample = rep(c("s1", "s2"), length.out = n))
  rownames(md) <- md$cell_id
  list(md = md, truth = g1)
}

# Fragment file (bgzip + tabix) with uniform fragments for n_cells cells on
# two chromosomes; the first n_g1 cells are G1 with decreasing coverage.
synth_fragments <- function(chr_len = 300000L, n_cells = 30L, n_g1 = 15L, seed = 5L) {
  cells <- sprintf("BC%02d-1", seq_len(n_cells))
  withr::with_seed(seed, {
    frags <- do.call(rbind, lapply(c("chr1", "chr2"), function(chr) {
      do.call(rbind, lapply(seq_len(n_cells), function(j) {
        n <- stats::rpois(1, 150 + 10 * (n_cells - j))
        s <- sort(sample.int(chr_len - 500L, n))
        data.frame(chr = chr, start = s - 1L, end = s - 1L + sample(80:300, n, TRUE),
                   barcode = cells[j], count = 1L)
      }))
    }))
  })
  frags <- frags[order(frags$chr, frags$start, frags$end), ]
  tsv <- tempfile(fileext = ".tsv")
  utils::write.table(frags, tsv, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = FALSE)
  gz <- Rsamtools::bgzip(tsv, dest = paste0(tsv, ".gz"))
  Rsamtools::indexTabix(gz, seq = 1, start = 2, end = 3, zeroBased = TRUE)
  md <- data.frame(row.names = cells,
                   is_G1 = seq_len(n_cells) <= n_g1,
                   coverage = rev(seq_len(n_cells)) + 0.5)
  # retained bins: 1 kb tiles minus a gap on chr1
  si <- GenomeInfoDb::Seqinfo(c("chr1", "chr2"), c(chr_len, chr_len))
  bins <- make_fixed_bins(si, 1000L, use_chr = 1:2)
  bins <- bins[!IRanges::overlapsAny(
    bins, GenomicRanges::GRanges("chr1", IRanges::IRanges(100001, 150000)))]
  list(frag_path = gz, md = md, bins = bins, frags = frags)
}

