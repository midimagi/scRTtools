# Synthetic bins x cells state matrix with a shared replication-timing order.
synth_states <- function(n_bin_chr = 120L, n_cells = 40L, seed = 6L) {
  withr::with_seed(seed, {
    chrs <- c("chr1", "chr2", "chr3")
    gr <- GenomicRanges::GRanges(
      rep(chrs, each = n_bin_chr),
      IRanges::IRanges(rep(seq(1L, by = 20000L, length.out = n_bin_chr), 3), width = 20000L)
    )
    names(gr) <- paste(GenomicRanges::seqnames(gr), BiocGenerics::start(gr), sep = "-")
    # smooth RT profile per chromosome
    rt <- unlist(lapply(chrs, function(i) {
      x <- cumsum(stats::rnorm(n_bin_chr))
      (x - min(x)) / diff(range(x))
    }))
    p <- sort(stats::runif(n_cells, 0.05, 0.95))
    st <- vapply(p, function(pi) {
      s <- ifelse(rt < pi, 4, 2)
      flip <- stats::runif(length(s)) < 0.03
      s[flip] <- 6 - s[flip]
      s[stats::runif(length(s)) < 0.05] <- NA
      s
    }, numeric(length(gr)))
    st[, 5] <- sample(c(2, 4), nrow(st), TRUE)                        # noise cell
    st[, 10] <- ifelse(as.character(GenomicRanges::seqnames(gr)) == "chr1", 2, 4)  # CNV-like
    st[, 1] <- 2                                                      # no replication
  })
  cells <- sprintf("cell%02d", seq_len(n_cells))
  dimnames(st) <- list(names(gr), cells)
  md <- data.frame(row.names = cells, cell_id = cells)
  list(gr = gr, st = st, md = md)
}
