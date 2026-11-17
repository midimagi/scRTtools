# Genome-bin counts with GC / mappability bias and an RT profile.
synth_pseudobulk <- function(n_chr = 2L, n_bin = 400L, n_g1 = 15L, n_s = 25L, seed = 21L) {
  withr::with_seed(seed, {
    chrs <- paste0("chr", seq_len(n_chr))
    gr <- GenomicRanges::GRanges(
      rep(chrs, each = n_bin),
      IRanges::IRanges(rep(seq(1L, by = 20000L, length.out = n_bin), n_chr), width = 20000L)
    )
    names(gr) <- paste(GenomicRanges::seqnames(gr), BiocGenerics::start(gr),
                       BiocGenerics::end(gr), sep = "-")
    n <- length(gr)
    gr$gc_frequency <- round(stats::runif(n, 0.35, 0.5), 3)
    gr$mappability <- stats::runif(n, 0.7, 1)
    gr$is_mappable <- TRUE
    gr$is_mappable[c(1:3, 120:125, 300, 401:404, 600:640)] <- FALSE
    gr$gc_frequency[7] <- NA
    rt <- unlist(lapply(chrs, function(i) {
      0.5 + 0.5 * sin(seq(0, 4 * pi, length.out = n_bin) + stats::runif(1, 0, pi))
    }))
    bias <- gr$mappability * (1 + 2 * (gr$gc_frequency - 0.42))
    bias[is.na(bias)] <- 1
    g1 <- vapply(seq_len(n_g1), function(j) stats::rpois(n, 8 * bias), numeric(n))
    p <- stats::runif(n_s, 0.2, 0.8)
    s <- vapply(p, function(pj) stats::rpois(n, 8 * bias * (1 + (rt > 1 - pj))), numeric(n))
    m <- cbind(g1, s)
  })
  cells <- sprintf("c%02d", seq_len(ncol(m)))
  dimnames(m) <- list(names(gr), cells)
  md <- data.frame(row.names = cells,
                   is_G1 = rep(c(TRUE, FALSE), c(n_g1, n_s)),
                   pass_qc = TRUE)
  md$pass_qc[n_g1 + 1:3] <- FALSE
  list(gr = gr, m = Matrix::Matrix(m, sparse = TRUE), md = md)
}

# Split synth_pseudobulk() into per-sample raw matrix files, as the S/G1 ratio
# functions read them: sample "sA" stores bare barcodes with shuffled rows,
# sample "sB" stores already prefixed cell ids.
synth_pseudobulk_files <- function(x = synth_pseudobulk(), dir = tempfile("pb_")) {
  dir.create(dir)
  n <- ncol(x$m)
  grp <- rep(c("sA", "sB"), length.out = n)
  md <- x$md
  rownames(md) <- paste(grp, colnames(x$m), sep = "-")
  # extra cells that are in the metadata but in no matrix are allowed
  md <- rbind(md, data.frame(row.names = "sC-zz", is_G1 = NA, pass_qc = TRUE))
  mA <- x$m[rev(seq_len(nrow(x$m))), grp == "sA"]
  mB <- x$m[, grp == "sB"]
  colnames(mB) <- paste0("sB-", colnames(mB))
  files <- c(sA = file.path(dir, "sA.rds"), sB = file.path(dir, "sB.rds"))
  saveRDS(mA, files[["sA"]])
  saveRDS(mB, files[["sB"]])
  list(files = files, mats = list(sA = mA, sB = mB), md = md, gr = x$gr, m = x$m, grp = grp)
}
