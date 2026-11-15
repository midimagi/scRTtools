# Synthetic inputs and original-code loaders for the RNA tests.

# calculate_sc2p.R + gex_tempo_with_cp10k.R in one environment, so the
# original gex_tempo() finds SC2P_gene_unit() the way the notebooks do.
source_orig_gex <- function() {
  for (pkg in c("dplyr", "tibble", "tidyr")) testthat::skip_if_not_installed(pkg)
  env <- new.env(parent = asNamespace("dplyr"))
  for (f in c("calculate_sc2p.R", "gex_tempo_with_cp10k.R")) {
    path <- file.path(orig_rcode_dir(), "analysis", f)
    testthat::skip_if_not(file.exists(path), paste("original", f, "not available"))
    sys.source(path, envir = env)
  }
  env
}

# The progression-marker script only exists in the scRT_learning tree.
source_orig_markers <- function() {
  for (pkg in c("dplyr", "mgcv")) testthat::skip_if_not_installed(pkg)
  path <- file.path(
    local_dir("SCRTTOOLS_ORIG_LEARNING_RCODE"),
    "processing", "identify_S_phase_progression_markers.R"
  )
  testthat::skip_if_not(file.exists(path), "original progression-marker script not available")
  env <- new.env(parent = asNamespace("dplyr"))
  env$s <- mgcv::s
  sys.source(path, envir = env)
  env
}

synth_gex <- function(n_genes = 60L, n_cells = 120L, seed = 21L) {
  withr::with_seed(seed, {
    progress <- stats::runif(n_cells, 0.02, 0.98)
    mu <- c(rep(0.2, 15), rep(1.5, 15), rep(6, 15), rep(25, 10), rep(0, 5))[seq_len(n_genes)]
    slope <- stats::rnorm(n_genes, 0, 0.8)
    lam <- outer(mu, progress, function(m, p) m) * exp(outer(slope, progress - 0.5))
    m <- matrix(stats::rpois(n_genes * n_cells, lam), n_genes, n_cells,
                dimnames = list(sprintf("gene%02d", seq_len(n_genes)),
                                sprintf("cell%03d", seq_len(n_cells))))
    # Zero-inflate, as in real single-cell counts.
    m[stats::runif(length(m)) < 0.3] <- 0L
    meta <- data.frame(rep_percent = progress, other = seq_len(n_cells),
                       row.names = colnames(m))
    meta <- meta[sample(n_cells), , drop = FALSE]
    list(counts = m, meta = meta)
  })
}

# Two groups (ctr / ama) of two samples each; metadata row names lack the
# sample prefix that the RNA cell IDs carry.
synth_markers <- function(n_genes = 36L, n_per_sample = 45L, seed = 31L) {
  withr::with_seed(seed, {
    samples <- c("scRT6-mESC-rep1", "scRT6-mESC-rep2",
                 "scRT6-mESC-ama-rep1", "scRT6-mESC-ama-rep2")
    n_cells <- n_per_sample * length(samples)
    sample <- rep(samples, each = n_per_sample)
    barcode <- sprintf("BC%04d-1", seq_len(n_cells))
    progress <- stats::runif(n_cells, 0, 1)
    shape <- rep(c("up", "down", "peak", "flat"), length.out = n_genes)
    base <- stats::runif(n_genes, 2, 12)
    lam <- vapply(progress, function(p) {
      base * ifelse(shape == "up", 0.3 + 1.6 * p,
             ifelse(shape == "down", 1.9 - 1.6 * p,
             ifelse(shape == "peak", 0.4 + 1.5 * exp(-((p - 0.5) / 0.15)^2), 1)))
    }, numeric(n_genes))
    lam <- lam * rep(c(1, 1.3, 0.9, 1.2), each = n_per_sample)[col(lam)]
    counts <- matrix(stats::rpois(length(lam), lam), n_genes, n_cells)
    counts[n_genes, ] <- 3L
    counts[n_genes - 1L, ] <- c(1L, rep(0L, n_cells - 1L))  # barely expressed
    dimnames(counts) <- list(sprintf("gene%02d", seq_len(n_genes)),
                             paste0(sample, "-", barcode))
    meta <- data.frame(
      sample = sample,
      is_G1 = stats::runif(n_cells) < 0.1,
      rep_percent = progress,
      pass_qc = stats::runif(n_cells) > 0.05,
      row.names = barcode,
      stringsAsFactors = FALSE
    )
    meta$pass_qc[3] <- NA
    meta$rep_percent[5] <- NA
    list(counts = methods::as(counts, "CsparseMatrix"), meta = meta, samples = samples)
  })
}
