# Loader for the 6.6 / 6.7 analysis scripts and synthetic regression data.
# The scripts are flat pipelines that cannot be sourced, so only the named
# top-level assignments are evaluated.

script_objects <- function(script, object_names) {
  testthat::skip_if_not_installed("data.table")
  path <- file.path(notebook_dir(), script)
  testthat::skip_if_not(file.exists(path), paste("script", script, "not available"))
  # data.table syntax inside the functions needs a data.table-aware scope.
  env <- new.env(parent = asNamespace("data.table"))
  env$ns <- splines::ns
  for (e in parse(path, keep.source = FALSE)) {
    if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") && is.name(e[[2]]) &&
        as.character(e[[2]]) %in% object_names) {
      eval(e, envir = env)
    }
  }
  missing <- setdiff(object_names, ls(env))
  if (length(missing)) stop("Not found in ", script, ": ", paste(missing, collapse = ", "))
  env
}

source_orig_66 <- function() {
  script_objects("6.6.Chek1_paired_DNA_RNA_analysis.R",
                 c("zscore", "make_bin", "ols_hc3", "base_formula", "make_base",
                   "meta_combine", "gene_design", "count_hc3"))
}

source_orig_67 <- function() {
  script_objects("6.7.expression_to_other_gene_replication.R",
                 c("zscore", "base_matrix", "fit_hc3", "combine_two"))
}

# Cells of two samples with the covariates the 6.6 / 6.7 models use.
synth_cells_hc3 <- function(n = 400L, seed = 51L) {
  withr::with_seed(seed, {
    d <- data.frame(
      cell_id = sprintf("cell%04d", seq_len(n)),
      sample = rep(c("scRT1", "scRT2"), length.out = n),
      progress = stats::runif(n, 0.1, 0.9),
      log_rna = stats::rnorm(n, 9, 0.4),
      log_dna = stats::rnorm(n, 11, 0.3),
      complexity = stats::rnorm(n),
      call_fraction = stats::runif(n, 0.8, 1),
      score_S = stats::rnorm(n),
      score_G2M = stats::rnorm(n),
      score_pluripotency = stats::rnorm(n),
      stringsAsFactors = FALSE
    )
    d$chek1_local <- as.numeric(stats::runif(n) < d$progress)
    d$progress_c <- d$progress - 0.5
    d$chek1_z <- stats::rnorm(n)
    d$local <- as.numeric(stats::runif(n) < d$progress)
    d$local[sample(n, 12)] <- NA
    d$expr <- 0.4 * d$local + 0.2 * d$chek1_z + 0.5 * d$progress +
      stats::rnorm(n, sd = 0.5 + d$progress)
    d$count <- stats::rpois(n, exp(d$log_rna - 9 + 0.3 * d$chek1_z + 0.4 * d$progress))
    d$expr[7] <- NA
    d
  })
}
