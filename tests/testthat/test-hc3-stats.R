test_that("fit_ols_hc3 matches fit_hc3 of 6.7", {
  o67 <- source_orig_67()
  d <- synth_cells_hc3()
  B <- full_rank_model_matrix(~ splines::ns(progress, df = 5) + log_rna + log_dna +
                                complexity + chek1_local, d)
  for (y in list(d$local, d$expr)) {
    old <- o67$fit_hc3(B, d$chek1_z, y)
    new <- fit_ols_hc3(cbind(B, exposure = d$chek1_z), y)
    expect_identical(new$term, "exposure")
    expect_identical(as.list(new[, names(old)]), lapply(old, unname))
  }
  # too few rows, rank deficiency, constant response
  few <- seq_len(ncol(B) + 15L)
  expect_identical(fit_ols_hc3(cbind(B, exposure = d$chek1_z)[few, ], d$expr[few])$status,
                   o67$fit_hc3(B[few, ], d$chek1_z[few], d$expr[few])$status)
  dup <- fit_ols_hc3(cbind(B, exposure = B[, 2]), d$expr)
  expect_identical(dup$status, o67$fit_hc3(B, B[, 2], d$expr)$status)
  expect_identical(dup$status, "rank_or_n")
  const <- fit_ols_hc3(cbind(B, exposure = d$chek1_z), rep(1, nrow(d)))
  old <- o67$fit_hc3(B, d$chek1_z, rep(1, nrow(d)))
  expect_identical(const$status, old$status)

  # one cell alone in a factor level: leverage 1
  d$lonely <- as.numeric(seq_len(nrow(d)) == 3L)
  B2 <- cbind(B, lonely = d$lonely)
  old <- o67$fit_hc3(B2, d$chek1_z, d$expr)
  new <- fit_ols_hc3(cbind(B2, exposure = d$chek1_z), d$expr)
  expect_identical(old$status, "high_leverage")
  expect_identical(as.list(new[, names(old)]), old)
})

test_that("fit_ols_hc3 with max_leverage = Inf matches ols_hc3 of 6.6", {
  o66 <- source_orig_66()
  d <- synth_cells_hc3()
  B <- full_rank_model_matrix(~ splines::ns(progress, df = 5) + log_rna + log_dna, d)
  X <- add_interaction_terms(B, d$chek1_z, d$local, d$progress_c)
  terms <- c("local_rep", "interaction")
  old <- as.data.frame(o66$ols_hc3(X, d$expr, terms))
  new <- fit_ols_hc3(X, d$expr, terms, max_leverage = Inf)
  expect_identical(new[c("term", "n", "status")], old[c("term", "n", "status")])
  expect_identical(new$beta, old$beta)
  expect_equal(new$se, old$se, tolerance = 1e-10)
  expect_equal(new$p, old$p, tolerance = 1e-9)
  expect_equal(new$max_hat, old$max_hat, tolerance = 1e-10)
  # the default agrees here because no leverage is extreme
  expect_identical(fit_ols_hc3(X, d$expr, terms), new)

  # leverage 1: 6.6 caps it, the default rejects the fit
  X2 <- cbind(X, lonely = as.numeric(seq_len(nrow(d)) == 3L))
  old <- as.data.frame(o66$ols_hc3(X2, d$expr, terms))
  new <- fit_ols_hc3(X2, d$expr, terms, max_leverage = Inf)
  expect_identical(new$status, old$status)
  expect_equal(new$se, old$se, tolerance = 1e-6)
  expect_identical(fit_ols_hc3(X2, d$expr, terms)$status, rep("high_leverage", 2))

  few <- seq_len(20)
  expect_identical(fit_ols_hc3(X[few, ], d$expr[few], terms)[c("term", "n", "status")],
                   as.data.frame(o66$ols_hc3(X[few, ], d$expr[few], terms))[c("term", "n", "status")])
  expect_error(fit_ols_hc3(X, d$expr, "nope"), "not found")
  expect_error(fit_ols_hc3(X, d$expr[-1]), "one value per row")
})

test_that("fit_poisson_hc3 matches count_hc3 of 6.6", {
  o66 <- source_orig_66()
  d <- synth_cells_hc3()
  B <- full_rank_model_matrix(~ splines::ns(progress, df = 4) + complexity, d)
  X <- add_interaction_terms(B, d$chek1_z, d$local, d$progress_c)
  off <- d$log_rna - 9
  old <- as.data.frame(o66$count_hc3(X, d$count, off))
  new <- fit_poisson_hc3(X, d$count, off)
  expect_identical(new$term, "interaction")
  expect_identical(new[names(old)], old)
  expect_identical(old$status, "ok")
  expect_identical(fit_poisson_hc3(X, d$count, off, "interaction"), new)
  two <- fit_poisson_hc3(X, d$count, off, c("regulator", "interaction"))
  expect_identical(two[2, -1], new[, -1], ignore_attr = TRUE)
  expect_identical(two[1, names(old)],
                   as.data.frame(o66$count_hc3(X, d$count, off, "regulator")),
                   ignore_attr = TRUE)

  # a dosage offset only shifts the local_rep main effect
  dose <- fit_poisson_hc3(X, d$count, off + log1p(ifelse(is.na(d$local), NA, d$local)))
  expect_equal(dose$beta, new$beta, tolerance = 1e-6)

  few <- seq_len(25)
  expect_identical(fit_poisson_hc3(X[few, ], d$count[few], off[few])[names(old)],
                   as.data.frame(o66$count_hc3(X[few, ], d$count[few], off[few])))
  # an extreme coefficient is flagged in both
  sep <- d$count
  sep[d$local %in% 1 & d$chek1_z > 1.2] <- 0L
  X3 <- cbind(B, interaction = as.numeric(d$local %in% 1 & d$chek1_z > 1.2))
  old <- as.data.frame(o66$count_hc3(X3, sep, off))
  expect_identical(fit_poisson_hc3(X3, sep, off)[names(old)], old)
  expect_identical(fit_poisson_hc3(X3, sep, off, max_abs_beta = 1)$status, "unstable")
})

make_fits <- function(seed = 61L, n = 40L) {
  withr::with_seed(seed, {
    tab <- data.frame(
      gene = rep(sprintf("g%02d", sample(n)), each = 2),
      term = "interaction",
      sample = rep(c("scRT2", "scRT1"), n),
      beta = stats::rnorm(2 * n, 0.1, 0.3),
      se = stats::runif(2 * n, 0.05, 0.3),
      status = "ok",
      stringsAsFactors = FALSE
    )
    tab$p <- 2 * stats::pnorm(-abs(tab$beta / tab$se))
    tab
  })
}

test_that("combine_two_samples matches combine_two of 6.7", {
  o67 <- source_orig_67()
  tab <- make_fits()
  tab[5, c("beta", "se", "p")] <- NA
  tab$status[5] <- "rank_or_n"
  tab$status[12] <- "high_leverage"
  tab$se[20] <- 0
  dt <- data.table::as.data.table(tab)
  old <- as.data.frame(dt[, o67$combine_two(.SD), by = c("gene", "term")])
  new <- combine_two_samples(tab, c("gene", "term"))
  expect_identical(new[names(old)], old)
  expect_identical(sum(new$status == "fit_failed"), 3L)
  expect_true(all(new$p[new$status == "fit_failed"] == 1))
  ok <- new$status == "ok"
  expect_equal(new$I2[ok], ifelse(new$Q[ok] > 0, pmax(0, (new$Q[ok] - 1) / new$Q[ok]), 0))
  expect_named(new, c("gene", "term", "beta", "se", "p", "replication_p", "same_direction",
                      "Q", "heterogeneity_p", "I2", "status",
                      "beta_scRT1", "se_scRT1", "p_scRT1", "beta_scRT2", "se_scRT2", "p_scRT2"))
})

test_that("combine_two_samples(failed = 'drop') matches meta_combine of 6.6", {
  o66 <- source_orig_66()
  tab <- make_fits()
  tab[5, c("beta", "se", "p")] <- NA
  tab$status[5] <- "rank_or_n"
  tab$se[20] <- 0
  old <- as.data.frame(o66$meta_combine(data.table::as.data.table(tab), c("gene", "term")))
  new <- combine_two_samples(tab, c("gene", "term"), failed = "drop")
  expect_identical(nrow(new), nrow(old))
  new <- new[match(old$gene, new$gene), ]
  rownames(new) <- NULL
  expect_identical(new[names(old)], old)

  # other sample names, single 'by' column, tibble-free data.table input
  tab2 <- tab
  tab2$sample <- ifelse(tab$sample == "scRT1", "a", "b")
  alt <- combine_two_samples(data.table::as.data.table(tab2), "gene", failed = "drop")
  expect_identical(alt$beta_a, combine_two_samples(tab, "gene", failed = "drop")$beta_scRT1)
  expect_identical(class(alt), "data.frame")

  expect_error(combine_two_samples(tab[-1, ], "gene"), "both samples")
  expect_error(combine_two_samples(rbind(tab, tab[1, ]), "gene"), "once per")
  three <- tab
  three$sample[1] <- "x"
  expect_error(combine_two_samples(three, "gene"), "exactly two samples")
})

test_that("design helpers match the originals", {
  o66 <- source_orig_66()
  o67 <- source_orig_67()
  d <- synth_cells_hc3()
  dt <- data.table::as.data.table(d)

  # 6.6: fixed-knot spline and progress bins; the formulas live in the script
  expect_identical(full_rank_model_matrix(o66$base_formula, d), o66$make_base(dt))
  f_bins <- ~ factor(floor(progress / .05)) + log_rna + log_dna + complexity + call_fraction +
    score_S + score_G2M + score_pluripotency + chek1_local
  expect_identical(full_rank_model_matrix(f_bins, d), o66$make_base(dt, "bins"))
  # one sample only: empty bins and collinear columns are dropped
  one <- d[d$sample == "scRT1" & d$progress < 0.5, ]
  rownames(one) <- NULL
  expect_identical(full_rank_model_matrix(f_bins, one),
                   o66$make_base(data.table::as.data.table(one), "bins"))
  expect_lt(ncol(full_rank_model_matrix(f_bins, one)), ncol(full_rank_model_matrix(f_bins, d)))

  # 6.7: reads the regulator's own replication state from a global matrix
  ii <- which(d$sample == "scRT2")
  o67$own_R <- rbind(Chek1 = d$chek1_local)
  ns <- splines::ns
  f67 <- ~ ns(progress, df = 5) + log_rna + log_dna + complexity + call_fraction +
    score_S + score_G2M + score_pluripotency + own
  d$own <- d$chek1_local
  d_ii <- d[ii, ]
  rownames(d_ii) <- NULL
  expect_identical(full_rank_model_matrix(f67, d_ii),
                   o67$base_matrix(dt[ii], "Chek1", ii))
  bad <- d
  bad$log_rna[2] <- NA
  expect_error(full_rank_model_matrix(f67, bad))

  B <- full_rank_model_matrix(~ log_rna, d)
  expect_identical(add_interaction_terms(B, d$chek1_z, d$local, d$progress_c),
                   o66$gene_design(B, dt, d$local))
  expect_error(add_interaction_terms(B, 1, d$local, d$progress_c), "one value per row")
})

test_that("permute_within_bins reproduces the permutation of 6.6", {
  d <- synth_cells_hc3()
  # the loop body of permute_pairs() in 6.6
  reference <- function(bin_width) {
    groups <- split(seq_len(nrow(d)), interaction(d$sample, floor(d$progress / bin_width), drop = TRUE))
    perm <- seq_len(nrow(d))
    for (gr in groups) if (length(gr) > 1L) perm[gr] <- sample(gr)
    perm
  }
  for (bw in c(0.025, 0.1)) {
    new <- withr::with_seed(7, list(permute_within_bins(d$sample, d$progress, bw),
                                    permute_within_bins(d$sample, d$progress, bw)))
    old <- withr::with_seed(7, list(reference(bw), reference(bw)))
    expect_identical(new, old)
    perm <- new[[1]]
    expect_identical(sort(perm), seq_len(nrow(d)))
    expect_identical(d$sample[perm], d$sample)
    expect_lt(max(abs(d$progress - d$progress[perm])), bw + 1e-10)
    expect_gt(mean(perm != seq_along(perm)), 0.3)
  }
  expect_error(permute_within_bins(d$sample, c(NA, d$progress[-1]), 0.1), "missing")
})

test_that("make_bin and .zscore match the originals", {
  o66 <- source_orig_66()
  o67 <- source_orig_67()
  pos <- c(1L, 20000L, 20001L, 45575327L, 3001L, 158974082L)
  expect_identical(make_bin("chr7", pos), o66$make_bin("chr7", pos))
  # bin ends such as 100000 or 1.1e+08: the original switches to scientific
  # notation unless scipen is raised, and those bins then match nothing
  round_pos <- c(99999L, 109996398L, 24985842L, 2999990L)
  expect_identical(make_bin("chr5", round_pos),
                   withr::with_options(list(scipen = 100), o66$make_bin("chr5", round_pos)))
  expect_identical(make_bin("chr5", 109996398L), "chr5-109980001-110000000")
  expect_match(withr::with_options(list(scipen = 0), o66$make_bin("chr5", 109996398L)), "e\\+")
  expect_identical(make_bin("chr1", c(5, NA)), c("chr1-1-20000", NA))
  expect_identical(make_bin("chr1", 250001, bin_width = 50000), "chr1-250001-300000")

  x <- c(1, 4, 2, NA, 9)
  expect_identical(.zscore(x), o67$zscore(x))
  expect_identical(.zscore(x, 0), o66$zscore(x))
  expect_identical(.zscore(rep(3, 4)), o67$zscore(rep(3, 4)))
  expect_identical(.zscore(rep(3, 4), 0), o66$zscore(rep(3, 4)))
})
