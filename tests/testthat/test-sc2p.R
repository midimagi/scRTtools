test_that("SC2P_gene_unit matches the original", {
  orig <- source_orig_gex()
  d <- synth_gex()
  for (thr in c(15, 4)) {
    new <- suppressWarnings(SC2P_gene_unit(d$counts, yy_thres = thr))
    old <- suppressWarnings(orig$SC2P_gene_unit(d$counts, yy_thres = thr))
    expect_identical(new, old)
  }
  # data.frame input, and a matrix with a single expressed gene
  expect_identical(
    suppressWarnings(SC2P_gene_unit(as.data.frame(d$counts))),
    suppressWarnings(orig$SC2P_gene_unit(as.data.frame(d$counts)))
  )
  one <- d$counts[c(20, 60), , drop = FALSE] + 1L
  one[2, ] <- 0L
  expect_identical(suppressWarnings(SC2P_gene_unit(one)),
                   suppressWarnings(orig$SC2P_gene_unit(one)))
})

test_that("SC2P_gene_unit output shape and edge cases", {
  d <- synth_gex()
  res <- suppressWarnings(SC2P_gene_unit(d$counts))
  expressed <- rownames(d$counts)[rowSums(d$counts) > 0]
  expect_identical(rownames(res), expressed)
  expect_named(res, c("gene", "phase_I_ratio", "phase_II_mean", "phase_II_varcoff"))
  expect_true(all(res$phase_I_ratio >= 0 & res$phase_I_ratio <= 1))
  expect_true(all(is.na(res$phase_II_mean) == (res$phase_I_ratio == 1)))

  expect_error(SC2P_gene_unit(d$counts * 0), "expressed genes")
  bad <- d$counts
  bad[, 1] <- 0L
  expect_error(SC2P_gene_unit(bad), "total count of zero")
  bad <- d$counts
  bad[1, 1] <- NA
  expect_error(SC2P_gene_unit(bad), "missing values")
  expect_error(SC2P_gene_unit(unname(d$counts)), "gene names")
  expect_error(SC2P_gene_unit("a"), "numeric matrix")
})

test_that(".sc2p_robust_poisson matches the original RobustPoi1", {
  orig <- source_orig_gex()
  withr::local_seed(4)
  for (x in list(rpois(200, 1.2), c(rep(0, 50), rpois(100, 3)), integer(0),
                 rep(0, 10), rep(2, 10), c(1, 1, 2, 3, 5))) {
    expect_identical(suppressWarnings(.sc2p_robust_poisson(x)),
                     suppressWarnings(orig$RobustPoi1(x)))
  }
})
