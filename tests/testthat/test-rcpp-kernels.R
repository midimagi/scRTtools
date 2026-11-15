# Smoke tests for the compiled HMM / changepoint kernels, plus an
# equivalence check against the original sourceCpp() builds when the
# original source tree is available.

make_hmm_inputs <- function(n = 300L, k = 10L, seed = 1L) {
  withr::with_seed(seed, {
    counts <- c(sample(0:3, n / 2, TRUE), sample(5:(k - 1), n / 2, TRUE))
    counts[sample.int(n, 10)] <- NA
    emis <- matrix(runif(2 * k), 2, k)
    emis <- emis / rowSums(emis)
  })
  list(
    counts = as.numeric(counts),
    trans = matrix(c(0.99, 0.01, 0.01, 0.99), 2, 2),
    emis = emis,
    init = c(0.5, 0.5)
  )
}

run_kernels <- function(env, x) {
  bw <- env$single_cell_baum_welch_cpp(x$counts, x$trans, x$emis, x$init)
  ex <- env$compute_baum_welch_expectations(
    x$counts, bw$alpha, bw$beta, bw$scalar, x$trans, x$emis
  )
  # initial_states is per position (NaN = unknown), as built in
  # initialize_single_cell_hmm(); it is not a length-2 prior.
  vit <- env$viterbi_decode_cpp(
    x$counts, x$trans, x$emis, rep(NaN, length(x$counts)), integer()
  )
  cps <- env$findchangepts_cpp_max2(
    withr::with_seed(2L, c(rnorm(100), rnorm(100, 3), rnorm(100))), 5L, 2L
  )
  list(bw = bw, ex = ex, vit = vit, cps = cps)
}

pkg_env <- function() {
  ns <- asNamespace("scRTtools")
  list(
    single_cell_baum_welch_cpp = ns$single_cell_baum_welch_cpp,
    compute_baum_welch_expectations = ns$compute_baum_welch_expectations,
    viterbi_decode_cpp = ns$viterbi_decode_cpp,
    findchangepts_cpp_max2 = ns$findchangepts_cpp_max2
  )
}

test_that("compiled kernels run and return expected shapes", {
  x <- make_hmm_inputs()
  res <- run_kernels(pkg_env(), x)
  expect_true(is.finite(res$bw$log_probability_seq))
  expect_length(res$vit, length(x$counts))
  expect_true(all(res$vit[!is.na(res$vit)] %in% c(2, 4)))
  expect_equal(dim(res$ex$trans_exp), c(2L, 2L))
  expect_length(res$cps, 2L)
})

test_that("compiled kernels match the original sourceCpp() builds", {
  src_dir <- file.path(orig_rcode_dir(), "src")
  skip_if_not(dir.exists(src_dir), "original rcode/src not available")
  skip_on_cran()
  orig <- new.env()
  for (f in c("single_cell_baum_welch.cpp", "viterbi_decode.cpp",
              "findchangepts_cpp_2cpts.cpp")) {
    Rcpp::sourceCpp(file.path(src_dir, f), env = orig,
                    cacheDir = file.path(tempdir(), "scRTtools_orig_cpp"))
  }
  x <- make_hmm_inputs()
  expect_identical(run_kernels(pkg_env(), x), run_kernels(as.list(orig), x))
})
