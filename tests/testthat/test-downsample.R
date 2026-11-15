synth_coverage <- function(seed = 7L) {
  withr::with_seed(seed, {
    data.frame(
      state = rep(c("ama", "ctr", "other"), c(300, 400, 50)),
      coverage = c(stats::rlnorm(300, log(80), 0.3), stats::rlnorm(400, log(240), 0.3),
                   stats::rlnorm(50, log(10), 0.3))
    )
  })
}

test_that("optimize_shrink_factor recovers the coverage ratio", {
  md <- synth_coverage()
  res <- suppressWarnings(optimize_shrink_factor(md))
  expect_named(res, c("optimum", "objective", "shrink_factor"))
  expect_equal(res$optimum, 1 / 3, tolerance = 0.05)
  expect_equal(res$shrink_factor, round(res$optimum, 2))
  res2 <- suppressWarnings(optimize_shrink_factor(md, offset = 0.08))
  expect_equal(res2$shrink_factor, round(res$optimum, 2) + 0.08)
  expect_error(optimize_shrink_factor(md, target = "x"), "target")
  expect_error(optimize_shrink_factor(md, coverage_col = "nope"), "nope")
})

test_that("optimize_shrink_factor matches the notebook objective", {
  skip_if_not_installed("dplyr")
  md <- synth_coverage()
  # notebook 1.2.1 code, verbatim apart from namespacing
  optimize_shrink <- function(shrink_factor, meta_data_quantile) {
    meta_tmp <- meta_data_quantile |>
      dplyr::mutate(coverage_scaled = dplyr::case_when(
        state == "ama" ~ coverage,
        state == "ctr" ~ coverage * shrink_factor
      ))
    ama_data <- meta_tmp |> dplyr::filter(state == "ama")
    ctr_data <- meta_tmp |> dplyr::filter(state == "ctr")
    stats::ks.test(ama_data$coverage_scaled, ctr_data$coverage_scaled)$statistic
  }
  old <- suppressWarnings(optimize(function(x) optimize_shrink(x, md),
                                   interval = c(0.1, 1), maximum = FALSE))
  new <- suppressWarnings(optimize_shrink_factor(md, offset = 0.08))
  expect_identical(new$optimum, old$minimum)
  expect_identical(new$objective, unname(old$objective))
  expect_identical(new$shrink_factor, round(old$minimum, 2) + 0.08)
})

test_that("downsample_fragments validates its arguments", {
  f <- tempfile(fileext = ".tsv.gz")
  writeLines("chr1\t1\t100\tBC-1\t1", f)
  expect_error(downsample_fragments("nope.tsv.gz", tempfile(), 0.5), "input")
  expect_error(downsample_fragments(f, f, 0.5), "differ")
  expect_error(downsample_fragments(f, tempfile(), 0), "keep_prob")
  expect_error(downsample_fragments(f, tempfile(), 0.5, seed = -1), "seed")
  expect_error(downsample_fragments(f, file.path(tempdir(), "no", "x.gz"), 0.5),
               "directory")
  expect_error(downsample_fragments(f, tempfile(), 0.5, bgzip = ""), "bgzip")
})

# End-to-end run; needs bgzip/tabix (set SCRTTOOLS_BGZIP / SCRTTOOLS_TABIX or
# have them on PATH) and a C compiler.
find_tool <- function(env, name) {
  p <- Sys.getenv(env, unset = Sys.which(name))
  if (nzchar(p) && file.exists(p)) p else ""
}

test_that("downsample_fragments subsamples reproducibly and matches the original", {
  skip_on_cran()
  bgzip <- find_tool("SCRTTOOLS_BGZIP", "bgzip")
  tabix <- find_tool("SCRTTOOLS_TABIX", "tabix")
  cc <- Sys.which("cc")
  skip_if(!nzchar(bgzip) || !nzchar(tabix) || !nzchar(cc), "bgzip/tabix/cc not available")

  x <- synth_fragments()
  out1 <- tempfile(fileext = ".tsv.gz")
  out2 <- tempfile(fileext = ".tsv.gz")
  run <- function(out, seed) {
    downsample_fragments(x$frag_path, out, keep_prob = 0.3, seed = seed, threads = 1,
                         bgzip = bgzip, tabix = tabix, cc = cc, verbose = FALSE)
  }
  expect_identical(run(out1, 1), out1)
  run(out2, 1)
  expect_true(file.exists(paste0(out1, ".tbi")))
  expect_true(file.exists(paste0(out1, ".downsample_info.txt")))
  l_in <- readLines(gzfile(x$frag_path))
  l1 <- readLines(gzfile(out1))
  expect_identical(l1, readLines(gzfile(out2)))
  expect_true(all(l1 %in% l_in))
  expect_equal(length(l1) / length(l_in), 0.3, tolerance = 0.05)
  expect_error(run(out1, 1), "exit status")          # exists, overwrite = FALSE

  orig <- file.path(local_dir("SCRTTOOLS_NOTEBOOKS_DS"), "downsample_fragments_fast.sh")
  skip_if_not(file.exists(orig), "original script not available")
  out3 <- tempfile(fileext = ".tsv.gz")
  status <- system2("bash", c(shQuote(orig), "-p", "0.3", "-s", "1", "-t", "1",
                              shQuote(x$frag_path), shQuote(out3)),
                    env = c(paste0("BGZIP_BIN=", bgzip), paste0("TABIX_BIN=", tabix),
                            paste0("CC_BIN=", cc)),
                    stdout = FALSE, stderr = FALSE)
  expect_equal(status, 0L)
  expect_identical(l1, readLines(gzfile(out3)))
})
