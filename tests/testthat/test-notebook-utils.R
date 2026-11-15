helpers_file <- test_path("fixtures", "notebook_helpers_2026-11-12.R")

load_notebook_helpers <- function() {
  env <- new.env(parent = globalenv())
  sys.source(helpers_file, envir = env)
  env
}

test_that("safe summaries handle non-finite values", {
  expect_equal(cor_safe(c(1, 2, NA, 4), c(2, 4, 6, Inf)), 1)
  expect_identical(cor_safe(c(1, NA), c(1, 2)), NA_real_)
  z <- zscore_safe(c(1, 2, 3, NA))
  expect_equal(z, c(-1, 0, 1, NA))
  expect_identical(zscore_safe(c(5, 5, NA)), rep(NA_real_, 3))
  expect_identical(zscore_safe(c(1, NA)), c(NA_real_, NA_real_))
  expect_identical(zscore_safe(numeric()), numeric())
  expect_identical(peak_bin(c(0, 25, 50, 75), c(0.1, 0.8, NaN, 0.3)), 25)
  expect_identical(peak_bin(c(0, 25, 50), c(1, 2, 2)), 25)
  expect_identical(peak_bin(c(0, 25), c(NA, NA)), NA_real_)
})

test_that("sample helpers prefix IDs and find exactly one file", {
  expect_identical(add_sample_prefix(c("A-1", "s1-B-1", NA), "s1"),
                   c("s1-A-1", "s1-B-1", NA))
  expect_identical(add_sample_prefix(factor("A"), "s1"), "s1-A")
  expect_error(add_sample_prefix(NULL, "s1"), "NULL for sample: s1")

  d <- withr::local_tempdir()
  dir.create(file.path(d, "a", "s1"), recursive = TRUE)
  dir.create(file.path(d, "merge"))
  dir.create(file.path(d, "s10"))
  file.create(file.path(d, "a", "s1", "m.rds"), file.path(d, "merge", "m.rds"),
              file.path(d, "s10", "m.rds"))
  expect_identical(find_sample_file(d, "s1", "m.rds"),
                   normalizePath(file.path(d, "a", "s1", "m.rds")))
  expect_error(find_sample_file(d, "s2", "m.rds"), "Cannot find m.rds")
  expect_error(find_sample_file(d, "s1", "M.rds"), "Cannot find")
  dir.create(file.path(d, "b", "s1"), recursive = TRUE)
  file.create(file.path(d, "b", "s1", "m.rds"))
  expect_error(find_sample_file(d, "s1", "m.rds"), "Found multiple m.rds")
})

test_that("notebook utilities match notebook_helpers.R", {
  h <- load_notebook_helpers()
  set.seed(1)
  x <- c(rnorm(50), NA, Inf, -Inf, NaN)
  y <- c(rnorm(50), 1, NA, 2, 3)
  pos <- seq(0, 100, length.out = length(x))
  cases <- list(x, y, c(1, NA), c(4, 4, 4), numeric(), c(NA_real_, NaN))
  for (v in cases) {
    expect_identical(zscore_safe(v), h$zscore_safe(v))
    # constant input: both versions warn (zero sd) and return NA
    expect_identical(suppressWarnings(cor_safe(v, rev(v))),
                     suppressWarnings(h$cor_safe(v, rev(v))))
    expect_identical(peak_bin(seq_along(v) * 1.5, v), h$peak_bin(seq_along(v) * 1.5, v))
  }
  expect_identical(cor_safe(x, y), h$cor_safe(x, y))
  expect_identical(peak_bin(pos, y), h$peak_bin(pos, y))

  ids <- list(c("AAAC-1", "s1-AAAG-1", NA), factor(c("x", "s1-y")), character(), 1:3)
  for (v in ids) {
    expect_identical(add_sample_prefix(v, "s1"), h$add_sample_prefix(v, "s1"))
  }
  expect_identical(tryCatch(add_sample_prefix(NULL, "s1"), error = conditionMessage),
                   tryCatch(h$add_sample_prefix(NULL, "s1"), error = conditionMessage))

  d <- withr::local_tempdir()
  dir.create(file.path(d, "x", "s1"), recursive = TRUE)
  dir.create(file.path(d, "merge"))
  file.create(file.path(d, "x", "s1", "m.rds"), file.path(d, "merge", "m.rds"))
  expect_identical(find_sample_file(d, "s1", "m.rds"), h$find_sample_file(d, "s1", "m.rds"))
  for (s in c("s2", "merge2")) {
    expect_identical(tryCatch(find_sample_file(d, s, "m.rds"), error = conditionMessage),
                     tryCatch(h$find_sample_file(d, s, "m.rds"), error = conditionMessage))
  }
  dir.create(file.path(d, "y", "s1"), recursive = TRUE)
  file.create(file.path(d, "y", "s1", "m.rds"))
  expect_identical(tryCatch(find_sample_file(d, "s1", "m.rds"), error = conditionMessage),
                   tryCatch(h$find_sample_file(d, "s1", "m.rds"), error = conditionMessage))
})
