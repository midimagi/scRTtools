test_that("pivot_state_wide and add_ma_columns match the original", {
  orig <- source_orig_fig6()
  d <- synth_state_long()
  for (na in c("drop", "zero", "keep")) {
    new <- pivot_state_wide(d, "iz_id", "usage", na = na)
    old <- orig$pivot_state_wide(d, "iz_id", "usage", na = na)
    expect_identical(new, plain_df(old))
    for (type in c("log2", "relative", "sum")) {
      expect_identical(add_ma_columns(new, type), plain_df(orig$add_ma_columns(old, type)))
    }
  }
  # an integer feature keeps its type when missing values become zero
  di <- d
  di$usage <- as.integer(round(di$usage))
  expect_identical(pivot_state_wide(di, "iz_id", "usage", na = "zero"),
                   plain_df(orig$pivot_state_wide(di, "iz_id", "usage", na = "zero")))
  expect_type(pivot_state_wide(di, "iz_id", "usage", na = "zero")$ctr, "integer")
  # two id columns, tibble input
  d2 <- d
  d2$chr <- sub(":.*", "", d2$iz_id)
  expect_identical(
    pivot_state_wide(dplyr::as_tibble(d2), c("chr", "iz_id"), "width"),
    plain_df(orig$pivot_state_wide(d2, c("chr", "iz_id"), "width"))
  )
  expect_error(pivot_state_wide(rbind(d, d[1, ]), "iz_id", "usage"), "must occur once")
  expect_error(pivot_state_wide(d, "iz_id", "nope"), "lacks column")
  expect_error(pivot_state_wide(d, "iz_id", "usage", states = "ctr"), "two different")
})

test_that("other state names give the same numbers", {
  d <- synth_state_long()
  d2 <- d
  d2$state <- ifelse(d$state == "ctr", "wt", "ko")
  st <- c("wt", "ko")
  ref <- pivot_state_wide(d, "iz_id", "usage", na = "zero")
  alt <- pivot_state_wide(d2, "iz_id", "usage", na = "zero", states = st)
  expect_identical(unname(ref), unname(alt))
  expect_identical(add_ma_columns(ref)$M, add_ma_columns(alt, states = st)$M)
  g_ref <- call_diff_groups(d, "iz_id", "usage", 1)
  g_alt <- call_diff_groups(d2, "iz_id", "usage", 1, states = st)
  expect_identical(levels(g_alt$iz_group), c("ko_up", "ko_down", "no_sig"))
  expect_identical(as.integer(g_ref$iz_group), as.integer(g_alt$iz_group))
})

test_that("add_rep_cluster and combine_rep_tempo_feature match the original", {
  orig <- source_orig_fig6()
  x <- data.frame(rep_percent = c(-1, 0, 10, 33, 33.5, 66, 70, 100, 101, NA), k = 1:10)
  expect_identical(add_rep_cluster(x, c(0, 33, 66, 100)),
                   orig$add_rep_cluster(x, c(0, 33, 66, 100)))
  x$rep_percent <- x$rep_percent / 100
  expect_identical(add_rep_cluster(x, c(0, 0.33, 0.66, 1)),
                   orig$add_rep_cluster(x, c(0, 0.33, 0.66, 1)))
  expect_error(add_rep_cluster(x, c(0, 1)), "four")

  withr::local_seed(5)
  mk <- function(n) {
    t25 <- runif(n, 10, 40)
    df <- data.frame(region_id = sprintf("chr1:%d-%d", seq_len(n) * 1000L, seq_len(n) * 1000L + 500L),
                     t10 = t25 - 5, t25 = t25, t50 = t25 + runif(n, 5, 20),
                     t75 = t25 + runif(n, 25, 40), rmse = runif(n))
    df$t75[2] <- NA
    df$t50[n] <- NA
    df
  }
  a <- mk(20)
  b <- mk(15)
  expect_identical(combine_rep_tempo_feature(a, b),
                   plain_df(orig$combine_rep_tempo_feature(a, b)))
})

test_that("call_diff_groups and summarise_by_diff_group match the original", {
  orig <- source_orig_fig6()
  d <- synth_state_long()
  d$usage[1] <- 0   # log2(0): infinite fold change
  for (thres in c(0.5, 1, 3)) {
    new <- call_diff_groups(d, "iz_id", "usage", thres)
    old <- orig$call_diff_groups(d, "iz_id", "usage", thres)
    expect_identical(new, plain_df(old))
    expect_identical(
      summarise_by_diff_group(d, "iz_id", "width", new),
      plain_df(orig$summarise_by_diff_group(d, "iz_id", "width", old))
    )
  }
  expect_s3_class(new$iz_group, "factor")
  expect_setequal(as.character(unique(new$iz_group)), c("ama_up", "ama_down", "no_sig"))
  expect_false(anyDuplicated(new$iz_id) > 0)
})

test_that("region_id_to_granges and keep_autosomes match the original", {
  orig <- source_orig_fig6()
  d <- data.frame(score = c(1.5, 2, 3, 4),
                  region_id = c("chr1:1-100", "chr19:50-80", "chrX:5-9", "chrM:1-2"),
                  label = letters[1:4], stringsAsFactors = FALSE)
  new <- region_id_to_granges(d)
  expect_identical(new, orig$region_id_to_granges(d))
  expect_identical(names(S4Vectors::mcols(new)), c("score", "region_id", "label"))
  kept <- keep_autosomes(new)
  expect_identical(kept, orig$keep_autosomes(new))
  expect_identical(GenomeInfoDb::seqlevels(kept), paste0("chr", 1:19))
  expect_identical(length(keep_autosomes(new, c("chr1", "chrX"))), 2L)
})

test_that("annotate_track_segments matches the original", {
  skip_if_not_installed("DNAcopy")
  orig <- source_orig_fig6()
  withr::local_seed(8)
  sl <- c(chr1 = 12050000L, chr2 = 8000000L)
  bins <- GenomicRanges::tileGenome(sl, tilewidth = 20000, cut.last.tile.in.chrom = TRUE)
  block <- rep(rep(c(1, -1, 0.6, -0.8), each = 100), length.out = length(bins))
  bins$score <- block + rnorm(length(bins), 0, 0.3)
  regions <- data.frame(
    region_id = c(sprintf("chr1:%d-%d", seq(1L, 11900001L, by = 610000L),
                          seq(1L, 11900001L, by = 610000L) + 99999L),
                  "chr2:450001-650000", "chr2:7900001-8000000", "chr3:1-1000"),
    k = 1, stringsAsFactors = FALSE)

  run <- function(f, track, align) {
    withr::with_seed(99, suppressWarnings(utils::capture.output(
      res <- f(track, regions, "pc1", "compartment", c("A", "B"), align))))
    res
  }
  new <- run(annotate_track_segments, bins, FALSE)
  expect_identical(new, run(orig$annotate_track_segments, bins, FALSE))
  expect_true(all(c("A", "B") %in% new$compartment))
  expect_true(is.na(new$pc1[nrow(new)]))
  expect_identical(new$compartment == "A", new$pc1 > 0)

  # tiles without signal need align_tiles = TRUE
  holes <- bins[-(301:340)]
  new <- run(annotate_track_segments, holes, TRUE)
  expect_identical(new, run(orig$annotate_track_segments, holes, TRUE))
  expect_error(run(annotate_track_segments, holes, FALSE), "align_tiles")
})

test_that("find_intersection and format_bp match the original", {
  orig <- source_orig_fig6()
  withr::local_seed(3)
  x <- seq(2, 98, length.out = 40)
  d <- data.frame(
    state = rep(c("ctr", "ama", "other"), each = 40),
    rep_percent = c(x, rev(x) + 1, x),
    length_velocity = c(1e4 + 300 * x, 3e4 - 150 * (rev(x) + 1) + rnorm(40, 0, 200), x)
  )
  new <- find_intersection(d, "ctr", "ama")
  expect_identical(new, orig$find_intersection(d, "ctr", "ama"))
  expect_true(new$x > 30 && new$x < 60)
  expect_null(find_intersection(d, "ctr", "other"))
  expect_null(orig$find_intersection(d, "ctr", "other"))

  v <- c(0, 999, 1000, 20000, 1499999, 1.5e6, 2.75e9, NA)
  expect_identical(format_bp(v), orig$format_bp(v))
  expect_identical(format_bp(numeric(0)), character(0))
})
