test_that("adjacent IZ pairs, firing matrix and co-firing statistics", {
  f <- iz_analysis_fixture()
  ap <- make_adjacent_pairs(f$diz$izones, f$fr)
  expect_named(ap, c("izones", "features", "pairs"))
  p <- ap$pairs
  expect_true(all(p$chr == ap$features$chr[p$right_index]))
  expect_true(all(p$right_index == p$left_index + 1L))
  fm <- build_firing_matrix(ap$izones, f$diz$tracks, f$x$md)
  expect_s4_class(fm$fired, "dgCMatrix")
  expect_identical(dim(fm$fired), c(length(f$diz$tracks), length(ap$izones)))
  expect_true(all(fm$fired@x == 1))
  fm2 <- build_firing_matrix(ap$izones, f$diz$tracks, f$x$md, max_track_width_bp = 1e5)
  expect_true(all(fm2$cell_info$n_retained_tracks <= fm$cell_info$n_retained_tracks))
  st <- assign_rep_stage(c(0, 0.25, 0.6, 1, NA))
  expect_identical(as.character(st), c("Q1", "Q2", "Q3", "Q4", "unknown"))
  cs <- summarise_cell_adjacency(fm$fired, p, ap$features)
  expect_identical(nrow(cs), nrow(fm$fired))
  expect_equal(cs$n11 + cs$n10 + cs$n01 + cs$n00, rep(nrow(p), nrow(cs)))
  stage <- assign_rep_stage(fm$cell_info$rep_percent)
  bs <- summarise_adjacency_by_stage(fm$fired, p, cbind(fm$cell_info, cs), stage)
  expect_identical(bs$stage_summary$rep_stage, c("all", "Q1", "Q2", "Q3", "Q4"))
  expect_identical(nrow(bs$pair_stats), 5L * nrow(p))
  pi <- make_pair_init(f$fr)
  expect_identical(as.character(pi$pair_id), p$pair_id)
  expect_identical(nrow(make_pair_init(f$fr, n_adj = 1000L)), 0L)
  cb <- summarise_pair_cofiring_bins(pi, f$diz$tracks, f$x$md, bin_num = 20L)
  expect_identical(levels(cb$pair_fired_info$pair_id), levels(pi$pair_id))
  expect_true(all(cb$pair_fired_info$n_all == length(unique(cb$pair_bins$bin_mid))))
})

test_that("adjacency functions match the 4.3.1 notebooks", {
  skip_on_cran()
  f <- iz_analysis_fixture()
  nb <- notebook_functions(
    "4.3.1.analyze_adjacent_izones_features_per_cell.ipynb",
    c("make_adjacent_pairs", "build_firing_matrix", "assign_rep_stage", "safe_ratio",
      "binary_pair_stats", "summarise_cell_adjacency", "summarise_pair_adjacency")
  )
  for (nm in c("mcols", "seqnames", "start", "end", "width", "findOverlaps", "queryHits",
               "subjectHits", "unlist")) {
    assign(nm, getExportedValue(if (nm %in% c("queryHits", "subjectHits", "mcols")) "S4Vectors"
                                else if (nm == "findOverlaps") "GenomicRanges"
                                else if (nm == "seqnames") "GenomeInfoDb" else "BiocGenerics", nm),
           envir = nb)
  }
  nb$`ranges<-` <- IRanges::`ranges<-`
  nb$IRanges <- IRanges::IRanges
  nb$rep_breaks <- c(0, 0.25, 0.5, 0.75, 1)
  nb$rep_labels <- c("Q1", "Q2", "Q3", "Q4")
  old_ap <- nb$make_adjacent_pairs(f$diz$izones, f$fr)
  ap <- make_adjacent_pairs(f$diz$izones, f$fr)
  expect_identical(ap, old_ap)
  for (w in list(NULL, 1e5)) {
    expect_identical(build_firing_matrix(ap$izones, f$diz$tracks, f$x$md, w),
                     nb$build_firing_matrix(ap$izones, f$diz$tracks, f$x$md, w))
  }
  fm <- build_firing_matrix(ap$izones, f$diz$tracks, f$x$md)
  rp <- c(fm$cell_info$rep_percent, NA, 1)
  expect_identical(assign_rep_stage(rp), nb$assign_rep_stage(rp))
  expect_identical(summarise_cell_adjacency(fm$fired, ap$pairs, ap$features),
                   nb$summarise_cell_adjacency(fm$fired, ap$pairs, ap$features))
  sel <- which(fm$cell_info$rep_percent > 0.4)
  expect_identical(summarise_pair_adjacency(fm$fired, ap$pairs, sel, "x"),
                   nb$summarise_pair_adjacency(fm$fired, ap$pairs, sel, "x"))

  # binned pair co-firing: notebook 4.3.1.analyze_adjacent_izones (dplyr code)
  skip_if_not_installed("dplyr")
  nb1 <- notebook_functions("4.3.1.analyze_adjacent_izones.ipynb", "make_pair_init",
                            parent = asNamespace("dplyr"))
  nb1$seqnames <- GenomeInfoDb::seqnames
  nb1$start <- BiocGenerics::start
  nb1$end <- BiocGenerics::end
  pair_raw <- nb1$make_pair_init(f$fr, n_adj = 0)
  expect_identical(make_pair_init(f$fr), pair_raw)
  `%>%` <- dplyr::`%>%`
  md <- f$x$md
  md$rep_percent[c(3, 9)] <- NA
  pair_summary_list <- list()
  for (cell_id in names(f$diz$tracks)) {
    putative_iz <- quiet_deprecated(
      S4Vectors::mcols(f$diz$tracks[[cell_id]]) %>% as.data.frame() %>%
        dplyr::arrange(iz_id) %>% dplyr::pull(iz_id) %>% stats::na.omit())
    pair_info <- pair_raw %>% dplyr::mutate(
      left_iz_state = ifelse(left_iz_id %in% putative_iz, 1, 0),
      right_iz_state = ifelse(right_iz_id %in% putative_iz, 1, 0)
    ) %>% dplyr::mutate(cell = cell_id, .after = pair_id)
    pair_summary_list[[cell_id]] <- pair_info
  }
  pair_summary_all <- dplyr::bind_rows(pair_summary_list)
  bin_num <- 25
  bin_breaks <- seq(0, 1, length.out = bin_num + 1)
  rep_percent_map <- md %>% dplyr::select(cell_id, rep_percent) %>% stats::na.omit() %>%
    dplyr::arrange(rep_percent) %>%
    dplyr::mutate(rep_bin = cut(rep_percent, breaks = bin_breaks, labels = FALSE,
                                include.lowest = TRUE),
                  bin_mid = (rep_bin - 0.5) / bin_num)
  pair_summary_all_bins <- pair_summary_all %>%
    dplyr::left_join(rep_percent_map, by = dplyr::join_by(cell == cell_id)) %>%
    dplyr::group_by(pair_id, bin_mid, left_iz_id, right_iz_id) %>%
    dplyr::summarise(rep_percent = stats::median(rep_percent),
                     left_iz_state = sum(left_iz_state),
                     right_iz_state = sum(right_iz_state), .groups = "drop") %>%
    dplyr::relocate(rep_percent, .after = bin_mid)
  pair_fired_info <- pair_summary_all_bins %>% dplyr::group_by(pair_id) %>%
    dplyr::summarise(
      n_both_no_fired = sum(left_iz_state == 0 & right_iz_state == 0),
      n_left_fired = sum(left_iz_state != 0 & right_iz_state == 0),
      n_right_fired = sum(left_iz_state == 0 & right_iz_state != 0),
      n_both_fired = sum(left_iz_state != 0 & right_iz_state != 0)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(n_all = (n_both_no_fired + n_left_fired + n_right_fired + n_both_fired),
                  both_no_fired_rate = n_both_no_fired / n_all,
                  left_fired_rate = n_left_fired / n_all,
                  right_fired_rate = n_right_fired / n_all,
                  both_fired_rate = n_both_fired / n_all)
  new <- summarise_pair_cofiring_bins(pair_raw, f$diz$tracks, md, bin_num = bin_num)
  expect_identical(new$pair_bins, as.data.frame(pair_summary_all_bins))
  expect_identical(new$pair_fired_info, as.data.frame(pair_fired_info))
})
