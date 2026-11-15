# Regression on real data: the same chain of steps is run on the example
# dataset with either the package functions or the original scripts.

regression_ref_file <- function() {
  file.path(local_dir("SCRTTOOLS_REFDATA"), "regression_scRT1_genome.rds")
}

regression_infer_args <- function() {
  # values of the main-pipeline notebooks
  list(hmm_aggr_width = 15L, init_ratio_lower = 1.2, final_ratio_lower = 1.4)
}

regression_impl_package <- function() {
  list(
    mapd = function(m, gr, md) calculate_mapd(m, gr, md, verbose = FALSE),
    g1 = function(md) find_g1_cells(md, seed = 1, verbose = FALSE),
    states = function(wc, g1w, bins) {
      assign_replication_states(wc, g1w, bins, infer_args = regression_infer_args(),
                                verbose = FALSE)
    },
    rep_percent = function(gs, bins, md) compute_rep_percent(gs, bins, md, verbose = FALSE),
    tracks = function(gs, bins, rp) define_replication_tracks(gs, bins, rp, verbose = FALSE),
    lengths = calc_rep_track_lengths,
    gaps = calc_rep_track_gaps,
    pseudobulk = function(m, gr, md) suppressWarnings(calculate_pseudobulk_rt(m, gr, md)),
    peaks = function(rt) find_rt_peaks(rt),
    izones = function(tr, pk) define_initiation_zones(tr, pk, min_usage = 2L, verbose = FALSE),
    firing_range = function(iz, gs, tr) calculate_firing_range(iz, gs, tr, verbose = FALSE),
    firing_order = function(iz, gs, md) calculate_firing_order(iz, gs, md),
    predict_order = function(iz, tr) predict_firing_order(iz, tr),
    rep_tempo = function(gs, regions, md) {
      rep_tempo(gs, regions, meta_data = md, nbin_rep = regression_nbin_rep())
    },
    fit_tempo = function(x) {
      lapply(regression_tempo_models(), function(model) {
        fit_rep_tempo(x, model_type = model, verbose = FALSE)
      })
    },
    sc2p = function(y) suppressWarnings(SC2P_gene_unit(y)),
    gex_tempo = function(counts, md) {
      lapply(regression_gex_methods(), function(method) {
        suppressWarnings(gex_tempo(counts, md, nbin_rep = regression_nbin_rep(),
                                   method = method))
      })
    },
    markers = function(counts, meta) {
      cfg <- regression_marker_config()
      meta <- prepare_metadata(meta, colnames(counts))
      norm <- normalize_sparse_counts(counts)
      res <- analyze_one_group(
        "ctr", meta, norm$log_normalized, norm$normalized,
        min_cells_per_group = cfg$min_cells_per_group,
        min_cell_fraction = cfg$min_cell_fraction, fdr_threshold = cfg$fdr_threshold,
        gam_basis_k = cfg$gam_basis_k, n_curve_points = cfg$n_curve_points, verbose = FALSE)
      list(meta_data = meta, library_size = norm$library_size, result = res)
    },
    compare = list(
      pivot = function(df, id, feature, na) pivot_state_wide(df, id, feature, na = na),
      ma = add_ma_columns,
      rep_cluster = add_rep_cluster,
      tempo_feature = combine_rep_tempo_feature,
      to_granges = region_id_to_granges,
      autosomes = keep_autosomes,
      segments = function(track, df) regression_segments(annotate_track_segments, track, df),
      diff_groups = call_diff_groups,
      by_diff_group = summarise_by_diff_group,
      intersection = find_intersection,
      format_bp = format_bp
    ),
    hc3 = list(
      make_bin = make_bin,
      zscore = .zscore,
      design = function(d) full_rank_model_matrix(regression_hc3_formula(), d),
      interactions = add_interaction_terms,
      ols_exposure = function(B, x, y) fit_ols_hc3(cbind(B, exposure = x), y),
      ols_terms = function(X, y, terms) fit_ols_hc3(X, y, terms, max_leverage = Inf),
      poisson = function(X, y, offset, term) fit_poisson_hc3(X, y, offset, term),
      combine_keep = function(tab) {
        combine_two_samples(tab, "term")[regression_combine_cols("keep", tab$sample)]
      },
      combine_drop = function(tab) {
        combine_two_samples(tab, "term", failed = "drop")[
          regression_combine_cols("drop", tab$sample)]
      },
      permute = permute_within_bins
    )
  )
}

regression_nbin_rep <- function() 5L

regression_tempo_models <- function() {
  stats::setNames(nm = c("sslogis", "gam_low", "gam_fast"))
}

regression_gex_methods <- function() stats::setNames(nm = c("sc2p", "cp10k"))

# Settings of the progression-marker script, with the group size, basis
# dimension and marker thresholds loosened to what 28 S-phase cells allow
# (chosen on the earlier chr19-only example, where at the script's 0.20 / 0.01
# no gene was called and the curve and clustering code would not run).
regression_marker_config <- function() {
  list(group_by_col = "state", use_sample_covariate = FALSE,
       min_cells_per_group = 10L, min_normalized_expression = 0.1,
       min_cell_fraction = 0.05, fdr_threshold = 0.5, aic_rule = "gam_lower",
       gam_basis_k = 5L, n_curve_points = 41L, n_clusters = 3L, workers = 1L)
}

# The example holds one sample and one condition. For the two-state and
# two-sample functions its S-phase cells are split into two pseudo-groups,
# alternating along rep_percent, so that both groups span S phase. The split
# has no biological meaning; it only gives both implementations the same
# two-group input.
regression_pseudo_groups <- function(md, labels) {
  group <- character(nrow(md))
  group[order(md$rep_percent, rownames(md))] <- rep(labels, length.out = nrow(md))
  stats::setNames(group, rownames(md))
}

regression_region_id <- function(gr) {
  sprintf("%s:%d-%d", as.character(GenomicRanges::seqnames(gr)),
          as.integer(BiocGenerics::start(gr)), as.integer(BiocGenerics::end(gr)))
}

# DNAcopy draws random numbers and prints its progress.
regression_segments <- function(f, track, df) {
  withr::with_seed(99, suppressWarnings(utils::capture.output(
    res <- f(track, df, "rt_segment", "rt_class", c("early", "late"), TRUE))))
  res
}

regression_hc3_formula <- function() ~ splines::ns(pos, df = 4) + gc

# Columns shared by combine_two_samples() and the 6.7 (`keep`) or 6.6
# (`drop`) function it replaces.
regression_combine_cols <- function(failed, sample) {
  if (failed == "keep") {
    return(c("term", "beta", "se", "p", "replication_p", "same_direction", "Q",
             "heterogeneity_p", "status"))
  }
  per_sample <- as.vector(outer(c("beta_", "se_", "p_"), sort(unique(sample)), paste0))
  c("term", "beta", "se", "p", per_sample, "same_direction", "replication_p",
    "heterogeneity_p", "I2")
}

regression_impl_original <- function() {
  serial <- BiocParallel::SerialParam()
  o_mapd <- source_orig("calculate_mapd.R")
  o_g1 <- source_orig("find_g1_cells.R")
  o_hmm <- source_orig_hmm(with_assign = TRUE)
  o_rp <- source_orig("compute_rep_percent.R")
  o_tr <- source_orig_file("processing/define_replication_tracks.R")
  o_len <- source_orig_file("analysis/calc_rep_track_lengths.R")
  o_gap <- source_orig_file("analysis/calc_rep_track_gaps.R")
  o_pb <- source_orig_file("processing/calculate_pseudobulk_rt.R", dplyr = TRUE)
  o_pk <- source_orig_file("processing/find_rt_peaks.R")
  o_iz <- source_orig_izones()
  o_tempo <- source_orig_tempo("tempo.R")
  o_fit <- source_orig_tempo("tempo_models1.R")
  o_gex <- source_orig_gex()
  o_mk <- source_orig_markers()
  o_fig6 <- source_orig_fig6()
  o66 <- source_orig_66()
  o67 <- source_orig_67()
  o66$base_formula <- regression_hc3_formula()
  list(
    mapd = function(m, gr, md) suppressMessages(o_mapd$calculate_mapd(m, gr, md, bpparam = serial)),
    g1 = function(md) withr::with_seed(1, suppressMessages(o_g1$find_g1_cells(md))),
    states = function(wc, g1w, bins) {
      o_hmm$assign_replication_states(wc, g1w, bins, bpparam = serial,
                                      infer_args = regression_infer_args(), verbose = FALSE)
    },
    rep_percent = function(gs, bins, md) {
      o_rp$compute_rep_percent(gs, bins, md, bpparam = serial, verbose = FALSE)
    },
    tracks = function(gs, bins, rp) {
      o_tr$define_replication_tracks(gs, bins, rp, bpparam = serial, verbose = FALSE)
    },
    lengths = o_len$calc_rep_track_lengths,
    gaps = o_gap$calc_rep_track_gaps,
    pseudobulk = function(m, gr, md) suppressWarnings(o_pb$calculate_pseudobulk_rt(m, gr, md)),
    peaks = function(rt) o_pk$find_rt_peaks(rt),
    izones = function(tr, pk) {
      o_iz$define_initiation_zones(tr, pk, min_usage = 2L, verbose = FALSE, bpparam = serial)
    },
    firing_range = function(iz, gs, tr) o_iz$calculate_firing_range(iz, gs, tr, verbose = FALSE),
    firing_order = function(iz, gs, md) o_iz$calculate_firing_order(iz, gs, md),
    predict_order = function(iz, tr) o_iz$predict_firing_order(iz, tr),
    rep_tempo = function(gs, regions, md) {
      o_tempo$rep_tempo(gs, regions, meta_data = md, nbin_rep = regression_nbin_rep())
    },
    fit_tempo = function(x) {
      lapply(regression_tempo_models(), function(model) {
        o_fit$fit_rep_tempo(x, model_type = model, n_boot = 0L, verbose = FALSE)
      })
    },
    sc2p = function(y) suppressWarnings(o_gex$SC2P_gene_unit(y)),
    gex_tempo = function(counts, md) {
      lapply(regression_gex_methods(), function(method) {
        suppressWarnings(o_gex$gex_tempo(counts, md, nbin_rep = regression_nbin_rep(),
                                         method = method))
      })
    },
    markers = function(counts, meta) {
      cfg <- regression_marker_config()
      meta_file <- tempfile(fileext = ".rds")
      on.exit(unlink(meta_file))
      saveRDS(meta, meta_file)
      meta <- o_mk$prepare_metadata(meta_file, colnames(counts), "state", 0.01, 0.99)
      norm <- o_mk$normalize_sparse_counts(counts)
      res <- suppressMessages(o_mk$analyze_one_group(
        "ctr", meta, norm$log_normalized, norm$normalized, cfg))
      list(meta_data = meta, library_size = norm$library_size, result = res)
    },
    # the originals return tibbles where the package returns data.frames
    compare = list(
      pivot = function(df, id, feature, na) {
        plain_df(o_fig6$pivot_state_wide(df, id, feature, na = na))
      },
      ma = function(df, type) plain_df(o_fig6$add_ma_columns(df, type)),
      rep_cluster = o_fig6$add_rep_cluster,
      tempo_feature = function(a, b) plain_df(o_fig6$combine_rep_tempo_feature(a, b)),
      to_granges = o_fig6$region_id_to_granges,
      autosomes = o_fig6$keep_autosomes,
      segments = function(track, df) {
        regression_segments(o_fig6$annotate_track_segments, track, df)
      },
      diff_groups = function(df, id, feature, thres) {
        plain_df(o_fig6$call_diff_groups(df, id, feature, thres))
      },
      by_diff_group = function(df, id, value, groups) {
        plain_df(o_fig6$summarise_by_diff_group(df, id, value, groups))
      },
      intersection = o_fig6$find_intersection,
      format_bp = o_fig6$format_bp
    ),
    hc3 = list(
      make_bin = o66$make_bin,
      zscore = o67$zscore,
      design = function(d) o66$make_base(data.table::as.data.table(d)),
      interactions = function(B, regulator, local_rep, progress_c) {
        o66$gene_design(B, list(progress_c = progress_c), local_rep, z = regulator)
      },
      ols_exposure = function(B, x, y) {
        # fit_hc3() leaves the coefficient name on `p`
        data.frame(term = "exposure", lapply(o67$fit_hc3(B, x, y), unname),
                   stringsAsFactors = FALSE)
      },
      ols_terms = function(X, y, terms) as.data.frame(o66$ols_hc3(X, y, terms)),
      poisson = function(X, y, offset, term) {
        data.frame(term = term, as.data.frame(o66$count_hc3(X, y, offset, term)),
                   stringsAsFactors = FALSE)
      },
      combine_keep = function(tab) {
        dt <- data.table::as.data.table(tab)
        old <- as.data.frame(dt[, o67$combine_two(.SD), by = "term"])
        old[regression_combine_cols("keep", tab$sample)]
      },
      combine_drop = function(tab) {
        old <- as.data.frame(o66$meta_combine(data.table::as.data.table(tab), "term"))
        # 6.6 writes its own sample names into the column names
        smp <- sort(unique(tab$sample))
        names(old) <- sub("_scRT1$", paste0("_", smp[1]), names(old))
        names(old) <- sub("_scRT2$", paste0("_", smp[2]), names(old))
        old[regression_combine_cols("drop", tab$sample)]
      },
      # the loop body of permute_pairs() in 6.6, which itself reads globals
      permute = function(group, progress, bin_width) {
        groups <- split(seq_along(progress),
                        interaction(group, floor(progress / bin_width), drop = TRUE))
        perm <- seq_along(progress)
        for (gr in groups) if (length(gr) > 1L) perm[gr] <- sample(gr)
        perm
      }
    )
  )
}

# Two-state comparison of the pseudo-groups: RT peaks (usage = tracks of the
# group over the peak, rt = replicated fraction at the peak centre), tempo
# features of the initiation zones and track length along S phase.
regression_group_compare <- function(impl, ex, states, md, out) {
  fn <- impl$compare
  res <- list()
  state_levels <- c("ctr", "ama")
  state <- regression_pseudo_groups(md, state_levels)
  cells_of <- function(s) names(state)[state == s]

  peaks <- out$rt_peaks
  region_id <- regression_region_id(peaks)
  bin <- GenomicRanges::findOverlaps(GenomicRanges::resize(peaks, 1L, fix = "center"),
                                     ex$bins, select = "first")
  tracks <- BiocGenerics::unlist(out$izones$tracks, use.names = FALSE)
  long <- do.call(rbind, lapply(state_levels, function(s) {
    usage <- GenomicRanges::countOverlaps(peaks, tracks[tracks$cell_id %in% cells_of(s)])
    data.frame(region_id = region_id, state = s,
               usage = ifelse(usage > 0, usage, NA),
               rt = rowMeans(states[bin, cells_of(s), drop = FALSE] == 4, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  long$state <- factor(long$state, levels = state_levels)

  res$wide <- lapply(stats::setNames(nm = c("drop", "zero", "keep")), function(na) {
    fn$pivot(long, "region_id", "usage", na)
  })
  res$ma <- lapply(stats::setNames(nm = c("log2", "relative", "sum")), function(type) {
    fn$ma(res$wide$drop, type)
  })
  res$diff_groups <- fn$diff_groups(long, "region_id", "usage", 0.5)
  res$by_diff_group <- fn$by_diff_group(long, "region_id", "rt", res$diff_groups)
  res$rep_cluster <- fn$rep_cluster(
    data.frame(cell_id = rownames(md), state = unname(state), rep_percent = md$rep_percent,
               stringsAsFactors = FALSE),
    c(0, 0.33, 0.66, 1))

  # sslogis tempo fit of the initiation zones within each group
  iz <- out$izones$izones
  iz$rt <- 0
  fits <- lapply(stats::setNames(nm = state_levels), function(s) {
    x <- impl$rep_tempo(states[, cells_of(s), drop = FALSE], iz, md[cells_of(s), , drop = FALSE])
    fit <- impl$fit_tempo(x)$sslogis
    stopifnot(nrow(fit) == length(iz))
    data.frame(region_id = regression_region_id(iz), t25 = fit$t25, t50 = fit$t50,
               t75 = fit$t75, stringsAsFactors = FALSE)
  })
  res$tempo_fits <- fits
  res$tempo_feature <- fn$tempo_feature(fits$ctr, fits$ama)
  res$tempo_granges <- fn$autosomes(fn$to_granges(res$tempo_feature))

  # early / late segments of the pseudobulk RT profile under each peak
  res$segments <- fn$segments(out$pseudobulk_rt, res$wide$keep)

  # median track length of every cell along S phase, and where the groups cross
  len <- out$track_lengths
  med <- tapply(len$width_bp, len$cell_id, stats::median)
  med <- med[intersect(rownames(md), names(med))]
  velocity <- data.frame(cell_id = names(med), state = unname(state[names(med)]),
                         rep_percent = md[names(med), "rep_percent"],
                         length_velocity = as.numeric(med), stringsAsFactors = FALSE)
  res$velocity <- velocity
  res$intersection <- fn$intersection(velocity, "ctr", "ama")
  res$length_labels <- fn$format_bp(velocity$length_velocity)
  res
}

# HC3 regressions with genes as observations, fitted in each pseudo-group and
# combined: replicated fraction at the gene against its expression (6.7
# style), the same with peak membership and its interaction (6.6 style), and
# a Poisson model of the counts.
regression_hc3 <- function(impl, ex, states, md, out) {
  fn <- impl$hc3
  res <- list()
  genes <- ex$genes
  sample <- regression_pseudo_groups(md, c("half1", "half2"))
  tss <- GenomicRanges::GRanges(genes$chr, IRanges::IRanges(genes$start, width = 1L))
  bin <- GenomicRanges::findOverlaps(tss, ex$bins, select = "first")
  in_peak <- GenomicRanges::countOverlaps(tss, out$rt_peaks, maxgap = 100000L) > 0
  rna <- ex$rna_counts[genes$gene, rownames(md), drop = FALSE]

  fits <- lapply(stats::setNames(nm = sort(unique(sample))), function(s) {
    cells <- names(sample)[sample == s]
    d <- data.frame(
      gene = genes$gene, pos = genes$start / 1e6,
      gc = ex$bins$gc_frequency[bin],
      rt = rowMeans(states[bin, cells, drop = FALSE] == 4, na.rm = TRUE),
      count = Matrix::rowSums(rna[, cells, drop = FALSE]),
      in_peak = as.numeric(in_peak), stringsAsFactors = FALSE)
    d <- d[!is.na(bin) & ex$bins$is_mappable[bin] %in% TRUE & is.finite(d$rt), ]
    rownames(d) <- NULL
    d$expr_z <- fn$zscore(log1p(d$count / sum(d$count) * 1e4))
    d$pos_c <- d$pos - mean(d$pos)
    B <- fn$design(d)
    X <- fn$interactions(B, d$expr_z, d$in_peak, d$pos_c)
    duplicate <- fn$ols_exposure(B, B[, 2], d$rt)
    duplicate$term <- "duplicate"
    poisson <- fn$poisson(cbind(B, rt = d$rt), d$count, rep(log(sum(d$count)), nrow(d)), "rt")
    poisson$term <- "poisson_rt"
    list(n_genes = nrow(d), design = B, interactions = X,
         ols_exposure = fn$ols_exposure(B, d$expr_z, d$rt),
         ols_terms = fn$ols_terms(X, d$rt, c("local_rep", "interaction")),
         duplicate = duplicate, poisson = poisson)
  })
  res$fits <- fits

  cols <- c("term", "beta", "se", "p", "n", "status")
  tab <- do.call(rbind, lapply(names(fits), function(s) {
    f <- fits[[s]]
    one <- rbind(f$ols_exposure[cols], f$ols_terms[cols], f$duplicate[cols], f$poisson[cols])
    one$sample <- s
    one
  }))
  rownames(tab) <- NULL
  res$fit_table <- tab
  by_term <- function(x) {
    x <- x[order(x$term), , drop = FALSE]
    rownames(x) <- NULL
    x
  }
  res$combined_keep <- by_term(fn$combine_keep(tab))
  res$combined_drop <- by_term(fn$combine_drop(tab))
  res$permutation <- withr::with_seed(7, fn$permute(unname(sample), md$rep_percent, 0.1))
  res
}

# Cell QC -> states -> tracks -> pseudobulk RT -> initiation zones ->
# replication tempo -> RNA -> group comparison -> HC3 statistics on the
# example data. Every step feeds on the previous results of the same
# implementation.
regression_pipeline <- function(impl, ex = scrt_example_data()) {
  out <- list()
  cells <- rownames(ex$meta)
  s_cells <- cells[!ex$meta$is_G1]
  bin_counts <- Matrix::Matrix(ex$bin_counts, sparse = TRUE)
  window_counts <- Matrix::Matrix(ex$window_counts[, s_cells], sparse = TRUE)

  out$mapd <- impl$mapd(bin_counts, ex$bins, ex$meta[, c("cell_id", "sample")])
  out$g1 <- impl$g1(out$mapd)

  out$states <- impl$states(window_counts, ex$g1_windows, ex$bins)
  out$rep_percent <- impl$rep_percent(out$states$genobins_rep_states, ex$bins,
                                      ex$meta[s_cells, c("cell_id", "sample")])
  states <- out$rep_percent$genobins_rep_states_filtered
  md <- out$rep_percent$meta_data
  keep <- intersect(colnames(states), rownames(md)[is.finite(md$rep_percent)])
  states <- states[, keep, drop = FALSE]
  md <- md[keep, , drop = FALSE]

  out$tracks <- impl$tracks(states, ex$bins, stats::setNames(md$rep_percent, keep))
  out$track_lengths <- impl$lengths(out$tracks)
  out$track_gaps <- impl$gaps(out$tracks)

  out$pseudobulk_rt <- impl$pseudobulk(bin_counts, ex$bins,
                                       ex$meta[, c("cell_id", "sample", "is_G1", "pass_qc")])
  out$rt_peaks <- impl$peaks(out$pseudobulk_rt)

  out$izones <- impl$izones(out$tracks, out$rt_peaks)
  iz <- out$izones$izones
  # fraction of cells replicated at the IZ centre, as the notebooks attach it
  centre <- GenomicRanges::resize(iz, 1L, fix = "center")
  bin <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(centre, ex$bins))
  iz$rt <- rowMeans(states[bin, , drop = FALSE] == 4, na.rm = TRUE)
  with_tracks <- intersect(keep, unique(BiocGenerics::unlist(out$izones$tracks,
                                                             use.names = FALSE)$cell_id))
  out$firing_range <- impl$firing_range(iz, states[, with_tracks, drop = FALSE],
                                        out$izones$tracks)
  # genome-wide a few IZ centres have no state call in any cell (rt = NaN);
  # as in notebook 4.3 those IZs are dropped for the two ordering steps (the
  # tracks' iz_id index the full IZ set, so firing_range keeps all IZs)
  iz_rt <- iz[is.finite(iz$rt)]
  out$firing_order <- impl$firing_order(iz_rt, states, md)
  out$predicted_order <- impl$predict_order(iz_rt, out$izones$tracks)

  # replication tempo of the initiation zones along S-phase progress
  out$rep_tempo <- impl$rep_tempo(states, iz, md)
  out$tempo_fit <- impl$fit_tempo(out$rep_tempo)

  # RNA of the same S-phase cells; genes seen in fewer than 3 of them carry
  # no information at this size
  rna <- ex$rna_counts[, keep, drop = FALSE]
  rna <- rna[Matrix::rowSums(rna > 0) >= 3, , drop = FALSE]
  out$sc2p <- impl$sc2p(as.matrix(rna))
  out$gex_tempo <- impl$gex_tempo(rna, md)

  # S-phase progression markers on all cells with RNA, using the recomputed
  # rep_percent. The original derives `state` from hard-coded sample-name
  # patterns that do not fit this sample, so the column is supplied.
  meta <- ex$meta[, c("cell_id", "sample", "is_G1", "rep_percent", "pass_qc")]
  meta[keep, "rep_percent"] <- md$rep_percent
  meta$state <- "ctr"
  out$markers <- impl$markers(ex$rna_counts, meta)

  # two-group comparison and HC3 statistics on two pseudo-groups of the S cells
  out$group_compare <- regression_group_compare(impl, ex, states, md, out)
  # 20 kb bin of every gene start; the 6.6 original writes round bin ends in
  # scientific notation, so this step is compared separately in the test
  out$gene_bins <- impl$hc3$make_bin(ex$genes$chr, ex$genes$start)
  out$hc3 <- regression_hc3(impl, ex, states, md, out)
  out
}
