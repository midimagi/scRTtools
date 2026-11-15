# Two-state (e.g. control vs treatment) comparison helpers
# (non-plotting functions of retreat_figure6_helpers.R, dplyr/tidyr rewritten
# in base R; the hard-coded "ctr"/"ama" became the `states` argument).

#' @keywords internal
#' @noRd
.check_states <- function(states) {
  if (!is.character(states) || length(states) != 2L || anyNA(states) ||
      states[1] == states[2]) {
    stop("'states' must be two different state names: c(reference, treatment).")
  }
  invisible(states)
}

#' @keywords internal
#' @noRd
.check_columns <- function(df, cols, what = "df") {
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols)) {
    stop("'", what, "' lacks column(s): ", paste(missing_cols, collapse = ", "))
  }
  invisible(TRUE)
}

#' Drop rows holding any NA and renumber the rows
#' @keywords internal
#' @noRd
.drop_na_rows <- function(df) {
  df <- df[stats::complete.cases(df), , drop = FALSE]
  rownames(df) <- NULL
  df
}

#' Long table to one column per state
#'
#' Spreads `feature` of a long table (one row per `id_col` and `state`) into
#' one column per state.
#'
#' @param df data.frame with the `id_col` column(s), a `state` column and
#'   `feature`.
#' @param id_col Name(s) of the identifier column(s).
#' @param feature Name of the value column.
#' @param na `"drop"`: remove rows with any `NA`; `"zero"`: remove rows where
#'   both states are `NA` and set the remaining `NA` to 0; `"keep"`: leave
#'   as is.
#' @param states The two state names, `c(reference, treatment)`.
#'
#' @return A data.frame with `id_col` and one column per state found in
#'   `df$state` (in order of first appearance), one row per identifier.
#'
#' @examples
#' long <- data.frame(id = rep(c("a", "b", "c"), 2),
#'                    state = rep(c("ctr", "ama"), each = 3),
#'                    usage = c(1, 2, NA, 2, 2, 5))
#' pivot_state_wide(long, "id", "usage")
#' pivot_state_wide(long, "id", "usage", na = "zero")
#' @export
pivot_state_wide <- function(df, id_col, feature, na = c("drop", "zero", "keep"),
                             states = c("ctr", "ama")) {
  na <- match.arg(na)
  .check_states(states)
  if (length(feature) != 1L) stop("'feature' must be one column name.")
  .check_columns(df, c(id_col, "state", feature))
  df <- as.data.frame(df)

  key <- do.call(paste, c(unname(as.list(df[id_col])), sep = "\r"))
  state <- as.character(df$state)
  if (anyDuplicated(paste(key, state, sep = "\r"))) {
    stop("Each combination of '", paste(id_col, collapse = "', '"),
         "' and state must occur once.")
  }
  first <- !duplicated(key)
  wide <- df[first, id_col, drop = FALSE]
  rownames(wide) <- NULL
  wide_key <- key[first]
  for (s in unique(state)) {
    in_state <- which(state == s)
    wide[[s]] <- df[[feature]][in_state[match(wide_key, key[in_state])]]
  }

  if (na == "drop") {
    wide <- .drop_na_rows(wide)
  } else if (na == "zero") {
    .check_columns(wide, states, "the wide table")
    wide <- wide[!(is.na(wide[[states[1]]]) & is.na(wide[[states[2]]])), , drop = FALSE]
    for (s in states) {
      # an integer feature stays integer, as with tidyr::replace_na()
      zero <- if (is.integer(wide[[s]])) 0L else 0
      wide[[s]][is.na(wide[[s]])] <- zero
    }
    rownames(wide) <- NULL
  }
  wide
}

#' Add A / M columns comparing two states
#'
#' @param df Wide data.frame with one column per state, e.g. from
#'   [pivot_state_wide()].
#' @param type `"log2"`: `A` = mean of the log2 values, `M` = log2 fold
#'   change (treatment / reference); `"relative"`: `A` = reference,
#'   `M = (treatment - reference) / reference`; `"sum"`: `A` = sum,
#'   `M = treatment - reference`.
#' @inheritParams pivot_state_wide
#' @return `df` with columns `A` and `M` added.
#' @examples
#' add_ma_columns(data.frame(id = 1:3, ctr = c(1, 2, 4), ama = c(2, 2, 1)))
#' @export
add_ma_columns <- function(df, type = c("log2", "relative", "sum"),
                           states = c("ctr", "ama")) {
  type <- match.arg(type)
  .check_states(states)
  .check_columns(df, states)
  ref <- df[[states[1]]]
  trt <- df[[states[2]]]
  if (type == "log2") {
    df$A <- (log2(ref) + log2(trt)) / 2
    df$M <- log2(trt) - log2(ref)
  } else if (type == "relative") {
    df$A <- ref
    df$M <- (trt - ref) / ref
  } else {
    df$A <- ref + trt
    df$M <- trt - ref
  }
  df
}

#' Label cells or regions as early / middle / late by replication progress
#'
#' @param df data.frame with a `rep_percent` column.
#' @param breaks Four increasing boundaries in the unit of `rep_percent`,
#'   e.g. `c(0, 33, 66, 100)` or `c(0, 0.33, 0.66, 1)`. Intervals are closed;
#'   a value on a boundary goes to the earlier class.
#' @return `df` with a factor column `rep_cluster` (levels `early`, `middle`,
#'   `late`; `NA` outside the breaks).
#' @examples
#' add_rep_cluster(data.frame(rep_percent = c(5, 33, 50, 99, 120)),
#'                 breaks = c(0, 33, 66, 100))
#' @export
add_rep_cluster <- function(df, breaks) {
  if (!is.numeric(breaks) || length(breaks) != 4L) {
    stop("'breaks' must be four numeric boundaries.")
  }
  .check_columns(df, "rep_percent")
  x <- df$rep_percent
  inside <- function(lo, hi) !is.na(x) & x >= lo & x <= hi
  cluster <- rep(NA_character_, length(x))
  cluster[inside(breaks[3], breaks[4])] <- "late"
  cluster[inside(breaks[2], breaks[3])] <- "middle"
  cluster[inside(breaks[1], breaks[2])] <- "early"
  df$rep_cluster <- factor(cluster, levels = c("early", "middle", "late"))
  df
}

#' Combine the replication-tempo fits of two states
#'
#' Stacks the per-region results of [fit_rep_tempo()] for two states and
#' derives `twidth = t75 - t25`.
#'
#' @param ctr_df,ama_df Fit tables of the reference and the treatment state,
#'   with columns `region_id`, `t25`, `t50`, `t75`.
#' @inheritParams pivot_state_wide
#' @return A data.frame with `region_id`, `t50`, `twidth` and the factor
#'   `state`; rows with `NA` are dropped.
#' @export
combine_rep_tempo_feature <- function(ctr_df, ama_df, states = c("ctr", "ama")) {
  .check_states(states)
  one <- function(df, state) {
    .check_columns(df, c("region_id", "t25", "t50", "t75"), "fit table")
    data.frame(region_id = df$region_id, t50 = df$t50,
               twidth = df$t75 - df$t25, state = state,
               stringsAsFactors = FALSE)
  }
  out <- rbind(one(ctr_df, states[1]), one(ama_df, states[2]))
  out$state <- factor(out$state, levels = states)
  .drop_na_rows(out)
}

#' Convert a table with a `region_id` column to GRanges
#'
#' @param df data.frame with a `region_id` column formatted `chr:start-end`.
#' @return A `GRanges` with all columns of `df` (including `region_id`) as
#'   metadata.
#' @examples
#' region_id_to_granges(data.frame(region_id = c("chr1:1-100", "chr2:50-80"),
#'                                 score = 1:2))
#' @export
region_id_to_granges <- function(df) {
  .check_columns(df, "region_id")
  df <- as.data.frame(df)
  id <- as.character(df$region_id)
  range <- sub("^[^:]*:", "", id)
  df$seqnames <- sub(":.*", "", id)
  df$start <- as.numeric(sub("-.*", "", range))
  df$end <- as.numeric(sub("^[^-]*-", "", range))
  GenomicRanges::makeGRangesFromDataFrame(
    df,
    keep.extra.columns = TRUE,
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end"
  )
}

#' Keep autosomal ranges only
#'
#' @param gr A `GRanges`.
#' @param chromosomes Chromosomes to keep; the mouse autosomes by default.
#'   The seqlevels of the result are exactly these.
#' @return The filtered `GRanges`.
#' @examples
#' gr <- GenomicRanges::GRanges(c("chr1", "chrX", "chr19"),
#'                              IRanges::IRanges(1, 10))
#' keep_autosomes(gr)
#' @export
keep_autosomes <- function(gr, chromosomes = paste0("chr", 1:19)) {
  gr <- gr[as.character(GenomicRanges::seqnames(gr)) %in% chromosomes]
  GenomeInfoDb::seqlevels(gr, pruning.mode = "coarse") <- chromosomes
  gr
}

#' Annotate regions with a segmented signal track
#'
#' Averages a signal track over genome tiles, segments the tile means with
#' `DNAcopy` (circular binary segmentation after outlier smoothing) and
#' writes each region's segment mean and a two-class label (segment mean
#' `> 0` or not) back to `df`. Typical use: A/B compartment or early/late
#' labels from an eigenvector or log-ratio track.
#'
#' `DNAcopy::segment()` uses random permutations; set a seed beforehand for
#' reproducible segments.
#'
#' @param track_gr `GRanges` with a numeric `score` column and seqlengths.
#' @param df data.frame with a `region_id` column (`chr:start-end`).
#' @param value_col,label_col Names of the columns to write.
#' @param labels Two labels: for segment mean `> 0` and otherwise.
#' @param align_tiles `TRUE`: tiles without signal get `NA`; `FALSE`: every
#'   tile must have signal.
#' @param tilewidth Tile width in bp.
#' @param verbose Let `DNAcopy::segment()` print its progress.
#'
#' @return `df` with `value_col` and `label_col` added (`NA` for regions
#'   overlapping no segment; when a region overlaps several segments the
#'   last one wins).
#' @export
annotate_track_segments <- function(track_gr, df, value_col, label_col, labels,
                                    align_tiles, tilewidth = 100000,
                                    verbose = FALSE) {
  rlang::check_installed("DNAcopy", reason = "to segment the signal track.")
  .check_columns(df, "region_id")
  if (length(labels) != 2L) stop("'labels' must hold two labels.")
  sq <- GenomeInfoDb::seqinfo(track_gr)
  tiles <- GenomicRanges::tileGenome(sq, tilewidth = tilewidth,
                                     cut.last.tile.in.chrom = TRUE)

  hits <- GenomicRanges::findOverlaps(track_gr, tiles)
  scores <- tapply(
    track_gr$score[S4Vectors::queryHits(hits)],
    S4Vectors::subjectHits(hits),
    function(x) mean(x, na.rm = TRUE)
  )
  if (align_tiles) {
    tile_score <- rep(NA_real_, length(tiles))
    tile_score[as.integer(names(scores))] <- as.numeric(scores)
  } else {
    if (length(scores) != length(tiles)) {
      stop("Some tiles have no signal; use align_tiles = TRUE.")
    }
    tile_score <- as.numeric(scores)
  }

  cna_object <- DNAcopy::CNA(
    genomdat = tile_score,
    chrom = as.factor(GenomicRanges::seqnames(tiles)),
    maploc = BiocGenerics::start(tiles),
    data.type = "logratio"
  )
  smoothed_cna <- DNAcopy::smooth.CNA(cna_object)
  segmentation_result <- DNAcopy::segment(smoothed_cna, verbose = as.integer(isTRUE(verbose)))
  segments <- segmentation_result$output

  id <- as.character(df$region_id)
  fit_gr <- GenomicRanges::GRanges(
    seqnames = sub(":.*", "", id),
    ranges = IRanges::IRanges(
      start = as.integer(sub("-.*", "", sub(".*:", "", id))),
      end = as.integer(sub(".*-", "", id))
    )
  )
  segment_gr <- GenomicRanges::GRanges(
    seqnames = segments$chrom,
    ranges = IRanges::IRanges(segments$loc.start, segments$loc.end)
  )
  hits <- suppressWarnings(GenomicRanges::findOverlaps(fit_gr, segment_gr))
  query <- S4Vectors::queryHits(hits)
  subject <- S4Vectors::subjectHits(hits)

  df[[value_col]] <- NA_real_
  df[[label_col]] <- NA_character_
  df[[value_col]][query] <- segments$seg.mean[subject]
  df[[label_col]][query] <- ifelse(segments$seg.mean > 0, labels[1], labels[2])[subject]
  df
}

#' Call regions that differ between two states
#'
#' Regions measured in both states are classified by the log2 fold change
#' `M` (treatment / reference): `<treatment>_up` when `M > thres`,
#' `<treatment>_down` when `M < -thres`, otherwise `no_sig`. Regions with a
#' value in only one state are appended as `<treatment>_up` (treatment only)
#' or `<treatment>_down` (reference only).
#'
#' @inheritParams pivot_state_wide
#' @param thres Threshold on the absolute log2 fold change.
#' @return A data.frame with `id_col` and the factor `iz_group` (levels
#'   `<treatment>_up`, `<treatment>_down`, `no_sig`).
#' @examples
#' long <- data.frame(id = rep(c("a", "b", "c"), 2),
#'                    state = rep(c("ctr", "ama"), each = 3),
#'                    usage = c(1, 2, NA, 4, 2, 5))
#' call_diff_groups(long, "id", "usage", thres = 1)
#' @export
call_diff_groups <- function(df, id_col, feature, thres, states = c("ctr", "ama")) {
  .check_states(states)
  up <- paste0(states[2], "_up")
  down <- paste0(states[2], "_down")
  group_levels <- c(up, down, "no_sig")

  both <- add_ma_columns(
    pivot_state_wide(df, id_col, feature, na = "drop", states = states),
    "log2", states = states
  )
  group <- rep("no_sig", nrow(both))
  group[!is.na(both$M) & both$M < -thres] <- down
  group[!is.na(both$M) & both$M > thres] <- up
  diff_1 <- both[, id_col, drop = FALSE]
  diff_1$iz_group <- factor(group, levels = group_levels)

  either <- pivot_state_wide(df, id_col, feature, na = "keep", states = states)
  .check_columns(either, states, "the wide table")
  ref_na <- is.na(either[[states[1]]])
  trt_na <- is.na(either[[states[2]]])
  group <- rep(NA_character_, nrow(either))
  group[!ref_na & trt_na] <- down
  group[ref_na & !trt_na] <- up
  diff_2 <- either[!is.na(group), id_col, drop = FALSE]
  diff_2$iz_group <- factor(group[!is.na(group)], levels = group_levels)

  out <- rbind(diff_1, diff_2)
  rownames(out) <- NULL
  out
}

#' Per-state values of regions with their difference group
#'
#' @inheritParams pivot_state_wide
#' @param value_col Name of the value column to spread by state.
#' @param diff_groups Output of [call_diff_groups()] (unique `id_col`).
#' @return The wide table of `value_col` (regions measured in both states)
#'   with the columns of `diff_groups` joined; rows with `NA` are dropped.
#' @export
summarise_by_diff_group <- function(df, id_col, value_col, diff_groups,
                                    states = c("ctr", "ama")) {
  if (length(id_col) != 1L) stop("'id_col' must be one column name here.")
  wide <- pivot_state_wide(df, id_col, value_col, na = "drop", states = states)
  .drop_na_rows(.left_join_by_id(wide, as.data.frame(diff_groups), by = id_col))
}

#' Crossing point of the velocity curves of two states
#'
#' Interpolates both curves on a common grid of 1000 points over the shared
#' `rep_percent` range and returns the first sign change of their
#' difference.
#'
#' @param data data.frame with columns `state`, `rep_percent` and
#'   `length_velocity` (bp), e.g. the prediction table of
#'   [fit_track_length_models()].
#' @param state1,state2 The two states to compare.
#' @return A one-row data.frame with `x` (`rep_percent` at the crossing) and
#'   `y` (velocity there, in kb), or `NULL` when the curves do not cross.
#' @examples
#' d <- data.frame(state = rep(c("ctr", "ama"), each = 11),
#'                 rep_percent = rep(seq(0, 100, 10), 2),
#'                 length_velocity = c(seq(1e4, 3e4, length.out = 11),
#'                                     seq(2.5e4, 1.5e4, length.out = 11)))
#' find_intersection(d, "ctr", "ama")
#' @export
find_intersection <- function(data, state1, state2) {
  .check_columns(data, c("state", "rep_percent", "length_velocity"), "data")
  one <- function(s) {
    d <- data[!is.na(data$state) & data$state == s, , drop = FALSE]
    d[order(d$rep_percent), , drop = FALSE]
  }
  df1 <- one(state1)
  df2 <- one(state2)

  common_x <- seq(
    max(min(df1$rep_percent), min(df2$rep_percent)),
    min(max(df1$rep_percent), max(df2$rep_percent)),
    length.out = 1000
  )

  y1 <- stats::approx(df1$rep_percent, df1$length_velocity / 1000, xout = common_x)$y
  y2 <- stats::approx(df2$rep_percent, df2$length_velocity / 1000, xout = common_x)$y

  diff <- y1 - y2
  sign_changes <- which(diff[-1] * diff[-length(diff)] < 0)
  if (!length(sign_changes)) return(NULL)

  idx <- sign_changes[1]
  x_intersect <- common_x[idx] + (common_x[idx + 1] - common_x[idx]) *
    abs(diff[idx]) / (abs(diff[idx]) + abs(diff[idx + 1]))
  y_intersect <- stats::approx(common_x, y1, xout = x_intersect)$y
  data.frame(x = x_intersect, y = y_intersect)
}

#' Format base-pair numbers as readable labels
#'
#' @param x Numeric vector of lengths in bp.
#' @return Character vector such as `"20 kb"`, `"1.5 Mb"`; `NA` stays `NA`.
#' @examples
#' format_bp(c(500, 20000, 1.5e6, 2.7e9, NA))
#' @export
format_bp <- function(x) {
  vapply(x, function(val) {
    if (is.na(val)) return(NA_character_)
    if (val >= 1e9) {
      paste0(round(val / 1e9, 2), " Gb")
    } else if (val >= 1e6) {
      paste0(round(val / 1e6, 2), " Mb")
    } else if (val >= 1e3) {
      paste0(round(val / 1e3, 0), " kb")
    } else {
      paste0(val, " bp")
    }
  }, character(1), USE.NAMES = FALSE)
}
