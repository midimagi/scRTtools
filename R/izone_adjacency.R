# Co-firing of adjacent initiation zones.
# make_pair_init() and summarise_pair_cofiring_bins() come from notebook
# 4.3.1.analyze_adjacent_izones; make_adjacent_pairs(), build_firing_matrix(),
# assign_rep_stage(), summarise_cell_adjacency(), summarise_pair_adjacency()
# and summarise_adjacency_by_stage() from
# 4.3.1.analyze_adjacent_izones_features_per_cell (inline code). Plots are not
# included.

#' Pairs of IZs separated by `n_adj` IZs
#'
#' Order IZs by chromosome (seqlevel order), start and end, and pair each IZ
#' with the IZ `n_adj + 1` positions to its right on the same chromosome
#' (`n_adj = 0`: adjacent IZs).
#'
#' @param izones `GRanges` of IZs with an `iz_id` column.
#' @param n_adj Number of IZs between the two members of a pair.
#'
#' @return A `data.frame` with `pair_id` (factor, `"left_right"`), `chr`,
#'   `left_iz_id` and `right_iz_id` (character).
#'
#' @export
make_pair_init <- function(izones, n_adj = 0L) {
  if (!inherits(izones, "GRanges") || !"iz_id" %in% names(S4Vectors::mcols(izones))) {
    stop("'izones' must be a GRanges with an 'iz_id' column.")
  }
  n_adj <- as.integer(n_adj)
  izones <- izones[order(as.integer(GenomeInfoDb::seqnames(izones)),
                         BiocGenerics::start(izones), BiocGenerics::end(izones))]
  iz <- data.frame(
    iz_id = as.character(izones$iz_id),
    chr = as.character(GenomeInfoDb::seqnames(izones)),
    start = BiocGenerics::start(izones),
    end = BiocGenerics::end(izones),
    center = floor((BiocGenerics::start(izones) + BiocGenerics::end(izones)) / 2),
    stringsAsFactors = FALSE
  )
  N <- nrow(iz)
  step <- n_adj + 1L
  left <- if (N > step) seq_len(N - step) else integer()
  right <- left + step
  # keep pairs with both ends on the same chromosome
  if (length(left)) {
    keep <- iz$chr[left] == iz$chr[right]
    left <- left[keep]
    right <- right[keep]
  }
  if (length(left) == 0) {
    pairs <- data.frame(
      pair_id = character(), chr = character(), left_iz_id = character(),
      right_iz_id = character(), stringsAsFactors = FALSE
    )
    pairs$pair_id <- factor(pairs$pair_id, levels = character())
    return(pairs)
  }
  pairs <- data.frame(
    pair_id = paste(iz$iz_id[left], iz$iz_id[right], sep = "_"),
    chr = iz$chr[left],
    left_iz_id = iz$iz_id[left],
    right_iz_id = iz$iz_id[right],
    stringsAsFactors = FALSE
  )
  pairs$pair_id <- factor(pairs$pair_id, levels = pairs$pair_id)
  pairs
}

#' Co-firing of IZ pairs in replication-progress bins
#'
#' An IZ is fired in a cell when one of the cell's tracks is assigned to it
#' (`iz_id`). Cells are grouped into `bin_num` equal bins of `rep_percent`
#' (cells without `rep_percent` form one `NA` bin); per pair and bin the
#' numbers of cells with the left / right IZ fired are summed, and per pair
#' the bins are classified by which members fired at least once.
#'
#' @param pairs Output of [make_pair_init()].
#' @param tracks Named `GRangesList` of tracks with an `iz_id` column (e.g.
#'   `define_initiation_zones()$tracks`).
#' @param meta_data A `data.frame` with columns `cell_id` and `rep_percent`.
#' @param bin_num Number of `rep_percent` bins.
#'
#' @return A list with `pair_bins` (per pair and bin: `bin_mid`, median
#'   `rep_percent`, summed `left_iz_state` / `right_iz_state`) and
#'   `pair_fired_info` (per pair: numbers of bins with neither / only left /
#'   only right / both fired, `n_all` and the corresponding rates).
#'
#' @export
summarise_pair_cofiring_bins <- function(pairs, tracks, meta_data, bin_num = 100L) {
  need <- c("pair_id", "left_iz_id", "right_iz_id")
  if (!is.data.frame(pairs) || !all(need %in% names(pairs))) {
    stop("'pairs' must be the output of make_pair_init().")
  }
  if (!inherits(tracks, "GRangesList") || is.null(names(tracks))) {
    stop("'tracks' must be a named GRangesList.")
  }
  if (!is.data.frame(meta_data) || !all(c("cell_id", "rep_percent") %in% names(meta_data))) {
    stop("'meta_data' must contain columns 'cell_id' and 'rep_percent'.")
  }
  cells <- names(tracks)
  n_pairs <- nrow(pairs)

  # per-cell fired state of the left / right members (pairs x cells)
  left_state <- matrix(0, n_pairs, length(cells))
  right_state <- matrix(0, n_pairs, length(cells))
  for (k in seq_along(cells)) {
    putative_iz <- stats::na.omit(sort(S4Vectors::mcols(tracks[[k]])$iz_id))
    left_state[, k] <- ifelse(pairs$left_iz_id %in% putative_iz, 1, 0)
    right_state[, k] <- ifelse(pairs$right_iz_id %in% putative_iz, 1, 0)
  }

  # rep_percent bins (cells without rep_percent -> NA bin)
  map <- stats::na.omit(meta_data[, c("cell_id", "rep_percent")])
  map <- map[order(map$rep_percent), , drop = FALSE]
  bin_breaks <- seq(0, 1, length.out = bin_num + 1)
  map$rep_bin <- cut(map$rep_percent, breaks = bin_breaks, labels = FALSE,
                     include.lowest = TRUE)
  map$bin_mid <- (map$rep_bin - 0.5) / bin_num
  j <- match(cells, map$cell_id)
  cell_bin <- map$bin_mid[j]
  cell_rep <- map$rep_percent[j]

  bins <- sort(unique(cell_bin), na.last = TRUE)
  cell_groups <- lapply(bins, function(b) if (is.na(b)) which(is.na(cell_bin)) else
    which(!is.na(cell_bin) & cell_bin == b))
  bin_rep <- vapply(cell_groups, function(ii) stats::median(cell_rep[ii]), numeric(1))
  L <- vapply(cell_groups, function(ii) rowSums(left_state[, ii, drop = FALSE]),
              numeric(n_pairs))
  R <- vapply(cell_groups, function(ii) rowSums(right_state[, ii, drop = FALSE]),
              numeric(n_pairs))
  L <- matrix(L, nrow = n_pairs)
  R <- matrix(R, nrow = n_pairs)
  n_bins <- length(bins)

  pair_bins <- data.frame(
    pair_id = rep(pairs$pair_id, each = n_bins),
    bin_mid = rep(bins, times = n_pairs),
    rep_percent = rep(bin_rep, times = n_pairs),
    left_iz_id = rep(pairs$left_iz_id, each = n_bins),
    right_iz_id = rep(pairs$right_iz_id, each = n_bins),
    left_iz_state = as.vector(t(L)),
    right_iz_state = as.vector(t(R)),
    stringsAsFactors = FALSE
  )

  n_both_no_fired <- rowSums(L == 0 & R == 0)
  n_left_fired <- rowSums(L != 0 & R == 0)
  n_right_fired <- rowSums(L == 0 & R != 0)
  n_both_fired <- rowSums(L != 0 & R != 0)
  n_all <- n_both_no_fired + n_left_fired + n_right_fired + n_both_fired
  pair_fired_info <- data.frame(
    pair_id = pairs$pair_id,
    n_both_no_fired = as.integer(n_both_no_fired),
    n_left_fired = as.integer(n_left_fired),
    n_right_fired = as.integer(n_right_fired),
    n_both_fired = as.integer(n_both_fired)
  )
  pair_fired_info$n_all <- as.integer(n_all)
  pair_fired_info$both_no_fired_rate <- pair_fired_info$n_both_no_fired / pair_fired_info$n_all
  pair_fired_info$left_fired_rate <- pair_fired_info$n_left_fired / pair_fired_info$n_all
  pair_fired_info$right_fired_rate <- pair_fired_info$n_right_fired / pair_fired_info$n_all
  pair_fired_info$both_fired_rate <- pair_fired_info$n_both_fired / pair_fired_info$n_all
  pair_fired_info <- pair_fired_info[order(pair_fired_info$pair_id), , drop = FALSE]
  rownames(pair_fired_info) <- NULL

  list(pair_bins = pair_bins, pair_fired_info = pair_fired_info)
}

#' Adjacent IZ pairs with firing-range features
#'
#' Order IZs by chromosome (seqlevel order), start and end, attach the firing
#' range of each IZ (matched by `iz_id`) and build all adjacent pairs on the
#' same chromosome with distances and firing-range overlaps.
#'
#' A firing range is valid when both bounds are finite, `earliest >= 0`,
#' `latest <= 1` and `earliest <= latest`; `firing_range_overlap` is computed
#' for pairs with two valid ranges.
#'
#' @param izones `GRanges` of IZs with an `iz_id` column.
#' @param fire_ranges Output of [calculate_firing_range()] (`iz_id`,
#'   `earliest_rep_percent`, `latest_unrep_percent`, optionally
#'   `corroborated`), with the same coordinates for matching ids.
#'
#' @return A list with `izones` (ordered), `features` (one row per IZ) and
#'   `pairs` (one row per adjacent pair: indices into `features`, ids,
#'   coordinates, `center_distance_bp`, `gap_bp`, `overlap_bp`, left/right
#'   firing-range columns, `firing_range_overlap`,
#'   `earliest_rep_difference`).
#'
#' @export
make_adjacent_pairs <- function(izones, fire_ranges) {
  stopifnot(inherits(izones, "GRanges"), inherits(fire_ranges, "GRanges"))
  mc <- function(x) S4Vectors::mcols(x)
  required <- c("iz_id", "earliest_rep_percent", "latest_unrep_percent")
  if (!"iz_id" %in% names(mc(izones)) ||
      !all(required %in% names(mc(fire_ranges)))) stop("Missing IZ metadata.")
  ids <- as.character(izones$iz_id)
  fr_ids <- as.character(fire_ranges$iz_id)
  if (!length(ids) || anyNA(ids) || any(ids == "") || anyDuplicated(ids) ||
      anyNA(fr_ids) || any(fr_ids == "") || anyDuplicated(fr_ids)) {
    stop("IZ IDs must be nonmissing and unique.")
  }
  sn <- function(x) GenomeInfoDb::seqnames(x)
  st <- function(x) BiocGenerics::start(x)
  en <- function(x) BiocGenerics::end(x)
  izones <- izones[order(as.integer(sn(izones)), st(izones), en(izones))]
  iz <- data.frame(iz_id = as.character(izones$iz_id), chr = as.character(sn(izones)),
                   start = st(izones), end = en(izones),
                   center = floor((st(izones) + en(izones)) / 2))
  idx <- match(iz$iz_id, fr_ids)
  found <- !is.na(idx)
  if (any(iz$chr[found] != as.character(sn(fire_ranges))[idx[found]] |
          iz$start[found] != st(fire_ranges)[idx[found]] |
          iz$end[found] != en(fire_ranges)[idx[found]])) {
    stop("IZ coordinates disagree between input files for matching iz_id.")
  }
  iz$earliest_rep_percent <- fire_ranges$earliest_rep_percent[idx]
  iz$latest_unrep_percent <- fire_ranges$latest_unrep_percent[idx]
  iz$valid_firing_range <- is.finite(iz$earliest_rep_percent) &
    is.finite(iz$latest_unrep_percent) & iz$earliest_rep_percent >= 0 &
    iz$latest_unrep_percent <= 1 & iz$earliest_rep_percent <= iz$latest_unrep_percent
  if ("corroborated" %in% names(mc(fire_ranges))) iz$corroborated <- fire_ranges$corroborated[idx]
  left <- if (nrow(iz) > 1L) seq_len(nrow(iz) - 1L) else integer()
  left <- left[iz$chr[left] == iz$chr[left + 1L]]
  if (!length(left)) stop("No within-chromosome adjacent IZ pairs.")
  right <- left + 1L
  pairs <- data.frame(
    pair_id = paste(iz$iz_id[left], iz$iz_id[right], sep = "_"),
    chr = iz$chr[left], left_index = left, right_index = right,
    left_iz_id = iz$iz_id[left], right_iz_id = iz$iz_id[right],
    left_start = iz$start[left], left_end = iz$end[left],
    right_start = iz$start[right], right_end = iz$end[right],
    center_distance_bp = iz$center[right] - iz$center[left],
    gap_bp = pmax(0, iz$start[right] - iz$end[left] - 1),
    overlap_bp = pmax(0, pmin(iz$end[left], iz$end[right]) - iz$start[right] + 1)
  )
  for (nm in intersect(c("earliest_rep_percent", "latest_unrep_percent",
                         "valid_firing_range", "corroborated"), names(iz))) {
    pairs[[paste0("left_", nm)]] <- iz[[nm]][left]
    pairs[[paste0("right_", nm)]] <- iz[[nm]][right]
  }
  valid <- iz$valid_firing_range[left] & iz$valid_firing_range[right]
  pairs$firing_range_overlap <- NA_real_
  pairs$firing_range_overlap[valid] <- pmax(
    0,
    pmin(iz$latest_unrep_percent[left[valid]], iz$latest_unrep_percent[right[valid]]) -
      pmax(iz$earliest_rep_percent[left[valid]], iz$earliest_rep_percent[right[valid]])
  )
  pairs$earliest_rep_difference <- abs(iz$earliest_rep_percent[left] -
                                         iz$earliest_rep_percent[right])
  list(izones = izones, features = iz, pairs = pairs)
}

#' Cell x IZ firing matrix from track centres
#'
#' An IZ is fired in a cell when the centre of one of the cell's
#' IZ-assigned tracks (optionally only tracks shorter than
#' `max_track_width_bp`) falls inside it.
#'
#' @param izones `GRanges` of IZs with an `iz_id` column (e.g.
#'   `make_adjacent_pairs()$izones`).
#' @param tracks Named `GRangesList` of tracks with columns `cell_id`
#'   (matching the names), `rep_percent` and `iz_id`.
#' @param meta_data A `data.frame` with columns `cell_id` (unique) and
#'   `rep_percent`; missing values are taken from the tracks.
#' @param max_track_width_bp `NULL` or the exclusive maximum track width.
#'
#' @return A list with `fired` (sparse 0/1 `dgCMatrix`, cells x IZs) and
#'   `cell_info` (`cell_id`, `rep_percent`, `n_tracks`, `n_retained_tracks`,
#'   `n_fired_izones`, `n_tracks_hitting_iz`).
#'
#' @export
build_firing_matrix <- function(izones, tracks, meta_data, max_track_width_bp = NULL) {
  if (!inherits(tracks, "GRangesList") || !length(tracks)) stop("No track cells.")
  cell_ids <- names(tracks)
  if (is.null(cell_ids) || anyNA(cell_ids) || any(cell_ids == "") || anyDuplicated(cell_ids)) {
    stop("names(tracks) must contain unique, nonmissing cell IDs.")
  }
  if (!all(c("cell_id", "rep_percent") %in% names(meta_data)) ||
      anyNA(meta_data$cell_id) || anyDuplicated(meta_data$cell_id) ||
      !is.numeric(meta_data$rep_percent)) stop("Invalid cell metadata.")
  if (!is.null(max_track_width_bp) &&
      (length(max_track_width_bp) != 1L || !is.finite(max_track_width_bp) ||
       max_track_width_bp <= 0)) stop("Invalid max_track_width_bp.")
  all_tracks <- BiocGenerics::unlist(tracks, use.names = FALSE)
  md <- S4Vectors::mcols(all_tracks)
  if (!all(c("cell_id", "rep_percent", "iz_id") %in% names(md)) ||
      !is.numeric(md$rep_percent)) stop("Missing or invalid track metadata.")
  owner <- rep(seq_along(tracks), S4Vectors::elementNROWS(tracks))
  if (anyNA(md$cell_id) || any(as.character(md$cell_id) != cell_ids[owner])) {
    stop("Track cell_id does not agree with names(tracks).")
  }
  track_rep <- unique(data.frame(cell_id = as.character(md$cell_id), rep_percent = md$rep_percent))
  track_rep <- track_rep[!is.na(track_rep$rep_percent), , drop = FALSE]
  if (anyDuplicated(track_rep$cell_id)) stop("Inconsistent rep_percent within a cell.")
  cell_info <- data.frame(cell_id = cell_ids,
                          rep_percent = meta_data$rep_percent[match(cell_ids, meta_data$cell_id)])
  from_tracks <- track_rep$rep_percent[match(cell_ids, track_rep$cell_id)]
  both <- !is.na(cell_info$rep_percent) & !is.na(from_tracks)
  if (any(abs(cell_info$rep_percent[both] - from_tracks[both]) > 1e-8)) {
    stop("rep_percent differs between metadata and tracks.")
  }
  missing <- is.na(cell_info$rep_percent)
  cell_info$rep_percent[missing] <- from_tracks[missing]
  rp <- cell_info$rep_percent
  if (any(!is.na(rp) & (!is.finite(rp) | rp < 0 | rp > 1))) stop("rep_percent must be in [0, 1].")
  cell_info$n_tracks <- S4Vectors::elementNROWS(tracks)
  keep <- !is.na(md$iz_id)
  if (any(!as.character(md$iz_id[keep]) %in% as.character(izones$iz_id))) {
    stop("Tracks reference IZ IDs absent from gr_izones.")
  }
  if (!is.null(max_track_width_bp)) keep <- keep & BiocGenerics::width(all_tracks) < max_track_width_bp
  cell_info$n_retained_tracks <- tabulate(owner[keep], nbins = length(cell_ids))
  retained <- all_tracks[keep]
  retained_owner <- owner[keep]
  # Reuse the GRanges structure so zero retained tracks remain a valid empty object.
  centers <- retained
  midpoint <- floor((BiocGenerics::start(retained) + BiocGenerics::end(retained)) / 2)
  IRanges::ranges(centers) <- IRanges::IRanges(start = midpoint, end = midpoint)
  hits <- GenomicRanges::findOverlaps(centers, izones, ignore.strand = TRUE)
  # Multiple tracks hitting the same cell/IZ are collapsed to one event.
  fired <- Matrix::sparseMatrix(
    i = retained_owner[S4Vectors::queryHits(hits)], j = S4Vectors::subjectHits(hits),
    x = rep(1, length(hits)), dims = c(length(cell_ids), length(izones)),
    dimnames = list(cell_ids, as.character(izones$iz_id))
  )
  fired <- (fired > 0) * 1
  cell_info$n_fired_izones <- as.numeric(Matrix::rowSums(fired))
  cell_info$n_tracks_hitting_iz <- tabulate(retained_owner[unique(S4Vectors::queryHits(hits))],
                                            nbins = length(cell_ids))
  list(fired = fired, cell_info = cell_info)
}

#' Replication stage labels from replication percentage
#'
#' `cut(rep_percent, breaks, labels, right = FALSE, include.lowest = TRUE)`,
#' with `rep_percent == 1` in the last stage and `NA` as `"unknown"`.
#'
#' @param rep_percent Numeric vector of replicated fractions.
#' @param breaks,labels Stage breaks and labels (e.g. quartiles
#'   `Q1`-`Q4`).
#'
#' @return A factor with levels `c(labels, "unknown")`.
#'
#' @export
assign_rep_stage <- function(rep_percent,
                             breaks = c(0, 0.25, 0.5, 0.75, 1),
                             labels = c("Q1", "Q2", "Q3", "Q4")) {
  stage <- as.character(cut(rep_percent, breaks = breaks, labels = labels,
                            right = FALSE, include.lowest = TRUE))
  stage[!is.na(rep_percent) & rep_percent == 1] <- labels[length(labels)]
  stage[is.na(stage)] <- "unknown"
  factor(stage, levels = c(labels, "unknown"))
}

#' Ratio with non-finite results set to NA
#' @keywords internal
#' @noRd
.safe_ratio <- function(numerator, denominator) {
  out <- numerator / denominator
  out[!is.finite(out)] <- NA_real_
  out
}

#' Co-firing statistics from 2 x 2 counts
#' @keywords internal
#' @noRd
.binary_pair_stats <- function(n11, n10, n01, n00) {
  n <- n11 + n10 + n01 + n00
  p_left <- .safe_ratio(n11 + n10, n)
  p_right <- .safe_ratio(n11 + n01, n)
  p11 <- .safe_ratio(n11, n)
  expected <- p_left * p_right
  data.frame(
    n_observations = n, n11 = n11, n10 = n10, n01 = n01, n00 = n00,
    p_left = p_left, p_right = p_right, p11 = p11, p_expected = expected,
    p_excess = p11 - expected, cofire_lift = .safe_ratio(p11, expected),
    phi = .safe_ratio(p11 - expected,
                      sqrt(p_left * (1 - p_left) * p_right * (1 - p_right))),
    p_right_given_left = .safe_ratio(n11, n11 + n10),
    p_left_given_right = .safe_ratio(n11, n11 + n01),
    p_right_given_left_unfired = .safe_ratio(n01, n01 + n00)
  )
}

#' Per-cell co-firing of adjacent IZs
#'
#' For each cell, cross-tabulate the fired state of the left and right IZ of
#' every adjacent pair and derive co-firing statistics (`p11`, expected
#' `p_left * p_right`, excess, lift, phi, conditional probabilities). The
#' chromosome-shuffle expectation is the expected number of fired adjacent
#' pairs if the cell's fired IZs were placed at random within each
#' chromosome.
#'
#' @param fired Output of [build_firing_matrix()]`$fired`.
#' @param pairs,iz_features `make_adjacent_pairs()$pairs` and `$features`.
#'
#' @return A `data.frame`, one row per cell (rows of `fired`), with the 2 x 2
#'   counts, probabilities and `p_expected_chr_shuffle`,
#'   `cofire_lift_chr_shuffle`, `p_excess_chr_shuffle`.
#'
#' @export
summarise_cell_adjacency <- function(fired, pairs, iz_features) {
  left <- fired[, pairs$left_index, drop = FALSE]
  right <- fired[, pairs$right_index, drop = FALSE]
  n11 <- as.numeric(Matrix::rowSums(left * right))
  n10 <- as.numeric(Matrix::rowSums(left)) - n11
  n01 <- as.numeric(Matrix::rowSums(right)) - n11
  stats <- .binary_pair_stats(n11, n10, n01, nrow(pairs) - n11 - n10 - n01)
  expected_n11 <- numeric(nrow(fired))
  for (idx in split(seq_len(ncol(fired)), iz_features$chr)) {
    n <- length(idx)
    if (n < 2L) next
    k <- as.numeric(Matrix::rowSums(fired[, idx, drop = FALSE]))
    expected_n11 <- expected_n11 + (n - 1) * k * (k - 1) / (n * (n - 1))
  }
  stats$p_expected_chr_shuffle <- expected_n11 / nrow(pairs)
  stats$cofire_lift_chr_shuffle <- .safe_ratio(n11, expected_n11)
  stats$p_excess_chr_shuffle <- stats$p11 - stats$p_expected_chr_shuffle
  stats
}

#' Per-pair co-firing of adjacent IZs across selected cells
#'
#' @param fired Output of [build_firing_matrix()]`$fired`.
#' @param pairs `make_adjacent_pairs()$pairs`.
#' @param selected_cells Row indices (or names) of `fired` to use.
#' @param label Value of the `rep_stage` column.
#'
#' @return `pairs` with `rep_stage` and the co-firing statistics of each pair
#'   over the selected cells.
#'
#' @export
summarise_pair_adjacency <- function(fired, pairs, selected_cells, label) {
  left <- fired[selected_cells, pairs$left_index, drop = FALSE]
  right <- fired[selected_cells, pairs$right_index, drop = FALSE]
  n11 <- as.numeric(Matrix::colSums(left * right))
  n10 <- as.numeric(Matrix::colSums(left)) - n11
  n01 <- as.numeric(Matrix::colSums(right)) - n11
  stats <- .binary_pair_stats(n11, n10, n01, nrow(left) - n11 - n10 - n01)
  cbind(pairs, rep_stage = label, stats)
}

#' Adjacent-IZ co-firing by replication stage
#'
#' Run [summarise_pair_adjacency()] for all cells and for each stage, and
#' summarise cells and pairs per stage (means of finite values; median pair
#' phi).
#'
#' @param fired Output of [build_firing_matrix()]`$fired`.
#' @param pairs `make_adjacent_pairs()$pairs`.
#' @param cell_stats Per-cell table with `n_fired_izones`, `p11`, `phi`,
#'   `p_expected_chr_shuffle` (e.g. `cbind(cell_info,
#'   summarise_cell_adjacency(...))`), rows matching `fired`.
#' @param rep_stage Stage of each cell (e.g. [assign_rep_stage()]).
#' @param stages Stages to summarise besides `"all"`.
#'
#' @return A list with `pair_stats` (stacked per-pair statistics, `rep_stage`
#'   = `"all"` or a stage) and `stage_summary` (one row per group).
#'
#' @export
summarise_adjacency_by_stage <- function(fired, pairs, cell_stats, rep_stage,
                                         stages = c("Q1", "Q2", "Q3", "Q4")) {
  cell_groups <- c(list(all = seq_len(nrow(fired))),
                   stats::setNames(lapply(stages, function(x) which(rep_stage == x)), stages))
  pair_stats <- do.call(rbind, lapply(names(cell_groups), function(label) {
    summarise_pair_adjacency(fired, pairs, cell_groups[[label]], label)
  }))
  rownames(pair_stats) <- NULL
  finite_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
  finite_median <- function(x) if (any(is.finite(x))) stats::median(x[is.finite(x)]) else NA_real_
  stage_summary <- do.call(rbind, lapply(names(cell_groups), function(label) {
    cells <- cell_stats[cell_groups[[label]], , drop = FALSE]
    ps <- pair_stats[pair_stats$rep_stage == label, , drop = FALSE]
    data.frame(
      rep_stage = label, n_cells = nrow(cells), n_pairs = nrow(ps),
      mean_fired_izones = finite_mean(cells$n_fired_izones),
      mean_cell_p11 = finite_mean(cells$p11), mean_cell_phi = finite_mean(cells$phi),
      mean_cell_p_expected_chr_shuffle = finite_mean(cells$p_expected_chr_shuffle),
      median_pair_phi = finite_median(ps$phi),
      n_pairs_with_defined_phi = sum(is.finite(ps$phi))
    )
  }))
  list(pair_stats = pair_stats, stage_summary = stage_summary)
}
