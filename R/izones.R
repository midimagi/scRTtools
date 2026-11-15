# Initiation zones (IZs): definition, firing range / order and IZ velocity.
# define_initiation_zones() is migrated from
# rcode/processing/define_initiation_zones_fixed.R (the version used by all
# notebooks); calculate_firing_range(), calculate_firing_order(),
# predict_firing_order() and fit_all_izone_velocities() from rcode/analysis/
# (plot_firing_order_frequency() and write_izone_velocity_results() are not
# included).

#' Define single-cell replication initiation zones
#'
#' Cluster single-cell replication tracks into initiation zones (IZs) seeded by
#' replication-timing peaks (e.g. [find_rt_peaks()] on
#' [count_track_centers()] output). Per chromosome, tracks are processed from
#' short to long. Tracks flagged as bounded by a gap (`flag != 0`) or longer
#' than `max_track_length` are skipped. Each track, widened by `2.5 * bin_size`
#' on both sides, either seeds a new putative IZ, joins the single putative IZ
#' it overlaps (by more than `min_overlap`, plus one bin for tracks longer than
#' `short_track_th`; if two are overlapped the one with a compatible
#' track-centre interval is used), or is skipped (more than two overlaps, a
#' track of the same cell already in the IZ, or incompatible centres). Short
#' tracks update the putative IZ interval. Putative IZs supported by fewer
#' than `min_usage` tracks are removed, all tracks are reassigned once more,
#' and IZs with fewer than `min_usage` reassigned tracks are removed.
#'
#' Each final IZ spans the 5%-95% quantiles of its track centres padded by
#' half a bin on each side (the "fixed" version: `start = round(q5 -
#' bin_size / 2)`, `width = q95 - q5 + bin_size`). IZs are numbered by genomic
#' order.
#'
#' Lengths can be given as numbers (bp) or strings such as `"500Kb"`, `"1Mb"`,
#' `"0bp"`.
#'
#' @param tracks Named `GRangesList` of tracks (one element per cell, unique
#'   names) with a `flag` metadata column (e.g. [define_replication_tracks()]).
#' @param rt_peaks `GRanges` of seed peaks.
#' @param short_track_th Tracks up to this length define IZ widths.
#' @param max_track_length Longer tracks are not used.
#' @param min_overlap Minimum overlap between a widened track and an IZ.
#' @param min_usage Minimum number of supporting tracks of an IZ (default 3;
#'   lower it for small datasets).
#' @param bin_size Bin size of the tracks.
#' @param BPPARAM A `BiocParallel::BiocParallelParam`; chromosomes are
#'   processed in parallel.
#' @param verbose Logical; print progress messages. Default `FALSE`.
#'
#' @return A list with `tracks` (input tracks with an `iz_id` metadata column,
#'   `NA` if unassigned) and `izones` (sorted `GRanges` with `iz_id`, `q5`,
#'   `q25`, `q50`, `q75`, `q95` and `n_tracks`).
#'
#' @export
define_initiation_zones <- function(
    tracks,
    rt_peaks,
    short_track_th = "500Kb",
    max_track_length = "1Mb",
    min_overlap = "0bp",
    min_usage = 3L,
    bin_size = "20Kb",
    BPPARAM = BiocParallel::SerialParam(),
    verbose = FALSE
) {

  short_track_th <- .parse_coordinate(short_track_th)
  max_track_length <- .parse_coordinate(max_track_length)
  min_overlap <- .parse_coordinate(min_overlap)
  bin_size <- .parse_coordinate(bin_size)

  short_track_th_bins <- short_track_th / bin_size

  if (!inherits(tracks, "GRangesList")) {
    stop("`tracks` must be a GRangesList or CompressedGRangesList.")
  }

  if (!inherits(rt_peaks, "GRanges")) {
    stop("`rt_peaks` must be a GRanges object.")
  }

  if (!inherits(BPPARAM, "BiocParallelParam")) {
    stop("'BPPARAM' must be a BiocParallelParam object.")
  }

  if (is.null(names(tracks)) || anyNA(names(tracks)) || anyDuplicated(names(tracks))) {
    stop("`names(tracks)` must contain non-missing unique cell names.")
  }
  #min_usage = floor(as.numeric(min_usage) * length(names(tracks)))

  .vmsg(verbose, "Unlisting tracks.")

  gr_all <- BiocGenerics::unlist(tracks, use.names = FALSE)

  if (length(gr_all) == 0L) {
    S4Vectors::mcols(gr_all)$iz_id <- integer(0)
    return(list(
      tracks = IRanges::relist(gr_all, tracks),
      izones = GenomicRanges::GRanges()
    ))
  }

  if (!"flag" %in% colnames(S4Vectors::mcols(gr_all))) {
    stop("`mcols(tracks)` must contain column `flag`.")
  }

  cell_id <- rep(names(tracks), lengths(tracks))
  original_index <- seq_along(gr_all)

  S4Vectors::mcols(gr_all)$.cell_id <- cell_id
  S4Vectors::mcols(gr_all)$.original_index <- original_index
  S4Vectors::mcols(gr_all)$iz_id <- NA_integer_

  chrs <- intersect(GenomeInfoDb::seqlevels(gr_all), GenomeInfoDb::seqlevels(rt_peaks))

  .vmsg(verbose, "Processing ", length(chrs), " chromosome(s) in parallel.")

  process_one_chr <- function(chr) {
    gr_chr <- gr_all[as.character(GenomeInfoDb::seqnames(gr_all)) == chr]
    peaks_chr <- rt_peaks[as.character(GenomeInfoDb::seqnames(rt_peaks)) == chr]

    if (length(gr_chr) == 0L) {
      return(list(
        chr = chr,
        assignment = data.frame(
          original_index = integer(0),
          local_iz_id = integer(0)
        ),
        chr_izs = GenomicRanges::GRanges()
      ))
    }

    track_width <- BiocGenerics::width(gr_chr)
    width_bins <- track_width / bin_size
    center <- floor((BiocGenerics::start(gr_chr) + BiocGenerics::end(gr_chr)) / 2)

    ord <- order(width_bins, center, na.last = TRUE)
    gr_chr <- gr_chr[ord]
    track_width <- track_width[ord]
    width_bins <- width_bins[ord]
    center <- center[ord]

    n_tracks <- length(gr_chr)

    flag <- S4Vectors::mcols(gr_chr)$flag
    is_bounded_by_gap <- !is.na(flag) & flag != 0
    is_too_long <- BiocGenerics::width(gr_chr) > max_track_length

    blur_bp <- 2.5 * bin_size ## set resolution of centers
    blur_start <- BiocGenerics::start(gr_chr) - blur_bp
    blur_end <- BiocGenerics::end(gr_chr) + blur_bp
    center_left <- (blur_start + BiocGenerics::end(gr_chr)) / 2
    center_right <- (BiocGenerics::start(gr_chr) + blur_end) / 2

    if (length(peaks_chr) > 0L) {
      peak_center <- floor(
        (BiocGenerics::start(peaks_chr) + BiocGenerics::end(peaks_chr)) / 2
      )

      putative_IZs <- data.frame(
        start = peak_center - blur_bp,
        end = peak_center + blur_bp,
        peak_id = seq_along(peaks_chr),
        stringsAsFactors = FALSE
      )
    } else {
      putative_IZs <- data.frame(
        start = numeric(0),
        end = numeric(0),
        peak_id = integer(0)
      )
    }

    cluster_id1 <- rep(NA_integer_, n_tracks)
    
  ## 1st loop to define putative IZs
    for (t in seq_len(n_tracks)) {
      # skip track bounded by a gap, or too long
      if (is_bounded_by_gap[t] || is_too_long[t]) {next}

      if (nrow(putative_IZs) == 0L) {
        putative_IZs <- rbind(
          putative_IZs,
          data.frame(
            start = center_left[t],
            end = center_right[t],
            peak_id = NA_integer_
          )
        )

        cluster_id1[t] <- nrow(putative_IZs)
        next
      }

      left1 <- pmax(blur_start[t], putative_IZs$start)
      right1 <- pmin(blur_end[t], putative_IZs$end)

      overlap_threshold <- if (width_bins[t] <= short_track_th_bins) {
        min_overlap
      } else {
        min_overlap + bin_size
      }

      is_overlapping <- right1 - left1 > overlap_threshold
      num_overlaps <- sum(is_overlapping)

      if (num_overlaps > 2L) {next}

      if (num_overlaps == 0L) {
        putative_IZs <- rbind(
          putative_IZs,
          data.frame(
            start = center_left[t],
            end = center_right[t],
            peak_id = NA_integer_
          )
        )

        cluster_id1[t] <- nrow(putative_IZs)
        next
      }

      # ------------------------------------------------------------
      # If current track overlaps exactly 2 IZs, use blur-center
      # compatibility with previous tracks to decide which IZ it belongs to.
      # ------------------------------------------------------------
      if (num_overlaps == 2L) {
        candidate_clusters <- which(is_overlapping)

        center_compatible <- logical(length(candidate_clusters))
        names(center_compatible) <- candidate_clusters

        for (i in seq_along(candidate_clusters)) {
          cl <- candidate_clusters[i]

          is_in_candidate_cluster <- !is.na(cluster_id1) & cluster_id1 == cl
          clustered_idx <- which(is_in_candidate_cluster)

          if (length(clustered_idx) == 0L) {
            center_compatible[i] <- FALSE
            next
          }

          left2 <- pmax(center_left[t], center_left[clustered_idx])
          right2 <- pmin(center_right[t], center_right[clustered_idx])
          center_compatible[i] <- any(right2 - left2 >= 0, na.rm = TRUE)
        }

        if (sum(center_compatible, na.rm = TRUE) != 1L) {next}
        
        overlapped_cluster <- candidate_clusters[which(center_compatible)]

      } else {
        overlapped_cluster <- which(is_overlapping)
      }

      is_in_same_cluster <- !is.na(cluster_id1) & cluster_id1 == overlapped_cluster

      if (sum(is_in_same_cluster, na.rm = TRUE) == 0L) {
        cluster_id1[t] <- overlapped_cluster
        next
      }

      same_cell <- S4Vectors::mcols(gr_chr)$.cell_id[is_in_same_cluster] ==
        S4Vectors::mcols(gr_chr)$.cell_id[t]

      if (any(same_cell, na.rm = TRUE)) {next}

      clustered_idx <- which(is_in_same_cluster)

      left3 <- pmax(center_left[t], center_left[clustered_idx])
      right3 <- pmin(center_right[t], center_right[clustered_idx])

      if (!any(right3 - left3 >= 0, na.rm = TRUE)) {next}
      
  ## updata the interval of putative IZ by the blur center of a short track
      if (width_bins[t] <= short_track_th_bins) {
        idx2 <- c(clustered_idx, t)

        putative_IZs$start[overlapped_cluster] <- min(
          center_left[idx2],
          na.rm = TRUE
        )

        putative_IZs$end[overlapped_cluster] <- max(
          center_right[idx2],
          na.rm = TRUE
        )
      }
      cluster_id1[t] <- overlapped_cluster
    }

    if (nrow(putative_IZs) == 0L) {
      return(list(
        chr = chr,
        assignment = data.frame(
          original_index = S4Vectors::mcols(gr_chr)$.original_index,
          local_iz_id = rep(NA_integer_, n_tracks)
        ),
        chr_izs = GenomicRanges::GRanges()
      ))
    }

  ## filter and re-order putative IZs
    putative_IZs$old_id <- seq_len(nrow(putative_IZs))
    
    putative_IZs$n_tracks <- tabulate(
      cluster_id1[!is.na(cluster_id1)],
      nbins = nrow(putative_IZs)
    )
    
    keep_iz <- putative_IZs$n_tracks >= min_usage
    putative_IZs <- putative_IZs[keep_iz, , drop = FALSE]
    putative_IZs <- putative_IZs[order(putative_IZs$start), , drop = FALSE]
    
    if (nrow(putative_IZs) == 0L) {
      return(list(
        chr = chr,
        assignment = data.frame(
          original_index = S4Vectors::mcols(gr_chr)$.original_index,
          local_iz_id = rep(NA_integer_, n_tracks)
        ),
        chr_izs = GenomicRanges::GRanges()
      ))
    }

    ## 2nd loop to re-asign tracks to new IZ IDs
      n_putative_IZs <- nrow(putative_IZs)
    cluster_id2 <- rep(NA_integer_, n_tracks)
    
    for (t in seq_len(n_tracks)) {
        left1 <- pmax(blur_start[t], putative_IZs$start)
        right1 <- pmin(blur_end[t], putative_IZs$end)

        overlap_threshold <- if (width_bins[t] <= short_track_th_bins) {
          min_overlap
        } else {
          min_overlap + bin_size
        }

        is_overlapping <- right1 - left1 >= overlap_threshold
        num_overlaps <- sum(is_overlapping, na.rm = TRUE)

        if (num_overlaps == 0L || num_overlaps > 2L) { next }

        if (num_overlaps == 2L) {
          candidate_clusters <- which(is_overlapping)

          left2 <- pmax(center_left[t], putative_IZs$start[candidate_clusters])
          right2 <- pmin(center_right[t], putative_IZs$end[candidate_clusters])
          center_compatible <- right2 - left2 >= 0
          names(center_compatible) <- candidate_clusters

          if (sum(center_compatible, na.rm = TRUE) == 1L) {

              overlapped_cluster <- candidate_clusters[which(center_compatible)]
              cluster_id2[t] <- overlapped_cluster
          }
        } else { # here num_overlaps == 1L
          overlapped_cluster <- which(is_overlapping)
            
          left3 <- max(center_left[t], putative_IZs$start[overlapped_cluster])
          right3 <- min(center_right[t], putative_IZs$end[overlapped_cluster])
            
          if (right3 - left3 >= 0) { 
              cluster_id2[t] <- overlapped_cluster }
        }        
    }

    n_support <- tabulate(
      cluster_id2[!is.na(cluster_id2)],
      nbins = nrow(putative_IZs)
    )

    keep_iz <- n_support >= min_usage

    if (!any(keep_iz)) {
      return(list(
        chr = chr,
        assignment = data.frame(
          original_index = S4Vectors::mcols(gr_chr)$.original_index,
          local_iz_id = rep(NA_integer_, n_tracks)
        ),
        chr_izs = GenomicRanges::GRanges()
      ))
    }

    putative_IZs <- putative_IZs[keep_iz, , drop = FALSE]
    old_to_new <- rep(NA_integer_, length(keep_iz))
    old_to_new[which(keep_iz)] <- seq_len(sum(keep_iz))
    cluster_id3 <- old_to_new[cluster_id2]

    num_izs <- nrow(putative_IZs)
    chr_iz_list <- vector("list", num_izs)

    for (id in seq_len(num_izs)) {
      idx <- which(cluster_id3 == id)
      centers <- center[idx]

      qs <- as.numeric(stats::quantile(
        centers,
        probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
        na.rm = TRUE,
        names = FALSE,
        type = 7
      ))

      names(qs) <- c("q5", "q25", "q50", "q75", "q95")
      iz_width <- qs["q95"] - qs["q5"] + bin_size  #iz_width <- qs["q75"] - qs["q25"] + bin_size
      # q5 and q95 locate the centers of the two edge bins. Pad the
      # q5-q95 interval by half a bin on each side instead of centering the
      # full IZ width on q5, which shifts every IZ to the left.
      iz_start <- round(qs["q5"] - bin_size / 2)
      iz_start <- max(1L, as.integer(iz_start))
      iz_end <- iz_start + round(iz_width) - 1L
      
      iz_gr <- GenomicRanges::GRanges(
        seqnames = chr,
        ranges = IRanges::IRanges(start = iz_start, end = iz_end),
        strand = "*"
      )

      S4Vectors::mcols(iz_gr)$local_iz_id <- id
      S4Vectors::mcols(iz_gr)$`q5` <- qs["q5"]
      S4Vectors::mcols(iz_gr)$`q25` <- qs["q25"]
      S4Vectors::mcols(iz_gr)$`q50` <- qs["q50"]
      S4Vectors::mcols(iz_gr)$`q75` <- qs["q75"]
      S4Vectors::mcols(iz_gr)$`q95` <- qs["q95"]
      S4Vectors::mcols(iz_gr)$n_tracks <- length(idx)

      chr_iz_list[[id]] <- iz_gr
    }

    chr_izs <- do.call(c, chr_iz_list)
    GenomeInfoDb::seqlevels(chr_izs, pruning.mode = "coarse") <- chrs

    assignment <- data.frame(
      original_index = S4Vectors::mcols(gr_chr)$.original_index,
      local_iz_id = cluster_id3
    )

    list(
      chr = chr,
      assignment = assignment,
      chr_izs = chr_izs
    )
  }

  res <- BiocParallel::bplapply(
    chrs,
    process_one_chr,
    BPPARAM = BPPARAM
  )

  .vmsg(verbose, "Combining chromosome-level results.")

  izs_list <- lapply(res, `[[`, "chr_izs")
  n_izs <- vapply(izs_list, length, integer(1))
  offsets <- cumsum(c(0L, utils::head(n_izs, -1L)))

  assignment_list <- vector("list", length(res))

  for (i in seq_along(res)) {
    iz_gr <- izs_list[[i]]
    offset <- offsets[i]

    if (length(iz_gr) > 0L) {
      local_id <- S4Vectors::mcols(iz_gr)$local_iz_id
      global_id <- local_id + offset

      S4Vectors::mcols(iz_gr)$iz_id <- global_id
      S4Vectors::mcols(iz_gr)$local_iz_id <- NULL

      izs_list[[i]] <- iz_gr
    }

    assignment <- res[[i]]$assignment

    if (nrow(assignment) > 0L) {
      assignment$iz_id <- ifelse(
        is.na(assignment$local_iz_id),
        NA_integer_,
        assignment$local_iz_id + offset
      )
    } else {
      assignment$iz_id <- integer(0)
    }

    assignment_list[[i]] <- assignment[, c("original_index", "iz_id"), drop = FALSE]
  }

  izones <- do.call(c, izs_list)

  if (is.null(izones)) {
    izones <- GenomicRanges::GRanges()
  }

  assignment_df <- do.call(rbind, assignment_list)

  iz_id_all <- rep(NA_integer_, length(gr_all))

  if (!is.null(assignment_df) && nrow(assignment_df) > 0L) {
    iz_id_all[assignment_df$original_index] <- assignment_df$iz_id
  }

  S4Vectors::mcols(gr_all)$iz_id <- iz_id_all
  S4Vectors::mcols(gr_all)$.cell_id <- NULL
  S4Vectors::mcols(gr_all)$.original_index <- NULL

  if (length(izones) > 0L) {
      izones <- BiocGenerics::sort(izones)
      old_iz_id <- S4Vectors::mcols(izones)$iz_id
      new_iz_id <- seq_along(izones)
      
      id_map <- new_iz_id
      names(id_map) <- old_iz_id
      S4Vectors::mcols(izones)$iz_id <- new_iz_id
      
      keep_cols <- c("iz_id", "q5", "q25", "q50", "q75", "q95", "n_tracks")
      S4Vectors::mcols(izones) <- S4Vectors::mcols(izones)[, keep_cols, drop = FALSE]
      
      old_ids_on_tracks <- S4Vectors::mcols(gr_all)$iz_id
      remapped <- rep(NA_integer_, length(old_ids_on_tracks))
      non_na <- !is.na(old_ids_on_tracks)
      remapped[non_na] <- id_map[as.character(old_ids_on_tracks[non_na])]
      S4Vectors::mcols(gr_all)$iz_id <- remapped
      }
    
    S4Vectors::mcols(gr_all)$.cell_id <- NULL
    S4Vectors::mcols(gr_all)$.original_index <- NULL
    out_tracks <- IRanges::relist(gr_all, tracks)
    
    .vmsg(verbose, "Done. Called ", length(izones), " initiation zone(s).")
    
    list(
        tracks = out_tracks,
        izones = izones
        )
}

#' Firing range of initiation zones
#'
#' For each IZ with at least `min_rep_tracks` assigned tracks, the earliest
#' firing time is the smallest `rep_percent` among cells with a track assigned
#' to the IZ, and the latest unreplicated time is the largest `rep_percent`
#' among cells whose unreplicated segment (run of bins in state
#' `unrep_state`) covers the IZ centre, using only the segments that cover the
#' fewest IZ centres. Each bound is corroborated when the second-best value
#' is within 0.1 of it. IZs with `latest_unrep_percent <
#' earliest_rep_percent` get `NA` bounds and `corroborated = FALSE`.
#'
#' @param gr_izones `GRanges` of IZs with an `n_tracks` column (e.g.
#'   `define_initiation_zones()$izones`); IZ `i` must have `iz_id == i`.
#' @param genobins_rep_states Matrix of states, bins x cells, with row names
#'   `"chr-start-end"` and cell IDs as column names; every IZ centre must lie
#'   in exactly one bin.
#' @param replication_tracks `CompressedGRangesList` of tracks with metadata
#'   columns `iz_id`, `cell_id` and `rep_percent` (e.g.
#'   `define_initiation_zones()$tracks`); every matrix column must have a
#'   `rep_percent`.
#' @param unrep_state State of unreplicated bins.
#' @param min_rep_tracks Minimum `n_tracks` of an IZ.
#' @param verbose Logical; print progress messages. Default `FALSE`.
#'
#' @return `gr_izones` with columns `earliest_rep_percent`,
#'   `latest_unrep_percent` and `corroborated`.
#'
#' @export
calculate_firing_range <- function(
  gr_izones,
  genobins_rep_states,
  replication_tracks,
  unrep_state = 2,
  min_rep_tracks = 2,
  verbose = FALSE
) {
  if (!inherits(gr_izones, "GRanges")) {
    stop("'gr_izones' must be a GRanges object.")
  }

  if (!"n_tracks" %in% names(S4Vectors::mcols(gr_izones))) {
    stop("'mcols(gr_izones)' must contain column 'n_tracks'.")
  }

  if (!is.matrix(genobins_rep_states)) {
    stop("'genobins_rep_states' must be a matrix.")
  }

  if (is.null(rownames(genobins_rep_states))) {
    stop("'genobins_rep_states' must have row names in format 'chr-start-end'.")
  }

  if (is.null(colnames(genobins_rep_states))) {
    stop("'genobins_rep_states' must have cell IDs as column names.")
  }

  if (!inherits(replication_tracks, "CompressedGRangesList")) {
    stop("'replication_tracks' must be a CompressedGRangesList.")
  }

  rep_tracks_unlisted <- methods::slot(replication_tracks, "unlistData")
  track_mcols <- S4Vectors::mcols(rep_tracks_unlisted)

  required_track_cols <- c("iz_id", "cell_id", "rep_percent")
  missing_track_cols <- setdiff(required_track_cols, names(track_mcols))

  if (length(missing_track_cols) > 0L) {
    stop(
      "'mcols(replication_tracks@unlistData)' must contain column(s): ",
      paste(missing_track_cols, collapse = ", ")
    )
  }

  if (anyNA(track_mcols$cell_id) || any(as.character(track_mcols$cell_id) == "")) {
    stop(
      "'mcols(replication_tracks@unlistData)$cell_id' contains missing or empty values."
    )
  }

  if (!is.numeric(track_mcols$rep_percent)) {
    stop("'mcols(replication_tracks@unlistData)$rep_percent' must be numeric.")
  }

  track_cell_ids <- as.character(track_mcols$cell_id)
  track_rep_percent <- track_mcols$rep_percent

  cell_rep_df <- unique(data.frame(
    cell_id = track_cell_ids,
    rep_percent = track_rep_percent,
    stringsAsFactors = FALSE
  ))

  non_na_cell_rep_df <- cell_rep_df[
    !is.na(cell_rep_df$rep_percent),
    ,
    drop = FALSE
  ]

  duplicated_cell_rep <- split(
    non_na_cell_rep_df$rep_percent,
    non_na_cell_rep_df$cell_id
  )

  inconsistent_cells <- names(duplicated_cell_rep)[
    vapply(
      duplicated_cell_rep,
      function(x) length(unique(x)) > 1L,
      logical(1)
    )
  ]

  if (length(inconsistent_cells) > 0L) {
    stop(
      "Some cells have inconsistent 'rep_percent' values in ",
      "'mcols(replication_tracks@unlistData)': ",
      paste(inconsistent_cells, collapse = ", ")
    )
  }

  named_rep_percent <- tapply(
    non_na_cell_rep_df$rep_percent,
    non_na_cell_rep_df$cell_id,
    function(x) x[1]
  )

  missing_bin_cells <- setdiff(
    colnames(genobins_rep_states),
    names(named_rep_percent)
  )

  if (length(missing_bin_cells) > 0L) {
    stop(
      "Some cells in 'genobins_rep_states' are missing from ",
      "'mcols(replication_tracks@unlistData)$cell_id' or have NA 'rep_percent'."
    )
  }

  bin_ids <- rownames(genobins_rep_states)

  parsed_bins <- utils::strcapture(
    pattern = "^(.*)-([0-9]+)-([0-9]+)$",
    x = bin_ids,
    proto = data.frame(
      chr = character(),
      start = integer(),
      end = integer(),
      stringsAsFactors = FALSE
    )
  )

  if (
    anyNA(parsed_bins$chr) ||
      anyNA(parsed_bins$start) ||
      anyNA(parsed_bins$end)
  ) {
    stop("Failed to parse bin row names. Expected format: 'chr-start-end'.")
  }

  if (any(parsed_bins$start < 1L) || any(parsed_bins$end < 1L)) {
    stop("Bin start and end coordinates must be positive integers.")
  }

  if (any(parsed_bins$end < parsed_bins$start)) {
    stop("Some bins have end < start.")
  }

  bin_gr <- GenomicRanges::GRanges(
    seqnames = parsed_bins$chr,
    ranges = IRanges::IRanges(
      start = parsed_bins$start,
      end = parsed_bins$end
    )
  )

  iz_centers <- GenomicRanges::resize(
    gr_izones,
    width = 1L,
    fix = "center"
  )

  n_iz <- length(gr_izones)
  n_tracks <- S4Vectors::mcols(gr_izones)$n_tracks

  if (length(n_tracks) != n_iz) {
    stop("'mcols(gr_izones)$n_tracks' must have length equal to length(gr_izones).")
  }

  earliest_rep_percent <- rep(NA_real_, n_iz)
  latest_unrep_percent <- rep(NA_real_, n_iz)

  corroborated_earliest <- rep(FALSE, n_iz)
  corroborated_latest <- rep(FALSE, n_iz)

  rep_iz_ids_all <- track_mcols$iz_id

  if (length(rep_iz_ids_all) != length(rep_tracks_unlisted)) {
    stop("'iz_id' length does not match the number of unlisted replication tracks.")
  }

  if (length(track_cell_ids) != length(rep_tracks_unlisted)) {
    stop(
      "'cell_id' length does not match the number of unlisted replication tracks."
    )
  }

  assigned_rep_track_idx <- which(!is.na(rep_iz_ids_all))

  rep_iz_ids <- rep_iz_ids_all[assigned_rep_track_idx]

  ## Important:
  ## Use mcols(replication_tracks@unlistData)$cell_id directly.
  ## Do not use names(replication_tracks).
  rep_cell_ids_per_track <- track_cell_ids[assigned_rep_track_idx]

  if (length(rep_iz_ids) == 0L) {
    warning(
      "No replication tracks have non-NA 'iz_id'. ",
      "Earliest firing timing will be NA for all IZs."
    )
  }

  if (length(rep_iz_ids) > 0L) {
    if (any(rep_iz_ids < 1L | rep_iz_ids > n_iz)) {
      stop("Non-NA 'iz_id' values must be valid 1-based indices into 'gr_izones'.")
    }

    if (any(rep_iz_ids != round(rep_iz_ids))) {
      stop("Non-NA 'iz_id' values must be integer-valued 1-based indices.")
    }

    rep_iz_ids <- as.integer(rep_iz_ids)
  }

  all_chrs <- unique(as.character(GenomeInfoDb::seqnames(gr_izones)))

  for (chr in all_chrs) {
    .vmsg(verbose, "Processing ", chr)

    iz_idx_chr <- which(as.character(GenomeInfoDb::seqnames(gr_izones)) == chr)
    bin_idx_chr <- which(as.character(GenomeInfoDb::seqnames(bin_gr)) == chr)

    if (length(iz_idx_chr) == 0L || length(bin_idx_chr) == 0L) {
      next
    }

    iz_centers_chr <- iz_centers[iz_idx_chr]
    bin_gr_chr <- bin_gr[bin_idx_chr]

    hits_center_bin <- GenomicRanges::findOverlaps(
      query = iz_centers_chr,
      subject = bin_gr_chr,
      type = "within"
    )

    n_hits <- tabulate(
      S4Vectors::queryHits(hits_center_bin),
      nbins = length(iz_idx_chr)
    )

    if (any(n_hits == 0L)) {
      stop("Some IZ centers on ", chr, " do not overlap any genomic bin.")
    }

    if (any(n_hits > 1L)) {
      stop("Some IZ centers on ", chr, " overlap more than one genomic bin.")
    }

    chr_states <- genobins_rep_states[bin_idx_chr, , drop = FALSE]

    unrep_tracks_chr <- .find_unreplicated_segments(
      states = chr_states,
      bin_gr = bin_gr_chr,
      cell_ids = colnames(genobins_rep_states),
      unrep_state = unrep_state
    )

    if (nrow(unrep_tracks_chr) > 0L) {
      unrep_gr <- GenomicRanges::GRanges(
        seqnames = chr,
        ranges = IRanges::IRanges(
          start = unrep_tracks_chr$start,
          end = unrep_tracks_chr$end
        )
      )

      hits_unrep_iz <- GenomicRanges::findOverlaps(
        query = unrep_gr,
        subject = iz_centers_chr,
        type = "any"
      )

      n_izs_per_unrep_track <- tabulate(
        S4Vectors::queryHits(hits_unrep_iz),
        nbins = length(unrep_gr)
      )

      unrep_tracks_chr$n_overlapped_izs <- n_izs_per_unrep_track
    }

    for (ii in seq_along(iz_idx_chr)) {
      iz_global_idx <- iz_idx_chr[ii]

      if (
        is.na(n_tracks[iz_global_idx]) ||
          n_tracks[iz_global_idx] < min_rep_tracks
      ) {
        next
      }

      rep_track_idx <- which(rep_iz_ids == iz_global_idx)

      if (length(rep_track_idx) < n_tracks[iz_global_idx]) {
        stop(
          "IZ ", iz_global_idx,
          " has fewer non-NA assigned tracks in 'replication_tracks' than ",
          "mcols(gr_izones)$n_tracks indicates."
        )
      }

      rep_cells <- rep_cell_ids_per_track[rep_track_idx]
      rep_percent <- named_rep_percent[rep_cells]
      rep_percent <- sort(rep_percent[!is.na(rep_percent)])

      if (length(rep_percent) >= 1L) {
        earliest_rep_percent[iz_global_idx] <- rep_percent[1]

        if (length(rep_percent) >= 2L) {
          corroborated_earliest[iz_global_idx] <-
            rep_percent[2] - rep_percent[1] < 0.1
        }
      }

      if (nrow(unrep_tracks_chr) == 0L) {
        next
      }

      center_pos <- BiocGenerics::start(iz_centers_chr[ii])

      overlaps_this_iz <- (
        unrep_tracks_chr$start <= center_pos &
          unrep_tracks_chr$end >= center_pos
      )

      unrep_tracks_this_iz <- unrep_tracks_chr[
        overlaps_this_iz,
        ,
        drop = FALSE
      ]

      if (nrow(unrep_tracks_this_iz) == 0L) {
        next
      }

      min_num_izs <- min(unrep_tracks_this_iz$n_overlapped_izs)

      unrep_tracks_this_iz <- unrep_tracks_this_iz[
        unrep_tracks_this_iz$n_overlapped_izs == min_num_izs,
        ,
        drop = FALSE
      ]

      unrep_cells <- unrep_tracks_this_iz$cell_id
      unrep_percent <- named_rep_percent[unrep_cells]
      unrep_percent <- sort(unrep_percent[!is.na(unrep_percent)])

      if (length(unrep_percent) >= 1L) {
        latest_unrep_percent[iz_global_idx] <-
          unrep_percent[length(unrep_percent)]

        if (length(unrep_percent) > 1L) {
          corroborated_latest[iz_global_idx] <-
            unrep_percent[length(unrep_percent)] -
            unrep_percent[length(unrep_percent) - 1L] < 0.1
        }
      }
    }
  }

  is_imbalanced <- latest_unrep_percent < earliest_rep_percent
  is_imbalanced[is.na(is_imbalanced)] <- FALSE

  earliest_rep_percent[is_imbalanced] <- NA_real_
  latest_unrep_percent[is_imbalanced] <- NA_real_

  corroborated <- corroborated_earliest & corroborated_latest
  corroborated[is_imbalanced] <- FALSE

  IZ_range <- gr_izones
  S4Vectors::mcols(IZ_range)$earliest_rep_percent <- earliest_rep_percent
  S4Vectors::mcols(IZ_range)$latest_unrep_percent <- latest_unrep_percent
  S4Vectors::mcols(IZ_range)$corroborated <- corroborated

  return(IZ_range)
}


#' Runs of unreplicated bins per cell on one chromosome
#' @keywords internal
#' @noRd
.find_unreplicated_segments <- function(
  states,
  bin_gr,
  cell_ids,
  unrep_state = 2
) {
  out <- vector("list", length(cell_ids))

  for (j in seq_along(cell_ids)) {
    is_unrep <- states[, j] == unrep_state
    is_unrep[is.na(is_unrep)] <- FALSE

    r <- rle(is_unrep)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1L

    true_runs <- which(r$values)

    if (length(true_runs) == 0L) {
      out[[j]] <- NULL
      next
    }

    seg_start_idx <- starts[true_runs]
    seg_end_idx <- ends[true_runs]

    out[[j]] <- data.frame(
      cell_id = cell_ids[j],
      start = BiocGenerics::start(bin_gr)[seg_start_idx],
      end = BiocGenerics::end(bin_gr)[seg_end_idx],
      start_bin = seg_start_idx,
      end_bin = seg_end_idx,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, out)

  if (is.null(out)) {
    out <- data.frame(
      cell_id = character(),
      start = integer(),
      end = integer(),
      start_bin = integer(),
      end_bin = integer(),
      stringsAsFactors = FALSE
    )
  }

  rownames(out) <- NULL
  out
}

#' Firing-order consistency of initiation zones from replication states
#'
#' IZs are ordered by decreasing `rt` (earliest first). For each cell with
#' `rep_percent` in `rep_range` (inclusive), the state of the bin containing
#' each IZ centre is read (cells with less than 90% non-`NA` IZ states are set
#' to `NA`); if `K` IZs are fired (state `4`), the first `K` IZs in timing
#' order are expected to be fired. Observed (`4` fired, `2` unfired) and
#' expected states are cross-tabulated per IZ.
#'
#' @param izs `GRanges` of IZs with a numeric `rt` column (no `NA`).
#' @param genobins_rep_states Matrix of states, bins x cells, with row names
#'   `"chr-start-end"`; every IZ centre must lie in exactly one bin.
#' @param meta_data A `data.frame` with columns `cell_id` (covering all matrix
#'   columns) and `rep_percent` in `[0, 1]`.
#' @param rep_range Inclusive `rep_percent` range of the cells used.
#'
#' @return A numeric matrix, one row per IZ in timing order (row names
#'   `"chr:start-end"`), with the percentages
#'   `observed_fired_expected_fired`, `observed_unfired_expected_unfired`,
#'   `observed_fired_expected_unfired` and `observed_unfired_expected_fired`
#'   (`NA` if an IZ has no observation).
#'
#' @export
calculate_firing_order <- function(
  izs,
  genobins_rep_states,
  meta_data,
  rep_range = c(0, 1)
) {

  if (!inherits(izs, "GRanges")) {
    stop("'izs' must be a GRanges object.")
  }

  if (!is.matrix(genobins_rep_states)) {
    stop("'genobins_rep_states' must be a matrix.")
  }

  if (is.null(rownames(genobins_rep_states))) {
    stop("'genobins_rep_states' must have rownames in the format 'chr-start-end'.")
  }

  if (is.null(colnames(genobins_rep_states))) {
    stop("'genobins_rep_states' must have colnames as cell IDs.")
  }

  if (!is.data.frame(meta_data)) {
    stop("'meta_data' must be a data frame.")
  }

  if (!all(c("cell_id", "rep_percent") %in% colnames(meta_data))) {
    stop("'meta_data' must contain columns 'cell_id' and 'rep_percent'.")
  }

  if (anyDuplicated(meta_data$cell_id)) {
    stop("'meta_data$cell_id' contains duplicated cell IDs.")
  }

  if (!is.numeric(meta_data$rep_percent)) {
    stop("'meta_data$rep_percent' must be numeric.")
  }

  if (any(!is.na(meta_data$rep_percent) &
          (meta_data$rep_percent < 0 | meta_data$rep_percent > 1))) {
    stop("'meta_data$rep_percent' must be within [0, 1], excluding NA values.")
  }

  if (!is.numeric(rep_range) || length(rep_range) != 2L || anyNA(rep_range)) {
    stop("'rep_range' must be a numeric vector of length 2 without NA.")
  }

  if (any(rep_range < 0 | rep_range > 1)) {
    stop("'rep_range' must be within [0, 1].")
  }

  if (rep_range[1] > rep_range[2]) {
    stop("'rep_range[1]' must be <= 'rep_range[2]'.")
  }

  if (!all(colnames(genobins_rep_states) %in% meta_data$cell_id)) {
    missing_cells <- setdiff(colnames(genobins_rep_states), meta_data$cell_id)
    stop(
      "Some cells in 'genobins_rep_states' are missing from 'meta_data'. ",
      "First missing cell: ", missing_cells[1], "."
    )
  }

  meta_data <- meta_data[
    match(colnames(genobins_rep_states), meta_data$cell_id),
    ,
    drop = FALSE
  ]

  selected_cells <- !is.na(meta_data$rep_percent) &
    meta_data$rep_percent >= rep_range[1] &
    meta_data$rep_percent <= rep_range[2]

  if (!any(selected_cells)) {
    stop(
      "No cells found within the specified 'rep_range' after excluding ",
      "cells with NA 'rep_percent'."
    )
  }

  genobins_rep_states <- genobins_rep_states[, selected_cells, drop = FALSE]

  if (!"rt" %in% names(S4Vectors::mcols(izs))) {
    stop("'izs' must contain metadata column 'rt'.")
  }

  if (length(izs) == 0L) {
    stop("'izs' contains no initiation zones.")
  }

  if (nrow(genobins_rep_states) == 0L) {
    stop("'genobins_rep_states' contains no genomic bins.")
  }

  bin_ids <- rownames(genobins_rep_states)

  parsed_bins <- utils::strcapture(
    pattern = "^(.*)-([0-9]+)-([0-9]+)$",
    x = bin_ids,
    proto = data.frame(
      chr = character(),
      start = integer(),
      end = integer(),
      stringsAsFactors = FALSE
    )
  )

  if (
    anyNA(parsed_bins$chr) ||
      anyNA(parsed_bins$start) ||
      anyNA(parsed_bins$end)
  ) {
    stop(
      "Failed to parse rownames(genobins_rep_states). ",
      "Expected format: 'chr-start-end', for example 'chr1-100001-200000'."
    )
  }

  if (any(parsed_bins$start < 1L) || any(parsed_bins$end < 1L)) {
    stop("Bin start and end coordinates must be positive integers.")
  }

  if (any(parsed_bins$end < parsed_bins$start)) {
    stop("Some bins have end < start.")
  }

  genobin_gr <- GenomicRanges::GRanges(
    seqnames = parsed_bins$chr,
    ranges = IRanges::IRanges(
      start = parsed_bins$start,
      end = parsed_bins$end
    )
  )

  iz_centers <- GenomicRanges::resize(
    izs,
    width = 1L,
    fix = "center"
  )

  hits <- GenomicRanges::findOverlaps(
    query = iz_centers,
    subject = genobin_gr,
    type = "within",
    select = "all"
  )

  query_hits <- S4Vectors::queryHits(hits)
  subject_hits <- S4Vectors::subjectHits(hits)

  n_hits_per_iz <- tabulate(query_hits, nbins = length(izs))

  if (any(n_hits_per_iz == 0L)) {
    missing_izs <- which(n_hits_per_iz == 0L)
    stop(
      "Some IZ centers do not overlap any genomic bin. ",
      "First missing IZ index: ", missing_izs[1], "."
    )
  }

  if (any(n_hits_per_iz > 1L)) {
    duplicated_izs <- which(n_hits_per_iz > 1L)
    stop(
      "Some IZ centers overlap more than one genomic bin. ",
      "First ambiguous IZ index: ", duplicated_izs[1], "."
    )
  }

  bin_index_for_iz <- integer(length(izs))
  bin_index_for_iz[query_hits] <- subject_hits

  iz_state <- genobins_rep_states[bin_index_for_iz, , drop = FALSE]

  rt <- S4Vectors::mcols(izs)$rt

  if (anyNA(rt)) {
    stop("'mcols(izs)$rt' contains NA values.")
  }

  sort_order <- order(rt, decreasing = TRUE)

  iz_state <- iz_state[sort_order, , drop = FALSE]
  izs_ordered <- izs[sort_order]

  rownames(iz_state) <- paste0(
    as.character(GenomeInfoDb::seqnames(izs_ordered)),
    ":",
    BiocGenerics::start(izs_ordered),
    "-",
    BiocGenerics::end(izs_ordered)
  )

  completeness <- colSums(!is.na(iz_state)) / nrow(iz_state)

  low_quality_cells <- completeness < 0.9

  if (any(low_quality_cells)) {
    iz_state[, low_quality_cells] <- NA
  }

  num_fired_izs <- colSums(iz_state == 4, na.rm = TRUE)

  expected_iz_state <- matrix(
    FALSE,
    nrow = nrow(iz_state),
    ncol = ncol(iz_state),
    dimnames = dimnames(iz_state)
  )

  for (cell_idx in seq_len(ncol(iz_state))) {
    k <- num_fired_izs[cell_idx]

    if (!is.na(k) && k > 0L) {
      expected_iz_state[seq_len(k), cell_idx] <- TRUE
    }
  }

  observed_fired_expected_fired <- rowSums(
    iz_state == 4 & expected_iz_state,
    na.rm = TRUE
  )

  observed_unfired_expected_unfired <- rowSums(
    iz_state == 2 & !expected_iz_state,
    na.rm = TRUE
  )

  observed_fired_expected_unfired <- rowSums(
    iz_state == 4 & !expected_iz_state,
    na.rm = TRUE
  )

  observed_unfired_expected_fired <- rowSums(
    iz_state == 2 & expected_iz_state,
    na.rm = TRUE
  )

  firing_order_counts <- cbind(
    observed_fired_expected_fired = observed_fired_expected_fired,
    observed_unfired_expected_unfired = observed_unfired_expected_unfired,
    observed_fired_expected_unfired = observed_fired_expected_unfired,
    observed_unfired_expected_fired = observed_unfired_expected_fired
  )

  denominator <- rowSums(firing_order_counts)

  firing_order_frequency <- firing_order_counts / denominator * 100

  firing_order_frequency[denominator == 0, ] <- NA_real_

  rownames(firing_order_frequency) <- rownames(iz_state)

  firing_order_frequency
}

#' Firing-order consistency of initiation zones from track centres
#'
#' Like [calculate_firing_order()], but an IZ is fired in a cell when at least
#' one of the cell's track centres (`floor((start + end) / 2)`) falls inside
#' it. Cells are selected by `rep_range[1] <= rep_percent < rep_range[2]` and
#' tracks optionally filtered with `track_filter`.
#'
#' @param gr_izones `GRanges` of IZs with a numeric `rt` column (no `NA`).
#' @param rep_tracks Named `GRangesList` of tracks with metadata columns
#'   `cell_id` (equal to the element name), `rep_percent` and `num_bins`.
#' @param rep_range `rep_percent` range (left-closed, right-open).
#' @param track_filter Optional unquoted expression evaluated against each
#'   cell's track metadata (non-standard evaluation), e.g.
#'   `num_bins <= 10`.
#'
#' @return A numeric matrix as in [calculate_firing_order()].
#'
#' @export
predict_firing_order <- function(
  gr_izones,
  rep_tracks,
  rep_range = c(0, 1),
  track_filter = NULL
) {
  if (!inherits(gr_izones, "GRanges")) {
    stop("'gr_izones' must be a GenomicRanges::GRanges object.")
  }

  if (!inherits(rep_tracks, "GRangesList")) {
    stop("'rep_tracks' must be a GRangesList or CompressedGRangesList.")
  }

  if (length(gr_izones) == 0L) {
    stop("'gr_izones' contains no initiation zones.")
  }

  if (!"rt" %in% names(S4Vectors::mcols(gr_izones))) {
    stop("'gr_izones' must contain metadata column 'rt'.")
  }

  rt <- S4Vectors::mcols(gr_izones)$rt

  if (!is.numeric(rt)) {
    stop("'mcols(gr_izones)$rt' must be numeric.")
  }

  if (anyNA(rt)) {
    stop("'mcols(gr_izones)$rt' contains NA values.")
  }

  if (length(rep_tracks) == 0L) {
    stop("'rep_tracks' contains no cells.")
  }

  if (!is.numeric(rep_range) || length(rep_range) != 2L) {
    stop("'rep_range' must be a numeric vector of length 2.")
  }

  if (anyNA(rep_range)) {
    stop("'rep_range' must not contain NA values.")
  }

  if (rep_range[1] >= rep_range[2]) {
    stop("'rep_range[1]' must be smaller than 'rep_range[2]'.")
  }

  cell_ids <- names(rep_tracks)

  if (is.null(cell_ids) || anyNA(cell_ids) || any(cell_ids == "")) {
    stop("'names(rep_tracks)' must contain non-missing cell IDs.")
  }

  if (anyDuplicated(cell_ids)) {
    stop("'names(rep_tracks)' must contain unique cell IDs.")
  }

  rep_tracks_unlisted <- BiocGenerics::unlist(rep_tracks, use.names = FALSE)
  rep_mcols <- S4Vectors::mcols(rep_tracks_unlisted)

  required_track_cols <- c("cell_id", "rep_percent", "num_bins")
  missing_track_cols <- setdiff(required_track_cols, names(rep_mcols))

  if (length(missing_track_cols) > 0L) {
    stop(
      "'mcols(rep_tracks@unlistData)' must contain column(s): ",
      paste(missing_track_cols, collapse = ", ")
    )
  }

  if (length(rep_tracks_unlisted) > 0L) {
    if (anyNA(rep_mcols$cell_id) || any(as.character(rep_mcols$cell_id) == "")) {
      stop("'mcols(rep_tracks@unlistData)$cell_id' contains missing or empty values.")
    }

    if (!is.numeric(rep_mcols$rep_percent)) {
      stop("'mcols(rep_tracks@unlistData)$rep_percent' must be numeric.")
    }

    if (!is.numeric(rep_mcols$num_bins)) {
      stop("'mcols(rep_tracks@unlistData)$num_bins' must be numeric.")
    }
  }

  ## ------------------------------------------------------------
  ## Filter cells by rep_range using metadata columns:
  ## cell_id and rep_percent in rep_tracks@unlistData
  ## ------------------------------------------------------------

  rep_df <- unique(data.frame(
    cell_id = as.character(rep_mcols$cell_id),
    rep_percent = rep_mcols$rep_percent,
    stringsAsFactors = FALSE
  ))

  rep_percent_by_cell <- split(rep_df$rep_percent, rep_df$cell_id)

  invalid_rep_percent_cells <- names(rep_percent_by_cell)[
    vapply(
      rep_percent_by_cell,
      function(x) length(unique(x[!is.na(x)])) != 1L,
      logical(1)
    )
  ]

  if (length(invalid_rep_percent_cells) > 0L) {
    stop(
      "Each cell must have exactly one non-missing 'rep_percent'. ",
      "Problematic cell(s): ",
      paste(invalid_rep_percent_cells, collapse = ", ")
    )
  }

  cell_rep_percent <- vapply(
    rep_percent_by_cell,
    function(x) unique(x[!is.na(x)])[1],
    numeric(1)
  )

  selected_cell_ids <- names(cell_rep_percent)[
    cell_rep_percent >= rep_range[1] &
      cell_rep_percent < rep_range[2]
  ]

  selected_cell_ids <- intersect(cell_ids, selected_cell_ids)

  if (length(selected_cell_ids) == 0L) {
    stop(
      "No cells remain after filtering by 'rep_range = c(",
      rep_range[1], ", ", rep_range[2], ")'."
    )
  }

  rep_tracks <- rep_tracks[selected_cell_ids]
  cell_ids <- names(rep_tracks)

  filter_expr <- substitute(track_filter)
  has_track_filter <- !identical(filter_expr, quote(NULL))

  sort_order <- order(rt, decreasing = TRUE)
  izones_ordered <- gr_izones[sort_order]

  izone_names <- paste0(
    as.character(GenomicRanges::seqnames(izones_ordered)),
    ":",
    BiocGenerics::start(izones_ordered),
    "-",
    BiocGenerics::end(izones_ordered)
  )

  observed_fired <- matrix(
    FALSE,
    nrow = length(izones_ordered),
    ncol = length(rep_tracks),
    dimnames = list(izone_names, cell_ids)
  )

  for (i in seq_along(rep_tracks)) {
    gr <- rep_tracks[[i]]
    cell_id <- cell_ids[i]

    if (length(gr) == 0L) {
      next
    }

    gr_mcols <- S4Vectors::mcols(gr)

    if (!"cell_id" %in% names(gr_mcols)) {
      stop("'mcols(rep_tracks[[", i, "]])' must contain column 'cell_id'.")
    }

    gr_cell_ids <- unique(as.character(gr_mcols$cell_id))

    if (length(gr_cell_ids) != 1L || is.na(gr_cell_ids) || gr_cell_ids == "") {
      stop(
        "'mcols(rep_tracks[[", i, "]])$cell_id' must contain exactly one ",
        "non-missing cell ID."
      )
    }

    if (!identical(gr_cell_ids, cell_id)) {
      stop(
        "Cell ID mismatch in 'rep_tracks[[", i, "]]': names(rep_tracks)[", i,
        "] is '", cell_id, "', but mcols(... )$cell_id is '", gr_cell_ids, "'."
      )
    }

    if (has_track_filter) {
      mdat <- .mcols_df(S4Vectors::mcols(gr))

      keep <- tryCatch(
        eval(filter_expr, envir = mdat, enclos = parent.frame()),
        error = function(e) {
          stop(
            "Failed to evaluate 'track_filter' for cell '", cell_id, "'. ",
            "Please make sure all variables used in 'track_filter' are columns ",
            "in mcols(rep_tracks[[i]]). Original error: ",
            conditionMessage(e),
            call. = FALSE
          )
        }
      )

      if (!is.logical(keep)) {
        stop("'track_filter' must evaluate to a logical vector.")
      }

      if (length(keep) == 1L) {
        keep <- rep(keep, length(gr))
      }

      keep[is.na(keep)] <- FALSE
      gr <- gr[keep]

      if (length(gr) == 0L) {
        next
      }
    }

    midpoint <- floor((BiocGenerics::start(gr) + BiocGenerics::end(gr)) / 2)

    centers <- GenomicRanges::GRanges(
      seqnames = GenomicRanges::seqnames(gr),
      ranges = IRanges::IRanges(start = midpoint, end = midpoint),
      strand = "*",
      seqinfo = GenomeInfoDb::seqinfo(gr)
    )

    hits <- GenomicRanges::findOverlaps(
      query = izones_ordered,
      subject = centers,
      ignore.strand = TRUE
    )

    if (length(hits) > 0L) {
      observed_fired[unique(S4Vectors::queryHits(hits)), i] <- TRUE
    }
  }

  num_fired_izones <- colSums(observed_fired)

  expected_fired <- matrix(
    FALSE,
    nrow = nrow(observed_fired),
    ncol = ncol(observed_fired),
    dimnames = dimnames(observed_fired)
  )

  for (cell_idx in seq_len(ncol(observed_fired))) {
    k <- num_fired_izones[cell_idx]

    if (!is.na(k) && k > 0L) {
      expected_fired[seq_len(k), cell_idx] <- TRUE
    }
  }

  observed_fired_expected_fired <- rowSums(
    observed_fired & expected_fired
  )

  observed_unfired_expected_unfired <- rowSums(
    !observed_fired & !expected_fired
  )

  observed_fired_expected_unfired <- rowSums(
    observed_fired & !expected_fired
  )

  observed_unfired_expected_fired <- rowSums(
    !observed_fired & expected_fired
  )

  firing_order_counts <- cbind(
    observed_fired_expected_fired = observed_fired_expected_fired,
    observed_unfired_expected_unfired = observed_unfired_expected_unfired,
    observed_fired_expected_unfired = observed_fired_expected_unfired,
    observed_unfired_expected_fired = observed_unfired_expected_fired
  )

  denominator <- rowSums(firing_order_counts)

  firing_order_frequency <- firing_order_counts / denominator * 100
  firing_order_frequency[denominator == 0, ] <- NA_real_

  rownames(firing_order_frequency) <- rownames(observed_fired)

  firing_order_frequency
}

#' Central finite differences (one-sided at the ends)
#' @keywords internal
#' @noRd
.finite_difference <- function(x, y) {
    n <- length(x)
    c(
        (y[2L] - y[1L]) / (x[2L] - x[1L]),
        (y[3L:n] - y[1L:(n - 2L)]) / (x[3L:n] - x[1L:(n - 2L)]),
        (y[n] - y[n - 1L]) / (x[n] - x[n - 1L])
    )
}

#' GAM velocity fit of one IZ
#' @keywords internal
#' @noRd
.fit_one_izone_velocity <- function(d, iz_id, config) {
    d <- d[
        is.finite(d$rep_percent) & is.finite(d$mid_track_len) &
            d$rep_percent >= 0 & d$rep_percent <= 1 & d$mid_track_len > 0,
        c("rep_percent", "mid_track_len"), drop = FALSE
    ]
    n_unique <- length(unique(d$rep_percent))
    rep_range <- if (nrow(d)) range(d$rep_percent) else c(NA_real_, NA_real_)

    fail <- function(message) list(
        summary = data.frame(
            iz_id = as.integer(iz_id), n_points = nrow(d), n_unique_rep = n_unique,
            rep_min = rep_range[1L], rep_max = rep_range[2L], rep_span = diff(rep_range),
            gam_basis_k = NA_integer_, p_value = NA_real_, edf = NA_real_,
            f_statistic = NA_real_, adjusted_r_squared = NA_real_,
            deviance_explained = NA_real_, aic_gam = NA_real_, aic_null = NA_real_,
            early_length_bp = NA_real_, late_length_bp = NA_real_,
            early_to_late_change_bp = NA_real_, median_velocity_bp_per_pct = NA_real_,
            max_abs_velocity_bp_per_pct = NA_real_, fit_status = message
        ),
        curve = NULL
    )

    if (nrow(d) < config$min_points_per_izone) return(fail("insufficient_points"))
    if (n_unique < 5L || !is.finite(diff(rep_range)) || diff(rep_range) <= 0) {
        return(fail("insufficient_replication_range"))
    }

    d$log_length <- log1p(d$mid_track_len)
    if (!is.finite(stats::var(d$log_length)) || stats::var(d$log_length) <= 0) {
        return(fail("constant_length"))
    }

    tryCatch({
        k_use <- max(4L, min(config$gam_basis_k, n_unique - 1L, nrow(d) - 2L))
        # mgcv::gam() resolves s() from the formula environment; mgcv is only
        # suggested, so bind it there instead of attaching mgcv.
        form_env <- new.env(parent = environment())
        form_env$s <- mgcv::s
        gam_formula <- stats::as.formula("log_length ~ s(rep_percent, k = k_use)",
                                         env = form_env)
        gam_fit <- suppressWarnings(mgcv::gam(
            gam_formula,
            data = d, family = stats::gaussian(), method = "GCV.Cp"
        ))
        null_fit <- stats::lm(log_length ~ 1, data = d)
        gam_summary <- summary(gam_fit)
        smooth <- gam_summary$s.table

        grid <- seq(rep_range[1L], rep_range[2L], length.out = config$n_curve_points)
        fitted_length <- pmax(expm1(stats::predict(
            gam_fit, newdata = data.frame(rep_percent = grid), type = "response"
        )), 0)
        velocity <- .finite_difference(grid * 100, fitted_length)

        list(
            summary = data.frame(
                iz_id = as.integer(iz_id), n_points = nrow(d), n_unique_rep = n_unique,
                rep_min = rep_range[1L], rep_max = rep_range[2L], rep_span = diff(rep_range),
                gam_basis_k = k_use, p_value = smooth[1L, ncol(smooth)],
                edf = smooth[1L, "edf"],
                f_statistic = if ("F" %in% colnames(smooth)) smooth[1L, "F"] else NA_real_,
                adjusted_r_squared = unname(gam_summary$r.sq),
                deviance_explained = unname(gam_summary$dev.expl),
                aic_gam = stats::AIC(gam_fit), aic_null = stats::AIC(null_fit),
                early_length_bp = fitted_length[1L], late_length_bp = fitted_length[length(grid)],
                early_to_late_change_bp = fitted_length[length(grid)] - fitted_length[1L],
                median_velocity_bp_per_pct = stats::median(velocity),
                max_abs_velocity_bp_per_pct = max(abs(velocity)), fit_status = "ok"
            ),
            curve = data.frame(
                iz_id = as.integer(iz_id), rep_percent = grid,
                rep_progress_pct = grid * 100, fitted_track_len_bp = fitted_length,
                velocity_bp_per_pct = velocity
            )
        )
    }, error = function(e) fail(conditionMessage(e)))
}

#' Track-length velocity of each initiation zone
#'
#' For each IZ, fit a GAM `log1p(mid_track_len) ~ s(rep_percent)`
#' (`mgcv::gam`, GCV, `k = max(4, min(gam_basis_k, n_unique - 1, n - 2))`)
#' to its tracks (finite values, `rep_percent` in `[0, 1]`, positive length)
#' and compare it with an intercept-only model. The fitted length curve on
#' `n_curve_points` grid points gives the velocity (bp per percent, central
#' finite differences). IZs with fewer than `min_points_per_izone` tracks or
#' fewer than 5 distinct `rep_percent` values are not fitted. GAM p-values are
#' BH-adjusted (non-finite p-values count as 1); a velocity model passes if
#' the fit is ok, `fdr < fdr_threshold` and the GAM has the lower AIC.
#'
#' @param df A `data.frame` with columns `iz_id`, `rep_percent` (fraction)
#'   and `mid_track_len` (bp), one row per track.
#' @param izone_ids IZ ids to fit (also those without rows).
#' @param min_points_per_izone Minimum number of tracks per IZ.
#' @param gam_basis_k Maximum basis dimension (e.g. `10`).
#' @param n_curve_points Number of curve points, at least 3 (e.g. `101`).
#' @param fdr_threshold FDR threshold (e.g. `0.05`).
#' @param BPPARAM A `BiocParallel::BiocParallelParam` used over IZs.
#'
#' @return A list with `summary` (one row per IZ: fit statistics, `fdr`,
#'   `delta_aic_null_minus_gam`, `aic_pass`, `velocity_model_pass`,
#'   `velocity_direction`, `fit_status`), `curves` (fitted length and velocity
#'   per grid point) and `config` (the parameters used).
#'
#' @export
fit_all_izone_velocities <- function(df,
                                     izone_ids,
                                     min_points_per_izone,
                                     gam_basis_k = 10L,
                                     n_curve_points = 101L,
                                     fdr_threshold = 0.05,
                                     BPPARAM = BiocParallel::SerialParam()) {
    rlang::check_installed("mgcv", reason = "to fit IZ velocity GAMs.")
    if (!inherits(BPPARAM, "BiocParallelParam")) {
        stop("'BPPARAM' must be a BiocParallelParam object.")
    }
    config <- list(
        min_points_per_izone = min_points_per_izone,
        gam_basis_k = gam_basis_k,
        n_curve_points = n_curve_points,
        fdr_threshold = fdr_threshold
    )
    stopifnot(
        all(c("iz_id", "rep_percent", "mid_track_len") %in% names(df)),
        config$n_curve_points >= 3L
    )
    izone_ids <- as.integer(izone_ids)
    groups <- split(
        as.data.frame(df),
        factor(as.integer(df$iz_id), levels = izone_ids),
        drop = FALSE
    )
    fit_index <- function(i) .fit_one_izone_velocity(groups[[i]], izone_ids[i], config)
    fits <- BiocParallel::bplapply(seq_along(groups), fit_index, BPPARAM = BPPARAM)

    model_summary <- do.call(rbind, lapply(fits, `[[`, "summary"))
    p_for_adjustment <- model_summary$p_value
    p_for_adjustment[!is.finite(p_for_adjustment)] <- 1
    model_summary$fdr <- stats::p.adjust(p_for_adjustment, method = "BH")
    model_summary$delta_aic_null_minus_gam <- model_summary$aic_null - model_summary$aic_gam
    model_summary$aic_pass <- model_summary$delta_aic_null_minus_gam > 0
    model_summary$velocity_model_pass <-
        model_summary$fit_status == "ok" &
        model_summary$fdr < config$fdr_threshold & model_summary$aic_pass
    model_summary$velocity_direction <- ifelse(
        !model_summary$velocity_model_pass, "not_supported",
        ifelse(model_summary$early_to_late_change_bp > 0, "increasing", "decreasing")
    )

    curves <- do.call(rbind, Filter(Negate(is.null), lapply(fits, `[[`, "curve")))
    if (is.null(curves)) curves <- data.frame()
    if (nrow(curves)) {
        curves <- merge(
            curves,
            model_summary[c("iz_id", "fdr", "delta_aic_null_minus_gam", "velocity_model_pass")],
            by = "iz_id", all.x = TRUE, sort = FALSE
        )
        curves <- curves[order(curves$iz_id, curves$rep_percent), ]
        rownames(curves) <- NULL
    }

    list(summary = model_summary[order(model_summary$iz_id), ], curves = curves, config = config)
}
