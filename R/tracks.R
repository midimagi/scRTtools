# Replication tracks and their lengths / gaps.
# define_replication_tracks() is migrated from
# rcode/processing/define_replication_tracks.R, calc_rep_track_lengths() and
# calc_rep_track_gaps() and count_track_centers() from rcode/analysis/.
# merge_track_centers() (no notebook caller) is not included.

#' Identify replication tracks from bin-level replication states
#'
#' Identify contiguous segments with replication state `cn` (typically `4`)
#' in each cell of a bins x cells state matrix. R translation of the MATLAB
#' functions `define_replication_tracks()` and `list_cn_segments()`.
#'
#' Within each chromosome and cell, bins that are `NA` in all cells or not
#' mappable are treated as `NA`. A new segment starts at the first bin, when
#' the state changes, at the first non-`NA` bin after an `NA` bin with no
#' previous state in the segment, and after more than `max_na_run`
#' consecutive `NA` bins. Trailing `NA` bins belong to the preceding segment.
#' Segments whose `NA` fraction exceeds `max_na_fraction` are dropped, and
#' only segments with state `cn` are returned.
#'
#' The `flag` column (MATLAB fourth column) is `1` if a segment neighbours a
#' segment without a state (all `NA`) or is adjacent to another segment with
#' state `cn`.
#'
#' @param genobins_matrix Numeric matrix of replication states (`2`, `4` or
#'   `NA`), bins x cells, with row names equal to `names(genome_granges)` and
#'   unique column names (e.g. `compute_rep_percent()$genobins_rep_states_filtered`).
#' @param genome_granges Named `GRanges`, one range per row, with an
#'   `is_mappable` metadata column (`NA` counts as not mappable).
#' @param rep_percent Named numeric vector of replicated fractions; names
#'   must include all column names of `genobins_matrix` (no `NA`).
#' @param cn State to extract.
#' @param max_na_run Non-negative integer; a new segment starts after more
#'   than this many consecutive `NA` bins.
#' @param max_na_fraction Number in `[0, 1]`; maximum `NA` fraction of a kept
#'   segment.
#' @param BPPARAM A `BiocParallel::BiocParallelParam`; chromosomes are
#'   processed in parallel.
#' @param verbose Logical; print progress messages. Default `FALSE`.
#'
#' @return A `GRangesList` named by cell; each element holds the cell's
#'   segments with metadata columns `cell_id`, `rep_percent`, `cn`,
#'   `num_bins`, `flag` and `na_fraction`.
#'
#' @export
define_replication_tracks <- function(
    genobins_matrix,
    genome_granges,
    rep_percent,
    cn = 4,
    max_na_run = 5L,
    max_na_fraction = 0.5,
    BPPARAM = BiocParallel::SerialParam(),
    verbose = FALSE
) {

  stopifnot(is.matrix(genobins_matrix))
  stopifnot(length(genome_granges) == nrow(genobins_matrix))
  stopifnot(!is.null(rownames(genobins_matrix)))
  stopifnot(!is.null(names(genome_granges)))
  stopifnot(all(rownames(genobins_matrix) == names(genome_granges)))

  cell_names <- colnames(genobins_matrix)

  stopifnot(
    !is.null(cell_names),
    !anyNA(cell_names),
    !anyDuplicated(cell_names)
  )

  stopifnot(
    is.numeric(rep_percent),
    !is.null(names(rep_percent)),
    all(cell_names %in% names(rep_percent)),
    !anyNA(rep_percent[cell_names])
  )

  if (!"is_mappable" %in% colnames(S4Vectors::mcols(genome_granges))) {
    stop("mcols(genome_granges) must contain an 'is_mappable' column.")
  }

  stopifnot(
    length(cn) == 1L,
    is.numeric(cn),
    !is.na(cn)
  )

  stopifnot(
    length(max_na_run) == 1L,
    is.numeric(max_na_run),
    !is.na(max_na_run),
    max_na_run >= 0,
    max_na_run == as.integer(max_na_run)
  )
  max_na_run <- as.integer(max_na_run)
  if (!inherits(BPPARAM, "BiocParallelParam")) {
    stop("'BPPARAM' must be a BiocParallelParam object.")
  }

  stopifnot(
    length(max_na_fraction) == 1L,
    is.numeric(max_na_fraction),
    !is.na(max_na_fraction),
    max_na_fraction >= 0,
    max_na_fraction <= 1
  )

  num_cells <- ncol(genobins_matrix)
  seqs <- as.character(GenomeInfoDb::seqnames(genome_granges))
  chromosomes <- unique(seqs)

  if (verbose) {
    message(
      "[define_replication_tracks] Input: ",
      nrow(genobins_matrix), " bins x ",
      num_cells, " cells across ",
      length(chromosomes), " chromosomes."
    )
    message(
      "[define_replication_tracks] Extracting CN = ",
      cn, " tracks using ",
      class(BPPARAM)[1], "."
    )
  }

  chr_results <- BiocParallel::bplapply(
    chromosomes,
    function(chr) {
      chr_idx <- which(seqs == chr)

      mat_chr <- genobins_matrix[chr_idx, , drop = FALSE]
      gr_chr <- genome_granges[chr_idx]

      .define_cn_segments(
        replication_state_chr = mat_chr,
        genome_granges_chr = gr_chr,
        rep_percent = rep_percent,
        cn = cn,
        max_na_run = max_na_run,
        max_na_fraction = max_na_fraction,
        cell_names = cell_names
      )
    },
    BPPARAM = BPPARAM
  )

  names(chr_results) <- chromosomes

  tracks_by_cell <- vector("list", num_cells)
  names(tracks_by_cell) <- cell_names

  for (i in seq_len(num_cells)) {
    tracks_by_cell[[i]] <- GenomicRanges::GRanges()
  }

  for (chr in chromosomes) {
    chr_tracks <- chr_results[[chr]]

    for (i in seq_len(num_cells)) {
      if (length(chr_tracks[[i]]) > 0L) {
        tracks_by_cell[[i]] <- c(tracks_by_cell[[i]], chr_tracks[[i]])
      }
    }
  }

  tracks_by_cell <- GenomicRanges::GRangesList(tracks_by_cell)

  if (verbose) {
    n_tracks <- sum(lengths(tracks_by_cell))
    message(
      "[define_replication_tracks] Done. Defined ",
      n_tracks, " replicated segments."
    )
  }

  tracks_by_cell
}


#' Segments of one chromosome (MATLAB list_cn_segments)
#' @keywords internal
#' @noRd
.define_cn_segments <- function(
    replication_state_chr,
    genome_granges_chr,
    rep_percent,
    cn,
    max_na_run,
    max_na_fraction,
    cell_names
) {
  replication_state_chr <- as.matrix(replication_state_chr)

  num_cells <- ncol(replication_state_chr)
  num_bins <- nrow(replication_state_chr)

  stopifnot(
    num_bins > 0L,
    num_cells > 0L,
    length(genome_granges_chr) == num_bins,
    length(cell_names) == num_cells
  )

  is_mappable <- S4Vectors::mcols(genome_granges_chr)$is_mappable

  if (is.null(is_mappable)) {
    is_mappable <- rep(TRUE, num_bins)
  }

  stopifnot(length(is_mappable) == num_bins)

  is_mappable <- as.logical(is_mappable)
  is_mappable[is.na(is_mappable)] <- FALSE

  is_gap <- rowSums(is.na(replication_state_chr)) == num_cells
  is_gap <- is_gap | !is_mappable

  replicated_segments <- vector("list", num_cells)
  names(replicated_segments) <- cell_names

  for (barcode in seq_len(num_cells)) {
    x <- replication_state_chr[, barcode]
    x[is_gap] <- NA

    segment_start_flag <- rep(FALSE, num_bins)
    segment_start_flag[1] <- TRUE

    last_non_na_cn <- NA_real_
    na_counter <- 0L
    need_new_segment_after_long_na <- FALSE

    for (row in seq_len(num_bins)) {
      current_cn <- x[row]

      if (is.na(current_cn)) {
        na_counter <- na_counter + 1L

        if (na_counter > max_na_run) {
          need_new_segment_after_long_na <- TRUE
          last_non_na_cn <- NA_real_
        }

        next
      }

      if (need_new_segment_after_long_na) {
        segment_start_flag[row] <- TRUE
        last_non_na_cn <- current_cn
        need_new_segment_after_long_na <- FALSE
        na_counter <- 0L
        next
      }

      na_counter <- 0L

      if (is.na(last_non_na_cn)) {
        last_non_na_cn <- current_cn

        if (row > 1L && is.na(x[row - 1L])) {
          segment_start_flag[row] <- TRUE
        }

      } else if (current_cn != last_non_na_cn) {
        segment_start_flag[row] <- TRUE
        last_non_na_cn <- current_cn
      }
    }

    segment_starts <- which(segment_start_flag)
    segment_ends <- c(segment_starts[-1] - 1L, num_bins)

    segments <- data.frame(
      bin_start = segment_starts,
      bin_end = segment_ends,
      cn = NA_real_,
      flag = 0L,
      stringsAsFactors = FALSE
    )

    for (i in seq_len(nrow(segments))) {
      idx <- segments$bin_start[i]:segments$bin_end[i]
      segment_cn <- unique(x[idx][!is.na(x[idx])])

      if (length(segment_cn) == 1L) {
        segments$cn[i] <- segment_cn
      }
    }

    for (i in seq_len(nrow(segments))) {
      neighbors <- c(i - 1L, i + 1L)
      neighbors <- neighbors[
        neighbors >= 1L & neighbors <= nrow(segments)
      ]

      segments$flag[i] <- as.integer(
        !is.na(segments$cn[i]) &&
          length(neighbors) > 0L &&
          any(is.na(segments$cn[neighbors]))
      )
    }

    if (nrow(segments) >= 2L) {
      for (i in seq_len(nrow(segments) - 1L)) {
        if (
          !is.na(segments$cn[i]) &&
            !is.na(segments$cn[i + 1L]) &&
            segments$cn[i] == cn &&
            segments$cn[i + 1L] == cn
        ) {
          segments$flag[i] <- 1L
          segments$flag[i + 1L] <- 1L
        }
      }
    }

    segments$num_bins <- segments$bin_end - segments$bin_start + 1L
    segments$na_fraction <- NA_real_

    for (i in seq_len(nrow(segments))) {
      idx <- segments$bin_start[i]:segments$bin_end[i]
      segments$na_fraction[i] <- mean(is.na(x[idx]))
    }

    has_too_many_NAs <- segments$na_fraction > max_na_fraction
    segments <- segments[!has_too_many_NAs, , drop = FALSE]

    keep <- !is.na(segments$cn) & segments$cn == cn
    segments <- segments[keep, , drop = FALSE]

    if (nrow(segments) == 0L) {
      replicated_segments[[barcode]] <- GenomicRanges::GRanges()
      next
    }

    gr <- GenomicRanges::GRanges(
      seqnames = GenomeInfoDb::seqnames(
        genome_granges_chr[segments$bin_start]
      ),
      ranges = IRanges::IRanges(
        start = BiocGenerics::start(
          genome_granges_chr[segments$bin_start]
        ),
        end = BiocGenerics::end(
          genome_granges_chr[segments$bin_end]
        )
      ),
      strand = "*"
    )

    cell_id <- cell_names[barcode]

    S4Vectors::mcols(gr)$cell_id <- cell_id
    S4Vectors::mcols(gr)$rep_percent <- unname(rep_percent[cell_id])
    S4Vectors::mcols(gr)$cn <- segments$cn
    S4Vectors::mcols(gr)$num_bins <- segments$num_bins
    S4Vectors::mcols(gr)$flag <- segments$flag
    S4Vectors::mcols(gr)$na_fraction <- segments$na_fraction

    replicated_segments[[barcode]] <- gr
  }

  replicated_segments
}

#' Table of replication-track lengths
#'
#' Flatten a track `GRangesList` into one row per track.
#'
#' @param rep_tracks Named `GRangesList` (one element per cell), e.g. from
#'   [define_replication_tracks()].
#'
#' @return A `data.frame` with columns `chr`, `start`, `end`, `width_bp`
#'   followed by the metadata columns of the tracks. If there are no tracks,
#'   an empty `data.frame` with columns `cell`, `chr`, `start`, `end` and
#'   `width_bp` (note: `cell`, not `cell_id`).
#'
#' @export
calc_rep_track_lengths <- function(rep_tracks) {
  if (!inherits(rep_tracks, "GRangesList")) {
    stop("'rep_tracks' must be a GRangesList or CompressedGRangesList object.")
  }

  if (is.null(names(rep_tracks)) || anyNA(names(rep_tracks))) {
    stop("names(rep_tracks) must contain cell names.")
  }

  length_list <- vector("list", length(rep_tracks))
  names(length_list) <- names(rep_tracks)

  for (i in seq_along(rep_tracks)) {
    gr <- rep_tracks[[i]]

    if (length(gr) == 0L) {
      next
    }

    mdat <- .mcols_df(S4Vectors::mcols(gr))

    df <- data.frame(
      chr = as.character(GenomicRanges::seqnames(gr)),
      start = BiocGenerics::start(gr),
      end = BiocGenerics::end(gr),
      width_bp = BiocGenerics::width(gr),
      stringsAsFactors = FALSE
    )

    if (ncol(mdat) > 0L) {
      df <- cbind(df, mdat)
    }

    length_list[[i]] <- df
  }

  length_df <- do.call(rbind, length_list)

  if (is.null(length_df)) {
    length_df <- data.frame(
      cell = character(),
      chr = character(),
      start = integer(),
      end = integer(),
      width_bp = integer(),
      stringsAsFactors = FALSE
    )
  }

  rownames(length_df) <- NULL

  length_df
}

#' Gaps between adjacent replication tracks
#'
#' For each cell, sort the tracks and compute the distance between
#' consecutive tracks on the same chromosome (`start[i + 1] - end[i] - 1`).
#'
#' @param rep_tracks Named `GRangesList` (one element per cell), e.g. from
#'   [define_replication_tracks()], with metadata columns `cell_id` and
#'   `rep_percent`.
#'
#' @return A `data.frame` with columns `cell_id`, `rep_percent`, `chr` and
#'   `gap_bp`, one row per adjacent pair.
#'
#' @export
calc_rep_track_gaps <- function(rep_tracks) {
  if (!inherits(rep_tracks, "GRangesList")) {
    stop("'rep_tracks' must be a GRangesList or CompressedGRangesList object.")
  }

  if (is.null(names(rep_tracks)) || anyNA(names(rep_tracks))) {
    stop("names(rep_tracks) must contain cell names.")
  }

  gap_list <- vector("list", length(rep_tracks))
  names(gap_list) <- names(rep_tracks)

  for (i in seq_along(rep_tracks)) {
    gr <- rep_tracks[[i]]
    cell_id <- S4Vectors::mcols(gr)$cell_id[1]
    rep_percent <- S4Vectors::mcols(gr)$rep_percent[1]

    if (length(gr) < 2L) {
      next
    }

    ## Sort tracks by chromosome and genomic coordinate
    gr <- BiocGenerics::sort(gr)

    seqs <- as.character(GenomicRanges::seqnames(gr))
    starts <- BiocGenerics::start(gr)
    ends <- BiocGenerics::end(gr)

    ## Only calculate gaps between adjacent tracks on the same chromosome
    same_chr <- seqs[-1] == seqs[-length(seqs)]

    if (!any(same_chr)) {
      next
    }

    gap_bp <- starts[-1] - ends[-length(ends)] - 1L

    ## Important: chr and gap_bp must be filtered by the same same_chr vector
    gap_list[[i]] <- data.frame(
      cell_id = cell_id,
      rep_percent = rep_percent,
      chr = seqs[-1][same_chr],
      gap_bp = as.numeric(gap_bp[same_chr]),
      stringsAsFactors = FALSE
    )
  }

  gap_df <- do.call(rbind, gap_list)

  if (is.null(gap_df)) {
    gap_df <- data.frame(
      cell_id = character(),
      rep_percent = numeric(),
      chr = character(),
      gap_bp = numeric(),
      stringsAsFactors = FALSE
    )
  }

  rownames(gap_df) <- NULL

  gap_df
}

#' Count replication-track centres in genomic bins
#'
#' Optionally filter each cell's tracks with an expression evaluated on their
#' metadata columns, then count, for every range of `granges`, the tracks
#' whose centre `floor((start + end) / 2)` falls inside it.
#'
#' @param tracks Named `GRangesList` of tracks (one element per cell), e.g.
#'   from [define_replication_tracks()].
#' @param granges `GRanges` of the bins to count in.
#' @param track_filter Optional unquoted expression evaluated against
#'   `mcols(tracks[[cell]])` (non-standard evaluation), e.g.
#'   `(flag == 0) & (num_bins > 50)` or `!is.na(iz_id)`. It must give a
#'   logical scalar or one value per track; `NA` counts as `FALSE`.
#'
#' @return `granges`, sorted, with an integer `score` column (number of track
#'   centres; zero-count ranges are kept).
#'
#' @export
count_track_centers <- function(
    tracks,
    granges,
    track_filter = NULL
) {
  if (!inherits(tracks, "GRangesList")) {
    stop("`tracks` must be a GRangesList or CompressedGRangesList.")
  }

  if (!inherits(granges, "GRanges")) {
    stop("`granges` must be a GenomicRanges::GRanges object.")
  }

  out <- granges
  S4Vectors::mcols(out)$score <- integer(length(out))

  filter_expr <- substitute(track_filter)
  has_track_filter <- !identical(filter_expr, quote(NULL))

  gr_list <- vector("list", length(tracks))

  track_names <- names(tracks)

  for (i in seq_along(tracks)) {
    gr <- tracks[[i]]

    if (length(gr) == 0L) {
      next
    }

    if (has_track_filter) {
      mdat <- .mcols_df(S4Vectors::mcols(gr))

      keep <- tryCatch(
        eval(filter_expr, envir = mdat, enclos = parent.frame()),
        error = function(e) {
          cell_msg <- if (!is.null(track_names) && !is.na(track_names[i])) {
            paste0(" for cell '", track_names[i], "'")
          } else {
            paste0(" for element ", i)
          }

          stop(
            "Failed to evaluate 'track_filter'", cell_msg, ". ",
            "Please make sure all variables used in 'track_filter' are columns ",
            "in mcols(tracks[[i]]). Original error: ",
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

      if (length(keep) != length(gr)) {
        stop(
          "'track_filter' must evaluate to length 1 or length of each cell's GRanges."
        )
      }

      keep[is.na(keep)] <- FALSE
      gr <- gr[keep]

      if (length(gr) == 0L) {
        next
      }
    }

    gr_list[[i]] <- gr
  }

  gr_list <- gr_list[!vapply(gr_list, is.null, logical(1))]

  if (length(gr_list) == 0L) {
    return(BiocGenerics::sort(out))
  }

  gr <- do.call(c, gr_list)

  if (length(gr) == 0L) {
    return(BiocGenerics::sort(out))
  }

  midpoint <- floor((BiocGenerics::start(gr) + BiocGenerics::end(gr)) / 2)

  centers <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(gr),
    ranges = IRanges::IRanges(
      start = midpoint,
      end = midpoint
    ),
    strand = "*"
  )

  GenomeInfoDb::seqinfo(centers) <- GenomeInfoDb::seqinfo(gr)

  hits <- GenomicRanges::findOverlaps(
    query = out,
    subject = centers,
    ignore.strand = TRUE
  )

  if (length(hits) > 0L) {
    score <- tabulate(
      S4Vectors::queryHits(hits),
      nbins = length(out)
    )

    S4Vectors::mcols(out)$score <- as.integer(score)
  }

  BiocGenerics::sort(out)
}
