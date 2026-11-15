# Genome binning: fixed-width windows, GC fraction and mappability.
# Migrated from rcode/processing/make_fixed_bins.R (notebook 1.1).

#' Define fixed-width genomic windows
#'
#' Create fixed-width genomic windows for selected chromosomes using chromosome
#' information from a `Seqinfo` object. Chromosome names are normalised to the
#' `"chr"`-prefixed style, both in `seqinfo` and in `use_chr`.
#'
#' @param seqinfo A `GenomeInfoDb::Seqinfo` object containing chromosome names
#'   and lengths.
#' @param bin_width A single positive integer specifying the window size (bp).
#'   We use `20000L`.
#' @param use_chr Chromosomes to use, with or without the `"chr"` prefix.
#'   Defaults to `c(1:19, "X")` (mm10 autosomes + X); use `1:19` for autosomes only.
#'
#' @return A `GRanges` of contiguous, non-overlapping windows
#'   `[1 + (k-1)*bin_width, min(k*bin_width, seqlength)]` per chromosome, with
#'   no metadata columns and with `seqinfo` set to the (chr-normalised) input
#'   `seqinfo`. The last window of each chromosome may be shorter.
#'
#' @examples
#' si <- GenomeInfoDb::Seqinfo(c("chr1", "chr2"), c(105000L, 60000L))
#' make_fixed_bins(si, bin_width = 20000L, use_chr = 1:2)
#'
#' @export
make_fixed_bins <- function(seqinfo,
                            bin_width,
                            use_chr = c(1:19, "X")) {
  if (!methods::is(seqinfo, "Seqinfo")) {
    stop("'seqinfo' must be a Seqinfo object.")
  }

  if (!is.numeric(bin_width) || length(bin_width) != 1L ||
      is.na(bin_width) || bin_width <= 0) {
    stop("'bin_width' must be a single positive integer.")
  }
  bin_width <- as.integer(bin_width)

  ## Normalize seqlevels in Seqinfo to chr-style
  old_seqlevels <- GenomeInfoDb::seqlevels(seqinfo)
  new_seqlevels <- .normalize_chr(old_seqlevels)

  if (anyDuplicated(new_seqlevels)) {
    stop("Duplicated chromosome names found in 'seqinfo' after normalization.")
  }

  GenomeInfoDb::seqlevels(seqinfo) <- new_seqlevels

  chr_names <- .normalize_chr(use_chr)

  if (!all(chr_names %in% GenomeInfoDb::seqlevels(seqinfo))) {
    missing_chr <- chr_names[!chr_names %in% GenomeInfoDb::seqlevels(seqinfo)]
    stop(
      "Not all requested chromosomes are present in 'seqinfo': ",
      paste(missing_chr, collapse = ", ")
    )
  }

  seqlens <- GenomeInfoDb::seqlengths(seqinfo)[chr_names]

  if (anyNA(seqlens) || any(seqlens <= 0)) {
    stop("Requested chromosomes must have positive and non-missing seqlengths.")
  }

  seqlens <- as.integer(seqlens)
  names(seqlens) <- chr_names

  n_win <- as.integer(ceiling(seqlens / bin_width))

  seqnames_vec <- rep.int(chr_names, times = n_win)
  window_id <- sequence(n_win)

  starts <- 1L + (window_id - 1L) * bin_width
  chr_end <- rep.int(seqlens, times = n_win)
  ends <- pmin.int(starts + bin_width - 1L, chr_end)

  gr <- GenomicRanges::GRanges(
    seqnames = seqnames_vec,
    ranges = IRanges::IRanges(start = starts, end = ends)
  )

  ## Assign the Seqinfo object to returned GRanges
  GenomeInfoDb::seqinfo(gr) <- seqinfo

  gr
}

#' Compute GC fraction and mappability for fixed genome windows
#'
#' Compute the per-window GC fraction from a reference FASTA and the per-window
#' mappable fraction from a mappability track (e.g. a GenMap bedGraph).
#' Windows whose coordinates fall outside the reference sequence are retained
#' and get `NA` GC fraction.
#'
#' A base is counted as mappable when it lies in an interval with
#' `score >= th`; the mappable fraction of a window is the number of such bases
#' divided by the window width, and `is_mappable` is `mappability > window_prop`
#' (strictly greater).
#'
#' @param genome_bins A `GRanges` of fixed windows, e.g. from
#'   [make_fixed_bins()].
#' @param genome_mappability A bedGraph file path, or a `GRanges` with a
#'   numeric `score` metadata column. Importing a file requires
#'   \pkg{rtracklayer}.
#' @param th Numeric scalar in `(0, 1]`. Mappability score threshold.
#' @param window_prop Numeric scalar in `[0, 1]`. A window is mappable when its
#'   mappable fraction is strictly greater than this value.
#' @param ref_path Character scalar. Path to the reference genome FASTA
#'   (requires \pkg{Biostrings}). Sequence names may be with or without the
#'   `"chr"` prefix.
#' @param BPPARAM A `BiocParallel::BiocParallelParam`; GC content is computed
#'   per chromosome in parallel. We used
#'   `BiocParallel::MulticoreParam(workers = 8L)`.
#' @param verbose Logical; print progress messages. Default `FALSE`.
#'
#' @return `genome_bins` with added metadata columns:
#'   \describe{
#'     \item{`bin_id`}{character, `"<chr>-<start>-<end>"`.}
#'     \item{`gc_frequency`}{numeric GC fraction in `[0, 1]` (`NA` for windows
#'       outside the reference).}
#'     \item{`mappability`}{numeric mappable fraction in `[0, 1]`.}
#'     \item{`is_mappable`}{logical, `mappability > window_prop`.}
#'   }
#'   If no mappability interval survives chromosome filtering or the `th`
#'   threshold, the column `frac_mappable = 0` is added instead of
#'   `mappability` and all windows
#'   are `is_mappable = FALSE`.
#'
#' @examples
#' if (requireNamespace("Biostrings", quietly = TRUE)) {
#'   fa <- tempfile(fileext = ".fa")
#'   seqs <- Biostrings::DNAStringSet(c(
#'     chr1 = paste(rep("ACGTGGCCAT", 100), collapse = "")
#'   ))
#'   Biostrings::writeXStringSet(seqs, fa)
#'   si <- GenomeInfoDb::Seqinfo("chr1", 1000L)
#'   bins <- make_fixed_bins(si, 200L, use_chr = 1)
#'   map <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 501), c(500, 1000)),
#'                                 score = c(1, 0.5))
#'   GenomeInfoDb::seqinfo(map) <- si
#'   compute_bin_mappability(bins, map, ref_path = fa)
#' }
#'
#' @export
compute_bin_mappability <- function(genome_bins,
                                    genome_mappability,
                                    th = 0.8,
                                    window_prop = 0.7,
                                    ref_path,
                                    BPPARAM = BiocParallel::SerialParam(),
                                    verbose = FALSE) {
  rlang::check_installed("Biostrings", reason = "to compute GC fraction.")

  vmsg <- function(...) .vmsg(verbose, ...)

  vmsg("Validating input parameters...")

  .validate_bins(genome_bins, "genome_bins")

  if (!is.numeric(th) || length(th) != 1L || is.na(th) || th <= 0 || th > 1) {
    stop("'th' must be a single numeric value in (0, 1].")
  }
  if (!is.numeric(window_prop) || length(window_prop) != 1L ||
      is.na(window_prop) || window_prop < 0 || window_prop > 1) {
    stop("'window_prop' must be a single numeric value in [0, 1].")
  }

  if (missing(ref_path) || !is.character(ref_path) || length(ref_path) != 1L ||
      !file.exists(ref_path)) {
    stop("'ref_path' must be a valid FASTA file path.")
  }

  if (!inherits(BPPARAM, "BiocParallelParam")) {
    stop("'BPPARAM' must be a BiocParallelParam object.")
  }

  vmsg("Input validation complete.")
  vmsg("Number of genome bins: ", length(genome_bins))
  vmsg("Mappability threshold: ", th)
  vmsg("Window proportion threshold: ", window_prop)

  import_mappability <- function(x) {
    vmsg("Importing mappability data...")
    if (inherits(x, "GRanges")) {
      gr <- x
      vmsg("Using provided GRanges object for mappability.")
    } else if (is.character(x) && length(x) == 1L && file.exists(x)) {
      rlang::check_installed("rtracklayer",
                             reason = "to import a mappability bedGraph.")
      gr <- rtracklayer::import(x)
      vmsg("Imported mappability from file: ", x)
    } else {
      stop("'genome_mappability' must be a bedGraph file path or a GRanges object.")
    }

    if (!"score" %in% names(S4Vectors::mcols(gr))) {
      stop("'genome_mappability' must contain a metadata column named 'score'.")
    }

    score <- as.numeric(S4Vectors::mcols(gr)$score)
    if (all(is.na(score))) {
      stop("All mappability scores are NA.")
    }
    S4Vectors::mcols(gr)$score <- score
    vmsg("Mappability data imported successfully. Total intervals: ", length(gr))
    gr
  }

  out <- genome_bins
  S4Vectors::mcols(out)$bin_id <- paste(GenomeInfoDb::seqnames(out),
                                        BiocGenerics::start(out),
                                        BiocGenerics::end(out), sep = "-")

  ## -----------------------------
  ## GC fraction per window
  ## -----------------------------
  vmsg("Computing GC fraction per window...")

  vmsg("Loading reference genome from: ", ref_path)
  reference <- Biostrings::readDNAStringSet(ref_path)
  if (length(reference) == 0L) {
    stop("No sequences found in 'ref_path'.")
  }
  if (is.null(names(reference)) || anyNA(names(reference)) ||
      any(names(reference) == "")) {
    stop("Reference FASTA must contain non-empty sequence names.")
  }

  vmsg("Reference genome loaded. Number of sequences: ", length(reference))

  ref_norm <- .normalize_chr(names(reference))
  if (anyDuplicated(ref_norm)) {
    stop("Duplicated chromosome names found in reference after normalization.")
  }
  ref_idx_map <- stats::setNames(seq_along(reference), ref_norm)

  out_chr <- as.character(GenomeInfoDb::seqnames(out))
  out_chr_norm <- .normalize_chr(out_chr)
  out_start <- GenomicRanges::start(out)
  out_end <- GenomicRanges::end(out)

  frac_gc <- rep(NA_real_, length(out))

  has_ref <- out_chr_norm %in% names(ref_idx_map)

  ref_len <- rep(NA_integer_, length(out))
  if (any(has_ref)) {
    ref_len[has_ref] <- Biostrings::width(reference[ref_idx_map[out_chr_norm[has_ref]]])
  }

  valid_gc <- has_ref &
    !is.na(out_start) &
    !is.na(out_end) &
    out_start >= 1L &
    out_end >= out_start &
    !is.na(ref_len) &
    out_end <= ref_len

  valid_idx <- which(valid_gc)

  vmsg("Total bins: ", length(out))
  vmsg("Valid bins for GC calculation: ", length(valid_idx))
  vmsg("Invalid bins (will have NA GC): ", length(out) - length(valid_idx))

  if (length(valid_idx) > 0L) {
    valid_split <- split(valid_idx, out_chr_norm[valid_idx])
    vmsg("Computing GC content for ", length(valid_split), " chromosomes in parallel...")

    gc_res <- BiocParallel::bplapply(
      X = names(valid_split),
      FUN = function(chr_norm) {
        idx <- valid_split[[chr_norm]]
        ref_seq <- reference[[ref_idx_map[[chr_norm]]]]

        views <- IRanges::Views(
          subject = ref_seq,
          start = out_start[idx],
          end = out_end[idx]
        )

        gc_mat <- Biostrings::letterFrequency(
          x = views,
          letters = c("G", "C"),
          as.prob = TRUE
        )

        list(
          idx = idx,
          frac_gc = rowSums(gc_mat)
        )
      },
      BPPARAM = BPPARAM
    )

    for (res in gc_res) {
      frac_gc[res$idx] <- res$frac_gc
    }
  }

  S4Vectors::mcols(out)$gc_frequency <- frac_gc
  vmsg("GC fraction calculation complete.")

  ## -----------------------------
  ## Mappability
  ## -----------------------------
  vmsg("Computing mappability...")

  map_gr <- import_mappability(genome_mappability)

  win_seqlevels <- GenomeInfoDb::seqlevels(out)
  win_norm <- .normalize_chr(win_seqlevels)

  if (anyDuplicated(win_norm)) {
    stop("Duplicated chromosome names found in 'genome_bins' after normalization.")
  }

  norm_to_win <- stats::setNames(win_seqlevels, win_norm)

  map_chr_norm <- .normalize_chr(as.character(GenomicRanges::seqnames(map_gr)))
  keep_chr <- map_chr_norm %in% names(norm_to_win)

  map_gr <- map_gr[keep_chr]
  map_chr_norm <- map_chr_norm[keep_chr]

  vmsg("Mappability intervals after chromosome filtering: ", length(map_gr))

  if (length(map_gr) == 0L) {
    vmsg("No mappability intervals remain after chromosome filtering. Returning all bins as unmappable.")
    # NOTE (kept from original): column name differs from the main path
    # ('frac_mappable' vs 'mappability'). Reported as a suspected bug.
    S4Vectors::mcols(out)$frac_mappable <- 0
    S4Vectors::mcols(out)$is_mappable <- FALSE
    return(out)
  }

  mapped_seqnames <- unname(norm_to_win[map_chr_norm])

  map_gr <- GenomicRanges::GRanges(
    seqnames = mapped_seqnames,
    ranges = GenomicRanges::ranges(map_gr),
    strand = "*",
    score = as.numeric(S4Vectors::mcols(map_gr)$score),
    seqinfo = GenomicRanges::seqinfo(out)
  )

  seqlens <- GenomeInfoDb::seqlengths(out)
  idx_has_len <- !is.na(seqlens[as.character(GenomicRanges::seqnames(map_gr))])

  if (any(idx_has_len)) {
    chr_len <- seqlens[as.character(GenomicRanges::seqnames(map_gr[idx_has_len]))]
    if (any(GenomicRanges::start(map_gr[idx_has_len]) < 1L)) {
      stop("Imported mappability intervals contain coordinates < 1.")
    }
    if (any(GenomicRanges::end(map_gr[idx_has_len]) > chr_len)) {
      stop("Some imported mappability intervals exceed chromosome lengths in seqinfo(genome_bins).")
    }
  }

  keep_gr <- map_gr[
    !is.na(S4Vectors::mcols(map_gr)$score) &
      S4Vectors::mcols(map_gr)$score >= th
  ]

  vmsg("Mappability intervals meeting threshold (score >= ", th, "): ", length(keep_gr))

  if (length(keep_gr) == 0L) {
    vmsg("No intervals meet mappability threshold. Returning all bins as unmappable.")
    S4Vectors::mcols(out)$frac_mappable <- 0
    S4Vectors::mcols(out)$is_mappable <- FALSE
    return(out)
  }

  keep_gr <- GenomicRanges::reduce(keep_gr, ignore.strand = TRUE)
  vmsg("Reduced mappability regions: ", length(keep_gr))

  vmsg("Calculating overlaps between bins and mappability regions...")
  hits <- GenomicRanges::findOverlaps(out, keep_gr, ignore.strand = TRUE)

  retained_bp <- numeric(length(out))
  window_width <- IRanges::width(out)

  if (length(hits) > 0L) {
    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)

    ov <- IRanges::pintersect(
      GenomicRanges::ranges(out)[qh],
      GenomicRanges::ranges(keep_gr)[sh]
    )
    ovw <- IRanges::width(ov)

    retained_by_window <- rowsum(ovw, group = qh, reorder = FALSE)
    retained_bp[as.integer(rownames(retained_by_window))] <- retained_by_window[, 1]
  }

  frac_mappable <- retained_bp / window_width

  S4Vectors::mcols(out)$mappability <- frac_mappable
  S4Vectors::mcols(out)$is_mappable <- frac_mappable > window_prop

  n_mappable <- sum(frac_mappable > window_prop, na.rm = TRUE)
  n_unmappable <- sum(frac_mappable <= window_prop, na.rm = TRUE)
  n_na <- sum(is.na(frac_mappable))

  vmsg("Mappability calculation complete.")
  vmsg("Results summary:")
  vmsg("  Mappable bins: ", n_mappable)
  vmsg("  Unmappable bins: ", n_unmappable)
  vmsg("  NA bins: ", n_na)

  out
}
