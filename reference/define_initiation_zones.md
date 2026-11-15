# Define single-cell replication initiation zones

Cluster single-cell replication tracks into initiation zones (IZs)
seeded by replication-timing peaks (e.g.
[`find_rt_peaks()`](https://midimagi.github.io/scRTtools/reference/find_rt_peaks.md)
on
[`count_track_centers()`](https://midimagi.github.io/scRTtools/reference/count_track_centers.md)
output). Per chromosome, tracks are processed from short to long. Tracks
flagged as bounded by a gap (`flag != 0`) or longer than
`max_track_length` are skipped. Each track, widened by `2.5 * bin_size`
on both sides, either seeds a new putative IZ, joins the single putative
IZ it overlaps (by more than `min_overlap`, plus one bin for tracks
longer than `short_track_th`; if two are overlapped the one with a
compatible track-centre interval is used), or is skipped (more than two
overlaps, a track of the same cell already in the IZ, or incompatible
centres). Short tracks update the putative IZ interval. Putative IZs
supported by fewer than `min_usage` tracks are removed, all tracks are
reassigned once more, and IZs with fewer than `min_usage` reassigned
tracks are removed.

## Usage

``` r
define_initiation_zones(
  tracks,
  rt_peaks,
  short_track_th = "500Kb",
  max_track_length = "1Mb",
  min_overlap = "0bp",
  min_usage = 3L,
  bin_size = "20Kb",
  BPPARAM = BiocParallel::SerialParam(),
  verbose = FALSE
)
```

## Arguments

- tracks:

  Named `GRangesList` of tracks (one element per cell, unique names)
  with a `flag` metadata column (e.g.
  [`define_replication_tracks()`](https://midimagi.github.io/scRTtools/reference/define_replication_tracks.md)).

- rt_peaks:

  `GRanges` of seed peaks.

- short_track_th:

  Tracks up to this length define IZ widths.

- max_track_length:

  Longer tracks are not used.

- min_overlap:

  Minimum overlap between a widened track and an IZ.

- min_usage:

  Minimum number of supporting tracks of an IZ (default 3; lower it for
  small datasets).

- bin_size:

  Bin size of the tracks.

- BPPARAM:

  A
  [`BiocParallel::BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html);
  chromosomes are processed in parallel.

- verbose:

  Logical; print progress messages. Default `FALSE`.

## Value

A list with `tracks` (input tracks with an `iz_id` metadata column, `NA`
if unassigned) and `izones` (sorted `GRanges` with `iz_id`, `q5`, `q25`,
`q50`, `q75`, `q95` and `n_tracks`).

## Details

Each final IZ spans the 5%-95% quantiles of its track centres padded by
half a bin on each side (the "fixed" version:
`start = round(q5 - bin_size / 2)`, `width = q95 - q5 + bin_size`). IZs
are numbered by genomic order.

Lengths can be given as numbers (bp) or strings such as `"500Kb"`,
`"1Mb"`, `"0bp"`.
