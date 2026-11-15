# Annotate regions with a segmented signal track

Averages a signal track over genome tiles, segments the tile means with
`DNAcopy` (circular binary segmentation after outlier smoothing) and
writes each region's segment mean and a two-class label (segment mean
`> 0` or not) back to `df`. Typical use: A/B compartment or early/late
labels from an eigenvector or log-ratio track.

## Usage

``` r
annotate_track_segments(
  track_gr,
  df,
  value_col,
  label_col,
  labels,
  align_tiles,
  tilewidth = 1e+05,
  verbose = FALSE
)
```

## Arguments

- track_gr:

  `GRanges` with a numeric `score` column and seqlengths.

- df:

  data.frame with a `region_id` column (`chr:start-end`).

- value_col, label_col:

  Names of the columns to write.

- labels:

  Two labels: for segment mean `> 0` and otherwise.

- align_tiles:

  `TRUE`: tiles without signal get `NA`; `FALSE`: every tile must have
  signal.

- tilewidth:

  Tile width in bp.

- verbose:

  Let
  [`DNAcopy::segment()`](https://rdrr.io/pkg/DNAcopy/man/segment.html)
  print its progress.

## Value

`df` with `value_col` and `label_col` added (`NA` for regions
overlapping no segment; when a region overlaps several segments the last
one wins).

## Details

[`DNAcopy::segment()`](https://rdrr.io/pkg/DNAcopy/man/segment.html)
uses random permutations; set a seed beforehand for reproducible
segments.
