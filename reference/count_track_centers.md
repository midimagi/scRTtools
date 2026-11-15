# Count replication-track centres in genomic bins

Optionally filter each cell's tracks with an expression evaluated on
their metadata columns, then count, for every range of `granges`, the
tracks whose centre `floor((start + end) / 2)` falls inside it.

## Usage

``` r
count_track_centers(tracks, granges, track_filter = NULL)
```

## Arguments

- tracks:

  Named `GRangesList` of tracks (one element per cell), e.g. from
  [`define_replication_tracks()`](https://midimagi.github.io/scRTtools/reference/define_replication_tracks.md).

- granges:

  `GRanges` of the bins to count in.

- track_filter:

  Optional unquoted expression evaluated against `mcols(tracks[[cell]])`
  (non-standard evaluation), e.g. `(flag == 0) & (num_bins > 50)` or
  `!is.na(iz_id)`. It must give a logical scalar or one value per track;
  `NA` counts as `FALSE`.

## Value

`granges`, sorted, with an integer `score` column (number of track
centres; zero-count ranges are kept).
