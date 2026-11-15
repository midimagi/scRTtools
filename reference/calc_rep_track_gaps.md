# Gaps between adjacent replication tracks

For each cell, sort the tracks and compute the distance between
consecutive tracks on the same chromosome (`start[i + 1] - end[i] - 1`).

## Usage

``` r
calc_rep_track_gaps(rep_tracks)
```

## Arguments

- rep_tracks:

  Named `GRangesList` (one element per cell), e.g. from
  [`define_replication_tracks()`](https://midimagi.github.io/scRTtools/reference/define_replication_tracks.md),
  with metadata columns `cell_id` and `rep_percent`.

## Value

A `data.frame` with columns `cell_id`, `rep_percent`, `chr` and
`gap_bp`, one row per adjacent pair.
