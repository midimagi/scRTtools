# Table of replication-track lengths

Flatten a track `GRangesList` into one row per track.

## Usage

``` r
calc_rep_track_lengths(rep_tracks)
```

## Arguments

- rep_tracks:

  Named `GRangesList` (one element per cell), e.g. from
  [`define_replication_tracks()`](https://midimagi.github.io/scRTtools/reference/define_replication_tracks.md).

## Value

A `data.frame` with columns `chr`, `start`, `end`, `width_bp` followed
by the metadata columns of the tracks. If there are no tracks, an empty
`data.frame` with columns `cell`, `chr`, `start`, `end` and `width_bp`
(note: `cell`, not `cell_id`).
