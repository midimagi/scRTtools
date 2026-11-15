# Cell x IZ firing matrix from track centres

An IZ is fired in a cell when the centre of one of the cell's
IZ-assigned tracks (optionally only tracks shorter than
`max_track_width_bp`) falls inside it.

## Usage

``` r
build_firing_matrix(izones, tracks, meta_data, max_track_width_bp = NULL)
```

## Arguments

- izones:

  `GRanges` of IZs with an `iz_id` column (e.g.
  `make_adjacent_pairs()$izones`).

- tracks:

  Named `GRangesList` of tracks with columns `cell_id` (matching the
  names), `rep_percent` and `iz_id`.

- meta_data:

  A `data.frame` with columns `cell_id` (unique) and `rep_percent`;
  missing values are taken from the tracks.

- max_track_width_bp:

  `NULL` or the exclusive maximum track width.

## Value

A list with `fired` (sparse 0/1 `dgCMatrix`, cells x IZs) and
`cell_info` (`cell_id`, `rep_percent`, `n_tracks`, `n_retained_tracks`,
`n_fired_izones`, `n_tracks_hitting_iz`).
