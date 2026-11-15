# Firing range of initiation zones

For each IZ with at least `min_rep_tracks` assigned tracks, the earliest
firing time is the smallest `rep_percent` among cells with a track
assigned to the IZ, and the latest unreplicated time is the largest
`rep_percent` among cells whose unreplicated segment (run of bins in
state `unrep_state`) covers the IZ centre, using only the segments that
cover the fewest IZ centres. Each bound is corroborated when the
second-best value is within 0.1 of it. IZs with
`latest_unrep_percent < earliest_rep_percent` get `NA` bounds and
`corroborated = FALSE`.

## Usage

``` r
calculate_firing_range(
  gr_izones,
  genobins_rep_states,
  replication_tracks,
  unrep_state = 2,
  min_rep_tracks = 2,
  verbose = FALSE
)
```

## Arguments

- gr_izones:

  `GRanges` of IZs with an `n_tracks` column (e.g.
  `define_initiation_zones()$izones`); IZ `i` must have `iz_id == i`.

- genobins_rep_states:

  Matrix of states, bins x cells, with row names `"chr-start-end"` and
  cell IDs as column names; every IZ centre must lie in exactly one bin.

- replication_tracks:

  `CompressedGRangesList` of tracks with metadata columns `iz_id`,
  `cell_id` and `rep_percent` (e.g. `define_initiation_zones()$tracks`);
  every matrix column must have a `rep_percent`.

- unrep_state:

  State of unreplicated bins.

- min_rep_tracks:

  Minimum `n_tracks` of an IZ.

- verbose:

  Logical; print progress messages. Default `FALSE`.

## Value

`gr_izones` with columns `earliest_rep_percent`, `latest_unrep_percent`
and `corroborated`.
