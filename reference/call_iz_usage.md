# Cell x IZ usage matrix from assigned tracks

Cell x IZ usage matrix from assigned tracks

## Usage

``` r
call_iz_usage(
  track_intervals = NULL,
  cell_ids = NULL,
  iz_ids = NULL,
  method = c("assigned_track", "existing"),
  existing_usage = NULL
)
```

## Arguments

- track_intervals:

  `prepare_iz_data()$track_intervals`.

- cell_ids, iz_ids:

  Row and column axes.

- method:

  `"assigned_track"` (1 where a cell has a track assigned to the IZ) or
  `"existing"` (`existing_usage != 0`).

- existing_usage:

  Matrix used with `method = "existing"`.

## Value

A sparse matrix, cells x IZs.
