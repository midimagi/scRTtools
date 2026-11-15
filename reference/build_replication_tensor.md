# Dense cell x IZ x position replication tensor

`1` where a cell's merged track covers position `center + rel_bp` of the
IZ (`rel_bp` from `-window_bp` to `window_bp` by `bin_size`), else `0`.

## Usage

``` r
build_replication_tensor(
  track_intervals,
  cell_ids,
  iz_metadata,
  iz_ids,
  window_bp = 500000L,
  bin_size = 20000L,
  max_elements = 5e+07
)
```

## Arguments

- track_intervals:

  `prepare_iz_data()$track_intervals`.

- cell_ids:

  Cells (first dimension).

- iz_metadata:

  `prepare_iz_data()$iz_metadata`.

- iz_ids:

  IZs (second dimension).

- window_bp, bin_size:

  Window half-width and step (bp).

- max_elements:

  Maximum number of tensor elements.

## Value

A numeric array, cells x IZs x positions.
