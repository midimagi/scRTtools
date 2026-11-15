# Sample cells evenly across replication progress

Sample cells evenly across replication progress

## Usage

``` r
sample_cells_by_progress(
  cell_metadata,
  n_cells = 250L,
  n_bins = 10L,
  seed = 1L
)
```

## Arguments

- cell_metadata:

  Table with `cell_id` and `replication_progress`.

- n_cells:

  Number of cells to keep (all cells if there are fewer).

- n_bins:

  Number of progress bins sampled from.

- seed:

  Seed.

## Value

Character vector of cell ids.
