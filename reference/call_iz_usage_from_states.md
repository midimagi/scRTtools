# Cell x IZ usage matrix from replication states

Read the state of the genome bin containing each IZ centre: `1` if
replicated, `0` if unreplicated, `NA` otherwise.

## Usage

``` r
call_iz_usage_from_states(
  genomic_states,
  cell_ids,
  iz_metadata,
  replicated_state = 4L,
  unreplicated_state = 2L,
  bin_size = 20000L
)
```

## Arguments

- genomic_states:

  Matrix of states, bins x cells, with row names `"chr-start-end"` of
  `bin_size` bins (1-based starts).

- cell_ids:

  Cells to use (columns of `genomic_states`).

- iz_metadata:

  `prepare_iz_data()$iz_metadata` (`IZ_id`, `chr`, `center`).

- replicated_state, unreplicated_state:

  State codes.

- bin_size:

  Bin size (bp).

## Value

An integer matrix, cells x IZs.
