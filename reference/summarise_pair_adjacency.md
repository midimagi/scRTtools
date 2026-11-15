# Per-pair co-firing of adjacent IZs across selected cells

Per-pair co-firing of adjacent IZs across selected cells

## Usage

``` r
summarise_pair_adjacency(fired, pairs, selected_cells, label)
```

## Arguments

- fired:

  Output of
  [`build_firing_matrix()`](https://midimagi.github.io/scRTtools/reference/build_firing_matrix.md)`$fired`.

- pairs:

  `make_adjacent_pairs()$pairs`.

- selected_cells:

  Row indices (or names) of `fired` to use.

- label:

  Value of the `rep_stage` column.

## Value

`pairs` with `rep_stage` and the co-firing statistics of each pair over
the selected cells.
