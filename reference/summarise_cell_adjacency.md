# Per-cell co-firing of adjacent IZs

For each cell, cross-tabulate the fired state of the left and right IZ
of every adjacent pair and derive co-firing statistics (`p11`, expected
`p_left * p_right`, excess, lift, phi, conditional probabilities). The
chromosome-shuffle expectation is the expected number of fired adjacent
pairs if the cell's fired IZs were placed at random within each
chromosome.

## Usage

``` r
summarise_cell_adjacency(fired, pairs, iz_features)
```

## Arguments

- fired:

  Output of
  [`build_firing_matrix()`](https://midimagi.github.io/scRTtools/reference/build_firing_matrix.md)`$fired`.

- pairs, iz_features:

  `make_adjacent_pairs()$pairs` and `$features`.

## Value

A `data.frame`, one row per cell (rows of `fired`), with the 2 x 2
counts, probabilities and `p_expected_chr_shuffle`,
`cofire_lift_chr_shuffle`, `p_excess_chr_shuffle`.
