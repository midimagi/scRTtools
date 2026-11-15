# Pairs of IZs separated by `n_adj` IZs

Order IZs by chromosome (seqlevel order), start and end, and pair each
IZ with the IZ `n_adj + 1` positions to its right on the same chromosome
(`n_adj = 0`: adjacent IZs).

## Usage

``` r
make_pair_init(izones, n_adj = 0L)
```

## Arguments

- izones:

  `GRanges` of IZs with an `iz_id` column.

- n_adj:

  Number of IZs between the two members of a pair.

## Value

A `data.frame` with `pair_id` (factor, `"left_right"`), `chr`,
`left_iz_id` and `right_iz_id` (character).
