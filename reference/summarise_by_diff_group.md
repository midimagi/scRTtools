# Per-state values of regions with their difference group

Per-state values of regions with their difference group

## Usage

``` r
summarise_by_diff_group(
  df,
  id_col,
  value_col,
  diff_groups,
  states = c("ctr", "ama")
)
```

## Arguments

- df:

  data.frame with the `id_col` column(s), a `state` column and
  `feature`.

- id_col:

  Name(s) of the identifier column(s).

- value_col:

  Name of the value column to spread by state.

- diff_groups:

  Output of
  [`call_diff_groups()`](https://midimagi.github.io/scRTtools/reference/call_diff_groups.md)
  (unique `id_col`).

- states:

  The two state names, `c(reference, treatment)`.

## Value

The wide table of `value_col` (regions measured in both states) with the
columns of `diff_groups` joined; rows with `NA` are dropped.
