# Call regions that differ between two states

Regions measured in both states are classified by the log2 fold change
`M` (treatment / reference): `<treatment>_up` when `M > thres`,
`<treatment>_down` when `M < -thres`, otherwise `no_sig`. Regions with a
value in only one state are appended as `<treatment>_up` (treatment
only) or `<treatment>_down` (reference only).

## Usage

``` r
call_diff_groups(df, id_col, feature, thres, states = c("ctr", "ama"))
```

## Arguments

- df:

  data.frame with the `id_col` column(s), a `state` column and
  `feature`.

- id_col:

  Name(s) of the identifier column(s).

- feature:

  Name of the value column.

- thres:

  Threshold on the absolute log2 fold change.

- states:

  The two state names, `c(reference, treatment)`.

## Value

A data.frame with `id_col` and the factor `iz_group` (levels
`<treatment>_up`, `<treatment>_down`, `no_sig`).

## Examples

``` r
long <- data.frame(id = rep(c("a", "b", "c"), 2),
                   state = rep(c("ctr", "ama"), each = 3),
                   usage = c(1, 2, NA, 4, 2, 5))
call_diff_groups(long, "id", "usage", thres = 1)
#>   id iz_group
#> 1  a   ama_up
#> 2  b   no_sig
#> 3  c   ama_up
```
