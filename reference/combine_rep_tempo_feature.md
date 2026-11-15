# Combine the replication-tempo fits of two states

Stacks the per-region results of
[`fit_rep_tempo()`](https://midimagi.github.io/scRTtools/reference/fit_rep_tempo.md)
for two states and derives `twidth = t75 - t25`.

## Usage

``` r
combine_rep_tempo_feature(ctr_df, ama_df, states = c("ctr", "ama"))
```

## Arguments

- ctr_df, ama_df:

  Fit tables of the reference and the treatment state, with columns
  `region_id`, `t25`, `t50`, `t75`.

- states:

  The two state names, `c(reference, treatment)`.

## Value

A data.frame with `region_id`, `t50`, `twidth` and the factor `state`;
rows with `NA` are dropped.
