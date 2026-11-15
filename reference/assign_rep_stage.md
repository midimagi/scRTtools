# Replication stage labels from replication percentage

`cut(rep_percent, breaks, labels, right = FALSE, include.lowest = TRUE)`,
with `rep_percent == 1` in the last stage and `NA` as `"unknown"`.

## Usage

``` r
assign_rep_stage(
  rep_percent,
  breaks = c(0, 0.25, 0.5, 0.75, 1),
  labels = c("Q1", "Q2", "Q3", "Q4")
)
```

## Arguments

- rep_percent:

  Numeric vector of replicated fractions.

- breaks, labels:

  Stage breaks and labels (e.g. quartiles `Q1`-`Q4`).

## Value

A factor with levels `c(labels, "unknown")`.
