# Correlate heterogeneity metrics with epigenomic features

Correlate heterogeneity metrics with epigenomic features

## Usage

``` r
correlate_epigenomic_features(
  d,
  metrics,
  features,
  methods = c("spearman", "pearson")
)
```

## Arguments

- d:

  Per-IZ table with metric and feature columns.

- metrics, features:

  Column names.

- methods:

  Correlation methods.

## Value

A `data.frame` with one row per metric, feature and method
(`correlation`, `p_value`, `n`, BH `FDR`; `NA` with fewer than 10
values).
