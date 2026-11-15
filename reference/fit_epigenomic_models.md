# Linear models of heterogeneity metrics on epigenomic features

For each metric, fit `scale(metric) ~ features + covariates` with
features transformed by signed `log1p` and scaled, and `width` (log1p),
`S50` and `mean_usage` scaled.

## Usage

``` r
fit_epigenomic_models(
  d,
  metrics,
  features,
  covariates = c("width", "S50", "mean_usage", "chr")
)
```

## Arguments

- d:

  Per-IZ table with the metrics, features, `width`, `S50`, `mean_usage`
  and the covariates.

- metrics, features:

  Column names.

- covariates:

  Covariate column names.

## Value

A `data.frame` of feature coefficients per metric with BH `FDR`.
