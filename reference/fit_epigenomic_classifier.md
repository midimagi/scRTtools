# Logistic classifier of high versus low heterogeneity

IZs in the bottom / top quartile of the composite heterogeneity (mean of
the scaled metrics) are labelled 0 / 1 and modelled by logistic
regression on the transformed features plus `width`, `S50`, `mean_usage`
and `chr`.

## Usage

``` r
fit_epigenomic_classifier(d, metrics, features)
```

## Arguments

- d:

  Per-IZ table with the metrics, features, `width`, `S50`, `mean_usage`
  and the covariates.

- metrics, features:

  Column names.

## Value

A list with `model`, `coefficients` (with odds ratios, training AUC and
BH `FDR`) and `labeled_data`.
