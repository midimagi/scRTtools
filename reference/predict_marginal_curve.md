# Predict the fitted expression curve of one gene

Fits the GAM of
[`build_gene_models()`](https://midimagi.github.io/scRTtools/reference/build_gene_models.md)
and predicts on `rep_grid`. With a sample covariate the prediction is
averaged over samples.

## Usage

``` r
predict_marginal_curve(
  y,
  model_data,
  rep_grid,
  gam_basis_k = 10L,
  use_sample_covariate = FALSE
)
```

## Arguments

- y:

  Numeric expression vector, one value per row of `model_data`.

- model_data:

  data.frame with `rep_percent` and a factor `sample`.

- rep_grid:

  Replication-progress values to predict at.

- gam_basis_k:

  Upper bound of the smooth basis dimension.

- use_sample_covariate:

  Add `sample` as a covariate when it has more than one level.

## Value

Numeric vector, one value per element of `rep_grid`.
