# Fit the progression GAM and its null model for one gene

Fits `expression ~ s(rep_percent, k)` (optionally `+ sample`) with
[`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) (Gaussian, GCV) and
the matching null model without the smooth.
`k = max(4, min(gam_basis_k, n_unique_rep - 1))`.

## Usage

``` r
build_gene_models(
  y,
  model_data,
  gam_basis_k = 10L,
  use_sample_covariate = FALSE
)
```

## Arguments

- y:

  Numeric expression vector, one value per row of `model_data`.

- model_data:

  data.frame with `rep_percent` and a factor `sample`.

- gam_basis_k:

  Upper bound of the smooth basis dimension.

- use_sample_covariate:

  Add `sample` as a covariate when it has more than one level.

## Value

A list with `gam`, `null`, `model_data` (with `expression` added) and
`sample_is_used`. Errors when there are fewer than five unique
`rep_percent` values or the expression is constant.
