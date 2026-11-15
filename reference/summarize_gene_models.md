# Summarise the models of one gene

Summarise the models of one gene

## Usage

``` r
summarize_gene_models(model_list)
```

## Arguments

- model_list:

  Output of
  [`build_gene_models()`](https://midimagi.github.io/scRTtools/reference/build_gene_models.md).

## Value

A one-row data.frame: smooth-term `p_value`, `edf`, `f_statistic`,
`adjusted_r_squared`, `deviance_explained`, `aic_gam`, `aic_null` and
`sample_covariate_used`.
