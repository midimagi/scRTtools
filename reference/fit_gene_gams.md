# Fit progression GAMs for every gene

Applies
[`build_gene_models()`](https://midimagi.github.io/scRTtools/reference/build_gene_models.md)
and
[`summarize_gene_models()`](https://midimagi.github.io/scRTtools/reference/summarize_gene_models.md)
to each row of `group_log_expr`. A gene whose fit fails gets `NA`
statistics and the error message in `fit_status`.

## Usage

``` r
fit_gene_gams(
  group_log_expr,
  group_normalized_expr,
  model_data,
  min_normalized_expression = 0.1,
  gam_basis_k = 10L,
  use_sample_covariate = FALSE,
  BPPARAM = BiocParallel::SerialParam(),
  verbose = FALSE
)
```

## Arguments

- group_log_expr, group_normalized_expr:

  Gene-by-cell matrices of one group, columns aligned with the rows of
  `model_data`.

- model_data:

  data.frame from
  [`prepare_group_inputs()`](https://midimagi.github.io/scRTtools/reference/prepare_group_inputs.md).

- min_normalized_expression:

  A gene counts as expressed in a cell when its normalised expression is
  at least this value.

- gam_basis_k:

  Upper bound of the smooth basis dimension.

- use_sample_covariate:

  Add `sample` as a covariate when it has more than one level.

- BPPARAM:

  A `BiocParallel` back-end; serial by default.

- verbose:

  Print a progress message. Default `FALSE`.

## Value

A data.frame with one row per gene: the columns of
[`summarize_gene_models()`](https://midimagi.github.io/scRTtools/reference/summarize_gene_models.md),
`fit_status`, `gene` and `expressed_fraction`.
