# Fitted expression curves of selected genes and their clustering

Predicts each gene's GAM curve on an even grid over the observed
replication progress, z-scores the curves, clusters them (Ward's method
on `1 - Pearson correlation`) and summarises each curve.

## Usage

``` r
compute_gene_curves(
  group_log_expr,
  model_data,
  genes,
  gam_basis_k = 10L,
  use_sample_covariate = FALSE,
  n_curve_points = 101L,
  n_clusters = 3L
)
```

## Arguments

- group_log_expr:

  Gene-by-cell log-normalised matrix of one group.

- model_data:

  data.frame from
  [`prepare_group_inputs()`](https://midimagi.github.io/scRTtools/reference/prepare_group_inputs.md).

- genes:

  Genes to compute curves for; unknown genes are ignored.

- gam_basis_k:

  Upper bound of the smooth basis dimension.

- use_sample_covariate:

  Add `sample` as a covariate when it has more than one level.

- n_curve_points:

  Number of grid points.

- n_clusters:

  Number of curve clusters (capped at the number of curves).

## Value

A list with `curves` and `curve_z` (gene-by-grid matrices), `rep_grid`,
`clusters`, `hclust` (`NULL` for a single curve), `gene_order` and
`curve_summary` (`gene`, `cluster`, `peak_rep_percent`,
`late_minus_early`). Genes whose curve is non-finite or constant are
dropped. With no usable gene name the elements are `NULL` / empty.
