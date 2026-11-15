# Identify S-phase progression markers in one group of cells

Runs
[`prepare_group_inputs()`](https://midimagi.github.io/scRTtools/reference/prepare_group_inputs.md),
[`fit_gene_gams()`](https://midimagi.github.io/scRTtools/reference/fit_gene_gams.md),
[`classify_progression_markers()`](https://midimagi.github.io/scRTtools/reference/classify_progression_markers.md)
and
[`compute_gene_curves()`](https://midimagi.github.io/scRTtools/reference/compute_gene_curves.md)
(on the marker genes) for one analysis group.

## Usage

``` r
analyze_one_group(
  group_name,
  meta_data,
  log_expr,
  normalized_expr,
  group_by_col = "state",
  min_cells_per_group = 50L,
  min_normalized_expression = 0.1,
  min_cell_fraction = 0.2,
  fdr_threshold = 0.01,
  aic_rule = c("gam_lower", "paper_literal_gam_greater"),
  gam_basis_k = 10L,
  use_sample_covariate = FALSE,
  n_curve_points = 101L,
  n_clusters = 3L,
  BPPARAM = BiocParallel::SerialParam(),
  verbose = FALSE
)
```

## Arguments

- group_name:

  Value of `group_by_col` to analyse.

- meta_data:

  Output of
  [`prepare_metadata()`](https://midimagi.github.io/scRTtools/reference/prepare_metadata.md).

- log_expr, normalized_expr:

  Gene-by-cell matrices from
  [`normalize_sparse_counts()`](https://midimagi.github.io/scRTtools/reference/normalize_sparse_counts.md).

- group_by_col:

  Metadata column defining the analysis groups.

- min_cells_per_group:

  Minimum number of eligible cells; fewer is an error.

- min_normalized_expression:

  A gene counts as expressed in a cell when its normalised expression is
  at least this value.

- min_cell_fraction:

  Minimum fraction of cells expressing the gene.

- fdr_threshold:

  FDR cut-off for the smooth term.

- aic_rule:

  `"gam_lower"`: the GAM must have the lower AIC;
  `"paper_literal_gam_greater"`: the opposite.

- gam_basis_k:

  Upper bound of the smooth basis dimension.

- use_sample_covariate:

  Add `sample` as a covariate when it has more than one level.

- n_curve_points:

  Number of grid points.

- n_clusters:

  Number of curve clusters (capped at the number of curves).

- BPPARAM:

  A `BiocParallel` back-end; serial by default.

- verbose:

  Print a progress message. Default `FALSE`.

## Value

A list with `analysis_group`, `gam_results` (all genes), `marker_table`
(markers with a valid curve, with `cluster`, `peak_rep_percent`,
`late_minus_early`, sorted by them and `fdr`), `group_meta`, `cells`,
and the `curves`, `curve_z`, `rep_grid`, `clusters`, `hclust`,
`gene_order` elements of
[`compute_gene_curves()`](https://midimagi.github.io/scRTtools/reference/compute_gene_curves.md).
