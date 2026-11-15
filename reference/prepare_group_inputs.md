# Collect the cells and matrices of one analysis group

Collect the cells and matrices of one analysis group

## Usage

``` r
prepare_group_inputs(
  group_name,
  meta_data,
  log_expr,
  normalized_expr,
  group_by_col = "state",
  min_cells_per_group = 50L
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

## Value

A list with `group_meta` (sorted by `rep_percent`), `model_data`
(`rep_percent`, `sample`), `log_expr`, `normalized_expr` and `cells`.
