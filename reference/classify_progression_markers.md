# Classify genes as S-phase progression markers

Adds BH-adjusted `fdr` (non-finite p-values count as 1), the AIC
difference `delta_aic_null_minus_gam` and the flags `aic_pass`,
`expression_pass` (`expressed_fraction >= min_cell_fraction`), `dynamic`
(`fdr < fdr_threshold` and `aic_pass`) and `marker` (`dynamic` and
`expression_pass`).

## Usage

``` r
classify_progression_markers(
  fit_table,
  group_name,
  n_cells,
  min_cell_fraction = 0.2,
  fdr_threshold = 0.01,
  aic_rule = c("gam_lower", "paper_literal_gam_greater")
)
```

## Arguments

- fit_table:

  Output of
  [`fit_gene_gams()`](https://midimagi.github.io/scRTtools/reference/fit_gene_gams.md).

- group_name, n_cells:

  Group label and cell number, stored in the `analysis_group` and
  `n_cells` columns.

- min_cell_fraction:

  Minimum fraction of cells expressing the gene.

- fdr_threshold:

  FDR cut-off for the smooth term.

- aic_rule:

  `"gam_lower"`: the GAM must have the lower AIC;
  `"paper_literal_gam_greater"`: the opposite.

## Value

`fit_table` with the added columns, `gene`, `analysis_group`, `n_cells`
first, sorted by `fdr` then decreasing AIC difference.
