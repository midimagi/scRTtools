# Summarise gene-expression dynamics along replication progress

Cells with finite progress are first ordered by increasing replication
progress. They are then grouped either into consecutive groups of
`ncells` cells or into `nbin_rep` equal-width progress intervals. Within
each non-empty group, method = "sc2p" applies
[`SC2P_gene_unit()`](https://midimagi.github.io/scRTtools/reference/SC2P_gene_unit.md)
with warnings suppressed and returns Phase I ratios and Phase II means.
Alternatively, method = "cp10k" normalises the entire matrix once by
each cell's total counts across all input genes and multiplies by
10,000. It then takes the arithmetic mean across cells in each
replication bin, including zero gene expression. No log transformation
or phase selection is applied. Supply the full gene count matrix before
selecting genes for downstream analysis so library sizes reflect all
genes. Sparse input remains sparse during normalisation.

## Usage

``` r
gex_tempo(
  gex_matrix,
  meta_data,
  var_use = "rep_percent",
  ncells = NULL,
  nbin_rep = NULL,
  yy_thres = 15,
  method = c("sc2p", "cp10k")
)
```

## Arguments

- gex_matrix:

  Numeric matrix, data.frame, or Matrix-like object with genes
  (`gene_id`) in rows and cells (`cell_id`) in columns. Values must be
  finite, non-negative integer counts. Row and column names must be
  unique.

- meta_data:

  Required data.frame with unique cell IDs as row names, containing
  every matrix cell and the replication-progress column.

- var_use:

  Column of meta_data containing 0–1 replication progress. Defaults to
  "rep_percent".

- ncells:

  Optional positive integer giving the number of consecutive
  progress-sorted cells per group. The final smaller group is retained.

- nbin_rep:

  Optional positive integer giving the number of equal-width
  replication-progress intervals. Empty intervals are retained as `NA`
  columns. `ncells` and `nbin_rep` are mutually exclusive. If no
  grouping parameter is supplied, each cell forms its own group.

- yy_thres:

  Count threshold passed to
  [`SC2P_gene_unit()`](https://midimagi.github.io/scRTtools/reference/SC2P_gene_unit.md);
  ignored for cp10k.

- method:

  Either "sc2p" (default) or "cp10k".

## Value

For sc2p, a named list with `df_phase_I` and `df_phase_II`. Both are
data frames with genes in rows and replication bins in columns. The list
has `rep_bins`, `cell_info`, `progress_source`, and `grouping`
attributes. In non-empty bins, genes with all-zero counts have Phase I
ratio 1 and Phase II mean NA. Empty bins retain NA for both summaries.
For cp10k, a list containing `df_cp10k_mean`, a gene-by-bin data.frame,
with the same grouping and cell metadata plus a `method` attribute.
Zero-library cells are excluded from CP10K means; bins with no positive
libraries retain NA. Otherwise all-zero genes have mean 0. The cp10k
rep_bins attribute also includes `n_cells_used`, the mean denominator.
All input counts are validated for cp10k, including cells with NA
progress. SC2P is not run when method = "cp10k".

## Details

Replication progress is read only from `meta_data[[var_use]]`, matched
by metadata row names to matrix column names. Values must be in 0–1, as
in
[`rep_tempo()`](https://midimagi.github.io/scRTtools/reference/rep_tempo.md).
Cells with NA/NaN progress are excluded; missing metadata cell IDs,
infinite/out-of-range progress, or no usable progress cause an error.
Matrix attributes and column ranks are never used as progress.

## Examples

``` r
set.seed(1)
counts <- matrix(rpois(30 * 40, 2), nrow = 30,
                 dimnames = list(paste0("g", 1:30), paste0("c", 1:40)))
meta <- data.frame(rep_percent = runif(40), row.names = colnames(counts))
res <- gex_tempo(counts, meta, ncells = 10, method = "cp10k")
dim(res$df_cp10k_mean)
#> [1] 30  4
attr(res, "rep_bins")
#>        rep_bin rep_percent   rep_min  rep_max n_cells n_cells_used
#> 1 rep_bin_0001    9.129881  1.854641 28.08001      10           10
#> 2 rep_bin_0002   37.614283 31.423464 44.56445      10           10
#> 3 rep_bin_0003   62.450985 46.938107 70.61567      10           10
#> 4 rep_bin_0004   87.315584 73.617154 99.62718      10           10
```
