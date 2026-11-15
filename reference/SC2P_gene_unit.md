# Gene-level SC2P Phase I / Phase II summaries

For every expressed gene, counts not above `yy_thres` are fitted with a
robust zero-inflated Poisson (Phase I, the background). A cell is in
Phase II for a gene when its count reaches the 90% quantile of that
Poisson (threshold 1 when the fit is undefined). The function reports
the fraction of Phase I cells and, among Phase II cells, the mean and
coefficient of variation of counts-per-10,000 expression (library size =
the cell's total over the expressed genes of `Y`).

## Usage

``` r
SC2P_gene_unit(Y, yy_thres = 15)
```

## Arguments

- Y:

  Numeric matrix or data.frame of counts, genes in rows and cells in
  columns, with unique row and column names and no missing values. Genes
  with zero total are dropped; cells with zero total are an error.

- yy_thres:

  Counts above this value are excluded from the Phase I fit.

## Value

A data.frame with one row per expressed gene (row names = gene) and
columns `gene`, `phase_I_ratio`, `phase_II_mean`, `phase_II_varcoff`.
Genes with no Phase II cell have `NA` Phase II summaries; with a single
Phase II cell the coefficient of variation is `NA`.

## Examples

``` r
set.seed(1)
counts <- matrix(rpois(40 * 60, lambda = rep(c(0.3, 2), each = 20)),
                 nrow = 40,
                 dimnames = list(paste0("g", 1:40), paste0("c", 1:60)))
head(suppressWarnings(SC2P_gene_unit(counts)))
#>    gene phase_I_ratio phase_II_mean phase_II_varcoff
#> g1   g1     0.7000000     261.91088        0.4170683
#> g2   g2     0.6833333     261.34765        0.4034421
#> g3   g3     0.6000000     259.62834        0.2826473
#> g4   g4     0.9333333     437.84608        0.2317330
#> g5   g5     0.6500000     261.13258        0.4914219
#> g6   g6     0.0000000      78.23481        1.4845847
```
