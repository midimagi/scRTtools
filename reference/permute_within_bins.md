# Permute cells within sample and replication-progress bins

Returns a permutation of the cell indices that only exchanges cells of
the same group whose progress falls in the same bin of width
`bin_width`. Use it to break the DNA-RNA pairing while keeping sample
and progress (almost) fixed: `rna_columns[perm]` against the original
DNA columns. Set a seed beforehand for reproducible permutations.

## Usage

``` r
permute_within_bins(group, progress, bin_width)
```

## Arguments

- group:

  Sample (or other stratum) of each cell.

- progress:

  Replication progress of each cell.

- bin_width:

  Width of the progress bins.

## Value

Integer vector `perm` of the same length; cell `i` receives the data of
cell `perm[i]`.

## Examples

``` r
set.seed(1)
progress <- runif(20)
perm <- permute_within_bins(rep(c("a", "b"), 10), progress, 0.25)
max(abs(progress - progress[perm])) < 0.25
#> [1] TRUE
```
