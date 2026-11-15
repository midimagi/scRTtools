# Library-size normalise a sparse count matrix

Library-size normalise a sparse count matrix

## Usage

``` r
normalize_sparse_counts(counts, scale_factor = 10000)
```

## Arguments

- counts:

  Sparse gene-by-cell count matrix (`dgCMatrix`).

- scale_factor:

  Counts are scaled to this total per cell.

## Value

A list with `normalized` (counts per `scale_factor`), `log_normalized`
(`log1p` of it) and `library_size`. Cells with a zero or non-finite
library size are an error.

## Examples

``` r
m <- Matrix::rsparsematrix(20, 10, 0.4, rand.x = function(n) rpois(n, 3) + 1)
dimnames(m) <- list(paste0("g", 1:20), paste0("c", 1:10))
str(normalize_sparse_counts(m), max.level = 1)
#> List of 3
#>  $ normalized    :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>  $ log_normalized:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>  $ library_size  : Named num [1:10] 29 40 42 35 26 15 44 32 33 23
#>   ..- attr(*, "names")= chr [1:10] "c1" "c2" "c3" "c4" ...
```
