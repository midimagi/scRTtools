# Row-wise z-scores, dropping constant rows

Row-wise z-scores, dropping constant rows

## Usage

``` r
scale_rows(x)
```

## Arguments

- x:

  Numeric matrix.

## Value

A list with `matrix` (z-scored rows with a finite, positive standard
deviation) and `valid` (logical, which rows of `x` were kept).

## Examples

``` r
scale_rows(rbind(a = c(1, 2, 3), b = c(5, 5, 5)))
#> $matrix
#>   [,1] [,2] [,3]
#> a   -1    0    1
#> 
#> $valid
#>     a     b 
#>  TRUE FALSE 
#> 
```
