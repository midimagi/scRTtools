# Format base-pair numbers as readable labels

Format base-pair numbers as readable labels

## Usage

``` r
format_bp(x)
```

## Arguments

- x:

  Numeric vector of lengths in bp.

## Value

Character vector such as `"20 kb"`, `"1.5 Mb"`; `NA` stays `NA`.

## Examples

``` r
format_bp(c(500, 20000, 1.5e6, 2.7e9, NA))
#> [1] "500 bp" "20 kb"  "1.5 Mb" "2.7 Gb" NA      
```
