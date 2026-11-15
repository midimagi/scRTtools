# Model matrix without linearly dependent columns

Builds `stats::model.matrix(formula, data)` and keeps a full-rank subset
of its columns (QR pivoting; original column order is kept). Useful when
factor levels or spline bases are empty or collinear within one sample.

## Usage

``` r
full_rank_model_matrix(formula, data)
```

## Arguments

- formula:

  One-sided model formula. Functions used in it (e.g.
  [`splines::ns`](https://rdrr.io/r/splines/ns.html)) must be visible
  from the formula's environment.

- data:

  data.frame with the variables of `formula`; missing or non-finite
  values are an error.

## Value

A numeric matrix with one row per row of `data`.

## Examples

``` r
d <- data.frame(x = rnorm(30), g = gl(3, 10))
d$x2 <- 2 * d$x
colnames(full_rank_model_matrix(~ x + x2 + g, d))
#> [1] "(Intercept)" "x"           "g2"          "g3"         
```
