# Poisson regression with HC3 robust standard errors

Fits a log-link Poisson model with an offset by
[`stats::glm.fit`](https://rdrr.io/r/stats/glm.html) and reports
selected coefficients with HC3 sandwich standard errors, which tolerate
a variance different from the mean, and normal p-values.

## Usage

``` r
fit_poisson_hc3(
  X,
  y,
  offset,
  terms = NULL,
  maxit = 100,
  min_rcond = 1e-12,
  max_abs_beta = 20
)
```

## Arguments

- X:

  Numeric design matrix (include the intercept column yourself).

- y:

  Count response.

- offset:

  Offset on the log scale, e.g. `log(library_size / 1e4)`.

- terms:

  Column names of `X` to report; `NULL` reports the last column.

- maxit:

  Maximum number of IWLS iterations.

- min_rcond, max_abs_beta:

  Thresholds of the `"ill_conditioned"` and `"unstable"` checks.

## Value

A data.frame with one row per term: `term`, `beta`, `se`, `p`, `n`,
`status`.

## Details

Status codes: `"rank_or_n"` (`n <= k + 20` or rank-deficient `X`),
`"nonconverged"`, `"ill_conditioned"` (reciprocal condition number of
`X' diag(mu) X` below `min_rcond`), `"unstable"` (standard error zero or
not finite, or `|beta| > max_abs_beta`), `"ok"`.

## Examples

``` r
set.seed(1)
X <- cbind(`(Intercept)` = 1, dose = rnorm(300))
depth <- runif(300, 5e3, 2e4)
y <- rpois(300, depth / 1e4 * exp(0.2 + 0.3 * X[, "dose"]))
fit_poisson_hc3(X, y, offset = log(depth / 1e4))
#>   term      beta        se            p   n status
#> 1 dose 0.3161844 0.0443064 9.586725e-13 300     ok
```
