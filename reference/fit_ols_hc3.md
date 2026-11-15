# Ordinary least squares with HC3 robust standard errors

Fits `y ~ X` by QR and reports selected coefficients with
heteroscedasticity-consistent (HC3) standard errors and t-based p-values
(`n - k` degrees of freedom). Rows with a non-finite `y` or `X` are
dropped first.

## Usage

``` r
fit_ols_hc3(X, y, terms = NULL, max_leverage = 0.99)
```

## Arguments

- X:

  Numeric design matrix (include the intercept column yourself).

- y:

  Numeric response, one value per row of `X`.

- terms:

  Column names of `X` to report; `NULL` reports the last column.

- max_leverage:

  Largest leverage accepted. `Inf` never rejects a fit on leverage.

## Value

A data.frame with one row per term: `term`, `beta`, `se`, `p`, `n`,
`max_hat`, `status`.

## Details

The fit is not attempted (`status = "rank_or_n"`) when `n <= k + 20` or
`X` is rank deficient. When the largest leverage reaches `max_leverage`
the result is `NA` with `status = "high_leverage"`; otherwise leverages
are capped at `1 - 1e-8`. A term whose standard error is zero or not
finite gets `status = "zero_se"`.

## Examples

``` r
set.seed(1)
X <- cbind(`(Intercept)` = 1, age = rnorm(200), dose = rnorm(200))
y <- 0.5 * X[, "dose"] + rnorm(200, sd = 1 + abs(X[, "age"]))
fit_ols_hc3(X, y, terms = c("age", "dose"))
#>   term      beta        se           p   n    max_hat status
#> 1  age 0.2153581 0.2220164 0.333231020 200 0.04744096     ok
#> 2 dose 0.5017496 0.1584040 0.001782137 200 0.04744096     ok
```
