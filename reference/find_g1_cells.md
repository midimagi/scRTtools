# Identify G1-phase cells from coverage and scaled MAPD

Classify cells as G1 using the relationship between coverage and scaled
MAPD. A two-component linear model is fitted with a hard-assignment EM
(internal `.find_g1_coefficients_em()`); the component with the lower
predicted MAPD at the maximum coverage is taken as G1, and cells whose
residual from that line is `< 0.05` are classified as G1. With fewer
than 4 cells a single linear regression is used instead (with a
warning).

## Usage

``` r
find_g1_cells(meta_data, by_sample = FALSE, seed = NULL, verbose = FALSE)
```

## Arguments

- meta_data:

  A `data.frame` with columns `cell_id`, `coverage`, `scaled_mapd` and
  `sample` (e.g. output of
  [`calculate_mapd()`](https://midimagi.github.io/scRTtools/reference/calculate_mapd.md));
  row names must equal `cell_id`.

- by_sample:

  Logical; fit the model separately for each `sample`. We use `FALSE`.

- seed:

  Optional integer. The EM breaks likelihood ties and re-seeds
  degenerate assignments at random; if `seed` is given, this runs inside
  `withr::with_seed(seed, ...)`. With `NULL` (default) the current RNG
  state is used.

- verbose:

  Logical; print EM iteration messages. Default `FALSE`.

## Value

`meta_data` **restricted to rows with non-`NA` `coverage` and
`scaled_mapd`**, with an added/updated logical column `is_G1`.

## Examples

``` r
set.seed(1)
cov <- runif(200, 0.5, 3)
g1 <- rep(c(TRUE, FALSE), c(120, 80))
smapd <- ifelse(g1, 0.2 - 0.02 * cov, 0.5 - 0.02 * cov) + rnorm(200, 0, 0.01)
md <- data.frame(cell_id = paste0("c", 1:200), coverage = cov,
                 scaled_mapd = smapd, sample = "s1")
rownames(md) <- md$cell_id
res <- find_g1_cells(md, seed = 1, verbose = FALSE)
table(res$is_G1, g1)
#>        g1
#>         FALSE TRUE
#>   FALSE    80    0
#>   TRUE      0  120
```
