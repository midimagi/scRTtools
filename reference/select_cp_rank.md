# Choose the rank of a non-negative CP decomposition

Fit
[`fit_nonnegative_cp()`](https://midimagi.github.io/scRTtools/reference/fit_nonnegative_cp.md)
`n_init` times for every rank (seed `seed + 1000 * rank + run`). Among
ranks whose factor stability reaches `stability_threshold` (all ranks if
none does), the smallest rank whose mean error is within 10% of the
error range from the minimum is chosen.

## Usage

``` r
select_cp_rank(
  y,
  ranks = 2:6,
  n_init = 2L,
  seed = 1L,
  max_iter = 30L,
  tol = 1e-04,
  stability_threshold = 0.8,
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- y:

  Non-negative three-way array.

- ranks:

  Candidate ranks.

- n_init:

  Number of initialisations per rank.

- seed:

  Base seed.

- max_iter, tol:

  Passed to
  [`fit_nonnegative_cp()`](https://midimagi.github.io/scRTtools/reference/fit_nonnegative_cp.md).

- stability_threshold:

  Minimum stability.

- BPPARAM:

  A
  [`BiocParallel::BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  used over fits.

## Value

A list with `summary`, `chosen_rank`, `best_fits` (per rank), `all_fits`
and `jobs`.
