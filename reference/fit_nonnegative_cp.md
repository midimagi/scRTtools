# Non-negative CP decomposition of a three-way tensor

Multiplicative-update non-negative CANDECOMP/PARAFAC with random uniform
initialisation. The relative error is checked every `check_every`
iterations; the fit stops when it changes by less than `tol`.

## Usage

``` r
fit_nonnegative_cp(
  y,
  rank,
  seed = 1L,
  max_iter = 50L,
  tol = 1e-04,
  check_every = 5L,
  epsilon = 1e-10
)
```

## Arguments

- y:

  Non-negative array with three dimensions and no `NA`.

- rank:

  Number of components.

- seed:

  Seed of the initialisation.

- max_iter:

  Maximum number of iterations.

- tol:

  Convergence tolerance.

- check_every:

  Error-check interval.

- epsilon:

  Lower bound of denominators.

## Value

A list with unit-norm factor matrices `A`, `B`, `C`, `lambda`, `error`
(relative), `explained`, `iterations`, `converged`, `history` and
`rank`.
