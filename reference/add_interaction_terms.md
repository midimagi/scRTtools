# Add regulator and local-replication interaction terms to a design matrix

Appends the columns of the per-gene model of the paired DNA-RNA
analysis: `regulator`, `regulator_progress`, `local_rep`,
`local_progress`, `local_progress2` and `interaction`
(`local_rep * regulator`).

## Usage

``` r
add_interaction_terms(B, regulator, local_rep, progress_c)
```

## Arguments

- B:

  Base design matrix (nuisance covariates), one row per cell.

- regulator:

  Regulator expression per cell (e.g. z-scored).

- local_rep:

  Local replication state per cell (0 = unreplicated, 1 = replicated,
  `NA` = missing).

- progress_c:

  Centred replication progress per cell.

## Value

`B` with the six columns appended.
