# Remove the replication-progress trend from a signal matrix

Fit the expected value of each column as a function of replication
progress, either as the mean within progress bins (`"bin"`) or with a
natural spline of `spline_df` degrees of freedom (`"spline"`, falling
back to `"bin"` with fewer than 4 cells or 3 distinct progress values),
and return the residuals.

## Usage

``` r
correct_signal_matrix(
  x,
  progress,
  method = c("bin", "spline"),
  n_progress_bins = 10L,
  spline_df = 4L
)
```

## Arguments

- x:

  Numeric matrix, cells x features.

- progress:

  Replication progress of the cells.

- method:

  `"bin"` or `"spline"`.

- n_progress_bins:

  Number of progress bins.

- spline_df:

  Spline degrees of freedom.

## Value

A list with `residual`, `fitted` and `method`.
