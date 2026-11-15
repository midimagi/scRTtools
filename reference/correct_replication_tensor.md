# Remove the replication-progress trend from a replication tensor

Apply
[`correct_signal_matrix()`](https://midimagi.github.io/scRTtools/reference/correct_signal_matrix.md)
to the tensor flattened to cells x (IZ, position).

## Usage

``` r
correct_replication_tensor(
  x,
  progress,
  method = c("bin", "spline"),
  n_progress_bins = 10L,
  spline_df = 4L,
  keep_fitted = FALSE
)
```

## Arguments

- x:

  Array, cells x IZs x positions (e.g.
  [`build_replication_tensor()`](https://midimagi.github.io/scRTtools/reference/build_replication_tensor.md)).

- progress:

  Replication progress of the cells.

- method:

  `"bin"` or `"spline"`.

- n_progress_bins:

  Number of progress bins.

- spline_df:

  Spline degrees of freedom.

- keep_fitted:

  Logical; also return the fitted tensor.

## Value

A list with `residual`, `fitted` (or `NULL`) and `method`.
