# Tables of CP factors

Tables of CP factors

## Usage

``` r
cp_factor_tables(
  fit,
  cell_ids,
  iz_ids,
  relative_bp,
  cell_metadata,
  iz_summary,
  signed_track = TRUE
)
```

## Arguments

- fit:

  Output of
  [`fit_nonnegative_cp()`](https://midimagi.github.io/scRTtools/reference/fit_nonnegative_cp.md).

- cell_ids, iz_ids, relative_bp:

  Labels of the three tensor modes.

- cell_metadata:

  Table with `cell_id` and `replication_progress`.

- iz_summary:

  Per-IZ table with `IZ_id`.

- signed_track:

  Logical; the tensor came from
  [`split_signed_residual()`](https://midimagi.github.io/scRTtools/reference/split_signed_residual.md)
  (track factors = positive minus negative part).

## Value

A list with `cell`, `iz`, `track`, `progress` (Spearman correlation of
each cell factor with progress) and `top_iz` (top 5% IZs per factor).
