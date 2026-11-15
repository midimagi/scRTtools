# Track-length velocity of each initiation zone

For each IZ, fit a GAM `log1p(mid_track_len) ~ s(rep_percent)`
([`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html), GCV,
`k = max(4, min(gam_basis_k, n_unique - 1, n - 2))`) to its tracks
(finite values, `rep_percent` in `[0, 1]`, positive length) and compare
it with an intercept-only model. The fitted length curve on
`n_curve_points` grid points gives the velocity (bp per percent, central
finite differences). IZs with fewer than `min_points_per_izone` tracks
or fewer than 5 distinct `rep_percent` values are not fitted. GAM
p-values are BH-adjusted (non-finite p-values count as 1); a velocity
model passes if the fit is ok, `fdr < fdr_threshold` and the GAM has the
lower AIC.

## Usage

``` r
fit_all_izone_velocities(
  df,
  izone_ids,
  min_points_per_izone,
  gam_basis_k = 10L,
  n_curve_points = 101L,
  fdr_threshold = 0.05,
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- df:

  A `data.frame` with columns `iz_id`, `rep_percent` (fraction) and
  `mid_track_len` (bp), one row per track.

- izone_ids:

  IZ ids to fit (also those without rows).

- min_points_per_izone:

  Minimum number of tracks per IZ.

- gam_basis_k:

  Maximum basis dimension (e.g. `10`).

- n_curve_points:

  Number of curve points, at least 3 (e.g. `101`).

- fdr_threshold:

  FDR threshold (e.g. `0.05`).

- BPPARAM:

  A
  [`BiocParallel::BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  used over IZs.

## Value

A list with `summary` (one row per IZ: fit statistics, `fdr`,
`delta_aic_null_minus_gam`, `aic_pass`, `velocity_model_pass`,
`velocity_direction`, `fit_status`), `curves` (fitted length and
velocity per grid point) and `config` (the parameters used).
