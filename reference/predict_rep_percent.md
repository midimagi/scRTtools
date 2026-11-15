# Replication percentage needed to reach given track lengths

Invert a fitted track-length model: for each length, the replication
percentage at which the predicted median track length equals it, and the
average length velocity between consecutive lengths,
`diff(track_length) / diff(rep_percent)`.

## Usage

``` r
predict_rep_percent(fit, track_length)
```

## Arguments

- fit:

  Output of
  [`fit_track_velocity_linear()`](https://midimagi.github.io/scRTtools/reference/fit_track_velocity_linear.md)
  or
  [`fit_track_velocity_ode()`](https://midimagi.github.io/scRTtools/reference/fit_track_velocity_ode.md).

- track_length:

  Positive numeric vector of track lengths (bp).

## Value

A `data.frame` with `track_length`, `rep_percent` and `length_velocity`
(`NA` in the first row).

## Details

For the linear model the right-hand root of the quadratic
`b2 x^2 + b1 x + b0 = log(length)` is returned (`NA` if there is none).
For the ODE model the root is searched with `stats::uniroot(tol = 1e-8)`
between the turning point (clamped to the data range) and the largest
observed `rep_percent` (`NA` if the length is out of range).
