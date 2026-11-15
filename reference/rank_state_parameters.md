# Rank evaluated replication-state parameters

Keep the successful rows of
[`evaluate_state_parameters()`](https://midimagi.github.io/scRTtools/reference/evaluate_state_parameters.md)
and order them by decreasing `target_prop_all_cells`,
`target_prop_non_missing` and `pass_qc_prop`, then by increasing
`distance_to_current` (sum of absolute differences to `reference`) and
finally by the parameter columns named in `reference` (increasing). `NA`
values sort last.

## Usage

``` r
rank_state_parameters(results, reference)
```

## Arguments

- results:

  Output of
  [`evaluate_state_parameters()`](https://midimagi.github.io/scRTtools/reference/evaluate_state_parameters.md).

- reference:

  Named numeric vector (or list) of the current parameter values. We use
  `c(init_ratio_lower = 1.2, final_ratio_lower = 1.4)`.

## Value

`results` restricted to `status == "ok"`, with a `distance_to_current`
column, in ranked order. An error is raised if no row succeeded.
