# Split a signed residual tensor into non-negative parts

Split a signed residual tensor into non-negative parts

## Usage

``` r
split_signed_residual(residual)
```

## Arguments

- residual:

  Array, cells x IZs x positions (e.g.
  `correct_replication_tensor()$residual`).

## Value

An array with the third dimension doubled: positive parts followed by
negative parts.
