# Genome-bin identifier of a position

Genome-bin identifier of a position

## Usage

``` r
make_bin(chr, pos, bin_width = 20000)
```

## Arguments

- chr:

  Chromosome name(s).

- pos:

  1-based position(s).

- bin_width:

  Bin width in bp.

## Value

Character vector `"<chr>-<start>-<end>"` of the fixed-width bin holding
each position, always in plain (never scientific) notation.

## Examples

``` r
make_bin("chr1", c(1, 20000, 20001, 99999))
#> [1] "chr1-1-20000"      "chr1-1-20000"      "chr1-20001-40000" 
#> [4] "chr1-80001-100000"
```
