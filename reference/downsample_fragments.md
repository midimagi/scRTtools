# Randomly downsample a 10x fragment file

Stream a coordinate-sorted, bgzipped fragment file once and keep each
fragment record independently with probability `keep_prob` (SplitMix64
RNG seeded with `seed`; header lines starting with `#` are always kept),
then write BGZF output and build a tabix BED index. This is an R wrapper
around the packaged shell script
`system.file("scripts", "downsample_fragments_fast.sh", package = "scRTtools")`,
which compiles a small C sampler with `cc` at run time.

## Usage

``` r
downsample_fragments(
  input,
  output,
  keep_prob,
  seed = 20260821,
  threads = 16L,
  compression_level = 1L,
  overwrite = FALSE,
  bgzip = Sys.which("bgzip"),
  tabix = Sys.which("tabix"),
  cc = Sys.which("cc"),
  verbose = FALSE
)
```

## Arguments

- input:

  Path to the input `.tsv.gz` fragment file.

- output:

  Path to the output `.tsv.gz` file (must differ from `input`; its
  directory must exist and be writable).

- keep_prob:

  Keep probability in `(0, 1]`, e.g.
  `optimize_shrink_factor()$shrink_factor`.

- seed:

  Non-negative integer seed of the sampler.

- threads:

  Positive integer; bgzip threads.

- compression_level:

  Integer 0-9; bgzip compression level.

- overwrite:

  Logical; replace an existing output and index.

- bgzip, tabix, cc:

  Paths to the `bgzip`, `tabix` and C compiler executables. Default:
  found on `PATH`.

- verbose:

  Logical; show the script's progress output. Default `FALSE`.

## Value

Invisibly, `output`.

## Details

Besides `output`, the script writes `<output>.tbi` and
`<output>.downsample_info.txt` (parameters and kept/total record
counts). Results are reproducible for a given `seed` and input.
