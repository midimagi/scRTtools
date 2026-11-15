# Per-sample cell IDs and input files

Helpers for merging several samples.

## Usage

``` r
add_sample_prefix(cell_ids, sample_name)

find_sample_file(data_dir, sample_name, filename)
```

## Arguments

- cell_ids:

  Character vector of cell IDs (coerced with
  [`as.character()`](https://rdrr.io/r/base/character.html)). `NULL` is
  an error.

- sample_name:

  Sample name (a single string).

- data_dir:

  Directory to search.

- filename:

  File name to look for (exact match, case-sensitive).

## Value

`add_sample_prefix()`: a character vector as long as `cell_ids`.
`find_sample_file()`: the normalised absolute path of the single match.

## Details

`add_sample_prefix()` prepends `"<sample_name>-"` to cell IDs that do
not already start with it, so merged samples get unique cell IDs.

`find_sample_file()` searches `data_dir` recursively for a file called
`filename` whose immediate parent directory is exactly `sample_name`.
Files in other directories (e.g. a `merge/` output folder) and partial
name matches are ignored. It is an error if there is no match or more
than one.

## Examples

``` r
add_sample_prefix(c("AAAC-1", "s1-AAAG-1"), "s1")
#> [1] "s1-AAAC-1" "s1-AAAG-1"

d <- tempfile()
dir.create(file.path(d, "s1"), recursive = TRUE)
dir.create(file.path(d, "merge"))
invisible(file.create(file.path(d, "s1", "meta_data.rds"),
                      file.path(d, "merge", "meta_data.rds")))
find_sample_file(d, "s1", "meta_data.rds")
#> [1] "/tmp/RtmpEXCQa9/file2b5f147a221e/s1/meta_data.rds"
```
