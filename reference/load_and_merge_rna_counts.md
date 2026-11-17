# Merge RNA count matrices of several samples

Takes one Seurat object (or gene-by-cell count matrix) per sample,
extracts the raw counts of `assay`, prefixes cell IDs with `"<sample>-"`
(unless already prefixed), keeps the genes common to all samples and
column-binds the matrices. To start from files, read them first, e.g.
`lapply(stats::setNames(files, samples), readRDS)`.

## Usage

``` r
load_and_merge_rna_counts(rna_objects, assay = "rna")
```

## Arguments

- rna_objects:

  Named list of Seurat objects or gene-by-cell count matrices, one per
  sample; the names are the sample names.

- assay:

  Name of the RNA assay (used for Seurat objects).

## Value

A gene-by-cell count matrix (sparse, as stored in the objects).
