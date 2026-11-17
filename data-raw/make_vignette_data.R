# Build the larger example dataset used by vignettes/figures.Rmd.
#
# Source: sample scRT1 of the scRT results directory (mm10, 20 kb bins),
# genome-wide (chr1-chr19) for 500 cells: 100 cells called G1 and 400 cells
# called non-G1 by the original analysis. Unlike data-raw/make_example_data.R,
# cells are NOT restricted to pass_qc, so the vignette can show how many cells
# fail the replication-state QC. Within each group the cells are picked evenly
# along their coverage (no random numbers), so the rep_percent distribution and
# the QC failure rate stay representative of the sample. All cells have RNA.
#
# The output (about 49 MB) is NOT part of the package or of the git
# repository: upload it as an asset of a GitHub release and update
# `.vignette_data` (url, md5) in R/example_data.R; scrt_example_data("vignette")
# downloads it from there.
#
# Run from the package root:
#   conda run -n new_sc_envs Rscript data-raw/make_vignette_data.R
#
# Environment variables (required):
#   SCRTTOOLS_RESULTS  the scRT results directory
#   SCRTTOOLS_GTF      genes.gtf of refdata-cellranger-arc-mm10-2020-A-2.0.0
#   SCRTTOOLS_BULK_RT  bulk Repli-seq RT bigWig (Loess smoothed, mm10)
# Optional:
#   SCRTTOOLS_VIGNETTE_OUT  output file (default ../scRTtools_release/scRT1_vignette.rds)

t_start <- proc.time()[["elapsed"]]
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Matrix)
})

results_dir <- Sys.getenv("SCRTTOOLS_RESULTS")
gtf_file <- Sys.getenv("SCRTTOOLS_GTF")
bulk_rt_file <- Sys.getenv("SCRTTOOLS_BULK_RT")
if (!nzchar(results_dir) || !nzchar(gtf_file) || !nzchar(bulk_rt_file)) {
  stop("set SCRTTOOLS_RESULTS, SCRTTOOLS_GTF and SCRTTOOLS_BULK_RT")
}
sample_name <- "scRT1"
n_g1 <- 100L
n_s <- 400L
out_file <- Sys.getenv("SCRTTOOLS_VIGNETTE_OUT",
                       file.path("..", "scRTtools_release", "scRT1_vignette.rds"))

sample_dir <- file.path(results_dir, "1.data", sample_name)
merge_dir <- file.path(results_dir, "1.data", "merge")

# ---- cells -----------------------------------------------------------------
meta <- as.data.frame(readRDS(file.path(merge_dir, "meta_data.rds")))
rownames(meta) <- meta$cell_id
rna <- readRDS(file.path(merge_dir, "merge_rna_counts_filtered.csv"))  # an RDS despite its name
message("rna: ", class(rna)[1], " ", paste(dim(rna), collapse = " x "))

usable <- meta$sample == sample_name & !(meta$has_cnv %in% TRUE) &
  !is.na(meta$is_G1) & meta$cell_id %in% colnames(rna)

pick_evenly <- function(ids, n) ids[unique(round(seq(1, length(ids), length.out = n)))]
pick_by_coverage <- function(pool, n) {
  pick_evenly(pool$cell_id[order(pool$coverage, pool$cell_id)], n)
}
g1_cells <- pick_by_coverage(meta[usable & meta$is_G1, ], n_g1)
s_cells <- pick_by_coverage(meta[usable & !meta$is_G1, ], n_s)
cells <- c(g1_cells, s_cells)
stopifnot(length(cells) == n_g1 + n_s, !anyDuplicated(cells))

# ---- genome bins and raw counts ---------------------------------------------
bins <- readRDS(file.path(results_dir, "0.mm10", "genobins_20kb_granges.rds"))
chromosomes <- seqlevels(bins)
message("bins: ", length(bins), " on ", paste(chromosomes, collapse = ","))

subset_matrix <- function(file, row_ids, col_ids) {
  x <- readRDS(file)
  message(basename(file), ": ", class(x)[1], " ", paste(dim(x), collapse = " x "))
  # per-sample files name cells by bare barcode, merged files by "<sample>-<barcode>"
  lookup <- col_ids
  if (!all(lookup %in% colnames(x))) lookup <- sub(paste0("^", sample_name, "-"), "", col_ids)
  if (!all(row_ids %in% rownames(x)) || !all(lookup %in% colnames(x))) {
    stop(basename(file), ": requested rows or cells are missing")
  }
  x <- x[row_ids, lookup, drop = FALSE]
  colnames(x) <- col_ids
  if (inherits(x, "Matrix") || inherits(x, "IterableMatrix")) x <- as.matrix(x)
  storage.mode(x) <- "integer"
  x
}

bin_counts <- subset_matrix(file.path(sample_dir, "genobins_20kb_raw_matrix.rds"),
                            names(bins), cells)
gc()
g1_windows <- readRDS(file.path(sample_dir, "g1bins_granges.rds"))
stopifnot(all(seqlevels(g1_windows) %in% chromosomes))
window_counts <- subset_matrix(file.path(sample_dir, "g1bins_raw_matrix.rds"),
                               names(g1_windows), cells)
gc()

# ---- RNA counts of the genes on the bin chromosomes -------------------------
gtf <- data.table::fread(
  cmd = paste("awk -F '\t'", shQuote('$3 == "gene"'), shQuote(gtf_file)),
  header = FALSE, sep = "\t", quote = "", showProgress = FALSE
)
genes <- data.frame(
  gene = sub('.*gene_name "([^"]+)".*', "\\1", gtf$V9),
  chr = gtf$V1, start = as.integer(gtf$V4), end = as.integer(gtf$V5), strand = gtf$V7,
  stringsAsFactors = FALSE
)
genes <- genes[genes$chr %in% chromosomes, ]
genes <- genes[!duplicated(genes$gene) & genes$gene %in% rownames(rna), ]
genes <- genes[order(match(genes$chr, chromosomes), genes$start, genes$gene), ]
rownames(genes) <- NULL
rna_counts <- Matrix::Matrix(as.matrix(rna[genes$gene, cells, drop = FALSE]), sparse = TRUE)
rna_counts <- methods::as(rna_counts, "generalMatrix")
rm(rna); gc()

# ---- bulk Repli-seq RT on the 20 kb bins -------------------------------------
bulk <- rtracklayer::import(bulk_rt_file)
bulk <- keepSeqlevels(bulk, intersect(seqlevels(bulk), chromosomes), pruning.mode = "coarse")
hits <- findOverlaps(bins, bulk)
bulk_rt <- granges(bins)
mcols(bulk_rt) <- NULL
bulk_rt$score <- NA_real_
agg <- tapply(score(bulk)[subjectHits(hits)], queryHits(hits), mean)
bulk_rt$score[as.integer(names(agg))] <- as.numeric(agg)

# ---- metadata ---------------------------------------------------------------
keep_cols <- c("cell_id", "sample", "nCount_rna", "nFeature_rna", "coverage", "mapd",
               "scaled_mapd", "is_G1", "rep_percent", "sigma_chr", "has_cnv", "pass_qc")
cell_meta <- meta[cells, keep_cols]
cell_meta$has_rna <- TRUE

example <- list(
  bins = bins,
  bin_counts = bin_counts,
  g1_windows = g1_windows,
  window_counts = window_counts,
  rna_counts = rna_counts,
  genes = genes,
  bulk_rt = bulk_rt,
  meta = cell_meta,
  info = list(sample = sample_name, chromosomes = chromosomes, genome = "mm10",
              bin_width = 20000L, bulk_rt = basename(bulk_rt_file),
              created = as.character(Sys.Date()))
)
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(example, out_file, compress = "xz")

message("\nWrote ", out_file, " (", round(file.size(out_file) / 1024^2, 1), " MB, md5 ",
        unname(tools::md5sum(out_file)), ")")
f <- tempfile()
for (nm in names(example)) {
  saveRDS(example[[nm]], f, compress = "xz")
  x <- example[[nm]]
  message(sprintf("  %-14s %-12s %-16s %8.0f kB on disk", nm, class(x)[1],
                  paste(if (is.null(dim(x))) length(x) else dim(x), collapse = " x "),
                  file.size(f) / 1024))
}
print(table(G1 = cell_meta$is_G1, pass_qc = cell_meta$pass_qc, useNA = "ifany"))
print(table(G1 = cell_meta$is_G1, rep_percent_NA = is.na(cell_meta$rep_percent)))
print(summary(cell_meta$rep_percent[!cell_meta$is_G1]))
print(summary(bulk_rt$score))
message(sprintf("elapsed: %.1f s", proc.time()[["elapsed"]] - t_start))
