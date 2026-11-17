# scRTtools: main figures from a 500-cell example

> **The code of this vignette was not run when this copy was built.** It
> takes 20-35 minutes (mostly
> [`assign_replication_states()`](https://midimagi.github.io/scRTtools/reference/assign_replication_states.md)
> on 500 cells) and needs ggplot2, hexbin, robustbase, irlba and uwot.
> Build it with `SCRTTOOLS_BUILD_FIGURES=true R CMD build scRTtools`, or
> knit this file after setting that environment variable.

This vignette recomputes the main figures of a scRT analysis with the
package functions, from raw counts to replication tracks, initiation
zones and gene expression. It uses the larger example dataset:
`scrt_example_data("vignette")` holds genome-wide 20 kb counts of 500
cells of the mouse ESC multiome sample `scRT1` (100 cells called G1 and
400 called non-G1 in the original analysis, not filtered on QC), their
RNA counts and a bulk Repli-seq replication-timing track. These data (49
MB) are not installed with the package; the first call of
`scrt_example_data("vignette")` downloads them from a GitHub release and
caches them. Because it has about a tenth of the cells of the full
sample, the curves are noisier and the counts (cells, initiation zones)
are smaller than in the full analysis, but the shapes are the same.

The package has no plotting functions; the figures are drawn with
ggplot2 (plus hexbin, irlba and uwot for single plots), which the
vignette needs in addition to the package. Each step section names the
package functions it uses.

| Section | Figure | Main functions |
|----|----|----|
| 1 | MAPD vs coverage, G1 calls | [`calculate_mapd()`](https://midimagi.github.io/scRTtools/reference/calculate_mapd.md), [`find_g1_cells()`](https://midimagi.github.io/scRTtools/reference/find_g1_cells.md) |
| 2 | QC tables | [`assign_replication_states()`](https://midimagi.github.io/scRTtools/reference/assign_replication_states.md), [`compute_rep_percent()`](https://midimagi.github.io/scRTtools/reference/compute_rep_percent.md) |
| 3 | rep_percent of S cells | [`compute_rep_percent()`](https://midimagi.github.io/scRTtools/reference/compute_rep_percent.md) |
| 4 | Pseudobulk vs bulk RT | [`aggregate_pseudobulk_s_g1_counts()`](https://midimagi.github.io/scRTtools/reference/aggregate_pseudobulk_s_g1_counts.md), [`calculate_pseudobulk_rt_g1_ratio()`](https://midimagi.github.io/scRTtools/reference/calculate_pseudobulk_rt_g1_ratio.md) |
| 5 | Single-cell replication states | [`compute_rep_percent()`](https://midimagi.github.io/scRTtools/reference/compute_rep_percent.md) |
| 6 | UMAP of replication states | (irlba, uwot) |
| 7 | Track length and its velocity | [`define_replication_tracks()`](https://midimagi.github.io/scRTtools/reference/define_replication_tracks.md), [`summarize_track_lengths()`](https://midimagi.github.io/scRTtools/reference/summarize_track_lengths.md), [`fit_track_velocity_linear()`](https://midimagi.github.io/scRTtools/reference/fit_track_velocity_linear.md) |
| 8 | Initiation zones: width, distance, RT | [`count_track_centers()`](https://midimagi.github.io/scRTtools/reference/count_track_centers.md), [`define_initiation_zones()`](https://midimagi.github.io/scRTtools/reference/define_initiation_zones.md) |
| 9 | Gene replication vs expression | [`rep_tempo()`](https://midimagi.github.io/scRTtools/reference/rep_tempo.md), [`fit_rep_tempo()`](https://midimagi.github.io/scRTtools/reference/fit_rep_tempo.md), [`gex_tempo()`](https://midimagi.github.io/scRTtools/reference/gex_tempo.md) |

``` r

library(scRTtools)
library(GenomicRanges)
library(ggplot2)

ex <- scrt_example_data("vignette")
names(ex)
table(called_G1 = ex$meta$is_G1)
```

All figures share one ggplot2 theme (grey panel, white grid):

``` r

theme_scrt <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(panel.background = element_rect(fill = "grey90", colour = NA),
          panel.grid.major = element_line(colour = "white", linewidth = 0.5),
          panel.grid.minor = element_line(colour = "white", linewidth = 0.25),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 0.5),
          axis.text = element_text(colour = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.key = element_blank())
}
col_g1 <- "orangered"
col_s <- "dodgerblue"
```

## 1. Cell-cycle calls from MAPD and coverage

Replicating (S-phase) cells have a mix of one- and two-copy regions, so
their read counts vary more between neighbouring bins than those of G1
cells. The MAPD (median absolute pairwise difference of log counts of
adjacent bins) captures this variation; it also falls with sequencing
depth, so
[`find_g1_cells()`](https://midimagi.github.io/scRTtools/reference/find_g1_cells.md)
scales it by coverage and calls the low-MAPD cells G1.

``` r

bin_counts <- Matrix::Matrix(ex$bin_counts, sparse = TRUE)
mapd <- calculate_mapd(bin_counts, ex$bins, ex$meta[, c("cell_id", "sample")])
g1 <- find_g1_cells(mapd, seed = 1)
table(G1 = g1$is_G1)
```

``` r

ggplot(g1, aes(coverage, scaled_mapd, colour = is_G1)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_colour_manual(values = c(`TRUE` = col_g1, `FALSE` = col_s), name = "is_G1") +
  labs(x = "Coverage", y = "Scaled MAPD") +
  theme_scrt()
```

## 2. Replication states and QC

Replication states are inferred for every cell, G1 cells included: on
variable-width windows that hold about the same number of G1 reads, a
hidden Markov model labels each window as unreplicated (2 copies) or
replicated (4 copies).
[`compute_rep_percent()`](https://midimagi.github.io/scRTtools/reference/compute_rep_percent.md)
then summarizes each cell (fraction replicated, `rep_percent`) and flags
cells whose chromosomes disagree too much (`pass_qc`). This is the
slowest step of the vignette.

``` r

window_counts <- Matrix::Matrix(ex$window_counts, sparse = TRUE)
res_states <- assign_replication_states(
  window_counts, ex$g1_windows, ex$bins,
  BPPARAM = BiocParallel::MulticoreParam(n_workers),
  chunk_size = ceiling(ncol(window_counts) / n_workers)  # one block of cells per worker
)
res_rp <- compute_rep_percent(res_states$genobins_rep_states, ex$bins,
                              ex$meta[, c("cell_id", "sample")])
md <- res_rp$meta_data
md$is_G1 <- g1[rownames(md), "is_G1"]
md$pass_qc <- md$pass_qc %in% TRUE
```

G1 cells should fail the replication-state QC or get no `rep_percent` at
all, and most S cells should pass:

``` r

plot_table <- function(x, y, xlab, ylab) {
  tab <- as.data.frame(table(x = x, y = y))
  ggplot(tab, aes(x, y, fill = Freq)) +
    geom_tile(colour = "grey70") +
    geom_text(aes(label = Freq), size = 5) +
    scale_fill_gradient(low = "lightyellow", high = "orangered", name = "Count") +
    labs(x = xlab, y = ylab) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank())
}
p1 <- plot_table(md$pass_qc, md$is_G1, "pass_qc", "is_G1")
p2 <- plot_table(md$is_G1, is.na(md$rep_percent), "is_G1", "is.na(rep_percent)")
gridExtra_free <- function(...) {   # side by side without extra packages
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  for (i in seq_along(list(...))) {
    print(list(...)[[i]], vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i))
  }
}
gridExtra_free(p1, p2)
```

From here on, “S cells” are the cells called non-G1 that pass QC.

``` r

s_cells <- rownames(md)[!md$is_G1 & md$pass_qc & is.finite(md$rep_percent)]
length(s_cells)
states <- res_rp$genobins_rep_states_filtered[, s_cells, drop = FALSE]
md_s <- md[s_cells, , drop = FALSE]
```

## 3. Replication progress of S cells

``` r

ggplot(md_s, aes(rep_percent)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, boundary = 0,
                 fill = "steelblue", colour = "white", alpha = 0.8) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "Replication percentage", y = "Density",
       title = "Replication percentage per cell") +
  theme_scrt()
```

## 4. Pseudobulk replication timing vs bulk Repli-seq

Pooling the reads of all S cells and dividing by the pooled G1 reads
gives a pseudobulk replication-timing profile: early-replicating regions
have more copies in S phase. It should agree with bulk Repli-seq of the
same cell type. Both tracks are averaged in 50 kb windows before
comparing them.

``` r

rt_meta <- md[, c("cell_id", "sample", "is_G1", "rep_percent", "pass_qc")]
rt_counts <- aggregate_pseudobulk_s_g1_counts(bin_counts, rt_meta, names(ex$bins))
attr(rt_counts, "selected_cells")
rt <- calculate_pseudobulk_rt_g1_ratio(rt_counts, ex$bins)
```

``` r

windows <- tileGenome(seqinfo(ex$bins), tilewidth = 50000L,
                      cut.last.tile.in.chrom = TRUE)
window_mean <- function(gr) {
  gr <- gr[is.finite(gr$score)]
  hits <- findOverlaps(windows, gr)
  out <- rep(NA_real_, length(windows))
  agg <- tapply(gr$score[subjectHits(hits)], queryHits(hits), mean)
  out[as.integer(names(agg))] <- agg
  out
}
df_rt <- data.frame(bulk = window_mean(ex$bulk_rt), pseudo = window_mean(rt))
df_rt <- df_rt[stats::complete.cases(df_rt), ]
r <- cor(df_rt$bulk, df_rt$pseudo)
ggplot(df_rt, aes(bulk, pseudo)) +
  geom_hex(bins = 100) +
  scale_fill_viridis_c(trans = "log10", name = "Windows") +
  coord_cartesian(ylim = c(-2.5, 2.5)) +   # a few outlier windows lie beyond
  labs(title = "Bulk RT vs pseudobulk RT",
       subtitle = sprintf("S/G1 ratio, 50 kb windows, Pearson r = %.3f", r),
       x = "Bulk RT (Loess smoothed)", y = "Pseudobulk RT (z-score)") +
  theme_bw(base_size = 13)
```

## 5. Single-cell replication states

Each row is one S cell, sorted by `rep_percent`; replicated regions (4
copies, blue) widen from early- to late-replicating cells.

``` r

in_region <- which(as.character(seqnames(ex$bins)) == "chr4" &
                   start(ex$bins) >= 1e6 & end(ex$bins) <= 50e6)
cell_order <- order(md_s$rep_percent)
df_states <- data.frame(
  pos = rep(start(ex$bins)[in_region] / 1e6, times = length(cell_order)),
  cell = rep(seq_along(cell_order), each = length(in_region)),
  state = factor(as.vector(states[in_region, cell_order]), levels = c(2, 4))
)
y_breaks <- unique(round(seq(1, length(cell_order), length.out = 6)))
ggplot(df_states, aes(pos, cell, fill = state)) +
  geom_raster() +
  scale_fill_manual(values = c(`2` = "lightyellow", `4` = "royalblue"),
                    na.value = "grey90", name = "State") +
  scale_y_reverse(breaks = y_breaks, expand = c(0, 0),
                  labels = paste0(round(100 * md_s$rep_percent[cell_order][y_breaks]), "%")) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "chr4 coordinate (Mb)", y = "% Replicated", title = "chr4:1-50 Mb") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
```

## 6. UMAP of replication states

The replication states alone order the S cells along S phase. As in the
original analysis, the binary state matrix is normalized like an ATAC
matrix (TF-IDF, log), reduced with a truncated SVD and embedded with
UMAP; bins without a call count as unreplicated.

``` r

x <- states
x[is.na(x)] <- 2
x <- Matrix::Matrix(x == 4, sparse = TRUE) * 1
x <- x[Matrix::rowSums(x) > 0, , drop = FALSE]
tfidf <- Matrix::t(Matrix::t(x) / Matrix::colSums(x)) / Matrix::rowMeans(x)
tfidf <- log1p(1e4 * tfidf)
set.seed(2026)
svd <- irlba::irlba(scale(Matrix::t(tfidf), center = TRUE, scale = FALSE), nv = 30)
pcs <- svd$u %*% diag(svd$d)
umap <- uwot::umap(pcs, n_neighbors = 15, seed = 2026)
df_umap <- data.frame(UMAP_1 = umap[, 1], UMAP_2 = umap[, 2],
                      rep_percent = md_s$rep_percent,
                      rna_expr = ex$meta[s_cells, "nCount_rna"] /
                        mean(ex$meta[s_cells, "nCount_rna"]))
```

``` r

p_rep <- ggplot(df_umap, aes(UMAP_1, UMAP_2, colour = rep_percent)) +
  geom_point(size = 1) +
  scale_colour_gradient(low = "grey85", high = "midnightblue") +
  labs(title = "Coloured by rep_percent") + theme_scrt(11)
p_rna <- ggplot(df_umap, aes(UMAP_1, UMAP_2, colour = rna_expr)) +
  geom_point(size = 1) +
  scale_colour_gradient(low = "grey90", high = "firebrick", limits = c(0, 2),
                        oob = scales::squish) +
  labs(title = "Coloured by relative RNA depth") + theme_scrt(11)
gridExtra_free(p_rep, p_rna)
```

The RNA depth of a cell (its RNA counts over the mean of the S cells)
rises along the same trajectory: cells late in S phase have more RNA.

## 7. Replication track length and its velocity

A replication track is a run of replicated bins in one cell. Tracks grow
as forks move outwards from initiation zones, so the median track length
of a cell increases with its replication progress, slowly in early S and
fast when neighbouring tracks merge in late S.
[`fit_track_velocity_linear()`](https://midimagi.github.io/scRTtools/reference/fit_track_velocity_linear.md)
fits `log(length) ~ rep_percent + rep_percent^2` robustly and derives
the growth rate (“velocity”) of the track length per percentage point.

``` r

tracks <- define_replication_tracks(states, ex$bins,
                                    setNames(md_s$rep_percent, s_cells))
track_lengths <- calc_rep_track_lengths(tracks)
median_lengths <- summarize_track_lengths(track_lengths)
head(median_lengths, 3)
velocity <- fit_track_velocity_linear(median_lengths, n_grid = 200, fit_seed = 1)
velocity$validation_metrics
```

``` r

p_len <- ggplot(median_lengths, aes(rep_percent, median_length / 1e6)) +
  geom_point(colour = col_s, alpha = 0.6, size = 1) +
  geom_line(data = velocity$prediction, aes(rep_percent, lmrob_length / 1e6),
            colour = col_g1, linewidth = 1) +
  labs(x = "Replication percentage (%)", y = "Median length of tracks (Mb)",
       title = "Median track length per cell") + theme_scrt(11)
p_vel <- ggplot(velocity$prediction, aes(rep_percent, length_velocity / 1e3)) +
  geom_line(colour = col_s, linewidth = 1) +
  labs(x = "Replication percentage (%)",
       y = "Track-length velocity (kb / percentage point)",
       title = "Growth rate of the track length") + theme_scrt(11)
gridExtra_free(p_len, p_vel)
```

## 8. Initiation zones

Short tracks are centred on the initiation zones (IZs) that fired them.
The centres of short, gap-free tracks are counted per bin; bins hit by
several cells seed the IZs, and
[`define_initiation_zones()`](https://midimagi.github.io/scRTtools/reference/define_initiation_zones.md)
then clusters the tracks of all cells around these seeds. An IZ must be
supported by tracks of at least `min_usage` cells.

``` r

centers <- count_track_centers(tracks, ex$bins,
                               track_filter = (num_bins <= 50) & (flag == 0))
table(pmin(centers$score, 5))
seeds <- centers[centers$score >= 3]
res_iz <- define_initiation_zones(
  tracks, seeds, max_track_length = "1Mb", min_usage = 3L, bin_size = "20Kb",
  BPPARAM = BiocParallel::MulticoreParam(n_workers)
)
iz <- res_iz$izones
length(iz)
```

For each IZ we report its width, the distance to the next IZ on the same
chromosome, and its replication timing (the mean pseudobulk RT of
section 4 over the IZ). With about 300 S cells only the more frequently
used IZs reach `min_usage`, so there are fewer IZs than in the full
sample and the distance between neighbouring IZs is larger; widths and
timing are comparable.

``` r

next_start <- c(start(iz)[-1], NA)
same_chr <- c(as.character(seqnames(iz))[-1] == as.character(seqnames(iz))[-length(iz)], FALSE)
iz_distance <- (next_start - end(iz) - 1)[same_chr]
hits <- findOverlaps(iz, rt[is.finite(rt$score)])
iz$rt <- NA_real_
agg <- tapply(rt[is.finite(rt$score)]$score[subjectHits(hits)], queryHits(hits), mean)
iz$rt[as.integer(names(agg))] <- agg
summary(width(iz))
```

``` r

iz_hist <- function(df, x, label, xlab, log_x = TRUE, density = FALSE) {
  med <- stats::median(df[[x]], na.rm = TRUE)
  p <- ggplot(df, aes(.data[[x]])) +
    geom_histogram(aes(y = after_stat(density)), bins = 40, boundary = 0,
                   fill = "steelblue", colour = "white", alpha = 0.8) +
    geom_vline(xintercept = med, linetype = "dashed", colour = "grey30") +
    labs(x = xlab, y = "Density",
         subtitle = paste0(nrow(df), " IZs, median ", label, " = ",
                           if (log_x) format_bp(med) else round(med, 2))) +
    theme_scrt(11) + theme(plot.margin = margin(5, 12, 5, 5))
  if (density) p <- p + geom_density(colour = col_g1, linewidth = 1)
  if (log_x) p <- p + scale_x_log10(labels = format_bp)
  p
}
p_w <- iz_hist(data.frame(width = width(iz)), "width", "width", "Width of IZs")
p_d <- iz_hist(data.frame(distance = iz_distance[iz_distance > 0]), "distance",
               "distance", "Distance between IZs")
p_rt <- iz_hist(data.frame(rt = iz$rt[is.finite(iz$rt)]), "rt", "RT",
                "Pseudobulk RT of IZs", log_x = FALSE, density = TRUE)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 3)))
for (i in 1:3) print(list(p_w, p_d, p_rt)[[i]],
                     vp = grid::viewport(layout.pos.row = 1, layout.pos.col = i))
```

## 9. Gene replication vs gene expression

For every gene, the S cells are split into bins of replication progress
and
[`rep_tempo()`](https://midimagi.github.io/scRTtools/reference/rep_tempo.md)
reports the fraction of cells of each bin in which the gene is
replicated.
[`fit_rep_tempo()`](https://midimagi.github.io/scRTtools/reference/fit_rep_tempo.md)
fits a logistic curve per gene; its `t50` (the progress at which half of
the cells have replicated the gene) orders the genes from early to late
replicating.

``` r

genes <- ex$genes
rna <- ex$rna_counts[, s_cells]
genes <- genes[genes$gene %in% rownames(rna)[Matrix::rowSums(rna > 0) >= 20], ]
gene_gr <- GRanges(genes$chr, IRanges(genes$start, genes$end), gene = genes$gene)
names(gene_gr) <- genes$gene
nbin <- 20L
gene_tempo <- rep_tempo(states, gene_gr, meta_data = md_s, nbin_rep = nbin)
gene_fit <- fit_rep_tempo(gene_tempo, model_type = "sslogis",
                          BPPARAM = BiocParallel::MulticoreParam(n_workers))
table(gene_fit$fit_status)
```

[`gex_tempo()`](https://midimagi.github.io/scRTtools/reference/gex_tempo.md)
bins the same cells in the same way and summarizes the expression of
each gene per bin; with `method = "sc2p"` it separates the fraction of
cells expressing the gene (phase I) from the mean expression in the
expressing cells (phase II).

``` r

gex <- suppressWarnings(gex_tempo(as.matrix(rna[genes$gene, ]), md_s,
                                  nbin_rep = nbin, method = "sc2p"))
names(gex)
```

Both heatmaps use the same gene order (by `t50`, earliest at the top)
and z-scores per gene. Replication (left) moves as a front through S
phase, whereas expression (right) has no such front.

``` r

rep_bins <- attr(gene_tempo, "rep_bins")   # progress bins, in % (0-100)
filled <- rep_bins$n_cells > 0
# genes whose logistic fit passed QC; "partial_thresholds" only misses t10/t90
ok <- gene_fit$fit_status %in% c("ok", "partial_thresholds") & is.finite(gene_fit$t50)
gene_order <- gene_fit$region_id[ok][order(gene_fit$t50[ok])]
gene_order <- gene_order[gene_order %in% rownames(gex$df_phase_II)]
heat_long <- function(mat, title) {
  mat <- as.matrix(mat[gene_order, filled, drop = FALSE])
  z <- t(apply(mat, 1, zscore_safe))
  data.frame(progress = rep((rep_bins$rep_min + rep_bins$rep_max)[filled] / 2,
                            each = nrow(z)),
             gene = rep(rev(seq_len(nrow(z))), times = ncol(z)),
             z = as.vector(z), panel = title)
}
df_heat <- rbind(heat_long(gene_tempo, "Gene replication"),
                 heat_long(gex$df_phase_II, "Phase II mean expression"))
bin_width <- (rep_bins$rep_max - rep_bins$rep_min)[1]
ggplot(df_heat, aes(progress, gene, fill = z)) +
  geom_tile(width = bin_width, height = 1) +
  facet_wrap(~panel) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-3, 3),
                       oob = scales::squish, na.value = "grey80", name = "Z-score") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Replication progress of cell (%)", y = "Genes (ordered by t50)") +
  theme_bw(base_size = 13) +
  theme(strip.text = element_text(face = "bold"))
```

## Session

``` r

sessionInfo()
```
