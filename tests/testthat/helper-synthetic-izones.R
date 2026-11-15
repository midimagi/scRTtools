# Replication states generated from fixed origins: origin o fires at time
# t_o; in a cell at replication progress p, bins within (p - t_o) * speed bins
# of a fired origin are replicated (state 4), others 2; a few NA bins.
synth_izone_data <- function(n_chr = 2L, n_bin = 500L, n_cells = 80L, seed = 41L) {
  withr::with_seed(seed, {
    chrs <- paste0("chr", seq_len(n_chr))
    gr <- GenomicRanges::GRanges(
      rep(chrs, each = n_bin),
      IRanges::IRanges(rep(seq(1L, by = 20000L, length.out = n_bin), n_chr), width = 20000L)
    )
    names(gr) <- paste(GenomicRanges::seqnames(gr), BiocGenerics::start(gr),
                       BiocGenerics::end(gr), sep = "-")
    gr$is_mappable <- TRUE
    origins <- do.call(rbind, lapply(seq_along(chrs), function(k) {
      pos <- sort(sample(seq(15, n_bin - 15, by = 30), 12))
      data.frame(chr = chrs[k], bin = pos + (k - 1L) * n_bin, t = stats::runif(12, 0, 0.7))
    }))
    p <- sort(stats::runif(n_cells, 0.05, 0.95))
    st <- vapply(p, function(pj) {
      s <- rep(2, length(gr))
      for (o in seq_len(nrow(origins))) {
        if (origins$t[o] < pj) {
          w <- round((pj - origins$t[o]) * 40 + stats::runif(1, -1, 1))
          chr_rng <- range(which(as.character(GenomicRanges::seqnames(gr)) == origins$chr[o]))
          b <- origins$bin[o] + sample(-1:1, 1)
          s[max(chr_rng[1], b - w):min(chr_rng[2], b + w)] <- 4
        }
      }
      s[stats::runif(length(s)) < 0.01] <- NA
      s
    }, numeric(length(gr)))
  })
  cells <- sprintf("cell%03d", seq_len(n_cells))
  dimnames(st) <- list(names(gr), cells)
  rp <- stats::setNames(colSums(st == 4, na.rm = TRUE) / colSums(!is.na(st)), cells)
  tracks <- define_replication_tracks(st, gr, rp, verbose = FALSE)
  peaks <- gr[origins$bin[-c(3, 15)]]              # two origins without a seed peak
  S4Vectors::mcols(peaks) <- NULL
  list(gr = gr, st = st, rp = rp, tracks = tracks, peaks = peaks, origins = origins,
       md = data.frame(cell_id = cells, rep_percent = unname(rp)))
}

source_orig_izones <- function() {
  env <- source_orig_file("processing/helper_functions.R")
  env$seqlevels <- GenomeInfoDb::seqlevels
  env$relist <- IRanges::relist
  env$unlist <- BiocGenerics::unlist
  for (f in c("processing/define_initiation_zones_fixed.R", "analysis/calculate_firing_range.R",
              "analysis/calculate_firing_order.R", "analysis/predict_firing_order.R",
              "analysis/izone_velocity_model.R")) {
    path <- file.path(orig_rcode_dir(), f)
    testthat::skip_if_not(file.exists(path), paste("original", f, "not available"))
    sys.source(path, envir = env)
  }
  env$s <- mgcv::s
  env
}

# IZs with firing range and an rt column (fraction of cells replicated at the
# IZ centre bin).
izone_fixture <- function(x = synth_izone_data()) {
  diz <- define_initiation_zones(x$tracks, x$peaks, min_usage = 5L, verbose = FALSE)
  iz <- diz$izones
  ctr <- GenomicRanges::resize(iz, 1L, fix = "center")
  bin <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(ctr, x$gr))
  iz$rt <- rowMeans(x$st[bin, , drop = FALSE] == 4, na.rm = TRUE)
  list(x = x, diz = diz, iz = iz)
}

# Function definitions taken verbatim from a notebook's code cells.
notebook_dir <- function() {
  local_dir("SCRTTOOLS_NOTEBOOKS")
}

notebook_functions <- function(notebook, fun_names, parent = globalenv()) {
  testthat::skip_if_not_installed("jsonlite")
  path <- file.path(notebook_dir(), notebook)
  testthat::skip_if_not(file.exists(path), paste("notebook", notebook, "not available"))
  nb <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  code <- unlist(lapply(nb$cells, function(c) {
    if (identical(c$cell_type, "code")) paste(unlist(c$source), collapse = "") else NULL
  }))
  env <- new.env(parent = parent)
  for (chunk in code) {
    exprs <- tryCatch(parse(text = chunk, keep.source = FALSE), error = function(e) NULL)
    for (e in exprs) {
      if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") && is.name(e[[2]]) &&
          as.character(e[[2]]) %in% fun_names) {
        eval(e, envir = env)
      }
    }
  }
  missing <- setdiff(fun_names, ls(env))
  if (length(missing)) stop("Not found in notebook: ", paste(missing, collapse = ", "))
  env
}

source_orig_heterogeneity <- function() {
  testthat::skip_if_not_installed("dplyr")
  env <- new.env(parent = asNamespace("dplyr"))
  for (f in c("iz_heterogeneity_phase1_3.R", "iz_heterogeneity_state_usage.R")) {
    path <- file.path(notebook_dir(), f)
    testthat::skip_if_not(file.exists(path), paste(f, "not available"))
    sys.source(path, envir = env, toplevel.env = env)
  }
  env
}

# Assigned tracks with firing ranges (inputs of the adjacency / heterogeneity
# functions).
iz_analysis_fixture <- function() {
  f <- izone_fixture()
  fr <- calculate_firing_range(f$iz, f$x$st, f$diz$tracks, verbose = FALSE)
  c(f, list(fr = fr))
}
