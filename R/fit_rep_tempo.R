# Migrated from rcode/analysis/tempo_models1.R (the version sourced by
# notebooks 4.8 and 6.2). Changes: n_boot defaults to 0 (was 200), and the
# parallel backend argument is BPPARAM = SerialParam() (was bpparam = NULL).
# The dyn.load/.C OpenBLAS thread setter was dropped (R CMD check NOTE); worker
# threads are limited through environment variables and RhpcBLASctl.

# Limit nested numerical-library threads in child processes only.
# OpenBLAS exposes a Fortran ABI setter accepting an integer pointer, which
# can be called with .C(); the C ABI (without trailing "_") cannot.
.tm_worker_threads <- function() {
    Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
               MKL_NUM_THREADS = "1", BLIS_NUM_THREADS = "1")
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
        RhpcBLASctl::blas_set_num_threads(1L)
        RhpcBLASctl::omp_set_num_threads(1L)
    }
    invisible(NULL)
}

.tm_ids <- function(x, n, label) {
    if (is.null(x) || length(x) != n || anyNA(x) || any(!nzchar(x)) ||
        anyDuplicated(x)) stop(label, " must be unique, non-empty IDs.")
}
.tm_integer <- function(x, lower, label) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < lower || x != floor(x) || x >= .Machine$integer.max) {
        stop(label, " must be an integer >= ", lower, ".")
    }
}
.tm_try <- function(expr) {
    messages <- character()
    value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
        messages <<- c(messages, conditionMessage(w))
        invokeRestart("muffleWarning")
    }), error = function(e) {
        messages <<- c(messages, conditionMessage(e))
        NULL
    })
    list(value = value, message = paste(unique(messages), collapse = " | "))
}

# Same three-parameter NLS and fixed anchors as 4.6, including plinear fallback.
.tm_sslogis <- function(d, anchors) {
    augmented <- data.frame(time = (50 - d$x) / 5, percentage = d$y)
    if (anchors) augmented <- rbind(augmented,
        data.frame(time = c(-10, 10), percentage = c(1, 0)))
    messages <- character()
    for (algorithm in c("port", "plinear")) {
        attempt <- .tm_try({
            if (algorithm == "port") {
                stats::nls(percentage ~ SSlogis(time, Asym, xmid, scal),
                    data = augmented, start = c(Asym = 1, xmid = 0, scal = -0.5),
                    algorithm = "port", control = stats::nls.control(maxiter = 100L))
            } else {
                stats::nls(percentage ~ SSlogis(time, 1, xmid, scal),
                    data = augmented, start = c(xmid = 0, scal = -0.5),
                    algorithm = "plinear",
                    control = stats::nls.control(maxiter = 100L, tol = 1e-4,
                                                warnOnly = TRUE))
            }
        })
        messages <- c(messages, attempt$message)
        model <- attempt$value
        if (is.null(model)) next
        pars <- stats::coef(model)
        if (algorithm == "plinear") {
            pars <- c(Asym = unname(pars[".lin"]), pars[c("xmid", "scal")])
        }
        if (!isTRUE(model$convInfo$isConv) || any(!is.finite(pars)) ||
            pars["Asym"] <= 0 || pars["scal"] >= 0) next
        pars <- c(Asym = unname(pars["Asym"]), xmid = 50 - 5 * unname(pars["xmid"]),
                  scal = -5 * unname(pars["scal"]))
        predict_curve <- function(x) {
            unname(pars["Asym"] * stats::plogis((x - pars["xmid"]) / pars["scal"]))
        }
        derivative <- function(x) {
            p <- stats::plogis((x - pars["xmid"]) / pars["scal"])
            unname(pars["Asym"] * p * (1 - p) / pars["scal"])
        }
        return(list(predict = predict_curve, derivative = derivative, pars = pars,
                    algorithm = algorithm, edf = 3, gam_k = NA_integer_,
                    message = paste(messages[nzchar(messages)], collapse = " | ")))
    }
    stop("SSlogis failed to converge to an increasing curve. ",
         paste(messages[nzchar(messages)], collapse = " | "))
}

# Adapt 4.4's Gaussian GAM / GCV.Cp / adaptive k to regional fractions.
# The response is y in 0--1; log1p(track length) is not used for fractions.
.tm_gam_low <- function(d, gam_k, gam_method, gam_fixed_df) {
    k_use <- min(gam_k, length(unique(d$x)) - 1L, nrow(d) - 2L)
    if (k_use < 3L) stop("Insufficient distinct progress values for GAM.")
    s <- mgcv::s
    fit <- mgcv::gam(y ~ s(x, k = k_use, fx = gam_fixed_df), data = d,
                     family = stats::gaussian(), method = gam_method)
    if (!isTRUE(fit$converged)) stop("GAM did not converge.")
    predict_curve <- function(x) as.numeric(stats::predict(
        fit, newdata = data.frame(x = x), type = "response"))
    derivative <- function(x) {
        h <- 0.001
        lo <- pmax(0, x - h)
        hi <- pmin(100, x + h)
        (predict_curve(hi) - predict_curve(lo)) / (hi - lo)
    }
    algorithm <- if (gam_fixed_df) paste0("mgcv_gaussian_fixed_df_k", k_use) else
        paste0("mgcv_gaussian_", gam_method)
    list(predict = predict_curve, derivative = derivative, pars = NULL,
         algorithm = algorithm, edf = sum(fit$edf),
         gam_k = k_use, message = "")
}

.tm_gam_fast <- function(d, gam_k, gam_method, gam_fixed_df) {
    fit <- .tm_gam_low(d, gam_k, gam_method, gam_fixed_df)
    fit$algorithm <- sub("^mgcv_gaussian_", "mgcv_gaussian_fast_", fit$algorithm)
    fit
}

.tm_fit <- function(d, model_type, anchors, gam_k, gam_method, gam_fixed_df) {
    if (model_type == "sslogis") return(.tm_sslogis(d, anchors))
    if (model_type == "gam_low") return(.tm_gam_low(d, gam_k, gam_method, gam_fixed_df))
    if (model_type == "gam_fast") return(.tm_gam_fast(d, gam_k, gam_method, gam_fixed_df))
    stop("Unsupported model_type: ", model_type)
}

# Never substitute the nearest boundary for a threshold that was not crossed.
.tm_features_low <- function(fit, model_type, domain, grid_points, bootstrap_only = FALSE) {
    targets <- if (bootstrap_only) c(t10 = .1, t50 = .5, t90 = .9) else
        c(t10 = .1, t25 = .25, t50 = .5, t75 = .75, t90 = .9)
    grid <- seq(domain[1L], domain[2L], length.out = grid_points)
    y <- fit$predict(grid)
    if (any(!is.finite(y))) stop("Non-finite curve predictions.")
    if (!bootstrap_only) {
        full_grid <- seq(0, 100, length.out = grid_points)
        full_y <- fit$predict(full_grid)
        if (any(!is.finite(full_y))) stop("Non-finite curve predictions.")
    }
    times <- stats::setNames(rep(NA_real_, length(targets)), names(targets))
    status <- stats::setNames(rep("not_crossed", length(targets)), names(targets))
    for (key in names(targets)) {
        p <- targets[key]
        if (model_type == "sslogis") {
            a <- fit$pars["Asym"]
            if (p >= a) {
                status[key] <- "above_asymptote"
                next
            }
            root <- unname(fit$pars["xmid"] + fit$pars["scal"] * stats::qlogis(p / a))
            if (root < domain[1L]) status[key] <- "before_observed_range" else
            if (root > domain[2L]) status[key] <- "after_observed_range" else {
                times[key] <- root
                status[key] <- "ok"
            }
        } else {
            delta <- y - p
            exact <- which(abs(delta) < 1e-10)
            crossing <- which(delta[-length(delta)] * delta[-1L] < 0)
            roots <- c(grid[exact], vapply(crossing, function(j) {
                stats::uniroot(function(x) fit$predict(x) - p,
                               c(grid[j], grid[j + 1L]), tol = 1e-7)$root
            }, numeric(1)))
            roots <- sort(unique(round(roots, 7)))
            if (!length(roots)) next
            if (length(roots) > 1L) {
                status[key] <- "multiple_crossings"
            } else if (fit$derivative(roots) <= 1e-8) {
                status[key] <- "non_increasing_crossing"
            } else {
                times[key] <- roots
                status[key] <- "ok"
            }
        }
    }
    window <- unname(times["t90"] - times["t10"])
    if (is.finite(window) && window < 0) window <- NA_real_
    # Bootstrap consumes only t50 and t90-t10; avoid unused dense predictions,
    # derivatives, AUC, and t25/t75 root searches for every replicate.
    if (bootstrap_only) return(list(times = times, window = window))
    slopes <- fit$derivative(grid)
    if (model_type == "sslogis") {
        peak_x <- min(domain[2L], max(domain[1L], fit$pars["xmid"]))
        max_slope <- fit$derivative(peak_x)
    } else max_slope <- max(slopes)
    auc <- sum(diff(full_grid) * (full_y[-1L] + full_y[-length(full_y)]) / 2)
    list(times = times, status = status, window = window,
         max_slope = max_slope, auc = auc, auc_normalized = auc / 100,
         non_monotonic = any(slopes < -1e-6),
         out_of_bounds = any(full_y < -1e-6 | full_y > 1 + 1e-6),
         grid = full_grid, curve = full_y)
}

# Fast GAM feature extraction follows 4.9's dense-grid strategy: predict a
# combined grid once and derive crossings directly from adjacent grid points.
# This avoids repeated predict.gam() calls inside uniroot(). At the default
# grid_points=1001, crossing interpolation is resolved to about 0.1 progress
# percentage point; gam_low retains the exact uniroot implementation.
.tm_grid_derivative <- function(grid, y) {
    interval_slope <- diff(y) / diff(grid)
    if (length(interval_slope) == 1L) return(rep(interval_slope, 2L))
    c(interval_slope[1L],
      (interval_slope[-1L] + interval_slope[-length(interval_slope)]) / 2,
      interval_slope[length(interval_slope)])
}

.tm_grid_roots <- function(grid, y, target) {
    delta <- y - target
    exact <- which(abs(delta) < 1e-10)
    crossing <- which(delta[-length(delta)] * delta[-1L] < 0)
    roots <- grid[exact]
    if (length(crossing)) {
        interpolated <- grid[crossing] +
            (target - y[crossing]) * (grid[crossing + 1L] - grid[crossing]) /
            (y[crossing + 1L] - y[crossing])
        roots <- c(roots, interpolated)
    }
    sort(unique(round(roots, 7)))
}

.tm_features_fast <- function(fit, domain, grid_points, bootstrap_only = FALSE) {
    targets <- if (bootstrap_only) c(t10 = .1, t50 = .5, t90 = .9) else
        c(t10 = .1, t25 = .25, t50 = .5, t75 = .75, t90 = .9)
    grid <- seq(domain[1L], domain[2L], length.out = grid_points)
    if (bootstrap_only) {
        y <- fit$predict(grid)
        slopes <- .tm_grid_derivative(grid, y)
    } else {
        full_grid <- seq(0, 100, length.out = grid_points)
        h <- 0.001
        lo <- pmax(0, grid - h)
        hi <- pmin(100, grid + h)
        n <- length(grid)
        combined <- fit$predict(c(grid, full_grid, lo, hi))
        y <- combined[seq_len(n)]
        full_y <- combined[n + seq_len(n)]
        y_lo <- combined[2L * n + seq_len(n)]
        y_hi <- combined[3L * n + seq_len(n)]
        slopes <- (y_hi - y_lo) / (hi - lo)
        if (any(!is.finite(full_y))) stop("Non-finite curve predictions.")
    }
    if (any(!is.finite(y)) || any(!is.finite(slopes))) {
        stop("Non-finite curve predictions or derivatives.")
    }
    times <- stats::setNames(rep(NA_real_, length(targets)), names(targets))
    status <- stats::setNames(rep("not_crossed", length(targets)), names(targets))
    for (key in names(targets)) {
        roots <- .tm_grid_roots(grid, y, targets[key])
        if (!length(roots)) next
        if (length(roots) > 1L) {
            status[key] <- "multiple_crossings"
        } else {
            root_slope <- stats::approx(grid, slopes, xout = roots,
                                        rule = 2, ties = "ordered")$y
            if (!is.finite(root_slope) || root_slope <= 1e-8) {
                status[key] <- "non_increasing_crossing"
            } else {
                times[key] <- roots
                status[key] <- "ok"
            }
        }
    }
    window <- unname(times["t90"] - times["t10"])
    if (is.finite(window) && window < 0) window <- NA_real_
    if (bootstrap_only) return(list(times = times, window = window))
    auc <- sum(diff(full_grid) * (full_y[-1L] + full_y[-length(full_y)]) / 2)
    list(times = times, status = status, window = window,
         max_slope = max(slopes), auc = auc, auc_normalized = auc / 100,
         non_monotonic = any(slopes < -1e-6),
         out_of_bounds = any(full_y < -1e-6 | full_y > 1 + 1e-6),
         grid = full_grid, curve = full_y)
}

.tm_features <- function(fit, model_type, domain, grid_points, bootstrap_only = FALSE) {
    if (model_type == "gam_fast") {
        .tm_features_fast(fit, domain, grid_points, bootstrap_only)
    } else {
        .tm_features_low(fit, model_type, domain, grid_points, bootstrap_only)
    }
}

#' Fit a replication tempo curve to every region
#'
#' Fits one curve per row of a [rep_tempo()] table (regional replication
#' percentage against genome-wide replication progress) and extracts timing
#' features from it.
#'
#' `"sslogis"` is an unweighted three-parameter logistic NLS
#' (`Asym * plogis((x - xmid) / scal)`) with two optional fixed end-point
#' observations (0 at 0 percent, 1 at 100 percent) and a `plinear` fallback.
#' `"gam_low"` and `"gam_fast"` fit the same Gaussian `mgcv` GAM `y ~ s(x)`
#' and differ only in how thresholds are extracted: `"gam_low"` refines them
#' with `uniroot`, `"gam_fast"` interpolates on a dense prediction grid.
#' GAM fits are neither bounded nor monotonic; QC flags record this.
#'
#' Thresholds (`t10` ... `t90`) are real increasing crossings inside each
#' region's observed progress range; thresholds that are not crossed, or are
#' crossed more than once, are `NA`. They are in genome-wide percentage
#' points, as is `window_10_90 = t90 - t10`.
#'
#' The bootstrap resamples the observed (progress, percentage) points with
#' replacement and refits the model. It is a group-level interval, not a
#' cell-level one. A CI is reported only when at least
#' `max(20, ceiling(0.7 * n_boot))` replicates give a valid estimate. Row `i`
#' uses seed `seed + i - 1`, so serial and parallel runs agree. The caller's
#' RNG state is restored.
#'
#' @param df_rep_tempo Output of [rep_tempo()]: numeric regions by progress
#'   groups, values in 0-100, with unique row names.
#' @param rep_percent Optional sorted progress values in 0-100, one per
#'   column, for a plain matrix without the `rep_bins` attribute.
#' @param model_type `"sslogis"` (default), `"gam_low"` or `"gam_fast"`.
#' @param n_boot Number of bootstrap replicates; 0 (default) disables the
#'   bootstrap. 200 replicates take hours genome-wide.
#' @param conf_level Confidence level of the percentile intervals.
#' @param seed Non-negative integer seed for the bootstrap.
#' @param gam_k Maximum basis dimension of the GAM smooth (reduced
#'   automatically for regions with few points).
#' @param gam_method Smoothing-parameter method, `"GCV.Cp"` or `"REML"`.
#' @param gam_fixed_df If `TRUE`, fit an unpenalised regression spline with
#'   fixed degrees of freedom (`fx = TRUE`).
#' @param anchors Add the two fixed end-point observations to the SSlogis
#'   fit.
#' @param grid_points Number of grid points on 0-100 used for curves, AUC and
#'   grid-based thresholds.
#' @param curve_ids Region IDs whose fitted curves are returned in the
#'   `curves` attribute; `NULL` means all regions.
#' @param verbose Print progress messages. Default `FALSE`.
#' @param BPPARAM A `BiocParallel::BiocParallelParam`; regions are fitted
#'   independently. Worker processes are limited to one BLAS/OpenMP thread.
#'
#' @return A data.frame with one row per region: `region_id`, `model_type`,
#'   `fit_status`, `t10`/`t25`/`t50`/`t75`/`t90` with their status columns,
#'   `window_10_90`, `max_slope`, `auc`, `auc_normalized`, fit diagnostics
#'   (`rmse`, `r_squared`, `deviance`, residual summaries), SSlogis
#'   parameters, QC flags and bootstrap intervals. Attributes: `fitted` and
#'   `residuals` (percent, at the input progress values), `curves`,
#'   `rep_bins`, `range_info`, `settings` and `units`.
#' @seealso [rep_tempo()], [subset_rep_tempo()]
#' @export
fit_rep_tempo <- function(df_rep_tempo, rep_percent = NULL,
                          model_type = c("sslogis", "gam_low", "gam_fast"),
                          n_boot = 0L, conf_level = 0.95, seed = 1L,
                          gam_k = 10L, gam_method = c("GCV.Cp", "REML"),
                          gam_fixed_df = FALSE,
                          anchors = TRUE, grid_points = 1001L,
                          curve_ids = NULL, verbose = FALSE,
                          BPPARAM = BiocParallel::SerialParam()) {
    if (!is.null(BPPARAM) && !inherits(BPPARAM, "BiocParallelParam")) {
        stop("'BPPARAM' must be a BiocParallelParam object.")
    }
    # A serial backend runs the plain loop of the original (bpparam = NULL).
    bpparam <- if (is.null(BPPARAM) || inherits(BPPARAM, "SerialParam")) NULL else BPPARAM
    model_type <- match.arg(model_type)
    gam_method <- match.arg(gam_method)
    .tm_integer(n_boot, 0, "n_boot")
    .tm_integer(gam_k, 3, "gam_k")
    .tm_integer(grid_points, 101, "grid_points")
    .tm_integer(seed, 0, "seed")
    if (!is.numeric(conf_level) || length(conf_level) != 1L ||
        !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
        stop("conf_level must be in (0,1).")
    }
    if (!is.logical(anchors) || length(anchors) != 1L || is.na(anchors) ||
        !is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop("anchors and verbose must be TRUE or FALSE.")
    }
    if (!is.logical(gam_fixed_df) || length(gam_fixed_df) != 1L || is.na(gam_fixed_df)) {
        stop("gam_fixed_df must be TRUE or FALSE.")
    }
    if (model_type %in% c("gam_low", "gam_fast") &&
        !requireNamespace("mgcv", quietly = TRUE)) {
        stop("Package 'mgcv' is required for GAM.")
    }
    if (!(is.matrix(df_rep_tempo) || is.data.frame(df_rep_tempo))) {
        stop("df_rep_tempo must be a numeric data.frame or matrix.")
    }
    y <- as.matrix(df_rep_tempo)
    if (is.data.frame(df_rep_tempo) && !nrow(y) &&
        all(vapply(df_rep_tempo, is.numeric, logical(1)))) storage.mode(y) <- "double"
    if (!is.numeric(y) || !ncol(y) ||
        any(!is.na(y) & (!is.finite(y) | y < 0 | y > 100))) {
        stop("df_rep_tempo must contain percentages in 0--100 or NA.")
    }
    ids <- if (!nrow(y)) character() else rownames(y)
    .tm_ids(ids, nrow(y), "Region row names")
    bin_ids <- colnames(y)
    if (is.null(bin_ids)) bin_ids <- sprintf("rep_bin_%04d", seq_len(ncol(y)))
    .tm_ids(bin_ids, ncol(y), "Column names")
    dimnames(y) <- list(ids, bin_ids)
    bins <- attr(df_rep_tempo, "rep_bins")
    if (is.null(rep_percent)) {
        if (!is.data.frame(bins) || !all(c("rep_bin", "rep_percent") %in% names(bins))) {
            stop("Supply rep_percent or valid rep_bins metadata.")
        }
        .tm_ids(bins$rep_bin, nrow(bins), "rep_bins IDs")
        ix <- match(bin_ids, bins$rep_bin)
        if (anyNA(ix)) stop("rep_bins metadata does not match matrix columns.")
        bins <- bins[ix, , drop = FALSE]
        rep_percent <- bins$rep_percent
    } else bins <- data.frame(rep_bin = bin_ids, rep_percent = rep_percent)
    if (!is.numeric(rep_percent) || length(rep_percent) != ncol(y) ||
        any(!is.finite(rep_percent) | rep_percent < 0 | rep_percent > 100) ||
        is.unsorted(rep_percent)) stop("rep_percent must be sorted, finite, in 0--100.")
    counts <- attr(df_rep_tempo, "n_valid_cells")
    if (!is.null(counts)) {
        if (!is.matrix(counts) || !is.numeric(counts)) stop("Invalid n_valid_cells matrix.")
        .tm_ids(if (!nrow(counts)) character() else rownames(counts),
                nrow(counts), "n_valid_cells row names")
        .tm_ids(colnames(counts), ncol(counts), "n_valid_cells column names")
        ri <- match(ids, rownames(counts)); ci <- match(bin_ids, colnames(counts))
        if (anyNA(ri) || anyNA(ci)) stop("n_valid_cells IDs do not match the input.")
        counts <- counts[ri, ci, drop = FALSE]
        if (any(!is.finite(counts) | counts < 0 | counts != floor(counts)) ||
            any(is.finite(y) & counts == 0)) stop("Invalid n_valid_cells counts.")
    }
    if (is.null(curve_ids)) curve_ids <- ids
    if (!is.character(curve_ids) || anyNA(curve_ids) || anyDuplicated(curve_ids) ||
        any(!curve_ids %in% ids)) stop("curve_ids must be unique input region IDs.")
    if (n_boot > 0L || !is.null(bpparam)) {
        old_rng_kind <- RNGkind()
        had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
        if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
        on.exit({
            do.call(RNGkind, as.list(old_rng_kind))
            if (had_seed) assign(".Random.seed", old_seed, envir = globalenv()) else
                if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
                    rm(".Random.seed", envir = globalenv())
        }, add = TRUE)
    }
    if (!is.null(bpparam) && length(ids)) {
        # Fork workers share the read-only input; each result contains only one
        # region. Do not ship unrelated genome-wide objects or fitted models.
        controller_pid <- Sys.getpid()
        worker_ready <- FALSE
        fit_one <- function(i) {
            if (Sys.getpid() != controller_pid && !worker_ready) {
                .tm_worker_threads()
                worker_ready <<- TRUE
            }
            part <- as.data.frame(y[i, , drop = FALSE])
            if (!is.null(counts)) attr(part, "n_valid_cells") <- counts[i, , drop = FALSE]
            fit_rep_tempo(
                part, rep_percent = rep_percent, model_type = model_type,
                n_boot = n_boot, conf_level = conf_level,
                seed = as.integer((as.double(seed) + i - 1) %% .Machine$integer.max),
                gam_k = gam_k, gam_method = gam_method, gam_fixed_df = gam_fixed_df,
                anchors = anchors,
                grid_points = grid_points, curve_ids = intersect(ids[i], curve_ids),
                verbose = FALSE, BPPARAM = NULL
            )
        }
        if (verbose) message(model_type, ": ", length(ids), " regions, ",
                             BiocParallel::bpworkers(bpparam), " workers")
        results <- BiocParallel::bplapply(seq_along(ids), fit_one, BPPARAM = bpparam)
        # Combine columns once instead of repeatedly appending data.frames.
        out <- as.data.frame(stats::setNames(lapply(names(results[[1L]]), function(key) {
            unlist(lapply(results, function(z) z[[key]]), use.names = FALSE)
        }), names(results[[1L]])), stringsAsFactors = FALSE)
        rownames(out) <- ids
        attr(out, "fitted") <- do.call(rbind, lapply(results, attr, "fitted"))
        attr(out, "residuals") <- do.call(rbind, lapply(results, attr, "residuals"))
        curves <- matrix(NA_real_, length(curve_ids), grid_points,
                         dimnames = list(curve_ids, NULL))
        for (j in seq_along(curve_ids)) {
            curves[j, ] <- attr(results[[match(curve_ids[j], ids)]],
                                "curves")$fraction[1L, ]
        }
        attr(out, "curves") <- list(rep_percent = seq(0, 100, length.out = grid_points),
                                    fraction = curves)
        attr(out, "rep_bins") <- bins
        attr(out, "range_info") <- attr(df_rep_tempo, "range_info")
        settings <- attr(results[[1L]], "settings")
        settings$seed <- seed
        attr(out, "settings") <- settings
        attr(out, "units") <- attr(results[[1L]], "units")
        attr(out, "execution") <- list(backend = class(bpparam)[1L],
                                       workers = BiocParallel::bpworkers(bpparam))
        return(out)
    }
    out <- data.frame(region_id = ids, stringsAsFactors = FALSE)
    numeric_cols <- c("t10", "t25", "t50", "t75", "t90", "window_10_90",
        "max_slope", "auc", "fit_quality", "t50_ci_low", "t50_ci_high",
        "window_10_90_ci_low", "window_10_90_ci_high", "auc_normalized",
        "midpoint", "slope", "dispersion", "replication_window",
        "Asym", "xmid", "scal", "rmse", "r_squared", "deviance",
        "residual_mean", "residual_sd", "residual_lag1_cor", "edf", "gam_k",
        "n_cells", "rep_min", "rep_max")
    for (key in numeric_cols) out[[key]] <- rep(NA_real_, length(ids))
    out$model_type <- rep(model_type, length(ids))
    out$fit_status <- rep("no_data", length(ids))
    out$fit_algorithm <- rep(NA_character_, length(ids))
    out$fit_message <- rep("", length(ids))
    for (key in c("t10", "t25", "t50", "t75", "t90")) {
        out[[paste0(key, "_status")]] <- rep("no_data", length(ids))
    }
    for (key in c("n_points", "n_unique_progress", "bootstrap_success",
                  "bootstrap_t50_valid", "bootstrap_window_valid")) out[[key]] <- integer(length(ids))
    for (key in c("non_monotonic", "out_of_bounds", "auc_extrapolated")) {
        out[[key]] <- rep(NA, length(ids))
    }
    out$bootstrap_requested <- rep(as.integer(n_boot), length(ids))
    for (key in c("t50_ci_status", "window_10_90_ci_status")) {
        out[[key]] <- rep(if (n_boot) "not_fitted" else "disabled", length(ids))
    }
    fitted <- residuals <- matrix(NA_real_, nrow(y), ncol(y), dimnames = dimnames(y))
    curve_grid <- seq(0, 100, length.out = grid_points)
    curves <- matrix(NA_real_, length(curve_ids), grid_points,
                     dimnames = list(curve_ids, NULL))
    min_points <- if (model_type %in% c("gam_low", "gam_fast")) 5L else 3L
    min_boot <- max(20L, ceiling(.7 * n_boot))
    for (i in seq_along(ids)) {
        if (n_boot > 0L) {
            set.seed(as.integer((as.double(seed) + i - 1) %% .Machine$integer.max),
                     kind = "Mersenne-Twister", normal.kind = "Inversion",
                     sample.kind = "Rejection")
        }
        if (verbose && (i == 1L || i %% 250L == 0L || i == length(ids))) {
            message(model_type, ": region ", i, "/", length(ids))
        }
        keep <- is.finite(y[i, ])
        if (!any(keep)) next
        d <- data.frame(x = rep_percent[keep], y = unname(y[i, keep]) / 100)
        out$n_points[i] <- nrow(d)
        out$n_unique_progress[i] <- length(unique(d$x))
        out$n_cells[i] <- if (is.null(counts)) NA_real_ else sum(counts[i, keep])
        out$rep_min[i] <- min(d$x); out$rep_max[i] <- max(d$x)
        failure <- if (length(unique(d$x)) < min_points) "insufficient_points" else
                   if (diff(range(d$y)) < 1e-10) "constant_response" else NULL
        attempt <- if (is.null(failure)) .tm_try({
            fit <- .tm_fit(d, model_type, anchors, gam_k, gam_method, gam_fixed_df)
            features <- .tm_features(fit, model_type, range(d$x), grid_points)
            list(fit = fit, features = features)
        }) else list(value = NULL, message = failure)
        if (is.null(attempt$value)) {
            out$fit_status[i] <- if (is.null(failure)) "fit_failed" else failure
            out$fit_message[i] <- attempt$message
            for (key in c("t10", "t25", "t50", "t75", "t90")) {
                out[[paste0(key, "_status")]][i] <- out$fit_status[i]
            }
            next
        }
        fit <- attempt$value$fit; f <- attempt$value$features
        for (key in names(f$times)) {
            out[[key]][i] <- f$times[key]
            out[[paste0(key, "_status")]][i] <- f$status[key]
        }
        out$window_10_90[i] <- f$window
        for (key in c("max_slope", "auc", "auc_normalized",
                      "non_monotonic", "out_of_bounds")) out[[key]][i] <- f[[key]]
        out$auc_extrapolated[i] <- min(d$x) > 0 || max(d$x) < 100
        out$fit_status[i] <- if (f$non_monotonic || f$out_of_bounds) "curve_qc_flag" else
                            if (any(f$status != "ok")) "partial_thresholds" else "ok"
        out$fit_algorithm[i] <- fit$algorithm
        out$fit_message[i] <- paste(Filter(nzchar, c(attempt$message, fit$message)),
                                    collapse = " | ")
        out$edf[i] <- fit$edf; out$gam_k[i] <- fit$gam_k
        if (model_type == "sslogis") {
            for (key in names(fit$pars)) out[[key]][i] <- fit$pars[key]
            out$midpoint[i] <- out$t50[i]
            out$slope[i] <- 1 / fit$pars["scal"]
            out$dispersion[i] <- fit$pars["scal"]
            out$replication_window[i] <- f$window
        }
        fitted[i, ] <- 100 * fit$predict(rep_percent)
        residuals[i, keep] <- y[i, keep] - fitted[i, keep]
        error <- residuals[i, keep] / 100
        sse <- sum(error^2); sst <- sum((d$y - mean(d$y))^2)
        out$rmse[i] <- 100 * sqrt(mean(error^2))
        out$r_squared[i] <- out$fit_quality[i] <- if (sst > 0) 1 - sse / sst else NA_real_
        out$deviance[i] <- sse
        out$residual_mean[i] <- 100 * mean(error)
        out$residual_sd[i] <- 100 * stats::sd(error)
        if (length(error) > 2L && stats::sd(error[-1L]) > 0 &&
            stats::sd(error[-length(error)]) > 0) {
            out$residual_lag1_cor[i] <- stats::cor(error[-1L], error[-length(error)])
        }
        curve_row <- match(ids[i], curve_ids)
        if (!is.na(curve_row)) curves[curve_row, ] <- f$curve
        if (!n_boot) next
        boot <- matrix(NA_real_, n_boot, 2L)
        for (b in seq_len(n_boot)) {
            bd <- d[sample.int(nrow(d), nrow(d), replace = TRUE), , drop = FALSE]
            if (length(unique(bd$x)) < min_points || diff(range(bd$y)) < 1e-10) next
            bf <- .tm_try({
                bfit <- .tm_fit(bd, model_type, anchors, gam_k, gam_method, gam_fixed_df)
                # Restrict each replicate to its own observed support.
                .tm_features(bfit, model_type, range(bd$x), grid_points, bootstrap_only = TRUE)
            })$value
            if (is.null(bf)) next
            out$bootstrap_success[i] <- out$bootstrap_success[i] + 1L
            boot[b, ] <- c(bf$times["t50"], bf$window)
        }
        valid_n <- colSums(is.finite(boot))
        out$bootstrap_t50_valid[i] <- valid_n[1L]
        out$bootstrap_window_valid[i] <- valid_n[2L]
        for (j in 1:2) {
            key <- c("t50", "window_10_90")[j]
            ci_status <- if (!is.finite(out[[key]][i])) "estimate_unavailable" else
                         if (valid_n[j] < min_boot) "insufficient_bootstrap" else "ok"
            out[[paste0(key, "_ci_status")]][i] <- ci_status
            if (ci_status == "ok") {
                ci <- stats::quantile(boot[, j], c((1 - conf_level) / 2,
                                                  (1 + conf_level) / 2), na.rm = TRUE)
                out[[paste0(key, "_ci_low")]][i] <- ci[1L]
                out[[paste0(key, "_ci_high")]][i] <- ci[2L]
            }
        }
    }
    first <- c("region_id", "t10", "t25", "t50", "t75", "t90", "window_10_90",
               "max_slope", "auc", "model_type", "fit_quality", "t50_ci_low", "t50_ci_high")
    out <- out[, c(first, setdiff(names(out), first)), drop = FALSE]
    rownames(out) <- ids
    attr(out, "fitted") <- fitted
    attr(out, "residuals") <- residuals
    attr(out, "curves") <- list(rep_percent = curve_grid, fraction = curves)
    attr(out, "rep_bins") <- bins
    attr(out, "range_info") <- attr(df_rep_tempo, "range_info")
    attr(out, "settings") <- list(model_type = model_type, n_boot = n_boot,
        conf_level = conf_level, seed = seed, anchors = model_type == "sslogis" && anchors,
        gam_k = gam_k, gam_method = gam_method, gam_fixed_df = gam_fixed_df,
        gam_smoothing = if (!model_type %in% c("gam_low", "gam_fast")) NA_character_ else
            if (gam_fixed_df) "fixed degrees of freedom; no smoothing-parameter estimation" else
                paste0("penalized smooth; method=", gam_method),
        gam_feature_extraction = if (model_type == "gam_fast")
            "batched prediction; grid-linear crossing interpolation" else
            if (model_type == "gam_low") "repeated prediction; uniroot crossing refinement" else
                NA_character_,
        grid_points = grid_points,
        bootstrap = "pairs of observed rep_bin points; equal weights; percentile CI",
        bootstrap_min_valid = min_boot,
        bootstrap_rng = "Mersenne-Twister; per-row seed + row_index - 1 (mod integer.max)",
        threshold_domain = "observed range",
        auc_domain = c(0, 100), fit_quality = "observation-only R-squared",
        response = "mean regional replicated fraction, not whole-region completion counts")
    attr(out, "units") <- list(time = "genome-wide percent (0--100), not hours",
        window = "genome-wide percentage points",
        max_slope = "regional fraction per genome-wide percentage point",
        auc = "integral of regional fraction over genome-wide percent 0--100",
        rmse = "regional percentage points", deviance = "Gaussian SSE in fraction units")
    attr(out, "execution") <- list(backend = "serial", workers = 1L)
    out
}
