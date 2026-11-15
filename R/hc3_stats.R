# Per-sample regression with HC3 standard errors, two-sample combination and
# design helpers (from jupyter/6.6.Chek1_paired_DNA_RNA_analysis.R and
# jupyter/6.7.expression_to_other_gene_replication.R; decisions H1-H12).

#' Z-score with a fixed value for constant input
#' @param x Numeric vector.
#' @param constant Value returned for every element when the standard
#'   deviation is zero or not finite.
#' @keywords internal
#' @noRd
.zscore <- function(x, constant = NA_real_) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(constant, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

#' Resolve `terms` to column indices of a design matrix
#' @keywords internal
#' @noRd
.hc3_terms <- function(X, terms) {
  if (!is.matrix(X) || !is.numeric(X)) stop("'X' must be a numeric matrix.")
  if (is.null(terms)) {
    j <- ncol(X)
    label <- if (is.null(colnames(X))) as.character(j) else colnames(X)[j]
    return(list(index = j, label = label))
  }
  j <- match(terms, colnames(X))
  if (anyNA(j)) {
    stop("Term(s) not found in colnames(X): ", paste(terms[is.na(j)], collapse = ", "))
  }
  list(index = j, label = as.character(terms))
}

#' Ordinary least squares with HC3 robust standard errors
#'
#' Fits `y ~ X` by QR and reports selected coefficients with
#' heteroscedasticity-consistent (HC3) standard errors and t-based p-values
#' (`n - k` degrees of freedom). Rows with a non-finite `y` or `X` are
#' dropped first.
#'
#' The fit is not attempted (`status = "rank_or_n"`) when `n <= k + 20` or
#' `X` is rank deficient. When the largest leverage reaches `max_leverage`
#' the result is `NA` with `status = "high_leverage"`; otherwise leverages
#' are capped at `1 - 1e-8`. A term whose standard error is zero or not
#' finite gets `status = "zero_se"`.
#'
#' @param X Numeric design matrix (include the intercept column yourself).
#' @param y Numeric response, one value per row of `X`.
#' @param terms Column names of `X` to report; `NULL` reports the last
#'   column.
#' @param max_leverage Largest leverage accepted. `Inf` never rejects a fit
#'   on leverage.
#'
#' @return A data.frame with one row per term: `term`, `beta`, `se`, `p`,
#'   `n`, `max_hat`, `status`.
#'
#' @examples
#' set.seed(1)
#' X <- cbind(`(Intercept)` = 1, age = rnorm(200), dose = rnorm(200))
#' y <- 0.5 * X[, "dose"] + rnorm(200, sd = 1 + abs(X[, "age"]))
#' fit_ols_hc3(X, y, terms = c("age", "dose"))
#' @export
fit_ols_hc3 <- function(X, y, terms = NULL, max_leverage = 0.99) {
  sel <- .hc3_terms(X, terms)
  if (length(y) != nrow(X)) stop("'y' must have one value per row of 'X'.")
  ok <- is.finite(y) & rowSums(!is.finite(X)) == 0L
  X <- X[ok, , drop = FALSE]
  y <- y[ok]
  n <- length(y)
  k <- ncol(X)
  bad <- function(status) {
    data.frame(term = sel$label, beta = NA_real_, se = NA_real_, p = NA_real_,
               n = n, max_hat = NA_real_, status = status,
               stringsAsFactors = FALSE)
  }
  q <- qr(X)
  if (n <= k + 20L || q$rank < k) return(bad("rank_or_n"))

  b <- qr.coef(q, y)
  e <- as.numeric(y - X %*% b)
  h <- rowSums(qr.Q(q)^2)
  max_hat <- max(h)
  if (max_hat >= max_leverage) return(bad("high_leverage"))
  h <- pmin(h, 1 - 1e-8)

  inv0 <- tcrossprod(backsolve(qr.R(q), diag(k)))
  inv <- inv0
  inv[q$pivot, q$pivot] <- inv0

  out <- bad("ok")
  out$max_hat <- max_hat
  for (i in seq_along(sel$index)) {
    j <- sel$index[i]
    v <- as.numeric(X %*% inv[, j])
    se <- sqrt(sum((v * e / (1 - h))^2))
    if (!is.finite(se) || se <= 0) {
      out$status[i] <- "zero_se"
      out$max_hat[i] <- NA_real_
      next
    }
    out$beta[i] <- unname(b[j])
    out$se[i] <- se
    out$p[i] <- 2 * stats::pt(-abs(b[j] / se), df = n - k)
  }
  out
}

#' Poisson regression with HC3 robust standard errors
#'
#' Fits a log-link Poisson model with an offset by `stats::glm.fit` and
#' reports selected coefficients with HC3 sandwich standard errors, which
#' tolerate a variance different from the mean, and normal p-values.
#'
#' Status codes: `"rank_or_n"` (`n <= k + 20` or rank-deficient `X`),
#' `"nonconverged"`, `"ill_conditioned"` (reciprocal condition number of
#' `X' diag(mu) X` below `min_rcond`), `"unstable"` (standard error zero or
#' not finite, or `|beta| > max_abs_beta`), `"ok"`.
#'
#' @inheritParams fit_ols_hc3
#' @param y Count response.
#' @param offset Offset on the log scale, e.g. `log(library_size / 1e4)`.
#' @param maxit Maximum number of IWLS iterations.
#' @param min_rcond,max_abs_beta Thresholds of the `"ill_conditioned"` and
#'   `"unstable"` checks.
#'
#' @return A data.frame with one row per term: `term`, `beta`, `se`, `p`,
#'   `n`, `status`.
#'
#' @examples
#' set.seed(1)
#' X <- cbind(`(Intercept)` = 1, dose = rnorm(300))
#' depth <- runif(300, 5e3, 2e4)
#' y <- rpois(300, depth / 1e4 * exp(0.2 + 0.3 * X[, "dose"]))
#' fit_poisson_hc3(X, y, offset = log(depth / 1e4))
#' @export
fit_poisson_hc3 <- function(X, y, offset, terms = NULL, maxit = 100,
                            min_rcond = 1e-12, max_abs_beta = 20) {
  sel <- .hc3_terms(X, terms)
  if (length(y) != nrow(X) || length(offset) != nrow(X)) {
    stop("'y' and 'offset' must have one value per row of 'X'.")
  }
  ok <- is.finite(y) & is.finite(offset) & rowSums(!is.finite(X)) == 0L
  X <- X[ok, , drop = FALSE]
  y <- y[ok]
  offset <- offset[ok]
  bad <- function(status) {
    data.frame(term = sel$label, beta = NA_real_, se = NA_real_, p = NA_real_,
               n = length(y), status = status, stringsAsFactors = FALSE)
  }
  if (nrow(X) <= ncol(X) + 20 || qr(X)$rank < ncol(X)) return(bad("rank_or_n"))
  fit <- suppressWarnings(stats::glm.fit(
    X, y, offset = offset, family = stats::poisson(),
    control = stats::glm.control(maxit = maxit)
  ))
  if (!isTRUE(fit$converged) || any(!is.finite(fit$coefficients))) {
    return(bad("nonconverged"))
  }
  mu <- fit$fitted.values
  bread <- crossprod(X, X * mu)
  if (!is.finite(rcond(bread)) || rcond(bread) < min_rcond) return(bad("ill_conditioned"))
  inv <- solve(bread)
  h <- pmin(mu * rowSums((X %*% inv) * X), 1 - 1e-8)

  out <- bad("ok")
  for (i in seq_along(sel$index)) {
    j <- sel$index[i]
    se <- sqrt(sum(((X %*% inv[, j]) * (y - mu) / (1 - h))^2))
    beta <- unname(fit$coefficients[j])
    if (!is.finite(se) || se <= 0 || abs(beta) > max_abs_beta) {
      out$status[i] <- "unstable"
      next
    }
    out$beta[i] <- beta
    out$se[i] <- se
    out$p[i] <- 2 * stats::pnorm(-abs(beta / se))
  }
  out
}

#' Combine the estimates of two samples
#'
#' Fixed-effect (inverse-variance) combination of per-sample estimates, one
#' pair per combination of the `by` columns, with a normal p-value, Cochran's
#' Q heterogeneity test (1 df) and `I2`. `replication_p` is the larger of the
#' two per-sample p-values when both estimates have the same non-zero sign
#' and 1 otherwise: evidence that both samples support the effect.
#'
#' A pair is failed when a standard error is missing or not positive or,
#' if `tab` has a `status` column, when a status is not `"ok"`. Combining two
#' samples does not add biological replicates.
#'
#' @param tab data.frame of per-sample fits with columns `beta`, `se`, `p`,
#'   the `by` columns, the sample column and optionally `status`, e.g. the
#'   stacked output of [fit_ols_hc3()]. Exactly two samples, each once per
#'   combination of `by`.
#' @param by Columns identifying one test.
#' @param sample_col Column holding the sample name.
#' @param failed `"keep"`: failed pairs stay in the result with `beta`/`se`
#'   `NA`, `p = 1`, `replication_p = 1` and `status = "fit_failed"`, so they
#'   count in later multiple-testing corrections; `"drop"` removes them.
#'
#' @return A data.frame with one row per test (in order of first
#'   appearance): the `by` columns, `beta`, `se`, `p`, `replication_p`,
#'   `same_direction`, `Q`, `heterogeneity_p`, `I2`, `status`, and
#'   `beta_<sample>`, `se_<sample>`, `p_<sample>` for both samples (sorted by
#'   name).
#'
#' @examples
#' fits <- data.frame(gene = rep(c("a", "b"), each = 2),
#'                    sample = rep(c("s1", "s2"), 2),
#'                    beta = c(0.5, 0.4, 0.3, -0.2), se = c(0.1, 0.12, 0.1, 0.1),
#'                    p = c(1e-6, 1e-3, 3e-3, 0.05))
#' combine_two_samples(fits, by = "gene")
#' @export
combine_two_samples <- function(tab, by, sample_col = "sample",
                                failed = c("keep", "drop")) {
  failed <- match.arg(failed)
  tab <- as.data.frame(tab)
  .check_columns(tab, c(by, sample_col, "beta", "se", "p"), "tab")
  sample_names <- sort(unique(as.character(tab[[sample_col]])))
  if (length(sample_names) != 2L) {
    stop("'tab' must hold exactly two samples; found ", length(sample_names), ".")
  }
  key <- do.call(paste, c(unname(as.list(tab[by])), sep = "\r"))
  smp <- as.character(tab[[sample_col]])
  if (anyDuplicated(paste(key, smp, sep = "\r"))) {
    stop("Each sample must occur once per combination of the 'by' columns.")
  }
  first <- which(!duplicated(key))
  i1 <- which(smp == sample_names[1])
  i2 <- which(smp == sample_names[2])
  i1 <- i1[match(key[first], key[i1])]
  i2 <- i2[match(key[first], key[i2])]
  if (anyNA(i1) || anyNA(i2)) {
    stop("Every combination of the 'by' columns needs both samples.")
  }

  n <- length(first)
  b1 <- tab$beta[i1]; b2 <- tab$beta[i2]
  s1 <- tab$se[i1]; s2 <- tab$se[i2]
  p1 <- tab$p[i1]; p2 <- tab$p[i2]
  ok <- is.finite(s1) & is.finite(s2) & s1 > 0 & s2 > 0
  if ("status" %in% names(tab)) {
    ok <- ok & tab$status[i1] %in% "ok" & tab$status[i2] %in% "ok"
  }

  beta <- se <- Q <- rep(NA_real_, n)
  for (g in which(ok)) {
    bb <- c(b1[g], b2[g])
    w <- 1 / c(s1[g], s2[g])^2
    beta[g] <- sum(w * bb) / sum(w)
    se[g] <- sqrt(1 / sum(w))
    Q[g] <- sum(w * (bb - beta[g])^2)
  }
  same <- ok & ((b1 > 0 & b2 > 0) | (b1 < 0 & b2 < 0))
  same[is.na(same)] <- FALSE

  out <- tab[first, by, drop = FALSE]
  out$beta <- beta
  out$se <- se
  out$p <- ifelse(ok, 2 * stats::pnorm(-abs(beta / se)), 1)
  out$replication_p <- ifelse(same, pmax(p1, p2), 1)
  out$same_direction <- same
  out$Q <- Q
  out$heterogeneity_p <- stats::pchisq(Q, df = 1, lower.tail = FALSE)
  out$I2 <- ifelse(ok, ifelse(Q > 0, pmax(0, (Q - 1) / Q), 0), NA_real_)
  out$status <- ifelse(ok, "ok", "fit_failed")
  out[[paste0("beta_", sample_names[1])]] <- b1
  out[[paste0("se_", sample_names[1])]] <- s1
  out[[paste0("p_", sample_names[1])]] <- p1
  out[[paste0("beta_", sample_names[2])]] <- b2
  out[[paste0("se_", sample_names[2])]] <- s2
  out[[paste0("p_", sample_names[2])]] <- p2
  if (failed == "drop") out <- out[ok, , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Model matrix without linearly dependent columns
#'
#' Builds `stats::model.matrix(formula, data)` and keeps a full-rank subset
#' of its columns (QR pivoting; original column order is kept). Useful when
#' factor levels or spline bases are empty or collinear within one sample.
#'
#' @param formula One-sided model formula. Functions used in it (e.g.
#'   `splines::ns`) must be visible from the formula's environment.
#' @param data data.frame with the variables of `formula`; missing or
#'   non-finite values are an error.
#'
#' @return A numeric matrix with one row per row of `data`.
#'
#' @examples
#' d <- data.frame(x = rnorm(30), g = gl(3, 10))
#' d$x2 <- 2 * d$x
#' colnames(full_rank_model_matrix(~ x + x2 + g, d))
#' @export
full_rank_model_matrix <- function(formula, data) {
  b <- stats::model.matrix(formula, data, na.action = stats::na.fail)
  if (nrow(b) != nrow(data) || !all(is.finite(b))) {
    stop("The model matrix has missing or non-finite values.")
  }
  q <- qr(b)
  b[, sort(q$pivot[seq_len(q$rank)]), drop = FALSE]
}

#' Add regulator and local-replication interaction terms to a design matrix
#'
#' Appends the columns of the per-gene model of the paired DNA-RNA analysis:
#' `regulator`, `regulator_progress`, `local_rep`, `local_progress`,
#' `local_progress2` and `interaction` (`local_rep * regulator`).
#'
#' @param B Base design matrix (nuisance covariates), one row per cell.
#' @param regulator Regulator expression per cell (e.g. z-scored).
#' @param local_rep Local replication state per cell (0 = unreplicated,
#'   1 = replicated, `NA` = missing).
#' @param progress_c Centred replication progress per cell.
#'
#' @return `B` with the six columns appended.
#' @export
add_interaction_terms <- function(B, regulator, local_rep, progress_c) {
  n <- nrow(B)
  if (length(regulator) != n || length(local_rep) != n || length(progress_c) != n) {
    stop("'regulator', 'local_rep' and 'progress_c' need one value per row of 'B'.")
  }
  cbind(
    B,
    regulator = regulator,
    regulator_progress = regulator * progress_c,
    local_rep = local_rep,
    local_progress = local_rep * progress_c,
    local_progress2 = local_rep * progress_c^2,
    interaction = local_rep * regulator
  )
}

#' Permute cells within sample and replication-progress bins
#'
#' Returns a permutation of the cell indices that only exchanges cells of
#' the same group whose progress falls in the same bin of width `bin_width`.
#' Use it to break the DNA-RNA pairing while keeping sample and progress
#' (almost) fixed: `rna_columns[perm]` against the original DNA columns.
#' Set a seed beforehand for reproducible permutations.
#'
#' @param group Sample (or other stratum) of each cell.
#' @param progress Replication progress of each cell.
#' @param bin_width Width of the progress bins.
#'
#' @return Integer vector `perm` of the same length; cell `i` receives the
#'   data of cell `perm[i]`.
#'
#' @examples
#' set.seed(1)
#' progress <- runif(20)
#' perm <- permute_within_bins(rep(c("a", "b"), 10), progress, 0.25)
#' max(abs(progress - progress[perm])) < 0.25
#' @export
permute_within_bins <- function(group, progress, bin_width) {
  if (length(group) != length(progress)) {
    stop("'group' and 'progress' must have the same length.")
  }
  if (anyNA(group) || any(!is.finite(progress))) {
    stop("'group' and 'progress' must not be missing.")
  }
  n <- length(progress)
  groups <- split(seq_len(n),
                  interaction(group, floor(progress / bin_width), drop = TRUE))
  perm <- seq_len(n)
  for (gr in groups) if (length(gr) > 1L) perm[gr] <- sample(gr)
  perm
}

#' Genome-bin identifier of a position
#'
#' @param chr Chromosome name(s).
#' @param pos 1-based position(s).
#' @param bin_width Bin width in bp.
#' @return Character vector `"<chr>-<start>-<end>"` of the fixed-width bin
#'   holding each position, always in plain (never scientific) notation.
#' @examples
#' make_bin("chr1", c(1, 20000, 20001, 99999))
#' @export
make_bin <- function(chr, pos, bin_width = 20000) {
  index <- floor((pos - 1) / bin_width)
  plain <- function(x) format(x, scientific = FALSE, trim = TRUE)
  out <- paste0(chr, "-", plain(index * bin_width + 1), "-",
                plain((index + 1) * bin_width))
  out[is.na(pos)] <- NA_character_
  out
}
