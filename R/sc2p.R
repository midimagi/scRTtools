# Gene-level SC2P phase summaries (from rcode/analysis/calculate_sc2p.R).

#' @keywords internal
#' @noRd
.check_sc2p_counts <- function(Y) {
    if (!(is.matrix(Y) || is.data.frame(Y)) || length(dim(Y)) != 2L) {
        stop("Y must be a numeric matrix or data.frame.", call. = FALSE)
    }
    if (nrow(Y) == 0L || ncol(Y) == 0L) {
        stop("Y must contain at least one gene and one cell.", call. = FALSE)
    }

    numeric_input <- if (is.data.frame(Y)) {
        all(vapply(Y, is.numeric, FUN.VALUE = logical(1)))
    } else {
        is.numeric(Y)
    }
    if (!numeric_input) {
        stop("All values in Y must be numeric counts.", call. = FALSE)
    }
    if (anyNA(Y)) {
        stop("Y must not contain missing values.", call. = FALSE)
    }

    gene_names <- rownames(Y)
    cell_names <- colnames(Y)
    if (is.null(gene_names) || anyNA(gene_names) || any(gene_names == "") ||
        anyDuplicated(gene_names)) {
        stop("Y must have unique, non-empty gene names as row names.", call. = FALSE)
    }
    if (is.null(cell_names) || anyNA(cell_names) || any(cell_names == "") ||
        anyDuplicated(cell_names)) {
        stop("Y must have unique, non-empty cell names as column names.", call. = FALSE)
    }

    invisible(TRUE)
}

#' Robust zero-inflated Poisson fit used for SC2P Phase I
#'
#' Regresses `log(frequency) + log(x!)` on the non-zero count values `x`
#' (weights `1 / exp(x)`) to estimate the Poisson mean, then derives the
#' zero-inflation proportion from the observed number of zeros.
#'
#' @param x Numeric vector of counts for one gene.
#' @return Named numeric vector `c(p0, mu)`; both `NA` when `x` is empty, all
#'   zero, or the regression is not estimable.
#' @keywords internal
#' @noRd
.sc2p_robust_poisson <- function(x) {
    tt <- table(x)
    n0 <- sum(x == 0)

    # After SC2P_gene_unit() removes counts above the threshold, a highly
    # expressed gene can have no observations (or only zeros) left. The fit is
    # undefined then; NA lets the caller fall back to a threshold of 1.
    if (length(x) == 0L || n0 == length(x)) {
        return(c(p0 = NA_real_, mu = NA_real_))
    }

    if (names(tt)[1] == "0") {
        xx <- as.numeric(names(tt))[-1]
        tt <- as.vector(tt)[-1]
    } else {
        xx <- as.numeric(names(tt))
        tt <- as.vector(tt)
    }

    tt <- log(tt) + log(gamma(1 + xx))
    beta <- stats::lm(tt ~ xx, weights = 1 / exp(xx))$coef
    if (is.na(as.numeric(beta)[2]) || is.na(as.numeric(beta)[1])) {
        return(c(p0 = NA_real_, mu = NA_real_))
    }

    mu_hat <- exp(beta[2])
    p0 <- (n0 - exp(beta[1])) /
        (exp(beta[1] + mu_hat) + n0 - exp(beta[1]))
    if (any(p0 < 0)) {
        warning("Proportion of zero inflation is negative", call. = FALSE)
    }
    p0 <- pmax(p0, 0)

    c(p0 = unname(p0), mu = unname(mu_hat))
}

#' Gene-level SC2P Phase I / Phase II summaries
#'
#' For every expressed gene, counts not above `yy_thres` are fitted with a
#' robust zero-inflated Poisson (Phase I, the background). A cell is in
#' Phase II for a gene when its count reaches the 90% quantile of that
#' Poisson (threshold 1 when the fit is undefined). The function reports the
#' fraction of Phase I cells and, among Phase II cells, the mean and
#' coefficient of variation of counts-per-10,000 expression (library size =
#' the cell's total over the expressed genes of `Y`).
#'
#' @param Y Numeric matrix or data.frame of counts, genes in rows and cells in
#'   columns, with unique row and column names and no missing values. Genes
#'   with zero total are dropped; cells with zero total are an error.
#' @param yy_thres Counts above this value are excluded from the Phase I fit.
#'
#' @return A data.frame with one row per expressed gene (row names = gene) and
#'   columns `gene`, `phase_I_ratio`, `phase_II_mean`, `phase_II_varcoff`.
#'   Genes with no Phase II cell have `NA` Phase II summaries; with a single
#'   Phase II cell the coefficient of variation is `NA`.
#'
#' @examples
#' set.seed(1)
#' counts <- matrix(rpois(40 * 60, lambda = rep(c(0.3, 2), each = 20)),
#'                  nrow = 40,
#'                  dimnames = list(paste0("g", 1:40), paste0("c", 1:60)))
#' head(suppressWarnings(SC2P_gene_unit(counts)))
#' @export
SC2P_gene_unit <- function(Y, yy_thres = 15) {
    .check_sc2p_counts(Y)

    Y <- Y[rowSums(Y) > 0, , drop = FALSE]
    if (nrow(Y) == 0L) {
        stop("Y does not contain any expressed genes.", call. = FALSE)
    }
    if (any(colSums(Y) == 0)) {
        stop("Y contains one or more cells with a total count of zero.", call. = FALSE)
    }
    Y <- as.matrix(Y)

    par1 <- apply(Y, 1, function(yy) {
        yy <- yy[yy <= yy_thres]
        .sc2p_robust_poisson(yy)
    })

    thres <- stats::qpois(0.9, par1[2L, ])
    thres[is.na(thres)] <- 1

    # thres has one value per gene (row), so it recycles down the columns.
    phase_II <- !(Y < thres)
    cp10k <- sweep(Y * 1e4, 2L, apply(Y, 2L, sum), "/")

    n_genes <- nrow(Y)
    phase_I_ratio <- numeric(n_genes)
    phase_II_mean <- rep(NA_real_, n_genes)
    phase_II_varcoff <- rep(NA_real_, n_genes)
    for (g in seq_len(n_genes)) {
        in_phase_II <- phase_II[g, ]
        phase_I_ratio[g] <- 1 - mean(in_phase_II + 0)
        if (any(in_phase_II)) {
            x <- cp10k[g, in_phase_II]
            phase_II_mean[g] <- mean(x)
            phase_II_varcoff[g] <- stats::sd(x) / phase_II_mean[g]
        }
    }

    data.frame(
        gene = rownames(Y),
        phase_I_ratio = phase_I_ratio,
        phase_II_mean = phase_II_mean,
        phase_II_varcoff = phase_II_varcoff,
        row.names = rownames(Y),
        stringsAsFactors = FALSE
    )
}
