# Per-track table resembling calc_rep_track_lengths() output, with a
# replication-dependent median length (convex in rep_percent) and outliers.
synth_track_table <- function(n_cells = 160L, states = c("ctr", "ama"), seed = 31L,
                              ode = FALSE) {
  withr::with_seed(seed, {
    do.call(rbind, lapply(seq_along(states), function(k) {
      rp <- stats::runif(n_cells, 0.05, 0.95)
      do.call(rbind, lapply(seq_len(n_cells), function(i) {
        x <- rp[i] * 100
        mu <- if (ode) {
          12 - 0.01 * (x - 50) + 0.02 / 0.04 * expm1(0.04 * (x - 50))
        } else {
          11.5 + 0.004 * x + (0.00035 + 0.00005 * k) * x^2
        }
        if (i %% 37 == 0) mu <- mu + 2                      # outlier cell
        n_tr <- sample(5:15, 1)
        data.frame(
          chr = "chr1", start = 1L, end = 2L,
          width_bp = round(exp(mu + stats::rnorm(n_tr, 0, 0.4))),
          cell_id = sprintf("%s_c%03d", states[k], i),
          rep_percent = rp[i], state = states[k],
          stringsAsFactors = FALSE
        )
      }))
    }))
  })
}

# Notebook 2.3 ODE code (inline functions), verbatim apart from namespacing;
# rep_center is looked up in the returned environment.
notebook_ode_env <- function(df_true_length) {
  env <- new.env()
  env$rep_center <- stats::median(df_true_length$rep_percent, na.rm = TRUE)
  local(envir = env, {
    lmrob <- robustbase::lmrob
    nlrob <- robustbase::nlrob
    ode_log_length <- function(rep_percent, log_L_ref, alpha, log_beta, log_gamma) {
      beta <- exp(log_beta)
      gamma <- exp(log_gamma)
      rep_delta <- rep_percent - rep_center
      log_L_ref + alpha * rep_delta + beta / gamma * expm1(gamma * rep_delta)
    }
    get_ode_start <- function(data) {
      start_data <- transform(data, rep_delta = rep_percent - rep_center)
      fit_quadratic_start <- lmrob(log_length ~ rep_delta + I(rep_delta^2), data = start_data)
      start_gamma <- 0.03
      quadratic_coef <- coef(fit_quadratic_start)
      start_beta <- max(2 * quadratic_coef["I(rep_delta^2)"] / start_gamma, 1e-04)
      start_alpha <- quadratic_coef["rep_delta"] - start_beta
      list(log_L_ref = unname(quadratic_coef["(Intercept)"]), alpha = unname(start_alpha),
           log_beta = log(unname(start_beta)), log_gamma = log(start_gamma))
    }
    fit_ode_model <- function(data, start = NULL) {
      if (is.null(start)) start <- get_ode_start(data)
      nlrob(log_length ~ ode_log_length(rep_percent, log_L_ref, alpha, log_beta, log_gamma),
            data = data, start = start, method = "M", doCov = TRUE, model = TRUE,
            control = nls.control(maxiter = 500), trace = FALSE)
    }
    ode_predict_from_coef <- function(rep_percent, model_coef) {
      ode_log_length(rep_percent = rep_percent, log_L_ref = model_coef["log_L_ref"],
                     alpha = model_coef["alpha"], log_beta = model_coef["log_beta"],
                     log_gamma = model_coef["log_gamma"])
    }
    ode_gradient <- function(rep_percent, model_coef) {
      beta <- exp(model_coef["log_beta"])
      gamma <- exp(model_coef["log_gamma"])
      rep_delta <- rep_percent - rep_center
      exp_term <- exp(gamma * rep_delta)
      gradient <- cbind(log_L_ref = rep(1, length(rep_percent)), alpha = rep_delta,
                        log_beta = beta / gamma * expm1(gamma * rep_delta),
                        log_gamma = beta * (rep_delta * exp_term - expm1(gamma * rep_delta) / gamma))
      gradient[, names(model_coef), drop = FALSE]
    }
  })
  env
}
