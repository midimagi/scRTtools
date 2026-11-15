# Single-cell replication-state inference (two-state Poisson HMM).
# Migrated from rcode/processing/infer_single_cell_states.R. The three
# Rcpp::sourceCpp() calls of the original are replaced by the kernels compiled
# in src/ (single_cell_baum_welch_cpp, compute_baum_welch_expectations,
# viterbi_decode_cpp, findchangepts_cpp_max2). The algorithm is unchanged.

#' Infer replication states of one cell
#'
#' Infer per-window replication states (2 = unreplicated, 4 = replicated) of a
#' single cell from read counts in G1-derived windows, following the steps of
#' the original MATLAB-derived pipeline:
#'
#' 1. **CNV filtering**: counts are aggregated by `cnv_aggr_width` windows;
#'    two-change-point segmentation (linear statistic, BIC penalty) plus
#'    Poisson tests (BH-adjusted) iteratively mask segments whose mean differs
#'    from the background by at least `cnv_fc_cutoff`-fold with FDR below
#'    `cnv_fdr_cutoff`. Masked windows become `NaN` in `filtered_counts`.
#' 2. **Autocorrelation filter**: if the lag-20 autocorrelation of counts
#'    aggregated by `autocorr_aggr_width` exceeds `autocorr_threshold`, all
#'    states are `NA`.
#' 3. **Initialisation**: a two-component Poisson mixture on counts aggregated
#'    by `hmm_aggr_width`; if the ratio of the two means is outside
#'    `[init_ratio_lower, init_ratio_upper]`, all states are set to `2`.
#' 4. **Baum-Welch training** (up to `maxiter` iterations, convergence
#'    tolerance `tolerance`, switch probabilities at least `min_switch_prob`);
#'    a cell that does not converge gets all states `NA`.
#' 5. **Viterbi decoding**, restarting at `breakpoints` (chromosome starts).
#' 6. **Post-filter**: if the mean-count ratio between state-4 and state-2
#'    segments is outside `[final_ratio_lower, final_ratio_upper]`, all
#'    states are `NA`.
#'
#' @param counts Numeric vector of read counts per window (names are kept).
#' @param breakpoints Numeric vector of 1-based indices where a new chromosome
#'   starts.
#' @param autocorr_aggr_width Aggregation width (windows) for the
#'   autocorrelation filter.
#' @param autocorr_threshold Maximum lag-20 autocorrelation.
#' @param cnv_aggr_width Aggregation width for CNV filtering.
#' @param cnv_fc_cutoff Minimum fold change of a CNV segment (must be > 1).
#' @param cnv_fdr_cutoff FDR cutoff for CNV segments, in (0, 1).
#' @param min_segment_length Minimum segment length (aggregated windows) for
#'   CNV segmentation (must be >= 2).
#' @param hmm_aggr_width Aggregation width for HMM initialisation.
#' @param init_ratio_lower,init_ratio_upper Accepted range of the ratio of the
#'   two initial Poisson means.
#' @param min_switch_prob Minimum state-switch probability in `[0, 0.5)`.
#' @param tolerance Convergence tolerance of Baum-Welch.
#' @param maxiter Maximum number of Baum-Welch iterations.
#' @param final_ratio_lower,final_ratio_upper Accepted range of the final
#'   state-4 / state-2 read-depth ratio.
#'
#' @return A list with `filtered_counts` (numeric, `NaN` at missing or
#'   CNV-masked windows) and `states` (numeric, `2`, `4` or `NA`), both of
#'   length `length(counts)`.
#'
#' @seealso [assign_replication_states()] to run this over all cells.
#' @export
infer_single_cell_states <- function(
    counts,
    breakpoints,
    autocorr_aggr_width = 15L,
    autocorr_threshold = 0.20,
    cnv_aggr_width = 15L,  
    cnv_fc_cutoff = 1.5,
    cnv_fdr_cutoff = 0.01,
    min_segment_length = 2L,
    hmm_aggr_width = 10L, 
    init_ratio_lower = 1.5, 
    init_ratio_upper = 2.5,  
    min_switch_prob = 1e-4, 
    tolerance = 1e-6,  
    maxiter = 500L, 
    final_ratio_lower = 1.5, 
    final_ratio_upper = 2.2  
) {
  # ---------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------
  if (!is.numeric(counts) || !is.null(dim(counts))) {
    stop("'counts' must be a numeric vector.")
  }

  if (!is.numeric(breakpoints) || !is.null(dim(breakpoints))) {
    stop("'breakpoints' must be a numeric vector.")
  }

  if (!is.numeric(autocorr_aggr_width) ||
      length(autocorr_aggr_width) != 1L ||
      autocorr_aggr_width < 2) {
    stop("'autocorr_aggr_width' must be a single integer >= 2.")
  }

  if (!is.numeric(autocorr_threshold) ||
      length(autocorr_threshold) != 1L ||
      autocorr_threshold < 0 ||
      autocorr_threshold > 1) {
    stop("'autocorr_threshold' must be a single numeric value in [0, 1].")
  }

  if (!is.numeric(cnv_aggr_width) ||
      length(cnv_aggr_width) != 1L ||
      cnv_aggr_width < 1) {
    stop("'cnv_aggr_width' must be a single integer >= 1.")
  }

  if (!is.numeric(cnv_fc_cutoff) ||
      length(cnv_fc_cutoff) != 1L ||
      cnv_fc_cutoff <= 0) {
    stop("'cnv_fc_cutoff' must be a single positive numeric value.")
  }

  if (!is.numeric(cnv_fdr_cutoff) ||
      length(cnv_fdr_cutoff) != 1L ||
      cnv_fdr_cutoff < 0 ||
      cnv_fdr_cutoff > 1) {
    stop("'cnv_fdr_cutoff' must be a single numeric value in [0, 1].")
  }

  if (!is.numeric(min_segment_length) ||
      length(min_segment_length) != 1L ||
      min_segment_length < 1) {
    stop("'min_segment_length' must be a single integer >= 1.")
  }

  if (!is.numeric(hmm_aggr_width) ||
      length(hmm_aggr_width) != 1L ||
      hmm_aggr_width < 1) {
    stop("'hmm_aggr_width' must be a single integer >= 1.")
  }

  if (!is.numeric(init_ratio_lower) ||
      length(init_ratio_lower) != 1L ||
      init_ratio_lower <= 0) {
    stop("'init_ratio_lower' must be a single positive numeric value.")
  }

  if (!is.numeric(init_ratio_upper) ||
      length(init_ratio_upper) != 1L ||
      init_ratio_upper <= init_ratio_lower) {
    stop("'init_ratio_upper' must be a single numeric value greater than 'init_ratio_lower'.")
  }

  if (!is.numeric(min_switch_prob) ||
      length(min_switch_prob) != 1L ||
      min_switch_prob < 0 ||
      min_switch_prob >= 0.5) {
    stop("'min_switch_prob' must be a single numeric value in [0, 0.5).")
  }

  if (!is.numeric(tolerance) ||
      length(tolerance) != 1L ||
      tolerance <= 0) {
    stop("'tolerance' must be a single positive numeric value.")
  }

  if (!is.numeric(maxiter) ||
      length(maxiter) != 1L ||
      maxiter < 1) {
    stop("'maxiter' must be a single integer >= 1.")
  }

  if (!is.numeric(final_ratio_lower) ||
      length(final_ratio_lower) != 1L ||
      final_ratio_lower <= 0) {
    stop("'final_ratio_lower' must be a single positive numeric value.")
  }

  if (!is.numeric(final_ratio_upper) ||
      length(final_ratio_upper) != 1L ||
      final_ratio_upper <= final_ratio_lower) {
    stop("'final_ratio_upper' must be a single numeric value greater than 'final_ratio_lower'.")
  }

  # ---------------------------------------------------------------
  # Normalize argument types
  # ---------------------------------------------------------------
  counts_names <- names(counts)
  counts <- as.numeric(counts)
  names(counts) <- counts_names

  breakpoints <- as.numeric(breakpoints)

  autocorr_aggr_width <- as.integer(autocorr_aggr_width)
  cnv_aggr_width <- as.integer(cnv_aggr_width)
  min_segment_length <- as.integer(min_segment_length)
  hmm_aggr_width <- as.integer(hmm_aggr_width)
  maxiter <- as.integer(maxiter)

  # Initialize states vector
  states <- rep(NA_real_, length(counts))
  names(states) <- counts_names

  # ---------------------------------------------------------------
  # Step 1: CNV filtering
  # ---------------------------------------------------------------
  filtered_counts <- .filter_single_cell_cnv(
    counts,
    cnv_aggr_width = cnv_aggr_width,
    cnv_fc_cutoff = cnv_fc_cutoff,
    cnv_fdr_cutoff = cnv_fdr_cutoff,
    min_segment_length = min_segment_length
  )

  # ---------------------------------------------------------------
  # Step 2: Autocorrelation quality filter
  # ---------------------------------------------------------------
  non_nan_values <- filtered_counts[!is.na(filtered_counts)]

  if (length(non_nan_values) == 0L) {
    return(list(
      filtered_counts = filtered_counts,
      states = states
    ))
  }

  seq_agg <- .aggregate_counts(non_nan_values, autocorr_aggr_width)

  if (length(seq_agg) >= 21L) {
    autocorr_agg <- as.numeric(
      stats::acf(seq_agg, plot = FALSE, lag.max = 20L)$acf
    )

    if (autocorr_agg[21L] > autocorr_threshold) {
      return(list(
        filtered_counts = filtered_counts,
        states = states
      ))
    }
  }

  # ---------------------------------------------------------------
  # Step 3: HMM initialization
  # ---------------------------------------------------------------
  init <- .initialize_single_cell_hmm(
    count_data = filtered_counts,
    hmm_aggr_width = hmm_aggr_width,
    init_ratio_lower = init_ratio_lower,
    init_ratio_upper = init_ratio_upper
  )

  initial_emissions <- init$initial_emissions
  initial_transitions <- init$initial_transitions
  initial_states <- init$initial_states

  if (length(initial_states) == 0L) {
    states[] <- 2
    return(list(
      filtered_counts = filtered_counts,
      states = states
    ))
  }

  # ---------------------------------------------------------------
  # Step 4: HMM training
  # ---------------------------------------------------------------
  fit <- .train_single_cell_hmm(
    count_data = filtered_counts,
    initial_emissions = initial_emissions,
    initial_transitions = initial_transitions,
    initial_states = initial_states,
    min_switch_prob = min_switch_prob,
    tolerance = tolerance,
    maxiter = maxiter
  )

  emissions_matrix <- fit$emissions_matrix
  transition_matrix <- fit$transition_matrix
  converged <- isTRUE(fit$converged)

  if (!converged) {
    states[] <- NA_real_
    return(list(
      filtered_counts = filtered_counts,
      states = states
    ))
  }

  # ---------------------------------------------------------------
  # Step 5: Viterbi decoding
  # ---------------------------------------------------------------
  states <- .decode_single_cell_hmm(
    count_data = filtered_counts,
    transition_matrix = transition_matrix,
    emissions_matrix = emissions_matrix,
    initial_states = initial_states,
    breakpoints = breakpoints
  )

  names(states) <- counts_names

  # ---------------------------------------------------------------
  # Step 6: Post-filter based on count ratio between states
  # ---------------------------------------------------------------
  read_count_ratio <- .calculate_read_count_ratio(
    filtered_counts,
    states
  )

  if (is.na(read_count_ratio) ||
      read_count_ratio < final_ratio_lower ||
      read_count_ratio > final_ratio_upper) {
    states[] <- NA_real_
  }

  # ---------------------------------------------------------------
  # Return results
  # ---------------------------------------------------------------
  list(
    filtered_counts = filtered_counts,
    states = states
  )
}

#' Initial HMM parameters from a two-component Poisson mixture
#' @keywords internal
#' @noRd
.initialize_single_cell_hmm <- function(count_data,
                                       hmm_aggr_width = 10L,
                                       init_ratio_lower = 1.2,
                                       init_ratio_upper = 2.5) {
  # ---------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------
  if (!is.numeric(count_data) || !is.null(dim(count_data))) {
    stop("'count_data' must be a numeric vector.")
  }

  if (!is.numeric(hmm_aggr_width) ||
      length(hmm_aggr_width) != 1L ||
      hmm_aggr_width < 1) {
    stop("'hmm_aggr_width' must be a single integer >= 1.")
  }

  if (!is.numeric(init_ratio_lower) ||
      length(init_ratio_lower) != 1L ||
      init_ratio_lower <= 0) {
    stop("'init_ratio_lower' must be a single positive numeric value.")
  }

  if (!is.numeric(init_ratio_upper) ||
      length(init_ratio_upper) != 1L ||
      init_ratio_upper <= init_ratio_lower) {
    stop("'init_ratio_upper' must be a single numeric value greater than 'init_ratio_lower'.")
  }

  hmm_aggr_width <- as.integer(hmm_aggr_width)

  count_names <- names(count_data)

  # ---------------------------------------------------------------
  # Index NaNs / NAs
  # ---------------------------------------------------------------
  missing_values <- is.na(count_data)
  non_nan_values <- count_data[!missing_values]

  if (length(non_nan_values) == 0L) {
    return(list(
      initial_emissions = numeric(0),
      initial_transitions = numeric(0),
      initial_states = numeric(0)
    ))
  }

  # ---------------------------------------------------------------
  # Aggregate
  # ---------------------------------------------------------------
  seq <- .aggregate_counts(non_nan_values, hmm_aggr_width)

  if (length(seq) == 0L || all(is.na(seq))) {
    return(list(
      initial_emissions = numeric(0),
      initial_transitions = numeric(0),
      initial_states = numeric(0)
    ))
  }

  # ---------------------------------------------------------------
  # Fit two-component Poisson mixture
  # ---------------------------------------------------------------
  mu <- .fit_two_component_poisson_mixture(seq)

  if (length(mu) == 0L || any(!is.finite(mu))) {
    return(list(
      initial_emissions = numeric(0),
      initial_transitions = numeric(0),
      initial_states = numeric(0)
    ))
  }

  mu <- sort(mu)

  if (length(mu) < 2L || mu[1] <= 0) {
    return(list(
      initial_emissions = numeric(0),
      initial_transitions = numeric(0),
      initial_states = numeric(0)
    ))
  }

  mu_ratio <- mu[2] / mu[1]

  if (mu_ratio < init_ratio_lower || mu_ratio > init_ratio_upper) {
    return(list(
      initial_emissions = numeric(0),
      initial_transitions = numeric(0),
      initial_states = numeric(0)
    ))
  }

  # ---------------------------------------------------------------
  # Assign initial states on aggregated bins
  # ---------------------------------------------------------------
  probability_of_state <- matrix(NaN, nrow = length(seq), ncol = 2L)

  for (m in 1:2) {
    probability_of_state[, m] <- stats::dpois(seq, lambda = mu[m])
  }

  states <- rep(NaN, length(non_nan_values))
  max_probability <- apply(probability_of_state, 1L, max)

  for (state in seq_along(mu)) {
    index <- which(max_probability == probability_of_state[, state])

    if (length(index) > 0L) {
      starts <- (index - 1L) * hmm_aggr_width + 1L
      ends <- pmin(starts + hmm_aggr_width - 1L, length(non_nan_values))

      for (row in seq_along(starts)) {
        states[starts[row]:ends[row]] <- state
      }
    }
  }

  # ---------------------------------------------------------------
  # Restore states to original vector length
  # ---------------------------------------------------------------
  initial_states <- rep(NaN, length(missing_values))
  initial_states[!missing_values] <- states
  names(initial_states) <- count_names

  # ---------------------------------------------------------------
  # Initialize emission matrix
  # ---------------------------------------------------------------
  max_count <- as.integer(max(count_data, na.rm = TRUE))
  initial_emissions <- matrix(NaN, nrow = 2L, ncol = max_count + 1L)

  for (state in 1:2) {
    state_mask <- !is.na(initial_states) & initial_states == state

    if (!any(state_mask)) {
      lambda_state <- mean(count_data, na.rm = TRUE)
    } else {
      lambda_state <- mean(count_data[state_mask], na.rm = TRUE)
    }

    if (!is.finite(lambda_state) || lambda_state < 0) {
      lambda_state <- mean(count_data, na.rm = TRUE)
    }

    initial_emissions[state, ] <- stats::dpois(
      0:max_count,
      lambda = lambda_state
    )
  }

  # ---------------------------------------------------------------
  # Initialize transition matrix
  # ---------------------------------------------------------------
  initial_transitions <- matrix(0, nrow = 2L, ncol = 2L)

  if (length(initial_states) >= 2L) {
    from_states <- initial_states[-length(initial_states)]
    to_states <- initial_states[-1L]

    for (state1 in 1:2) {
      for (state2 in 1:2) {
        initial_transitions[state1, state2] <- sum(
          from_states == state1 & to_states == state2,
          na.rm = TRUE
        )
      }
    }
  }

  # Normalize transition rows
  row_sums <- rowSums(initial_transitions)

  zero_rows <- row_sums == 0 | !is.finite(row_sums)

  if (any(zero_rows)) {
    idx <- which(zero_rows)
    initial_transitions[idx, ] <- 0
    initial_transitions[cbind(idx, idx)] <- 1
    row_sums <- rowSums(initial_transitions)
  }

  initial_transitions <- sweep(
    initial_transitions,
    1L,
    row_sums,
    "/"
  )

  initial_transitions[is.na(initial_transitions)] <- 0

  # ---------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------
  list(
    initial_emissions = initial_emissions,
    initial_transitions = initial_transitions,
    initial_states = initial_states
  )
}

#' Baum-Welch training (E-step in C++)
#' @keywords internal
#' @noRd
.train_single_cell_hmm <- function(count_data,
                                  initial_emissions,
                                  initial_transitions,
                                  initial_states,
                                  min_switch_prob = 0,
                                  tolerance = 1e-6,
                                  maxiter = 500L) {

  # ---------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------
  if (!is.numeric(count_data) || !is.null(dim(count_data))) {
    stop("'count_data' must be a numeric vector.")
  }

  observed_counts <- count_data[!is.na(count_data)]
  if (length(observed_counts) == 0L) {
    stop("'count_data' must contain at least one non-missing value.")
  }

  if (any(observed_counts < 0) || any(observed_counts != floor(observed_counts))) {
    stop("'count_data' must contain non-negative integer-like values.")
  }

  if (!is.numeric(min_switch_prob) ||
      length(min_switch_prob) != 1L ||
      min_switch_prob < 0 ||
      min_switch_prob >= 0.5) {
    stop("'min_switch_prob' must be a single numeric value in [0, 0.5).")
  }

  if (!is.numeric(tolerance) ||
      length(tolerance) != 1L ||
      tolerance <= 0) {
    stop("'tolerance' must be a single positive numeric value.")
  }

  if (!is.numeric(maxiter) ||
      length(maxiter) != 1L ||
      maxiter < 1) {
    stop("'maxiter' must be a single integer >= 1.")
  }

  maxiter <- as.integer(maxiter)

  num_states <- 2L
  num_emissions <- as.integer(max(count_data, na.rm = TRUE)) + 1L

  if (!is.matrix(initial_transitions) ||
      !all(dim(initial_transitions) == c(num_states, num_states))) {
    stop("'initial_transitions' must be a 2 x 2 matrix.")
  }

  if (!is.matrix(initial_emissions) ||
      nrow(initial_emissions) != num_states ||
      ncol(initial_emissions) != num_emissions) {
    stop("'initial_emissions' has incompatible dimensions.")
  }

  # ---------------------------------------------------------------
  # Compute initial state probabilities from 'initial_states'
  # ---------------------------------------------------------------
  valid_states <- initial_states[!is.na(initial_states)]

  if (length(valid_states) > 0) {
    init_probs <- c(
      sum(valid_states == 1) / length(valid_states),
      sum(valid_states == 2) / length(valid_states)
    )
  } else {
    init_probs <- c(0.5, 0.5)
  }

  # Fallback if init_probs become invalid
  if (any(!is.finite(init_probs)) || sum(init_probs) <= 0) {
    init_probs <- c(0.5, 0.5)
  } else {
    init_probs <- init_probs / sum(init_probs)
  }

  # ---------------------------------------------------------------
  # Initialize parameters
  # ---------------------------------------------------------------
  transition_matrix <- initial_transitions
  emissions_matrix <- initial_emissions

  transition_matrix[is.na(transition_matrix)] <- 0
  emissions_matrix[is.na(emissions_matrix)] <- 0

  transition_matrix <- sweep(
    transition_matrix,
    1,
    rowSums(transition_matrix),
    "/"
  )

  emissions_matrix <- sweep(
    emissions_matrix,
    1,
    rowSums(emissions_matrix),
    "/"
  )

  transition_matrix[is.na(transition_matrix)] <- 0
  emissions_matrix[is.na(emissions_matrix)] <- 0

  # Apply minimum switch probability to the initial transition matrix
  if (min_switch_prob > 0) {
    transition_matrix[1, 2] <- max(transition_matrix[1, 2], min_switch_prob)
    transition_matrix[2, 1] <- max(transition_matrix[2, 1], min_switch_prob)

    transition_matrix[1, 1] <- 1 - transition_matrix[1, 2]
    transition_matrix[2, 2] <- 1 - transition_matrix[2, 1]
  }

  converged <- FALSE
  log_likelihood <- 1

  # ---------------------------------------------------------------
  # Baum-Welch iterations
  # ---------------------------------------------------------------
  for (iteration in seq_len(maxiter)) {

    previous_iter <- list(
      log_likelihood = log_likelihood,
      emissions = emissions_matrix,
      transitions = transition_matrix
    )

    # ---------------------------------------------------------------
    # E-step: forward-backward algorithm
    # ---------------------------------------------------------------
    bw <- single_cell_baum_welch_cpp(
      count_data = count_data,
      transition_matrix = transition_matrix,
      emissions_matrix = emissions_matrix,
      init_probs = init_probs
    )

    log_likelihood <- bw$log_probability_seq

    # ---------------------------------------------------------------
    # M-step: compute expected counts
    # ---------------------------------------------------------------
    exp_obj <- compute_baum_welch_expectations(
      count_data = count_data,
      alpha = bw$alpha,
      beta = bw$beta,
      scalar = bw$scalar,
      transition_matrix = transition_matrix,
      emissions_matrix = emissions_matrix
    )

    transition_matrix <- exp_obj$trans_exp
    emissions_matrix  <- exp_obj$emiss_exp

    # ---------------------------------------------------------------
    # Normalize rows to obtain probabilities
    # ---------------------------------------------------------------
    transition_row_sums <- rowSums(transition_matrix)
    emission_row_sums   <- rowSums(emissions_matrix)

    zero_transition_rows <- transition_row_sums == 0 | !is.finite(transition_row_sums)
    zero_emission_rows   <- emission_row_sums == 0 | !is.finite(emission_row_sums)

    if (any(zero_transition_rows)) {
      idx <- which(zero_transition_rows)
      transition_matrix[idx, ] <- 0
      transition_matrix[cbind(idx, idx)] <- 1
      transition_row_sums <- rowSums(transition_matrix)
    }

    if (any(zero_emission_rows)) {
      idx <- which(zero_emission_rows)
      emissions_matrix[idx, ] <- 1 / num_emissions
      emission_row_sums <- rowSums(emissions_matrix)
    }

    transition_matrix <- sweep(transition_matrix, 1, transition_row_sums, "/")
    emissions_matrix  <- sweep(emissions_matrix,  1, emission_row_sums,  "/")

    transition_matrix[is.na(transition_matrix)] <- 0
    emissions_matrix[is.na(emissions_matrix)] <- 0

    # ---------------------------------------------------------------
    # Enforce minimum switch probability
    # ---------------------------------------------------------------
    if (min_switch_prob > 0) {
      transition_matrix[1, 2] <- max(transition_matrix[1, 2], min_switch_prob)
      transition_matrix[2, 1] <- max(transition_matrix[2, 1], min_switch_prob)

      transition_matrix[1, 1] <- 1 - transition_matrix[1, 2]
      transition_matrix[2, 2] <- 1 - transition_matrix[2, 1]
    }

    # ---------------------------------------------------------------
    # Check for illegal values after M-step
    # ---------------------------------------------------------------
    if (any(!is.finite(transition_matrix)) ||
        any(!is.finite(emissions_matrix))) {
      converged <- FALSE
      break
    }

    # ---------------------------------------------------------------
    # Check convergence
    # ---------------------------------------------------------------
    likelihood_converged <- FALSE

    if (is.finite(log_likelihood) &&
        is.finite(previous_iter$log_likelihood)) {
      rel_change <- abs(log_likelihood - previous_iter$log_likelihood) /
        (1 + abs(previous_iter$log_likelihood))

      likelihood_converged <- rel_change < tolerance
    }

    trans_diff <- tryCatch(
      norm(transition_matrix - previous_iter$transitions, type = "I"),
      error = function(e) NA_real_
    )

    transitions_converged <- if (is.finite(trans_diff)) {
      (trans_diff / num_states) < tolerance
    } else {
      FALSE
    }

    emiss_diff <- tryCatch(
      norm(emissions_matrix - previous_iter$emissions, type = "I"),
      error = function(e) NA_real_
    )

    emissions_converged <- if (is.finite(emiss_diff)) {
      (emiss_diff / num_emissions) < tolerance
    } else {
      FALSE
    }

    if (likelihood_converged &&
        transitions_converged &&
        emissions_converged) {
      converged <- TRUE
      break
    }
  }

  list(
    emissions_matrix = emissions_matrix,
    transition_matrix = transition_matrix,
    converged = converged
  )
}

#' Viterbi decoding (C++), restarting at breakpoints
#' @keywords internal
#' @noRd
.decode_single_cell_hmm <- function(count_data,
                                   transition_matrix,
                                   emissions_matrix,
                                   initial_states,
                                   breakpoints) {
  # ---------------------------------------------------------------
  # Input validation (unchanged from original)
  # ---------------------------------------------------------------
  if (!is.numeric(count_data) || !is.null(dim(count_data))) {
    stop("'count_data' must be a numeric vector.")
  }

  observed_idx <- which(!is.na(count_data))
  if (length(observed_idx) == 0L) {
    stop("'count_data' must contain at least one non-missing value.")
  }

  num_states <- 2L
  num_emissions <- as.integer(max(count_data, na.rm = TRUE)) + 1L

  if (!is.matrix(transition_matrix) ||
      !all(dim(transition_matrix) == c(num_states, num_states))) {
    stop("'transition_matrix' must be a 2 x 2 matrix.")
  }

  if (!is.matrix(emissions_matrix) ||
      nrow(emissions_matrix) != num_states ||
      ncol(emissions_matrix) != num_emissions) {
    stop("'emissions_matrix' has incompatible dimensions.")
  }

  # Ensure breakpoints are integers within range
  if (!is.numeric(breakpoints) && !is.integer(breakpoints)) {
    stop("'breakpoints' must be a numeric vector of indices.")
  }
  breakpoints <- as.integer(breakpoints)
  breakpoints <- breakpoints[breakpoints >= 1 & breakpoints <= length(count_data)]

  # ---------------------------------------------------------------
  # Call Rcpp core
  # ---------------------------------------------------------------
  path <- viterbi_decode_cpp(
    count_data = count_data,
    transition_matrix = transition_matrix,
    emissions_matrix = emissions_matrix,
    initial_states = initial_states,
    breakpoints = breakpoints
  )

  # Return path (already contains 2/4 states)
  return(path)
}

#' Mask CNV segments; returns a vector of length(counts)
#' @keywords internal
#' @noRd
.filter_single_cell_cnv <- function(counts,
                                   cnv_aggr_width = 15L,
                                   cnv_fc_cutoff = 1.5,
                                   cnv_fdr_cutoff = 0.01,
                                   min_segment_length = 2L) {

  # ---------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------
  if (!is.numeric(counts) || !is.null(dim(counts))) {
    stop("'counts' must be a numeric vector.")
  }

  if (!is.numeric(cnv_aggr_width) ||
      length(cnv_aggr_width) != 1L ||
      cnv_aggr_width < 1) {
    stop("'cnv_aggr_width' must be a single integer >= 1.")
  }

  if (!is.numeric(cnv_fc_cutoff) ||
      length(cnv_fc_cutoff) != 1L ||
      cnv_fc_cutoff <= 1) {
    stop("'cnv_fc_cutoff' must be a single numeric value > 1.")
  }

  if (!is.numeric(cnv_fdr_cutoff) ||
      length(cnv_fdr_cutoff) != 1L ||
      cnv_fdr_cutoff <= 0 ||
      cnv_fdr_cutoff >= 1) {
    stop("'cnv_fdr_cutoff' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(min_segment_length) ||
      length(min_segment_length) != 1L ||
      min_segment_length < 1) {
    stop("'min_segment_length' must be a single integer >= 1.")
  }

  cnv_aggr_width <- as.integer(cnv_aggr_width)
  min_segment_length <- as.integer(min_segment_length)

  # ---------------------------------------------------------------
  # Index NaNs / NAs
  # ---------------------------------------------------------------
  missing_values <- is.na(counts)
  non_nan_values <- counts[!missing_values]

  if (length(non_nan_values) == 0L) {
    return(rep(NaN, length(counts)))
  }

  # ---------------------------------------------------------------
  # Aggregate
  # ---------------------------------------------------------------
  seq <- .aggregate_counts(non_nan_values, cnv_aggr_width)

  if (length(seq) == 0L) {
    filtered_counts <- rep(NaN, length(missing_values))
    filtered_counts[!missing_values] <- non_nan_values
    return(filtered_counts)
  }

  # ---------------------------------------------------------------
  # Identify CNVs via Poisson tests
  # ---------------------------------------------------------------
  is_not_cnv <- rep(TRUE, length(seq))

  repeat {
    index <- which(is_not_cnv)

    if (length(index) < 4L) {
      break
    }

    counts_to_filter <- seq[index]
    counts_to_filter <- .fillmissing_nearest(
      counts_to_filter,
      endvalues = "nearest"
    )

    if (anyNA(counts_to_filter) || any(!is.finite(counts_to_filter))) {
      break
    }

    changepoints <- .findchangepts_2cpts(
      counts_to_filter,
      MaxNumChanges = 2L,
      Statistic = "linear",
      MinSegmentLength = min_segment_length,
      Penalty = "BIC"
    )

    if (length(changepoints) == 0L) {
      break
    }

    # -------------------------------------------------------------
    # Convert change points to non-overlapping regions
    # -------------------------------------------------------------
    n_counts <- length(counts_to_filter)

    changepoints <- sort(unique(as.integer(changepoints)))
    changepoints <- changepoints[changepoints >= 1L & changepoints < n_counts]

    region_start <- c(1L, changepoints + 1L)
    region_end   <- c(changepoints, n_counts)

    regions <- cbind(region_start, region_end)

    if (nrow(regions) == 0L) {
      break
    }

    csum <- c(0, cumsum(counts_to_filter))
    total_sum <- csum[length(csum)]

    log2_fc_cutoff <- log2(cnv_fc_cutoff)

    region_stats <- vector("list", nrow(regions))
    keep_region <- logical(nrow(regions))

    for (k in seq_len(nrow(regions))) {
      s <- regions[k, 1]
      e <- regions[k, 2]

      seg_len <- e - s + 1L
      bg_len <- n_counts - seg_len

      if (seg_len < min_segment_length || bg_len < min_segment_length) {
        next
      }

      seg_sum <- csum[e + 1L] - csum[s]
      bg_sum <- total_sum - seg_sum

      mu_seg <- seg_sum / seg_len
      mu_bg  <- bg_sum / bg_len

      if (!is.finite(mu_seg) || !is.finite(mu_bg) || mu_bg <= 0) {
        next
      }

      ratio <- mu_seg / mu_bg
      log2_ratio <- log2(ratio)

      if (!is.finite(log2_ratio) || abs(log2_ratio) < log2_fc_cutoff) {
        next
      }

      pt <- stats::poisson.test(
        x = c(seg_sum, bg_sum),
        T = c(seg_len, bg_len),
        alternative = "two.sided"
      )

      keep_region[k] <- TRUE
      region_stats[[k]] <- data.frame(
        start = s,
        end = e,
        width = seg_len,
        log2_ratio = log2_ratio,
        p_value = pt$p.value,
        stringsAsFactors = FALSE
      )
    }

    if (!any(keep_region)) {
      break
    }

    region_stats <- do.call(rbind, region_stats[keep_region])

    region_stats$p_adj <- stats::p.adjust(
      region_stats$p_value,
      method = "BH"
    )

    sig_regions <- region_stats[
      region_stats$p_adj < cnv_fdr_cutoff,
      ,
      drop = FALSE
    ]

    if (nrow(sig_regions) == 0L) {
      break
    }

    # -------------------------------------------------------------
    # Select the strongest candidate region:
    # smallest adjusted p-value,
    # then largest absolute effect size,
    # then longer width
    # -------------------------------------------------------------
    best_idx <- order(
      sig_regions$p_adj,
      -abs(sig_regions$log2_ratio),
      -sig_regions$width
    )[1L]

    best_region <- sig_regions[best_idx, , drop = FALSE]

    # -------------------------------------------------------------
    # Map back to original aggregated indices and mark as CNV
    # -------------------------------------------------------------
    cnv_region_filtered <- best_region$start:best_region$end
    cnv_region_original <- index[cnv_region_filtered]

    is_not_cnv[cnv_region_original] <- FALSE
  }

  cnv_windows <- which(!is_not_cnv)

  # ---------------------------------------------------------------
  # Disaggregate CNV windows back to original non-NA positions
  # ---------------------------------------------------------------
  if (length(cnv_windows) > 0L) {
    cnv_starts <- (cnv_windows - 1L) * cnv_aggr_width + 1L
    cnv_ends <- pmin(
      cnv_starts + cnv_aggr_width - 1L,
      length(non_nan_values)
    )

    cnv_windows <- cbind(cnv_starts, cnv_ends)
  } else {
    cnv_windows <- matrix(integer(0), ncol = 2L)
  }

  filtered_non_nan_counts <- non_nan_values

  if (nrow(cnv_windows) > 0L) {
    for (row in seq_len(nrow(cnv_windows))) {
      filtered_non_nan_counts[
        cnv_windows[row, 1]:cnv_windows[row, 2]
      ] <- NaN
    }
  }

  # ---------------------------------------------------------------
  # Add back masked windows
  # ---------------------------------------------------------------
  filtered_counts <- rep(NaN, length(missing_values))
  filtered_counts[!missing_values] <- filtered_non_nan_counts
  names(filtered_counts) <- names(counts)
  filtered_counts
}

#' Sum consecutive blocks of rows (NA removed)
#' @keywords internal
#' @noRd
.aggregate_counts <- function(counts, windows_to_aggregate) {
  if (length(windows_to_aggregate) != 1L ||
      !is.numeric(windows_to_aggregate) ||
      is.na(windows_to_aggregate) ||
      windows_to_aggregate <= 0) {
    stop("'windows_to_aggregate' must be one positive integer.")
  }

  windows_to_aggregate <- as.integer(windows_to_aggregate)

  # make sure counts is a matrix
  if (is.vector(counts)) {
    counts <- matrix(counts, ncol = 1L)
  } else {
    counts <- as.matrix(counts)
  }

  n_rows <- nrow(counts)
  n_cols <- ncol(counts)

  aggregated_height <- floor((n_rows - 1L) / windows_to_aggregate) + 1L
  aggregated <- matrix(0, nrow = aggregated_height, ncol = n_cols)

  index <- seq.int(1L, n_rows, by = windows_to_aggregate)

  for (row in seq_along(index)) {
    low <- index[row]
    high <- min(low + windows_to_aggregate - 1L, n_rows)
    aggregated[row, ] <- colSums(counts[low:high, , drop = FALSE], na.rm = TRUE)
  }

  colnames(aggregated) <- colnames(counts)
  aggregated
}

#' Fill NA with the nearest value (MATLAB fillmissing "nearest")
#' @keywords internal
#' @noRd
.fillmissing_nearest <- function(x, endvalues = NA) {
  if (!is.vector(x)) stop("x must be a vector")
  
  # Indices and values of non-missing elements
  idx_ok <- which(!is.na(x))
  val_ok <- x[idx_ok]
  n <- length(x)
  
  # If all are NA, return unchanged (or handle endvalues? but no valid values)
  if (length(idx_ok) == 0) return(x)
  
  result <- x
  
  for (i in seq_len(n)) {
    if (is.na(x[i])) {
      # Find left and right nearest indices
      left_idx <- max(idx_ok[idx_ok < i])
      right_idx <- min(idx_ok[idx_ok > i])
      
      # Both sides exist
      if (!is.infinite(left_idx) && !is.infinite(right_idx)) {
        left_dist <- i - left_idx
        right_dist <- right_idx - i
        if (left_dist < right_dist) {
          result[i] <- x[left_idx]
        } else if (right_dist < left_dist) {
          result[i] <- x[right_idx]
        } else {  # equal distances -> take right (MATLAB default)
          result[i] <- x[right_idx]
        }
      }
      # Only left exists (end on the right side)
      else if (!is.infinite(left_idx) && is.infinite(right_idx)) {
        if (!is.na(endvalues) && endvalues == "nearest") {
          result[i] <- x[left_idx]   # nearest is the only side
        } else {
          result[i] <- NA  # keep NA
        }
      }
      # Only right exists (end on the left side)
      else if (is.infinite(left_idx) && !is.infinite(right_idx)) {
        if (!is.na(endvalues) && endvalues == "nearest") {
          result[i] <- x[right_idx]  # nearest is the only side
        } else {
          result[i] <- NA
        }
      }
      # No valid values (should not happen because we already checked length>0)
      else {
        result[i] <- NA
      }
    }
  }
  return(result)
}

#' Sorted means of a two-component Poisson mixture (L-BFGS-B)
#' @keywords internal
#' @noRd
.fit_two_component_poisson_mixture <- function(seq) {
  seq <- as.numeric(seq)
  seq <- seq[is.finite(seq)]

  if (length(seq) == 0L) {
    return(numeric(0))
  }

  mixed_poisson <- function(seq, proportion1, logmu1, logmu2) {
    proportion1 * stats::dpois(seq, lambda = exp(logmu1)) +
      (1 - proportion1) * stats::dpois(seq, lambda = exp(logmu2))
  }

  med <- stats::median(seq, na.rm = TRUE)
  initial_guess <- log(c(med, 2 * med))

  # MATLAB equivalent:
  # if any(~isfinite(initial_guess)) || any(initial_guess == 0)
  if (any(!is.finite(initial_guess)) || any(initial_guess == 0)) {
    return(numeric(0))
  }

  neg_loglik <- function(par) {
    proportion1 <- par[1]
    logmu1 <- par[2]
    logmu2 <- par[3]

    prob <- mixed_poisson(seq, proportion1, logmu1, logmu2)
    prob[!is.finite(prob) | prob <= 0] <- .Machine$double.xmin

    -sum(log(prob))
  }

  fit <- try(
    stats::optim(
      par = c(0.5, initial_guess),
      fn = neg_loglik,
      method = "L-BFGS-B",
      lower = c(0, 0, 0),
      upper = c(1, Inf, Inf),
      control = list(maxit = 1000)
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    return(numeric(0))
  }

  mu <- sort(exp(fit$par[2:3]), decreasing = FALSE)
  mu
}

#' At most two change points (linear statistic), C++ backend
#' @keywords internal
#' @noRd
.findchangepts_2cpts <- function(
  x,
  MaxNumChanges = 2L,
  Statistic = "linear",
  MinSegmentLength = 2L,
  Penalty = c("BIC", "AIC", "none")
) {
  Penalty <- match.arg(Penalty)

  # Input validation
  if (!is.numeric(x) || !is.null(dim(x))) {
    stop("'x' must be a numeric vector.")
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop("'x' cannot contain NA/NaN/Inf.")
  }
  if (!identical(Statistic, "linear")) {
    stop("Only Statistic = 'linear' is currently supported.")
  }
  if (length(MaxNumChanges) != 1L || is.na(MaxNumChanges)) {
    stop("'MaxNumChanges' must be a single non-missing integer.")
  }
  if (length(MinSegmentLength) != 1L || is.na(MinSegmentLength)) {
    stop("'MinSegmentLength' must be a single non-missing integer.")
  }

  x <- as.numeric(x)
  n <- length(x)
  MaxNumChanges <- as.integer(MaxNumChanges)
  MinSegmentLength <- as.integer(MinSegmentLength)

  if (MaxNumChanges != 2L) {
    stop("This specialized fast version only supports MaxNumChanges = 2L.")
  }
  if (MinSegmentLength < 2L) {
    stop("'MinSegmentLength' must be >= 2.")
  }
  if (n < 2L * MinSegmentLength) {
    return(integer(0))
  }

  penalty_code <- switch(
    Penalty,
    "none" = 0L,
    "AIC"  = 1L,
    "BIC"  = 2L
  )

  cps <- findchangepts_cpp_max2(
    x = x,
    MinSegmentLength = MinSegmentLength,
    penalty_code = penalty_code
  )

  sort(as.integer(cps))
}

#' Mean read depth of state-4 over state-2 segments
#' @keywords internal
#' @noRd
.calculate_read_count_ratio <- function(counts, states) {
  if (!is.numeric(counts) || !is.null(dim(counts))) {
    stop("'counts' must be a numeric vector.")
  }

  if (!is.numeric(states) || !is.null(dim(states))) {
    stop("'states' must be a numeric vector.")
  }

  if (length(counts) != length(states)) {
    stop("'counts' and 'states' must have the same length.")
  }

  if (all(is.na(states))) {
    return(NaN)
  }

  non_na_idx <- which(!is.na(states))
  non_na_states <- states[non_na_idx]

  if (all(non_na_states == 2) || all(non_na_states == 4)) {
    return(NaN)
  }

  # Identify segments of consistent state
  is_breakpoint <- abs(diff(non_na_states)) == 2
  segment_starts <- c(1L, non_na_idx[which(is_breakpoint)] + 1L)
  segment_ends <- c(segment_starts[-1L] - 1L, length(states))

  copy_number <- rep(NaN, length(segment_starts))
  read_depth <- rep(NaN, length(segment_starts))

  # Calculate mean read depth in each segment
  for (i in seq_along(segment_starts)) {
    start <- segment_starts[i]
    end <- segment_ends[i]

    segment_states <- states[start:end]
    is_assigned_state <- !is.na(segment_states)

    segment_cn <- unique(segment_states[is_assigned_state])
    if (length(segment_cn) != 1L) {
      stop("Each segment must contain exactly one assigned state.")
    }

    copy_number[i] <- segment_cn

    segment_counts <- counts[start:end]
    segment_counts <- segment_counts[is_assigned_state]
    read_depth[i] <- mean(segment_counts, na.rm = TRUE)
  }

  mean(read_depth[copy_number == 4], na.rm = TRUE) /
    mean(read_depth[copy_number == 2], na.rm = TRUE)
}
