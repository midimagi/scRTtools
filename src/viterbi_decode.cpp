#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector viterbi_decode_cpp(NumericVector count_data,
                                 NumericMatrix transition_matrix,
                                 NumericMatrix emissions_matrix,
                                 NumericVector initial_states,
                                 IntegerVector breakpoints) {
  int n = count_data.size();
  int num_states = 2;

  // ---------------------------------------------------------------
  // Precompute log matrices with a tiny offset to avoid log(0)
  // ---------------------------------------------------------------
  const double MIN_PROB = 1e-300;
  NumericMatrix logA(num_states, num_states);
  for (int i = 0; i < num_states; ++i) {
    for (int j = 0; j < num_states; ++j) {
      double val = transition_matrix(i, j);
      if (val <= 0) val = MIN_PROB;
      logA(i, j) = log(val);
    }
  }

  int num_emissions = emissions_matrix.ncol();
  NumericMatrix logE(num_states, num_emissions);
  for (int i = 0; i < num_states; ++i) {
    for (int j = 0; j < num_emissions; ++j) {
      double val = emissions_matrix(i, j);
      if (val <= 0) val = MIN_PROB;
      logE(i, j) = log(val);
    }
  }

  // ---------------------------------------------------------------
  // Compute initial state probabilities from 'initial_states'
  // ---------------------------------------------------------------
  double init_prob[2] = {0.0, 0.0};
  int valid_count = 0;
  for (int i = 0; i < n; ++i) {
    if (!NumericVector::is_na(initial_states[i])) {
      int s = static_cast<int>(initial_states[i]);
      if (s >= 1 && s <= num_states) {
        init_prob[s - 1] += 1.0;
        valid_count++;
      }
    }
  }
  if (valid_count > 0) {
    init_prob[0] /= valid_count;
    init_prob[1] /= valid_count;
  } else {
    init_prob[0] = 0.5;
    init_prob[1] = 0.5;
  }

  // ---------------------------------------------------------------
  // Precompute next valid index for each position (for backtracking)
  // ---------------------------------------------------------------
  IntegerVector next_valid(n, NA_INTEGER);
  int last_valid = -1;
  for (int i = n - 1; i >= 0; --i) {
    if (!NumericVector::is_na(count_data[i])) {
      last_valid = i;
    }
    next_valid[i] = last_valid;
  }

  // ---------------------------------------------------------------
  // Find first valid position
  // ---------------------------------------------------------------
  int first_obs = -1;
  for (int i = 0; i < n; ++i) {
    if (!NumericVector::is_na(count_data[i])) {
      first_obs = i;
      break;
    }
  }
  if (first_obs == -1) {
    return NumericVector::create();  // should never reach here due to R validation
  }

  // ---------------------------------------------------------------
  // Viterbi forward recursion
  // ---------------------------------------------------------------
  // V matrix: current log probabilities for each state
  NumericMatrix V(num_states, n);
  std::fill(V.begin(), V.end(), R_NegInf);
  // Ptr matrix: store the best previous state for backtracking
  IntegerMatrix Ptr(num_states, n);
  std::fill(Ptr.begin(), Ptr.end(), NA_INTEGER);

  // Initialize at first valid position
  int obs = static_cast<int>(count_data[first_obs]);
  if (obs < 0 || obs >= num_emissions) {
    stop("Observation value out of range of emissions matrix.");
  }
  for (int s = 0; s < num_states; ++s) {
    V(s, first_obs) = log(init_prob[s]) + logE(s, obs);
  }

  // Track prior V (for breakpoint restarts, we need the V of previous *valid* position)
  // We will compute V for all valid positions sequentially.
  int prev_valid_pos = first_obs;
  NumericVector prior_V(2);
  prior_V[0] = V(0, first_obs);
  prior_V[1] = V(1, first_obs);

  // Helper lambda to reset prior_V at a breakpoint
  auto reset_prior_at_breakpoint = [&](int pos) {
    int obs_bp = static_cast<int>(count_data[pos]);
    if (obs_bp < 0 || obs_bp >= num_emissions) {
      stop("Observation value out of range of emissions matrix at breakpoint.");
    }
    for (int s = 0; s < num_states; ++s) {
      prior_V[s] = log(init_prob[s]) + logE(s, obs_bp);
    }
  };

  // Iterate over all positions
  for (int pos = first_obs + 1; pos < n; ++pos) {
    if (NumericVector::is_na(count_data[pos])) {
      continue;  // skip NA positions, they remain -Inf in V
    }

    // Check if this position is a breakpoint
    bool is_break = false;
    for (int bp : breakpoints) {
      if (bp == pos + 1) {  // R uses 1-based indexing for breakpoints
        is_break = true;
        break;
      }
    }

    // If breakpoint, reset prior_V using initial probabilities and current emission
    if (is_break) {
      reset_prior_at_breakpoint(pos);
    }

    // Current observation
    obs = static_cast<int>(count_data[pos]);
    if (obs < 0 || obs >= num_emissions) {
      stop("Observation value out of range of emissions matrix.");
    }

    // Compute Viterbi step for both states
    for (int curr = 0; curr < num_states; ++curr) {
      double max_prob = R_NegInf;
      int best_prev = NA_INTEGER;
      for (int prev = 0; prev < num_states; ++prev) {
        double prob = prior_V[prev] + logA(prev, curr);
        if (prob > max_prob) {
          max_prob = prob;
          best_prev = prev;
        }
      }
      V(curr, pos) = max_prob + logE(curr, obs);
      Ptr(curr, pos) = best_prev + 1;  // store 1-based state index for R compatibility
    }

    // Update prior_V for the next valid position
    prior_V[0] = V(0, pos);
    prior_V[1] = V(1, pos);
    prev_valid_pos = pos;
  }

  // ---------------------------------------------------------------
  // Backtracking
  // ---------------------------------------------------------------
  NumericVector path(n, NA_REAL);

  // Find last valid position
  int last_obs = -1;
  for (int i = n - 1; i >= 0; --i) {
    if (!NumericVector::is_na(count_data[i])) {
      last_obs = i;
      break;
    }
  }

  // Determine best final state
  int best_final = (V(0, last_obs) >= V(1, last_obs)) ? 0 : 1;
  path[last_obs] = best_final + 1;  // store as 1 or 2

  // Backtrack from last_obs backwards
  int curr_pos = last_obs;
  int curr_state = best_final;

  while (curr_pos > first_obs) {
    // Find previous valid position
    int prev_pos = -1;
    for (int i = curr_pos - 1; i >= 0; --i) {
      if (!NumericVector::is_na(count_data[i])) {
        prev_pos = i;
        break;
      }
    }
    if (prev_pos == -1) break;

    // Retrieve best previous state from Ptr
    int prev_state = Ptr(curr_state, curr_pos) - 1;  // convert back to 0-based
    if (prev_state < 0 || prev_state >= num_states) {
      // Degenerate case: breakpoint? For robustness, if Ptr invalid, choose max V
      double v0 = V(0, prev_pos);
      double v1 = V(1, prev_pos);
      prev_state = (v0 >= v1) ? 0 : 1;
    }

    path[prev_pos] = prev_state + 1;

    curr_pos = prev_pos;
    curr_state = prev_state;
  }

  // ---------------------------------------------------------------
  // Convert state labels to copy-number (1->2, 2->4)
  // ---------------------------------------------------------------
  for (int i = 0; i < n; ++i) {
    if (!NumericVector::is_na(path[i])) {
      if (path[i] == 1) path[i] = 2;
      else if (path[i] == 2) path[i] = 4;
    }
  }

  return path;
}