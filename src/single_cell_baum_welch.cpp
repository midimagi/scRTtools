#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List single_cell_baum_welch_cpp(NumericVector count_data,
                                NumericMatrix transition_matrix,
                                NumericMatrix emissions_matrix,
                                NumericVector init_probs) {
  // ---------------------------------------------------------------
  // Input validation
  // ---------------------------------------------------------------
  int n = count_data.size();
  int num_states = transition_matrix.nrow();

  if (num_states != 2) {
    stop("This implementation is optimized for exactly 2 states.");
  }
  if (emissions_matrix.nrow() != 2) {
    stop("'emissions_matrix' must have 2 rows.");
  }
  if (init_probs.size() != 2) {
    stop("'init_probs' must be a numeric vector of length 2.");
  }

  // Extract transition probabilities as scalars for fast access
  double a11 = transition_matrix(0, 0);
  double a12 = transition_matrix(0, 1);
  double a21 = transition_matrix(1, 0);
  double a22 = transition_matrix(1, 1);

  // ---------------------------------------------------------------
  // Precompute previous and next valid observation indices
  // ---------------------------------------------------------------
  IntegerVector prev_valid(n, NA_INTEGER);
  IntegerVector next_valid(n, NA_INTEGER);

  int last_valid = -1;
  for (int i = 0; i < n; ++i) {
    if (!NumericVector::is_na(count_data[i])) {
      last_valid = i;
    }
    prev_valid[i] = last_valid;
  }

  last_valid = -1;
  for (int i = n - 1; i >= 0; --i) {
    if (!NumericVector::is_na(count_data[i])) {
      last_valid = i;
    }
    next_valid[i] = last_valid;
  }

  // Find first and last non‑missing positions
  int first_obs = -1, last_obs = -1;
  for (int i = 0; i < n; ++i) {
    if (!NumericVector::is_na(count_data[i])) {
      first_obs = i;
      break;
    }
  }
  for (int i = n - 1; i >= 0; --i) {
    if (!NumericVector::is_na(count_data[i])) {
      last_obs = i;
      break;
    }
  }
  if (first_obs == -1) {
    stop("'count_data' contains no non‑missing values.");
  }

  // ---------------------------------------------------------------
  // Allocate matrices for alpha, beta, and scaling factors
  // ---------------------------------------------------------------
  NumericMatrix alpha(num_states, n);
  NumericMatrix beta(num_states, n);
  NumericVector scalar(n, NA_REAL);

  // ---------------------------------------------------------------
  // Forward procedure
  // ---------------------------------------------------------------
  // Initialization at first valid observation
  int obs = count_data[first_obs];
  double e1 = emissions_matrix(0, obs);
  double e2 = emissions_matrix(1, obs);

  double a1 = init_probs[0] * e1;
  double a2 = init_probs[1] * e2;
  double sum_a = a1 + a2;

  alpha(0, first_obs) = a1 / sum_a;
  alpha(1, first_obs) = a2 / sum_a;
  scalar[first_obs] = sum_a;

  // Recursion for remaining positions
  for (int pos = first_obs + 1; pos < n; ++pos) {
    if (NumericVector::is_na(count_data[pos])) {
      continue;
    }
    int prev_pos = prev_valid[pos - 1];
    if (prev_pos == -1) continue;  // should not happen

    double prev_a1 = alpha(0, prev_pos);
    double prev_a2 = alpha(1, prev_pos);

    obs = count_data[pos];
    e1 = emissions_matrix(0, obs);
    e2 = emissions_matrix(1, obs);

    // Unrolled computation for 2 states
    double new_a1 = e1 * (prev_a1 * a11 + prev_a2 * a21);
    double new_a2 = e2 * (prev_a1 * a12 + prev_a2 * a22);

    sum_a = new_a1 + new_a2;
    alpha(0, pos) = new_a1 / sum_a;
    alpha(1, pos) = new_a2 / sum_a;
    scalar[pos] = sum_a;
  }

  // ---------------------------------------------------------------
  // Backward procedure
  // ---------------------------------------------------------------
  // Initialization at last valid observation
  beta(0, last_obs) = 1.0;
  beta(1, last_obs) = 1.0;

  // Backward recursion
  for (int pos = last_obs - 1; pos >= 0; --pos) {
    if (NumericVector::is_na(count_data[pos])) {
      continue;
    }
    int next_pos = next_valid[pos + 1];
    if (next_pos == -1) continue;

    double s_next = scalar[next_pos];
    obs = count_data[next_pos];
    double e1_next = emissions_matrix(0, obs);
    double e2_next = emissions_matrix(1, obs);

    double beta1_next = beta(0, next_pos);
    double beta2_next = beta(1, next_pos);

    // Unrolled formula for 2 states:
    // beta_i = sum_j (a_ij * e_j(next) * beta_j(next)) / scalar[next]
    double new_beta1 = (a11 * e1_next * beta1_next + a12 * e2_next * beta2_next) / s_next;
    double new_beta2 = (a21 * e1_next * beta1_next + a22 * e2_next * beta2_next) / s_next;

    beta(0, pos) = new_beta1;
    beta(1, pos) = new_beta2;
  }

  // ---------------------------------------------------------------
  // Compute log‑likelihood
  // ---------------------------------------------------------------
  double loglik = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!NumericVector::is_na(scalar[i])) {
      loglik += log(scalar[i]);
    }
  }

  // ---------------------------------------------------------------
  // Return results
  // ---------------------------------------------------------------
  return List::create(
    Named("log_probability_seq") = loglik,
    Named("alpha") = alpha,
    Named("beta") = beta,
    Named("scalar") = scalar
  );
}

// [[Rcpp::export]]
List compute_baum_welch_expectations(NumericVector count_data,
                                     NumericMatrix alpha,
                                     NumericMatrix beta,
                                     NumericVector scalar,
                                     NumericMatrix transition_matrix,
                                     NumericMatrix emissions_matrix) {
  int n = count_data.size();
  int num_states = 2;
  int num_emissions = emissions_matrix.ncol();

  // Precompute valid positions and mapping to next/prev valid
  IntegerVector prev_valid(n, NA_INTEGER);
  IntegerVector next_valid(n, NA_INTEGER);
  int last_valid = -1;
  for (int i = 0; i < n; ++i) {
    if (!NumericVector::is_na(count_data[i])) last_valid = i;
    prev_valid[i] = last_valid;
  }
  last_valid = -1;
  for (int i = n - 1; i >= 0; --i) {
    if (!NumericVector::is_na(count_data[i])) last_valid = i;
    next_valid[i] = last_valid;
  }

  // Extract transition probabilities (log‑space not needed here)
  double a11 = transition_matrix(0, 0);
  double a12 = transition_matrix(0, 1);
  double a21 = transition_matrix(1, 0);
  double a22 = transition_matrix(1, 1);

  // Transition expectations (2x2 matrix)
  NumericMatrix trans_exp(num_states, num_states);

  // Emission expectations (2 x num_emissions matrix)
  NumericMatrix emiss_exp(num_states, num_emissions);

  // Accumulate transition expectations
  for (int pos = 0; pos < n - 1; ++pos) {
    if (NumericVector::is_na(count_data[pos])) continue;
    int next_pos = next_valid[pos + 1];
    if (next_pos == -1 || next_pos >= n) continue;

    int obs_next = count_data[next_pos];
    double e1_next = emissions_matrix(0, obs_next);
    double e2_next = emissions_matrix(1, obs_next);

    double alpha1 = alpha(0, pos);
    double alpha2 = alpha(1, pos);
    double beta1_next = beta(0, next_pos);
    double beta2_next = beta(1, next_pos);
    double s_next = scalar[next_pos];

    // xi(i,j) = alpha_i * a_ij * e_j(next) * beta_j_next / s_next
    double xi11 = alpha1 * a11 * e1_next * beta1_next / s_next;
    double xi12 = alpha1 * a12 * e2_next * beta2_next / s_next;
    double xi21 = alpha2 * a21 * e1_next * beta1_next / s_next;
    double xi22 = alpha2 * a22 * e2_next * beta2_next / s_next;

    trans_exp(0, 0) += xi11;
    trans_exp(0, 1) += xi12;
    trans_exp(1, 0) += xi21;
    trans_exp(1, 1) += xi22;
  }

  // Accumulate emission expectations
  for (int pos = 0; pos < n; ++pos) {
    if (NumericVector::is_na(count_data[pos])) continue;
    int obs = count_data[pos];
    double gamma1 = alpha(0, pos) * beta(0, pos);
    double gamma2 = alpha(1, pos) * beta(1, pos);

    emiss_exp(0, obs) += gamma1;
    emiss_exp(1, obs) += gamma2;
  }

  return List::create(
    Named("trans_exp") = trans_exp,
    Named("emiss_exp") = emiss_exp
  );
}