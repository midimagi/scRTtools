#include <Rcpp.h>
#include <vector>
#include <limits>
#include <cmath>

using namespace Rcpp;

// penalty_code:
// 0 = none
// 1 = AIC
// 2 = BIC

inline size_t mat_idx(const int i, const int j, const int n) {
  return static_cast<size_t>(i) + static_cast<size_t>(n) * static_cast<size_t>(j);
}

// [[Rcpp::export]]
IntegerVector findchangepts_cpp_max2(
    NumericVector x,
    int MinSegmentLength,
    int penalty_code
) {
  const int n = x.size();

  if (n < 2 * MinSegmentLength) {
    return IntegerVector(0);
  }

  // Prefix sums
  NumericVector ct(n + 1), cy(n + 1), ct2(n + 1), cy2(n + 1), cty(n + 1);

  for (int i = 0; i < n; ++i) {
    const double t = static_cast<double>(i + 1);
    const double xi = x[i];

    ct[i + 1]  = ct[i]  + t;
    cy[i + 1]  = cy[i]  + xi;
    ct2[i + 1] = ct2[i] + t * t;
    cy2[i + 1] = cy2[i] + xi * xi;
    cty[i + 1] = cty[i] + t * xi;
  }

  const double inf = std::numeric_limits<double>::infinity();
  const double eps = std::numeric_limits<double>::epsilon();

  // RSS matrix in column-major layout
  std::vector<double> rss_mat(static_cast<size_t>(n) * static_cast<size_t>(n), inf);

  // Compute RSS(i, j) for all valid segments x[(i+1):(j+1)] in 0-based indexing
  for (int j = MinSegmentLength - 1; j < n; ++j) {
    if ((j & 127) == 0) {
      Rcpp::checkUserInterrupt();
    }

    const int i_max = j - MinSegmentLength + 1;

    const double ct_j1  = ct[j + 1];
    const double cy_j1  = cy[j + 1];
    const double ct2_j1 = ct2[j + 1];
    const double cy2_j1 = cy2[j + 1];
    const double cty_j1 = cty[j + 1];

    for (int i = 0; i <= i_max; ++i) {
      const double m = static_cast<double>(j - i + 1);

      const double sx  = ct_j1  - ct[i];
      const double sy  = cy_j1  - cy[i];
      const double sxx = ct2_j1 - ct2[i];
      const double syy = cy2_j1 - cy2[i];
      const double sxy = cty_j1 - cty[i];

      const double Sxx = sxx - sx * sx / m;
      const double Syy = syy - sy * sy / m;
      const double Sxy = sxy - sx * sy / m;

      double rss;
      if (!std::isfinite(Sxx) || Sxx <= eps) {
        rss = Syy;
      } else {
        rss = Syy - (Sxy * Sxy) / Sxx;
      }

      if (!std::isfinite(rss)) {
        rss = inf;
      } else if (rss < 0.0) {
        rss = 0.0;
      }

      rss_mat[mat_idx(i, j, n)] = rss;
    }
  }

  // ------------------------------------------------------------
  // Specialized DP for MaxNumChanges = 2
  //
  // k = 0  => 1 segment
  // k = 1  => 2 segments
  // k = 2  => 3 segments
  // ------------------------------------------------------------

  // One segment
  const double rss1 = rss_mat[mat_idx(0, n - 1, n)];

  // Two segments
  std::vector<double> dp2(n, inf);
  std::vector<int> prev2(n, NA_INTEGER);

  for (int j = 2 * MinSegmentLength - 1; j < n; ++j) {
    double best_val = inf;
    int best_i = NA_INTEGER;

    const int i_min = MinSegmentLength - 1;
    const int i_max = j - MinSegmentLength;

    for (int i = i_min; i <= i_max; ++i) {
      const double left  = rss_mat[mat_idx(0, i, n)];
      const double right = rss_mat[mat_idx(i + 1, j, n)];
      const double cand  = left + right;

      if (cand < best_val) {
        best_val = cand;
        best_i = i;
      }
    }

    if (std::isfinite(best_val)) {
      dp2[j] = best_val;
      prev2[j] = best_i;
    }
  }

  // Three segments
  std::vector<double> dp3(n, inf);
  std::vector<int> prev3(n, NA_INTEGER);

  if (n >= 3 * MinSegmentLength) {
    for (int j = 3 * MinSegmentLength - 1; j < n; ++j) {
      if ((j & 127) == 0) {
        Rcpp::checkUserInterrupt();
      }

      double best_val = inf;
      int best_i = NA_INTEGER;

      const int i_min = 2 * MinSegmentLength - 1;
      const int i_max = j - MinSegmentLength;

      for (int i = i_min; i <= i_max; ++i) {
        const double left = dp2[i];
        if (!std::isfinite(left)) {
          continue;
        }

        const double right = rss_mat[mat_idx(i + 1, j, n)];
        const double cand  = left + right;

        if (cand < best_val) {
          best_val = cand;
          best_i = i;
        }
      }

      if (std::isfinite(best_val)) {
        dp3[j] = best_val;
        prev3[j] = best_i;
      }
    }
  }

  // ------------------------------------------------------------
  // Model selection among k = 0, 1, 2
  // ------------------------------------------------------------
  double crit0 = inf, crit1 = inf, crit2 = inf;

  if (std::isfinite(rss1)) {
    if (penalty_code == 0) {
      crit0 = rss1;
    } else if (penalty_code == 1) {
      crit0 = n * std::log((rss1 + 1e-12) / n) + 2.0 * 2.0;
    } else {
      crit0 = n * std::log((rss1 + 1e-12) / n) + std::log((double)n) * 2.0;
    }
  }

  if (std::isfinite(dp2[n - 1])) {
    const double rss2 = dp2[n - 1];
    const double p = 5.0;  // 2*2 + 1
    if (penalty_code == 0) {
      crit1 = rss2;
    } else if (penalty_code == 1) {
      crit1 = n * std::log((rss2 + 1e-12) / n) + 2.0 * p;
    } else {
      crit1 = n * std::log((rss2 + 1e-12) / n) + std::log((double)n) * p;
    }
  }

  if (std::isfinite(dp3[n - 1])) {
    const double rss3 = dp3[n - 1];
    const double p = 8.0;  // 2*3 + 2
    if (penalty_code == 0) {
      crit2 = rss3;
    } else if (penalty_code == 1) {
      crit2 = n * std::log((rss3 + 1e-12) / n) + 2.0 * p;
    } else {
      crit2 = n * std::log((rss3 + 1e-12) / n) + std::log((double)n) * p;
    }
  }

  int best_k = 0;
  double best_crit = crit0;

  if (crit1 < best_crit) {
    best_crit = crit1;
    best_k = 1;
  }
  if (crit2 < best_crit) {
    best_crit = crit2;
    best_k = 2;
  }

  if (best_k == 0) {
    return IntegerVector(0);
  }

  if (best_k == 1) {
    const int cp1 = prev2[n - 1];
    if (cp1 == NA_INTEGER || cp1 < 0 || cp1 >= (n - 1)) {
      return IntegerVector(0);
    }
    return IntegerVector::create(cp1 + 1);
  }

  // best_k == 2
  const int cp2 = prev3[n - 1];
  if (cp2 == NA_INTEGER || cp2 < 0 || cp2 >= (n - 1)) {
    return IntegerVector(0);
  }

  const int cp1 = prev2[cp2];
  if (cp1 == NA_INTEGER || cp1 < 0 || cp1 >= cp2) {
    return IntegerVector(0);
  }

  return IntegerVector::create(cp1 + 1, cp2 + 1);
}