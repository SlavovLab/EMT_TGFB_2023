// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <random>
#include <vector>
#include <cmath>

using namespace Rcpp;

// ------------------------------------------------------------------
// Helpers
// ------------------------------------------------------------------
void build_complete_pair(const NumericVector& x,
                         const NumericVector& y,
                         std::vector<double>& x_out,
                         std::vector<double>& y_out) {
  int n = x.size();
  x_out.clear();
  y_out.clear();
  x_out.reserve(n);
  y_out.reserve(n);

  for (int i = 0; i < n; ++i) {
    if (!NumericVector::is_na(x[i]) && !NumericVector::is_na(y[i])) {
      x_out.push_back(x[i]);
      y_out.push_back(y[i]);
    }
  }
}

double cor_pearson_std(const std::vector<double>& x,
                       const std::vector<double>& y) {
  int n = x.size();
  if (y.size() != n || n < 2) return NA_REAL;

  double sumx = 0.0, sumy = 0.0;
  for (int i = 0; i < n; ++i) {
    sumx += x[i];
    sumy += y[i];
  }
  double meanx = sumx / n;
  double meany = sumy / n;

  double num = 0.0, sxx = 0.0, syy = 0.0;
  for (int i = 0; i < n; ++i) {
    double dx = x[i] - meanx;
    double dy = y[i] - meany;
    num += dx * dy;
    sxx += dx * dx;
    syy += dy * dy;
  }
  if (sxx <= 0.0 || syy <= 0.0) return NA_REAL;
  return num / std::sqrt(sxx * syy);
}

double cor_perm_x_std(const std::vector<double>& x,
                      const std::vector<double>& y,
                      const std::vector<int>& perm) {
  int n = x.size();
  if (y.size() != n || (int)perm.size() != n) return NA_REAL;

  double sumx = 0.0, sumy = 0.0;
  for (int i = 0; i < n; ++i) {
    double xi = x[perm[i]];
    double yi = y[i];
    sumx += xi;
    sumy += yi;
  }
  double meanx = sumx / n;
  double meany = sumy / n;

  double num = 0.0, sxx = 0.0, syy = 0.0;
  for (int i = 0; i < n; ++i) {
    double xi = x[perm[i]];
    double yi = y[i];
    double dx = xi - meanx;
    double dy = yi - meany;
    num += dx * dy;
    sxx += dx * dx;
    syy += dy * dy;
  }
  if (sxx <= 0.0 || syy <= 0.0) return NA_REAL;
  return num / std::sqrt(sxx * syy);
}

double cor_perm_y_std(const std::vector<double>& x,
                      const std::vector<double>& y,
                      const std::vector<int>& perm) {
  int n = x.size();
  if (y.size() != n || (int)perm.size() != n) return NA_REAL;

  double sumx = 0.0, sumy = 0.0;
  for (int i = 0; i < n; ++i) {
    double xi = x[i];
    double yi = y[perm[i]];
    sumx += xi;
    sumy += yi;
  }
  double meanx = sumx / n;
  double meany = sumy / n;

  double num = 0.0, sxx = 0.0, syy = 0.0;
  for (int i = 0; i < n; ++i) {
    double xi = x[i];
    double yi = y[perm[i]];
    double dx = xi - meanx;
    double dy = yi - meany;
    num += dx * dy;
    sxx += dx * dx;
    syy += dy * dy;
  }
  if (sxx <= 0.0 || syy <= 0.0) return NA_REAL;
  return num / std::sqrt(sxx * syy);
}

void shuffle_perm(std::vector<int>& perm, std::mt19937& rng) {
  int n = perm.size();
  for (int i = n - 1; i > 0; --i) {
    std::uniform_int_distribution<int> dist(0, i);
    int j = dist(rng);
    std::swap(perm[i], perm[j]);
  }
}

// ------------------------------------------------------------------
// Main exported function
// ------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List perm_diff_cor_both_cpp(Rcpp::NumericVector x1,
                                  Rcpp::NumericVector y1,
                                  Rcpp::NumericVector x2,
                                  Rcpp::NumericVector y2,
                                  int N              = 2000,
                                  bool do_trick      = true,
                                  bool return_nulls  = true,
                                  int seed           = 12345) {
  if (N <= 0) {
    N = 1;
  }

  // If lengths somehow differ, just return NA rather than erroring.
  if (x1.size() != y1.size() || x2.size() != y2.size()) {
    int n1 = x1.size();
    int n2 = x2.size();
    return Rcpp::List::create(
      Rcpp::_["r1_obs"]   = NA_REAL,
      Rcpp::_["r2_obs"]   = NA_REAL,
      Rcpp::_["D_obs"]    = NA_REAL,
      Rcpp::_["n1_eff"]   = n1,
      Rcpp::_["n2_eff"]   = n2,
      Rcpp::_["R1"]       = NumericVector(0),
      Rcpp::_["R2"]       = NumericVector(0),
      Rcpp::_["D_null"]   = NumericVector(0),
      Rcpp::_["p_direct"] = NA_REAL,
      Rcpp::_["p_trick"]  = NA_REAL
    );
  }

  std::vector<double> x1c, y1c, x2c, y2c;
  build_complete_pair(x1, y1, x1c, y1c);
  build_complete_pair(x2, y2, x2c, y2c);

  int n1_eff = x1c.size();
  int n2_eff = x2c.size();

  if (n1_eff < 4 || n2_eff < 4) {
    return Rcpp::List::create(
      Rcpp::_["r1_obs"]   = NA_REAL,
      Rcpp::_["r2_obs"]   = NA_REAL,
      Rcpp::_["D_obs"]    = NA_REAL,
      Rcpp::_["n1_eff"]   = n1_eff,
      Rcpp::_["n2_eff"]   = n2_eff,
      Rcpp::_["R1"]       = NumericVector(0),
      Rcpp::_["R2"]       = NumericVector(0),
      Rcpp::_["D_null"]   = NumericVector(0),
      Rcpp::_["p_direct"] = NA_REAL,
      Rcpp::_["p_trick"]  = NA_REAL
    );
  }

  double r1_obs = cor_pearson_std(x1c, y1c);
  double r2_obs = cor_pearson_std(x2c, y2c);

  if (NumericVector::is_na(r1_obs) || NumericVector::is_na(r2_obs)) {
    return Rcpp::List::create(
      Rcpp::_["r1_obs"]   = r1_obs,
      Rcpp::_["r2_obs"]   = r2_obs,
      Rcpp::_["D_obs"]    = NA_REAL,
      Rcpp::_["n1_eff"]   = n1_eff,
      Rcpp::_["n2_eff"]   = n2_eff,
      Rcpp::_["R1"]       = NumericVector(0),
      Rcpp::_["R2"]       = NumericVector(0),
      Rcpp::_["D_null"]   = NumericVector(0),
      Rcpp::_["p_direct"] = NA_REAL,
      Rcpp::_["p_trick"]  = NA_REAL
    );
  }

  double D_obs = std::fabs(std::fabs(r1_obs) - std::fabs(r2_obs));

  // allocate only if needed
  NumericVector R1_vec, R2_vec, D_null_vec;
  if (return_nulls || do_trick) {
    R1_vec = NumericVector(N);
    R2_vec = NumericVector(N);
    if (return_nulls) {
      D_null_vec = NumericVector(N);
    }
  }

  std::vector<int> perm1(n1_eff), perm2(n2_eff);
  for (int i = 0; i < n1_eff; ++i) perm1[i] = i;
  for (int i = 0; i < n2_eff; ++i) perm2[i] = i;

  std::mt19937 rng(seed);
  int count_direct = 0;

  // generate null draws
  for (int b = 0; b < N; ++b) {
    shuffle_perm(perm1, rng);
    shuffle_perm(perm2, rng);

    double r1_b = cor_perm_x_std(x1c, y1c, perm1);
    double r2_b = cor_perm_y_std(x2c, y2c, perm2);
    double D_b  = std::fabs(std::fabs(r1_b) - std::fabs(r2_b));

    if (return_nulls || do_trick) {
      R1_vec[b] = r1_b;
      R2_vec[b] = r2_b;
      if (return_nulls) {
        D_null_vec[b] = D_b;
      }
    }

    if (D_b >= D_obs) {
      count_direct++;
    }
  }

  double p_direct = static_cast<double>(count_direct) /
                    static_cast<double>(N);

  double p_trick = NA_REAL;
  if (do_trick) {
    long long count_trick = 0;
    for (int i = 0; i < N; ++i) {
      double a = std::fabs(R1_vec[i]);
      for (int j = 0; j < N; ++j) {
        double b = std::fabs(R2_vec[j]);
        double D = std::fabs(a - b);
        if (D >= D_obs) {
          count_trick++;
        }
      }
    }
    double denom = static_cast<double>(N) * static_cast<double>(N);
    p_trick = static_cast<double>(count_trick) / denom;
  }

  // If we don't want nulls, return empty vectors
  if (!return_nulls) {
    R1_vec     = NumericVector(0);
    R2_vec     = NumericVector(0);
    D_null_vec = NumericVector(0);
  }

  return Rcpp::List::create(
    Rcpp::_["r1_obs"]   = r1_obs,
    Rcpp::_["r2_obs"]   = r2_obs,
    Rcpp::_["D_obs"]    = D_obs,
    Rcpp::_["n1_eff"]   = n1_eff,
    Rcpp::_["n2_eff"]   = n2_eff,
    Rcpp::_["R1"]       = R1_vec,
    Rcpp::_["R2"]       = R2_vec,
    Rcpp::_["D_null"]   = D_null_vec,
    Rcpp::_["p_direct"] = p_direct,
    Rcpp::_["p_trick"]  = p_trick
  );
}
