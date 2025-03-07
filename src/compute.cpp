#include <Rcpp.h>
using namespace Rcpp;

//' Compute the generalized binomial test statistic
//' 
//' @noRd
// [[Rcpp::export]]
double compute_GBT(const std::vector<int> &s, const std::vector<double> &p) {
  const int n = s.size();
  std::vector<double> f(n + 1);
  std::vector<double> g(n + 1);
  f[0] = 1 - p[0];
  f[1] = p[0];
  for (int i = 1; i < n; ++i) {
    g[0] = (1 - p[i]) * f[0];
    for (int j = 1; j < (i + 1); ++j) {
      g[j] = (1 - p[i]) * f[j] + p[i] * f[j - 1];
    }
    g[i + 1] = p[i] * f[i];
    for (int j = 0; j <= (i + 1); ++j) {
      f[j] = g[j];
    }
  }
  int sum_s = 0;
  for (int i = 0; i < n; ++i) {
    sum_s += s[i];
  }
  double sum_f = 0;
  for (int i = sum_s; i < (n + 1); ++i) {
    sum_f += f[i];
  }
  return sum_f;
}

//' Compute the M4 statistic
//' 
//' @noRd
// [[Rcpp::export]]
double compute_M4(const std::vector<int> &s_1,
                  const std::vector<int> &s_0,
                  const std::vector<double> &p_1,
                  const std::vector<double> &p_0,
                  const std::vector<double> &p_n) {
  const int n = s_1.size();
  std::vector<std::vector<double>> f(n + 1, std::vector<double>(n + 1));
  std::vector<std::vector<double>> F(n + 1, std::vector<double>(n + 1));
  std::vector<std::vector<double>> g(n + 1, std::vector<double>(n + 1));
  f[0][0] = p_n[0];
  f[0][1] = p_0[0];
  f[1][0] = p_1[0];
  for (int i = 1; i < n; ++i) {
    g[0][0] = p_n[i] * f[0][0];
    for (int j = 1; j < (i + 1); ++j) {
      g[0][j] = p_0[i] * f[0][j - 1] + p_n[i] * f[0][j];
      g[j][0] = p_1[i] * f[j - 1][0] + p_n[i] * f[j][0];
      for (int k = 1; k <= (i - j + 1); ++k) {
        g[j][k] = p_1[i] * f[j - 1][k] + p_0[i] * f[j][k - 1] +
          p_n[i] * f[j][k];
      }
    }
    g[0][i + 1] = p_0[i] * f[0][i];
    g[i + 1][0] = p_1[i] * f[i][0];
    for (int j = 0; j <= (i + 1); ++j) {
      for (int k = 0; k <= (i - j + 1); ++k) {
        f[j][k] = g[j][k];
      }
    }
  }
  for (int i = 0; i < (n + 1); ++i) {
    for (int j = 0; j < (n - i + 1); ++j) {
      for (int k = i; k < (n - j + 1); ++k) {
        for (int l = j; l < (n - k + 1); ++l) {
          F[i][j] += f[k][l];
        }
      }
    }
  }
  int sum_s_1 = 0;
  int sum_s_0 = 0;
  double sum_f = 0;
  for (int i = 0; i < n; ++i) {
    sum_s_1 += s_1[i];
    sum_s_0 += s_0[i];
  }
  for (int i = 0; i < (n + 1); ++i) {
    for (int j = 0; j < (n - i + 1); ++j) {
      if (F[i][j] <= F[sum_s_1][sum_s_0]) {
        sum_f += f[i][j];
      }
    }
  }
  return sum_f;
}
