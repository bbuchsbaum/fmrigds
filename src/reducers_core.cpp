// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using arma::mat;
using arma::vec;
using arma::uword;

inline bool finite_(double x) { return std::isfinite(x); }

// Runtime control of OpenMP threads
// [[Rcpp::export]]
void set_omp_threads(const int n) {
#ifdef _OPENMP
  if (n > 0) {
    omp_set_num_threads(n);
  }
#else
  (void)n; // suppress unused warning when OpenMP is not enabled
#endif
}

// ----- Fixed-Effects meta (beta/var) -----

// tail: 0 = two.sided, 1 = less, 2 = greater
// [[Rcpp::export]]
Rcpp::List meta_fe_cpp(const arma::mat& beta, const arma::mat& var,
                       const int min_subj = 2,
                       const double eps = 1e-12,
                       const int tail = 0) {
  const uword S = beta.n_rows;
  const uword B = beta.n_cols;
  NumericVector mu(B), var_mu(B), se_mu(B), z(B), p(B), Q(B), I2(B);

  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (long b = 0; b < static_cast<long>(B); ++b) {
    double sw = 0.0, wy = 0.0, sum_w2 = 0.0;
    int k = 0;

    for (uword i = 0; i < S; ++i) {
      double y = beta(i, b);
      double v = var(i, b);
      if (finite_(y) && finite_(v) && v > 0.0) {
        double w = 1.0 / std::max(v, eps);
        sw += w; wy += w * y; sum_w2 += w * w; ++k;
      }
    }
    if (k < min_subj || sw <= 0.0) {
      mu[b] = var_mu[b] = se_mu[b] = z[b] = p[b] = Q[b] = I2[b] = NA_REAL;
      continue;
    }

    double mu_b = wy / sw;
    double var_b = 1.0 / sw;
    double se_b  = std::sqrt(var_b);

    double Qb = 0.0;
    for (uword i = 0; i < S; ++i) {
      double y = beta(i, b);
      double v = var(i, b);
      if (finite_(y) && finite_(v) && v > 0.0) {
        double w = 1.0 / std::max(v, eps);
        double r = y - mu_b;
        Qb += w * r * r;
      }
    }
    double I2b = (Qb > eps) ? std::max(0.0, (Qb - (k - 1.0)) / Qb) : 0.0;
    double zb = (se_b > 0.0) ? (mu_b / se_b) : NA_REAL;
    double pb = NA_REAL;
    if (finite_(zb)) {
      if (tail == 0)      pb = 2.0 * R::pnorm(-std::fabs(zb), 0.0, 1.0, 1, 0);
      else if (tail == 1) pb = R::pnorm(zb, 0.0, 1.0, 1, 0);
      else                pb = R::pnorm(zb, 0.0, 1.0, 0, 0);
    }

    mu[b]     = mu_b;
    var_mu[b] = var_b;
    se_mu[b]  = se_b;
    z[b]      = zb;
    p[b]      = pb;
    Q[b]      = Qb;
    I2[b]     = I2b;
  }

  return Rcpp::List::create(
    _["beta_g"] = mu,
    _["var_g"]  = var_mu,
    _["se_g"]   = se_mu,
    _["z_g"]    = z,
    _["p_g"]    = p,
    _["Q"]      = Q,
    _["I2"]     = I2
  );
}

// ----- Random-Effects meta (DL tau2) -----

// [[Rcpp::export]]
Rcpp::List meta_re_dl_cpp(const arma::mat& beta, const arma::mat& var,
                          const int min_subj = 2,
                          const double eps = 1e-12,
                          const int tail = 0) {
  const uword S = beta.n_rows;
  const uword B = beta.n_cols;
  NumericVector mu(B), var_mu(B), se_mu(B), z(B), p(B), Q(B), I2(B), tau2(B);

  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (long b = 0; b < static_cast<long>(B); ++b) {
    double sw = 0.0, wy = 0.0, sum_w2 = 0.0;
    int k = 0;

    for (uword i = 0; i < S; ++i) {
      double y = beta(i, b);
      double v = var(i, b);
      if (finite_(y) && finite_(v) && v > 0.0) {
        double w = 1.0 / std::max(v, eps);
        sw += w; wy += w * y; sum_w2 += w * w; ++k;
      }
    }
    if (k < min_subj || sw <= 0.0) {
      mu[b] = var_mu[b] = se_mu[b] = z[b] = p[b] = Q[b] = I2[b] = tau2[b] = NA_REAL;
      continue;
    }
    double mu_fe = wy / sw;

    double Qb = 0.0;
    for (uword i = 0; i < S; ++i) {
      double y = beta(i, b);
      double v = var(i, b);
      if (finite_(y) && finite_(v) && v > 0.0) {
        double w = 1.0 / std::max(v, eps);
        double r = y - mu_fe;
        Qb += w * r * r;
      }
    }
    double C = sw - (sum_w2 / sw);
    double t2 = std::max(0.0, (Qb - (k - 1.0)) / std::max(C, eps));

    double sws = 0.0, wys = 0.0;
    for (uword i = 0; i < S; ++i) {
      double y = beta(i, b);
      double v = var(i, b);
      if (finite_(y) && finite_(v) && v > 0.0) {
        double w = 1.0 / (std::max(v, eps) + t2);
        sws += w; wys += w * y;
      }
    }
    if (sws <= 0.0) {
      mu[b] = var_mu[b] = se_mu[b] = z[b] = p[b] = Q[b] = I2[b] = tau2[b] = NA_REAL;
      continue;
    }
    double mu_re = wys / sws;
    double var_re = 1.0 / sws;
    double se_re  = std::sqrt(var_re);

    double I2b = (Qb > eps) ? std::max(0.0, (Qb - (k - 1.0)) / Qb) : 0.0;

    double zb = (se_re > 0.0) ? (mu_re / se_re) : NA_REAL;
    double pb = NA_REAL;
    if (finite_(zb)) {
      if (tail == 0)      pb = 2.0 * R::pnorm(-std::fabs(zb), 0.0, 1.0, 1, 0);
      else if (tail == 1) pb = R::pnorm(zb, 0.0, 1.0, 1, 0);
      else                pb = R::pnorm(zb, 0.0, 1.0, 0, 0);
    }

    mu[b]     = mu_re;
    var_mu[b] = var_re;
    se_mu[b]  = se_re;
    z[b]      = zb;
    p[b]      = pb;
    Q[b]      = Qb;
    I2[b]     = I2b;
    tau2[b]   = t2;
  }

  return Rcpp::List::create(
    _["beta_g"] = mu,
    _["var_g"]  = var_mu,
    _["se_g"]   = se_mu,
    _["z_g"]    = z,
    _["p_g"]    = p,
    _["tau2"]   = tau2,
    _["Q"]      = Q,
    _["I2"]     = I2
  );
}

// ----- Fixed-Effects meta-regression (WLS) -----

// [[Rcpp::export]]
Rcpp::List meta_fe_reg_cpp(const arma::mat& beta, const arma::mat& var,
                           const arma::mat& X,
                           const int min_subj = 2,
                           const double eps = 1e-12) {
  const uword S = beta.n_rows;
  const uword B = beta.n_cols;
  const uword p = X.n_cols;

  arma::mat coef(p, B, arma::fill::value(NA_REAL));
  arma::mat se  (p, B, arma::fill::value(NA_REAL));
  NumericVector Q(B), df_res(B);

  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (long b = 0; b < static_cast<long>(B); ++b) {
    arma::mat G(p, p, arma::fill::zeros);
    arma::vec g(p, arma::fill::zeros);
    int k = 0;

    for (uword i = 0; i < S; ++i) {
      double y = beta(i, b);
      double v = var(i, b);
      if (finite_(y) && finite_(v) && v > 0.0) {
        double w = 1.0 / std::max(v, eps);
        arma::rowvec xi = X.row(i);
        if (!xi.is_finite()) continue;
        G  += w * (xi.t() * xi);
        g  += w * (xi.t() * y);
        ++k;
      }
    }
    if (k < std::max<int>(min_subj, p + 1)) { Q[b] = df_res[b] = NA_REAL; continue; }
    arma::mat A;
    bool ok = arma::inv_sympd(A, G);
    if (!ok) { Q[b] = df_res[b] = NA_REAL; continue; }
    arma::vec bh = A * g;
    double Qb = 0.0;
    for (uword i = 0; i < S; ++i) {
      double y = beta(i, b);
      double v = var(i, b);
      if (finite_(y) && finite_(v) && v > 0.0) {
        double w = 1.0 / std::max(v, eps);
        double r = y - arma::dot(X.row(i), bh);
        Qb += w * r * r;
      }
    }
    coef.col(b) = bh;
    for (uword j = 0; j < p; ++j) se(j, b) = std::sqrt(std::max(A(j, j), 0.0));
    Q[b] = Qb; df_res[b] = static_cast<double>(k - static_cast<int>(p));
  }

  return Rcpp::List::create(
    _["coef"]   = coef,
    _["se_coef"] = se,
    _["Q"]      = Q,
    _["df_res"]  = df_res
  );
}

// ----- Random-Effects meta-regression (DL tau2) -----

// [[Rcpp::export]]
Rcpp::List meta_re_reg_dl_cpp(const arma::mat& beta, const arma::mat& var,
                              const arma::mat& X,
                              const int min_subj = 2,
                              const double eps = 1e-12) {
  const uword S = beta.n_rows;
  const uword B = beta.n_cols;
  const uword p = X.n_cols;

  arma::mat coef(p, B, arma::fill::value(NA_REAL));
  arma::mat se  (p, B, arma::fill::value(NA_REAL));
  NumericVector tau2(B, NA_REAL), Q(B, NA_REAL), df_res(B, NA_REAL);

  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (long b = 0; b < static_cast<long>(B); ++b) {
    arma::mat G(p, p, arma::fill::zeros);
    arma::vec g(p, arma::fill::zeros);
    std::vector<double> w_save; w_save.resize(S, 0.0);
    int k = 0;

    for (uword i = 0; i < S; ++i) {
      double y = beta(i, b);
      double v = var(i, b);
      if (finite_(y) && finite_(v) && v > 0.0) {
        double w = 1.0 / std::max(v, eps);
        w_save[i] = w;
        arma::rowvec xi = X.row(i);
        if (!xi.is_finite()) { w_save[i] = 0.0; continue; }
        G += w * (xi.t() * xi);
        g += w * (xi.t() * y);
        ++k;
      } else {
        w_save[i] = 0.0;
      }
    }
    if (k < std::max<int>(min_subj, p + 1)) { tau2[b] = Q[b] = df_res[b] = NA_REAL; continue; }
    if (arma::rank(G) < p) { tau2[b] = Q[b] = df_res[b] = NA_REAL; continue; }
    arma::mat A;
    bool ok = arma::inv_sympd(A, G);
    if (!ok) { tau2[b] = Q[b] = df_res[b] = NA_REAL; continue; }
    arma::vec bh_fe = A * g;

    double Qb = 0.0, sw = 0.0, trace_correction = 0.0;
    for (uword i = 0; i < S; ++i) {
      double w = w_save[i];
      if (w <= 0.0) continue;
      sw += w;
      double r = beta(i, b) - arma::dot(X.row(i), bh_fe);
      Qb += w * r * r;
      double xAx = arma::as_scalar(X.row(i) * A * X.row(i).t());
      trace_correction += w * w * xAx;
    }
    // P = W - W X (X' W X)^-1 X' W, hence
    // tr(P) = sum(w_i) - sum(w_i^2 x_i' A x_i).
    double C = sw - trace_correction;
    double dfb = static_cast<double>(k - static_cast<int>(p));
    if (!finite_(C) || C <= 0.0) {
      tau2[b] = Q[b] = df_res[b] = NA_REAL;
      continue;
    }
    double t2 = std::max(0.0, (Qb - dfb) / C);

    arma::mat Gs(p, p, arma::fill::zeros);
    arma::vec gs(p, arma::fill::zeros);
    for (uword i = 0; i < S; ++i) {
      double y = beta(i, b);
      double v = var(i, b);
      if (finite_(y) && finite_(v) && v > 0.0) {
        double w = 1.0 / (std::max(v, eps) + t2);
        arma::rowvec xi = X.row(i);
        if (!xi.is_finite()) continue;
        Gs += w * (xi.t() * xi);
        gs += w * (xi.t() * y);
      }
    }
    arma::mat As;
    ok = arma::inv_sympd(As, Gs);
    if (!ok) { tau2[b] = Q[b] = df_res[b] = NA_REAL; continue; }
    arma::vec bh = As * gs;
    coef.col(b) = bh;
    for (uword j = 0; j < p; ++j) se(j, b) = std::sqrt(std::max(As(j, j), 0.0));
    tau2[b] = t2; Q[b] = Qb; df_res[b] = dfb;
  }

  return Rcpp::List::create(
    _["coef"]    = coef,
    _["se_coef"] = se,
    _["tau2"]    = tau2,
    _["Q"]       = Q,
    _["df_res"]  = df_res
  );
}

// ----- Stouffer (Z-combiner) -----

// [[Rcpp::export]]
Rcpp::List stouffer_combine_cpp(const arma::mat& z,
                                Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue,
                                const int min_subj = 1) {
  const uword S = z.n_rows, B = z.n_cols;
  NumericVector Z(B), p(B);

  std::vector<double> w(S, 1.0);
  if (weights.isNotNull()) {
    NumericVector ww(weights);
    if (ww.size() == S) {
      for (uword i = 0; i < S; ++i) w[i] = ww[i];
    }
  }

  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (long b = 0; b < static_cast<long>(B); ++b) {
    double num = 0.0, den2 = 0.0; int k = 0;
    for (uword i = 0; i < S; ++i) {
      double zi = z(i, b);
      if (finite_(zi)) { num += w[i] * zi; den2 += w[i] * w[i]; ++k; }
    }
    if (k < min_subj || den2 <= 0.0) { Z[b] = p[b] = NA_REAL; continue; }
    double Zb = num / std::sqrt(den2);
    Z[b] = Zb;
    p[b] = 2.0 * R::pnorm(-std::fabs(Zb), 0.0, 1.0, 1, 0);
  }
  return Rcpp::List::create(_["z_g"] = Z, _["p_g"] = p);
}

// ----- Fisher (p-combiner) -----

inline double clamp_p(double x) {
  if (!finite_(x)) return NA_REAL;
  if (x <= 0.0)    return 1e-300;
  if (x >= 1.0)    return 1.0 - 1e-16;
  return x;
}

// [[Rcpp::export]]
Rcpp::List fisher_combine_cpp(const arma::mat& pmat, const int min_subj = 1) {
  const uword S = pmat.n_rows, B = pmat.n_cols;
  NumericVector stat(B), df(B), p(B);

  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (long b = 0; b < static_cast<long>(B); ++b) {
    double X2 = 0.0; int k = 0;
    for (uword i = 0; i < S; ++i) {
      double pi = clamp_p(pmat(i, b));
      if (finite_(pi)) { X2 += -2.0 * std::log(pi); ++k; }
    }
    if (k < min_subj) { stat[b]=df[b]=p[b]=NA_REAL; continue; }
    double dfb = 2.0 * k;
    double pb  = R::pchisq(X2, dfb, /*lower_tail*/0, /*log_p*/0);
    stat[b] = X2; df[b] = dfb; p[b] = pb;
  }
  return Rcpp::List::create(_["chi2"] = stat, _["df"] = df, _["p_g"] = p);
}

// ----- Lancaster (weighted Fisher) -----

// [[Rcpp::export]]
Rcpp::List lancaster_combine_cpp(const arma::mat& pmat,
                                 const Rcpp::IntegerVector dfw,
                                 const int min_subj = 1) {
  const uword S = pmat.n_rows, B = pmat.n_cols;
  NumericVector stat(B), df(B), p(B);

  if (static_cast<uword>(dfw.size()) != S)
    stop("dfw length must equal number of subjects (rows of p).");

  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (long b = 0; b < static_cast<long>(B); ++b) {
    double X2 = 0.0; int df_sum = 0; int k = 0;
    for (uword i = 0; i < S; ++i) {
      double pi = clamp_p(pmat(i, b));
      if (finite_(pi)) {
        int wi = std::max(dfw[i], 1);
        double chi_i = R::qchisq(1.0 - pi, 2.0 * wi, /*lower_tail*/1, /*log_p*/0);
        X2 += chi_i; df_sum += 2 * wi; ++k;
      }
    }
    if (k < min_subj || df_sum <= 0) { stat[b]=df[b]=p[b]=NA_REAL; continue; }
    double pb = R::pchisq(X2, (double)df_sum, /*lower_tail*/0, /*log_p*/0);
    stat[b] = X2; df[b] = (double)df_sum; p[b] = pb;
  }
  return Rcpp::List::create(_["chi2"] = stat, _["df"] = df, _["p_g"] = p);
}

// ----- Permutation tests ----------------------------------------------------

inline double perm_tail_p(const double obs, const double null_stat, const int tail) {
  if (!finite_(obs) || !finite_(null_stat)) return 0.0;
  if (tail == 1) return null_stat <= obs ? 1.0 : 0.0;       // less
  if (tail == 2) return null_stat >= obs ? 1.0 : 0.0;       // greater
  return std::fabs(null_stat) >= std::fabs(obs) ? 1.0 : 0.0; // two.sided
}

inline double mean_from_sums(const double sum, const int n) {
  return n > 0 ? sum / static_cast<double>(n) : NA_REAL;
}

inline double t_onesample_from_sums(const double sum, const double sumsq, const int n) {
  if (n < 2) return NA_REAL;
  const double mean = sum / static_cast<double>(n);
  double var = (sumsq - static_cast<double>(n) * mean * mean) / static_cast<double>(n - 1);
  if (var < 0.0 && var > -1e-10) var = 0.0;
  if (var <= 0.0 || !finite_(var)) return NA_REAL;
  return mean / std::sqrt(var / static_cast<double>(n));
}

inline double t_twosample_from_sums(const double sum0, const double sumsq0, const int n0,
                                    const double sum1, const double sumsq1, const int n1,
                                    const bool welch) {
  if (n0 < 2 || n1 < 2) return NA_REAL;
  const double m0 = sum0 / static_cast<double>(n0);
  const double m1 = sum1 / static_cast<double>(n1);
  double v0 = (sumsq0 - static_cast<double>(n0) * m0 * m0) / static_cast<double>(n0 - 1);
  double v1 = (sumsq1 - static_cast<double>(n1) * m1 * m1) / static_cast<double>(n1 - 1);
  if (v0 < 0.0 && v0 > -1e-10) v0 = 0.0;
  if (v1 < 0.0 && v1 > -1e-10) v1 = 0.0;
  if (!finite_(v0) || !finite_(v1)) return NA_REAL;
  double se2;
  if (welch) {
    se2 = v0 / static_cast<double>(n0) + v1 / static_cast<double>(n1);
  } else {
    const double sp2 = ((n0 - 1.0) * v0 + (n1 - 1.0) * v1) / static_cast<double>(n0 + n1 - 2);
    se2 = sp2 * (1.0 / static_cast<double>(n0) + 1.0 / static_cast<double>(n1));
  }
  if (se2 <= 0.0 || !finite_(se2)) return NA_REAL;
  return (m1 - m0) / std::sqrt(se2);
}

// sign_mat is n_perm x subjects with entries -1/+1.
// [[Rcpp::export]]
Rcpp::List perm_onesample_t_cpp(const arma::mat& beta,
                                const arma::imat& sign_mat,
                                const int tail = 0,
                                const int min_subj = 2) {
  const uword S = beta.n_rows, B = beta.n_cols;
  const uword P = sign_mat.n_rows;
  if (sign_mat.n_cols != S) stop("sign_mat column count must equal number of subjects.");

  NumericVector estimate(B), se(B), stat(B), df(B), p_perm(B), p_fwer(B);
  NumericVector max_null(P);

  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (long pp = 0; pp < static_cast<long>(P); ++pp) {
    double max_abs = 0.0;
    for (uword b = 0; b < B; ++b) {
      double sum = 0.0, sumsq = 0.0;
      int n = 0;
      for (uword i = 0; i < S; ++i) {
        const double y = beta(i, b);
        if (!finite_(y)) continue;
        const double yp = static_cast<double>(sign_mat(pp, i)) * y;
        sum += yp;
        sumsq += yp * yp;
        ++n;
      }
      const double tp = t_onesample_from_sums(sum, sumsq, n);
      if (finite_(tp)) max_abs = std::max(max_abs, std::fabs(tp));
    }
    max_null[pp] = max_abs;
  }

  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (long b = 0; b < static_cast<long>(B); ++b) {
    double sum = 0.0, sumsq = 0.0;
    int n = 0;
    for (uword i = 0; i < S; ++i) {
      const double y = beta(i, b);
      if (!finite_(y)) continue;
      sum += y;
      sumsq += y * y;
      ++n;
    }
    if (n < min_subj) {
      estimate[b] = se[b] = stat[b] = df[b] = p_perm[b] = p_fwer[b] = NA_REAL;
      continue;
    }
    const double mean = mean_from_sums(sum, n);
    double var = (sumsq - static_cast<double>(n) * mean * mean) / static_cast<double>(n - 1);
    if (var < 0.0 && var > -1e-10) var = 0.0;
    const double seb = (var > 0.0 && finite_(var)) ? std::sqrt(var / static_cast<double>(n)) : NA_REAL;
    const double tobs = t_onesample_from_sums(sum, sumsq, n);
    estimate[b] = mean;
    se[b] = seb;
    stat[b] = tobs;
    df[b] = static_cast<double>(n - 1);
    if (!finite_(tobs)) {
      p_perm[b] = p_fwer[b] = NA_REAL;
      continue;
    }
    double count = 0.0, count_fwer = 0.0;
    for (uword pp = 0; pp < P; ++pp) {
      double psum = 0.0, psumsq = 0.0;
      int pn = 0;
      for (uword i = 0; i < S; ++i) {
        const double y = beta(i, b);
        if (!finite_(y)) continue;
        const double yp = static_cast<double>(sign_mat(pp, i)) * y;
        psum += yp;
        psumsq += yp * yp;
        ++pn;
      }
      const double tnull = t_onesample_from_sums(psum, psumsq, pn);
      count += perm_tail_p(tobs, tnull, tail);
      if (max_null[pp] >= std::fabs(tobs)) count_fwer += 1.0;
    }
    p_perm[b] = (count + 1.0) / (static_cast<double>(P) + 1.0);
    p_fwer[b] = (count_fwer + 1.0) / (static_cast<double>(P) + 1.0);
  }

  return Rcpp::List::create(
    _["beta_g"] = estimate,
    _["se_g"] = se,
    _["t_g"] = stat,
    _["df"] = df,
    _["p_perm"] = p_perm,
    _["p_fwer"] = p_fwer,
    _["max_t_null"] = max_null
  );
}

// group_mat is n_perm x subjects with entries 0/1.
// Observed effect is mean(group == 1) - mean(group == 0).
// [[Rcpp::export]]
Rcpp::List perm_twosample_t_cpp(const arma::mat& beta,
                                const arma::imat& group_mat,
                                const int tail = 0,
                                const bool welch = true,
                                const int min_group = 2) {
  const uword S = beta.n_rows, B = beta.n_cols;
  const uword P = group_mat.n_rows;
  if (group_mat.n_cols != S) stop("group_mat column count must equal number of subjects.");

  NumericVector estimate(B), se(B), stat(B), df(B), p_perm(B), p_fwer(B);
  NumericVector max_null(P);

  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (long pp = 0; pp < static_cast<long>(P); ++pp) {
    double max_abs = 0.0;
    for (uword b = 0; b < B; ++b) {
      double sum0 = 0.0, sum1 = 0.0, sumsq0 = 0.0, sumsq1 = 0.0;
      int n0 = 0, n1 = 0;
      for (uword i = 0; i < S; ++i) {
        const double y = beta(i, b);
        if (!finite_(y)) continue;
        if (group_mat(pp, i) == 1) {
          sum1 += y; sumsq1 += y * y; ++n1;
        } else {
          sum0 += y; sumsq0 += y * y; ++n0;
        }
      }
      const double tp = t_twosample_from_sums(sum0, sumsq0, n0, sum1, sumsq1, n1, welch);
      if (finite_(tp)) max_abs = std::max(max_abs, std::fabs(tp));
    }
    max_null[pp] = max_abs;
  }

  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (long b = 0; b < static_cast<long>(B); ++b) {
    double sum0 = 0.0, sum1 = 0.0, sumsq0 = 0.0, sumsq1 = 0.0;
    int n0 = 0, n1 = 0;
    for (uword i = 0; i < S; ++i) {
      const double y = beta(i, b);
      if (!finite_(y)) continue;
      if (group_mat(0, i) == 1) {
        sum1 += y; sumsq1 += y * y; ++n1;
      } else {
        sum0 += y; sumsq0 += y * y; ++n0;
      }
    }
    if (n0 < min_group || n1 < min_group) {
      estimate[b] = se[b] = stat[b] = df[b] = p_perm[b] = p_fwer[b] = NA_REAL;
      continue;
    }
    const double m0 = sum0 / static_cast<double>(n0);
    const double m1 = sum1 / static_cast<double>(n1);
    double v0 = (sumsq0 - static_cast<double>(n0) * m0 * m0) / static_cast<double>(n0 - 1);
    double v1 = (sumsq1 - static_cast<double>(n1) * m1 * m1) / static_cast<double>(n1 - 1);
    if (v0 < 0.0 && v0 > -1e-10) v0 = 0.0;
    if (v1 < 0.0 && v1 > -1e-10) v1 = 0.0;
    double se2;
    double dfb;
    if (welch) {
      const double a = v0 / static_cast<double>(n0);
      const double c = v1 / static_cast<double>(n1);
      se2 = a + c;
      const double den = (a * a) / static_cast<double>(n0 - 1) + (c * c) / static_cast<double>(n1 - 1);
      dfb = den > 0.0 ? (se2 * se2) / den : NA_REAL;
    } else {
      const double sp2 = ((n0 - 1.0) * v0 + (n1 - 1.0) * v1) / static_cast<double>(n0 + n1 - 2);
      se2 = sp2 * (1.0 / static_cast<double>(n0) + 1.0 / static_cast<double>(n1));
      dfb = static_cast<double>(n0 + n1 - 2);
    }
    const double tobs = t_twosample_from_sums(sum0, sumsq0, n0, sum1, sumsq1, n1, welch);
    estimate[b] = m1 - m0;
    se[b] = se2 > 0.0 ? std::sqrt(se2) : NA_REAL;
    stat[b] = tobs;
    df[b] = dfb;
    if (!finite_(tobs)) {
      p_perm[b] = p_fwer[b] = NA_REAL;
      continue;
    }
    double count = 0.0, count_fwer = 0.0;
    for (uword pp = 1; pp < P; ++pp) {
      double psum0 = 0.0, psum1 = 0.0, psumsq0 = 0.0, psumsq1 = 0.0;
      int pn0 = 0, pn1 = 0;
      for (uword i = 0; i < S; ++i) {
        const double y = beta(i, b);
        if (!finite_(y)) continue;
        if (group_mat(pp, i) == 1) {
          psum1 += y; psumsq1 += y * y; ++pn1;
        } else {
          psum0 += y; psumsq0 += y * y; ++pn0;
        }
      }
      const double tnull = t_twosample_from_sums(psum0, psumsq0, pn0, psum1, psumsq1, pn1, welch);
      count += perm_tail_p(tobs, tnull, tail);
      if (max_null[pp] >= std::fabs(tobs)) count_fwer += 1.0;
    }
    p_perm[b] = (count + 1.0) / static_cast<double>(P);
    p_fwer[b] = (count_fwer + 1.0) / static_cast<double>(P);
  }

  return Rcpp::List::create(
    _["beta_g"] = estimate,
    _["se_g"] = se,
    _["t_g"] = stat,
    _["df"] = df,
    _["p_perm"] = p_perm,
    _["p_fwer"] = p_fwer,
    _["max_t_null"] = max_null
  );
}
