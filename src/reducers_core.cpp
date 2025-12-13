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
  NumericVector tau2(B), Q(B), df_res(B);

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
    arma::mat A;
    bool ok = arma::inv_sympd(A, G);
    if (!ok) { tau2[b] = Q[b] = df_res[b] = NA_REAL; continue; }
    arma::vec bh_fe = A * g;

    double Qb = 0.0, sw = 0.0, trH = 0.0;
    for (uword i = 0; i < S; ++i) {
      double w = w_save[i];
      if (w <= 0.0) continue;
      sw += w;
      double r = beta(i, b) - arma::dot(X.row(i), bh_fe);
      Qb += w * r * r;
      double hi = w * arma::as_scalar(X.row(i) * A * X.row(i).t());
      trH += hi;
    }
    double C = sw - trH;
    double dfb = static_cast<double>(k - static_cast<int>(p));
    double t2 = std::max(0.0, (Qb - dfb) / std::max(C, eps));

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
