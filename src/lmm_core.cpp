// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]

#include <RcppArmadillo.h>

using namespace Rcpp;
using arma::mat;
using arma::vec;
using arma::rowvec;
using arma::uword;

namespace {

inline bool finite_lmm(double x) { return std::isfinite(x); }

inline arma::mat symmetrize_pd_input(const arma::mat& x) {
  return 0.5 * (x + x.t());
}

struct LmmRiStats {
  bool ok;
  double logdetV;
  double logdetA;
  arma::mat Ainv;
  arma::mat XtVinvY;
  arma::rowvec yVinvY;
  int n_obs;
  int df_res;
};

struct LmmLowRankStats {
  bool ok;
  double logdetV;
  double logdetA;
  arma::mat Ainv;
  arma::mat XtVinvY;
  arma::rowvec yVinvY;
  int n_obs;
  int df_res;
};

struct LmmKnownVarStats {
  bool ok;
  double logdetV;
  double logdetA;
  arma::mat Ainv;
  arma::vec XtVinvY;
  double yVinvY;
  int n_obs;
  int df_res;
};

LmmRiStats lmm_ri_stats(const arma::mat& Y, const arma::mat& X, const int n_repeat, const double lambda) {
  const uword n_obs = Y.n_rows;
  const uword B = Y.n_cols;
  const uword p = X.n_cols;
  if (n_repeat <= 0 || (n_obs % static_cast<uword>(n_repeat)) != 0) {
    return {false, NA_REAL, NA_REAL, arma::mat(), arma::mat(), arma::rowvec(), static_cast<int>(n_obs), NA_INTEGER};
  }
  const uword n_subject = n_obs / static_cast<uword>(n_repeat);
  const double lam = std::max(lambda, 0.0);
  const double a = (lam <= 0.0) ? 0.0 : lam / (1.0 + lam * static_cast<double>(n_repeat));

  arma::mat XtVinvX(p, p, arma::fill::zeros);
  arma::mat XtVinvY(p, B, arma::fill::zeros);
  arma::rowvec yVinvY(B, arma::fill::zeros);

  for (uword s = 0; s < n_subject; ++s) {
    const uword start = s * static_cast<uword>(n_repeat);
    const uword end = start + static_cast<uword>(n_repeat) - 1;
    arma::mat Xb = X.rows(start, end);
    arma::mat Yb = Y.rows(start, end);
    arma::rowvec sumX = arma::sum(Xb, 0);
    arma::rowvec sumY = arma::sum(Yb, 0);
    arma::mat VinvXb = Xb.each_row() - a * sumX;
    arma::mat VinvYb = Yb.each_row() - a * sumY;
    XtVinvX += Xb.t() * VinvXb;
    XtVinvY += Xb.t() * VinvYb;
    yVinvY += arma::sum(Yb % VinvYb, 0);
  }

  arma::mat XtVinvX_sym = symmetrize_pd_input(XtVinvX);
  arma::mat Ainv;
  if (!arma::inv_sympd(Ainv, XtVinvX_sym)) {
    return {false, NA_REAL, NA_REAL, arma::mat(), arma::mat(), arma::rowvec(), static_cast<int>(n_obs), NA_INTEGER};
  }

  double logdetA = NA_REAL, signA = NA_REAL;
  arma::log_det(logdetA, signA, XtVinvX_sym);
  if (!finite_lmm(logdetA) || signA <= 0.0) {
    return {false, NA_REAL, NA_REAL, arma::mat(), arma::mat(), arma::rowvec(), static_cast<int>(n_obs), NA_INTEGER};
  }

  const double logdetV = static_cast<double>(n_subject) * std::log1p(lam * static_cast<double>(n_repeat));
  return {true, logdetV, logdetA, Ainv, XtVinvY, yVinvY, static_cast<int>(n_obs), static_cast<int>(n_obs - p)};
}

LmmLowRankStats lmm_lowrank_stats(const arma::mat& Y, const arma::mat& X, const arma::mat& U_block) {
  const uword n_obs = Y.n_rows;
  const uword B = Y.n_cols;
  const uword p = X.n_cols;
  const uword n_repeat = U_block.n_rows;
  const uword q = U_block.n_cols;

  if (n_repeat == 0 || (n_obs % n_repeat) != 0) {
    return {false, NA_REAL, NA_REAL, arma::mat(), arma::mat(), arma::rowvec(), static_cast<int>(n_obs), NA_INTEGER};
  }

  const uword n_subject = n_obs / n_repeat;
  arma::mat XtVinvX(p, p, arma::fill::zeros);
  arma::mat XtVinvY(p, B, arma::fill::zeros);
  arma::rowvec yVinvY(B, arma::fill::zeros);

  arma::mat Kinv;
  double logdetK = 0.0;
  if (q > 0) {
    arma::mat K = arma::eye(q, q) + U_block.t() * U_block;
    arma::mat K_sym = symmetrize_pd_input(K);
    if (!arma::inv_sympd(Kinv, K_sym)) {
      return {false, NA_REAL, NA_REAL, arma::mat(), arma::mat(), arma::rowvec(), static_cast<int>(n_obs), NA_INTEGER};
    }
    double signK = NA_REAL;
    arma::log_det(logdetK, signK, K_sym);
    if (!finite_lmm(logdetK) || signK <= 0.0) {
      return {false, NA_REAL, NA_REAL, arma::mat(), arma::mat(), arma::rowvec(), static_cast<int>(n_obs), NA_INTEGER};
    }
  }

  for (uword s = 0; s < n_subject; ++s) {
    const uword start = s * n_repeat;
    const uword end = start + n_repeat - 1;
    arma::mat Xb = X.rows(start, end);
    arma::mat Yb = Y.rows(start, end);
    arma::mat VinvXb = Xb;
    arma::mat VinvYb = Yb;

    if (q > 0) {
      arma::mat UtX = U_block.t() * Xb;
      arma::mat UtY = U_block.t() * Yb;
      VinvXb -= U_block * (Kinv * UtX);
      VinvYb -= U_block * (Kinv * UtY);
    }

    XtVinvX += Xb.t() * VinvXb;
    XtVinvY += Xb.t() * VinvYb;
    yVinvY += arma::sum(Yb % VinvYb, 0);
  }

  arma::mat XtVinvX_sym = symmetrize_pd_input(XtVinvX);
  arma::mat Ainv;
  if (!arma::inv_sympd(Ainv, XtVinvX_sym)) {
    return {false, NA_REAL, NA_REAL, arma::mat(), arma::mat(), arma::rowvec(), static_cast<int>(n_obs), NA_INTEGER};
  }

  double logdetA = NA_REAL, signA = NA_REAL;
  arma::log_det(logdetA, signA, XtVinvX_sym);
  if (!finite_lmm(logdetA) || signA <= 0.0) {
    return {false, NA_REAL, NA_REAL, arma::mat(), arma::mat(), arma::rowvec(), static_cast<int>(n_obs), NA_INTEGER};
  }

  const double logdetV = static_cast<double>(n_subject) * logdetK;
  return {true, logdetV, logdetA, Ainv, XtVinvY, yVinvY, static_cast<int>(n_obs), static_cast<int>(n_obs - p)};
}

LmmKnownVarStats lmm_knownvar_stats(const arma::vec& y,
                                    const arma::vec& sampling_var,
                                    const arma::mat& X,
                                    const arma::mat& U_block,
                                    const double residual_var) {
  const uword n_obs = y.n_elem;
  const uword p = X.n_cols;
  const uword n_repeat = U_block.n_rows;
  const uword q = U_block.n_cols;

  if (sampling_var.n_elem != n_obs || X.n_rows != n_obs || n_repeat == 0 ||
      (n_obs % n_repeat) != 0 || !sampling_var.is_finite() ||
      arma::any(sampling_var <= 0.0) || !finite_lmm(residual_var) ||
      residual_var < 0.0) {
    return {false, NA_REAL, NA_REAL, arma::mat(), arma::vec(), NA_REAL,
            static_cast<int>(n_obs), NA_INTEGER};
  }

  const uword n_subject = n_obs / n_repeat;
  arma::mat XtVinvX(p, p, arma::fill::zeros);
  arma::vec XtVinvY(p, arma::fill::zeros);
  double yVinvY = 0.0;
  double logdetV = 0.0;

  for (uword s = 0; s < n_subject; ++s) {
    const uword start = s * n_repeat;
    const uword end = start + n_repeat - 1;
    arma::mat Xb = X.rows(start, end);
    arma::vec yb = y.subvec(start, end);
    arma::vec db = sampling_var.subvec(start, end) + residual_var;
    if (!db.is_finite() || arma::any(db <= 0.0)) {
      return {false, NA_REAL, NA_REAL, arma::mat(), arma::vec(), NA_REAL,
              static_cast<int>(n_obs), NA_INTEGER};
    }

    arma::vec dinv = 1.0 / db;
    arma::mat VinvXb = Xb.each_col() % dinv;
    arma::vec Vinvyb = yb % dinv;
    logdetV += arma::accu(arma::log(db));

    if (q > 0) {
      arma::mat DinvU = U_block.each_col() % dinv;
      arma::mat K = arma::eye(q, q) + U_block.t() * DinvU;
      arma::mat K_sym = symmetrize_pd_input(K);
      arma::mat Kinv;
      if (!arma::inv_sympd(Kinv, K_sym)) {
        return {false, NA_REAL, NA_REAL, arma::mat(), arma::vec(), NA_REAL,
                static_cast<int>(n_obs), NA_INTEGER};
      }
      double logdetK = NA_REAL, signK = NA_REAL;
      arma::log_det(logdetK, signK, K_sym);
      if (!finite_lmm(logdetK) || signK <= 0.0) {
        return {false, NA_REAL, NA_REAL, arma::mat(), arma::vec(), NA_REAL,
                static_cast<int>(n_obs), NA_INTEGER};
      }
      VinvXb -= DinvU * (Kinv * (U_block.t() * VinvXb));
      Vinvyb -= DinvU * (Kinv * (U_block.t() * Vinvyb));
      logdetV += logdetK;
    }

    XtVinvX += Xb.t() * VinvXb;
    XtVinvY += Xb.t() * Vinvyb;
    yVinvY += arma::dot(yb, Vinvyb);
  }

  arma::mat XtVinvX_sym = symmetrize_pd_input(XtVinvX);
  arma::mat Ainv;
  if (!arma::inv_sympd(Ainv, XtVinvX_sym)) {
    return {false, NA_REAL, NA_REAL, arma::mat(), arma::vec(), NA_REAL,
            static_cast<int>(n_obs), NA_INTEGER};
  }
  double logdetA = NA_REAL, signA = NA_REAL;
  arma::log_det(logdetA, signA, XtVinvX_sym);
  if (!finite_lmm(logdetA) || signA <= 0.0 || !finite_lmm(logdetV)) {
    return {false, NA_REAL, NA_REAL, arma::mat(), arma::vec(), NA_REAL,
            static_cast<int>(n_obs), NA_INTEGER};
  }

  return {true, logdetV, logdetA, Ainv, XtVinvY, yVinvY,
          static_cast<int>(n_obs), static_cast<int>(n_obs - p)};
}

} // namespace

// [[Rcpp::export]]
double lmm_ri_objective_cpp(const arma::mat& Y,
                            const arma::mat& X,
                            const int n_repeat,
                            const double lambda,
                            const std::string fit = "REML") {
  LmmRiStats stats = lmm_ri_stats(Y, X, n_repeat, lambda);
  if (!stats.ok) return R_PosInf;

  arma::mat beta = stats.Ainv * stats.XtVinvY;
  arma::rowvec quad = arma::sum(beta % stats.XtVinvY, 0);
  arma::rowvec rss = stats.yVinvY - quad;
  const double eps = 1e-12;

  if (fit == "ML") {
    const double n = static_cast<double>(stats.n_obs);
    rss.transform([&](double x) { return std::max(x / n, eps); });
    return static_cast<double>(Y.n_cols) * stats.logdetV + n * arma::accu(arma::log(rss));
  }

  const double df = static_cast<double>(std::max(stats.df_res, 1));
  rss.transform([&](double x) { return std::max(x / df, eps); });
  return static_cast<double>(Y.n_cols) * (stats.logdetV + stats.logdetA) + df * arma::accu(arma::log(rss));
}

// [[Rcpp::export]]
Rcpp::List lmm_ri_fit_cpp(const arma::mat& Y,
                          const arma::mat& X,
                          const int n_repeat,
                          const double lambda,
                          const std::string fit = "REML") {
  LmmRiStats stats = lmm_ri_stats(Y, X, n_repeat, lambda);
  if (!stats.ok) stop("lmm_ri_fit_cpp failed: design matrix is singular or invalid");

  const uword B = Y.n_cols;
  const uword p = X.n_cols;
  arma::mat coef = stats.Ainv * stats.XtVinvY;
  arma::rowvec quad = arma::sum(coef % stats.XtVinvY, 0);
  arma::rowvec rss = stats.yVinvY - quad;
  const double scale_df = (fit == "ML") ? static_cast<double>(stats.n_obs) : static_cast<double>(std::max(stats.df_res, 1));
  const double ll_df = (fit == "ML") ? static_cast<double>(stats.n_obs) : static_cast<double>(std::max(stats.df_res, 1));
  const double eps = 1e-12;

  arma::rowvec sigma2 = rss;
  sigma2.transform([&](double x) { return std::max(x / scale_df, eps); });

  arma::mat se(p, B, arma::fill::zeros);
  arma::vec diagA = stats.Ainv.diag();
  for (uword b = 0; b < B; ++b) {
    for (uword j = 0; j < p; ++j) {
      se(j, b) = std::sqrt(std::max(diagA(j) * sigma2(b), 0.0));
    }
  }

  NumericVector df_res(B), vc_intercept(B), vc_resid(B), logLik(B), lambda_out(B), converged(B);
  const double log_2pi = std::log(2.0 * arma::datum::pi);
  for (uword b = 0; b < B; ++b) {
    const double s2 = sigma2(b);
    df_res[b] = stats.df_res;
    vc_resid[b] = s2;
    vc_intercept[b] = std::max(lambda, 0.0) * s2;
    lambda_out[b] = std::max(lambda, 0.0);
    converged[b] = 1.0;
    if (fit == "ML") {
      logLik[b] = -0.5 * (ll_df * (log_2pi + std::log(s2) + 1.0) + stats.logdetV);
    } else {
      logLik[b] = -0.5 * (ll_df * (log_2pi + std::log(s2) + 1.0) + stats.logdetV + stats.logdetA);
    }
  }

  return Rcpp::List::create(
    _["coef"] = coef,
    _["se_coef"] = se,
    _["sigma2"] = sigma2,
    _["vc_intercept"] = vc_intercept,
    _["vc_resid"] = vc_resid,
    _["df_res"] = df_res,
    _["logLik"] = logLik,
    _["lambda"] = lambda_out,
    _["converged"] = converged
  );
}

// [[Rcpp::export]]
double lmm_lowrank_objective_cpp(const arma::mat& Y,
                                 const arma::mat& X,
                                 const arma::mat& U_block,
                                 const std::string fit = "REML") {
  LmmLowRankStats stats = lmm_lowrank_stats(Y, X, U_block);
  if (!stats.ok) return R_PosInf;

  arma::mat beta = stats.Ainv * stats.XtVinvY;
  arma::rowvec quad = arma::sum(beta % stats.XtVinvY, 0);
  arma::rowvec rss = stats.yVinvY - quad;
  const double eps = 1e-12;

  if (fit == "ML") {
    const double n = static_cast<double>(stats.n_obs);
    rss.transform([&](double x) { return std::max(x / n, eps); });
    return static_cast<double>(Y.n_cols) * stats.logdetV + n * arma::accu(arma::log(rss));
  }

  const double df = static_cast<double>(std::max(stats.df_res, 1));
  rss.transform([&](double x) { return std::max(x / df, eps); });
  return static_cast<double>(Y.n_cols) * (stats.logdetV + stats.logdetA) + df * arma::accu(arma::log(rss));
}

// [[Rcpp::export]]
Rcpp::List lmm_lowrank_fit_cpp(const arma::mat& Y,
                               const arma::mat& X,
                               const arma::mat& U_block,
                               const std::string fit = "REML") {
  LmmLowRankStats stats = lmm_lowrank_stats(Y, X, U_block);
  if (!stats.ok) stop("lmm_lowrank_fit_cpp failed: design matrix is singular or invalid");

  const uword B = Y.n_cols;
  const uword p = X.n_cols;
  arma::mat coef = stats.Ainv * stats.XtVinvY;
  arma::rowvec quad = arma::sum(coef % stats.XtVinvY, 0);
  arma::rowvec rss = stats.yVinvY - quad;
  const double scale_df = (fit == "ML") ? static_cast<double>(stats.n_obs) : static_cast<double>(std::max(stats.df_res, 1));
  const double ll_df = (fit == "ML") ? static_cast<double>(stats.n_obs) : static_cast<double>(std::max(stats.df_res, 1));
  const double eps = 1e-12;

  arma::rowvec sigma2 = rss;
  sigma2.transform([&](double x) { return std::max(x / scale_df, eps); });

  arma::mat se(p, B, arma::fill::zeros);
  arma::vec diagA = stats.Ainv.diag();
  for (uword b = 0; b < B; ++b) {
    for (uword j = 0; j < p; ++j) {
      se(j, b) = std::sqrt(std::max(diagA(j) * sigma2(b), 0.0));
    }
  }

  NumericVector df_res(B), logLik(B), converged(B);
  const double log_2pi = std::log(2.0 * arma::datum::pi);
  for (uword b = 0; b < B; ++b) {
    const double s2 = sigma2(b);
    df_res[b] = stats.df_res;
    converged[b] = 1.0;
    if (fit == "ML") {
      logLik[b] = -0.5 * (ll_df * (log_2pi + std::log(s2) + 1.0) + stats.logdetV);
    } else {
      logLik[b] = -0.5 * (ll_df * (log_2pi + std::log(s2) + 1.0) + stats.logdetV + stats.logdetA);
    }
  }

  return Rcpp::List::create(
    _["coef"] = coef,
    _["se_coef"] = se,
    _["sigma2"] = sigma2,
    _["df_res"] = df_res,
    _["logLik"] = logLik,
    _["converged"] = converged
  );
}

// [[Rcpp::export]]
double lmm_knownvar_objective_cpp(const arma::vec& y,
                                  const arma::vec& sampling_var,
                                  const arma::mat& X,
                                  const arma::mat& U_block,
                                  const double residual_var,
                                  const std::string fit = "REML") {
  LmmKnownVarStats stats = lmm_knownvar_stats(
    y, sampling_var, X, U_block, residual_var
  );
  if (!stats.ok) return R_PosInf;

  arma::vec beta = stats.Ainv * stats.XtVinvY;
  const double rss = std::max(
    stats.yVinvY - arma::dot(beta, stats.XtVinvY), 0.0
  );
  if (fit == "ML") return stats.logdetV + rss;
  return stats.logdetV + stats.logdetA + rss;
}

// [[Rcpp::export]]
Rcpp::List lmm_knownvar_fit_cpp(const arma::vec& y,
                                const arma::vec& sampling_var,
                                const arma::mat& X,
                                const arma::mat& U_block,
                                const double residual_var,
                                const std::string fit = "REML") {
  LmmKnownVarStats stats = lmm_knownvar_stats(
    y, sampling_var, X, U_block, residual_var
  );
  if (!stats.ok) {
    stop("lmm_knownvar_fit_cpp failed: covariance or design matrix is invalid");
  }

  arma::vec coef = stats.Ainv * stats.XtVinvY;
  const double rss = std::max(
    stats.yVinvY - arma::dot(coef, stats.XtVinvY), 0.0
  );
  arma::vec se = arma::sqrt(arma::clamp(stats.Ainv.diag(), 0.0, arma::datum::inf));
  const double log_2pi = std::log(2.0 * arma::datum::pi);
  const double ll = (fit == "ML")
    ? -0.5 * (static_cast<double>(stats.n_obs) * log_2pi + stats.logdetV + rss)
    : -0.5 * (static_cast<double>(stats.df_res) * log_2pi + stats.logdetV + stats.logdetA + rss);

  return Rcpp::List::create(
    _["coef"] = coef,
    _["se_coef"] = se,
    _["df_res"] = stats.df_res,
    _["logLik"] = ll
  );
}
