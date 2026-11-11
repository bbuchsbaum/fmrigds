// Subject-level sign-flip permutation inference with spatial statistics.
//
// One sign matrix drives every family-wise null in a single pass over the
// permutations: max-|t| (or max-t / max-(-t)), max-TFCE and the largest
// cluster mass and extent. The observed statistic is always computed from
// the unflipped data; when skip_first is true, row 1 of the sign matrix is
// the identity flip and is excluded from the null set.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

#include <algorithm>
#include <cmath>
#include <vector>

#include "nt_spatial_core.h"

using namespace Rcpp;
using arma::uword;

namespace {

// Same arithmetic as t_onesample_from_sums() in reducers_core.cpp, so that
// the t statistics (and hence p_perm / p_fwer) agree bit for bit with
// perm_onesample_t_cpp().
inline double sp_t_from_sums(const double sum, const double sumsq, const int n) {
  if (n < 2) return NA_REAL;
  const double mean = sum / static_cast<double>(n);
  double var = (sumsq - static_cast<double>(n) * mean * mean) / static_cast<double>(n - 1);
  if (var < 0.0 && var > -1e-10) var = 0.0;
  if (var <= 0.0 || !std::isfinite(var)) return NA_REAL;
  return mean / std::sqrt(var / static_cast<double>(n));
}

inline double sp_tail_extreme(const double t, const int tail) {
  if (tail == 1) return -t;     // less
  if (tail == 2) return t;      // greater
  return std::fabs(t);          // two.sided
}

inline double sp_tail_hit(const double obs, const double null_stat, const int tail) {
  if (!std::isfinite(obs) || !std::isfinite(null_stat)) return 0.0;
  if (tail == 1) return null_stat <= obs ? 1.0 : 0.0;
  if (tail == 2) return null_stat >= obs ? 1.0 : 0.0;
  return std::fabs(null_stat) >= std::fabs(obs) ? 1.0 : 0.0;
}

// Signed normal equivalent of t with df degrees of freedom (same one-tail
// probability), on the log scale so that it stays finite far in the tail.
// Matches neurothresh's .ss_t_to_z() and fmrigds' .z_from_t().
inline double sp_z_from_t(const double t, const double df) {
  if (!std::isfinite(t) || !(df > 0.0)) return NA_REAL;
  if (t == 0.0) return 0.0;
  const double lp = R::pt(std::fabs(t), df, /*lower*/ 0, /*log*/ 1);
  const double z = R::qnorm(lp, 0.0, 1.0, /*lower*/ 0, /*log*/ 1);
  return t > 0.0 ? z : -z;
}

struct Settings {
  int tail;
  int min_subj;
  bool do_tfce;
  bool do_cluster;
  bool use_z;
  double H, E, dh, u;
};

// Statistic maps for one sign vector (signs == nullptr: observed data).
// t: one-sample t per sample (NA if fewer than min_subj finite subjects or
// zero variance); x: spatial statistic (t or z), NaN where t is NA.
//
// For one-sided null maps (full == false) values on the untested side only
// need to stay on that side (they never reach a positive threshold), so the
// t-to-z conversion is skipped there.
void stat_maps(const arma::mat& beta, const int* signs, const uword P_stride,
               const Settings& s, double* t, double* x, const bool full = true) {
  const uword S = beta.n_rows, V = beta.n_cols;
  for (uword b = 0; b < V; ++b) {
    const double* col = beta.colptr(b);
    double sum = 0.0, sumsq = 0.0;
    int n = 0;
    for (uword i = 0; i < S; ++i) {
      const double y = col[i];
      if (!std::isfinite(y)) continue;
      const double yp = signs ? static_cast<double>(signs[i * P_stride]) * y : y;
      sum += yp;
      sumsq += yp * yp;
      ++n;
    }
    double tv = NA_REAL;
    if (n >= s.min_subj) tv = sp_t_from_sums(sum, sumsq, n);
    t[b] = tv;
    if (!std::isfinite(tv)) {
      x[b] = R_NaN;
    } else if (s.use_z) {
      if (!full && ((s.tail == 2 && tv <= 0.0) || (s.tail == 1 && tv >= 0.0))) {
        x[b] = tv;
      } else {
        x[b] = sp_z_from_t(tv, static_cast<double>(n - 1));
      }
    } else {
      x[b] = tv;
    }
  }
}

struct Work {
  ntcore::TfceWork tw;
  ntcore::ClusterWork cw;
  std::vector<double> neg, tf;
  void resize(int n) {
    tw.resize(n);
    cw.resize(n);
    neg.resize(n);
    tf.resize(n);
  }
};

// Family-wise null statistics of one statistic map x.
void spatial_max(const ntcore::Graph& g, const double* x, const Settings& s,
                 Work& w, double& tfce_max, double& mass_max, double& extent_max) {
  const int n = g.n;
  tfce_max = 0.0;
  mass_max = 0.0;
  extent_max = 0.0;
  const bool pos = (s.tail == 2 || s.tail == 0);
  const bool negt = (s.tail == 1 || s.tail == 0);
  if (negt) {
    for (int v = 0; v < n; ++v) w.neg[v] = -x[v];
  }
  if (s.do_tfce) {
    if (pos) tfce_max = std::max(tfce_max, ntcore::tfce(g, x, s.H, s.E, s.dh, w.tf.data(), w.tw));
    if (negt) tfce_max = std::max(tfce_max, ntcore::tfce(g, w.neg.data(), s.H, s.E, s.dh, w.tf.data(), w.tw));
  }
  if (s.do_cluster) {
    if (pos) {
      ntcore::ClusterSummary cs = ntcore::cluster_label(g, x, s.u, w.cw);
      mass_max = std::max(mass_max, cs.max_mass);
      extent_max = std::max(extent_max, static_cast<double>(cs.max_size));
    }
    if (negt) {
      ntcore::ClusterSummary cs = ntcore::cluster_label(g, w.neg.data(), s.u, w.cw);
      mass_max = std::max(mass_max, cs.max_mass);
      extent_max = std::max(extent_max, static_cast<double>(cs.max_size));
    }
  }
}

List cluster_list(const ntcore::Graph& g, const double* x, double u, Work& w) {
  IntegerVector labels(g.n);
  std::vector<int> sizes, peak;
  std::vector<double> mass;
  ntcore::cluster_label(g, x, u, w.cw, labels.begin(), &sizes, &mass, &peak);
  IntegerVector peak1(peak.size());
  for (std::size_t c = 0; c < peak.size(); ++c) peak1[c] = peak[c] + 1;
  return List::create(
    _["labels"] = labels,
    _["sizes"] = IntegerVector(sizes.begin(), sizes.end()),
    _["mass"] = NumericVector(mass.begin(), mass.end()),
    _["peak"] = peak1
  );
}

}  // namespace

// beta: subjects x samples; samples are in-mask voxels whose 1-based grid
//   indices are `voxels` on a grid of extent `dims`.
// sign_mat: n_perm x subjects (+1/-1).
// tail: 0 = two.sided, 1 = less, 2 = greater.
// stat_type: 0 = t, 1 = z-equivalent of t (per-sample df = n - 1).
// threads: OpenMP team size for this call (<= 0: current default).
// [[Rcpp::export]]
Rcpp::List perm_onesample_spatial_cpp(const arma::mat& beta,
                                      const arma::imat& sign_mat,
                                      const Rcpp::IntegerVector& dims,
                                      const Rcpp::IntegerVector& voxels,
                                      const int connectivity,
                                      const int tail,
                                      const int min_subj,
                                      const bool skip_first,
                                      const bool do_tfce,
                                      const bool do_cluster,
                                      const double H,
                                      const double E,
                                      const double dh,
                                      const double cluster_threshold,
                                      const int stat_type,
                                      const int threads = 0) {
  const uword S = beta.n_rows, V = beta.n_cols;
  const uword P = sign_mat.n_rows;
  if (sign_mat.n_cols != S) stop("sign_mat column count must equal number of subjects.");
  if (static_cast<uword>(voxels.size()) != V) stop("length(voxels) must equal the number of samples.");
  if (dims.size() != 3) stop("dims must have length 3.");
  if (do_tfce && (!(H > 0.0) || !(E > 0.0) || !(dh > 0.0))) stop("H, E and dh must be positive.");
  if (do_cluster && !std::isfinite(cluster_threshold)) stop("cluster threshold must be finite.");
  const uword first = (skip_first && P > 0) ? 1 : 0;
  const double n_null = static_cast<double>(P - first);

  ntcore::Graph g;
  {
    int d[3] = {dims[0], dims[1], dims[2]};
    std::vector<int> vox(V);
    for (uword v = 0; v < V; ++v) vox[v] = voxels[v] - 1;
    try {
      g = ntcore::build_grid_graph(d, vox, connectivity);
    } catch (const std::exception& e) {
      stop(e.what());
    }
  }

  Settings s;
  s.tail = tail;
  s.min_subj = min_subj;
  s.do_tfce = do_tfce;
  s.do_cluster = do_cluster;
  s.use_z = (stat_type == 1);
  s.H = H;
  s.E = E;
  s.dh = dh;
  s.u = cluster_threshold;

  // Observed maps.
  std::vector<double> t_obs(V), x_obs(V);
  stat_maps(beta, nullptr, 0, s, t_obs.data(), x_obs.data());
  Work w0;
  w0.resize(static_cast<int>(V));
  NumericVector tfce_pos(V), tfce_neg(V);
  if (do_tfce) {
    if (tail == 2 || tail == 0) ntcore::tfce(g, x_obs.data(), H, E, dh, tfce_pos.begin(), w0.tw);
    if (tail == 1 || tail == 0) {
      for (uword v = 0; v < V; ++v) w0.neg[v] = -x_obs[v];
      ntcore::tfce(g, w0.neg.data(), H, E, dh, tfce_neg.begin(), w0.tw);
    }
  }
  List clus_pos, clus_neg;
  if (do_cluster) {
    if (tail == 2 || tail == 0) clus_pos = cluster_list(g, x_obs.data(), cluster_threshold, w0);
    if (tail == 1 || tail == 0) {
      for (uword v = 0; v < V; ++v) w0.neg[v] = -x_obs[v];
      clus_neg = cluster_list(g, w0.neg.data(), cluster_threshold, w0);
    }
  }

  // Null distributions.
  NumericVector max_t(P, NA_REAL), max_tfce(P, NA_REAL), max_mass(P, NA_REAL),
      max_extent(P, NA_REAL);
  std::vector<double> count_unc(V, 0.0);
  const int* sign_ptr = sign_mat.memptr();  // column-major: (pp, i) at pp + i * P

  int n_threads = 1;
#ifdef _OPENMP
  n_threads = threads > 0 ? threads : omp_get_max_threads();
#endif
  (void)n_threads;

  #ifdef _OPENMP
  #pragma omp parallel num_threads(n_threads)
  #endif
  {
    std::vector<double> t(V), x(V), cnt(V, 0.0);
    Work w;
    w.resize(static_cast<int>(V));
    #ifdef _OPENMP
    #pragma omp for schedule(dynamic, 1)
    #endif
    for (long pp = static_cast<long>(first); pp < static_cast<long>(P); ++pp) {
      stat_maps(beta, sign_ptr + pp, P, s, t.data(), x.data(), /*full*/ false);
      double mt = R_NegInf;
      for (uword b = 0; b < V; ++b) {
        const double tv = t[b];
        if (std::isfinite(tv)) mt = std::max(mt, sp_tail_extreme(tv, tail));
        cnt[b] += sp_tail_hit(t_obs[b], tv, tail);
      }
      max_t[pp] = mt;
      if (do_tfce || do_cluster) {
        double tm = 0.0, mm = 0.0, em = 0.0;
        spatial_max(g, x.data(), s, w, tm, mm, em);
        if (do_tfce) max_tfce[pp] = tm;
        if (do_cluster) {
          max_mass[pp] = mm;
          max_extent[pp] = em;
        }
      }
    }
    #ifdef _OPENMP
    #pragma omp critical
    #endif
    {
      for (uword b = 0; b < V; ++b) count_unc[b] += cnt[b];
    }
  }

  // Voxelwise outputs.
  NumericVector estimate(V), se(V), stat(V), df(V), p_perm(V), p_fwer(V), xs(V);
  std::vector<double> null_t;
  null_t.reserve(P - first);
  for (uword pp = first; pp < P; ++pp) null_t.push_back(max_t[pp]);
  std::sort(null_t.begin(), null_t.end());
  for (uword b = 0; b < V; ++b) {
    const double* col = beta.colptr(b);
    double sum = 0.0, sumsq = 0.0;
    int n = 0;
    for (uword i = 0; i < S; ++i) {
      const double y = col[i];
      if (!std::isfinite(y)) continue;
      sum += y;
      sumsq += y * y;
      ++n;
    }
    if (n < min_subj) {
      estimate[b] = se[b] = stat[b] = df[b] = p_perm[b] = p_fwer[b] = xs[b] = NA_REAL;
      continue;
    }
    const double mean = sum / static_cast<double>(n);
    double var = (sumsq - static_cast<double>(n) * mean * mean) / static_cast<double>(n - 1);
    if (var < 0.0 && var > -1e-10) var = 0.0;
    estimate[b] = mean;
    se[b] = (var > 0.0 && std::isfinite(var)) ? std::sqrt(var / static_cast<double>(n)) : NA_REAL;
    stat[b] = t_obs[b];
    df[b] = static_cast<double>(n - 1);
    xs[b] = std::isfinite(x_obs[b]) ? x_obs[b] : NA_REAL;
    if (!std::isfinite(t_obs[b])) {
      p_perm[b] = p_fwer[b] = NA_REAL;
      continue;
    }
    const double obs_ext = sp_tail_extreme(t_obs[b], tail);
    // #{null max >= obs_ext}: the null maxima are sorted ascending.
    const auto it = std::lower_bound(null_t.begin(), null_t.end(), obs_ext);
    const double count_fwer = static_cast<double>(null_t.end() - it);
    p_perm[b] = (count_unc[b] + 1.0) / (n_null + 1.0);
    p_fwer[b] = (count_fwer + 1.0) / (n_null + 1.0);
  }

  return List::create(
    _["beta_g"] = estimate,
    _["se_g"] = se,
    _["t_g"] = stat,
    _["df"] = df,
    _["p_perm"] = p_perm,
    _["p_fwer"] = p_fwer,
    _["stat"] = xs,
    _["tfce_pos"] = tfce_pos,
    _["tfce_neg"] = tfce_neg,
    _["clusters_pos"] = clus_pos,
    _["clusters_neg"] = clus_neg,
    _["max_t_null"] = max_t,
    _["max_tfce_null"] = max_tfce,
    _["max_mass_null"] = max_mass,
    _["max_extent_null"] = max_extent
  );
}
