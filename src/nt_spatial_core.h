// nt_spatial_core.h -- vendored from neurothresh (inst/include/neurothresh/spatial_core.h,
// neurothresh commit 40ae4f3). Do not edit here; change neurothresh and re-copy.
// fmrigds tests check agreement with neurothresh when it is installed.
//
// Header-only C++ core for spatial permutation statistics on masked 3D
// volumes: neighbour graphs (6/18/26 connectivity), threshold-free cluster
// enhancement (TFCE) and cluster labelling with extent and mass.
//
// The functions depend only on the C++ standard library (no R or Rcpp
// headers) and allocate no R memory, so they can be called from OpenMP
// worker threads. Each thread must own its workspace objects.
//
// Semantics are those of the reference R implementation in neurothresh
// (`.tfce_from_vec_r()`, `.find_clusters_r()`):
//
// * TFCE heights are h_i = min(i * dh, h_max) for i = 1, ..., n with
//   n = (int)(h_max / dh + 1e-10), exactly as `seq(0, h_max, by = dh)`
//   minus the h = 0 step. A voxel with value z receives
//   sum_{i : h_i <= z} e_i(v)^E * h_i^H * dh, where e_i(v) is the number of
//   voxels in the connected component of {z >= h_i} that contains v.
// * Non-finite values never exceed a threshold (they act as 0).
// * Cluster labels are numbered 1..K in order of each cluster's smallest
//   node index (node order is the linear voxel order of the mask), matching
//   the breadth-first R labeller.
//
// fmrigds vendors a copy of this file (src/nt_spatial_core.h); keep the two
// in sync. The SPATIAL_CORE_VERSION macro records the revision.

#ifndef NEUROTHRESH_SPATIAL_CORE_H
#define NEUROTHRESH_SPATIAL_CORE_H

#define NEUROTHRESH_SPATIAL_CORE_VERSION 1

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace ntcore {

// Compressed sparse row neighbour graph over n nodes.
struct Graph {
  int n = 0;
  std::vector<int> ptr;  // length n + 1
  std::vector<int> idx;  // neighbour node indices
};

// Build the neighbour graph of the in-mask voxels of a 3D grid.
//
// dims: grid extents (x fastest, column-major as in R).
// voxels: 0-based linear voxel index of each node (node order defines the
//   label order of cluster_label()); must be unique and in range.
// connectivity: 6 (faces), 18 (faces + edges) or 26 (faces + edges + corners).
inline Graph build_grid_graph(const int dims[3], const std::vector<int>& voxels,
                              int connectivity) {
  if (connectivity != 6 && connectivity != 18 && connectivity != 26) {
    throw std::invalid_argument("connectivity must be 6, 18 or 26");
  }
  const long nx = dims[0], ny = dims[1], nz = dims[2];
  const long nvox = nx * ny * nz;
  const int n = static_cast<int>(voxels.size());
  std::vector<int> node_of(static_cast<std::size_t>(nvox), -1);
  for (int v = 0; v < n; ++v) {
    const long lin = voxels[v];
    if (lin < 0 || lin >= nvox) throw std::out_of_range("voxel index out of range");
    if (node_of[lin] >= 0) throw std::invalid_argument("duplicate voxel index");
    node_of[lin] = v;
  }

  int off[26][3];
  int n_off = 0;
  for (int dk = -1; dk <= 1; ++dk) {
    for (int dj = -1; dj <= 1; ++dj) {
      for (int di = -1; di <= 1; ++di) {
        const int m = std::abs(di) + std::abs(dj) + std::abs(dk);
        if (m == 0) continue;
        if (connectivity == 6 && m > 1) continue;
        if (connectivity == 18 && m > 2) continue;
        off[n_off][0] = di;
        off[n_off][1] = dj;
        off[n_off][2] = dk;
        ++n_off;
      }
    }
  }

  Graph g;
  g.n = n;
  g.ptr.assign(static_cast<std::size_t>(n) + 1, 0);
  g.idx.reserve(static_cast<std::size_t>(n) * static_cast<std::size_t>(n_off));
  for (int v = 0; v < n; ++v) {
    const long lin = voxels[v];
    const long i = lin % nx;
    const long j = (lin / nx) % ny;
    const long k = lin / (nx * ny);
    for (int o = 0; o < n_off; ++o) {
      const long ii = i + off[o][0], jj = j + off[o][1], kk = k + off[o][2];
      if (ii < 0 || ii >= nx || jj < 0 || jj >= ny || kk < 0 || kk >= nz) continue;
      const int u = node_of[ii + nx * (jj + ny * kk)];
      if (u >= 0) g.idx.push_back(u);
    }
    g.ptr[v + 1] = static_cast<int>(g.idx.size());
  }
  return g;
}

// Union-find with union by size and path halving.
struct UnionFind {
  std::vector<int> parent;
  std::vector<int> size;
  void resize(int n) {
    parent.resize(n);
    size.resize(n);
  }
  int find(int x) {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  }
  void make(int x) {
    parent[x] = x;
    size[x] = 1;
  }
  int unite(int a, int b) {
    a = find(a);
    b = find(b);
    if (a == b) return a;
    if (size[a] < size[b]) std::swap(a, b);
    parent[b] = a;
    size[a] += size[b];
    return a;
  }
};

// Per-thread scratch space for tfce().
struct TfceWork {
  UnionFind uf;
  std::vector<int> level;   // highest height index reached by each node
  std::vector<int> order;   // nodes sorted by level (descending)
  std::vector<int> start;   // bucket offsets
  std::vector<int> stamp;   // per-root cache stamp
  std::vector<double> cval; // per-root cached contribution
  std::vector<char> on;     // node already added
  void resize(int n) {
    uf.resize(n);
    level.resize(n);
    order.resize(n);
    stamp.resize(n);
    cval.resize(n);
    on.resize(n);
  }
};

inline double tfce_height(int i, double dh, double hmax) {
  const double h = static_cast<double>(i) * dh;
  return h > hmax ? hmax : h;
}

// Positive-tail TFCE of x (length g.n) written to out (length g.n). To get
// the negative tail pass -x. Returns the maximum TFCE value (0 when no node
// is positive).
inline double tfce(const Graph& g, const double* x, double H, double E, double dh,
                   double* out, TfceWork& w) {
  const int n = g.n;
  std::fill(out, out + n, 0.0);
  if (n == 0) return 0.0;
  double hmax = 0.0;
  for (int v = 0; v < n; ++v) {
    if (std::isfinite(x[v]) && x[v] > hmax) hmax = x[v];
  }
  if (!(hmax > 0.0)) return 0.0;
  const int nlev = static_cast<int>(hmax / dh + 1e-10);
  if (nlev < 1) return 0.0;
  if (static_cast<int>(w.level.size()) < n) w.resize(n);

  // Highest height index each node reaches (0 = none).
  w.start.assign(static_cast<std::size_t>(nlev) + 2, 0);
  for (int v = 0; v < n; ++v) {
    const double z = x[v];
    int k = 0;
    if (std::isfinite(z) && z > 0.0) {
      double kd = std::floor(z / dh);
      k = kd > nlev ? nlev : static_cast<int>(kd);
      while (k < nlev && tfce_height(k + 1, dh, hmax) <= z) ++k;
      while (k > 0 && tfce_height(k, dh, hmax) > z) --k;
    }
    w.level[v] = k;
    w.on[v] = 0;
    w.stamp[v] = 0;
    ++w.start[k];
  }
  // Bucket offsets so that order[] lists nodes by descending level.
  {
    int pos = 0;
    for (int k = nlev; k >= 0; --k) {
      const int c = w.start[k];
      w.start[k] = pos;
      pos += c;
    }
    w.start[nlev + 1] = 0;  // unused sentinel
  }
  {
    std::vector<int> fill(w.start.begin(), w.start.end());
    for (int v = 0; v < n; ++v) w.order[fill[w.level[v]]++] = v;
  }

  // Sweep heights from the top down, adding nodes and merging components.
  int n_active = 0;
  for (int i = nlev; i >= 1; --i) {
    // Nodes of bucket i occupy order[start[i], start[i - 1]).
    const int b0 = w.start[i];
    const int end_i = w.start[i - 1];
    for (int p = b0; p < end_i; ++p) {
      const int v = w.order[p];
      w.uf.make(v);
      w.on[v] = 1;
      for (int q = g.ptr[v]; q < g.ptr[v + 1]; ++q) {
        const int u = g.idx[q];
        if (w.on[u]) w.uf.unite(v, u);
      }
    }
    n_active = end_i;
    const double h = tfce_height(i, dh, hmax);
    const double hH = std::pow(h, H);
    for (int p = 0; p < n_active; ++p) {
      const int v = w.order[p];
      const int r = w.uf.find(v);
      if (w.stamp[r] != i) {
        w.stamp[r] = i;
        w.cval[r] = std::pow(static_cast<double>(w.uf.size[r]), E) * hH * dh;
      }
      out[v] += w.cval[r];
    }
  }
  double mx = 0.0;
  for (int v = 0; v < n; ++v) {
    if (out[v] > mx) mx = out[v];
  }
  return mx;
}

// Per-thread scratch space for cluster_label().
struct ClusterWork {
  UnionFind uf;
  std::vector<char> on;
  std::vector<int> root_label;
  void resize(int n) {
    uf.resize(n);
    on.resize(n);
    root_label.resize(n);
  }
};

struct ClusterSummary {
  int n_clusters = 0;
  int max_size = 0;
  double max_mass = 0.0;
};

// Label the connected components of {v : x[v] >= u} (u finite).
//
// labels (optional, length g.n): 1..K in order of each cluster's smallest
//   node index, 0 outside clusters.
// sizes, mass, peak (optional): per-cluster voxel count, excess mass
//   sum(x - u), and node index of the maximum (first node on ties).
// Returns the number of clusters plus the largest size and mass.
inline ClusterSummary cluster_label(const Graph& g, const double* x, double u,
                                    ClusterWork& w, int* labels = nullptr,
                                    std::vector<int>* sizes = nullptr,
                                    std::vector<double>* mass = nullptr,
                                    std::vector<int>* peak = nullptr) {
  const int n = g.n;
  ClusterSummary s;
  if (static_cast<int>(w.on.size()) < n) w.resize(n);
  for (int v = 0; v < n; ++v) {
    w.on[v] = (std::isfinite(x[v]) && x[v] >= u) ? 1 : 0;
  }
  for (int v = 0; v < n; ++v) {
    if (!w.on[v]) continue;
    w.uf.make(v);
    for (int q = g.ptr[v]; q < g.ptr[v + 1]; ++q) {
      const int nb = g.idx[q];
      if (nb < v && w.on[nb]) w.uf.unite(v, nb);
    }
  }
  std::vector<int> lsize;
  std::vector<double> lmass;
  std::vector<int> lpeak;
  int K = 0;
  // Assign labels in order of smallest member: scan nodes in order and give
  // a new label the first time a root is met.
  for (int v = 0; v < n; ++v) {
    if (w.on[v]) w.root_label[w.uf.find(v)] = 0;
  }
  for (int v = 0; v < n; ++v) {
    if (!w.on[v]) {
      if (labels) labels[v] = 0;
      continue;
    }
    const int r = w.uf.find(v);
    int lab = w.root_label[r];
    if (lab == 0) {
      lab = ++K;
      w.root_label[r] = lab;
      lsize.push_back(0);
      lmass.push_back(0.0);
      lpeak.push_back(v);
    }
    if (labels) labels[v] = lab;
    lsize[lab - 1] += 1;
    lmass[lab - 1] += x[v] - u;
    if (x[v] > x[lpeak[lab - 1]]) lpeak[lab - 1] = v;
  }
  s.n_clusters = K;
  for (int c = 0; c < K; ++c) {
    if (lsize[c] > s.max_size) s.max_size = lsize[c];
    if (lmass[c] > s.max_mass) s.max_mass = lmass[c];
  }
  if (sizes) sizes->swap(lsize);
  if (mass) mass->swap(lmass);
  if (peak) peak->swap(lpeak);
  return s;
}

}  // namespace ntcore

#endif  // NEUROTHRESH_SPATIAL_CORE_H
