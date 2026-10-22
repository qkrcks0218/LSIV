// siv_kliep.cpp
//
// KLIEP density-ratio engine, PORTED FROM UDID_Package/src/udid_kliep.cpp.
//
// The manuscript's own replication code estimates the density ratio of
// eq. (15),
//     R_Y(y, x) = f(y | A=0, Z=1, x) / f(y | A=0, Z=0, x),
// DIRECTLY and nonparametrically -- `Simulation_F_Both/1.Sim_Conti_SepIV.R`
// uses `densityratio::kliep()` with a cross-validated Gaussian-kernel
// bandwidth, cross-fitted, rather than dividing two separately estimated
// conditional densities.  That is the strategy Section 5.1 of the manuscript
// recommends.  This file provides the same estimator inside SIV, reusing the
// compiled engine already written for UDID rather than adding a dependency on
// `densratio` / `densityratio`.
//
// Ported verbatim except for:
//   * the UDID-specific EIF grid helpers (rowsums_grid_cpp,
//     sens_rowsums_mroot_cpp, mboot_sd_cpp) are NOT ported -- they belong to
//     the UDiD influence function, not to this package;
//   * every exported symbol is prefixed `siv_` for namespace hygiene, so that
//     both packages can be loaded at once;
//   * `kliep_optimize()` is exposed to R (`siv_kliep_optimize_cpp`) so that the
//     median-heuristic path does not need `densratio:::KLIEP_optimize_alpha`.
//
// The bandwidth search in siv_kliep_bw_cpp() IS the hyper-parameter tuning:
// coordinate descent over (sigma_Y, sigma_X) maximising the K-fold
// cross-validated KLIEP log-likelihood.
//
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cmath>

using namespace Rcpp;
using namespace arma;

//' Fast Gaussian kernel density ratio evaluation
//'
//' Computes  r(x_i) = sum_c w_c * exp(-||x_i - center_c||^2 / (2*sigma^2))
//' for every row x_i of \code{x}.
//'
//' Drop-in replacement for:
//'   as.vector(densratio:::compute_kernel_Gaussian(x, centers, sigma) %*% kernel_weights)
//'
//' @param x              (N_eval x p)    evaluation points
//' @param centers        (N_centers x p) kernel centres from KLIEP fit
//' @param sigma          scalar          bandwidth
//' @param kernel_weights (N_centers)     weights from KLIEP fit
//' @param n_threads      int             OpenMP threads (0 = all; ignored without OpenMP)
//' @return numeric vector of length N_eval
// [[Rcpp::export]]
arma::vec siv_gaussian_kernel_ratio(
    const arma::mat& x,
    const arma::mat& centers,
    double           sigma,
    const arma::vec& kernel_weights,
    int              n_threads = 0)
{
  const int    n_eval    = (int)x.n_rows;
  const int    n_centers = (int)centers.n_rows;
  const int    n_cols    = (int)x.n_cols;
  const double inv_2s2   = 1.0 / (2.0 * sigma * sigma);

  arma::vec result(n_eval, arma::fill::zeros);

#ifdef _OPENMP
  if (n_threads > 0) omp_set_num_threads(n_threads);
  #pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < n_eval; i++) {
    double r = 0.0;
    for (int c = 0; c < n_centers; c++) {
      double dist2 = 0.0;
      for (int p = 0; p < n_cols; p++) {
        double d = x(i, p) - centers(c, p);
        dist2 += d * d;
      }
      r += kernel_weights[c] * std::exp(-dist2 * inv_2s2);
    }
    result[i] = r;
  }

  return result;
}


//' Fast Gaussian kernel matrix with per-dimension bandwidths
//'
//' Computes  K[i,c] = exp( - sum_j (x[i,j] - center[c,j])^2 / (2 * sigma_j^2) )
//' for every row i and every centre c.
//'
//' @param x        (N x p) evaluation points
//' @param centers  (M x p) kernel centres
//' @param sigma_vec (p)    per-dimension bandwidths
//' @param n_threads int    OpenMP threads (0 = all; ignored without OpenMP)
//' @return (N x M) kernel matrix
// [[Rcpp::export]]
arma::mat siv_gaussian_kernel_matrix_bw(
    const arma::mat& x,
    const arma::mat& centers,
    const arma::vec& sigma_vec,
    int              n_threads = 0)
{
  const int n_eval    = (int)x.n_rows;
  const int n_centers = (int)centers.n_rows;
  const int n_cols    = (int)x.n_cols;

  // precompute 1/(2*sigma_j^2) for each dimension
  arma::vec inv_2s2(n_cols);
  for (int p = 0; p < n_cols; p++) {
    inv_2s2[p] = 1.0 / (2.0 * sigma_vec[p] * sigma_vec[p]);
  }

  arma::mat result(n_eval, n_centers);

#ifdef _OPENMP
  if (n_threads > 0) omp_set_num_threads(n_threads);
  #pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < n_eval; i++) {
    for (int c = 0; c < n_centers; c++) {
      double dist2 = 0.0;
      for (int p = 0; p < n_cols; p++) {
        double d = x(i, p) - centers(c, p);
        dist2 += d * d * inv_2s2[p];
      }
      result(i, c) = std::exp(-dist2);
    }
  }

  return result;
}


//' Fast Gaussian kernel density ratio evaluation with per-dimension bandwidths
//'
//' Computes  r(x_i) = sum_c w_c * exp( - sum_j (x[i,j] - center[c,j])^2 / (2 sigma_j^2) )
//'
//' @param x              (N_eval x p)    evaluation points
//' @param centers        (N_centers x p) kernel centres
//' @param sigma_vec      (p)             per-dimension bandwidths
//' @param kernel_weights (N_centers)     weights from KLIEP fit
//' @param n_threads      int             OpenMP threads (0 = all; ignored without OpenMP)
//' @return numeric vector of length N_eval
// [[Rcpp::export]]
arma::vec siv_gaussian_kernel_ratio_bw(
    const arma::mat& x,
    const arma::mat& centers,
    const arma::vec& sigma_vec,
    const arma::vec& kernel_weights,
    int              n_threads = 0)
{
  const int n_eval    = (int)x.n_rows;
  const int n_centers = (int)centers.n_rows;
  const int n_cols    = (int)x.n_cols;

  arma::vec inv_2s2(n_cols);
  for (int p = 0; p < n_cols; p++) {
    inv_2s2[p] = 1.0 / (2.0 * sigma_vec[p] * sigma_vec[p]);
  }

  arma::vec result(n_eval, arma::fill::zeros);

#ifdef _OPENMP
  if (n_threads > 0) omp_set_num_threads(n_threads);
  #pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < n_eval; i++) {
    double r = 0.0;
    for (int c = 0; c < n_centers; c++) {
      double dist2 = 0.0;
      for (int p = 0; p < n_cols; p++) {
        double d = x(i, p) - centers(c, p);
        dist2 += d * d * inv_2s2[p];
      }
      r += kernel_weights[c] * std::exp(-dist2);
    }
    result[i] = r;
  }

  return result;
}


// ================================================================
//  Full KLIEP implementation in C++ for speed
// ================================================================

// Compute the kernel matrix (N x M) with per-dimension bandwidths
static arma::mat compute_phi(const arma::mat& x,
                              const arma::mat& centers,
                              const arma::vec& inv_2s2)
{
  const int N = (int)x.n_rows;
  const int M = (int)centers.n_rows;
  const int p = (int)x.n_cols;
  arma::mat phi(N, M);

  for (int i = 0; i < N; i++) {
    for (int c = 0; c < M; c++) {
      double d2 = 0.0;
      for (int j = 0; j < p; j++) {
        double d = x(i, j) - centers(c, j);
        d2 += d * d * inv_2s2[j];
      }
      phi(i, c) = std::exp(-d2);
    }
  }
  return phi;
}

// KLIEP score: mean(log(phi * alpha)), returns -1e30 on failure
static double kliep_score(const arma::mat& phi, const arma::vec& alpha)
{
  arma::vec vals = phi * alpha;
  double s = 0.0;
  int n = (int)vals.n_elem;
  for (int i = 0; i < n; i++) {
    if (vals[i] <= 0.0) return -1e30;
    s += std::log(vals[i]);
  }
  return s / n;
}

// Project alpha: alpha += (1 - b'alpha)*c, clip negatives, normalise by b'alpha
static void project_alpha(arma::vec& alpha, const arma::vec& b, const arma::vec& c_vec)
{
  double ba = arma::dot(b, alpha);
  alpha += (1.0 - ba) * c_vec;
  // clip negatives
  for (arma::uword i = 0; i < alpha.n_elem; i++) {
    if (alpha[i] < 0.0) alpha[i] = 0.0;
  }
  ba = arma::dot(b, alpha);
  if (ba > 0.0) alpha /= ba;
}

// KLIEP optimise alpha (matches densratio:::KLIEP_optimize_alpha)
static arma::vec kliep_optimize(const arma::mat& phi_x1,
                                 const arma::vec& mean_phi_x2)
{
  const int M = (int)phi_x1.n_cols;
  arma::vec b = mean_phi_x2;
  double bb = arma::dot(b, b);
  if (bb < 1e-30) return arma::vec(M, arma::fill::zeros);
  arma::vec c_vec = b / bb;

  arma::vec alpha(M, arma::fill::ones);
  project_alpha(alpha, b, c_vec);
  double score = kliep_score(phi_x1, alpha);

  double epsilons[] = {1000, 100, 10, 1, 0.1, 0.01, 0.001};
  const int max_iter = 100;

  for (int e = 0; e < 7; e++) {
    double eps = epsilons[e];
    for (int it = 0; it < max_iter; it++) {
      // gradient step:  alpha_new = alpha + eps * A^T * (1 / (A * alpha))
      arma::vec Aa = phi_x1 * alpha;
      arma::vec inv_Aa(Aa.n_elem);
      for (arma::uword i = 0; i < Aa.n_elem; i++) {
        inv_Aa[i] = (Aa[i] > 1e-30) ? 1.0 / Aa[i] : 1e30;
      }
      arma::vec alpha_new = alpha + eps * (phi_x1.t() * inv_Aa);
      project_alpha(alpha_new, b, c_vec);
      double score_new = kliep_score(phi_x1, alpha_new);
      if (score_new <= score) break;
      alpha = alpha_new;
      score = score_new;
    }
  }
  return alpha;
}

// CV score for a given sigma_vec: build kernel matrices, run K-fold KLIEP
static double kliep_cv_score(const arma::mat& x1,
                              const arma::mat& x2,
                              const arma::mat& centers,
                              const arma::vec& inv_2s2,
                              const arma::ivec& cv_split,
                              int fold)
{
  arma::mat phi_x1 = compute_phi(x1, centers, inv_2s2);
  arma::mat phi_x2 = compute_phi(x2, centers, inv_2s2);

  // check for degenerate kernels
  if (!phi_x1.is_finite() || !phi_x2.is_finite()) return -1e30;

  arma::vec mean_phi_x2 = arma::mean(phi_x2, 0).t();  // column means
  if (arma::accu(arma::abs(mean_phi_x2)) < 1e-30) return -1e30;

  double total = 0.0;
  for (int k = 1; k <= fold; k++) {
    // extract train/test indices
    arma::uvec train_idx = arma::find(cv_split != k);
    arma::uvec test_idx  = arma::find(cv_split == k);
    if (train_idx.n_elem == 0 || test_idx.n_elem == 0) return -1e30;

    arma::vec alpha = kliep_optimize(phi_x1.rows(train_idx), mean_phi_x2);
    double s = kliep_score(phi_x1.rows(test_idx), alpha);
    if (s <= -1e29) return -1e30;
    total += s;
  }
  return total / fold;
}

// Single-coordinate search: mirrors .search_one() in R
static double search_one_coord(double sigma_init,
                                const arma::mat& x1,
                                const arma::mat& x2,
                                const arma::mat& centers,
                                arma::vec& sigma_vec,  // modified in-place
                                const arma::uvec& coord_idx,
                                const arma::ivec& cv_split,
                                int fold)
{
  auto set_and_score = [&](double s) -> double {
    arma::vec inv_2s2(sigma_vec.n_elem);
    for (arma::uword j = 0; j < sigma_vec.n_elem; j++) {
      double sv = sigma_vec[j];
      inv_2s2[j] = 1.0 / (2.0 * sv * sv);
    }
    // temporarily set the coordinates
    for (arma::uword ci = 0; ci < coord_idx.n_elem; ci++) {
      double tmp = 1.0 / (2.0 * s * s);
      inv_2s2[coord_idx[ci]] = tmp;
    }
    return kliep_cv_score(x1, x2, centers, inv_2s2, cv_split, fold);
  };

  double sigma = sigma_init;
  double score = set_and_score(sigma);
  if (!std::isfinite(score)) score = -1e30;

  double digit_powers[] = {10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001};
  for (int dp = 0; dp < 7; dp++) {
    double step = digit_powers[dp];
    // try increasing
    for (int i = 0; i < 9; i++) {
      double sigma_new = sigma + step;
      if (sigma_new <= 0) break;
      double score_new = set_and_score(sigma_new);
      if (!std::isfinite(score_new) || score_new <= score) break;
      score = score_new;
      sigma = sigma_new;
    }
    // try decreasing
    for (int i = 0; i < 9; i++) {
      double sigma_new = sigma - step;
      if (sigma_new <= 1e-6) break;
      double score_new = set_and_score(sigma_new);
      if (!std::isfinite(score_new) || score_new <= score) break;
      score = score_new;
      sigma = sigma_new;
    }
  }

  // update sigma_vec with the found value
  for (arma::uword ci = 0; ci < coord_idx.n_elem; ci++) {
    sigma_vec[coord_idx[ci]] = sigma;
  }
  return sigma;
}


//' KLIEP weight optimisation for fixed kernel matrices
//'
//' Exposed so that the median-heuristic path can optimise the KLIEP weights
//' without depending on `densratio:::KLIEP_optimize_alpha`.
//'
//' @param phi_x1 (n1 x M) kernel matrix of the numerator sample
//' @param phi_x2 (n2 x M) kernel matrix of the denominator sample
//' @return numeric vector of length M of kernel weights
//' @keywords internal
// [[Rcpp::export]]
arma::vec siv_kliep_optimize_cpp(const arma::mat& phi_x1,
                                 const arma::mat& phi_x2)
{
  arma::vec mean_phi_x2 = arma::mean(phi_x2, 0).t();
  return kliep_optimize(phi_x1, mean_phi_x2);
}


//' Full KLIEP with per-dimension bandwidth selection (all in C++)
//'
//' Runs coordinate-descent sigma search + KLIEP weight optimisation.
//' Drop-in replacement for KLIEP_bw() in R but much faster.
//'
//' @param x1         (n1 x p) numerator sample
//' @param x2         (n2 x p) denominator sample
//' @param centers    (M x p)  kernel centres
//' @param sigma_init (p)      starting bandwidths
//' @param n_y_cols   int      number of Y columns (get sigma_Y)
//' @param cv_split   (n1)     1-based fold label per row of x1, built in R so
//'   that the search is reproducible under set.seed()
//' @param fold       int      CV folds
//' @param n_cd_iter  int      coordinate descent iterations
//' @return List with kernel_weights (vec), sigma (vec), centers (mat)
// [[Rcpp::export]]
Rcpp::List siv_kliep_bw_cpp(const arma::mat& x1,
                         const arma::mat& x2,
                         const arma::mat& centers,
                         arma::vec sigma_init,
                         int n_y_cols,
                         const arma::ivec& cv_split,
                         int fold = 5,
                         int n_cd_iter = 3)
{
  const int p = (int)x1.n_cols;
  arma::vec sigma_vec = sigma_init;

  // The CV split is built in R (see siv_kliep()) rather than with
  // arma::randperm(), which draws from Armadillo's own generator and would
  // therefore ignore set.seed().  Cross-fitted nuisance estimates have to be
  // reproducible, so every random draw in this engine comes from R's RNG.
  const int nx1 = (int)x1.n_rows;
  if ((int)cv_split.n_elem != nx1)
    Rcpp::stop("siv_kliep_bw_cpp: length(cv_split) must equal nrow(x1).");

  // coordinate indices for Y and X
  arma::uvec idx_Y(n_y_cols);
  for (int j = 0; j < n_y_cols; j++) idx_Y[j] = j;
  arma::uvec idx_X(p - n_y_cols);
  for (int j = n_y_cols; j < p; j++) idx_X[j - n_y_cols] = j;

  // coordinate descent
  double best_score = -1e30;
  for (int cd = 0; cd < n_cd_iter; cd++) {
    // optimise sigma_Y
    search_one_coord(sigma_vec[0], x1, x2, centers, sigma_vec,
                     idx_Y, cv_split, fold);

    // optimise sigma_X (if X columns exist)
    if (p > n_y_cols) {
      search_one_coord(sigma_vec[n_y_cols], x1, x2, centers, sigma_vec,
                       idx_X, cv_split, fold);
    }

    // check convergence
    arma::vec inv_2s2(p);
    for (int j = 0; j < p; j++) {
      inv_2s2[j] = 1.0 / (2.0 * sigma_vec[j] * sigma_vec[j]);
    }
    double cur = kliep_cv_score(x1, x2, centers, inv_2s2, cv_split, fold);
    if (cur <= best_score) break;
    best_score = cur;
  }

  // final weight optimisation with chosen sigma
  arma::vec inv_2s2(p);
  for (int j = 0; j < p; j++) {
    inv_2s2[j] = 1.0 / (2.0 * sigma_vec[j] * sigma_vec[j]);
  }
  arma::mat phi_x1 = compute_phi(x1, centers, inv_2s2);
  arma::mat phi_x2 = compute_phi(x2, centers, inv_2s2);
  arma::vec mean_phi_x2 = arma::mean(phi_x2, 0).t();
  arma::vec alpha = kliep_optimize(phi_x1, mean_phi_x2);

  return Rcpp::List::create(
    Rcpp::Named("kernel_weights") = Rcpp::wrap(alpha),
    Rcpp::Named("sigma")          = Rcpp::wrap(sigma_vec),
    Rcpp::Named("centers")        = Rcpp::wrap(centers)
  );
}
