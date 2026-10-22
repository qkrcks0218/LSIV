// siv_fastpath.cpp
//
// Compiled fast paths for the two hot spots of the proposed estimator, both of
// which are dense n x M matrix recursions over an outcome grid.  Profiling
// SIV(type = "continuous", crossfit = 2) at n = 1000 attributed 51% of the
// total run time to siv_fixed_point() (24% of it inside siv_Psi()) and 21% to
// siv_eif(); nothing else exceeded 4%.
//
//   siv_contraction_cpp()    the normalised contraction of Theorem 3.2,
//                            eq. (19).  The R version allocates roughly eight
//                            n x M temporaries per iteration (`g * RY`, the
//                            three products in `num`, `num * S`, `pA1_zM * RY`,
//                            the quotient, `abs(gn - g) * Wm`); at n = 1000 and
//                            M = 578 that is ~37 MB of allocation per sweep,
//                            repeated for every iteration.  The loop below
//                            keeps two n x M buffers for the whole solve and
//                            allocates nothing inside the iteration.
//
//   siv_delta_moments_cpp()  the four delta-weighted grid integrals that
//                            siv_eif() needs for each level of Z.  The R
//                            version calls int_z() six times, each a
//                            rowSums(V * f * W) over the full n x M grid; this
//                            computes all four in a single sweep, reusing
//                            delta[i,m]*fY[i,m]*w[m].
//
// SEMANTICS ARE IDENTICAL TO THE R REFERENCE IMPLEMENTATIONS, which are
// retained in R/siv_core.R and are used whenever the compiled code is not
// available (e.g. when the package sources are sourced directly rather than
// installed).  Numerical agreement is verified in tests/testthat/test-rcpp.R;
// the two differ only by floating-point summation order.
//
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;

//' Contraction solver for the fixed point of Theorem 3.2 (compiled)
//'
//' Compiled equivalent of the `"contraction"` branch of [siv_fixed_point()].
//' Not intended to be called directly; [siv_fixed_point()] dispatches to it
//' when the compiled code is available and falls back to the R implementation
//' otherwise.
//'
//' @param g0 `n x M` matrix of starting values, already normalised to integrate
//'   to one against `weights`.
//' @param pA1_zm,pA1_zM Length-`n` vectors of treatment probabilities at the
//'   minimising / maximising level of the instrument.
//' @param fY_zm,fY_zM `n x M` conditional-density matrices at those levels.
//' @param RY `n x M` density ratio of eq. (15).
//' @param weights Length-`M` vector of integration weights.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance on the weighted L1 change.
//' @return A list with `g` (the fixed point), `iter` and `converged`.
//' @keywords internal
// [[Rcpp::export]]
List siv_contraction_cpp(NumericMatrix g0,
                         NumericVector pA1_zm, NumericVector pA1_zM,
                         NumericMatrix fY_zm, NumericMatrix fY_zM,
                         NumericMatrix RY, NumericVector weights,
                         int maxit, double tol) {
  const int n = g0.nrow(), M = g0.ncol();
  if (fY_zm.nrow() != n || fY_zm.ncol() != M ||
      fY_zM.nrow() != n || fY_zM.ncol() != M ||
      RY.nrow()    != n || RY.ncol()    != M)
    stop("siv_contraction_cpp: non-conformable matrices.");
  if (weights.size() != M) stop("siv_contraction_cpp: length(weights) != ncol.");
  if (pA1_zm.size() != n || pA1_zM.size() != n)
    stop("siv_contraction_cpp: length(pA1_z*) != nrow.");

  NumericMatrix g(clone(g0));
  NumericMatrix gn(n, M);

  const double *w = REAL(weights);
  const double *pzm = REAL(pA1_zm), *pzM = REAL(pA1_zM);
  const double *fzm = REAL(fY_zm), *fzM = REAL(fY_zM), *r = REAL(RY);
  double *gp = REAL(g), *gnp = REAL(gn);

  std::vector<double> S(n), den(n), dif(n);
  std::vector<int> bad(n);
  LogicalVector converged(n, false);

  int iter = 0;
  for (int j = 0; j < maxit; ++j) {
    iter = j + 1;

    // S_i = sum_m g[i,m] * RY[i,m] * w[m]
    std::fill(S.begin(), S.end(), 0.0);
    for (int m = 0; m < M; ++m) {
      const double wm = w[m];
      const std::size_t off = (std::size_t)m * n;
      for (int i = 0; i < n; ++i) S[i] += gp[off + i] * r[off + i] * wm;
    }

    // gn = pmax(Psi(g), 0);  den_i = sum_m gn[i,m] * w[m]
    std::fill(den.begin(), den.end(), 0.0);
    for (int m = 0; m < M; ++m) {
      const double wm = w[m];
      const std::size_t off = (std::size_t)m * n;
      for (int i = 0; i < n; ++i) {
        const std::size_t k = off + i;
        const double num = pzm[i] * gp[k] + (1.0 - pzm[i]) * fzm[k]
                                          - (1.0 - pzM[i]) * fzM[k];
        double v = (num * S[i]) / (pzM[i] * r[k]);
        if (!(v > 0.0)) v = 0.0;          // pmax(., 0); also maps NaN to 0
        gnp[k] = v;
        den[i] += v * wm;
      }
    }

    // Units whose update degenerated keep their previous (already normalised)
    // value and are reported as non-converged -- exactly the R semantics.
    for (int i = 0; i < n; ++i) {
      bad[i] = (!R_finite(den[i]) || den[i] <= 0.0) ? 1 : 0;
      if (bad[i]) den[i] = 1.0;
    }

    // gn <- gn / den ; dif_i = sum_m |gn - g| * w[m] ; g <- gn
    std::fill(dif.begin(), dif.end(), 0.0);
    for (int m = 0; m < M; ++m) {
      const double wm = w[m];
      const std::size_t off = (std::size_t)m * n;
      for (int i = 0; i < n; ++i) {
        const std::size_t k = off + i;
        const double v = bad[i] ? gp[k] : gnp[k] / den[i];
        dif[i] += std::fabs(v - gp[k]) * wm;
        gnp[k] = v;
      }
    }
    std::memcpy(gp, gnp, sizeof(double) * (std::size_t)n * (std::size_t)M);

    bool all_ok = true;
    for (int i = 0; i < n; ++i) {
      converged[i] = (dif[i] < tol) && !bad[i];
      if (!converged[i]) all_ok = false;
    }
    if (all_ok) break;
  }

  return List::create(_["g"] = g, _["iter"] = iter, _["converged"] = converged);
}

//' Delta-weighted grid integrals for the influence function (compiled)
//'
//' Computes, in a single sweep over the `n x M` grid, the four integrals
//' \eqn{\int \delta f}, \eqn{\int y \delta f}, \eqn{\int \delta^2 f} and
//' \eqn{\int y \delta^2 f} against the supplied integration weights.  Compiled
//' equivalent of the repeated `int_z()` calls in [siv_eif()].
//'
//' @param delta `n x M` matrix of odds-ratio values.
//' @param grid Length-`M` outcome grid.
//' @param fY `n x M` conditional-density matrix.
//' @param weights Length-`M` vector of integration weights.
//' @return An `n x 4` matrix with columns
//'   `denom`, `y`, `alpha`, `alphay`.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix siv_delta_moments_cpp(NumericMatrix delta, NumericVector grid,
                                    NumericMatrix fY, NumericVector weights) {
  const int n = delta.nrow(), M = delta.ncol();
  if (fY.nrow() != n || fY.ncol() != M)
    stop("siv_delta_moments_cpp: non-conformable matrices.");
  if (grid.size() != M || weights.size() != M)
    stop("siv_delta_moments_cpp: length(grid)/length(weights) != ncol.");

  NumericMatrix out(n, 4);
  double *o = REAL(out);
  const double *d = REAL(delta), *f = REAL(fY);
  const double *y = REAL(grid), *w = REAL(weights);

  for (int m = 0; m < M; ++m) {
    const double wm = w[m], ym = y[m];
    const std::size_t off = (std::size_t)m * n;
    for (int i = 0; i < n; ++i) {
      const std::size_t k = off + i;
      const double dv = d[k];
      const double base = dv * f[k] * wm;      // delta * f * w
      const double base2 = base * dv;          // delta^2 * f * w
      o[i]             += base;
      o[n + i]         += ym * base;
      o[2 * n + i]     += base2;
      o[3 * n + i]     += ym * base2;
    }
  }
  colnames(out) = CharacterVector::create("denom", "y", "alpha", "alphay");
  return out;
}
