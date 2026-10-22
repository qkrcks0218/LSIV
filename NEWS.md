# SIV 0.1.0.9000 (development)

* \NTC{The continuous-outcome defaults are now `outcome_density = "np_gss_cv"`
  and `ratio_bandwidth = "cv_separate"`. Both np/gss and KLIEP/uLSIF mixtures
  use MSPE alone, also known as ISE. Their quadratic scores are minimized
  directly on the simplex, including its endpoints, with no penalty or KL
  averaging and no loss-selection option. `mspe_weights` records the
  MSPE weights; `mse_weights` remains an alias for existing analysis code.
  Help text now renders change-marked text and learner names in plain font.}

* \NTC{`QTT_bounds` now supplies known bounds for both potential outcomes.
  QTT inference applies these bounds to exterior treated-quantile ranks,
  accepted counterfactual cells, and empty-inversion fallbacks. For example,
  `QTT_bounds = c(0, 52)` keeps QTT intervals within `[-52, 52]`. Binary
  outcomes automatically use bounds `c(0, 1)`; other outcome types use only
  explicitly supplied bounds.}

* \NTC{QTT inference now handles discrete CDF jumps with two one-sided
  quantile tests and includes complete accepted intervals for continuous
  outcomes. Treated-quantile confidence intervals use binomial order-statistic
  ranks. The default `QTT_c1 = 0.02 * (1 - level)` equals 0.001 at the default
  `level = 0.95`, and remains user configurable. `QTT_support` can specify a
  complete discrete treatment-free support independently of the density model.
  `QTT_numBB` is retained for compatibility and no longer affects inference.
  Empty inversions use specified support or an unbounded interval, and
  infinite confidence limits are preserved across cross-fitting repetitions.}

* \NTC{`SIV()` now uses median-adjusted repeated cross-fitting. The default is
  three independently randomized estimates with two folds each. The point
  estimate is their median, and its variance is the median of each split-specific
  variance plus the squared deviation from the median estimate. The public
  QTT estimate is likewise the median across repetitions. Its lower and upper
  confidence limits are the respective medians of the repetition-specific
  Berger---Boos limits, following Supplement A.8.4. When requested, the
  falsification test likewise reports the median of the repetition-specific
  multiplier-bootstrap p-values and bases its rejection decision on that
  median. The public
  controls are `Overall.Estimation_SL.list`,
  `Overall.Estimation_num.crossfit`, and `Overall.Estimation_num.fold`.
  `y_num.grid = 101L` replaces `ratio_num.grid`, and `y_ref = NULL` selects the
  untreated-outcome mean for continuous outcomes and zero otherwise.}

* \NTC{`SIV()` now has a compact public interface. Outcome, nuisance,
  density-ratio, contraction, QTT, and falsification controls use consistent
  prefixes. The nuisance-input and learner switches, `estimate_ratio`,
  `lower_density`, `QTT.n_grid`, and the joint falsification-cell option were
  removed. Super Learner is always used for the instrument and treatment
  nuisances, and also for binary or discrete outcome mass functions. Their
  density ratio is formed automatically. Continuous outcomes offer
  `"kliep_ulsif"`, `"kliep"`, `"ulsif"`, and `"logistic"`, with bandwidth
  selection controlled separately by `ratio_bandwidth`. QTT inversion now uses
  the main outcome grid. The density and ratio floors both default to `1e-6`.}

* For a continuous outcome, `siv_grid()` now expands the integration grid by
  `expand * sd(Y[A == 0])` rather than by the absolute amount `expand`, so the
  grid range is invariant to the units in which `Y` is recorded.  The default
  `expand = 0.1` is unchanged in value but is now read as a fraction of a
  standard deviation.

* `SIV()` now accepts `estimand = c("ATT", "QTT")`.  QTT point estimates use
  the one-step counterfactual CDF, and confidence intervals use the
  \NTC{Berger---Boos test inversion in Supplement A.6.3.}
* All package dependencies, including `SuperLearner` and its additional learner
  libraries, are required through `Imports`. `Install.R` installs every missing
  hard dependency before doing a clean local source build.
* Source builds are now portable between macOS and Windows: C++17 is declared
  explicitly, compiled artifacts are excluded from source distributions, and
  the bundled installation script uses a clean rebuild before installation.
* **BUG FIX (changes every reported `SIV()` standard error).** `siv_eif()`
  divided by `E{alpha | A=0,z,X}` twice when forming `beta`: it used
  `beta = gamma / denom` where the `gamma` returned by `siv_fixed_point()` is
  eq. (18) and already contains that division.  The manuscript's replication
  code uses `beta = odds / denom` with `odds = Pr(A=1|z,X)/Pr(A=0|z,X)`.  The
  weight `alpha * beta` was therefore too small by a factor of `denom` (range
  1.04-2.31), which broke the identity that `f(Y0=y|X)` reconstructed from the
  `Z=0` and `Z=1` arms must agree (residual 9e-04 before, 1e-25 after).  The
  ATT estimate stayed consistent either way -- the mis-scaled term is a
  mean-zero augmentation -- but the influence function, and hence every
  standard error and confidence interval, was wrong.  `SIV()` now reproduces
  `Simulation_F_Both/1.Sim_Conti_SepIV.R` exactly (ATT to 4e-13, SE to 2e-14).
* `siv_kliep()` / `siv_kliep_predict()` -- the ported engine, exported.

* `SIV_LATE()` -- efficient-influence-function estimator of the local average
  treatment effect, with cross-fitting and influence-function inference.  This
  targets a DIFFERENT estimand from the rest of the package: the complier
  average causal effect under monotonicity, not the ATT under separability.
* `siv_bound_prob()` -- the probability-bounding guard, exposed and documented.
* Compiled fast paths (`src/siv_fastpath.cpp`, Rcpp/RcppArmadillo) for the
  contraction of Theorem 3.2 and for the grid integrals of the influence
  function.  The R reference implementations are retained and are used when the
  compiled code is unavailable; `options(SIV.use_cpp = FALSE)` forces the R
  path.  Measured speedup of `siv_fixed_point()`: 4.2x-5.5x for a continuous
  outcome; of `SIV()` end to end: 2.6x-3.5x.
* Simulation and validation DGPs were removed from the installed package API.
  They now live with the simulation scripts in
  `Simulation_Version2_based_on_SIV_Package/00_common.R`; package tests use
  private, test-only fixtures.

# SIV 0.1.0

* Added `ratio_ensemble.loss`, defaulting to cross-fitted integrated squared
  error (`"ise"`), for data-driven KLIEP/uLSIF mixture weights. The former
  fixed average of the ISE and log-score weight vectors is no longer used.
  The optional `"logscore"` criterion now evaluates the log score of the
  density mixture.

* Repeated-cross-fitting median adjustment now excludes non-finite estimates
  parameter by parameter. Variance adjustment uses only repetitions having
  both a finite estimate and a finite variance, and QTT endpoint medians omit
  non-finite endpoints.

* First release.
* `SIV()` -- the proposed separable-IV estimator of the ATT under the
  logit-separable binary treatment choice model, with an explicit `type`
  argument selecting the binary/discrete or continuous formulation.
* `SIV_Liu()` -- IPW, outcome-regression and doubly robust estimators of
  Liu et al. (2021), with the stacked-estimating-equation sandwich covariance.
* `SIV_Sun()` -- the doubly robust estimator of Sun, Miao and Wickramarachchi
  (2025), mapped to the ATT, with the sandwich covariance.
* `SIV_Sun2021()` -- the MNAR-with-IV estimators of Sun et al. (2021).
* `siv_basis_spec()` -- configurable spline bases for the selection-bias
  (log odds-ratio) function, always normalised at a reference outcome value.
* `siv_sandwich()` -- the stacked-estimating-equation sandwich covariance with
  numerical diagnostics.
