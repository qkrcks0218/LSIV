# Logit-Separable Instrumental Variable Approach

This repository contains the `SIV` R package implementing the methodologies in [Park and Tchetgen Tchetgen (2026)](https://arxiv.org/abs/2602.02265). The package provides nonparametric inference for the average and quantile treatment effects on the treated under a logit-separable binary treatment choice model, together with a falsification test of the model's observed data implications.

The package is currently in development. 

## Installation

To install the package from a local copy, set your working directory to the package source directory containing `DESCRIPTION` and run the commands below.

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_local(".", dependencies = TRUE, upgrade = "never")
library(SIV)
```

The package requires R 4.0.0 or later and a C++17 compiler for source installation. On Windows, install the [Rtools version matching your R installation](https://cran.r-project.org/bin/windows/Rtools/).

## Methodology

Consider independent observations $O=(Y,A,Z,X)$, where $`A\in\{0,1\}`$ is a binary treatment indicator, $`Z\in\{0,1\}`$ is a binary instrumental variable (IV), $X$ is a vector of pretreatment covariates, and $Y$ is the observed outcome. Let $Y^{(a)}$ be the potential outcome under treatment $a$. The primary target is the **average treatment effect on the treated (ATT)**,

$$
\tau=E\left[Y^{(1)}-Y^{(0)}\mid A=1\right].
$$

Under consistency, $Y=Y^{(A)}$, the treated mean $E[Y^{(1)}\mid A=1]=E[Y\mid A=1]$ is directly identified. The main challenge is identifying and estimating the counterfactual mean $E[Y^{(0)}\mid A=1]$ in the presence of unmeasured confounding.

### Logit-Separable Treatment Mechanism

In addition to consistency, common support, exclusion, IV independence $Y^{(0)}\perp Z\mid X$, and relevance, the framework assumes a logit-separable treatment mechanism. Define the treatment odds conditional on the treatment-free potential outcome by

$$
\pi(y,z,x)=\frac{\Pr(A=1\mid Y^{(0)}=y,Z=z,X=x)}{\Pr(A=0\mid Y^{(0)}=y,Z=z,X=x)}.
$$

The key restriction is

$$
\pi(y,z,x)=\alpha(y,x)\beta(z,x),\qquad \alpha(y_R,x)=1,
$$

where $y_R$ is a reference outcome value in the common support. Here, $\beta(z,x)=\pi(y_R,z,x)$ is the baseline odds and $\alpha(y,x)=\pi(y,z,x)/\pi(y_R,z,x)$ is the generalized odds ratio relating treatment to the treatment-free potential outcome.

Equivalently,

$$
logit \Pr(A=1\mid Y^{(0)}=y,Z=z,X=x)=\log\alpha(y,x)+\log\beta(z,x).
$$

Thus, the instrument shifts the treatment log odds by an amount that can depend on $X$ but is constant across values of $Y^{(0)}$. The function $\alpha$ measures the degree of unmeasured confounding on the odds ratio scale and is shared across instrument levels. Both components are estimated without imposing a finite-dimensional model for their functional forms.

### ATT Identification

The method uses nuisance functions defined from the observed data law,

```math
\theta=\left\{f(Z\mid X),\ f(A\mid Z,X),\ f(Y\mid A=0,Z,X)\right\}.
```

Under the identifying assumptions, these functions determine $\alpha$ and $\beta$ through a normalized fixed-point equation. The identification and convergence results require the compatibility and regularity conditions stated in Park and Tchetgen Tchetgen (2026).

Once $\alpha$ is recovered, the counterfactual regression

$$
\mu(z,x)=E[Y^{(0)}\mid A=1,Z=z,X=x]
$$

is identified by

$$
\mu(z,x)=\frac{E[Y\alpha(Y,X)\mid A=0,Z=z,X=x]}{E[\alpha(Y,X)\mid A=0,Z=z,X=x]}.
$$

The ATT is then identified as

$$
\tau=E[Y-\mu(Z,X)\mid A=1].
$$

### Efficient Influence Function

Let $\omega(Y^{(0)},Z,X)$ be an auxiliary function satisfying

$$
E[\omega\mid Y^{(0)},X]=E[\omega\mid Z,X]=0,
$$

and write $\mu_\omega(z,x)=E[\omega(Y^{(0)},Z,X)\mid A=1,Z=z,X=x]$. The function $\omega$ solves

```math
E\left[\{Y^{(0)}-\mu(Z,X)\}-\{\omega(Y^{(0)},Z,X)-\mu_\omega(Z,X)\}\mid A=1,Y^{(0)},X\right]=0.
```

For a binary instrument, this integral equation has a unique solution under the conditions in the paper and admits a closed-form representation. Define

```math
\begin{aligned}
\nu(O)={}&(1-A)\pi(Y,Z,X)\{Y-\mu(Z,X)\}+A\mu(Z,X)\\
&-(1-A)\pi(Y,Z,X)\{\omega(Y,Z,X)-\mu_\omega(Z,X)\}\\
&-A\mu_\omega(Z,X)-(1-A)\omega(Y,Z,X).
\end{aligned}
```

The **efficient influence function (EIF)** for the ATT is

$$
IF(O)=\frac{AY-\nu(O)-A\tau}{\Pr(A=1)}.
$$

For one cross-fitting repetition, the estimator is

```math
\widehat\tau=\frac{\mathbb{P}_N\{AY-\widehat\nu(O)\}}{\mathbb{P}_N(A)},
```

where $\mathbb{P}_N$ denotes the empirical average and each observation's nuisance estimates are fitted using the other folds. Under the paper's admissibility, boundedness, relevance, and nuisance estimation rate conditions, the estimator is asymptotically normal and attains the semiparametric efficiency bound. Its influence function supplies standard errors and Wald confidence intervals.

### Nuisance Function Estimation

The default procedure uses **two-fold cross-fitting**. Nuisance functions are estimated on the training observations and evaluated on the held-out fold. `SIV()` repeats this procedure three times by default and reports the median point estimate. For the ATT, the variance adjustment combines the repetition-specific variance estimates with the squared deviations from the median estimate. The numbers of repetitions and folds are controlled by `Overall.Estimation_num.crossfit` and `Overall.Estimation_num.fold`.

**Instrument and treatment probabilities.** The conditional probabilities $\Pr(Z=1\mid X)$ and $\Pr(A=1\mid Z,X)$ are estimated with the Super Learner wrapper `MySL()`. Available learner groups include GLM, lasso/ridge, earth, GAM, XGBoost, polynomial splines, random forests, gradient boosting machines, and a one-layer neural network. `Overall.Estimation_SL.list` selects the groups; the default is `1`, corresponding to GLM.

**Continuous outcomes.** The default estimates two density components.

1. **Baseline conditional density** $f(Y\mid A=0,Z=0,X)$. With `outcome_density = "np_gss_cv"`, kernel conditional density estimates from `np` and smoothing-spline conditional density estimates from `gss` are combined using nonnegative weights summing to one. Inner validation chooses the weights using a quadratic density score, the package's mean squared prediction error (MSPE) criterion. The `"np_normal"` and `"np_cv"` options use only `np`, with normal-reference or likelihood cross-validated bandwidths, respectively.

2. **Conditional density ratio** $R_Y(y,x)=f(y\mid A=0,Z=1,X=x)/f(y\mid A=0,Z=0,X=x)$. With `ratio_method = "kliep_ulsif"`, KLIEP and uLSIF estimate density ratios, reconstruct the numerator conditional density from the baseline density, and combine the reconstructed densities using inner-validation MSPE weights. The default `ratio_bandwidth = "cv_separate"` selects separate outcome and standardized-covariate bandwidths by cross-validation. Standalone `"kliep"`, `"ulsif"`, and `"logistic"` ratio methods are also available.

The recovered conditional densities are normalized on an outcome grid before the fixed-point and EIF calculations. A bound-preserving adjustment enforces the model's necessary density-ratio inequality while maintaining normalization. This adjustment does not establish that an arbitrary fitted nuisance triple is compatible with the causal model.

**Binary and discrete outcomes.** Conditional mass functions are estimated using `MySL()`, and their ratio is computed directly. For binary outcomes, the outcome model is fitted among all untreated observations with $Z$ included as a predictor, and the fixed point has a closed-form solution. Set `y_type = "binary"`, `"discrete"`, or `"continuous"` explicitly for the corresponding outcome type.

### Quantile Treatment Effects on the Treated

For $q\in(0,1)$, the **quantile treatment effect on the treated (QTT)** is

$$
\tau_q=F_1^{-1}(q)-F_0^{-1}(q),\qquad F_a(t)=\Pr(Y^{(a)}\leq t\mid A=1),
$$

where $`F_a^{-1}(q)=\inf\{t\in\mathbb{R}\mid F_a(t)\geq q\}`$. Include `"QTT"` in `estimand` and select quantile indices with `QTT_quantile`.

The procedure combines a binomial order-statistic confidence interval for the observed treated quantile with counterfactual CDF tests, following the confidence-set maximization principle of Berger and Boos (1994) and the counterfactual moment approach of Lee et al. (2025). For continuous outcomes, it tests $F_0(t)=q$. For discrete outcomes, it tests both necessary restrictions $F_0(t-)\leq q\leq F_0(t)$, allowing for a jump across $q$. Accepted counterfactual quantiles are projected through the treated-quantile interval to construct a QTT interval.

The default adjustment is `QTT_c1 = 0.02 * (1 - level)`, giving 0.001 for a 95% interval. `QTT_support` supplies a complete discrete counterfactual support, while `QTT_bounds` supplies known bounds for both potential outcomes. Binary outcomes automatically use bounds `c(0, 1)`. When the acceptance set is empty, the package projects the specified support or bounds; without either, the fallback interval is unbounded.

With repeated cross-fitting, the reported lower and upper QTT confidence limits are the respective medians of the repetition-specific endpoints. Use `confint()` to extract these intervals. QTT inference uses test inversion, and the corresponding entries returned by `siv_se()` are `NA`.

### Falsification Test

The logit-separable IV model implies the observable inequality

$$
\frac{f(y,A=0\mid Z=1,x)-f(y,A=0\mid Z=0,x)}{\Pr(A=0\mid Z=1,x)-\Pr(A=0\mid Z=0,x)}\geq0
$$

for every outcome and covariate value in the common support. Setting `falsification = TRUE` tests this implication using cross-fitted influence function estimates of integrated moments and a Gaussian multiplier bootstrap. The test uses the unrestricted nuisance estimates before the density-ratio adjustment.

In practice, the moment dictionary uses outcome bins and covariate cells formed from one covariate at a time. `falsification_num.bin` controls their counts, `falsification_categorical` identifies categorical covariates, and `falsification_num.boot` controls the bootstrap draws. Rejection provides evidence against the model's necessary implication; non-rejection does not establish the IV assumptions. With repeated cross-fitting, the implementation reports the median of the repetition-specific p-values. Set `Overall.Estimation_num.crossfit = 1` to obtain a single repetition's test.

## Example usage

The following self-contained example generates data satisfying logit separability and requests the ATT, three QTTs, and the falsification test in a single fit. The continuous-outcome density ensembles use cross-validation, so this example can take several minutes.

```r
library(SIV)

###########################################
# Super Learner hyperparameters
###########################################

SL.list <- 1
# 1: GLM
# 2: lasso/ridge
# 3: earth
# 4: GAM
# 5: XGBoost
# 6: polynomial spline
# 7: random forest
# 8: gradient boosting machine
# 9: one-layer neural network

###########################################
# Data
###########################################

simulate_data <- function(n = 600, binary = FALSE) {
  X <- cbind(X1 = runif(n, -1, 1))
  W <- rowSums(X)
  Z <- rbinom(n, 1, plogis(0.1 * W))

  if (binary) {
    Y0 <- rbinom(n, 1, plogis(-0.25 + 0.2 * W))
    Y1 <- rbinom(n, 1, plogis(0.25 + 0.2 * W))
  } else {
    Y0 <- runif(n, -1, 1)
    Y1 <- Y0 + 1
  }

  A <- rbinom(n, 1, plogis(-1.5 + 3 * Z + 0.5 * tanh(Y0) + 0.1 * W))
  Y <- ifelse(A == 1, Y1, Y0)
  list(Y = Y, A = A, Z = Z, X = X)
}

set.seed(2026)
Data <- simulate_data()

###########################################
# SIV: continuous outcomes
###########################################

set.seed(2027)
SIV.Result <- SIV(
  Data$Y, Data$A, Data$Z, Data$X,
  y_type = "continuous",
  estimand = c("ATT", "QTT"),
  Overall.Estimation_SL.list = SL.list,
  Overall.Estimation_num.crossfit = 1,
  Overall.Estimation_num.fold = 2,
  outcome_density = "np_gss_cv",
  ratio_method = "kliep_ulsif",
  ratio_bandwidth = "cv_separate",
  y_grid = seq(-1, 1, length.out = 101),
  y_ref = 0,
  QTT_quantile = c(0.25, 0.5, 0.75),
  QTT_bounds = c(-1, 2),
  falsification = TRUE,
  falsification_num.boot = 399
)

coef(SIV.Result)
#       ATT  QTT_0.25   QTT_0.5  QTT_0.75
# 1.0005798 0.9137543 0.9933969 1.1540166

siv_se(SIV.Result)                # ATT standard error; NA for QTTs
confint(SIV.Result)
#              2.5 %   97.5 %
# ATT      0.7622993 1.238860
# QTT_0.25 0.5016478 1.319589
# QTT_0.5  0.4942908 1.480655
# QTT_0.75 0.4368687 1.574523

summary(SIV.Result)

###########################################
# QTT and falsification results
###########################################

SIV.Result$qtt$estimate
SIV.Result$qtt$ci
SIV.Result$falsification$p_value
# [1] 0.935

SIV.Result$falsification$reject
# [1] FALSE
```

The continuous example has $Y^{(0)}\in[-1,1]$ and $Y^{(1)}=Y^{(0)}+1$, so the true ATT and QTTs are all one. The counterfactual integration grid uses its known support, while `QTT_bounds = c(-1, 2)` supplies bounds containing both potential outcomes. With an unknown support, omit `y_grid` and `y_ref` to construct the grid and select the reference value automatically.

This example uses one cross-fitting repetition. To use the default median adjustment, set `Overall.Estimation_num.crossfit = 3`. `Overall.Estimation_SL.list = 1:9` includes all nine learner groups. The `g_floor_frac` and `min_relevance` entries of `diagnostics` summarize fixed-point flooring and instrument relevance; inspect them together with convergence and the fitted nuisance functions.

For binary outcomes, the probability models and closed-form fixed point are selected automatically.

```r
set.seed(2028)
Data.Binary <- simulate_data(binary = TRUE)

set.seed(2029)
SIV.Result.Binary <- SIV(
  Data.Binary$Y, Data.Binary$A, Data.Binary$Z, Data.Binary$X,
  y_type = "binary",
  Overall.Estimation_SL.list = SL.list,
  Overall.Estimation_num.crossfit = 3,
  Overall.Estimation_num.fold = 2
)

coef(SIV.Result.Binary)
siv_se(SIV.Result.Binary)
confint(SIV.Result.Binary)
```

For a discrete outcome with known support, supply the full support to QTT inversion. For example, an outcome measuring weeks worked with potential outcomes bounded between 0 and 52 uses `QTT_support = 0:52` and `QTT_bounds = c(0, 52)`. These arguments also allow discrete QTT inversion when the nuisance density model is fitted with `y_type = "continuous"`.

### Additional Estimators and Replication Files

The package also provides comparison estimators.

| Function | Method |
|---|---|
| `SIV_Liu()` | IPW, outcome-regression, and doubly robust estimators of [Liu et al. (2020)](https://doi.org/10.5705/ss.202017.0196) |
| `SIV_Sun()` | Doubly robust estimator of [Sun, Miao, and Wickramarachchi](https://doi.org/10.5705/ss.202023.0383), mapped to the ATT |
| `SIV_Sun2021()` | Instrumental variable estimator for nonignorable missing data of [Sun et al. (2018)](https://doi.org/10.5705/ss.202016.0324) |
| `SIV_LATE()` | Influence function estimator of the local average treatment effect under its corresponding IV assumptions |

`SIV_LATE()` targets the effect among compliers. Comparisons with the ATT require interpreting the difference in target populations and assumptions.

The [simulation scripts](ReplicationFiles/Simulation) and [Job Corps analysis](ReplicationFiles/Data) are included in `ReplicationFiles/`. The application directory contains `JC.csv` and scripts for fitting the methods, summarizing density ratios, and aggregating results. Batch-fitting scripts take a batch number as a command-line argument and read supporting files from their own directory.

## References

Park, C., and Tchetgen Tchetgen, E. J. (2026). **Nonparametric Inference with an Instrumental Variable under a Separable Binary Treatment Choice Model.** *arXiv*. [Paper](https://arxiv.org/abs/2602.02265).

Berger, R. L., and Boos, D. D. (1994). P Values Maximized Over a Confidence Set for the Nuisance Parameter. *Journal of the American Statistical Association*, 89, 1012-1016. [Article](https://doi.org/10.1080/01621459.1994.10476836).

Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., and Robins, J. (2018). Double/debiased machine learning for treatment and structural parameters. *The Econometrics Journal*, 21, C1-C68. [Article](https://doi.org/10.1111/ectj.12097).

Gu, C. (2014). Smoothing Spline ANOVA Models: R Package gss. *Journal of Statistical Software*, 58(5), 1-25. [Article](https://doi.org/10.18637/jss.v058.i05).

Hall, P., Racine, J. S., and Li, Q. (2004). Cross-validation and the estimation of conditional probability densities. *Journal of the American Statistical Association*, 99, 1015-1026. [Article](https://doi.org/10.1198/016214504000000548).

Hayfield, T., and Racine, J. S. (2008). Nonparametric Econometrics: The np Package. *Journal of Statistical Software*, 27(5), 1-32. [Article](https://doi.org/10.18637/jss.v027.i05).

Kanamori, T., Hido, S., and Sugiyama, M. (2009). A Least-Squares Approach to Direct Importance Estimation. *Journal of Machine Learning Research*, 10, 1391-1445. [Article](https://www.jmlr.org/papers/v10/kanamori09a.html).

Lee, Y., Yu, M., Liu, J., Park, C., Zhang, Y., Robins, J. M., and Tchetgen Tchetgen, E. J. (2025). Inference on Nonlinear Counterfactual Functionals under a Multiplicative IV Model. *arXiv*. [Paper](https://arxiv.org/abs/2507.15612).

Sugiyama, M., Nakajima, S., Kashima, H., von Buenau, P., and Kawanabe, M. (2007). Direct Importance Estimation with Model Selection and Its Application to Covariate Shift Adaptation. *Advances in Neural Information Processing Systems*, 20. [Paper](https://papers.nips.cc/paper/2007/hash/be83ab3ecd0db773eb2dc1b0a17836a1-Abstract.html).

van der Laan, M. J., Polley, E. C., and Hubbard, A. E. (2007). Super Learner. *Statistical Applications in Genetics and Molecular Biology*, 6(1). [Article](https://doi.org/10.2202/1544-6115.1309).

## AI Disclosure

Parts of the SIV package and this README were prepared with assistance from Claude (Anthropic) and ChatGPT (OpenAI).
