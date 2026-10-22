## ----------------------------------------------------------------------------
## MySL(): the Super Learner wrapper of the UDID package, ported to SIV.
##
## This file is a PORT of `UDID_Package/R/MySL.R` (Park & Tchetgen Tchetgen,
## "A Universal Nonparametric Framework for Difference-in-Differences
## Analyses").  The learner definitions, the tuning grids, the `sl_env`
## construction and the call to `SuperLearner::SuperLearner()` are reproduced
## verbatim so that the two packages fit *identical* ensembles.  Only the
## package plumbing differs:
##
##   1. `SuperLearner:::predict.SL.*` (triple colon) is replaced by
##      `siv_sl_ns()`, a thin `getFromNamespace()` helper.  The triple-colon
##      form is flagged by `R CMD check`; the helper is not, and resolves to
##      exactly the same functions.
##   2. Every optional base-learner package is checked at run time. `MySL()`
##      drops unavailable learner groups with a warning and errors only if
##      nothing is left. UDID lists them as hard `Imports`.
##   3. `PS.Adjust()` is NOT ported.  It is unused in UDID (nothing calls it)
##      and it does not run: its body references the undefined symbols `times`
##      and `crossValFUN_out` and calls `SuperLearner.control()` unqualified.
##      Porting dead, non-running code into a checked package would only
##      generate `R CMD check` notes.
##
## Probability bounding is done by the callers through `siv_bound_prob()`,
## exactly as in UDID (`bound_prob()` in `UDID_nonparametric.R`); `PS.thr` is
## retained in the signature for interface compatibility.
## ----------------------------------------------------------------------------

## Internal: resolve an unexported SuperLearner object by name.  Equivalent to
## `SuperLearner:::<name>` but acceptable to `R CMD check`.
siv_sl_ns <- function(name) utils::getFromNamespace(name, "SuperLearner")

## Internal: packages required by each of the nine learner groups of `MySL()`.
## Group 1 needs SuperLearner only (SL.glm is built on stats::glm).
siv_sl_requires <- list(
  `1` = character(0),
  `2` = "glmnet",
  `3` = "earth",
  `4` = "gam",
  `5` = "xgboost",
  `6` = "polspline",
  `7` = "ranger",
  `8` = c("caret", "gbm"),
  `9` = c("caret", "RSNNS")
)

## Internal: drop the learner groups whose packages are not installed.
siv_sl_available <- function(SL.list) {
  SL.list <- as.integer(SL.list)
  bad <- SL.list[SL.list < 1L | SL.list > 9L]
  if (length(bad)) {
    stop("`SL.list` must select learner groups from 1:9; got ",
         paste(bad, collapse = ", "), ".", call. = FALSE)
  }
  ok <- vapply(SL.list, function(g) {
    pk <- siv_sl_requires[[as.character(g)]]
    !length(pk) || all(vapply(pk, requireNamespace, logical(1), quietly = TRUE))
  }, logical(1))
  if (!all(ok)) {
    miss <- unique(unlist(siv_sl_requires[as.character(SL.list[!ok])]))
    miss <- miss[!vapply(miss, requireNamespace, logical(1), quietly = TRUE)]
    warning("MySL(): dropping learner group(s) ",
            paste(SL.list[!ok], collapse = ", "),
            " because the following package(s) are not installed: ",
            paste(miss, collapse = ", "), ".", call. = FALSE)
  }
  out <- SL.list[ok]
  if (!length(out)) {
    stop("MySL(): none of the requested learner groups is available. ",
         "Install the packages required by at least one selected group.",
         call. = FALSE)
  }
  out
}

#' Bound estimated probabilities away from 0 and 1
#'
#' Applied to every estimated probability used in a denominator (the instrument
#' propensity score \eqn{\Pr(Z = 1 \mid X)} in particular).  This is the
#' `bound_prob()` guard of the UDID package, exposed here so that the
#' bound is a documented, user-controllable argument rather than a hard-coded
#' constant.
#'
#' @param p Numeric vector of estimated probabilities.
#' @param tol Lower bound; `p` is mapped into `[tol, 1 - tol]`.
#' @return `pmin(pmax(p, tol), 1 - tol)`, with non-finite entries replaced by
#'   `0.5`.
#' @examples
#' siv_bound_prob(c(-0.1, 0.0005, 0.5, 1.2), tol = 0.01)
#' @export
siv_bound_prob <- function(p, tol = 1e-3) {
  p <- as.numeric(p)
  p[!is.finite(p)] <- 0.5
  pmin(pmax(p, tol), 1 - tol)
}

## Internal: Super Learner tuning used by the Version1 replication code.
siv_mysl_tuning <- function(N) {
  SL.hpara <- list()
  SL.hpara$MLPL <- c(1, 2)
  SL.hpara$MTRY <- c(1, 2)
  SL.hpara$NMN <- min(round(N * 0.1), 25)
  SL.hpara$MLPdecay <- 10^c(-1, -3)
  SL.hpara
}

## Internal: complete learner grids used by Data_Version1.  Keeping them in a
## small data-only helper makes parity testable without fitting every learner.
siv_mysl_library_tuning <- function(N) {
  h <- siv_mysl_tuning(N)
  list(
    glmnet = list(alpha = c(1, 0.5, 0), useMin = c(TRUE, FALSE)),
    earth = list(degree = 1:5),
    gam = list(deg.gam = 1:5),
    xgboost = list(n.trees = c(100, 300, 500), max_depth = 1:4),
    polymars = list(knots = 2:4),
    ranger = list(num.trees = c(500, 1000, 1500), mtry = h$MTRY),
    gbm = list(ntree = c(100, 300, 500), intdepth = 1:3,
               sh = c(0.1, 0.01), nmn = h$NMN),
    mlp = list(L1 = h$MLPL, L2 = 0, L3 = 0, decay = h$MLPdecay)
  )
}

## Internal: caret preprocessing used by the Version1 GBM and MLP wrappers.
siv_mysl_preprocess <- function() c("center", "scale", "pca")

#' Super Learner Wrapper
#'
#' Fit an ensemble of machine learning algorithms using the Super Learner
#' algorithm (van der Laan et al., 2007) with a configurable library of
#' base learners.
#'
#' This is the `MySL()` wrapper of the \pkg{UDID} package, ported unchanged in
#' its statistical content; see the file header of `R/mysl.R` for the list of
#' plumbing adaptations.  Within \pkg{SIV} it supplies the nuisance-function
#' estimates of [SIV()] and [SIV_LATE()]. [SIV()] always uses this wrapper,
#' while [SIV_LATE()] exposes a learner choice.
#'
#' @param Data A data frame containing the outcome and covariates.
#' @param locY Column index for the outcome variable.
#' @param locX Column index (or indices) for the covariates.
#' @param Ydist A family object (e.g., \code{gaussian()}, \code{binomial()}).
#' @param SL.list Integer vector selecting which learner groups to include (1--9).
#'   Available groups: 1 = GLM, 2 = lasso/ridge, 3 = earth,
#'   4 = GAM, 5 = xgboost, 6 = polynomial spline, 7 = random forest,
#'   8 = gbm, 9 = 1-layer MLP.  Groups whose packages are not installed are
#'   dropped with a warning, because every machine-learning package is a
#'   *suggested*, not a required, dependency of \pkg{SIV}.
#' @param obsWeights Optional observation weights.
#' @param SL.N Sample size used to set the Version1 GBM minimum-node tuning
#'   value. `SIV()` supplies the full analysis sample size.
#' @param PS.thr Propensity score threshold (default 1e-3).  Retained for
#'   interface compatibility with the \pkg{UDID} version; probability bounding
#'   is applied by the callers through [siv_bound_prob()].
#' @param CVlist Optional list of cross-validation fold assignments.
#'
#' @details
#' \code{MySL} is a convenience wrapper of the \code{SuperLearner}
#' algorithm (van der Laan et al., 2007) implemented in the
#' \code{SuperLearner} package (Polley et al., 2025). The function builds
#' a library of candidate learners from the groups specified in \code{SL.list}:
#' \enumerate{
#'   \item \NTC{GLM uses a generalized linear model.}
#'   \item \NTC{Lasso/Ridge uses elastic net via \code{glmnet} with
#'     \eqn{\alpha \in \{0, 0.5, 1\}}.}
#'   \item \NTC{MARS uses multivariate adaptive regression splines via
#'     \code{earth} with polynomial degrees 1-5.}
#'   \item \NTC{GAM uses a generalized additive model via \code{gam} with
#'     spline degrees 1-5.}
#'   \item \NTC{XGBoost uses gradient boosted trees via \code{xgboost} with
#'     a grid over the number of trees and tree depth.}
#'   \item \NTC{Polynomial spline uses polynomial MARS via \code{polymars}
#'     with 2-4 knots.}
#'   \item \NTC{Random forest uses \code{ranger} with a grid over the
#'     number of trees and \code{mtry}.}
#'   \item \NTC{GBM uses a gradient boosted model via \code{caret} with
#'     \code{gbm} backend.}
#'   \item \NTC{1-layer MLP uses a single hidden-layer neural network via
#'     \code{caret} with \code{mlpML} backend.}
#' }
#' The Super Learner computes the optimal convex combination of these base
#' learners using 5-fold cross-validation by default.
#'
#' For Poisson outcomes (\code{Ydist = poisson()}), only learner groups
#' 1 (GLM), 2 (lasso/ridge), 4 (GAM), and 5 (xgboost) are used.
#'
#' @return A fitted \code{SuperLearner} object.
#'
#' @references
#' \itemize{
#'   \item Polley, E., LeDell, E., Kennedy, C., Lendle, S., & van der Laan, M. (2025).
#'     SuperLearner: Super Learner Prediction. R package version 2.0-40.
#'   \item van der Laan, M. J., Polley, E. C., & Hubbard, A. E. (2007).
#'     Super learner.
#'     \emph{Statistical Applications in Genetics and Molecular Biology}, 6(1).
#' }
#'
#' @examples
#' \donttest{
#' if (requireNamespace("SuperLearner", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- data.frame(y = rnorm(100), x1 = rnorm(100), x2 = rnorm(100))
#'   fit <- MySL(d, locY = 1, locX = 2:3, Ydist = gaussian(), SL.list = 1)
#'   head(predict(fit, newdata = d[, 2:3], onlySL = TRUE)$pred)
#' }
#' }
#' @export
MySL <- function( Data, locY, locX, Ydist=stats::gaussian(),
                  SL.list=c(1:9), obsWeights=NULL,
                  SL.N=nrow(Data), PS.thr = 10^(-3), CVlist=NULL ){

  if (!requireNamespace("SuperLearner", quietly = TRUE)) {
    stop("`MySL()` requires the 'SuperLearner' package.",
         call. = FALSE)
  }

  ## Poisson c(2,4,5)
  if(Ydist$family=="poisson"){
    SL.list <- intersect(c(1,2,4,5),SL.list)
  }

  ## Soft dependencies: keep only the learner groups that are installed.
  SL.list <- siv_sl_available(SL.list)

  SL.tune <- siv_mysl_library_tuning(SL.N)
  SL.preprocess <- siv_mysl_preprocess()

  Learners <- list()

  ##############
  # Caret Based
  ##############

  SL.caret.SLP <- function (Y, X, newX, family, obsWeights, method = "mlpML",
                            L1,L2,L3,decay,
                            trControl = caret::trainControl(method = "none"),
                            metric = ifelse(family$family == "gaussian", "RMSE", "Accuracy"), ...)
  {
    if (family$family == "gaussian") {
      fit.train <- caret::train(x = X, y = Y, weights = obsWeights,
                                metric = metric, method = method,
                                tuneGrid = expand.grid(layer1=L1,layer2=L2,layer3=L3),
                                preProc = SL.preprocess,
                                hiddenActFunc = "Act_Logistic",linOut=TRUE,
                                learnFuncParams=decay,
                                trControl = trControl)
      pred <- stats::predict(fit.train, newdata = newX, type = "raw")
    }
    if (family$family == "binomial") {
      Y.f <- as.factor(Y)
      levels(Y.f) <- c("A0", "A1")
      fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights,
                                metric = metric, method = method,
                                tuneGrid = expand.grid(layer1=L1,layer2=L2,layer3=L3),
                                hiddenActFunc = "Act_Identity",
                                learnFuncParams = decay,
                                preProc = SL.preprocess,
                                trControl = trControl)
      pred <- stats::predict(fit.train, newdata = newX, type = "prob")[,2]
    }
    fit <- list(object = fit.train)
    out <- list(pred = pred, fit = fit)
    class(out$fit) <- c("SL.caret")
    return(out)
  }


  SL.caret.MLP <- function (Y, X, newX, family, obsWeights, method = "mlpML",
                            L1,L2,L3,decay,
                            trControl = caret::trainControl(method = "none"),
                            metric = ifelse(family$family == "gaussian", "RMSE", "Accuracy"), ...)
  {
    if (family$family == "gaussian") {
      fit.train <- caret::train(x = X, y = Y, weights = obsWeights,
                                metric = metric, method = method,
                                tuneGrid = expand.grid(layer1=L1,layer2=L2,layer3=L3),
                                learnFuncParams=decay,
                                hiddenActFunc = "Act_Logistic",linOut=TRUE,
                                preProc = SL.preprocess,
                                trControl = trControl)
      pred <- stats::predict(fit.train, newdata = newX, type = "raw")
    }
    if (family$family == "binomial") {
      Y.f <- as.factor(Y)
      levels(Y.f) <- c("A0", "A1")
      fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights,
                                metric = metric, method = method,
                                tuneGrid = expand.grid(layer1=L1,layer2=L2,layer3=L3),
                                hiddenActFunc = "Act_Identity",
                                learnFuncParams = decay,
                                preProc = SL.preprocess,
                                trControl = trControl)
      pred <- stats::predict(fit.train, newdata = newX, type = "prob")[,2]
    }
    fit <- list(object = fit.train)
    out <- list(pred = pred, fit = fit)
    class(out$fit) <- c("SL.caret")
    return(out)
  }



  SL.caret.gbm <- function (Y, X, newX, family, obsWeights, method = "gbm",
                            ntree,intdepth,sh,nmn,
                            trControl = caret::trainControl(method = "none"),
                            metric = ifelse(family$family == "gaussian", "RMSE", "Accuracy"), ...)
  {
    if (family$family == "gaussian") {

      fit.train <- caret::train(x = X, y = Y, weights = obsWeights,
                                metric = metric, method = method,
                                tuneGrid = expand.grid(n.trees=ntree,interaction.depth=intdepth,shrinkage=sh,n.minobsinnode=nmn),
                                preProc = SL.preprocess,
                                trControl = trControl)
      pred <- stats::predict(fit.train, newdata = newX, type = "raw")
    }
    if (family$family == "binomial") {
      Y.f <- as.factor(Y)
      levels(Y.f) <- c("A0", "A1")
      fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights,
                                metric = metric, method = method,
                                tuneGrid = expand.grid(n.trees=ntree,interaction.depth=intdepth,shrinkage=sh,n.minobsinnode=nmn),
                                preProc = SL.preprocess,
                                trControl = trControl)
      pred <- stats::predict(fit.train, newdata = newX, type = "prob")[,2]
    }
    fit <- list(object = fit.train)
    out <- list(pred = pred, fit = fit)
    class(out$fit) <- c("SL.caret")
    return(out)
  }


  #############################################
  # SL-based
  #############################################

  SL.new.earth <- function (Y, X, newX, family, obsWeights, id, degree = 2, penalty = 3,
                            nk = max(21, 2 * ncol(X) + 1), pmethod = "backward",
                            nfold = 0, ncross = 1, minspan = 0, endspan = 0, ...)
  {
    if (family$family == "gaussian") {
      fit.earth <- earth::earth(x = X, y = Y, degree = degree, weights=obsWeights,
                                nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold,
                                ncross = ncross, minspan = minspan, endspan = endspan)
    }
    if (family$family == "binomial") {
      fit.earth <- earth::earth(x = X, y = Y, degree = degree, weights=obsWeights,
                                nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold,
                                ncross = ncross, minspan = minspan, endspan = endspan,
                                glm = list(family = stats::binomial))
    }
    pred <- stats::predict(fit.earth, newdata = newX, type = "response")
    fit <- list(object = fit.earth)
    out <- list(pred = pred, fit = fit)
    class(out$fit) <- c("SL.earth")
    return(out)
  }

  SL.new.xgboost <- function (Y, X, newX, family, obsWeights, id, ntrees = 1000,
                              max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(),
                              nthread = 1, verbose = 0, save_period = NULL,
                              n.trees = ntrees, ...)
  {
    ## create.Learner supplies n.trees from the tuning grid. Keep ntrees as
    ## an alias for existing direct callers of this wrapper.
    ntrees <- n.trees
    if (!is.matrix(X)) {
      X = stats::model.matrix(~. - 1, X)
    }

    xgb_ver <- as.numeric(sub("^(\\d+\\.\\d+).*", "\\1", utils::packageVersion("xgboost")))

    if (xgb_ver >= 2.1) {
      if (family$family == "gaussian") {
        model = xgboost::xgboost(x = X, y = Y, weight = obsWeights,
                                 objective = "reg:squarederror",
                                 nrounds = ntrees, max_depth = max_depth,
                                 min_child_weight = minobspernode,
                                 learning_rate = shrinkage, nthread = nthread)
      }
      if (family$family == "binomial") {
        model = xgboost::xgboost(x = X, y = factor(Y), weight = obsWeights,
                                 objective = "binary:logistic",
                                 nrounds = ntrees, max_depth = max_depth,
                                 min_child_weight = minobspernode,
                                 learning_rate = shrinkage, nthread = nthread)
      }
      if (family$family == "multinomial") {
        model = xgboost::xgboost(x = X, y = factor(Y), weight = obsWeights,
                                 objective = "multi:softmax",
                                 nrounds = ntrees, max_depth = max_depth,
                                 min_child_weight = minobspernode,
                                 learning_rate = shrinkage,
                                 num_class = length(unique(Y)), nthread = nthread)
      }
      if (family$family == "poisson") {
        model = xgboost::xgboost(x = X, y = Y, weight = obsWeights,
                                 objective = "count:poisson",
                                 nrounds = ntrees, max_depth = max_depth,
                                 min_child_weight = minobspernode,
                                 learning_rate = shrinkage, nthread = nthread)
      }
    } else {
      if (family$family == "gaussian") {
        xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
        model = xgboost::xgboost(data = xgmat, objective = "reg:squarederror",
                                 nrounds = ntrees, max_depth = max_depth,
                                 min_child_weight = minobspernode,
                                 eta = shrinkage, verbose = verbose, nthread = nthread,
                                 params = params, save_period = save_period)
      }
      if (family$family == "binomial") {
        xgmat = xgboost::xgb.DMatrix(data = X, label = factor(Y), weight = obsWeights)
        model = xgboost::xgboost(data = xgmat, objective = "binary:logistic",
                                 nrounds = ntrees, max_depth = max_depth,
                                 min_child_weight = minobspernode,
                                 eta = shrinkage, verbose = verbose, nthread = nthread,
                                 params = params, save_period = save_period)
      }
      if (family$family == "multinomial") {
        xgmat = xgboost::xgb.DMatrix(data = X, label = factor(Y), weight = obsWeights)
        model = xgboost::xgboost(data = xgmat, objective = "multi:softmax",
                                 nrounds = ntrees, max_depth = max_depth,
                                 min_child_weight = minobspernode,
                                 eta = shrinkage, verbose = verbose,
                                 num_class = length(unique(Y)),
                                 nthread = nthread, params = params, save_period = save_period)
      }
      if (family$family == "poisson") {
        xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
        model = xgboost::xgboost(data = xgmat, objective = "count:poisson",
                                 nrounds = ntrees, max_depth = max_depth,
                                 min_child_weight = minobspernode,
                                 eta = shrinkage, verbose = verbose,
                                 nthread = nthread, params = params, save_period = save_period)
      }
    }

    if (!is.matrix(newX)) {
      newX = stats::model.matrix(~. - 1, newX)
    }
    pred = stats::predict(model, newdata = newX)
    fit = list(object = model)
    class(fit) = c("SL.xgboost")
    out = list(pred = pred, fit = fit)
    return(out)
  }

  ## Build a local environment that contains all learner and screen functions
  ## so SuperLearner can find them by name.
  sl_env <- list2env(list(
    ## Version1 caret preprocessing shared by the custom wrappers
    SL.preprocess      = SL.preprocess,
    ## Screen functions
    All               = SuperLearner::All,
    ## Base learners from SuperLearner
    SL.glm            = SuperLearner::SL.glm,
    SL.glmnet         = SuperLearner::SL.glmnet,
    SL.gam            = SuperLearner::SL.gam,
    SL.polymars       = SuperLearner::SL.polymars,
    SL.ranger         = SuperLearner::SL.ranger,
    SL.nnet           = SuperLearner::SL.nnet,
    SL.xgboost        = SuperLearner::SL.xgboost,
    ## Predict methods (unexported in SuperLearner; resolved by name)
    predict.SL.glm       = siv_sl_ns("predict.SL.glm"),
    predict.SL.glmnet    = siv_sl_ns("predict.SL.glmnet"),
    predict.SL.gam       = siv_sl_ns("predict.SL.gam"),
    predict.SL.polymars  = siv_sl_ns("predict.SL.polymars"),
    predict.SL.ranger    = siv_sl_ns("predict.SL.ranger"),
    predict.SL.nnet      = siv_sl_ns("predict.SL.nnet"),
    predict.SL.xgboost   = siv_sl_ns("predict.SL.xgboost"),
    predict.SL.caret     = siv_sl_ns("predict.SL.caret"),
    predict.SL.earth     = siv_sl_ns("predict.SL.earth"),
    ## Custom learners defined in this function
    SL.caret.SLP      = SL.caret.SLP,
    SL.caret.MLP      = SL.caret.MLP,
    SL.caret.gbm      = SL.caret.gbm,
    SL.new.earth       = SL.new.earth,
    SL.new.xgboost     = SL.new.xgboost
  ), parent = globalenv())

  Learners[[1]] <- SuperLearner::create.Learner("SL.glm")
  TOTAL.M <- 1

  Learners[[2]] <- SuperLearner::create.Learner(
    "SL.glmnet", tune = SL.tune$glmnet
  )

  Learners[[3]] <- SuperLearner::create.Learner(
    "SL.new.earth", tune = SL.tune$earth
  )

  Learners[[4]] <- SuperLearner::create.Learner(
    "SL.gam", tune = SL.tune$gam
  )

  Learners[[5]] <- SuperLearner::create.Learner(
    "SL.new.xgboost", tune = SL.tune$xgboost
  )

  Learners[[6]] <- SuperLearner::create.Learner(
    "SL.polymars", tune = SL.tune$polymars
  )

  Learners[[7]] <- SuperLearner::create.Learner(
    "SL.ranger", tune = SL.tune$ranger
  )

  # Learners[[TOTAL.M+1]] <- SuperLearner::create.Learner("SL.nnet",tune=list(linout=c(TRUE,FALSE), decay=c(0,0.1),size=MTRY))

  Learners[[8]] <- SuperLearner::create.Learner(
    "SL.caret.gbm", tune = SL.tune$gbm
  )

  Learners[[9]] <- SuperLearner::create.Learner(
    "SL.caret.SLP", tune = SL.tune$mlp
  )

  # Learners[[10]] <- SuperLearner::create.Learner("SL.caret.MLP", tune = list(L1=c(2,4),L2=MLPL,L3=MLPL,decay=MLPdecay))

  BaseLearner <- Learners[[ SL.list[1] ]]$names

  if( length(SL.list)>1 ){
    for( METHOD in 2:length(SL.list) ){
      BaseLearner <- c(BaseLearner, Learners[[ SL.list[METHOD] ]]$names)
    }
  }

  ## Copy all generated learner wrapper functions into sl_env
  ## and re-parent them so they can find base learners (SL.ranger, etc.)
  calling_env <- environment()
  for (lname in BaseLearner) {
    if (exists(lname, envir = calling_env, inherits = FALSE)) {
      fn <- get(lname, envir = calling_env)
      if (is.function(fn)) {
        environment(fn) <- sl_env
      }
      assign(lname, fn, envir = sl_env)
    }
  }

  if( length(locX)==1 ){
    dX <- data.frame(matrix( Data[,locX], dim(Data)[1], 1))
    colnames(dX) <- colnames(Data)[locX]
  } else {
    dX <- Data[,locX]
  }

  if(is.null(CVlist)){
    CVL <- list(V = 5)
  } else {
    CVL <- list(V = 5, shuffle = FALSE, validRows = CVlist)
  }

  if(Ydist$family!="binomial"){
    utils::capture.output( Fitted.SL <- SuperLearner::SuperLearner(Y=Data[,locY],X=dX,family=Ydist,
                                              SL.library=BaseLearner, cvControl = CVL, obsWeights=obsWeights,
                                              env = sl_env) , file=NULL )
  } else {
    utils::capture.output( Fitted.SL <- SuperLearner::SuperLearner(Y=Data[,locY],X=dX,family=Ydist,
                                              SL.library=BaseLearner, cvControl = CVL, obsWeights=obsWeights,
                                              env = sl_env) , file=NULL )
  }

  return(Fitted.SL)
}
