#' Therapist-cluster bootstrap for learner-based ATE inference
#'
#' Meta-learners and boosting provide **no native uncertainty
#' quantification** -- a property of the tool class, not an
#' implementation detail. Following the practice in the
#' methods-comparison literature (Li et al., 2026), we resample
#' whole therapists with replacement and refit.
#'
#' @param p Person table from [patients()].
#' @param tau_fun Function(p) -> tau_hat vector for that sample.
#' @param B Number of bootstrap replicates.
#' @return c(se, lo, hi) of the ATE across replicates.
#' @keywords internal
cluster_boot_ate <- function(p, tau_fun, B) {
  if (B <= 0) {  # explizit ohne Inferenz (z. B. teure Lerner im Gitter)
    return(c(se = NA_real_, lo = NA_real_, hi = NA_real_))
  }
  ids <- unique(p$therapist_id)
  ates <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    ziehung <- sample(ids, length(ids), replace = TRUE)
    pb <- do.call(rbind, lapply(ziehung, function(id)
      p[p$therapist_id == id, , drop = FALSE]))
    if (length(unique(pb$z)) < 2) next  # draw without both arms
    ates[b] <- mean(tau_fun(pb))
  }
  ates <- ates[!is.na(ates)]
  if (length(ates) < 10) {
    warning("cluster bootstrap: fewer than 10 valid draws ",
            "(both arms missing in resamples); returning NA.")
    return(c(se = NA_real_, lo = NA_real_, hi = NA_real_))
  }
  c(se = stats::sd(ates),
    lo = unname(stats::quantile(ates, 0.025)),
    hi = unname(stats::quantile(ates, 0.975)))
}

#' Estimate personalized effects (T-learner, random forest)
#'
#' The class of tools the field is currently being taught (see the
#' meta-learner tutorial literature): fit one outcome forest on the
#' treated, one on the controls, and take the difference of their
#' predictions as the personal effect. No shrinkage discipline, no
#' native inference -- ATE uncertainty comes from a therapist-cluster
#' bootstrap ([cluster_boot_ate()]).
#'
#' @param stream Data stream from [sim_stream()].
#' @param num_trees Trees per forest (default 500).
#' @param B Bootstrap replicates for the ATE interval (default 200).
#' @param seed Random seed (mandatory).
#' @return `data.frame` (`patient_id`, `x`, `z`, `tau_hat`) with
#'   attribute `ate` (`estimate`, `se`, `lo`, `hi`; cluster
#'   bootstrap).
#' @export
fit_cate_tlearner <- function(stream, num_trees = 500, B = 200, seed) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("fit_cate_tlearner() requires the 'ranger' package.",
         call. = FALSE)
  }
  if (missing(seed)) {
    stop("`seed` is mandatory (forests are stochastic).", call. = FALSE)
  }
  set.seed(seed)
  p <- patients(stream)
  feats <- colnames(cate_features(p))
  tau_fun <- function(d) {
    f <- stats::as.formula(paste("score_mean ~",
                                 paste(feats, collapse = " + ")))
    m1 <- ranger::ranger(f, data = d[d$z == 1, ],
                         num.trees = num_trees, seed = seed)
    m0 <- ranger::ranger(f, data = d[d$z == 0, ],
                         num.trees = num_trees, seed = seed)
    stats::predict(m1, d)$predictions -
      stats::predict(m0, d)$predictions
  }
  tau_hat <- tau_fun(p)
  out <- data.frame(patient_id = p$patient_id, x = p$x, z = p$z,
                    tau_hat = tau_hat)
  boot <- cluster_boot_ate(p, tau_fun, B)
  attr(out, "ate") <- c(estimate = mean(tau_hat), boot)
  out
}

#' Estimate personalized effects (S-learner, gradient boosting)
#'
#' The CatBoost/XGBoost class -- what applied matching studies
#' actually deploy: one boosted model of the outcome given features
#' plus treatment; the personal effect is the prediction difference
#' between z = 1 and z = 0. Field-typical defaults, deterministic
#' (`subsample = 1`, one thread). No shrinkage discipline toward
#' homogeneity, no native inference (therapist-cluster bootstrap).
#'
#' @param stream Data stream from [sim_stream()].
#' @param nrounds Boosting rounds (default 300).
#' @param max_depth Tree depth (default 3).
#' @param eta Learning rate (default 0.1).
#' @param B Bootstrap replicates for the ATE interval (default 200).
#' @param seed Random seed (mandatory).
#' @return `data.frame` (`patient_id`, `x`, `z`, `tau_hat`) with
#'   attribute `ate` (`estimate`, `se`, `lo`, `hi`; cluster
#'   bootstrap).
#' @export
fit_cate_sboost <- function(stream, nrounds = 300, max_depth = 3,
                            eta = 0.1, B = 200, seed) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("fit_cate_sboost() requires the 'xgboost' package.",
         call. = FALSE)
  }
  if (missing(seed)) {
    stop("`seed` is mandatory.", call. = FALSE)
  }
  set.seed(seed)
  p <- patients(stream)
  tau_fun <- function(d) {
    X <- cbind(cate_features(d), z = d$z)
    m <- suppressWarnings(xgboost::xgboost(
      data = X, label = d$score_mean, nrounds = nrounds,
      max_depth = max_depth, eta = eta, subsample = 1,
      nthread = 1, verbose = 0))
    X1 <- X; X1[, "z"] <- 1
    X0 <- X; X0[, "z"] <- 0
    stats::predict(m, X1) - stats::predict(m, X0)
  }
  tau_hat <- tau_fun(p)
  out <- data.frame(patient_id = p$patient_id, x = p$x, z = p$z,
                    tau_hat = tau_hat)
  boot <- cluster_boot_ate(p, tau_fun, B)
  attr(out, "ate") <- c(estimate = mean(tau_hat), boot)
  out
}

#' Estimate personalized effects (X-learner, random forest)
#'
#' Kuenzel et al.'s X-learner (full text in quellen/): T-learner
#' outcome models first, then *imputed* individual effects
#' (D1 = Y - m0(x) for the treated, D0 = m1(x) - Y for controls)
#' are themselves modeled by forests, and the two effect models are
#' blended with propensity weights. Designed to outperform the
#' T-learner under unbalanced arms. Same bootstrap inference as the
#' other learners ([cluster_boot_ate()]).
#'
#' @inheritParams fit_cate_tlearner
#' @return As [fit_cate_tlearner()].
#' @export
fit_cate_xlearner <- function(stream, num_trees = 500, B = 200, seed) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("fit_cate_xlearner() requires the 'ranger' package.",
         call. = FALSE)
  }
  if (missing(seed)) {
    stop("`seed` is mandatory (forests are stochastic).", call. = FALSE)
  }
  set.seed(seed)
  p <- patients(stream)
  feats <- colnames(cate_features(p))
  f1 <- stats::as.formula(paste("score_mean ~",
                                paste(feats, collapse = " + ")))
  f2 <- stats::as.formula(paste(".D ~", paste(feats, collapse = " + ")))
  tau_fun <- function(d) {
    d1 <- d[d$z == 1, ]; d0 <- d[d$z == 0, ]
    m1 <- ranger::ranger(f1, data = d1, num.trees = num_trees,
                         seed = seed)
    m0 <- ranger::ranger(f1, data = d0, num.trees = num_trees,
                         seed = seed)
    d1$.D <- d1$score_mean - stats::predict(m0, d1)$predictions
    d0$.D <- stats::predict(m1, d0)$predictions - d0$score_mean
    t1 <- ranger::ranger(f2, data = d1, num.trees = num_trees,
                         seed = seed)
    t0 <- ranger::ranger(f2, data = d0, num.trees = num_trees,
                         seed = seed)
    g <- mean(d$z)  # propensity weight (randomized setting)
    g * stats::predict(t0, d)$predictions +
      (1 - g) * stats::predict(t1, d)$predictions
  }
  tau_hat <- tau_fun(p)
  out <- data.frame(patient_id = p$patient_id, x = p$x, z = p$z,
                    tau_hat = tau_hat)
  attr(out, "ate") <- c(estimate = mean(tau_hat),
                        cluster_boot_ate(p, tau_fun, B))
  out
}

#' Estimate personalized effects (DR-learner, random forest)
#'
#' The doubly robust learner as taught to the field (see the
#' meta-learner tutorial; canonical treatment: Kennedy -- reference
#' to be verified before citing in prose): nuisance outcome models
#' feed a doubly robust pseudo-outcome, which is then regressed on
#' the features. Own-row overfitting of the pseudo-outcome is
#' mitigated by using out-of-bag predictions where a row was part of
#' the training data (a lightweight stand-in for full cross-fitting,
#' documented as such). Propensity is `mean(z)` (randomized setting).
#'
#' @inheritParams fit_cate_tlearner
#' @return As [fit_cate_tlearner()].
#' @export
fit_cate_drlearner <- function(stream, num_trees = 500, B = 200, seed) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("fit_cate_drlearner() requires the 'ranger' package.",
         call. = FALSE)
  }
  if (missing(seed)) {
    stop("`seed` is mandatory (forests are stochastic).", call. = FALSE)
  }
  set.seed(seed)
  p <- patients(stream)
  feats <- colnames(cate_features(p))
  f1 <- stats::as.formula(paste("score_mean ~",
                                paste(feats, collapse = " + ")))
  f3 <- stats::as.formula(paste(".psi ~", paste(feats, collapse = " + ")))
  tau_fun <- function(d) {
    i1 <- d$z == 1
    m1 <- ranger::ranger(f1, data = d[i1, ], num.trees = num_trees,
                         seed = seed)
    m0 <- ranger::ranger(f1, data = d[!i1, ], num.trees = num_trees,
                         seed = seed)
    n <- nrow(d)
    mu1 <- numeric(n); mu0 <- numeric(n)
    mu1[i1]  <- m1$predictions                               # OOB
    mu1[!i1] <- stats::predict(m1, d[!i1, ])$predictions
    mu0[!i1] <- m0$predictions                               # OOB
    mu0[i1]  <- stats::predict(m0, d[i1, ])$predictions
    e <- mean(d$z)
    mz <- ifelse(i1, mu1, mu0)
    d$.psi <- (d$z - e) / (e * (1 - e)) * (d$score_mean - mz) +
      mu1 - mu0
    td <- ranger::ranger(f3, data = d, num.trees = num_trees,
                         seed = seed)
    td$predictions  # OOB fuer Trainingszeilen
  }
  tau_hat <- tau_fun(p)
  out <- data.frame(patient_id = p$patient_id, x = p$x, z = p$z,
                    tau_hat = tau_hat)
  attr(out, "ate") <- c(estimate = mean(tau_hat),
                        cluster_boot_ate(p, tau_fun, B))
  out
}

#' Estimate personalized effects (R-learner, random forest)
#'
#' Nie & Wager's R-learner (full text in quellen/): local centering
#' of the outcome (out-of-bag), then the treatment-residual-weighted
#' regression of the pseudo-outcome -- the orthogonalization idea
#' whose importance Dandl et al. and Tabib & Larocque independently
#' confirm. Propensity is `mean(z)` (randomized setting), so the
#' treatment residual is z - mean(z).
#'
#' @inheritParams fit_cate_tlearner
#' @return As [fit_cate_tlearner()].
#' @export
fit_cate_rlearner <- function(stream, num_trees = 500, B = 200, seed) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("fit_cate_rlearner() requires the 'ranger' package.",
         call. = FALSE)
  }
  if (missing(seed)) {
    stop("`seed` is mandatory (forests are stochastic).", call. = FALSE)
  }
  set.seed(seed)
  p <- patients(stream)
  feats <- colnames(cate_features(p))
  f1 <- stats::as.formula(paste("score_mean ~",
                                paste(feats, collapse = " + ")))
  f4 <- stats::as.formula(paste(".po ~", paste(feats, collapse = " + ")))
  tau_fun <- function(d) {
    m <- ranger::ranger(f1, data = d, num.trees = num_trees,
                        seed = seed)
    ytil <- d$score_mean - m$predictions   # OOB-zentriert
    e <- mean(d$z)
    ztil <- d$z - e
    d$.po <- ytil / ztil
    tr <- ranger::ranger(f4, data = d, num.trees = num_trees,
                         seed = seed, case.weights = ztil^2)
    tr$predictions  # OOB
  }
  tau_hat <- tau_fun(p)
  out <- data.frame(patient_id = p$patient_id, x = p$x, z = p$z,
                    tau_hat = tau_hat)
  attr(out, "ate") <- c(estimate = mean(tau_hat),
                        cluster_boot_ate(p, tau_fun, B))
  out
}

#' Estimate personalized effects (PAI: linear T-learner)
#'
#' The field's own historical method -- the Personalized Advantage
#' Index tradition (DeRubeis school): one *linear* regression per
#' arm, and the difference of their predictions is the personal
#' effect. Implemented in its simplest field-typical form (all
#' candidate features as main effects, no variable selection --
#' selection variants can be added later). Deterministic given the
#' data; the seed feeds the bootstrap ([cluster_boot_ate()]).
#'
#' @inheritParams fit_cate_tlearner
#' @return As [fit_cate_tlearner()].
#' @export
fit_cate_pai <- function(stream, B = 200, seed) {
  if (missing(seed)) {
    stop("`seed` is mandatory (bootstrap is stochastic).", call. = FALSE)
  }
  set.seed(seed)
  p <- patients(stream)
  feats <- colnames(cate_features(p))
  f <- stats::as.formula(paste("score_mean ~",
                               paste(feats, collapse = " + ")))
  tau_fun <- function(d) {
    m1 <- stats::lm(f, data = d[d$z == 1, ])
    m0 <- stats::lm(f, data = d[d$z == 0, ])
    stats::predict(m1, d) - stats::predict(m0, d)
  }
  tau_hat <- tau_fun(p)
  out <- data.frame(patient_id = p$patient_id, x = p$x, z = p$z,
                    tau_hat = tau_hat)
  attr(out, "ate") <- c(estimate = mean(tau_hat),
                        cluster_boot_ate(p, tau_fun, B))
  out
}

#' Estimate personalized effects (model-based forest, MOB)
#'
#' The second tree philosophy (Zeileis/Hothorn/Seibold tradition;
#' the hybrid family that performed best overall in Tabib &
#' Larocque and anchors Dandl et al., full texts in quellen/):
#' instead of splitting on effect heterogeneity directly, each leaf
#' carries a *model* (here `lm(score_mean ~ z)`), and the forest of
#' personalized models yields a per-person treatment coefficient
#' (model4you::pmforest / pmodel). Bootstrap inference as in the
#' other learners.
#'
#' @inheritParams fit_cate_tlearner
#' @param num_trees Trees in the personalized-model forest
#'   (default 200 -- pmforest is costlier per tree than ranger).
#' @export
fit_cate_mob <- function(stream, num_trees = 200, B = 200, seed) {
  if (!requireNamespace("model4you", quietly = TRUE)) {
    stop("fit_cate_mob() requires the 'model4you' package.",
         call. = FALSE)
  }
  if (missing(seed)) {
    stop("`seed` is mandatory (forests are stochastic).", call. = FALSE)
  }
  set.seed(seed)
  p <- patients(stream)
  feats <- colnames(cate_features(p))
  tau_fun <- function(d) {
    dd <- d[, c("score_mean", "z", feats)]
    base <- stats::lm(score_mean ~ z, data = dd)
    frst <- model4you::pmforest(base, data = dd, ntree = num_trees,
                                perturb = list(replace = FALSE,
                                               fraction = 0.632))
    co <- model4you::pmodel(frst, fun = stats::coef)
    unname(co[, "z"])
  }
  tau_hat <- tau_fun(p)
  out <- data.frame(patient_id = p$patient_id, x = p$x, z = p$z,
                    tau_hat = tau_hat)
  attr(out, "ate") <- c(estimate = mean(tau_hat),
                        cluster_boot_ate(p, tau_fun, B))
  out
}
