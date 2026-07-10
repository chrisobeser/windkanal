#' Estimate personalized effects (mixed-effects boosting, GPBoost)
#'
#' Pruefling 15: gradient tree boosting combined with therapist
#' random intercepts (Sigrist, 2022, JMLR 23(232)) -- the
#' mixed-effects boosting branch of the Li et al. (2026) estimator
#' family. S-learner-style response surface: the boosting part
#' learns f(features, z), the GPModel carries the grouped random
#' effects. Individual effects are counterfactual differences of the
#' fixed-effect predictions (the random intercept cancels within a
#' therapist).
#'
#' Hyperparameters mirror the S-learner (xgboost) for comparability:
#' `nrounds = 300`, `learning_rate = 0.1`, `max_depth = 3`,
#' documented deliberately instead of tuned per world.
#'
#' @param stream Data stream from [sim_stream()].
#' @param nrounds Boosting rounds (default 300).
#' @param learning_rate Learning rate (default 0.1).
#' @param max_depth Tree depth (default 3).
#' @param B Bootstrap replicates for the ATE interval (default 200;
#'   therapist-cluster bootstrap as for the other learners).
#' @param seed Random seed (mandatory).
#' @return `data.frame` (`patient_id`, `x`, `z`, `tau_hat`) with
#'   attribute `ate` (`estimate`, `se`, `lo`, `hi`; cluster
#'   bootstrap).
#' @export
fit_cate_gpboost <- function(stream, nrounds = 300, learning_rate = 0.1,
                             max_depth = 3, B = 200, seed) {
  if (!requireNamespace("gpboost", quietly = TRUE)) {
    stop("fit_cate_gpboost() requires the 'gpboost' package.",
         call. = FALSE)
  }
  if (missing(seed)) stop("`seed` is mandatory.", call. = FALSE)
  set.seed(seed)
  p <- patients(stream)
  tau_fun <- function(d) {
    X <- cbind(cate_features(d), z = d$z)
    gp <- gpboost::GPModel(group_data = d$therapist_id,
                           likelihood = "gaussian")
    bst <- gpboost::gpboost(
      data = X, label = d$score_mean, gp_model = gp,
      nrounds = nrounds, verbose = 0,
      params = list(learning_rate = learning_rate,
                    max_depth = max_depth,
                    num_threads = 1))  # wie xgboost nthread=1:
                    # LightGBM-OpenMP x parallele Worker segfaultet sonst
    X1 <- X; X1[, "z"] <- 1
    X0 <- X; X0[, "z"] <- 0
    pr <- function(Xn) stats::predict(
      bst, data = Xn, group_data_pred = d$therapist_id,
      pred_latent = TRUE)$fixed_effect
    pr(X1) - pr(X0)
  }
  tau_hat <- tau_fun(p)
  out <- data.frame(patient_id = p$patient_id, x = p$x, z = p$z,
                    tau_hat = tau_hat)
  boot <- cluster_boot_ate(p, tau_fun, B)
  attr(out, "ate") <- c(estimate = mean(tau_hat), boot)
  out
}
